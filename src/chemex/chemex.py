"""The chemex module provides the entry point for the chemex script."""

import os
import sys
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

from chemex.cli import build_parser
from chemex.configuration.method_input import prepare_method_plan
from chemex.configuration.method_plan import (
    FormatOrigin,
    MethodPlan,
    StepPlan,
)
from chemex.configuration.methods import Method, Methods, Selection, read_method_plan
from chemex.configuration.parameters import ParameterConfigurationError, read_defaults
from chemex.containers.experiments import Experiments
from chemex.exceptions import ArtifactPublicationError, ChemExError
from chemex.experiments.builder import build_experiments
from chemex.messages import (
    print_logo,
    print_method_v1_deprecation_warning,
    print_no_data,
    print_reading_defaults,
    print_reading_methods,
    print_running_simulations,
    print_start_fit,
)
from chemex.optimize.fitting import invalidate_planned_outputs
from chemex.optimize.helper import execute_simulation
from chemex.optimize.method_plan_execution import execute_method_plan
from chemex.parameters.parameterization import (
    IncompleteParameterDependenciesError,
)
from chemex.parameters.sealed import InvalidConfigurationError
from chemex.run_info import (
    InputFile,
    capture_input_files,
    mark_failure_stage,
    write_run_info,
)
from chemex.runtime import (
    AnalysisSession,
    ExecutionSettings,
    ensure_plugins_registered,
)

_VERIFIED_PATH_ATTRIBUTES = (
    "diagnostics_path",
    "evidence_path",
    "samples_path",
    "failures_path",
    "restart_path",
    "outcome_path",
)


def _is_interrupted_failure(error: BaseException) -> bool:
    return isinstance(error, KeyboardInterrupt) or (
        isinstance(error, ChemExError)
        and error.__dict__.get("terminal") == "interrupted"
    )


def _copy_verified_paths(source: BaseException, target: BaseException) -> None:
    for name in _VERIFIED_PATH_ATTRIBUTES:
        value = getattr(source, name, None)
        if isinstance(value, Path):
            object.__setattr__(target, name, value)


def run_fit(
    args: Namespace,
    experiments: Experiments,
    session: AnalysisSession,
    *,
    input_files: Sequence[InputFile],
    argv: Sequence[str] | None = None,
    methods: Methods | MethodPlan | None = None,
) -> None:
    if methods is None:
        methods = _read_fit_methods(args)

    parameter_model = session.parameter_factory.sealed_parameter_model
    if parameter_model is None:
        raise RuntimeError("Native parameter model is unavailable")
    plan = prepare_method_plan(methods, parameter_model)
    starting_values = session.analysis_values.snapshot()
    run_info = write_run_info(
        args,
        parameter_model=parameter_model,
        starting_values=starting_values,
        execution=session.execution,
        input_files=input_files,
        argv=argv,
    )

    try:
        try:
            invalidate_planned_outputs(plan, args.output)
        except (Exception, KeyboardInterrupt) as error:
            mark_failure_stage(error, "output")
            raise

        resolved_values = session.resolve_current_values(experiments.param_ids)
        experiments.filter_from_values(resolved_values)

        print_start_fit()
        execute_method_plan(
            experiments,
            plan,
            args.output,
            args.plot,
            session=session,
            run_info=run_info,
        )
        run_info.write_outcome("complete", session.analysis_values.snapshot())
    except (Exception, KeyboardInterrupt) as error:
        mark_failure_stage(error, "deterministic_fit")
        try:
            run_info.write_outcome(
                "incomplete",
                session.analysis_values.snapshot(),
                failure=error,
                failure_stage="deterministic_fit",
            )
            object.__setattr__(error, "outcome_path", run_info.path / "outcome.toml")
            if run_info.restart_revision > 0:
                object.__setattr__(
                    error,
                    "restart_path",
                    run_info.path / "restart.toml",
                )
        except ArtifactPublicationError as finalization_error:
            if _is_interrupted_failure(error):
                _copy_verified_paths(error, finalization_error)
                if run_info.restart_revision > 0:
                    object.__setattr__(
                        finalization_error,
                        "restart_path",
                        run_info.path / "restart.toml",
                    )
                raise
            error.add_note(
                "ChemEx could not publish the incomplete run outcome: "
                f"{type(finalization_error).__name__}: {finalization_error}"
            )
        except (Exception, KeyboardInterrupt) as finalization_error:  # noqa: BLE001
            error.add_note(
                "ChemEx could not publish the incomplete run outcome: "
                f"{type(finalization_error).__name__}: {finalization_error}"
            )
        raise


def run_sim(
    args: Namespace,
    experiments: Experiments,
    session: AnalysisSession,
) -> None:
    print_running_simulations()

    path = args.output
    plot = args.plot == "normal"

    parameter_model = session.parameter_factory.sealed_parameter_model
    if parameter_model is None:
        raise RuntimeError("Native parameter model is unavailable")
    snapshot = session.analysis_values.snapshot()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    resolved_values = parameterization.resolve(
        parameterization.frame_from_snapshot(snapshot)
    )

    try:
        execute_simulation(
            experiments,
            path,
            parameter_values=resolved_values,
            parameter_model=parameter_model,
            parameterization=parameterization,
            plot=plot,
        )
    except (Exception, KeyboardInterrupt) as error:
        mark_failure_stage(error, "simulation")
        raise


def _read_fit_methods(args: Namespace) -> MethodPlan:
    if args.method is None:
        return MethodPlan(FormatOrigin.V1, (StepPlan(""),))
    print_reading_methods()
    plan = read_method_plan(args.method)
    if plan.format_origin is FormatOrigin.V1:
        print_method_v1_deprecation_warning()
    return plan


def run(
    args: Namespace,
    session: AnalysisSession | None = None,
    *,
    argv: Sequence[str] | None = None,
) -> None:
    """Run the fit or simulation."""
    if session is None:
        session = AnalysisSession.create()

    session.execution = ExecutionSettings.from_counts(
        workers=getattr(args, "workers", 1),
        native_threads=getattr(args, "native_threads", "auto"),
    )
    os.environ.update(session.execution.native_thread_env())

    input_files = capture_input_files(args) if args.commands == "fit" else ()

    # Parse kinetics model
    session.set_model(args.model)

    # Read experimental setup and data
    selection = Selection(args.include, args.exclude)
    experiments = build_experiments(args.experiments, selection, session=session)

    if not experiments:
        print_no_data()
        sys.exit()

    methods = _read_fit_methods(args) if args.commands == "fit" else None

    # Read initial values of fitting/fixed parameters
    print_reading_defaults()
    defaults = read_defaults(args.parameters)
    session.parameters.set_defaults(defaults)
    if not session.try_build_analysis_values():
        construction_error = getattr(
            session.parameter_factory,
            "native_construction_error",
            None,
        )
        msg = "Native parameter initialization failed"
        if construction_error is None:
            raise RuntimeError(msg)
        if isinstance(
            construction_error,
            (
                IncompleteParameterDependenciesError,
                InvalidConfigurationError,
            ),
        ):
            raise ParameterConfigurationError(
                tuple(args.parameters),
                str(construction_error),
            ) from construction_error
        raise RuntimeError(f"{msg}: {construction_error}") from construction_error

    if args.commands == "simulate":
        run_sim(args, experiments, session)
    else:
        run_fit(
            args,
            experiments,
            session,
            input_files=input_files,
            argv=argv,
            methods=methods,
        )


def main() -> None:
    """Do all the magic."""
    print_logo()
    ensure_plugins_registered()

    parser = build_parser()
    args = parser.parse_args()
    if args.analysis_command:
        run(args, argv=sys.argv)
    else:
        args.func(args)
