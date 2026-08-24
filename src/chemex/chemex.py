"""The chemex module provides the entry point for the chemex script."""

import os
import sys
from argparse import Namespace
from collections.abc import Sequence

from chemex.cli import build_parser
from chemex.configuration.method_plan import (
    FormatOrigin,
    MethodFormatError,
    MethodPlan,
    StepPlan,
)
from chemex.configuration.methods import Methods, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.experiments.builder import build_experiments
from chemex.messages import (
    console,
    print_logo,
    print_no_data,
    print_reading_defaults,
    print_reading_methods,
    print_running_simulations,
    print_start_fit,
)
from chemex.optimize.fitting import invalidate_planned_outputs, run_methods
from chemex.optimize.helper import execute_simulation
from chemex.run_info import write_run_info, write_run_outcome
from chemex.runtime import (
    AnalysisSession,
    ExecutionSettings,
    ensure_plugins_registered,
)


def run_fit(
    args: Namespace,
    experiments: Experiments,
    session: AnalysisSession,
    *,
    argv: Sequence[str] | None = None,
    methods: Methods | MethodPlan | None = None,
) -> None:
    if methods is None:
        methods = _read_fit_methods(args)

    write_run_info(args, experiments, argv=argv)

    try:
        invalidate_planned_outputs(methods, args.output)

        resolved_values = session.resolve_current_values(experiments.param_ids)
        experiments.filter_from_values(resolved_values)

        print_start_fit()
        run_methods(
            experiments,
            methods,
            args.output,
            args.plot,
            session=session,
        )
        write_run_outcome(args.output, "complete")
    except (Exception, KeyboardInterrupt) as error:
        try:
            write_run_outcome(args.output, "incomplete", failure=error)
        except (Exception, KeyboardInterrupt) as outcome_error:  # noqa: BLE001
            error.add_note(
                f"ChemEx could not publish the incomplete run outcome: {outcome_error}"
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

    session.parameters.fix_all()
    resolved_values = session.resolve_current_values(experiments.param_ids)

    execute_simulation(
        experiments,
        path,
        parameter_values=resolved_values,
        plot=plot,
    )


def _read_fit_methods(args: Namespace) -> MethodPlan:
    if args.method is None:
        return MethodPlan(FormatOrigin.V1, (StepPlan(""),))
    print_reading_methods()
    try:
        return read_method_plan(args.method)
    except MethodFormatError as error:
        console.print(f"[red] -- ERROR: {error}")
        sys.exit(1)


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

    # Parse kinetics model
    session.set_model(args.model)

    # Read experimental setup and data
    selection = Selection(args.include, args.exclude)
    experiments = build_experiments(args.experiments, selection, session=session)

    if not experiments:
        print_no_data()
        sys.exit()

    if experiments.parameter_store is not session.parameters:
        msg = "Experiments parameter store does not match the active session"
        raise ValueError(msg)

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
        raise RuntimeError(f"{msg}: {construction_error}") from construction_error

    if isinstance(methods, MethodPlan):
        session.validate_method_plan(methods)

    if args.commands == "simulate":
        run_sim(args, experiments, session)
    else:
        run_fit(args, experiments, session, argv=argv, methods=methods)


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
