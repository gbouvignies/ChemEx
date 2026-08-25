"""The fitting module contains the code for fitting the experimental data."""

import secrets
import shutil
from pathlib import Path

from chemex.configuration.method_execution import (
    mcmc_settings,
    normalize_methods_for_execution,
    step_names,
)
from chemex.configuration.method_plan import MethodPlan, StatisticsPlan
from chemex.configuration.methods import (
    Methods,
    Statistics,
)
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_fitmethod,
    print_mcmc_no_vary_warning,
    print_no_data,
    print_running_statistics,
    print_step_name,
)
from chemex.optimize.mcmc import run_native_mcmc
from chemex.optimize.native_deterministic import (
    NativeDeterministicFit,
    run_native_deterministic,
)
from chemex.optimize.resampling import run_native_resampling_statistics
from chemex.run_info import RunInfo, mark_failure_stage
from chemex.runtime import AnalysisSession

_CHEMEX_RESULT_PATHS = (
    "Parameters",
    "Data",
    "Plots",
    "Grid",
    "Groups",
    "All",
    "Components",
    "Statistics",
    "statistics.toml",
)


def invalidate_planned_outputs(methods: Methods | MethodPlan, path: Path) -> None:
    """Remove only ChemEx-owned results for every method step planned now."""
    names = step_names(methods)
    if len(names) > 1:
        output_root = path.resolve()
        step_roots = tuple(path / section for section in names)
        if any(
            not step_root.resolve().is_relative_to(output_root)
            for step_root in step_roots
        ):
            msg = "A planned method step root is outside the output directory"
            raise ValueError(msg)
    else:
        step_roots = (path,)
    for step_root in step_roots:
        for name in _CHEMEX_RESULT_PATHS:
            result_path = step_root / name
            if result_path.is_symlink() or result_path.is_file():
                result_path.unlink()
            elif result_path.is_dir():
                shutil.rmtree(result_path)


def _run_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: StatisticsPlan,
    fit: NativeDeterministicFit,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
    step_name: str = "DEFAULT",
) -> None:
    committed = session.analysis_values.snapshot()
    try:
        for name, request in (
            ("mc", statistics.mc),
            ("bs", statistics.bs),
            ("bsn", statistics.bsn),
        ):
            if request is None:
                continue
            root_seed = secrets.randbits(64) if request.seed is None else request.seed
            if run_info is not None:
                run_info.record_stochastic_operation(step_name, name, root_seed)
            run_native_resampling_statistics(
                experiments,
                path,
                Statistics(**{name: request.replicates}),
                fit,
                execution=session.execution,
                root_seed=root_seed,
            )
        mcmc = mcmc_settings(statistics.mcmc)
        if mcmc is not None:
            print_running_statistics("MCMC")
            run_native_mcmc(
                experiments,
                fit,
                mcmc,
                path,
                execution=session.execution,
                seed_recorder=(
                    None
                    if run_info is None
                    else lambda seed: run_info.record_stochastic_operation(
                        step_name,
                        "mcmc",
                        seed,
                    )
                ),
            )
    finally:
        if session.analysis_values.snapshot() != committed:
            raise RuntimeError("Native statistics mutated the committed central fit")


def _requests_only_mcmc(statistics: StatisticsPlan) -> bool:
    return (
        statistics.mcmc is not None
        and statistics.mc is None
        and statistics.bs is None
        and statistics.bsn is None
    )


def _run_requested_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: StatisticsPlan | None,
    fit: NativeDeterministicFit | None,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
    step_name: str = "DEFAULT",
) -> None:
    if statistics is None:
        return
    if fit is None:
        if _requests_only_mcmc(statistics):
            print_mcmc_no_vary_warning()
            return
        raise RuntimeError("Native resampling requires a committed deterministic fit")
    try:
        _run_native_statistics(
            experiments,
            path,
            statistics,
            fit,
            session=session,
            run_info=run_info,
            step_name=step_name,
        )
    except (Exception, KeyboardInterrupt) as error:
        mark_failure_stage(error, "statistics")
        raise


def run_methods(
    experiments: Experiments,
    methods: Methods | MethodPlan,
    path: Path,
    plot_level: str,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
) -> None:
    plan, operational = normalize_methods_for_execution(methods)
    effective_actions = plan.effective_role_actions()
    for index, step in enumerate(plan.steps, start=1):
        section = step.name
        method = operational[section]
        if section:
            print_step_name(section, index, len(plan.steps))

        # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
        experiments.select(method.selection)

        if not experiments:
            print_no_data()
            continue

        print_fitmethod(method.fitmethod)

        parameterization = session.compile_parameterization_from_actions(
            effective_actions[section],
            experiments.param_ids,
        )

        path_sect = path / section if len(plan.steps) > 1 else path

        fit = run_native_deterministic(
            experiments,
            method,
            path_sect,
            plot_level,
            session=session,
            parameterization=parameterization,
            search=step.search,
            run_info=run_info,
            step_name=section or "DEFAULT",
        )
        _run_requested_native_statistics(
            experiments,
            path_sect,
            step.statistics,
            fit,
            session=session,
            run_info=run_info,
            step_name=section or "DEFAULT",
        )
