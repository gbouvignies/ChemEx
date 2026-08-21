"""The fitting module contains the code for fitting the experimental data."""

import shutil
from pathlib import Path

from chemex.configuration.methods import Method, Methods, Statistics
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_fitmethod,
    print_grid_statistic_warning,
    print_group_name,
    print_mcmc_no_vary_warning,
    print_minimizing,
    print_no_data,
    print_running_statistics,
    print_step_name,
)
from chemex.optimize.grouping import create_groups
from chemex.optimize.helper import (
    execute_post_fit,
    execute_post_fit_groups,
)
from chemex.optimize.mcmc import run_native_mcmc
from chemex.optimize.minimizer import (
    minimize_with_report,
)
from chemex.optimize.native_deterministic import (
    NativeDeterministicFit,
    run_native_deterministic,
)
from chemex.optimize.resampling import (
    run_native_resampling_statistics,
    run_resampling_statistics,
)
from chemex.runtime import AnalysisSession, ExecutionSettings

_CHEMEX_RESULT_PATHS = (
    "Parameters",
    "Data",
    "Plots",
    "Grid",
    "Groups",
    "All",
    "Statistics",
    "statistics.toml",
)


def invalidate_planned_outputs(methods: Methods, path: Path) -> None:
    """Remove only ChemEx-owned results for every method step planned now."""
    if len(methods) > 1:
        output_root = path.resolve()
        step_roots = tuple(path / section for section in methods)
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


def _run_statistics(
    experiments: Experiments,
    path: Path,
    fitmethod: str,
    statistics: Statistics | None = None,
    *,
    execution: ExecutionSettings | None = None,
) -> None:
    if statistics is None:
        return

    run_resampling_statistics(
        experiments,
        path,
        fitmethod,
        statistics,
        execution=execution,
    )


def _fit_groups(
    experiments: Experiments,
    path: Path,
    plot: str,
    fitmethod: str,
    statistics: Statistics | None,
    *,
    execution: ExecutionSettings | None = None,
) -> None:
    parameter_store = experiments.parameter_store
    groups = create_groups(experiments)

    plot_flg = (plot == "normal" and len(groups) == 1) or plot == "all"

    print_minimizing()

    for group in groups:
        group_lmfit_params = parameter_store.build_lmfit_params(
            group.experiments.param_ids,
        )
        group_path = path / group.path

        if message := group.message:
            print_group_name(message)

        best_lmfit_params = minimize_with_report(
            group.experiments,
            group_lmfit_params,
            fitmethod,
        )

        parameter_store.update_from_parameters(best_lmfit_params)
        residuals = group.experiments.residuals(best_lmfit_params)
        nvarys = sum(
            parameter.vary and not parameter.expr
            for parameter in best_lmfit_params.values()
        )
        execute_post_fit(
            group.experiments,
            group_path,
            plot=plot_flg,
            residuals=residuals,
            nvarys=nvarys,
        )

        # Run Monte Carlo and/or bootstrap analysis
        _run_statistics(
            group.experiments,
            group_path,
            fitmethod,
            statistics,
            execution=execution,
        )

    if len(groups) > 1:
        params_lf = parameter_store.build_lmfit_params(experiments.param_ids)
        execute_post_fit_groups(
            experiments,
            path,
            plot,
            residuals=experiments.residuals(params_lf),
            nvarys=sum(
                parameter.vary and not parameter.expr
                for parameter in params_lf.values()
            ),
        )


def _run_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: Statistics,
    fit: NativeDeterministicFit,
    *,
    session: AnalysisSession,
) -> None:
    committed = session.analysis_values.snapshot()
    try:
        run_native_resampling_statistics(
            experiments,
            path,
            statistics,
            fit,
            execution=session.execution,
        )
        if statistics.mcmc is not None:
            print_running_statistics("MCMC")
            run_native_mcmc(
                experiments,
                fit,
                statistics.mcmc,
                path,
                execution=session.execution,
            )
    finally:
        if session.analysis_values.snapshot() != committed:
            raise RuntimeError("Native statistics mutated the committed central fit")


def _requests_only_mcmc(statistics: Statistics) -> bool:
    return (
        statistics.mcmc is not None
        and statistics.mc is None
        and statistics.bs is None
        and statistics.bsn is None
    )


def _run_requested_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: Statistics | None,
    fit: NativeDeterministicFit | None,
    *,
    session: AnalysisSession,
) -> None:
    if statistics is None:
        return
    if fit is None:
        if _requests_only_mcmc(statistics):
            print_mcmc_no_vary_warning()
            return
        raise RuntimeError("Native resampling requires a committed deterministic fit")
    _run_native_statistics(
        experiments,
        path,
        statistics,
        fit,
        session=session,
    )


def run_methods(
    experiments: Experiments,
    methods: Methods,
    path: Path,
    plot_level: str,
    *,
    session: AnalysisSession,
) -> None:
    for index, (section, method) in enumerate(methods.items(), start=1):
        if section:
            print_step_name(section, index, len(methods))

        # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
        experiments.select(method.selection)

        if not experiments:
            print_no_data()
            continue

        effective_method: Method = method
        if method.grid and method.statistics:
            print_grid_statistic_warning()
            effective_method = method.model_copy(update={"statistics": None})

        print_fitmethod(effective_method.fitmethod)

        # Update the parameter "vary" and "expr" status
        session.apply_current_parameter_roles(
            effective_method,
            experiments.param_ids,
        )

        path_sect = path / section if len(methods) > 1 else path

        fit = run_native_deterministic(
            experiments,
            effective_method,
            path_sect,
            plot_level,
            session=session,
        )
        _run_requested_native_statistics(
            experiments,
            path_sect,
            effective_method.statistics,
            fit,
            session=session,
        )
