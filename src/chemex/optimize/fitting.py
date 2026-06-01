"""The fitting module contains the code for fitting the experimental data."""

from pathlib import Path

from chemex.configuration.methods import Methods, Statistics
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_calculation_stopped_error,
    print_fitmethod,
    print_grid_statistic_warning,
    print_group_name,
    print_minimizing,
    print_no_data,
    print_running_statistics,
    print_step_name,
    print_value_error,
)
from chemex.optimize.gridding import run_grid
from chemex.optimize.grouping import create_groups
from chemex.optimize.helper import (
    execute_post_fit,
    execute_post_fit_groups,
)
from chemex.optimize.mcmc import run_mcmc
from chemex.optimize.minimizer import (
    minimize_with_report,
)
from chemex.optimize.resampling import run_resampling_statistics
from chemex.runtime import AnalysisSession


def _run_statistics(
    experiments: Experiments,
    path: Path,
    fitmethod: str,
    statistics: Statistics | None = None,
) -> None:
    if statistics is None:
        return

    run_resampling_statistics(experiments, path, fitmethod, statistics)
    parameter_store = experiments.parameter_store

    if statistics.mcmc is None:
        return

    print_running_statistics("MCMC")
    try:
        params_lf = parameter_store.build_lmfit_params(experiments.param_ids)
        run_mcmc(experiments, params_lf, statistics.mcmc, path)
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()


def _fit_groups(
    experiments: Experiments,
    path: Path,
    plot: str,
    fitmethod: str,
    statistics: Statistics | None,
    *,
    session: AnalysisSession | None = None,
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
        execute_post_fit(
            group.experiments,
            group_path,
            plot=plot_flg,
            session=session,
        )

        # Run Monte Carlo and/or bootstrap analysis
        _run_statistics(
            group.experiments,
            group_path,
            fitmethod,
            statistics,
        )

    if len(groups) > 1:
        execute_post_fit_groups(experiments, path, plot, session=session)


def run_methods(
    experiments: Experiments,
    methods: Methods,
    path: Path,
    plot_level: str,
    *,
    session: AnalysisSession | None = None,
) -> None:
    parameter_store = experiments.parameter_store

    for index, (section, method) in enumerate(methods.items(), start=1):
        if section:
            print_step_name(section, index, len(methods))

        # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
        experiments.select(method.selection)

        if not experiments:
            print_no_data()
            continue

        print_fitmethod(method.fitmethod)

        # Update the parameter "vary" and "expr" status
        parameter_store.set_parameter_status(method)

        if method.grid and method.statistics:
            print_grid_statistic_warning()
            method.statistics = None

        path_sect = path / section if len(methods) > 1 else path

        if method.grid:
            try:
                run_grid(
                    experiments,
                    method.grid,
                    path_sect,
                    plot_level,
                    method.fitmethod,
                    session=session,
                )
            except KeyboardInterrupt:
                print_calculation_stopped_error()
        else:
            _fit_groups(
                experiments,
                path_sect,
                plot_level,
                method.fitmethod,
                method.statistics,
                session=session,
            )
