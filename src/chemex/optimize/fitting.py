"""The fitting module contains the code for fitting the experimental data."""

from pathlib import Path

from rich.progress import track

from chemex.configuration.methods import Methods, Statistics
from chemex.containers.experiments import Experiments, generate_exp_for_statistics
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
    calculate_statistics,
    execute_post_fit,
    execute_post_fit_groups,
    print_header,
    print_values_stat,
)
from chemex.optimize.minimizer import (
    minimize,
    minimize_with_report,
)
from chemex.parameters import database
from chemex.runtime import AnalysisSession
from chemex.runtime.session import ParameterStore


def _get_parameter_store(session: AnalysisSession | None) -> ParameterStore:
    return session.parameters if session is not None else database


def _run_statistics(
    experiments: Experiments,
    path: Path,
    fitmethod: str,
    statistics: Statistics | None = None,
    *,
    session: AnalysisSession | None = None,
) -> None:
    if statistics is None:
        return

    methods = {
        "mc": {"message": "Monte Carlo", "filename": "monte_carlo.out"},
        "bs": {"message": "bootstrap", "filename": "bootstrap.out"},
        "bsn": {"message": "nucleus-based bootstrap", "filename": "bootstrap_ns.out"},
    }

    parameter_store = _get_parameter_store(session)
    params_lf = parameter_store.build_lmfit_params(experiments.param_ids)
    ids_vary = [param.name for param in params_lf.values() if param.vary]

    for statistic_name, iter_nb in statistics.model_dump().items():
        if iter_nb is None:
            continue

        method = methods[statistic_name]

        print_running_statistics(method["message"])

        with (path / method["filename"]).open(mode="w", encoding="utf-8") as fileout:
            fileout.write(print_header(ids_vary))

            try:
                for _ in track(range(iter_nb), total=iter_nb, description="   "):
                    exp_stat = generate_exp_for_statistics(experiments, statistic_name)
                    params_lf = parameter_store.build_lmfit_params(exp_stat.param_ids)
                    params_fit = minimize(exp_stat, params_lf, fitmethod)
                    stats = calculate_statistics(exp_stat, params_fit)
                    chisqr = stats.get("chisqr", 1e32)
                    fileout.write(
                        print_values_stat(params_fit, ids_vary, chisqr),
                    )
            except KeyboardInterrupt:
                print_calculation_stopped_error()
            except ValueError:
                print_value_error()
            finally:
                fileout.flush()


def _fit_groups(  # noqa: PLR0913
    experiments: Experiments,
    path: Path,
    plot: str,
    fitmethod: str,
    statistics: Statistics | None,
    *,
    session: AnalysisSession | None = None,
) -> None:
    parameter_store = _get_parameter_store(session)
    groups = create_groups(experiments, session=session)

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
            session=session,
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
    parameter_store = _get_parameter_store(session)

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
