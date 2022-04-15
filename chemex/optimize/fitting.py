"""The fitting module contains the code for fitting the experimental data."""
from __future__ import annotations

from pathlib import Path

from rich.progress import track

from chemex.configuration.methods import Methods
from chemex.configuration.methods import Statistics
from chemex.containers.experiments import Experiments
from chemex.containers.experiments import generate_exp_for_statistics
from chemex.messages import print_fitmethod
from chemex.messages import print_group_name
from chemex.messages import print_minimizing
from chemex.messages import print_no_data
from chemex.messages import print_running_statistics
from chemex.messages import print_step_name
from chemex.optimize.gridding import run_grid
from chemex.optimize.grouping import create_groups
from chemex.optimize.helper import calculate_statistics
from chemex.optimize.helper import execute_post_fit
from chemex.optimize.helper import execute_post_fit_groups
from chemex.optimize.helper import print_header
from chemex.optimize.helper import print_values_stat
from chemex.optimize.minimizer import minimize
from chemex.parameters import database


def _run_statistics(
    experiments: Experiments,
    path: Path,
    fitmethod: str | None = None,
    statistics: Statistics | None = None,
):
    if statistics is None:
        return

    methods = {
        "mc": {"message": "Monte Carlo", "filename": "monte_carlo.out"},
        "bs": {"message": "bootstrap", "filename": "bootstrap.out"},
        "bsn": {"message": "nucleus-based bootstrap", "filename": "bootstrap_ns.out"},
    }

    params_lf = database.build_lmfit_params(experiments.param_ids)
    ids_vary = [param.name for param in params_lf.values() if param.vary]

    for statistic_name, iter_nb in statistics.dict().items():

        if iter_nb is None:
            continue

        method = methods[statistic_name]

        print_running_statistics(method["message"])

        with open(path / method["filename"], "w") as fileout:

            fileout.write(print_header(ids_vary))

            for _ in track(range(iter_nb), total=iter_nb, description="   "):
                exp_stat = generate_exp_for_statistics(experiments, statistic_name)
                params_lf = database.build_lmfit_params(exp_stat.param_ids)
                params_fit = minimize(exp_stat, params_lf, fitmethod, verbose=False)
                chisqr = calculate_statistics(exp_stat, params_fit).get("chisqr", 1e32)
                fileout.write(print_values_stat(params_fit, ids_vary, chisqr))


def _fit_groups(
    experiments: Experiments,
    path: Path,
    plot: str,
    fitmethod: str,
    statistics: Statistics | None,
) -> None:
    groups = create_groups(experiments)

    plot_flg = (plot == "normal" and len(groups) == 1) or plot == "all"

    print_minimizing()

    for group in groups:

        group_lmfit_params = database.build_lmfit_params(group.experiments.param_ids)
        group_path = path / group.path

        if message := group.message:
            print_group_name(message)

        best_lmfit_params = minimize(
            group.experiments, group_lmfit_params, fitmethod, verbose=True
        )

        database.update_from_parameters(best_lmfit_params)
        execute_post_fit(group.experiments, group_path, plot_flg)

        # Run Monte Carlo and/or bootstrap analysis
        _run_statistics(
            group.experiments,
            group_path,
            fitmethod,
            statistics,
        )

    if len(groups) > 1:
        execute_post_fit_groups(experiments, path, plot)


def run_methods(
    experiments: Experiments, methods: Methods, path: Path, plot_level: str
) -> None:
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
        database.set_parameter_status(method)

        if method.grid and method.statistics:
            print(
                'Warning: "GRID" and "STATISTICS" options are mutually '
                'exclusive. Only the "GRID" calculation will be run.'
            )
            method.statistics = None

        path_sect = path / section if len(methods) > 1 else path

        if method.grid:
            run_grid(experiments, method.grid, path_sect, plot_level, method.fitmethod)
        else:
            _fit_groups(
                experiments, path_sect, plot_level, method.fitmethod, method.statistics
            )
