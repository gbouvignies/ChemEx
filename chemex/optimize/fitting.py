"""The fitting module contains the code for fitting the experimental data."""
import sys
from typing import Any
from typing import Union

import lmfit as lm
import tqdm

import chemex.optimize.methods
from chemex import helper as ch
from chemex.nmr.spin_system import SpinSystem
from chemex.optimize import gridding as cog
from chemex.optimize import grouping
from chemex.optimize import helper as coh
from chemex.parameters import helper as cph
from chemex.parameters import settings as cps


class Fit:
    def __init__(self, experiments, path, plot, defaults):
        self._experiments = experiments
        self._path = path
        self._method = {"": {}}
        self._plot = plot
        self._defaults = defaults

    def read_methods(self, filenames=None):
        if filenames is None:
            return
        self._method = chemex.optimize.methods.read_methods(filenames)

    def run_methods(self, params, path=None, plot=None, set_default=True):

        params_ = params.copy()

        if path is None:
            path = self._path

        if plot is None:
            plot = self._plot

        for index, (section, settings) in enumerate(self._method.items()):

            if section:
                ch.header2(f"{section.upper()}")

            # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
            selection = _get_selection(settings)
            self._experiments.select(selection)

            if not self._experiments:
                print("No data to fit...")
                continue

            # Set the fitting algorithm
            fitmethod = settings.get("fitmethod", "least_squares")
            print(f'\nFitting method -> "{fitmethod}"')

            # Update the parameter "vary" and "expr" status
            cps.set_status(params_, settings)
            if index == 0 and set_default:
                cps.set_values(params_, self._defaults)

            # Read grid search parameters
            grid = settings.get("grid")

            # Read statistics parameters
            statistics = settings.get("statistics")

            if grid and statistics:
                print(
                    'Warning: "GRID" and "STATISTICS" options are mutually '
                    'exclusive. Only the "GRID" calculation will be run.'
                )
                statistics = None

            path_sect = path / section.upper() if len(self._method) > 1 else path

            if grid:
                g_params = self._run_grid(grid, params_, path_sect, plot, fitmethod)
            else:
                g_params = self._fit_groups(
                    params_, path_sect, plot, fitmethod, statistics
                )

            params_.update(g_params)

        return params_

    def _fit_groups(self, params, path, plot, fitmethod, statistics=None):

        groups = grouping.create_groups(self._experiments, params)

        plot_flg = (plot == "normal" and len(groups) == 1) or plot == "all"

        params_list = []

        print("Minimizing...\n")

        for group in groups:

            g_exp = group["experiments"]
            g_par = group["params"]
            g_pat = path / group["path"]
            g_mes = group["message"]

            if g_mes:
                ch.header3(g_mes)

            best_par = minimize(g_exp, g_par, fitmethod, verbose=True)
            coh.post_fit(g_exp, best_par, g_pat, plot_flg)
            params_list.append(best_par)

            # Run Monte Carlo and/or bootstrap analysis
            self.run_statistics(best_par, g_pat, fitmethod, statistics)

        params_merged = cph.merge(params_list)

        if len(groups) > 1:
            ch.header3("All groups")
            coh.post_fit(
                self._experiments, params_merged, path / "All", plot != "nothing"
            )

        return params_merged

    def _run_grid(self, grid, params, path, plot, fitmethod):

        print("Running the grid search...\n")

        grid, params = cps.read_grid(grid, params)

        experiments = self._experiments
        groups = grouping.create_groups(experiments, params)

        grid_results = []
        for group in groups:
            print(f"{group['message']}")
            grid_result = cog.run_grid(group, grid, path, fitmethod)
            grid_results.append(grid_result)

        grids_combined = cog.combine_grids(grid, grid_results)
        grids_1d = cog.make_grids_nd(grid, grids_combined, params, 1)
        grids_2d = cog.make_grids_nd(grid, grids_combined, params, 2)
        cog.set_params_from_grid(grids_1d, params)
        cog.plot_grid_1d(grids_1d, params, path / "Grid")
        cog.plot_grid_2d(grids_2d, params, path / "Grid")

        if len(groups) > 1:
            ch.header3("All groups")
            coh.post_fit(self._experiments, params, path / "All", plot != "nothing")

        return params

    def run_statistics(self, params, path, fitmethod, statistics):

        if statistics is None:
            return

        methods = {
            "mc": {
                "message": "Monte Carlo",
                "filename": "monte_carlo.out",
                "args": [params],
            },
            "bs": {
                "message": "bootstrap",
                "filename": "bootstrap.out",
                "args": [],
            },
            "bsn": {
                "message": "nucleus-based bootstrap",
                "filename": "bootstrap_ns.out",
                "args": [],
            },
        }

        fnames_vary = [param.name for param in params.values() if param.vary]

        for method_name, iter_nb in statistics.items():

            method_name = method_name.lower()

            method = methods[method_name]

            ch.header3(f"Running {method['message']} simulations...")

            with open(path / method["filename"], "w") as fileout:

                fileout.write(coh.print_header(params, fnames_vary))

                for _ in tqdm.tqdm(range(iter_nb)):
                    exp_mc = self._experiments.statistics[method_name](*method["args"])
                    params_mc = params.copy()
                    params_mc = exp_mc.select_params(params)
                    params_mc = minimize(exp_mc, params_mc, fitmethod, verbose=False)
                    chisqr = coh.calculate_statistics(exp_mc, params_mc).get("chisqr")
                    fileout.write(coh.print_values_stat(params_mc, fnames_vary, chisqr))

            print()


def minimize(experiments, params, fitmethod=None, verbose=True):

    if fitmethod is None:
        fitmethod = "least_squares"

    kws = {
        "leastsq": {"factor": 0.1},
        "brute": {"keep": "all"},
        "least_squares": {"verbose": 2 if verbose else 0},
        "basinhopping": {"disp": verbose},
    }

    if fitmethod == "leastsq":
        experiments.verbose = verbose

    minimizer = lm.Minimizer(experiments.residuals, params)

    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
        result = minimizer.result
    except ValueError:
        result = minimizer.result
        sys.stderr.write("\n -- Got a ValueError: minimization stopped --\n")

    if verbose:
        print()

    experiments.verbose = False

    return result.params


def _get_selection(
    settings: dict[str, Any]
) -> dict[str, Union[list[SpinSystem], str, None]]:
    selection = {}
    for option in ("include", "exclude"):
        values = settings.get(option)
        if values is None or isinstance(values, str):
            selection[option] = values
            continue
        selection[option] = tuple(SpinSystem(value) for value in values)
    return selection
