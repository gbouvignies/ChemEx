"""The fitting module contains the code for fitting the experimental data."""
import sys

import lmfit as lm

from chemex import helper as ch
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
        self._method = cps.read_methods(filenames)

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
            selection = {key: settings.get(key) for key in ("include", "exclude")}
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

            path_sect = path / section.upper() if len(self._method) > 1 else path

            if grid:
                g_params = self._run_grid(grid, params_, path_sect, plot, fitmethod)
            else:
                g_params = self._fit_groups(params_, path_sect, plot, fitmethod)

            params_.update(g_params)

        return params_

    def _fit_groups(self, params, path, plot, fitmethod):

        groups = grouping.create_groups(self._experiments, params)

        plot_flg = (plot == "normal" and len(groups) == 1) or plot == "all"

        params_list = []

        print("Minimizing...\n")

        for group in groups:

            g_exp = group["experiments"]
            g_par = group["params"]
            g_pat = group["path"]
            g_mes = group["message"]

            if g_mes:
                print(f"{g_mes}")

            best_par = minimize(g_exp, g_par, fitmethod, verbose=True)
            coh.post_fit(g_exp, best_par, path / g_pat, plot_flg)
            params_list.append(best_par)

        params = cph.merge(params_list)

        if len(groups) > 1:
            print("\n\n-- All groups --")
            coh.post_fit(self._experiments, params, path / "All", plot != "nothing")

        return params

    def _run_grid(self, grid, params, path, plot, fitmethod):

        print("\nRunning the grid search...")

        grid, params = cps.read_grid(grid, params)

        experiments = self._experiments
        groups = grouping.create_groups(experiments, params)

        grid_results = []
        for group in groups:
            g_mes = group["message"]
            print(f"{g_mes}")
            grid_result = cog.run_grid(group, grid, path, fitmethod)
            grid_results.append(grid_result)

        params = cog.combine_grids(experiments, params, grid_results, path / "Grid")

        if len(groups) > 1:
            print("\n\n-- All groups --")
            coh.post_fit(self._experiments, params, path / "All", plot != "nothing")

        return params

    def mc_simulations(self, params, iter_nb, name=None):

        if iter_nb is None:
            return

        methods = {
            "mc": {
                "message": "Monte Carlo",
                "folder": "MonteCarlo",
                "attr": "monte_carlo",
                "args": [params],
            },
            "bs": {
                "message": "bootstrap",
                "folder": "Bootstrap",
                "attr": "bootstrap",
                "args": [],
            },
            "bsn": {
                "message": "nucleus-specific",
                "folder": "BootstrapNS",
                "attr": "bootstrap_ns",
                "args": [],
            },
        }

        method = methods[name]

        ch.header1(f"Running {method['message']} simulations")

        ndigits = len(str(iter_nb))

        experiments = self._experiments

        for index in range(1, iter_nb + 1):
            params_mc = params.copy()
            ch.header2(f"\nIteration {index} out of {iter_nb}")
            path = self._path / method["folder"] / f"{index:0{ndigits}}"
            self._experiments = getattr(experiments, method["attr"])(*method["args"])
            self.run_methods(params_mc, path, plot="nothing", set_default=False)

        self._experiments = experiments


def minimize(experiments, params, fitmethod=None, verbose=True):

    if fitmethod is None:
        fitmethod = "least_squares"

    kws = {
        "leastsq": {"factor": 0.1},
        "brute": {"keep": "all"},
        "least_squares": {"verbose": 2 if verbose else 0},
        "basinhopping": {"disp": verbose},
    }

    minimizer = lm.Minimizer(experiments.residuals, params)

    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
        result = minimizer.result
    except ValueError:
        result = minimizer.result
        sys.stderr.write("\n -- Got a ValueError: minimization stopped --\n")

    return result.params
