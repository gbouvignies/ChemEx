"""The fitting module contains the code for fitting the experimental data."""
import itertools as it
import pathlib
import sys

import lmfit as lm
import numpy as np
import scipy.stats as ss
import tqdm

import chemex.containers.plot as ccp
import chemex.helper as ch
import chemex.parameters.helper as cph
import chemex.parameters.name as cpn
import chemex.parameters.settings as cps


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

    def run_methods(self, params, path=None, plot=None):

        if path is None:
            path = self._path

        if plot is None:
            plot = self._plot

        params_ = params.copy()

        for index, (section, settings) in enumerate(self._method.items()):

            if section:
                ch.header2(f"\n{section.upper()}")

            # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
            selection = {key: settings.get(key) for key in ("include", "exclude")}
            self._experiments.select(selection)
            if not self._experiments:
                print("No data to fit...")
                continue

            # Set the fitting algorithm
            fitmethod = settings.get("fitmethod", "leastsq")
            print(f'\nFitting method -> "{fitmethod}"')

            # Update the parameter "vary" and "expr" status
            cps.set_status(params_, settings)
            if index == 0:
                cps.set_values(params_, self._defaults)

            # Red grid search parameters
            grid = settings.get("grid")

            path_sect = path / section.upper() if len(self._method) > 1 else path

            if grid:
                self._experiments.verbose = False
                g_params = self._run_grid(params_, grid, path_sect, plot, fitmethod)

            else:
                # Make cluster of data depending on indendent set of parameters
                self._experiments.verbose = True
                g_params = self._fit_groups(params_, path_sect, plot, fitmethod)

            params_.update(g_params)

        return params_

    def _run_grid(self, params, grid, path, plot, fitmethod):

        path_grid = path / "Grid"
        path_grid.mkdir(parents=True, exist_ok=True)

        params_select = self._experiments.select_params(params)

        grid_values, params_grid = cps.read_grid(grid, params_select)

        groups = self._create_groups(params_grid)

        print("\nRunning the grid search...\n")

        best_chi2 = 1e32
        best_params = None
        grid_n = np.prod([len(values) for values in grid_values.values()])

        with (path_grid / "grid.out").open("w") as fileout:

            names = [
                "# ",
                " ".join(
                    f"{str(params[fname].user_data['pname']):^17s}"
                    for fname in grid_values
                ),
                f"{'redchi2':^17s}",
                "\n",
            ]
            fileout.write("".join(names))

            redchi2 = []

            for values in tqdm.tqdm(it.product(*grid_values.values()), total=grid_n):

                for fname, value in zip(grid_values, values):
                    # Set parameters using grid values
                    for group in groups:
                        if fname in group["params"]:
                            group["params"][fname].value = value

                params = cph.merge(
                    _minimize(group["experiments"], group["params"], fitmethod)
                    for group in groups
                )

                statistics = _calculate_statistics(self._experiments, params)

                redchi2.append(statistics["redchi"])

                output = [
                    "  ",
                    " ".join(f"{value:17.8e}" for value in values),
                    f"{redchi2[-1]:17.8e}",
                    "\n",
                ]
                fileout.write("".join(output))
                fileout.flush()

                if redchi2[-1] < best_chi2:
                    best_chi2 = redchi2[-1]
                    best_params = params.copy()

        _post_fit(self._experiments, best_params, path, False)

        if plot != "nothing":
            _plot_grid(grid_values, redchi2, params, path_grid)

        return best_params

    def _fit_groups(self, params, path, plot, fitmethod):

        groups = self._create_groups(params)

        multi_groups = len(groups) > 1
        plot_group_flg = plot == "all" or (not multi_groups and plot == "normal")

        params_list = []
        for group in self._create_groups(params):
            g_mes, g_exp, g_par, g_pat = (
                group[key] for key in ("message", "experiments", "params", "path")
            )
            print(f"{g_mes}\nMinimizing...\n")
            g_par = _minimize(g_exp, g_par, fitmethod)
            _post_fit(g_exp, g_par, path / g_pat, plot_group_flg)
            params_list.append(g_par)
        params = cph.merge(params_list)

        if multi_groups:
            print("\n\n-- All clusters --")
            self._experiments.verbose = False
            _post_fit(self._experiments, params, path / "All", plot != "nothing")

        return params

    def mc_simulations(self, params, iter_nb, name="mc"):
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
            ch.header2(f"Iteration {index} out of {iter_nb}")
            path = self._path / method["folder"] / f"{index:0{ndigits}}"
            self._experiments = getattr(experiments, method["attr"])(*method["args"])
            self.run_methods(params_mc, path, plot="nothing")
        self._experiments = experiments

    def _create_groups(self, params):
        """Create clusters of datapoints that depend on disjoint sets of variables.

        For example, if the population of the minor state and the exchange
        rate are set to 'fix', chances are that the fit can be decomposed
        residue-specifically.
        """
        params_exp = self._experiments.select_params(params)
        pname_groups = self._group_pnames(params_exp)
        number = len(pname_groups)
        if number > 1:
            path = pathlib.Path("Clusters")
            message = f"\n\n-- Cluster {{index}}/{number} ({{name}}) --\n"
        else:
            path = pathlib.Path("")
            message = ""
        groups_dict = {}
        for index, pnames in enumerate(pname_groups, start=1):
            experiments = self._experiments.get_relevant_subset(pnames)
            name = experiments.get_cluster_name() if number > 1 else cpn.ParamName()
            group = {
                "experiments": experiments,
                "params": experiments.select_params(params_exp),
                "path": path / name.folder,
            }
            groups_dict[name] = group
        groups = []
        for index, (name, group) in enumerate(sorted(groups_dict.items()), start=1):
            group["message"] = message.format(index=index, name=name)
            groups.append(group)
        return groups

    def _group_pnames(self, params):
        pnames_vary = {
            name for name, param in params.items() if param.vary and not param.expr
        }
        groups = []
        for pnames in self._experiments.pname_sets:
            varies = pnames & pnames_vary
            found = False
            for group in groups:
                if varies & group:
                    group |= varies
                    found = True
                    break
            if not found and varies:
                groups.append(varies)
        return groups


def _minimize(experiments, params, fitmethod=None):
    if fitmethod is None:
        fitmethod = "leastsq"
    kws = {}
    if fitmethod == "leastsq":
        kws["factor"] = 0.1
    elif fitmethod == "brute":
        kws["keep"] = "all"
    elif fitmethod == "least_squares":
        kws["verbose"] = 2
        experiments.verbose = False
    elif fitmethod in ["basinhopping", "differential_evolution"]:
        kws["disp"] = True
        experiments.verbose = False

    minimizer = lm.Minimizer(experiments.residuals, params)
    try:
        result = minimizer.minimize(method=fitmethod, **kws)
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
        result = minimizer.result
    except ValueError:
        result = minimizer.result
        sys.stderr.write("\n -- Got a ValueError: minimization stopped --\n")
    return result.params


def _post_fit(experiments, params, path, plot=False):
    _print_chisqr(experiments, params)
    _write_files(experiments, params, path)
    if plot:
        ccp.write_plots(experiments, params, path)


def _write_files(experiments, params, path):
    """Write the results of the fit to output files."""
    print(f'\nWriting results -> "{path}/"')
    path.mkdir(parents=True, exist_ok=True)
    cps.write_par(params, path)
    experiments.write(params, path)
    _write_statistics(experiments, params, path=path)


def _write_statistics(experiments, params, path):
    """Write fitting statistics to a file."""
    statistics = _calculate_statistics(experiments, params)
    filename = path / "statistics.toml"
    with open(filename, "w") as f:
        f.write(f"number of data points          = {statistics['ndata']}\n")
        f.write(f"number of variables            = {statistics['nvarys']}\n")
        f.write(f"chi-square                     = {statistics['chisqr']: .5e}\n")
        f.write(f"reduced-chi-square             = {statistics['redchi']: .5e}\n")
        f.write(f"chi-squared test               = {statistics['pvalue']: .5e}\n")
        f.write(f"Kolmogorov-Smirnov test        = {statistics['ks_pvalue']: .5e}\n")
        f.write(f"Akaike Information Criterion   = {statistics['aic']: .5e}\n")
        f.write(f"Bayesian Information Criterion = {statistics['bic']: .5e}\n")


def _print_chisqr(experiments, params):
    statistics = _calculate_statistics(experiments, params)
    print("")
    print(f"Final Chi2        : {statistics['chisqr']:.3e}")
    print(f"Final Reduced Chi2: {statistics['redchi']:.3e}")


def _calculate_statistics(experiments, params):
    residuals = experiments.residuals(params)
    ndata = len(residuals)
    nvarys = len([param for param in params.values() if param.vary and not param.expr])
    chisqr = sum(residuals ** 2)
    redchi = chisqr / max(1, ndata - nvarys)
    _neg2_log_likel = ndata * np.log(chisqr / ndata)
    aic = _neg2_log_likel + 2 * nvarys
    bic = _neg2_log_likel + np.log(ndata) * nvarys
    _, ks_p_value = ss.kstest(residuals, "norm")
    pvalue = 1.0 - ss.chi2.cdf(chisqr, ndata - nvarys)
    return {
        "ndata": ndata,
        "nvarys": nvarys,
        "chisqr": chisqr,
        "redchi": redchi,
        "pvalue": pvalue,
        "ks_pvalue": ks_p_value,
        "aic": aic,
        "bic": bic,
    }


def _plot_grid(grid_values, redchi2, params, path):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square value per parameter and contour
    plots for all combination of two parameters.
    """

    shape = [len(values) for values in grid_values.values()]
    redchi2 = np.array(redchi2).reshape(shape)

    _plot_grid_1d(grid_values, redchi2, params, path)

    if len(grid_values) > 1:
        _plot_grid_2d(grid_values, redchi2, params, path)


def _plot_grid_1d(grid_values, redchi2, params, path):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path / "grid_1d.pdf") as pdf:
        for i, (fname, values) in enumerate(grid_values.items()):
            red_axis = tuple(a for a in range(len(grid_values)) if a != i)
            chi2_vector = np.minimum.reduce(redchi2, axis=red_axis)
            _fig, ax = plt.subplots(1, 1)
            ax.plot(values, chi2_vector, "o", ms=3)
            min_index = np.argmin(chi2_vector)
            ax.axvline(values[min_index], ls="dashed", color="r")
            ax.set_xlabel(str(params[fname].user_data["pname"]))
            ax.set_ylabel(r"$\chi^{2}$")
            pdf.savefig()


def _plot_grid_2d(grid_values, redchi2, params, path):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    from matplotlib import cm

    with PdfPages(path / "grid_2d.pdf") as pdf:
        combinations = it.combinations(enumerate(grid_values.items()), 2)
        for (i, (n1, v1)), (j, (n2, v2)) in combinations:
            fig, ax = plt.subplots(1, 1)
            grid_x, grid_y = np.meshgrid(v1, v2)
            red_axis = tuple(a for a in range(len(grid_values)) if a not in (i, j))
            chi2_matrix = np.minimum.reduce(redchi2, axis=red_axis).T
            cs = ax.pcolormesh(
                grid_x,
                grid_y,
                chi2_matrix,
                norm=LogNorm(),
                cmap=cm.Blues_r,
                shading="gouraud",
            )
            fig.colorbar(cs)
            min_indexes = np.unravel_index(np.argmin(chi2_matrix), chi2_matrix.shape)
            min_x = grid_x[min_indexes]
            min_y = grid_y[min_indexes]
            ax.axvline(min_x, ls="dashed", color="r")
            ax.axhline(min_y, ls="dashed", color="r")
            ax.plot(min_x, min_y, "rs", ms=3)
            ax.set_xlabel(str(params[n1].user_data["pname"]))
            ax.set_ylabel(str(params[n2].user_data["pname"]))
            pdf.savefig()
