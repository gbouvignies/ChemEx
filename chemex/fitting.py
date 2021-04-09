"""The fitting module contains the code for fitting the experimental data."""
import pathlib
import sys

import lmfit as lm
import numpy as np
import scipy.stats as ss

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

    def read_methods(self, filename):
        if filename is None:
            return
        self._method = ch.read_toml(filename)

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
            selection = {key: settings.pop(key, None) for key in ("include", "exclude")}
            self._experiments.select(selection)
            if not self._experiments:
                print("No data to fit...")
                continue

            # Set the fitting algorithm
            fitmethod = settings.pop("fitmethod", "leastsq")
            print(f'\nFitting method -> "{fitmethod}"')

            # Update the parameter "vary" and "expr" status
            cps.set_status(params_, settings)
            if index == 0:
                cps.set_values(params_, self._defaults)

            # Make cluster of data depending on indendent set of parameters
            groups = self._create_groups(params_)
            g_params_merged = self._fit_groups(path, plot, section, fitmethod, groups)
            params_.update(g_params_merged)

        return params_

    def _fit_groups(self, path, plot, section, fitmethod, groups):
        multi_groups = len(groups) > 1
        plot_group_flg = plot == "all" or (not multi_groups and plot == "normal")
        if len(self._method) > 1:
            path /= pathlib.Path(section.upper())
        params_list = []
        for group in groups:
            group_path = path / group["path"]
            print(group["message"], end="")
            params = _minimize(group["experiments"], group["params"], fitmethod)
            _write_files(group["experiments"], params, group_path)
            if plot_group_flg:
                ccp.write_plots(group["experiments"], params, group_path)
            params_list.append(params)
        params_merged = cph.merge(params_list)
        if multi_groups:
            print("\n\n-- All clusters --")
            _print_chisqr(self._experiments, params_merged)
            path_ = path / pathlib.Path("All")
            _write_files(self._experiments, params_merged, path_)
            if plot != "nothing":
                ccp.write_plots(self._experiments, params_merged, path_)
        return params_merged

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
    if fitmethod == "brute":
        kws["keep"] = "all"
    elif fitmethod == "least_squares":
        kws["verbose"] = 2
    elif fitmethod == "leastsq":
        kws["factor"] = 0.1
    elif fitmethod in ["basinhopping", "differential_evolution"]:
        kws["disp"] = True
    if fitmethod in ["leastsq", "shgo", "dual_annealing"]:
        experiments.verbose = True

    print("\nMinimizing...\n")
    minimizer = lm.Minimizer(experiments.residuals, params)
    try:
        result = minimizer.minimize(method=fitmethod, **kws)
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
        result = minimizer.result
    except ValueError:
        result = minimizer.result
        sys.stderr.write("\n -- Got a ValueError: minimization stopped --\n")
    _print_chisqr(experiments, result.params)
    return result.params


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
