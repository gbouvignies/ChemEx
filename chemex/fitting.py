"""The fitting module contains the code for fitting the experimental data."""
import copy
import sys

import lmfit as lm
import numpy as np
import scipy.stats as ss

import chemex.helper as ch
import chemex.nmr.helper as cnh
import chemex.parameters.helper as cph
import chemex.parameters.name as cpn
import chemex.parameters.settings as cps


def _pop_fitmethod(settings):
    fitmethod = settings.pop("fitmethod", "leastsq")
    print(f'\nFitting method -> "{fitmethod}"')
    return fitmethod


def _select(experiments, settings):
    selection = settings.pop("selection", None)
    if selection is None:
        return
    if isinstance(selection, list) and not isinstance(selection, str):
        selection = [cnh.SpinSystem(name) for name in selection]
    print("\nSelecting profiles...")
    experiments.select(selection)
    print(f"  - Profile(s): {len(experiments)}")


class Fit:
    def __init__(self, experiments, path, plot):
        self._experiments = experiments
        self._path = path
        self._method = {"Standard Calculation": {}}
        self._method_file = None
        self._result = None
        self._plot = plot

    def read_method(self, filename):
        if filename is None:
            self._method = {"": {}}
            self._method_file = None
        else:
            self._method = ch.read_toml(filename)
            self._method_file = filename

    def fit(self, params, method=None, path=None):
        if method is None:
            method = self._method
        if path is None:
            path = self._path
        params_ = copy.deepcopy(params)
        for section, settings in method.items():
            if section:
                ch.header2(f"\n{section.upper()}")
            fitmethod = _pop_fitmethod(settings)
            _select(self._experiments, settings)
            cps.set_param_status(params_, settings)
            if not self._experiments:
                print("No data to fit...")
                continue
            groups = self._cluster_data(params_)
            multi_groups = len(groups) > 1
            plot_group_flg = self._plot == "all" or (
                not multi_groups and self._plot == "normal"
            )
            section_path = section.upper() if len(method) > 1 else ""
            g_params_all = []
            for index, (g_name, (g_experiments, g_params)) in enumerate(groups.items()):
                group_path = path / section_path
                if multi_groups:
                    group_path = group_path / "Clusters" / str(g_name)
                    print(f"\n\n-- Cluster {index + 1}/{len(groups)} ({g_name}) --")
                g_params = _minimize(g_experiments, g_params, fitmethod)
                _write_files(g_experiments, g_params, group_path)
                if plot_group_flg:
                    _write_plots(g_experiments, g_params, group_path)
                g_params_all.append(g_params)
            g_params_merged = cph.merge(g_params_all)
            if multi_groups:
                print("\n\n-- All clusters --")
                _print_chisqr(self._experiments, g_params_merged)
                path_ = path / section_path / "All"
                _write_files(self._experiments, g_params_merged, path_)
                if self._plot != "nothing":
                    _write_plots(self._experiments, g_params_merged, path_)
            params_.update(g_params_merged)
        return params_

    def bootstrap(self, params, iter_nb):
        if iter_nb is None:
            return
        ch.header1("Running bootstrap simulations")
        ndigits = len(str(iter_nb))
        for index in range(1, iter_nb + 1):
            ch.header2(f"Iteration {index} out of {iter_nb}")
            path = self._path / "Bootstrap" / f"{index:0{ndigits}}"
            experiments = self._experiments.bootstrap()
            params_ = _minimize(experiments, params)
            _write_files(experiments, params_, path)

    def monte_carlo(self, params, iter_nb):
        if iter_nb is None:
            return
        ch.header1("Running Monte Carlo simulations")
        ndigits = len(str(iter_nb))
        for index in range(1, iter_nb + 1):
            ch.header2(f"Iteration {index} out of {iter_nb}")
            path = self._path / "MonteCarlo" / f"{index:0{ndigits}}"
            experiments = self._experiments.monte_carlo(params)
            params_ = _minimize(experiments, params)
            _write_files(experiments, params_, path)

    def _cluster_data(self, params):
        """Find clusters of datapoints that depend on disjoint sets of variables.

        For example, if the population of the minor state and the exchange
        rate are set to 'fix', chances are that the fit can be decomposed
        residue-specifically.
        """
        params_ = cph.merge(
            [{name: params[name] for name in self._experiments.params_default}]
        )
        par_name_groups = self._group_par_names(params_)
        if len(par_name_groups) <= 1:
            return {cpn.ParamName(): (self._experiments, params_)}
        groups = self._group_data(params_, par_name_groups)
        return {name: groups[name] for name in sorted(groups)}

    def _group_par_names(self, params):
        par_names_vary = {
            name for name, param in params.items() if param.vary and not param.expr
        }
        clusters = []
        for par_names in self._experiments.par_name_sets:
            varies = par_names & par_names_vary
            found = False
            for par_name_cluster in clusters:
                if varies & par_name_cluster:
                    par_name_cluster |= varies
                    found = True
                    break
            if not found and varies:
                clusters.append(varies)
        return clusters

    def _group_data(self, params, par_name_groups):
        clusters_ = {}
        for par_names in par_name_groups:
            cluster_experiments = self._experiments.get_relevant_subset(par_names)
            params_c = {
                name: params[name] for name in cluster_experiments.params_default
            }
            cluster_params = cph.merge([params_c])
            cluster_name = cluster_experiments.get_cluster_name()
            clusters_[cluster_name] = (cluster_experiments, cluster_params)
        return clusters_


def _minimize(experiments, params, fitmethod=None):
    if fitmethod is None:
        fitmethod = "least_squares"
    kws = {}
    if fitmethod == "brute":
        kws["keep"] = "all"
    elif fitmethod == "least_squares":
        kws["verbose"] = 2
    elif fitmethod in ["basinhopping", "differential_evolution"]:
        kws["disp"] = True
    if fitmethod in ["leastsq", "shgo", "dual_annealing"]:
        experiments.verbose = True

    print("\nMinimizing...")
    minimizer = lm.Minimizer(experiments.residuals, params)
    try:
        result = minimizer.minimize(method=fitmethod, **kws)
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
        result = minimizer.result
    except ValueError:
        result = minimizer.result
        sys.exit((result.params.pretty_print()))
    _print_chisqr(experiments, result.params)
    return result.params


def _write_files(experiments, params, path):
    """Write the results of the fit to output files.

    The files below are created and contain the following information:
      - parameters.toml: fitting parameters and their uncertainties
      - contstraints.fit: expression used for constraining parameters
      - *.dat: experimental and fitted data
      - statistics.fit: statistics for the fit

    """
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


def _write_plots(experiments, params, path):
    """Plot the experimental and fitted data."""
    print("\nPlotting data...")
    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_, params=params)
    except KeyboardInterrupt:
        print("  - Plotting cancelled")
