"""The fitting module contains the code for fitting the experimental data."""
import copy
import sys

import lmfit as lm
import numpy as np
import scipy.stats as ss

import chemex.containers.plot as ccp
import chemex.helper as ch
import chemex.parameters.helper as cph
import chemex.parameters.settings as cps


class Fit:
    def __init__(self, experiments, path, plot):
        self._experiments = experiments
        self._path = path
        self._method = {"": {}}
        self._result = None
        self._plot = plot

    def read_method(self, filename):
        if filename is None:
            return
        self._method = ch.read_toml(filename)

    def fit(self, params, path=None):

        if path is None:
            path = self._path

        params_ = copy.deepcopy(params)

        for index, (section, settings) in enumerate(self._method.items()):

            if section:
                ch.header2(f"\n{section.upper()}")

            # Set the fitting algorithm
            fitmethod = _pop_fitmethod(settings)

            # Select a subset of profiles based on "INCLUDE" and "EXCLUDE"
            selection = {key: settings.pop(key, None) for key in ("include", "exclude")}
            self._experiments.select(selection)

            # Update the parameter "vary" and "expr" status
            cps.set_status(params_, settings)

            # Hack to set back parameters without an "expr" to the user values
            if index == 0:
                cps.put_back_starting_values(params_)

            if not self._experiments:
                print("No data to fit...")
                continue

            # Make cluster of data depending on indendent set of parameters
            groups = self._cluster_data(params_)

            # Set section flags and path
            multi_groups = len(groups) > 1
            plot_group_flg = self._plot == "all" or (
                not multi_groups and self._plot == "normal"
            )
            section_path = section.upper() if len(self._method) > 1 else ""

            g_params_all = []
            for index, (g_name, (g_experiments, g_params)) in enumerate(groups.items()):
                group_path = path / section_path
                if multi_groups:
                    group_path = group_path / "Clusters" / str(g_name)
                    print(f"\n\n-- Cluster {index + 1}/{len(groups)} ({g_name}) --")
                g_params = _minimize(g_experiments, g_params, fitmethod)
                _write_files(g_experiments, g_params, group_path)
                if plot_group_flg:
                    ccp.write_plots(g_experiments, g_params, group_path)
                g_params_all.append(g_params)

            g_params_merged = cph.merge(g_params_all)

            if multi_groups:
                print("\n\n-- All clusters --")
                _print_chisqr(self._experiments, g_params_merged)
                path_ = path / section_path / "All"
                _write_files(self._experiments, g_params_merged, path_)
                if self._plot != "nothing":
                    ccp.write_plots(self._experiments, g_params_merged, path_)

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
        pnames_groups = self._group_pnames(params_)
        groups = self._group_data(params_, pnames_groups)
        return {name: groups[name] for name in sorted(groups)}

    def _group_pnames(self, params):
        pnames_vary = {
            name for name, param in params.items() if param.vary and not param.expr
        }
        clusters = []
        for pnames in self._experiments.pname_sets:
            varies = pnames & pnames_vary
            found = False
            for pname_cluster in clusters:
                if varies & pname_cluster:
                    pname_cluster |= varies
                    found = True
                    break
            if not found and varies:
                clusters.append(varies)
        return clusters

    def _group_data(self, params, par_name_groups):
        clusters_ = {}
        for pnames in par_name_groups:
            cluster_experiments = self._experiments.get_relevant_subset(pnames)
            params_c = {
                name: params[name] for name in cluster_experiments.params_default
            }
            cluster_params = cph.merge([params_c])
            cluster_name = cluster_experiments.get_cluster_name()
            clusters_[cluster_name] = (cluster_experiments, cluster_params)
        return clusters_


def _pop_fitmethod(settings):
    fitmethod = settings.pop("fitmethod", "leastsq")
    print(f'\nFitting method -> "{fitmethod}"')
    return fitmethod


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
        sys.exit(result.params.pretty_print())
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
