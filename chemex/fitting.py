"""The fitting module contains the code for fitting the experimental data."""
import copy
import shutil
import sys

import lmfit
from scipy import stats

import chemex.parameters.helper as cph
import chemex.parameters.name as cpn
from chemex import helper as ch
from chemex.parameters import settings as cpp


class Fit:
    def __init__(self, experiments, path, no_plot=False):
        self._experiments = experiments
        self._path = path
        self._no_plot = no_plot
        self._method = {"Standard Calculation": {}}
        self._method_file = None
        self._result = None

    def read_method(self, filename):
        self._method = ch.read_toml(filename)
        self._method_file = filename

    def fit(self, params, method=None, path=None, noplot=False):
        params_ = copy.deepcopy(params)
        fitmethod = None
        if method is None:
            method = self._method
        if path is None:
            path = self._path
        for name, settings in method.items():
            ch.header2(name.upper())
            cpp.set_param_status(params, settings)
            fitmethod = settings.get("fitmethod", "leastsq")
            print(f"\nFitting method: {fitmethod}")
            groups = self._cluster_data(params)
            for group_name, (experiments, group_params) in groups.items():
                group_path = path / name.upper()
                if group_name:
                    print(f"\n[{group_name}]")
                    folder = group_name.to_folder_name()
                    group_path = group_path / "Clusters" / folder
                group_result = _minimize(experiments, group_params, fitmethod)
                self._write_files(experiments, group_result, group_path, noplot)
                params_.update(group_result.params)
            if len(groups) > 1:
                print("[ALL]")
                path_ = path / name / "All"
                result = self._gather_result(params_, fitmethod)
                self._write_files(self._experiments, result, path_, noplot)
        result = self._gather_result(params_, fitmethod)
        return result

    def _write_files(self, experiments, result, path, no_plot=False):
        print("")
        print(f"Final Chi2        : {result.chisqr:.3e}")
        print(f"Final Reduced Chi2: {result.redchi:.3e}")
        path.mkdir(parents=True, exist_ok=True)
        _write(experiments, result, path, self._method_file)
        if not no_plot:
            _plot(experiments, result.params, path)

    def bootstrap(self, params, iter_nb):
        if iter_nb is None:
            return
        ch.header1("Running bootstrap simulations")
        method = {"": {}}
        experiments_orig = self._experiments
        for index in range(1, iter_nb + 1):
            ch.header2(f"Iteration {index} out of {iter_nb}")
            path = self._path / "Bootstrap" / str(index).zfill(len(str(iter_nb)))
            self._experiments = experiments_orig.bootstrap()
            self.fit(params, method, path, noplot=True)

    def monte_carlo(self, params, iter_nb):
        if iter_nb is None:
            return
        ch.header1("Running Monte Carlo simulations")
        method = {"": {}}
        experiments_orig = self._experiments
        for index in range(1, iter_nb + 1):
            ch.header2(f"Iteration {index} out of {iter_nb}")
            path = self._path / "Bootstrap" / str(index).zfill(len(str(iter_nb)))
            self._experiments = experiments_orig.monte_carlo(params)
            self.fit(params, method, path, noplot=True)

    def _cluster_data(self, params):
        """Find clusters of datapoints that depend on disjoint sets of variables.

        For example, if the population of the minor state and the exchange
        rate are set to 'fix', chances are that the fit can be decomposed
        residue-specifically.
        """
        par_name_groups = self._group_par_names(params)
        if len(par_name_groups) <= 1:
            return {"": (self._experiments, params)}
        groups = self._group_data(params, par_name_groups)
        return {name: groups[name] for name in sorted(groups)}

    def _group_data(self, params, par_name_groups):
        clusters_ = {}
        for par_names in par_name_groups:
            cluster_name = _get_cluster_name(par_names)
            cluster_experiments = self._experiments.get_relevant_subset(par_names)
            params_c = {
                name: params[name] for name in cluster_experiments.params_default
            }
            cluster_params = cph.merge([params_c])
            clusters_[cluster_name] = (cluster_experiments, cluster_params)
        return clusters_

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

    def _gather_result(self, params, fitmethod=None):
        minimizer = lmfit.Minimizer(self._experiments.residuals, params)
        result = minimizer.prepare_fit()
        result.residual = self._experiments.residuals(params, verbose=False)
        result.params = params
        result._calculate_statistics()
        result.method = fitmethod
        return result


def _minimize(experiments, params, fitmethod):
    kws = {}
    if fitmethod == "brute":
        kws["keep"] = "all"
    print("\nChi2 / Reduced Chi2:")
    minimizer = lmfit.Minimizer(experiments.residuals, params)
    try:
        result = minimizer.minimize(method=fitmethod, **kws)
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped\n")
        result = minimizer.result
    return result


def _get_result(func, params, fitmethod):
    minimizer = lmfit.Minimizer(func, params)
    result = minimizer.prepare_fit()
    result.residual = func(params, verbose=False)
    result.params = params
    result._calculate_statistics()
    result.method = fitmethod
    return result


def _get_cluster_name(par_names):
    par_names_ = par_names.copy()
    name = cpn.ParamName.from_full_name(par_names_.pop())
    while par_names_:
        name = name.intersection(cpn.ParamName.from_full_name(par_names_.pop()))
    return name


def _write_statistics(result, path):
    """Write fitting statistics to a file."""
    _, ks_p_value = stats.kstest(result.residual, "norm")
    chi2_p_value = 1.0 - stats.chi2.cdf(result.chisqr, result.nfree)
    filename = path / "statistics.toml"
    with open(filename, "w") as f:
        print(f"  - {filename}")
        f.write(f"# Number of data points: {result.ndata}\n")
        f.write(f"# Number of variables: {result.nvarys}\n")
        f.write(f"# Fitting method: {result.method}\n")
        f.write(f"# Number of function evaluations: {result.nfev}\n\n")
        f.write(f"chi-square         = {result.chisqr: .5e}\n")
        f.write(f"reduced-chi-square = {result.redchi: .5e}\n")
        f.write(f"chisq_p-value      = {chi2_p_value: .5e} # Chi-squared test\n")
        f.write(
            "ks_p-value         = {: .5e} # Kolmogorov-Smirnov test\n\n".format(
                ks_p_value
            )
        )
        f.write(f"Akaike Information Criterion = {result.aic: .5e}\n")
        f.write(f"Bayesian Information Criterion = {result.bic: .5e}\n")


def _write(experiments, result, path, method_file=None):
    """Write the results of the fit to output files.

    The files below are created and contain the following information:
      - parameters.fit: fitting parameters and their uncertainties
      - contstraints.fit: expression used for constraining parameters
      - *.dat: experimental and fitted data
      - statistics.fit: statistics for the fit

    """
    print("\nWriting Results...")
    if method_file:
        shutil.copyfile(method_file, path / "fitting-method.toml")
    cpp.write_par(result.params, path=path)
    cpp.write_constraints(result.params, path=path)
    experiments.write(path=path, params=result.params)
    _write_statistics(result, path=path)


def _plot(experiments, params, path):
    """Plot the experimental and fitted data."""
    print("\nPlotting data..")
    path_plots = path / "Plots"
    path_plots.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_plots, params=params)
    except KeyboardInterrupt:
        print("  - Plotting cancelled")
