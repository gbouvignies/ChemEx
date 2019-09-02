"""The fitting module contains the code for fitting the experimental data."""
import sys

import lmfit
from scipy import stats

import chemex.helper as ch
import chemex.parameters.helper as cph
import chemex.parameters.name as cpn
import chemex.parameters.settings as cpp


def run_fit(experiments, params, method_filename, cl_fitmethod):
    """Perform the fit."""
    ch.header1("Fit")
    fit_config = ch.read_toml(method_filename)
    if not fit_config:
        fit_config["Standard Calculation"] = {}
    result = None
    for section in fit_config:
        ch.header2(section)
        cpp.set_param_status(params, fit_config[section].items())
        clusters = find_independent_clusters(experiments, params)
        fitmethod = fit_config[section].get("fitmethod", cl_fitmethod)
        print(f"Fitting method: {fitmethod}\n")
        for name_c, (experiments_c, params_c) in clusters.items():
            if len(clusters) > 1:
                print(f"[{name_c}]")
            result_c = minimize(experiments_c.residuals, params_c, fitmethod)
            params.update(result_c.params)
            print("")
        result = get_result(experiments.residuals, params, fitmethod)
        print("")
        print(f"Final Chi2        : {result.chisqr:.3e}")
        print(f"Final Reduced Chi2: {result.redchi:.3e}")

    return result


def minimize(func, params, fitmethod):
    kws = {}
    if fitmethod == "brute":
        kws["keep"] = "all"
    print("Chi2 / Reduced Chi2:")
    minimizer = lmfit.Minimizer(func, params)
    try:
        result = minimizer.minimize(method=fitmethod, **kws)
    except KeyboardInterrupt:
        sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped\n")
        result = minimizer.result
    return result


def get_result(func, params, fitmethod):
    minimizer = lmfit.Minimizer(func, params)
    result = minimizer.prepare_fit()
    result.residual = func(params, verbose=False)
    result.params = params
    result._calculate_statistics()
    result.method = fitmethod
    return result


def find_independent_clusters(experiments, params):
    """Find clusters of datapoints that depend on disjoint sets of variables.

    For example, if the population of the minor state and the exchange
    rate are set to 'fix', chances are that the fit can be decomposed
    residue-specifically.
    """
    par_names_vary = {
        name for name, param in params.items() if param.vary and not param.expr
    }
    clusters = []
    for par_names in experiments.par_name_sets:
        varies = par_names & par_names_vary
        found = False
        for par_name_cluster in clusters:
            if varies & par_name_cluster:
                par_name_cluster |= varies
                found = True
                break
        if not found and varies:
            clusters.append(varies)
    clusters_ = {}
    for par_names in clusters:
        cluster_name = get_cluster_name(par_names)
        cluster_experiments = experiments.get_relevant_subset(par_names)
        params_c = {name: params[name] for name in cluster_experiments.params_default}
        cluster_params = cph.merge([params_c])
        clusters_[cluster_name] = (cluster_experiments, cluster_params)
    clusters_sorted = {key: clusters_[key] for key in sorted(clusters_)}
    return clusters_sorted


def get_cluster_name(par_names):
    par_names_ = par_names.copy()
    name = cpn.ParamName.from_full_name(par_names_.pop())
    while par_names_:
        name = name.intersection(cpn.ParamName.from_full_name(par_names_.pop()))
    return name


def write_statistics(result, path):
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
