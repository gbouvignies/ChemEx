"""The fitting module contains the code for fitting the experimental data."""
import functools
import sys

import lmfit
from scipy import stats

from chemex import datasets
from chemex import parameters
from chemex import util
from chemex.cli import FITMETHODS

ALL_FITMETHODS = lmfit.minimizer.SCALAR_METHODS
ALL_FITMETHODS.update(
    {
        "leastsq": "least-squares using Levenberg-Marquardt",
        "least_squares": "least-squares using Trust Region Reflective algorithm",
        "brute": "grid search using the brute force method",
        "basinhopping": "basinhopping",
        "ampgo": "Adaptive Memory Programming for Global Optimization (AMPGO)",
    }
)
ALLOWED_FITMETHODS = {
    name: desc for name, desc in ALL_FITMETHODS.items() if name in FITMETHODS
}


def run_fit(fit_filename, params, data, cl_fitmethod):
    """Perform the fit."""
    util.header1("Fit")

    fit_config = util.read_cfg_file(fit_filename)

    if not fit_config.sections():
        fit_config.add_section("Standard Calculation")

    for section in fit_config.sections():
        util.header2(section)
        items = fit_config.items(section)
        parameters.set_param_status(params, items)
        clusters = find_independent_clusters(data, params)
        fitmethod = fit_config.get(section, "fitmethod", fallback=cl_fitmethod)

        if fitmethod not in ALLOWED_FITMETHODS.keys():
            exit(
                "The fitting method '{}', as specified in section ['{}'],"
                "is invalid! Please choose from:\n  {}".format(
                    fitmethod, section, list(sorted(ALLOWED_FITMETHODS.keys()))
                )
            )

        print("Fitting method: {}\n".format(ALLOWED_FITMETHODS[fitmethod]))

        for c_name, c_data, c_params in clusters:
            if len(clusters) > 1:
                print(f"[{c_name}]")

            print("Chi2 / Reduced Chi2:")

            c_func = c_data.calculate_residuals
            c_minimizer = lmfit.Minimizer(c_func, c_params)

            try:
                if fitmethod == "brute":
                    c_result = c_minimizer.minimize(method=fitmethod, keep="all")
                else:
                    c_result = c_minimizer.minimize(method=fitmethod)

            except KeyboardInterrupt:
                sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped\n")
                c_result = c_minimizer.minimize(
                    params=c_minimizer.result.params, maxfev=1
                )

            for name, param in c_result.params.items():
                params[name] = param

            print("")

        if len(clusters) > 1:
            c_func = functools.partial(data.calculate_residuals, verbose=False)
            minimizer = lmfit.Minimizer(c_func, params)
            result = minimizer.minimize(maxfev=1)
            result.params = params
        else:
            result = c_result

        print(f"Final Chi2        : {result.chisqr:.3e}")
        print(f"Final Reduced Chi2: {result.redchi:.3e}")

    if result.method != "leastsq":
        print("\nWarning: uncertainties and covariance of fitting parameters are only")
        print("         calculated when using the 'leastsq' fitting method!")

    return result


def find_independent_clusters(data, params):
    """Find clusters of datapoints that depend on disjoint sets of variables.

    For example, if the population of the minor state and the exchange
    rate are set to 'fix', chances are that the fit can be decomposed
    residue-specifically.

    """
    clusters = []

    for profile in data:

        pnames = profile.params.keys()
        pnames_vary = {
            name for name in pnames if params[name].vary and not params[name].expr
        }

        for data_cluster, pnames_cluster, pnames_vary_cluster in clusters:

            if pnames_vary & pnames_vary_cluster:
                data_cluster.append(profile)
                pnames_cluster.extend(pnames)
                pnames_vary_cluster.update(pnames_vary)
                break

        else:

            data_cluster = datasets.DataSet(profile)
            pnames_cluster = list(pnames)
            pnames_vary_cluster = pnames_vary

            clusters.append((data_cluster, pnames_cluster, pnames_vary_cluster))

    clusters_ = []

    for data_cluster, pnames_cluster, pnames_vary_cluster in clusters:

        name_cluster = parameters.ParameterName.from_full_name(
            pnames_vary_cluster.pop()
        )

        while pnames_vary_cluster:
            name_cluster = name_cluster.intersection(
                parameters.ParameterName.from_full_name(pnames_vary_cluster.pop())
            )

        params_cluster = lmfit.Parameters()

        for pname in pnames_cluster:
            param = params[pname]
            params[pname]._delay_asteval = True
            params_cluster[pname] = param

        for param in params_cluster.values():
            param._delay_asteval = False

        clusters_.append((name_cluster, data_cluster, params_cluster))

    return sorted(clusters_)


def write_statistics(result, path):
    """Write fitting statistics to a file."""

    _, ks_p_value = stats.kstest(result.residual, "norm")
    chi2_p_value = 1.0 - stats.chi2.cdf(result.chisqr, result.nfree)

    filename = path / "statistics.fit"

    with open(filename, "w") as f:
        print(f"  * {filename}")

        f.write(f"# Number of data points: {result.ndata}\n")
        f.write(f"# Number of variables: {result.nvarys}\n")
        f.write(f"# Fitting method: {result.method}\n")
        f.write(f"# Number of function evaluations: {result.nfev}\n\n")
        f.write(f"chi-square         = {result.chisqr: .5e}\n")
        f.write(f"reduced chi-square = {result.redchi: .5e}\n")
        f.write(f"chisq_p-value      = {chi2_p_value: .5e} # Chi-squared test\n")
        f.write(
            "ks_p-value         = {: .5e} # Kolmogorov-Smirnov test\n\n".format(
                ks_p_value
            )
        )
        f.write(f"Akaike Information Criterion = {result.aic: .5e}\n")
        f.write(f"Bayesian Information Criterion = {result.bic: .5e}\n")
