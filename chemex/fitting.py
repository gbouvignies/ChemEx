"""The fitting module contains the code for fitting the experimental data."""

import sys

import lmfit

from chemex import datasets, parameters, util


def run_fit(fit_filename, params, data, cl_fitmethod):
    """Perform the fit."""
    util.header1("Fit")

    fit_config = util.read_cfg_file(fit_filename)

    if not fit_config.sections():
        fit_config.add_section('Standard Calculation')

    for section in fit_config.sections():
        util.header2(section)

        items = fit_config.items(section)

        parameters.set_param_status(params, items)

        independent_clusters = find_independent_clusters(data, params)
        independent_clusters_no = len(independent_clusters)

        fitmethod = fit_config.get(section, 'fitmethod', fallback=cl_fitmethod)
        allowed_fitmethods = lmfit.minimizer.SCALAR_METHODS
        allowed_fitmethods.update({
            'leastsq': 'least-squares using Levenberg-Marquardt',
            'least_squares': 'least-squares using Trust Region Reflective algorithm'
        })

        if fitmethod not in allowed_fitmethods.keys():
            exit(
                "The fitting method \'{}\', as specified in section \'{}\', is invalid! Please choose from:\n  {}".
                format(fitmethod, section, list(sorted(allowed_fitmethods.keys()))))

        print("Fitting method: {}\n".format(allowed_fitmethods[fitmethod]))

        for c_name, c_data, c_params in independent_clusters:
            if independent_clusters_no > 1:
                print("[{}]".format(c_name))

            print("Chi2 / Reduced Chi2:")

            func = c_data.calculate_residuals
            minimizer = lmfit.Minimizer(func, c_params)

            try:
                result = minimizer.minimize(method=fitmethod, params=c_params)

            except KeyboardInterrupt:
                result = minimizer.result
                sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped\n")

            for name, param in result.params.items():
                params[name] = param

            print('')

        print("Final Chi2        : {:.3e}".format(data.calculate_chisq(params)))
        print("Final Reduced Chi2: {:.3e}".format(data.calculate_redchi(params)))

        if result.method != 'leastsq':
            print("\nWarning: uncertainties and covariance of fitting parameters are only")
            print("         calculated when using the \'leastsq\' fitting method!")

    return params


def find_independent_clusters(data, params):
    """Find clusters of datapoints that depend on disjoint sets of variables.

    For example, if the population of the minor state and the exchange
    rate are set to 'fix', chances are that the fit can be decomposed
    residue-specifically.

    """
    clusters = []

    for profile in data:
        params_profile = lmfit.Parameters()
        for name, param in profile.default_params.items():
            params_profile[name] = lmfit.Parameter(name)

        names_vary_profile = []

        for name in profile.default_params:
            params_profile[name] = params[name]
            if params[name].vary and not params[name].expr:
                names_vary_profile.append(name)

        for name_cluster, data_cluster, params_cluster in clusters:
            names_vary_shared = [
                name for name, param in params_cluster.items()
                if param.vary and name in names_vary_profile
            ]

            if names_vary_shared:
                data_cluster.append(profile)
                for name in params_profile:
                    if name not in params_cluster:
                        params_cluster[name] = lmfit.Parameter(name)
                params_cluster.update(params_profile)
                for name in names_vary_shared:
                    name_cluster = name_cluster.intersection(
                        parameters.ParameterName.from_full_name(name))
                break

        else:
            data_cluster = datasets.DataSet(profile)
            params_cluster = params_profile

            name_cluster = parameters.ParameterName.from_full_name(names_vary_profile[0])

            for name in names_vary_profile:
                name_cluster = name_cluster.intersection(
                    parameters.ParameterName.from_full_name(name))

            clusters.append((name_cluster, data_cluster, params_cluster))

    return sorted(clusters)
