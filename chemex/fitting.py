from __future__ import print_function

import copy
import itertools
import sys

import lmfit

from chemex import parameters
from chemex import util
from chemex.experiments import datasets

product = itertools.product


def run_fit(fit_filename, params, data):
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

        for index, independent_cluster in enumerate(independent_clusters, 1):

            print(
                "\nChi2 / Reduced Chi2 (cluster {}/{}):"
                    .format(index, independent_clusters_no)
            )

            c_data, c_params = independent_cluster

            func = c_data.calculate_residuals

            minimizer = lmfit.Minimizer(func, c_params)

            try:
                minimizer.leastsq()
            except KeyboardInterrupt:
                sys.stderr.write("\n -- Keyboard Interrupt: minimization "
                                 "stopped\n")

            for name, param in c_params.items():
                params[name] = param

        print("\nFinal Chi2        : {:.3e}"
              .format(data.calculate_chisq(params)))
        print("Final Reduced Chi2: {:.3e}"
              .format(data.calculate_redchi(params)))

    return params


def find_independent_clusters(data, params):
    """
    Finds clusters of data points that depend on independent sets of variables.
    For example, if the population of the minor state and the exchange rate are
    set to 'fix', chances are that the fit can be decomposed
    residue-specifically.
    """

    clusters = []

    for profile in data:

        params_profile_all = lmfit.Parameters(
            (name, params[name]) for name in profile.map_names.values()
        )

        params_profile_fit = set([
                                     name for name, param in params_profile_all.items() if param.vary
                                     ])

        for data_cluster, params_cluster in clusters:

            params_cluster_fit = set([
                name for name, param in params_cluster.items() if param.vary
            ])

            if params_profile_fit.intersection(params_cluster_fit):
                data_cluster.append(profile)
                params_cluster.update(params_profile_all)
                break

        else:
            data_cluster = datasets.DataSet(profile)
            params_cluster = copy.deepcopy(params_profile_all)
            clusters.append((data_cluster, params_cluster))

    return clusters

