from __future__ import print_function

import itertools
import os.path
import ConfigParser
import sys

import scipy as sp
import scipy.optimize as opt

from chemex import utils
from chemex import chi2
from chemex import writing
from chemex.experiments import misc


product = itertools.product


def run_fit(fit_filename, par, par_indexes, par_fixed, data):
    fit_par_file = ConfigParser.SafeConfigParser()

    utils.header1("Fit")

    if fit_filename:

        if os.path.isfile(fit_filename):
            try:
                fit_par_file.read(fit_filename)
            except ConfigParser.MissingSectionHeaderError:
                exit(
                    'You are missing a section heading (default?) in {'
                    ':s}\n'.format(
                        fit_filename))
            except ConfigParser.ParsingError:
                exit(
                    'Having trouble reading your parameter file, have you '
                    'forgotten \'=\' signs?\n{:s}'
                    .format(sys.exc_info()[1]))
        else:
            exit("The file \'{}\' is empty or does not exist!\n".format(
                fit_filename))

    if not fit_par_file.sections():
        fit_par_file.add_section('Standard Calculation')

    for section in fit_par_file.sections():

        utils.header2(section)

        items = fit_par_file.items(section)
        par, par_indexes, par_fixed = fix_par(items, par, par_indexes,
                                              par_fixed)

        independent_clusters = find_independent_clusters(data, par,
                                                         par_indexes,
                                                         par_fixed)
        independent_clusters_no = len(independent_clusters)

        if independent_clusters_no > 1:

            par_err = list(par)

            for i, independent_cluster in enumerate(independent_clusters, 1):

                print('\nChi2 / Reduced Chi2 (cluster {}/{}):'
                      .format(i, independent_clusters_no))

                c_data, c_par, c_par_indexes = independent_cluster
                c_par, c_par_err, _c_reduced_chi2 = local_minimization(
                    c_par,
                    c_par_indexes,
                    par_fixed,
                    c_data,
                    verbose=True
                )

                for par_name in c_par_indexes:
                    par[par_indexes[par_name]] = c_par[c_par_indexes[par_name]]
                    par_err[par_indexes[par_name]] = c_par_err[
                        c_par_indexes[par_name]]

        else:
            print("\nChi2 / Reduced Chi2:")
            par, par_err, reduced_chi2 = local_minimization(par, par_indexes,
                                                            par_fixed, data)

        print("\nFinal Chi2        : {:.3e}".format(
            chi2.calc_chi2(par, par_indexes, par_fixed, data)))
        print("Final Reduced Chi2: {:.3e}".format(
            chi2.calc_reduced_chi2(par, par_indexes, par_fixed, data)))

    return par, par_err, par_indexes, par_fixed


def local_minimization(par, par_indexes, par_fixed, data, verbose=True):
    """
    Minimize the residuals using the Levenberg-Marquard algorithm.
    """

    func = chi2.make_calc_residuals(verbose=verbose)
    args = (par_indexes, par_fixed, data)

    try:
        out = opt.leastsq(func, par, args=args,
                          full_output=True,
                          ftol=1e-9,
                          xtol=1e-9,
                          maxfev=100000,
                          epsfcn=1e-10,
                          factor=0.1)
        par, pcov, _infodict, errmsg, ier = out

    except TypeError:
        sys.stderr.write(' -- Error encountered during minimization:\n')
        sys.stderr.write(' ----> {:s}\n'.format(sys.exc_info()[1]))
        sys.stderr.write(
            ' -- Check that all parameters are correctly initialized.\n')
        writing.dump_parameters(par, par_indexes, par_fixed, data)
        exit()

    if ier not in [1, 2, 3, 4]:
        print(''.join(('Optimal parameters not found: ', errmsg)))

    data_nb, par_nb = len(data), len(par)

    reduced_chi2 = chi2.calc_reduced_chi2(par, par_indexes, par_fixed, data)

    if (data_nb > par_nb) and pcov is not None:
        pcov = pcov * reduced_chi2
        par_err = sp.sqrt(sp.diag(pcov))

    else:
        par_err = par

    return par, par_err, reduced_chi2


def fix_par(items, par, par_indexes, par_fixed):
    """
    Fix (or not) fit variables according to what set in the protocol file.
    """

    fitted_pars = set(par_indexes)
    params_fix = set(par_fixed)

    options = {
        'fit': (params_fix, fitted_pars),
        'fix': (fitted_pars, params_fix)
    }

    for param_1, state in items:

        param_1_str = param_1.replace(' ', '').split(',')
        pool_start, pool_end = options[state]

        for param_2 in list(pool_start):

            param_2_str = [str(_) for _ in param_2]

            if set(param_1_str) <= set(param_2_str):
                pool_start.remove(param_2)
                pool_end.add(param_2)

    par_indexes_updated = dict()
    par_updated = list()
    par_fixed_updated = dict()

    for index, par_name in enumerate(fitted_pars):
        par_indexes_updated[par_name] = index
        par_updated.append(misc.get_par(par_name, par, par_indexes, par_fixed))

    par_updated = sp.array(par_updated)

    for par_name in params_fix:
        par_fixed_updated[par_name] = misc.get_par(par_name, par, par_indexes,
                                                   par_fixed)

    return par_updated, par_indexes_updated, par_fixed_updated


def get_params_fit(data_pt, params_fix):
    """Returns the fitted parameters a specific data point depends on."""

    data_pt_params = (
        (
            data_pt.get_fitting_parameter_names() |
            data_pt.get_fixed_parameter_names()
        ) - params_fix
    )

    return data_pt_params


def find_independent_clusters(data, par, par_indexes, par_fixed):
    """
    Finds clusters of data points that depend on independent sets of variables.
    For example, if the population of the minor state and the exchange rate are
    set to 'fix', chances are that the fit can be decomposed
    residue-specifically.
    """

    params_fix = set(par_fixed)

    clusters = []

    for data_pt in data:

        params_pt = get_params_fit(data_pt, params_fix)

        for data_cluster, params_cluster in clusters:

            if params_pt & params_cluster:
                data_cluster.append(data_pt)
                params_cluster.update(params_pt)
                break

        else:
            clusters.append(([data_pt], params_pt))

    clusters_final = list()

    for data_cluster, params_cluster in clusters:

        par_cluster = []
        par_indexes_cluster = {}

        for index, param in enumerate(params_cluster):
            par_cluster.append(par[par_indexes[param]])
            par_indexes_cluster[param] = index

        clusters_final.append((data_cluster, par_cluster, par_indexes_cluster))

    return clusters_final

