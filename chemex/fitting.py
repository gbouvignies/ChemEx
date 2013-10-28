"""
Created on Mar 31, 2011

@author: guillaume
"""

# Standard Libraries
import sys
import ConfigParser

# Specialized Libraries
import scipy as sc
import scipy.optimize as opt

# Local Libraries
from chemex.chi2 import make_calc_residuals, calc_reduced_chi2
from chemex.writing import dump_parameters
from chemex.experiments.misc import get_par


def run_fit(fit_filename, par, par_indexes, par_fixed, data):

    # fit_par_file
    fit_par_file = ConfigParser.SafeConfigParser()

    if fit_filename:
        fit_par_file.optionxform = str
        fit_par_file.read(fit_filename)

    else:
        fit_par_file.add_section('Standard Calculation')

    for section in fit_par_file.sections():

        print(''.join(['\n', section, ' ...']))

        items = fit_par_file.items(section)
        par, par_indexes, par_fixed = fix_par(items, par, par_indexes, par_fixed)

        independent_clusters = find_independent_clusters(data, par, par_indexes, par_fixed)
        independent_clusters_no = len(independent_clusters)

        if independent_clusters_no > 1:

            par_err = list(par)

            for i, independent_cluster in enumerate(independent_clusters):

                sys.stdout.write('\r    Fit cluster ({}/{})'.format(i + 1, independent_clusters_no))
                sys.stdout.flush()

                c_data, c_par, c_par_indexes = independent_cluster
                c_par, c_par_err, _c_reduced_chi2 = local_minimization(c_par, c_par_indexes, par_fixed, c_data, verbose=False)

                for par_name in c_par_indexes:
                    par[par_indexes[par_name]] = c_par[c_par_indexes[par_name]]
                    par_err[par_indexes[par_name]] = c_par_err[c_par_indexes[par_name]]

            reduced_chi2 = calc_reduced_chi2(par, par_indexes, par_fixed, data)
            print('    Reduced chi2: {:.2e}'.format(reduced_chi2))

        else:
            par, par_err, reduced_chi2 = local_minimization(par, par_indexes, par_fixed, data)

    return par, par_err, par_indexes, par_fixed, reduced_chi2


def local_minimization(par, par_indexes, par_fixed, data, verbose=True):
    """
    Minimize the residuals using the Levenberg-Marquard algorithm.
    """

    if verbose:
        print('\nMinimization:\n')


    func = make_calc_residuals(verbose=verbose)
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
        sys.stderr.write(' -- Check that all parameters are correctly initialized.\n')
        dump_parameters(par, par_indexes, par_fixed, data)
        exit()

    if ier not in [1, 2, 3, 4]:
        print(''.join(('Optimal parameters not found: ', errmsg)))

    data_nb, par_nb = len(data), len(par)

    reduced_chi2 = calc_reduced_chi2(par, par_indexes, par_fixed, data)

    if (data_nb > par_nb) and pcov is not None:
        pcov = pcov * reduced_chi2
        par_err = sc.sqrt(sc.diag(pcov))

    else:
        pcov = sc.inf
        par_err = par

    return par, par_err, reduced_chi2


def fix_par(items, par, par_indexes, par_fixed):
    """
    Fix (or not) fit variables according to what set in the protocol file.
    """

    fitted_pars = set(par_indexes)
    fixed_pars = set(par_fixed)

    options = {'fit': (fixed_pars, fitted_pars),
               'fix': (fitted_pars, fixed_pars)}

    for par_name_1, state in items:

        par_name_1_str = par_name_1.replace(' ', '').split(',')
        departure_pool, arrival_pool = options[state]

        for par_name_2 in list(departure_pool):

            par_name_2_str = [str(_) for _ in par_name_2]

            if set(par_name_1_str) <= set(par_name_2_str):
                departure_pool.remove(par_name_2)
                arrival_pool.add(par_name_2)

    updated_par_indexes = dict()
    updated_par = list()
    updated_par_fixed = dict()

    for index, par_name in enumerate(fitted_pars):

        updated_par_indexes[par_name] = index
        updated_par.append(get_par(par_name, par, par_indexes, par_fixed))

    updated_par = sc.array(updated_par)

    for par_name in fixed_pars:

        updated_par_fixed[par_name] = get_par(par_name, par, par_indexes, par_fixed)

    return updated_par, updated_par_indexes, updated_par_fixed


def find_independent_clusters(data, par, par_indexes, par_fixed):
    """
    Finds clusters of data points that depend on independent sets of variables.
    For example, if the population of the minor state and the exchange rate are
    set to 'fix', chances are that the fit can be decomposed residue-specifically.
    """

    fixed_par_set = set(par_fixed)

    data = [([data_point],
                 (data_point.get_fitting_parameter_names()
                  | data_point.get_fixed_parameter_names())
                 - fixed_par_set)
                for data_point in data]

    clusters = []

    try:

        for data_point in data:

            merged = False

            for cluster in clusters:

                par_name_point = data_point[1]
                par_name_cluster = cluster[1]

                if par_name_point & par_name_cluster:
                    cluster[0].extend(data_point[0])
                    cluster[1].update(data_point[1])
                    merged = True
                    break

            if not merged:
                clusters.append(data_point)

    except (KeyboardInterrupt):
        exit("\n -- ChemEx killed while clustering\n")

    final_clusters = list()

    for cluster in clusters:
        cluster_data, cluster_par_names = cluster
        cluster_par = list(par[par_indexes[par_name]] for par_name in cluster_par_names)
        cluster_par_indexes = dict((par_name, i) for i, par_name in enumerate(cluster_par_names))
        final_clusters.append((cluster_data, cluster_par, cluster_par_indexes))

    return final_clusters
