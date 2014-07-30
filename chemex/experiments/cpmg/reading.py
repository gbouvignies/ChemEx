"""
Created on Mar 30, 2011

@author: guillaume
"""

__updated__ = "2013-10-17"

import os
import scipy as sc

from chemex import tools


def read_data(cfg, working_dir, global_parameters, res_incl=None, res_excl=None):
    # Reads the path to get the intensities
    exp_data_dir = tools.normalize_path(working_dir, cfg.get('path', 'exp_data_dir'))

    data_points = list()

    experiment_name = name_experiment(global_parameters)

    for resonance_id, filename in cfg.items('data'):

        included = (
            (res_incl is not None and resonance_id in res_incl) or
            (res_excl is not None and resonance_id not in res_excl) or
            (res_incl is None and res_excl is None)
        )

        if not included:
            continue

        parameters = dict(global_parameters)

        parameters['experiment_name'] = experiment_name
        parameters['resonance_id'] = resonance_id

        abs_path_filename = os.path.join(exp_data_dir, filename)
        data_points += read_a_cpmg_profile(abs_path_filename, parameters)

    # Adjust the minimal uncertainty
    data_points = adjust_min_int_uncertainty(data_points)

    # Normalize intensities
    data_points = norm_int(data_points)

    return data_points


def name_experiment(global_parameters=dict()):
    if 'experiment_name' in global_parameters:
        name = global_parameters['experiment_name'].strip().replace(' ', '_')
    else:
        exp_type = global_parameters['experiment_type']
        h_larmor_frq = float(global_parameters['h_larmor_frq'])
        temperature = float(global_parameters['temperature'])

        name = '{:s}_{:.0f}MHz_{:.0f}C'.format(exp_type, h_larmor_frq, temperature).lower()

    return name


def read_a_cpmg_profile(filename, parameters):
    """Reads in the fuda file and spit out the intensities"""

    data = sc.loadtxt(filename, dtype=[('ncyc', '<f8'), ('intensity', '<f8'), ('intensity_err', '<f8')])

    uncertainty_from_duplicates = estimate_uncertainty_from_duplicates(data)

    data_points = list()

    exp_type = parameters['experiment_type'].replace('_cpmg', '')
    data_point = __import__(exp_type + '.data_point', globals(), locals(), ['DataPoint'], -1)

    intensity_ref = 1.0

    for ncyc, intensity_val, intensity_err in data:
        if ncyc == 0.0:
            intensity_ref = intensity_val

    parameters['profile_id'] = filename
    parameters['intensity_ref'] = intensity_ref

    for ncyc, intensity_val, intensity_err in data:
        parameters['ncyc'] = int(ncyc)

        # Used to keep reference points out of the bootstrapping
        parameters['reference'] = int(ncyc) == 0

        # Calculate r2 uncertainty from intensity uncertainty
        intensity_err = max([uncertainty_from_duplicates, intensity_err])

        data_points.append(data_point.DataPoint(intensity_val, intensity_err, parameters))

    return data_points


def estimate_uncertainty_from_duplicates(data):
    """Estimates uncertainty using duplicate measurements"""

    intensity_dict = dict()

    for ncyc, intensity_val, _ in data:
        intensity_dict.setdefault(ncyc, list()).append(intensity_val)

    intensity_std = [sc.std(duplicates)
                     for duplicates in intensity_dict.itervalues()
                     if len(duplicates) > 1]

    return sc.mean(intensity_std) if intensity_std else 0.0


def adjust_min_int_uncertainty(data_int):
    """
    Adjusts the uncertainty of data points to the maximum of
    either the present uncertainty or the median of all the uncertainties

    """

    int_err_list = list()

    for data_point in data_int:
        int_err_list.append(data_point.err)

    new_data_int = list()

    for data_point in data_int:
        data_point.err = max([data_point.err, sc.median(int_err_list)])
        new_data_int.append(data_point)

    return new_data_int


def norm_int(data_int):
    """Normalize intensities relative to the intensity of the reference plane"""

    new_data_int = list()

    for data_point in data_int:
        data_point.val /= data_point.par['intensity_ref']
        data_point.err /= data_point.par['intensity_ref']
        data_point.err = abs(data_point.err)
        new_data_int.append(data_point)

    return new_data_int
