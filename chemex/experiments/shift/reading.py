'''
Created on Mar 30, 2011

@author: guillaume
'''

import os
import scipy as sc

from chemex import tools


def read_data(cfg, working_dir, global_parameters, res_incl=None, res_excl=None):
    '''Read the shifts'''

    # Reads the path to get the shifts
    exp_data_dir = tools.normalize_path(working_dir, cfg.get('path', 'exp_data_dir'))

    data_points = list()

    experiment_name = name_experiment(global_parameters)

    for key, val in cfg.items('data'):

        if 'file' not in key:
            continue

        parameters = dict(global_parameters)

        parameters['experiment_name'] = experiment_name

        abs_path_filename = os.path.join(exp_data_dir, val)
        data_points += read_a_shift_file(abs_path_filename, parameters, res_incl, res_excl)

    return data_points


def name_experiment(global_parameters=None):
    if global_parameters is None:
        global_parameters = dict()

    if 'experiment_name' in global_parameters:
        name = global_parameters['experiment_name'].strip().replace(' ', '_')
    else:
        exp_type = global_parameters['experiment_type']
        h_larmor_frq = float(global_parameters['h_larmor_frq'])
        temperature = float(global_parameters['temperature'])

        name = '{:s}_{:.0f}MHz_{:.0f}C'.format(exp_type, h_larmor_frq, temperature)

    return name


def read_a_shift_file(filename, parameters, res_incl=None, res_excl=None):
    """Reads in the fuda file and spit out the intensities"""

    data = sc.loadtxt(filename, dtype=[('resonance_id', 'S10'), ('shift_ppb', 'f8'), ('shift_ppb_err', 'f8')])

    data_points = list()

    exp_type = parameters['experiment_type'].replace('_shift', '')
    data_point = __import__(exp_type + '.data_point', globals(), locals(), ['DataPoint'], -1)

    for resonance_id, shift_ppb, shift_ppb_err in data:

        included = (
            (res_incl is not None and resonance_id in res_incl) or
            (res_excl is not None and resonance_id not in res_excl) or
            (res_incl is None and res_excl is None)
        )

        if not included:
            continue

        parameters['resonance_id'] = resonance_id

        data_points.append(data_point.DataPoint(shift_ppb, shift_ppb_err, parameters))

    return data_points
