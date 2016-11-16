"""Read the experimental CEST data."""

import importlib
import os

import numpy as np

from chemex.experiments.cest import util


def read_profiles(path, profile_filenames, experiment_details, res_incl=None, res_excl=None):
    """Read the CEST profiles."""
    experiment_type = experiment_details['type'].split('.')
    experiment_details['name'] = name_experiment(experiment_details)
    experiment_module = importlib.import_module(
        '.'.join(['chemex.experiments',
                  experiment_type[0],
                  'profiles',
                  experiment_type[1]]))

    Profile = getattr(experiment_module, 'Profile')

    dtype = [
        ('b1_offsets', '<f8'),
        ('intensities', '<f8'),
        ('intensities_err', '<f8')
    ]

    profiles = []

    for profile_name, filename in profile_filenames.items():
        if (res_incl and profile_name not in res_incl) or (res_excl and profile_name in res_excl):
            continue
        full_path = os.path.join(path, filename)
        measurements = np.loadtxt(full_path, dtype=dtype)
        measurements.sort(order='b1_offsets')
        profile = Profile(profile_name, measurements, experiment_details)
        profiles.append(profile)

    error = experiment_details.get('error', 'file')
    if error not in {'file', 'auto'}:
        print('Warning: The \'error\' option should either be \'file\' or \'auto\'. Set to \'file\'')
        error = 'file'

    if error == 'auto':

        for profile in profiles:
            excl_reference = np.logical_not(profile.reference)
            error_value = util.estimate_noise(profile.val[excl_reference])
            profile.err = np.zeros_like(profile.err) + error_value

    ndata = sum(len(profile.val) for profile in profiles)

    return profiles, ndata


def name_experiment(experiment_details=None):
    """Generate a unique name for the experiment."""
    if experiment_details is None:
        experiment_details = dict()

    if 'name' in experiment_details:
        name = experiment_details['name'].strip().replace(' ', '_')

    else:
        exp_type = experiment_details['type'].replace('.', '_')
        h_larmor_frq = float(experiment_details['h_larmor_frq'])
        temperature = float(experiment_details['temperature'])
        b1_frq = float(experiment_details['b1_frq'])
        time_t1 = float(experiment_details['time_t1'])

        name = '{:s}_{:.0f}Hz_{:.0f}ms_{:.0f}MHz_{:.0f}C'.format(
            exp_type,
            b1_frq,
            time_t1 * 1e3,
            h_larmor_frq,
            temperature
        ).lower()

    return name
