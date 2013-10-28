"""
Created on Mar 30, 2011

@author: guillaume
"""

__updated__ = "2013-10-17"

# Standard libraries
import os

# Specialized libraries
import scipy as sc
import scipy.stats as st
import scipy.linalg as la
import scipy.signal as si
import scipy.interpolate as ip

import chemex.tools as tools


def read_data(cfg, working_dir, global_parameters, res_incl=None, res_excl=None):

    # Reads the path to get the intensities
    exp_data_dir = tools.normalize_path(working_dir, cfg.get('path', 'exp_data_dir'))

    data_points = list()

    experiment_name = name_experiment(global_parameters)

    for resonance_id, filename in cfg.items('data'):

        included = ((res_incl is not None and resonance_id in res_incl) or
                    (res_excl is not None and resonance_id not in res_excl) or
                    (res_incl is None and res_excl is None))

        if not included:
            continue

        parameters = dict(global_parameters)

        parameters['experiment_name'] = experiment_name
        parameters['resonance_id'] = resonance_id

        # Get the r2 values from the fuda files containing intensities
        abs_path_filename = os.path.join(exp_data_dir, filename)
        data_points += read_a_cest_profile(abs_path_filename, parameters)

    # Adjust the minimal uncertainty
    # data_points = adjust_min_int_uncertainty(data_points)

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

        if 'B1_frq' in global_parameters:
            B1_frq = float(global_parameters['B1_frq'])

        time_t1 = float(global_parameters['time_t1'])

        name = '{:s}_{:.0f}Hz_{:.0f}ms_{:.0f}MHz_{:.0f}C'.format(exp_type, B1_frq, time_t1 * 1e3, h_larmor_frq, temperature)

    return name

def read_a_cest_profile(filename, parameters):
    """Reads in the fuda file and spit out the intensities"""

    data = sc.loadtxt(filename, dtype=[('B1_offset', '<f8'), ('intensity', '<f8'), ('intensity_err', '<f8')])

    uncertainty = estimate_uncertainty(data)
#    data = find_subset(data)

    data_points = []

    intensity_ref = 1.0

    for B1_offset, intensity_val, intensity_err in data:
        if abs(B1_offset) >= 10000.0:
            intensity_ref = intensity_val

    parameters['intensity_ref'] = intensity_ref

    for B1_offset, intensity_val, intensity_err in data:

        parameters['B1_offset'] = B1_offset

        intensity_err = uncertainty

        exp_type = parameters['experiment_type'].replace('_cest', '')

        data_point = __import__(exp_type + '.data_point', globals(), locals(), ['DataPoint'], -1)

        data_points.append(data_point.DataPoint(intensity_val, intensity_err, parameters))

    return data_points


def estimate_uncertainty(data):
    """Estimates uncertainty using the baseline"""

    data.sort()
    int_list = sc.asarray([it for of, it , _er in data if abs(of) < 10000.0])

    return estimate_noise(int_list)


# def estimate_uncertainty(data):
#     """Estimates uncertainty using the baseline"""
#
#     data.sort()
#     off_list = sc.asarray([of for of, _it , _er in data if abs(of) < 10000.0])
#     int_list = sc.asarray([it for of, it , _er in data if abs(of) < 10000.0])
#
#     std_list = list()
#
#     N = 5
#
#     for i in range(len(int_list) - N):
#
#         x = sc.array(off_list[i:i + N])
#         y = sc.array(int_list[i:i + N])
#
#         n = 2
#
#         par = sc.polyfit(x, y, n)
#         y_calc = sum([p * (x ** (n - i)) for i, p in enumerate(par)])
#         std_list.append(sc.std(y - y_calc, ddof=n + 1))
#
#     return sc.median(std_list)


def adjust_min_int_uncertainty(data_int):
    """Adjusts the uncertainty of data points to the maximum of
    either the present uncertainty or the median of all the uncertainties
    """

    int_err = sc.median([data_point.err for data_point in data_int])

    new_data_int = list()

    for data_point in data_int:
        data_point.err = int_err
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


def estimate_noise(x):

    n = len(x)

    fda = [[1, -1],
           [1, -2, 1],
           [1, -3, 3, -1],
           [1, -4, 6, -4, 1],
           [1, -5, 10, -10, 5, -1],
           [1, -6, 15, -20, 15, -6, 1]]
    fda = [a_fda / la.norm(a_fda) for a_fda in fda]

    perc = sc.array([0.05] + list(sc.arange(0.1, 0.40, 0.025)))
    z = st.norm.ppf(1.0 - perc)

    sigma_est = []

    for fdai in fda:

        noisedata = si.convolve(x, fdai, mode='valid')
        ntrim = len(noisedata)

        if ntrim >= 2:

            noisedata.sort()

            p = 0.5 + sc.arange(1, ntrim + 1)
            p /= ntrim + 0.5

            q = []

            f = ip.interp1d(p, noisedata, 'linear')

            for a_perc, a_z in zip(perc, z):
                try:
                    val = (f(1.0 - a_perc) - f(a_perc)) / (2.0 * a_z)
                    q.append(val)
                except:
                    pass

            sigma_est.append(sc.median(q))

    noisevar = sc.median(sigma_est) ** 2
    noisevar /= (1.0 + 15.0 * (n + 1.225) ** -1.245)

    return sc.sqrt(noisevar)

