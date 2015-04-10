import os

import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
from matplotlib.ticker import MaxNLocator, NullFormatter
from matplotlib.backends.backend_pdf import PdfPages

from chemex.parsing import parse_assignment


dark_gray = '0.13'
red500 = '#F44336'
red200 = '#EF9A9A'


def sigma_estimator(x):
    """ Estimates standard deviation using median to exclude outliers. Up to
    50% can be bad """

    return sp.median([sp.median(abs(xi - sp.asarray(x))) for xi in x]) * 1.1926


def set_lim(values, scale):
    """Provides a range that contains all the value and adds a margin."""

    v_min, v_max = min(values), max(values)
    margin = (v_max - v_min) * scale
    v_min, v_max = v_min - margin, v_max + margin

    return v_min, v_max


def group_data(dataset):
    """Groups the data resonance specifically"""

    data_grouped = dict()

    for data_pt in dataset:
        id = data_pt.par['resonance_id']

        assignment = parse_assignment(id)
        index = int(assignment[0][0])

        data_grouped.setdefault((index, id), []).append(data_pt)

    return data_grouped


def compute_profiles(data_grouped, par, par_names, par_fixed):
    """Creates the arrays that will be used to plot one profile"""

    profiles = {}

    for (index, resonance_id), profile in data_grouped.items():

        profile_exp = []
        b1_offset_min = +1e16
        b1_offset_max = -1e16

        for data_pt in profile:

            b1_offset = data_pt.par['b1_offset']
            ppm_to_rads = data_pt.par['ppm_to_rads']
            carrier_ppm = data_pt.par['carrier']

            if abs(b1_offset) < 1e4:
                b1_ppm = (2.0 * sp.pi * b1_offset) / ppm_to_rads + carrier_ppm

                mag_cal = data_pt.cal
                mag_exp = data_pt.val
                mag_err = data_pt.err

                profile_exp.append([b1_ppm, mag_cal, mag_exp, mag_err])

                b1_offset_min = min(b1_offset_min, b1_offset)
                b1_offset_max = max(b1_offset_max, b1_offset)

        b1_ppm_exp, mag_cal, mag_exp, mag_err = zip(*sorted(profile_exp))

        profile_cal = []

        b1_offset_min, b1_offset_max = set_lim(
            [b1_offset_min, b1_offset_max], 0.02
        )

        data_pt = profile[0]

        for b1_offset in sp.linspace(b1_offset_min, b1_offset_max, 500):
            ppm_to_rads = data_pt.par['ppm_to_rads']
            carrier_ppm = data_pt.par['carrier']
            b1_ppm = (2.0 * sp.pi * b1_offset) / ppm_to_rads + carrier_ppm

            data_pt.update_b1_offset(b1_offset)
            data_pt.calc_val(par, par_names, par_fixed)

            profile_cal.append([b1_ppm, data_pt.cal])

        b1_ppm_fit, mag_fit = zip(*sorted(profile_cal))

        profiles.setdefault((index, resonance_id), []).append(
            [b1_ppm_exp, mag_cal, mag_exp, mag_err, b1_ppm_fit, mag_fit]
        )

    return profiles


def write_profile(resonance_id, b1_ppm_fit, mag_fit, file_txt):
    for b1_ppm_cal, mag_cal in zip(b1_ppm_fit, mag_fit):
        file_txt.write(
            "{:10s} {:8.3f} {:8.3f}\n".format(
                resonance_id.upper(), b1_ppm_cal, mag_cal
            )
        )


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    """Plot cest profiles and write a pdf file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.par['experiment_name']
        datasets.setdefault(experiment_name, []).append(data_point)

    for experiment_name, dataset in datasets.items():

        # ##### Matplotlib ######
        name_pdf = ''.join([experiment_name, '.pdf'])
        name_pdf = os.path.join(output_dir, name_pdf)

        name_txt = ''.join([experiment_name, '.fit'])
        name_txt = os.path.join(output_dir, name_txt)

        print("  * {} [.fit]".format(name_pdf))

        # #######################

        data_grouped = group_data(dataset)

        profiles = compute_profiles(
            data_grouped, par, par_names, par_fixed
        )

        with PdfPages(name_pdf) as file_pdf, open(name_txt, 'w') as file_txt:

            for (_index, resonance_id), profile in sorted(profiles.items()):
                b1_ppm, mag_cal, mag_exp, mag_err, b1_ppm_fit, mag_fit = \
                profile[0]

                write_profile(resonance_id, b1_ppm_fit, mag_fit, file_txt)

                ###### Matplotlib ######

                # fig = plt.figure(1)
                fig = plt.figure(1)

                gs = gsp.GridSpec(2, 1, height_ratios=[1, 4])

                ax1 = plt.subplot(gs[0])
                ax2 = plt.subplot(gs[1])

                ax1.axhline(0, color='black', alpha=0.87)
                ax2.axhline(0, color='black', alpha=0.87)

                ########################

                ax2.plot(b1_ppm_fit,
                         mag_fit,
                         linestyle='-',
                         color=red200,
                )

                ax2.plot(
                    b1_ppm,
                    mag_exp,
                    'o',
                    color=red500,
                )

                xmin, xmax = set_lim(b1_ppm_fit, 0.05)
                mags = list(mag_exp) + list(mag_fit)
                ymin, ymax = set_lim(mags, 0.10)

                ax2.set_xlim(xmin, xmax)
                ax2.set_ylim(ymin, ymax)

                ax2.invert_xaxis()

                ax2.xaxis.set_major_locator(MaxNLocator(9))
                # ax2.yaxis.set_major_locator(MaxNLocator(6))

                ax2.set_xlabel(r'$\mathregular{B_1 \ position \ (ppm)}$')
                ax2.set_ylabel(r'$\mathregular{I/I_0}$')

                ########################

                deltas = sp.asarray(mag_exp) - sp.asarray(mag_cal)
                max_val = max(sp.absolute(set_lim(deltas, 0.1))) + max(mag_err)
                power10 = int(sp.log10(max_val))
                deltas /= 10 ** power10
                mag_err = sp.array(mag_err) / 10 ** power10
                sigma = sigma_estimator(deltas)

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    1.0 * sigma * sp.asarray([-1.0, 1.0, 1.0, -1.0]),
                    fc='black',
                    alpha=0.12,
                    ec='none'
                )

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    2.0 * sigma * sp.asarray([-1.0, 1.0, 1.0, -1.0]),
                    fc='black',
                    alpha=0.12,
                    ec='none'
                )

                ax1.errorbar(
                    b1_ppm,
                    deltas,
                    mag_err,
                    fmt='o',
                    color=red500,
                )

                rmin, rmax = set_lim(deltas, 0.1)
                rmin = min([-3 * sigma, rmin - max(mag_err)])
                rmax = max([+3 * sigma, rmax + max(mag_err)])

                ax1.set_xlim(xmin, xmax)
                ax1.set_ylim(rmin, rmax)

                ax1.invert_xaxis()

                ax1.xaxis.set_major_locator(MaxNLocator(9))
                ax1.yaxis.set_major_locator(MaxNLocator(5))

                ax1.xaxis.set_major_formatter(NullFormatter())

                ax1.set_title('{:s}'.format(resonance_id.upper()))
                ax1.set_ylabel(r''.join([
                    r'$\mathregular{Resid. \ x10^{',
                    r'{:d}'.format(power10),
                    r'}}$'
                ]))

                ########################

                fig.tight_layout()

                ########################

                file_pdf.savefig()
                plt.close()

                ########################

    return
