import os

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator, NullFormatter

from chemex import peaks

# colors
dark_gray = '0.13'
red500 = '#F44336'
red200 = '#EF9A9A'

TWO_PI = 2.0 * np.pi


def sigma_estimator(x):
    """Estimates standard deviation using median to exclude outliers. Up to
    50% can be bad."""

    return np.median([np.median(abs(xi - np.asarray(x))) for xi in x]) * 1.1926


def set_lim(values, scale):
    """Provides a range that contains all the value and adds a margin."""

    v_min, v_max = min(values), max(values)
    margin = (v_max - v_min) * scale
    v_min, v_max = v_min - margin, v_max + margin

    return v_min, v_max


def group_data(dataset):
    """Groups the data resonance specifically"""

    data_grouped = dict()

    for profile in dataset:
        name = profile.profile_name
        peak = peaks.Peak(name)
        data_grouped[peak] = profile

    return data_grouped


def compute_profiles(data_grouped, params):
    """Creates the arrays that will be used to plot one profile"""

    profiles = {}

    for peak, profile in data_grouped.items():
        mask = profile.b1_offsets > -10000.0
        mask_ref = np.logical_not(mask)

        val_ref = np.mean(profile.val[mask_ref])

        b1_ppm_exp = profile.b1_offsets_to_ppm()[mask]
        mag_cal = profile.calculate_profile(params)[mask] / val_ref
        mag_exp = profile.val[mask] / val_ref
        mag_err = profile.err[mask] / np.absolute(val_ref)

        b1_offsets_min, b1_offsets_max = set_lim(profile.b1_offsets[mask], 0.02)
        b1_offsets = np.linspace(b1_offsets_min, b1_offsets_max, 500)

        b1_ppm_fit = profile.b1_offsets_to_ppm(b1_offsets)
        mag_fit = profile.calculate_profile(params, b1_offsets) / val_ref

        profiles[peak] = b1_ppm_exp, mag_cal, mag_exp, mag_err, b1_ppm_fit, mag_fit

    return profiles


def write_profile(name, b1_ppm_fit, mag_fit, file_txt):
    for b1_ppm_cal, mag_cal in zip(b1_ppm_fit, mag_fit):
        file_txt.write(
            "{:10s} {:8.3f} {:8.3f}\n".format(
                name.upper(), b1_ppm_cal, mag_cal
            )
        )


def plot_data(data, params, output_dir='./'):
    """Plot cest profiles and write a pdf file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.experiment_name
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

        profiles = compute_profiles(data_grouped, params)

        with PdfPages(name_pdf) as file_pdf, open(name_txt, 'w') as file_txt:

            for peak in sorted(profiles):
                b1_ppm, mag_cal, mag_exp, mag_err, b1_ppm_fit, mag_fit = profiles[peak]

                write_profile(peak.assignment, b1_ppm_fit, mag_fit, file_txt)

                ###### Matplotlib ######

                fig = plt.figure(1)

                gs = gsp.GridSpec(2, 1, height_ratios=[1, 4])

                ax1 = plt.subplot(gs[0])
                ax2 = plt.subplot(gs[1])

                ax1.axhline(0, color='black', alpha=0.87)
                ax2.axhline(0, color='black', alpha=0.87)

                ########################

                ax2.plot(
                    b1_ppm_fit,
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

                deltas = np.asarray(mag_exp) - np.asarray(mag_cal)
                max_val = max(np.absolute(set_lim(deltas, 0.1))) + max(mag_err)
                power10 = int(np.log10(max_val))
                deltas /= 10 ** power10
                mag_err = np.array(mag_err) / 10 ** power10
                sigma = sigma_estimator(deltas)

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    1.0 * sigma * np.asarray([-1.0, 1.0, 1.0, -1.0]),
                    fc='black',
                    alpha=0.12,
                    ec='none'
                )

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    2.0 * sigma * np.asarray([-1.0, 1.0, 1.0, -1.0]),
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
                    zorder=100,
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

                ax1.set_title('{:s}'.format(peak.assignment.upper()))
                ax1.set_ylabel(r''.join([
                    r'$\mathregular{Resid. \ x10^{',
                    r'{:d}'.format(power10),
                    r'}}$'
                ]))

                ########################

                fig.set_tight_layout(True)

                ########################

                file_pdf.savefig()
                plt.close()

                ########################

    return
