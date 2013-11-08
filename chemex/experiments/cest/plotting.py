"""
Created on Apr 1, 2011

@author: guillaume
"""

# Import
import os

# Specialized Libraries
from scipy import linspace, asarray, median, pi

import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp

from matplotlib.ticker import MaxNLocator, NullFormatter
from matplotlib.backends.backend_pdf import PdfPages

from chemex.tools import parse_assignment

# Constants
linewidth = 1.0


def sigma_estimator(x):
    """ Estimates standard deviation using median to exclude outliers. Up to 50% can be bad """
    return median([median(abs(xi - asarray(x))) for xi in x]) * 1.1926


def set_lim(x, scale):
    xmin, xmax = min(x), max(x)
    xmargin = (xmax - xmin) * scale
    xmin, xmax = xmin - xmargin, xmax + xmargin

    return xmin, xmax


def group_data(dataset):
    grouped_dataset = dict()

    for a_data in dataset:

        b1_offset = a_data.par['b1_offset']
        resonance_id = a_data.par['resonance_id']

        assignment = parse_assignment(resonance_id)
        index = assignment[0][0]

        if abs(b1_offset) < 10000.0:
            grouped_dataset.setdefault((index, resonance_id), []).append(a_data)

    return grouped_dataset


def make_val_for_plot(residue_dataset, par, par_names, par_fixed, fileout):
    values = []

    for point in residue_dataset:
        point.calc_val(par, par_names, par_fixed)

        b1_offset = point.par['b1_offset']
        ppm_to_rads = point.par['ppm_to_rads']
        b1_offset_ppm = (2.0 * pi * b1_offset) / ppm_to_rads

        values.append((b1_offset, b1_offset_ppm, point.val, point.err, point.cal))

    of, xe, ye, ee, yc = zip(*sorted(values))

    a_point = residue_dataset[0]

    values = []

    of_min, of_max = set_lim(of, 0.02)

    for b1_offset in linspace(of_min, of_max, 500):
        ppm_to_rads = a_point.par['ppm_to_rads']
        b1_offset_ppm = (2.0 * pi * b1_offset) / ppm_to_rads

        a_point.update_b1_offset(b1_offset)
        a_point.calc_val(par, par_names, par_fixed)

        values.append((b1_offset, b1_offset_ppm, a_point.cal))

        fileout.write(''.join([str(a_point), '\n']))

    of, xf, yf = zip(*sorted(values))

    return xe, ye, ee, yc, xf, yf


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    """Plot dispersion profiles and write a pdf file"""

    datasets = dict()

    for data_point in data:
        if 'cest' in data_point.par['experiment_type']:
            experiment_name = data_point.par['experiment_name']
            datasets.setdefault(experiment_name, []).append(data_point)

    for experiment_name, dataset in datasets.iteritems():

        ###### Matplotlib ######
        filename = ''.join([experiment_name, '.pdf'])
        filename = os.path.join(output_dir, filename)

        filename_calc = ''.join([experiment_name, '.fit'])
        filename_calc = os.path.join(output_dir, filename_calc)

        ########################

        grouped_dataset = group_data(dataset)

        pdf = PdfPages(filename)

        with open(filename_calc, 'w') as f:

            for (_index, resonance_id), residue_dataset in sorted(grouped_dataset.iteritems()):
                out = make_val_for_plot(residue_dataset, par, par_names, par_fixed, f)
                xe, ye, ee, yc, xf, yf = out

                ###### Matplotlib ######

                fig = plt.figure(1, linewidth=linewidth, figsize=(6, 6))

                gs = gsp.GridSpec(2, 1, height_ratios=[1, 4])

                ax1 = plt.subplot(gs[0])
                ax2 = plt.subplot(gs[1])

                ########################

                ax2.plot(
                    xf, yf, '-',
                    color='0.5',
                    linewidth=linewidth
                )

                ax2.plot(
                    xe, ye, 'ro',
                    linewidth=linewidth,
                    markersize=5.0,
                    markerfacecolor='w',
                    markeredgewidth=linewidth,
                    markeredgecolor='r'
                )

                ymin, ymax = set_lim(ye, 0.05)
                ax2.set_ylim(ymin, ymax)
                ax2.tick_params(length=3, top=True, right=False, labelsize=10)
                ax2.xaxis.set_major_locator(MaxNLocator(9))
                ax2.yaxis.set_major_locator(MaxNLocator(6))
                ax2.invert_xaxis()
                ax2.set_xlabel(r'$\mathregular{B_1}$' + ' offset (ppm)')
                ax2.set_ylabel(r'$\mathregular{I/I_0}$')

                ########################

                deltas = asarray(ye) - asarray(yc)
                sigma = sigma_estimator(deltas)

                xmin, xmax = set_lim(xf, 0.05)

                rmin, rmax = set_lim(deltas, 0.1)
                rmin = min([-4 * sigma, rmin - max(ee)])
                rmax = max([+4 * sigma, rmax + max(ee)])

                ax1.plot([xmin, xmax], [0.0, 0.0], 'k-', linewidth=linewidth)

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    2.0 * sigma * asarray([-1.0, 1.0, 1.0, -1.0]),
                    fc='#EEEEEE', ec='none'
                )

                ax1.fill(
                    (xmin, xmin, xmax, xmax),
                    1.0 * sigma * asarray([-1.0, 1.0, 1.0, -1.0]),
                    fc='#CCCCCC', ec='none'
                )

                ax1.errorbar(
                    xe, deltas, ee,
                    fmt='ro',
                    linewidth=linewidth,
                    markersize=5.0,
                    markerfacecolor='w',
                    markeredgewidth=linewidth,
                    markeredgecolor='r',
                    capsize=2.5
                )

                ax1.set_ylim(rmin, rmax)
                ax1.xaxis.set_major_locator(MaxNLocator(9))
                ax1.yaxis.set_major_locator(MaxNLocator(4))
                ax1.xaxis.set_major_formatter(NullFormatter())
                ax1.tick_params(length=3, top=False, right=False, labelsize=10)
                ax1.invert_xaxis()

                ########################

                ax2.set_xlim(xmin, xmax)
                ax1.set_xlim(xmin, xmax)
                ax1.set_title('{:s}'.format(resonance_id))
                ax1.set_ylabel('Residual')

                fig.tight_layout()

                ########################

                pdf.savefig()
                plt.close()

                ########################

        pdf.close()
    return
