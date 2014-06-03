"""
Created on Apr 1, 2011

@author: guillaume
"""

import os
import scipy as sc
import scipy.stats as stats

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

from chemex.tools import parse_assignment

# Constants
LINNEWIDTH = 1.0


def set_lim(values, scale):
    """Provides a range that contains all the value and adds a margin."""

    v_min, v_max = min(values), max(values)
    margin = (v_max - v_min) * scale
    v_min, v_max = v_min - margin, v_max + margin

    return v_min, v_max


def group_data(dataset):
    """Groups the data resonance specifically"""

    grouped_dataset = dict()

    for a_data in dataset:
        resonance_id = a_data.par['resonance_id']

        assignment = parse_assignment(resonance_id)
        index = assignment[0][0]

        grouped_dataset.setdefault((index, resonance_id), []).append(a_data)

    return grouped_dataset


def make_val_for_plot(residue_dataset, par, par_names, par_fixed, filename_out):
    """Creates the arrays that will be used to plot one profile"""

    mag_ref = sc.mean([point.val for point in residue_dataset if point.par['ncyc'] == 0])

    values = []

    for point in residue_dataset:

        ncyc = point.par['ncyc']

        if ncyc:
            time_t2 = point.par['time_t2']
            mag_exp = point.val
            mag_err = point.err

            frq = ncyc / time_t2

            r2_exp = -sc.log(point.val / mag_ref) / time_t2

            mag_dist = sc.random.normal(mag_exp, mag_err, 1000.0)
            r2_dist = -sc.log(mag_dist / mag_ref) / time_t2
            r2_err_down = abs(stats.scoreatpercentile(r2_dist, 15.9) - r2_exp)
            r2_err_up = abs(stats.scoreatpercentile(r2_dist, 84.1) - r2_exp)

            point.calc_val(par, par_names, par_fixed)
            r2_cal = -sc.log(point.cal / mag_ref) / time_t2

            values.append((frq, r2_exp, r2_err_down, r2_err_up, r2_cal))

    for a_xe, a_ye, a_ede, a_eue, a_yc in sorted(values):
        filename_out.write(
            "{:10s} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                residue_dataset[0].par['resonance_id'].upper(), a_xe, a_ye, a_ede, a_eue, a_yc
            )
        )

    xe, ye, ede, eue, yc = zip(*sorted(values))

    return xe, ye, ede, eue, yc


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    """Plot dispersion profiles and write a multi-page pdf file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.par['experiment_name']
        datasets.setdefault(experiment_name, list()).append(data_point)

    for experiment_name, dataset in datasets.iteritems():

        ###### Matplotlib ######

        filename = ''.join([experiment_name, '.pdf'])
        filename = os.path.join(output_dir, filename)

        filename_calc = ''.join([experiment_name, '.fit'])
        filename_calc = os.path.join(output_dir, filename_calc)

        print("     * {}".format(filename))

        ########################

        grouped_dataset = group_data(dataset)

        pdf = PdfPages(filename)

        with open(filename_calc, 'w') as f:

            for (_index, resonance_id), residue_dataset in sorted(grouped_dataset.iteritems()):
                out = make_val_for_plot(residue_dataset, par, par_names, par_fixed, f)
                xe, ye, ede, eue, yc = out

                ###### Matplotlib ######

                fig = plt.figure(1, linewidth=LINNEWIDTH, figsize=(6, 4.5))
                ax = fig.add_subplot(111)

                ########################

                ax.plot(
                    xe, yc, '-',
                    color='0.5',
                    linewidth=LINNEWIDTH
                )

                ax.errorbar(
                    xe, ye, yerr=[ede, eue], fmt='ro',
                    linewidth=LINNEWIDTH,
                    markersize=5.0,
                    markerfacecolor='w',
                    markeredgewidth=LINNEWIDTH,
                    markeredgecolor='r'
                )

                upper_vals = list(sc.array(ye) + sc.array(eue))
                lower_vals = list(sc.array(ye) - sc.array(ede))
                all_vals = upper_vals + lower_vals + list(yc)

                xmin, xmax = set_lim(xe, 0.05)
                ymin, ymax = set_lim(all_vals, 0.05)

                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)

                ax.tick_params(length=3, top=True, right=False, labelsize=10)

                ax.xaxis.set_major_locator(MaxNLocator(6))
                ax.yaxis.set_major_locator(MaxNLocator(6))

                ax.set_xlabel(r'$\mathregular{\nu_{CPMG}}$' + ' (Hz)')
                ax.set_ylabel(r'$\mathregular{R_{2,eff}}$ ($\mathregular{s^{-1}}$)')

                ax.set_title('{:s}'.format(resonance_id.upper()))

                fig.tight_layout()

                ########################

                pdf.savefig()
                plt.close()

                ########################

        pdf.close()

    return
