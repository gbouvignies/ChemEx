'''
Created on Apr 1, 2011

@author: guillaume
'''

# Standard Libraries
import os
import re
from sys import stdout, stderr

# Specialized Libraries
import scipy as sc
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

# Constants
linewidth = 1.5


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    '''Plot dispersion profiles and write a pdf file'''

    datasets = dict()

    for data_point in data:
        if 'cpmg' in data_point.par['experiment_type']:
            experiment_name = data_point.par['experiment_name']
            datasets.setdefault(experiment_name, list()).append(data_point)

    for experiment_name, dataset in datasets.iteritems():

        filename = ''.join([experiment_name, '.pdf'])
        filename = os.path.join(output_dir, filename)

        exp_ints = dict()
        cal_ints = dict()

        for a_data_point in dataset:
            resonance_id = a_data_point.par['resonance_id']
            ncyc = a_data_point.par['ncyc']
            time_t2 = a_data_point.par['time_t2']

            a_data_point.calc_val(par, par_names, par_fixed)

            exp_ints.setdefault(resonance_id, list()).append((ncyc, time_t2, a_data_point.val, a_data_point.err))
            cal_ints.setdefault(resonance_id, list()).append((ncyc, time_t2, a_data_point.cal))

        exp_r2s = dict()

        pdf = PdfPages(filename)

        index_id = [(int(re.findall(r'\d+', resonance_id)[0]), resonance_id)
                    for resonance_id, ints in exp_ints.iteritems()]

        for _, resonance_id in sorted(index_id):

            ints = exp_ints[resonance_id]
            ints.sort()
            int_ref = ints[0][2]

            exp_nus, exp_r2s, down_r2s, up_r2s = [], [], [], []

            for ncyc, time_t2, int_val, int_err in sorted(ints):
                if ncyc > 0:
                    r2_val = -sc.log(int_val / int_ref) / time_t2

                    int_val_dist = sc.random.normal(int_val, int_err, 1000.0)
                    r2_val_dist = -sc.log(int_val_dist / int_ref) / time_t2

                    r2_err_down = abs(stats.scoreatpercentile(r2_val_dist, 15.9) - r2_val)
                    r2_err_up = abs(stats.scoreatpercentile(r2_val_dist, 84.1) - r2_val)

                    exp_nus.append(ncyc / time_t2)
                    exp_r2s.append(r2_val)
                    down_r2s.append(r2_err_down)
                    up_r2s.append(r2_err_up)

            cal_nus, cal_r2s = [], []

            for ncyc, time_t2, int_cal in sorted(cal_ints[resonance_id]):
                stdout.flush()
                stderr.flush()

                if ncyc > 0:
                    cal_nus.append(ncyc / time_t2)
                    stdout.flush()
                    stderr.flush()

                    cal_r2s.append(-sc.log(int_cal / int_ref) / time_t2)
                    stdout.flush()
                    stderr.flush()

            fig = plt.figure(linewidth=linewidth)
            ax = fig.add_subplot(111)

            ax.plot(cal_nus, cal_r2s, '-', color='0.5', linewidth=linewidth)

            ax.errorbar(exp_nus, exp_r2s, [down_r2s, up_r2s],
                        fmt='ro',
                        linewidth=linewidth,
                        markerfacecolor='w',
                        markeredgewidth=linewidth,
                        markeredgecolor='r',
                        barsabove=False)

            ax.tick_params(length=2, top=False, right=False)

            xmin, xmax = min(cal_nus), max(cal_nus)
            xmargin = (xmax - xmin) * 0.05
            ax.set_xlim(max([0, xmin - xmargin]), xmax + xmargin)

            ymin = min([min(sc.array(exp_r2s) - sc.array(down_r2s)), min(cal_r2s)])
            ymax = max([max(sc.array(exp_r2s) + sc.array(up_r2s)), max(cal_r2s)])
            ymargin = (ymax - ymin) * 0.05
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(6))

            ax.set_xlim([0.0, ax.get_xlim()[1]])

            ax.set_xlabel(r'$\mathregular{\nu_{CPMG}}$' + ' (Hz)')
            ax.set_ylabel(r'$\mathregular{R_{2,eff} \; (s^{-1})}$')

            ax.set_title('{:s}'.format(resonance_id.upper()))

            pdf.savefig()

        pdf.close()
        plt.close()

    return
