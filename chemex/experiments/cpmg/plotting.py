import os

import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

from chemex.parsing import parse_assignment


dark_gray = '0.13'
red500 = '#F44336'
red200 = '#EF9A9A'


def set_lim(values, scale):
    """Provides a range that contains all the value and adds a margin."""

    v_min, v_max = min(values), max(values)
    margin = (v_max - v_min) * scale
    v_min, v_max = v_min - margin, v_max + margin

    return v_min, v_max


def group_data(data):
    """Groups the data resonance specifically"""

    data_grouped = dict()

    for data_pt in data:
        resonance_id = data_pt.par['resonance_id']

        assignment = parse_assignment(resonance_id)
        index = int(assignment[0][0])

        data_grouped.setdefault((index, resonance_id), []).append(data_pt)

    return data_grouped


def compute_profiles(data_grouped):
    profiles = {}
    r2_min = +1e16
    r2_max = -1e16

    for (index, resonance_id), profile in data_grouped.items():

        mag_ref = sp.mean(
            [data_pt.val for data_pt in profile if data_pt.par['ncyc'] == 0]
        )

        r2_profile = []

        for data_pt in profile:

            ncyc = data_pt.par['ncyc']
            time_t2 = data_pt.par['time_t2']

            frq = ncyc / time_t2

            if frq:
                mag_cal = data_pt.cal
                mag_exp = data_pt.val
                mag_err = data_pt.err
                mag_ens = sp.random.normal(mag_exp, mag_err, 10000)

                r2_cal = -sp.log(mag_cal / mag_ref) / time_t2
                r2_exp = -sp.log(mag_exp / mag_ref) / time_t2
                r2_ens = -sp.log(mag_ens / mag_ref) / time_t2
                r2_err = abs(sp.percentile(r2_ens, [15.9, 84.1]) - r2_exp)
                r2_erd, r2_eru = r2_err

                r2_profile.append([frq, r2_cal, r2_exp, r2_erd, r2_eru])

                r2_min = min(r2_min, r2_cal, r2_exp - r2_erd)
                r2_max = max(r2_max, r2_cal, r2_exp + r2_eru)

        r2_profile = zip(*sorted(r2_profile))
        profiles.setdefault((index, resonance_id), []).append(r2_profile)

    return profiles, r2_min, r2_max


def write_profile(id, r2_profile, file_txt):
    for frq, r2_cal, r2_exp, r2_erd, r2_eru in zip(*(r2_profile[0])):
        file_txt.write(
            "{:10s} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                id.upper(), frq, r2_cal, r2_exp, r2_erd, r2_eru
            )
        )


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    """Plot dispersion profiles and write a multi-page pdf file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.par['experiment_name']
        datasets.setdefault(experiment_name, list()).append(data_point)

    for experiment_name, dataset in datasets.items():

        # ##### Matplotlib ######

        name_pdf = ''.join([experiment_name, '.pdf'])
        name_pdf = os.path.join(output_dir, name_pdf)

        name_txt = ''.join([experiment_name, '.fit'])
        name_txt = os.path.join(output_dir, name_txt)

        print("  * {} [.fit]".format(name_pdf))

        # #######################

        data_grouped = group_data(dataset)
        profiles, r2_min, r2_max = compute_profiles(data_grouped)
        ymin, ymax = set_lim([r2_min, r2_max], 0.10)

        with PdfPages(name_pdf) as file_pdf, open(name_txt, 'w') as file_txt:

            for (_index, id), profile in sorted(profiles.items()):
                write_profile(id, profile, file_txt)

                ###### Matplotlib ######

                fig = plt.figure(1, frameon=True)
                ax = fig.add_subplot(111)

                ax.axhline(0, color='black', alpha=0.87)

                ########################

                frq, r2_cal, r2_exp, r2_erd, r2_eru = profile[0]

                ax.plot(
                    frq,
                    r2_cal,
                    linestyle='-',
                    color=red200,
                    zorder=2,
                )

                ax.errorbar(
                    frq,
                    r2_exp,
                    yerr=[r2_erd, r2_eru],
                    fmt='o',
                    color=red500,
                    zorder=3,
                )

                xmin, xmax = set_lim(frq, 0.10)

                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)

                ax.xaxis.set_major_locator(MaxNLocator(6))
                ax.yaxis.set_major_locator(MaxNLocator(6))

                ax.set_xlabel(r'$\mathregular{\nu_{CPMG} \ (Hz)}$')
                ax.set_ylabel(
                    r'$\mathregular{R_{2,eff} \ (s^{-1})}$')

                ax.set_title('{:s}'.format(id.upper()))

                fig.tight_layout()

                ########################

                file_pdf.savefig()
                plt.close()

                ########################

    return
