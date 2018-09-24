"""Plot the CPMG profiles."""

import os

import numpy as np
from matplotlib import gridspec as gsp, pyplot as plt, ticker
from matplotlib.backends import backend_pdf

from chemex.experiments import plotting


def compute_profiles(data_grouped, params):
    """Compute the CPMG profiles used for plotting."""
    profiles = {}

    for peak, profile in data_grouped.items():
        mask = profile.ncycs != 0
        mask_ref = np.logical_not(mask)

        mag_ref = np.mean(profile.val[mask_ref])

        nu_cpmg = profile.ncycs_to_nu_cpmgs()[mask]

        mag_cal = profile.calculate_profile(params)[mask]
        mag_exp = profile.val[mask]
        mag_err = profile.err[mask]

        # sort values according to increasing nu_cpmg values
        sorted_items = list(zip(*sorted(zip(nu_cpmg, mag_cal, mag_exp, mag_err))))
        nu_cpmg, mag_cal, mag_exp, mag_err = sorted_items

        mag_ens = mag_exp + mag_err * np.random.randn(1000, len(mag_exp))

        r2_cal = -np.log(mag_cal / mag_ref) / profile.exp_details["time_t2"]
        r2_exp = -np.log(mag_exp / mag_ref) / profile.exp_details["time_t2"]
        r2_ens = -np.log(mag_ens / mag_ref) / profile.exp_details["time_t2"]

        r2_erd, r2_eru = abs(np.percentile(r2_ens, [15.9, 84.1], axis=0) - r2_exp)

        profiles[peak] = nu_cpmg, r2_cal, r2_exp, r2_erd, r2_eru

    return profiles


def write_profile_fit(name, profile, file_txt):
    """Write the fitted CPMG profile."""

    file_txt.write("[{}]\n".format(name.upper()))

    file_txt.write("# {:>17s}   {:>17s}\n".format("NU_CPMG", "R2"))

    for nu_cpmg, r2_cal, _r2_exp, _r2_erd, _r2_eru in sorted(zip(*profile)):
        file_txt.write("  {0:17.8e} = {1:17.8e}\n".format(nu_cpmg, r2_cal))

    file_txt.write("\n")


def write_profile_exp(name, profile, file_txt):
    """Write the experimental CPMG profile."""

    file_txt.write("[{}]\n".format(name.upper()))

    file_txt.write(
        "# {:>17s}   {:>17s} {:>17s} {:>17s}\n".format(
            "NU_CPMG", "R2", "UNCERTAINTY_DOWN", "UNCERTAINTY_UP"
        )
    )

    for nu_cpmg, _r2_cal, r2_exp, r2_erd, r2_eru in sorted(zip(*profile)):
        file_txt.write(
            "  {0:17.8e} = {1:17.8e} {2:17.8e} {3:17.8e}\n".format(
                nu_cpmg, r2_exp, r2_erd, r2_eru
            )
        )

    file_txt.write("\n")


def plot_data(data, params, output_dir="./"):
    """Write experimental and fitted data to a file and plot the CPMG profiles.

    - *.fit: contains the experimental and fitted data
    - *.pdf: contains the plot of experimental and fitted data

    """
    datasets = dict()

    for data_point in data:
        experiment_name = data_point.experiment_name
        datasets.setdefault(experiment_name, []).append(data_point)

    for experiment_name, dataset in datasets.items():
        name_pdf = "".join([experiment_name, ".pdf"])
        name_pdf = os.path.join(output_dir, name_pdf)

        name_exp = "".join([experiment_name, ".exp"])
        name_exp = os.path.join(output_dir, name_exp)

        name_fit = "".join([experiment_name, ".fit"])
        name_fit = os.path.join(output_dir, name_fit)

        print(("  * {} [.fit]".format(name_pdf)))

        data_grouped = plotting.group_data(dataset)
        profiles = compute_profiles(data_grouped, params)

        with backend_pdf.PdfPages(name_pdf) as file_pdf, open(
            name_fit, "w"
        ) as file_fit, open(name_exp, "w") as file_exp:

            for peak in sorted(profiles):
                nu_cpmg, r2_cal, r2_exp, r2_erd, r2_eru = profiles[peak]
                write_profile_exp(peak.assignment, profiles[peak], file_exp)
                write_profile_fit(peak.assignment, profiles[peak], file_fit)

                r2_min = min(min(r2_cal), min(r2_exp - r2_erd))
                r2_max = max(max(r2_cal), max(r2_exp + r2_eru))
                ymin, ymax = plotting.set_lim([r2_min, r2_max], 0.10)

                # Matplotlib #
                gs = gsp.GridSpec(2, 1, height_ratios=[1, 4])

                ax1 = plt.subplot(gs[0])
                ax2 = plt.subplot(gs[1])

                ax1.axhline(0, color=plotting.palette["Black"]["Text"], linewidth=0.5)
                ax2.axhline(0, color=plotting.palette["Black"]["Text"], linewidth=0.5)

                ########################

                ax2.plot(
                    nu_cpmg,
                    r2_cal,
                    linestyle="-",
                    color=plotting.palette["Grey"]["700"],
                    zorder=2,
                )

                ax2.errorbar(
                    nu_cpmg,
                    r2_exp,
                    yerr=(r2_erd, r2_eru),
                    fmt="o",
                    markeredgecolor=plotting.palette["Red"]["500"],
                    ecolor=plotting.palette["Red"]["500"],
                    markerfacecolor="None",
                    zorder=3,
                )

                xmin, xmax = plotting.set_lim(nu_cpmg, 0.10)

                ax2.set_xlim(xmin, xmax)
                ax2.set_ylim(ymin, ymax)

                ax2.xaxis.set_major_locator(ticker.MaxNLocator(6))
                ax2.yaxis.set_major_locator(ticker.MaxNLocator(6))

                ax2.set_xlabel(r"$\nu_{CPMG} \ (Hz)$")
                ax2.set_ylabel(r"$R_{2,eff} \ (s^{-1})$")

                ax1.set_title("{:s}".format(peak.assignment.upper()))

                ax2.yaxis.set_ticks_position("left")
                ax2.xaxis.set_ticks_position("bottom")

                ########################

                deltas = np.asarray(r2_exp) - np.asarray(r2_cal)

                ax1.errorbar(
                    nu_cpmg,
                    deltas,
                    yerr=(r2_erd, r2_eru),
                    fmt="o",
                    markeredgecolor=plotting.palette["Red"]["500"],
                    ecolor=plotting.palette["Red"]["500"],
                    markerfacecolor="None",
                    zorder=100,
                )

                rmin, _rmax = plotting.set_lim(deltas - r2_erd, 0.2)
                _rmin, rmax = plotting.set_lim(deltas + r2_eru, 0.2)

                ax1.set_xlim(xmin, xmax)
                ax1.set_ylim(rmin, rmax)

                ax1.xaxis.set_major_locator(ticker.MaxNLocator(6))
                ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))

                ax1.xaxis.set_major_formatter(ticker.NullFormatter())

                ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")

                ax1.set_title("{:s}".format(peak.assignment.upper()))
                ax1.set_ylabel("Residual " + r"$(s^{-1})$")

                ########################

                for ax in (ax1, ax2):
                    ax.yaxis.set_ticks_position("left")
                    ax.xaxis.set_ticks_position("bottom")

                ########################

                file_pdf.savefig()
                plt.close()

                ########################

    return
