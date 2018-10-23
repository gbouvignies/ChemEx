"""The cest_profile module contains the code for handling CEST profiles."""
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import linalg
from scipy import signal
from scipy import stats

from chemex.experiments.base import plotting as pl
from chemex.experiments.base.base_profile import BaseProfile
from chemex.experiments.cest import plotting as plc
from chemex.spindynamics import constants

_EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t1": {"type": float},
    "b1_frq": {"type": float},
    "b1_inh": {"default": np.inf, "type": float},
    "b1_inh_res": {"default": 11, "type": int},
    "cn_label": {"default": False, "type": bool},
    "filter_offsets": {"default": 0.0, "type": float},
    "filter_bandwidths": {"default": 0.0, "type": float},
}


class ProfileCEST(BaseProfile):
    """CESTProfile class."""

    EXP_DETAILS = dict(**BaseProfile.EXP_DETAILS, **_EXP_DETAILS)
    DTYPE = [("offsets", "f8"), ("intensity", "f8"), ("error", "f8")]

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.carrier = self.exp_details["carrier"]
        self.time_t1 = self.exp_details["time_t1"]

        self.liouv.w1_i = 2 * np.pi * self.exp_details["b1_frq"]
        self.liouv.w1_i_inh = self.exp_details["b1_inh"]
        self.liouv.w1_i_inh_res = self.exp_details["b1_inh_res"]

        self.carriers_i = self.offsets_to_ppm()

        # Set the smaple labeling
        if self.exp_details["cn_label"]:
            symbol, nucleus = self.peak.symbols["i"], self.peak.nuclei["i"]
            j_values, j_weights = constants.get_multiplet(symbol, nucleus)
            self.liouv.j_eff_i = j_values
            self.liouv.j_eff_i_weights = j_weights

        self.plot_data = plc.plot_data
        self.dephasing = self.exp_details["b1_inh"] == np.inf

    @property
    def reference(self):
        return np.abs(self.data["offsets"]) >= 1.0e04

    def offsets_to_ppm(self, b1_offsets=None):
        """Convert B1 offset from Hz to ppm."""

        if b1_offsets is None:
            b1_offsets = self.data["offsets"]

        return 2.0 * np.pi * b1_offsets / self.ppms_i + self.carrier

    def print_profile(self, params=None):
        """Print the CEST profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.data["intensity"]

        iter_vals = list(
            zip(
                self.data["offsets"], self.data["intensity"], self.data["error"], values
            )
        )

        output.append(f"[{self.name}]")
        output.append(
            "# {'offset (Hz)':>12s}   {'intensity (exp)':>17s} "
            "{'uncertainty':>17s} {'intensity (calc)':>17s}"
        )

        for b1_offset, val, err, cal in iter_vals:
            line = f"  {b1_offset:12.2f} = {val:17.8e} {err:17.8e}"

            if params is not None:
                line += f"{cal:17.8e}"
            else:
                line += f"{'xxx':17s}"

            output.append(line)

        output.append("")
        output.append("")

        return "\n".join(output).upper()

    def estimate_noise(self):
        """Estimate the uncertainty in the CEST profile.

        Adapted from:
        https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise

        """

        data_sorted = np.sort(self.data[~self.reference], order="offsets")
        values = data_sorted["intensity"]
        size = len(values)

        fda = [
            [1, -1],
            [1, -2, 1],
            [1, -3, 3, -1],
            [1, -4, 6, -4, 1],
            [1, -5, 10, -10, 5, -1],
            [1, -6, 15, -20, 15, -6, 1],
        ]

        fda = [np.array(a_fda) / linalg.norm(a_fda) for a_fda in fda]

        percents = np.array([0.05] + list(np.arange(0.1, 0.40, 0.025)))
        percent_points = stats.norm.ppf(1.0 - percents)

        sigma_est = []

        for fdai in fda:
            noisedata = signal.convolve(values, fdai, mode="valid")
            ntrim = len(noisedata)

            if ntrim >= 2:
                noisedata.sort()

                xaxis = 0.5 + np.arange(1, ntrim + 1)
                xaxis /= ntrim + 0.5

                sigmas = []

                f_interpalated = interpolate.interp1d(xaxis, noisedata, "linear")

                for a_perc, a_z in zip(percents, percent_points):
                    try:
                        val = (
                            f_interpalated(1.0 - a_perc) - f_interpalated(a_perc)
                        ) / (2.0 * a_z)
                        sigmas.append(val)
                    except ValueError:
                        pass

                sigma_est.append(np.median(sigmas))

        noisevar = np.median(sigma_est) ** 2
        noisevar /= 1.0 + 15.0 * (size + 1.225) ** -1.245

        return np.sqrt(noisevar)

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether or not the point should be
        considered in the calculation."""

        cs_a = params[self.map_names["cs_i_a"]].value

        filter_offsets = np.asarray(self.exp_details["filter_offsets"]).reshape(-1)
        filter_bandwidths = np.asarray(self.exp_details["filter_bandwidths"]).reshape(
            -1
        )

        for offset, bandwidth in zip(filter_offsets, filter_bandwidths):
            nu_offsets = (
                (cs_a - self.carrier) * self.ppms_i / (2.0 * np.pi)
                - self.data["offsets"]
                + offset
            )

            self.mask = np.logical_and(self.mask, abs(nu_offsets) > bandwidth * 0.5)

    def get_plot_fig(self, params):

        profile_exp, profile_fit = self._get_plot_data(params)

        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 4]})

        fig.align_labels()
        fig.suptitle(f"{self.name.upper()}")

        ########################

        cs_values = [
            params[self.map_names[f"cs_i_{state}"]]
            for state in "abcd"
            if f"cs_i_{state}" in self.map_names
        ]

        lstyles = ("-", "--", "-.", ":")

        for ax_ in (ax1, ax2):
            for a_cs, lstyle in zip(cs_values, lstyles):
                ax_.axvline(
                    a_cs.value,
                    color=pl.PALETTE["Grey"]["400"],
                    linestyle=lstyle,
                    linewidth=0.75,
                )
            ax_.axhline(0, color="k", linewidth=0.5, zorder=0)

        ########################

        xmin, xmax = pl.set_lim(profile_fit["ppm"], 0.02)
        intensities = np.concatenate(
            (profile_exp["intensity"], profile_fit["intensity"])
        )
        ymin_ax2, ymax_ax2 = pl.set_lim(intensities, 0.10)

        ax2.set_xlim(xmin, xmax)
        ax2.set_ylim(ymin_ax2, ymax_ax2)

        ax2.invert_xaxis()

        ax2.set_xlabel(r"$B_1$ position (ppm)")
        ax2.set_ylabel(r"$I/I_0$")

        ndata_used = profile_exp[profile_exp["mask"]]
        ndata_unused = profile_exp[~profile_exp["mask"]]

        ax2.plot(
            profile_fit["ppm"],
            profile_fit["intensity"],
            linestyle="-",
            color=pl.PALETTE["Red"]["300"],
        )

        ax2.errorbar(
            ndata_used["ppm"],
            ndata_used["intensity"],
            ndata_used["error"],
            fmt=".",
            color=pl.PALETTE["Red"]["700"],
            zorder=3,
        )

        ax2.errorbar(
            ndata_unused["ppm"],
            ndata_unused["intensity"],
            ndata_unused["error"],
            fmt=".",
            color=pl.PALETTE["Red"]["100"],
            zorder=3,
        )

        ########################

        deltas = profile_exp["intensity"] - profile_exp["intensity_calc"]
        sigma = _sigma_estimator(deltas)

        deltas_used = deltas[profile_exp["mask"]]
        deltas_unused = deltas[~profile_exp["mask"]]

        ymin_ax1, ymax_ax1 = pl.set_lim(deltas, 0.1)
        ymin_ax1 = min([-3 * sigma, ymin_ax1 - max(profile_exp["error"])])
        ymax_ax1 = max([+3 * sigma, ymax_ax1 + max(profile_exp["error"])])

        ax1.set_ylabel("Residuals")
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin_ax1, ymax_ax1)
        ax1.invert_xaxis()
        ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y", useMathText=True)

        ax1.fill(
            (xmin, xmin, xmax, xmax),
            1.0 * sigma * np.asarray([-1.0, 1.0, 1.0, -1.0]),
            fc=(0, 0, 0, 0.12),
            ec="none",
        )

        ax1.fill(
            (xmin, xmin, xmax, xmax),
            2.0 * sigma * np.asarray([-1.0, 1.0, 1.0, -1.0]),
            fc=(0, 0, 0, 0.12),
            ec="none",
        )

        ax1.errorbar(
            ndata_used["ppm"],
            deltas_used,
            ndata_used["error"],
            fmt=".",
            color=pl.PALETTE["Red"]["700"],
            zorder=3,
        )

        ax1.errorbar(
            ndata_unused["ppm"],
            deltas_unused,
            ndata_unused["error"],
            fmt=".",
            color=pl.PALETTE["Red"]["100"],
            zorder=3,
        )

        #######################

        return fig

    def plot_data_string(self, params):
        """Write the fitted CEST profile."""

        profile_exp, profile_fit = self._get_plot_data(params)

        ndata_ = [
            f"[{self.name.upper()}]",
            f"# {'OFFSET':>17s}   {'INTENSITY':>17s} {'UNCERTAINTY':>17s}",
        ]
        ndata_.extend(
            [
                f"  {values['ppm']:17.8e} "
                f"= {values['intensity']:17.8e} {values['error']:17.8e}"
                for values in profile_exp
            ]
        )
        ndata_.append("")
        ndata_str = "\n".join(ndata_)

        nfit_ = [f"[{self.name.upper()}]", f"# {'OFFSET':>17s}   {'INTENSITY':>17s}"]
        nfit_.extend(
            [
                f"  {values['ppm']:17.8e} = {values['intensity']:17.8e}"
                for values in profile_fit
            ]
        )
        nfit_.append("")
        nfit_str = "\n".join(nfit_)

        return ndata_str, nfit_str

    def _get_plot_data(self, params):

        ndata, scale = self.normalize_profile(params)
        ndata = ndata[~self.reference]

        dt_exp = [
            ("ppm", "f8"),
            ("intensity", "f8"),
            ("intensity_calc", "f8"),
            ("error", "f8"),
            ("mask", bool),
        ]
        dt_fit = [("ppm", "f8"), ("intensity", "f8")]

        profile_exp = np.zeros(len(ndata), dtype=dt_exp)
        profile_fit = np.zeros(500, dtype=dt_fit)

        profile_exp["ppm"] = self.offsets_to_ppm()[~self.reference]
        profile_exp[["intensity", "intensity_calc", "error"]] = ndata[
            ["intensity", "intensity_calc", "error"]
        ]
        profile_exp["mask"] = self.mask[~self.reference]

        o_min, o_max = pl.set_lim(ndata["offsets"], 0.02)
        offsets = np.linspace(o_min, o_max, 500)
        ppms_fit = self.offsets_to_ppm(offsets)
        mag_fit = self.calculate_profile(params, offsets=offsets) * scale

        profile_fit["ppm"] = ppms_fit
        profile_fit["intensity"] = mag_fit

        return profile_exp, profile_fit


def _sigma_estimator(values):
    """Estimates standard deviation using median to exclude outliers.

    Up to 50% can be bad.

    Ref :
    Rousseeuw, Peter & Croux, Christophe. (1993). Alternatives to Median Absolute
    Deviation. Journal of the American Statistical Association. 88. 1273 - 1283.
    10.1080/01621459.1993.10476408.

    """
    _values = values.reshape(1, -1)
    return 1.1926 * np.median(np.median(abs(_values - _values.T), axis=0))
