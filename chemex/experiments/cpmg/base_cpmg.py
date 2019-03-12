"""The cpmg_profile module contains the code for handling CPMG profiles."""
import numpy as np
from matplotlib import pyplot as plt

from chemex.experiments.base import plotting as pl
from chemex.experiments.base.base_profile import BaseProfile
from chemex.experiments.cpmg import plotting

_EXP_DETAILS1 = {"time_t2": {"type": float}}

_EXP_DETAILS2 = {
    "carrier": {"type": float},
    "pw90": {"type": float},
    "time_equil": {"default": 0.0, "type": float},
}

np.random.seed(seed=0)


class ProfileCPMG1(BaseProfile):
    """CPMGProfile class."""

    RANDN = np.random.randn(10000, 1)
    EXP_DETAILS = dict(**BaseProfile.EXP_DETAILS, **_EXP_DETAILS1)
    DTYPE = [("ncycs", "i4"), ("intensity", "f8"), ("error", "f8")]

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.time_t2 = self.exp_details["time_t2"]

        # Set the delays
        ncycs = self.data["ncycs"][~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs)))
        self.tau_cps[-1.0] = 0.5 * self.time_t2
        self.delays = list(self.tau_cps.values())

        self.plot_data = plotting.plot_data

    @property
    def reference(self):
        return self.data["ncycs"] == 0

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""

        if ncycs is None:
            ncycs = self.data["ncycs"]

        # ncyc = -1 correspond to a single spin echo case
        ncycs[ncycs == -1] = 0.5

        return ncycs / self.time_t2

    def print_profile(self, params=None):
        """Print the CPMG profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.data["intensity"]

        iter_vals = list(
            zip(self.data["ncycs"], self.data["intensity"], self.data["error"], values)
        )

        output.append(f"[{self.name}]")
        output.append(
            "# {:>5s}   {:>17s} {:>17s} {:>17s}".format(
                "ncyc", "intensity (exp)", "uncertainty", "intensity (calc)"
            )
        )

        for ncyc, val, err, cal in iter_vals:
            line = f"  {ncyc:5.0f} = {val:17.8e} {err:17.8e}"

            if params is not None:
                line += f"{cal:17.8e}"
            else:
                line += "{:17s}".format("xxx")

            output.append(line)

        output.append("")
        output.append("")

        return "\n".join(output).upper()

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether the point should be
        considered in the calculation or not"""
        pass

    def get_variance_from_duplicates(self):
        """Estimate the variance of duplicate points"""

        groups = {}

        for ncyc, intensity in self.data[["ncycs", "intensity"]]:
            groups.setdefault(ncyc, []).append(intensity)

        variances = []
        for values in groups.values():
            n_duplicates = len(values)
            if n_duplicates > 1:
                variances.append(np.var(values, ddof=1))

        if groups:
            variance = np.mean(variances)
        else:
            variance = np.mean(self.data["error"])

        return variance

    def get_plot_fig(self, params):

        profile_exp, profile_fit = self._get_plot_data(params)

        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 4]})

        fig.align_labels()
        fig.suptitle(f"{self.name.upper()}")

        for ax_ in (ax1, ax2):
            ax_.axhline(0, color="k", linewidth=0.5)

        ########################

        ax2.plot(
            profile_fit["nu_cpmg"],
            profile_fit["r2"],
            linestyle="-",
            color=pl.PALETTE["Red"]["300"],
        )

        errors = profile_exp["error"].copy()
        errors[errors == np.inf] = 1e16

        ax2.errorbar(
            profile_exp["nu_cpmg"],
            profile_exp["r2"],
            yerr=abs(errors).T,
            fmt=".",
            color=pl.PALETTE["Red"]["700"],
            zorder=3,
        )

        xmin = 0.0
        xmax = max(profile_exp["nu_cpmg"]) + min(profile_exp["nu_cpmg"])

        r2s = np.concatenate(
            (profile_fit["r2"], np.ravel(profile_exp["r2"] + profile_exp["error"].T))
        )
        ymin_ax2, ymax_ax2 = pl.set_lim(r2s, 0.1)

        ax2.set_xlim(xmin, xmax)
        ax2.set_ylim(ymin_ax2, ymax_ax2)

        ax2.set_xlabel(r"$\nu_\mathregular{CPMG}$ (Hz)")
        ax2.set_ylabel(r"$R_{2,\mathregular{eff}}$ (s$^{-1}$)")

        ########################

        deltas = profile_exp["r2"] - profile_fit["r2"]

        ax1.errorbar(
            profile_exp["nu_cpmg"],
            deltas,
            yerr=abs(profile_exp["error"]).T,
            fmt=".",
            color=pl.PALETTE["Red"]["700"],
            zorder=3,
        )

        deltas_ = np.ravel(deltas + profile_exp["error"].T)
        ymin_ax1, ymax_ax1 = pl.set_lim(deltas_, 0.1)

        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin_ax1, ymax_ax1)

        ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y", useMathText=True)

        ax1.set_ylabel(r"Residuals (s$^{-1}$)")

        #######################

        return fig

    def plot_data_string(self, params):
        """Write the fitted CPMG profile."""

        p_exp, p_fit = self._get_plot_data(params)

        exp_str = [f"[{self.name.upper()}]"]
        exp_str += [
            f"# {'NU_CPMG':>17s}   {'R2':>17s} "
            f"{'UNCERTAINTY_DOWN':>17s} {'UNCERTAINTY_UP':>17s}"
        ]
        exp_str += [
            f"  {nu_cpmg:17.8e} = {r2:17.8e} {erd:17.8e} {eru:17.8e}"
            for nu_cpmg, r2, (erd, eru) in p_exp
        ]
        exp_str.append("")
        exp_str = "\n".join(exp_str)

        fit_str = [f"[{self.name.upper()}]"]
        fit_str += [f"# {'NU_CPMG':>17s}   {'R2':>17s}"]
        fit_str += [f"  {nu_cpmg:17.8e} = {r2:17.8e}" for nu_cpmg, r2 in p_fit]
        fit_str += ""
        fit_str = "\n".join(fit_str)

        return exp_str, fit_str

    def _get_plot_data(self, params):

        nu_cpmg = self.ncycs_to_nu_cpmgs()[~self.reference]
        ndata, _ = self.normalize_profile(params)
        ndata = ndata[~self.reference]

        dt_exp = [("nu_cpmg", "f8"), ("r2", "f8"), ("error", "f8", (2,))]
        dt_fit = [("nu_cpmg", "f8"), ("r2", "f8")]

        profile_exp = np.zeros(len(ndata), dtype=dt_exp)
        profile_fit = np.zeros(len(ndata), dtype=dt_fit)

        profile_exp["nu_cpmg"] = profile_fit["nu_cpmg"] = nu_cpmg

        profile_exp["r2"] = -np.log(ndata["intensity"]) / self.time_t2
        profile_fit["r2"] = -np.log(ndata["intensity_calc"]) / self.time_t2

        # MC simulation to propagate the error
        mag_ens = ndata["error"] * self.RANDN + ndata["intensity"]
        r2_ens = np.empty_like(mag_ens)
        neg_vals = mag_ens <= 0.0
        r2_ens[~neg_vals] = -np.log(mag_ens[~neg_vals]) / self.time_t2
        r2_ens[neg_vals] = np.inf
        r2_ens -= profile_exp["r2"]
        profile_exp["error"] = np.percentile(r2_ens, [15.9, 84.1], axis=0).transpose()

        profile_exp.sort(order="nu_cpmg")
        profile_fit.sort(order="nu_cpmg")

        return profile_exp, profile_fit


class ProfileCPMG2(ProfileCPMG1):
    """CPMGProfile class."""

    EXP_DETAILS = dict(**ProfileCPMG1.EXP_DETAILS, **_EXP_DETAILS2)

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.pw90 = self.exp_details["pw90"]
        self.time_eq = self.exp_details["time_equil"]

        self.liouv.carrier_i = self.exp_details["carrier"]
        self.liouv.w1_i = np.pi / (2.0 * self.pw90)

        # Set the delays in the experiments
        ncycs = self.data["ncycs"][~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs) - self.pw90))
        self.tau_cps[-1.0] = 0.5 * self.time_t2
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.delays = [self.t_neg, self.time_eq] + list(self.tau_cps.values())

        self.plot_data = plotting.plot_data
