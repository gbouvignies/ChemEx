"""The cest_profile module contains the code for handling CEST profiles."""
import copy

import numpy as np
from scipy import interpolate
from scipy import linalg
from scipy import signal
from scipy import stats

from chemex.experiments import base_intensity
from chemex.experiments.cest import plotting


class ProfileCEST(base_intensity.ProfileIntensity):
    """CESTProfile class."""

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.plot_data = plotting.plot_data
        self.mask = np.ones_like(self.data["offsets"], dtype=np.bool)

    @property
    def reference(self):
        return np.abs(self.data["offsets"]) >= 1.0e04

    def print_profile(self, params=None):
        """Print the CEST profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.data["intensities"]

        iter_vals = list(
            zip(
                self.data["offsets"],
                self.data["intensities"],
                self.data["errors"],
                values,
            )
        )

        output.append(f"[{self.name}]")
        output.append(
            "# {:>12s}   {:>17s} {:>17s} {:>17s}".format(
                "offset (Hz)", "intensity (exp)", "uncertainty", "intensity (calc)"
            )
        )

        for b1_offset, val, err, cal in iter_vals:
            line = f"  {b1_offset:12.2f} = {val:17.8e} {err:17.8e}"

            if params is not None:
                line += f"{cal:17.8e}"
            else:
                line += "{:17s}".format("xxx")

            output.append(line)

        output.append("")
        output.append("")

        return "\n".join(output).upper()

    def make_bs_profile(self):
        """Make a CEST profile for boostrap analysis."""
        indexes = np.array(range(len(self.data["intensities"])))
        pool1 = indexes[self.reference]
        pool2 = indexes[~self.reference]

        bs_indexes = []
        if pool1.size:
            bs_indexes.extend(np.random.choice(pool1, len(pool1)))
        bs_indexes.extend(np.random.choice(pool2, len(pool2)))

        bs_indexes = sorted(bs_indexes)

        profile = copy.copy(self)
        profile.data = self.data[bs_indexes]

        profile._calculate_unscaled_profile.cache_clear()

        return profile

    def estimate_noise(self):
        """Estimate the uncertainty in the CEST profile.

        # TODO: add reference to this method...

        """

        data_sorted = np.sort(self.data[~self.reference], order="offsets")
        values = data_sorted["intensities"]
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
