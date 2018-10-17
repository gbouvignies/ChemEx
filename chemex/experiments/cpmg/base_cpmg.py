"""The cpmg_profile module contains the code for handling CPMG profiles."""
import copy
import math
from functools import lru_cache

import numpy as np

from chemex.experiments.base.base_profile import BaseProfile
from chemex.experiments.cpmg import plotting


class ProfileCPMG(BaseProfile):
    """CPMGProfile class."""

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.plot_data = plotting.plot_data

    @property
    def reference(self):
        return self.data["ncycs"] == 0

    def print_profile(self, params=None):
        """Print the CPMG profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.data["intensities"]

        iter_vals = list(
            zip(
                self.data["ncycs"],
                self.data["intensities"],
                self.data["errors"],
                values,
            )
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
        considered in the calculation or not."""

        self.mask = np.ones_like(self.data["intensities"], dtype=np.bool)

    def make_bs_profile(self):
        """Make a CPMG profile for boostrap analysis."""
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
        profile.calculate_unscaled_profile = lru_cache(256)(
            profile._calculate_unscaled_profile
        )

        return profile

    def estimate_noise(self):
        """Estimate the uncertainty in CPMG data points (i.e., R2eff)."""

        intensity_dict = {}

        for ncyc, intensity in zip(self.data["ncycs"], self.data["intensities"]):
            intensity_dict.setdefault(ncyc, []).append(intensity)

        std_list = []
        for duplicates in intensity_dict.values():
            n_duplicates = len(duplicates)
            if n_duplicates > 1:
                std_list.append(np.std(duplicates, ddof=1) * _correction(n_duplicates))

        if std_list:
            error = np.mean(std_list)
        else:
            error = np.mean(self.data["errors"])

        return error


def _correction(n):
    """Calculate correction factor for noise estimate."""
    k = n // 2

    if n == 2 * k:
        factor = (
            math.sqrt(math.pi * (2 * k - 1) / 2.0)
            * math.factorial(2 * k - 2)
            / ((2 ** (2 * k - 2)) * (math.factorial(k - 1) ** 2))
        )
    else:
        factor = (
            math.sqrt(k / math.pi)
            * (2 ** (2 * k - 1))
            * (math.factorial(k - 1) ** 2)
            / math.factorial(2 * k - 1)
        )

    return factor
