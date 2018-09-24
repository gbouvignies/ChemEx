"""The cest_profile module contains the code for handling CEST profiles."""

import copy
from functools import lru_cache

import numpy as np

from chemex import peaks
from chemex.experiments import base_profile
from chemex.experiments.cest import plotting


class CESTProfile(base_profile.BaseProfile):
    """CESTProfile class."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, exp_details, model)

        self.plot_data = plotting.plot_data
        self.calculate_unscaled_profile_cached = lru_cache(8)(
            self.calculate_unscaled_profile
        )

        self.peak = peaks.Peak(self.profile_name)

        self.b1_offsets = measurements["b1_offsets"]
        self.val = measurements["intensities"]
        self.err = measurements["intensities_err"]

        self.mask = np.ones_like(self.b1_offsets, dtype=np.bool)
        self.reference = np.abs(self.b1_offsets) >= 1.0e+04

        self.map_names = {}
        self.basis = None

    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and back-calculated
        values."""

        values = self.calculate_profile(params)
        residuals = (self.val - values) / self.err

        return residuals[self.mask]

    def calculate_unscaled_profile(self, b1_offsets=None, **parvals):
        """Calculate the unscaled CEST profile."""
        pass

    def calculate_scale(self, cal):
        """Calculate the scaling factor."""

        norm = sum((cal / self.err) ** 2)

        if norm:
            scale = sum(cal * self.val / self.err ** 2) / norm
        else:
            scale = 0.0

        return scale

    def calculate_profile(self, params=None, b1_offsets=None):
        """Calculate the CEST profile."""
        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
        }

        values = self.calculate_unscaled_profile_cached(**kwargs)
        scale = self.calculate_scale(values)

        if b1_offsets is not None:
            values = self.calculate_unscaled_profile(b1_offsets=b1_offsets, **kwargs)

        return values * scale

    def print_profile(self, params=None):
        """Print the CEST profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.val

        iter_vals = list(zip(self.b1_offsets, self.val, self.err, values))

        output.append("[{}]".format(self.profile_name))
        output.append(
            "# {:>12s}   {:>17s} {:>17s} {:>17s}".format(
                "offset (Hz)", "intensity (exp)", "uncertainty", "intensity (calc)"
            )
        )

        for b1_offset, val, err, cal in iter_vals:
            line = "  {0:12.2f} = {1:17.8e} {2:17.8e}".format(b1_offset, val, err)

            if params is not None:
                line += "{:17.8e}".format(cal)
            else:
                line += "{:17s}".format("xxx")

            output.append(line)

        output.append("")
        output.append("")

        return "\n".join(output).upper()

    def make_mc_profile(self, params):
        """Make a CEST profile for boostrap analysis."""

        profile = copy.copy(self)
        profile.val = (
            self.calculate_profile(params) + np.random.randn(len(self.val)) * self.err
        )

        return profile

    def make_bs_profile(self):
        """Make a CEST profile for boostrap analysis."""
        indexes = np.array(range(len(self.val)))
        pool1 = indexes[self.reference]
        pool2 = indexes[~self.reference]

        bs_indexes = []
        if pool1.size:
            bs_indexes.extend(np.random.choice(pool1, len(pool1)))
        bs_indexes.extend(np.random.choice(pool2, len(pool2)))

        bs_indexes = sorted(bs_indexes)

        profile = copy.copy(self)
        profile.b1_offsets = profile.b1_offsets[bs_indexes]
        profile.val = profile.val[bs_indexes]
        profile.err = profile.err[bs_indexes]

        profile.calculate_unscaled_profile_cached = lru_cache(5)(
            profile.calculate_unscaled_profile
        )

        return profile
