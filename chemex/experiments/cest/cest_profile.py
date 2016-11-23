"""The cest_profile module contains the code for handling CEST profiles."""

import copy
from functools import lru_cache

import numpy as np
from scipy import stats

from chemex import constants, peaks
from chemex.experiments import base_profile
from chemex.experiments.cest import plotting


class CESTProfile(base_profile.BaseProfile):
    """CESTProfile class."""

    def __init__(self, profile_name, measurements, exp_details):
        self.profile_name = profile_name
        self.b1_offsets = measurements['b1_offsets']
        self.val = measurements['intensities']
        self.err = measurements['intensities_err']

        self.reference = (self.b1_offsets <= -1.0e+04) | (self.b1_offsets >= 1.0e+04)

        self.on_resonance_filter = base_profile.check_par(
            exp_details, 'on_resonance_filter', float, default=0.0)

        self.h_larmor_frq = base_profile.check_par(exp_details, 'h_larmor_frq', float)
        self.temperature = base_profile.check_par(exp_details, 'temperature', float)
        self.carrier = base_profile.check_par(exp_details, 'carrier', float)
        self.b1_frq = base_profile.check_par(exp_details, 'b1_frq', float)
        self.time_t1 = base_profile.check_par(exp_details, 'time_t1', float)
        self.b1_inh = base_profile.check_par(exp_details, 'b1_inh', float, np.inf)
        self.b1_inh_res = base_profile.check_par(exp_details, 'b1_inh_res', int, 11)
        self.p_total = base_profile.check_par(exp_details, 'p_total', float, required=False)
        self.l_total = base_profile.check_par(exp_details, 'l_total', float, required=False)
        self.model = base_profile.check_par(exp_details, 'model', default='2st.pb_kex')

        if self.b1_inh_res == 1 or self.b1_inh == np.inf:
            self.omega1s = 2.0 * np.pi * np.array([self.b1_frq])
        else:
            self.omega1s = 2.0 * np.pi * (
                np.linspace(-2.0, 2.0, self.b1_inh_res) * self.b1_inh + self.b1_frq)

        self.b1_weights = stats.norm.pdf(np.linspace(-2.0, 2.0, self.b1_inh_res))
        self.b1_weights /= sum(self.b1_weights)

        self.experiment_name = base_profile.check_par(exp_details, 'name')

        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance_i = self.peak.resonances[0]
        self.ppm_i = 2.0 * np.pi * self.h_larmor_frq * constants.xi_ratio[self.resonance_i['atom']]
        if len(self.peak.resonances) > 1:
            self.resonance_s = self.peak.resonances[1]
            self.ppm_s = 2.0 * np.pi * self.h_larmor_frq * constants.xi_ratio[self.resonance_s[
                'atom']]

        self.calculate_unscaled_profile_cached = lru_cache(5)(self.calculate_unscaled_profile)

        self.map_names = {}

    def calculate_unscaled_profile(self, *args, **kwargs):
        """Calculate the unscaled CEST profile."""
        pass

    def calculate_scale(self, cal):
        """Calculate the scaling factor."""
        scale = (sum(cal * self.val / self.err**2) / sum((cal / self.err)**2))

        return scale

    def calculate_profile(self, params):
        """Calculate the CEST profile."""
        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
        }

        values = self.calculate_unscaled_profile_cached(**kwargs)
        scale = self.calculate_scale(values)

        return values * scale

    def calculate_synthetic_profile(self, params=None, b1_offsets=None):
        """Calculate a synthetic CEST profile."""
        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
        }

        values = self.calculate_unscaled_profile_cached(**kwargs)
        scale = self.calculate_scale(values)

        if b1_offsets is not None:
            self.b1_offsets, b1_orig = b1_offsets, self.b1_offsets
            self.reference, ref_orig = [False for _ in b1_offsets], self.reference
            values = self.calculate_unscaled_profile(**kwargs)
            self.b1_offsets = b1_orig
            self.reference = ref_orig

        return values * scale

    def b1_offsets_to_ppm(self, b1_offsets=None):
        """Convert B1 offset from Hz to ppm."""
        if b1_offsets is None:
            b1_offsets = self.b1_offsets

        return 2.0 * np.pi * b1_offsets / self.ppm_i + self.carrier

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether or not the point should be
        considered in the calculation."""
        cs = params[self.map_names['cs_i_a']].value
        nu_offsets = ((cs - self.carrier) * self.ppm_i / (2.0 * np.pi) - self.b1_offsets)

        mask = abs(nu_offsets) > self.on_resonance_filter * 0.5

        self.b1_offsets = self.b1_offsets[mask]
        self.val = self.val[mask]
        self.err = self.err[mask]
        self.reference = self.reference[mask]

    def print_profile(self, params=None):
        """Print the CEST profile."""
        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.val

        iter_vals = list(zip(self.b1_offsets, self.val, self.err, values))

        output.append("[{}]".format(self.profile_name))
        output.append("# {:>12s}   {:>17s} {:>17s} {:>17s}".format(
            "offset (Hz)", "intensity (exp)", "uncertainty", "intensity (calc)"))

        for b1_offset, val, err, cal in iter_vals:
            line = ("  {0:12.2f} = {1:17.8e} {2:17.8e}".format(b1_offset, val, err))

            if params is not None:
                line += "{:17.8e}".format(cal)
            else:
                line += "{:17s}".format("xxx")

            output.append(line)

        output.append("")
        output.append("")

        return "\n".join(output).upper()

    def make_bs_profile(self):
        """Make a CEST profile for boostrap analysis."""
        indexes = np.array(range(len(self.val)))
        pool1 = indexes[self.reference]
        pool2 = indexes[np.logical_not(self.reference)]

        bs_indexes = []
        if pool1.size:
            bs_indexes.extend(np.random.choice(pool1, len(pool1)))
        bs_indexes.extend(np.random.choice(pool2, len(pool2)))

        bs_indexes = sorted(bs_indexes)

        profile = copy.deepcopy(self)
        profile.b1_offsets = profile.b1_offsets[bs_indexes]
        profile.val = profile.val[bs_indexes]
        profile.err = profile.err[bs_indexes]

        profile.calculate_unscaled_profile_cached = lru_cache(5)(profile.calculate_unscaled_profile)

        return profile
