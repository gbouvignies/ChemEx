"""HMQC-HSQC shifts

"""

import copy
from functools import lru_cache

import lmfit
import numpy as np

from chemex import constants, parameters, peaks
from chemex.bases import util
from chemex.experiments import base_profile
from chemex.experiments.shift import plotting

attributes_exp = {
    'h_larmor_frq': float,
    'temperature' : float,
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.val = np.array([measurements['shifts']])
        self.err = np.array([measurements['shifts_err']])

        self.h_larmor_frq = base_profile.check_par(exp_details, 'h_larmor_frq', float)
        self.temperature = base_profile.check_par(exp_details, 'temperature', float)

        self.experiment_name = base_profile.check_par(exp_details, 'name')

        self._calculate_profile_cached = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance_i, self.resonance_s = self.peak.resonances

        self.ppm_i = (self.h_larmor_frq * 2.0 * np.pi * constants.xi_ratio[self.resonance_i['atom']])
        self.ppm_s = (self.h_larmor_frq * 2.0 * np.pi * constants.xi_ratio[self.resonance_s['atom']])

        kwargs1 = {'temperature': self.temperature}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance_i['name']}
        kwargs3 = {'temperature': self.temperature, 'nuclei': self.resonance_s['name']}

        self.map_names = {
            'pb'     : parameters.ParameterName('pb', **kwargs1).to_full_name(),
            'kex_ab' : parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
            'dw_i_ab': parameters.ParameterName('dw_ab', **kwargs2).to_full_name(),
            'dw_s_ab': parameters.ParameterName('dw_ab', **kwargs3).to_full_name(),
        }

        self.default_params = self.create_default_parameters()

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
            (self.map_names['dw_s_ab'], 0.0, False, None, None, None),
        )

        return parameters

    def _calculate_profile(self, pb, kex_ab, dw_i_ab, dw_s_ab):
        """Calculate the intensity in presence of exchange after a CPMG block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_nxy : float
            Transverse relaxation rate of state A in /s.
        dr_nxy : float
            Transverse relaxation rate difference between states A and B in /s.

        Returns
        -------
        out : float
            Intensity after the CEST block
        """

        domega_i_ab = dw_i_ab * self.ppm_i
        domega_s_ab = dw_s_ab * self.ppm_s

        shift_sq = util.calculate_shift_ex_2st(pb, kex_ab, domega_i_ab)[0]
        shift_mq = 0.5 * (util.calculate_shift_ex_2st(pb, kex_ab, domega_i_ab + domega_s_ab)[0] +
                          util.calculate_shift_ex_2st(pb, kex_ab, domega_i_ab - domega_s_ab)[0])

        return np.asarray([(shift_sq - shift_mq) / self.ppm_i])

    def calculate_profile(self, params=None, **kwargs):
        kwargs_profile = {short_name: params[long_name].value for short_name, long_name in self.map_names.items()}
        values = self._calculate_profile_cached(**kwargs_profile)

        return values

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether the point should be considered
        in the calculation or not.

        Returns 'True' if the point should NOT be considered.
        """

        return False

    def print_profile(self, params=None):
        """Print the data point"""

        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.val

        for val, err, cal in zip(self.val, self.err, values):

            line = (
                "{0.profile_name:10s} "
                "{0.h_larmor_frq:8.1f} "
                "{0.temperature:5.1f} "
                "{1:15.8e} "
                "{2:15.8e} "
                    .format(self, val, err)
            )

            if params is not None:
                line += "{:15.8e}".format(cal)

            output.append(line)

        output.append("")

        return "\n".join(output).upper()

    def make_bs_profile(self):

        indexes = np.array(range(len(self.val)))

        bs_indexes = []
        bs_indexes.extend(np.random.choice(indexes, len(indexes)))

        bs_indexes = sorted(bs_indexes)

        profile = copy.deepcopy(self)
        profile.val = profile.val[bs_indexes]
        profile.err = profile.err[bs_indexes]

        profile._calculate_profile_cached = lru_cache(5)(profile.calculate_unscaled_profile)

        return profile
