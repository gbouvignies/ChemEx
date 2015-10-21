"""HMQC-HSQC shifts

"""
from __future__ import absolute_import

import lmfit
import numpy as np

from chemex import constants, parameters, peaks
from chemex.bases import util
from chemex.experiments import base_profile
from chemex.experiments.shift import plotting
from six.moves import zip

try:
    from functools import lru_cache
except ImportError:
    from chemex.lru_cache import lru_cache

check_par = base_profile.check_par
ParameterName = parameters.ParameterName

two_pi = 2.0 * np.pi

attributes_exp = {
    'h_larmor_frq_1': float,
    'h_larmor_frq_2': float,
    'temperature': float,
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.val = np.array([measurements['shifts']])
        self.err = np.array([measurements['shifts_err']])

        self.h_larmor_frq_1 = check_par(exp_details, 'h_larmor_frq_1', float)
        self.h_larmor_frq_2 = check_par(exp_details, 'h_larmor_frq_2', float)
        self.temperature = check_par(exp_details, 'temperature', float)

        self.experiment_name = check_par(exp_details, 'experiment_name')

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance = self.peak.resonances[0]

        self.ppm_to_rads_1 = (self.h_larmor_frq_1 * two_pi * constants.xi_ratio[self.resonance['atom']])
        self.ppm_to_rads_2 = (self.h_larmor_frq_2 * two_pi * constants.xi_ratio[self.resonance['atom']])

        kwargs1 = {'temperature': self.temperature}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance['name']}

        self.map_names = {
            'pb': ParameterName('pb', **kwargs1).to_full_name(),
            'kex_ab': ParameterName('kex_ab', **kwargs1).to_full_name(),
            'dw_i_ab': ParameterName('dw_ab', **kwargs2).to_full_name(),
        }

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
        )

        return parameters

    def _calculate_profile(self, pb, kex_ab, dw_i_ab):
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

        domega_i_ab_1 = dw_i_ab * self.ppm_to_rads_1
        domega_i_ab_2 = dw_i_ab * self.ppm_to_rads_2

        shift_sq_1 = util.calculate_shift_ex_2st(pb, kex_ab, domega_i_ab_1)[0] / self.ppm_to_rads_1
        shift_sq_2 = util.calculate_shift_ex_2st(pb, kex_ab, domega_i_ab_2)[0] / self.ppm_to_rads_2

        return np.asarray([shift_sq_2 - shift_sq_1])

    def calculate_profile(self, params=None, **kwargs):
        kwargs_profile = {short_name: params[long_name].value for short_name, long_name in self.map_names.items()}
        values = self._calculate_profile(**kwargs_profile)

        return values

    def calculate_residuals(self, params):
        """Calculates the residual between the experimental and
        back-calculated values.
        """

        values = self.calculate_profile(params)

        return (self.val - values) / self.err

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
                "{0.h_larmor_frq_1:8.1f} "
                "{0.h_larmor_frq_2:8.1f} "
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
