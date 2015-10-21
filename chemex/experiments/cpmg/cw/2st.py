"""Pure in-phase CPMG (2-state)

Analyzes chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]

Notes
-----
Off resonance effects are taken into account.

The calculation is designed specifically to analyze the experiment found in
the reference:

Journal of Physical Chemistry B (2008), 112, 5898-5904
"""
from __future__ import absolute_import

import functools
from functools import reduce

import lmfit
import numpy as np
from scipy import linalg

from chemex import constants, parameters, peaks
from chemex.bases import iph_2st, util as bases_util
from chemex.experiments import base_profile
from chemex.experiments.cpmg import plotting
from six.moves import zip

try:
    from functools import lru_cache
except ImportError:
    from chemex.lru_cache import lru_cache

reduce = functools.reduce
dot = np.dot
matrix_power = np.linalg.matrix_power
compute_liouvillian = iph_2st.compute_liouvillian
check_par = base_profile.check_par
compute_propagator_from_time_series = bases_util.compute_propagators_from_time_series
ParameterName = parameters.ParameterName

two_pi = 2.0 * np.pi

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'carrier': float,
    'pw': float,
    'time_t2': float,
    'time_equil': float,
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.ncycs = measurements['ncycs']
        self.val = measurements['intensities']
        self.err = measurements['intensities_err']

        self.h_larmor_frq = check_par(exp_details, 'h_larmor_frq', float)
        self.temperature = check_par(exp_details, 'temperature', float)
        self.carrier = check_par(exp_details, 'carrier', float)
        self.pw = check_par(exp_details, 'pw', float)
        self.time_t2 = check_par(exp_details, 'time_t2', float)
        self.time_equil = check_par(exp_details, 'time_equil', float)

        self.experiment_name = check_par(exp_details, 'experiment_name')
        self.tau_cp_list = np.array([self.time_t2 / (4.0 * ncyc) - self.pw if ncyc else 0.0 for ncyc in self.ncycs])

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance = self.peak.resonances[0]

        self.ppm_to_rads = (
            self.h_larmor_frq * two_pi *
            constants.xi_ratio[self.resonance['atom']]
        )

        kwargs1 = {'temperature': self.temperature}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance['name']}
        kwargs3 = {'temperature': self.temperature, 'nuclei': self.resonance['name'], 'h_larmor_frq': self.h_larmor_frq}

        self.map_names = {
            'pb': ParameterName('pb', **kwargs1).to_full_name(),
            'kex_ab': ParameterName('kex_ab', **kwargs1).to_full_name(),
            'cs_i_a': ParameterName('cs_a', **kwargs2).to_full_name(),
            'dw_i_ab': ParameterName('dw_ab', **kwargs2).to_full_name(),
            'lambda_i_a': ParameterName('lambda_a', **kwargs3).to_full_name(),
            'dlambda_i_ab': ParameterName('dlambda_ab', **kwargs3).to_full_name(),
            'rho_i_a': ParameterName('rho_a', **kwargs3).to_full_name(),
        }

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (self.map_names['cs_i_a'], 0.0, False, None, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
            (self.map_names['lambda_i_a'], 10.0, True, 0.0, None, None),
            (self.map_names['dlambda_i_ab'], 0.0, False, None, None, None),
            (self.map_names['rho_i_a'], 1.0, False, 0.0, None, None),
        )

        return parameters

    def _calculate_profile(self, pb, kex_ab, dw_i_ab, rho_i_a, lambda_i_a, dlambda_i_ab, cs_i_a):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex_ab : float
            Exchange rates between states A and B in /s.
        dw_i_ab : float
            Chemical shift difference between states A and B in rad/s.
        rho_i_a : float
            Longitudinal relaxation rate of states A in /s.
        lambda_i_a : float
            Transverse relaxation rate of state A in /s.
        dlambda_i_ab : float
            Transverse relaxation rate difference between states A and B in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block
        """

        domega_i_ab = dw_i_ab * self.ppm_to_rads
        omega_i_a = (cs_i_a - self.carrier) * self.ppm_to_rads

        pa = 1.0 - pb
        mag_eq = np.array([[0.0, 0.0, pa,
                            0.0, 0.0, pb]]).T

        w1 = 2.0 * np.pi / (4.0 * self.pw)

        # Calculation of the different liouvillians

        kwargs = {
            'pb': pb,
            'kex_ab': kex_ab,
            'rho_i_a': rho_i_a, 'lambda_i_a': lambda_i_a, 'omega_i_a': omega_i_a,
            'rho_i_b': rho_i_a, 'lambda_i_b': lambda_i_a + dlambda_i_ab, 'omega_i_b': omega_i_a + domega_i_ab
        }

        l_free = iph_2st.compute_liouvillian(**kwargs)
        l_pw1x = iph_2st.compute_liouvillian(omega1x_i=+w1, **kwargs)
        l_pw1y = iph_2st.compute_liouvillian(omega1y_i=+w1, **kwargs)
        l_mw1x = iph_2st.compute_liouvillian(omega1x_i=-w1, **kwargs)

        # Calculation of all the needed propagators

        p_90px = linalg.expm(l_pw1x * self.pw)
        p_90py = linalg.expm(l_pw1y * self.pw)
        p_90mx = linalg.expm(l_mw1x * self.pw)
        p_180pmx = 0.5 * (matrix_power(p_90px, 2) + matrix_power(p_90mx, 2))  # +/- phase cycling
        p_180py = matrix_power(p_90py, 2)

        t_neg = -2.0 * self.pw / np.pi
        time_series = [t_neg, self.time_equil]
        time_series.extend(self.tau_cp_list)

        p_free_list = compute_propagator_from_time_series(l_free, time_series)

        p_equil = p_free_list[self.time_equil]
        p_neg = p_free_list[t_neg]

        # Simulate the CPMG block as function of ncyc

        profile = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):

            if ncyc == 0:

                mag = reduce(dot, [p_equil, p_90px, p_180pmx, p_90px, mag_eq])

            else:

                p_free = p_free_list[tau_cp]
                p_cp = matrix_power(dot(dot(p_free, p_180py), p_free), int(ncyc))
                mag = reduce(dot, [p_equil, p_90px, p_neg, p_cp, p_180pmx, p_cp, p_neg, p_90px, mag_eq])

            profile.append(mag[2, 0])

        return np.asarray(profile)

    def calculate_profile(self, params, ncycs=None):

        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
            }

        values = self._calculate_profile(**kwargs)
        scale = self._calculate_scale(values)

        if ncycs is not None:
            self.ncycs, ncycs_orig = ncycs, self.ncycs
            values = self._calculate_profile.__wrapped__(**kwargs)
            self.ncycs = ncycs_orig

        return values * scale

    def _calculate_scale(self, cal):

        scale = (
            sum(cal * self.val / self.err ** 2) /
            sum((cal / self.err) ** 2)
        )

        return scale

    def calculate_residuals(self, params):
        """Calculates the residual between the experimental and
        back-calculated values.
        """

        values = self.calculate_profile(params)

        return (self.val - values) / self.err

    def ncycs_to_nu_cpmgs(self, ncycs=None):

        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.time_t2

    def filter_points(self, params):
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

        iter_vals = list(zip(self.ncycs, self.val, self.err, values))

        for ncyc, val, err, cal in iter_vals:

            line = (
                "{0.profile_name:10s} "
                "{0.h_larmor_frq:8.1f} "
                "{0.time_t2:8.1e} "
                "{0.temperature:5.1f} "
                "{1:15.8e} "
                "{2:15.8e} "
                "{3:15.8e} "
                    .format(self, ncyc, val, err)
            )

            if params is not None:
                line += "{:15.8e}".format(cal)

            output.append(line)

        output.append("")

        return "\n".join(output).upper()
