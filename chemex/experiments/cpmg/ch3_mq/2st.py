"""1H-13C(methyl) - Multiple Quantum CPMG (2-state)

Analyzes HyCx methyl group multiple quantum CPMG measured on site-specific
13CH3-labeled methyl groups in a highly deuterated background.  This is a
simplified basis set, which assumes you are on-resonance for 13C (ie, off-
resonance effects are not taken into account) as described in the reference:

[HxCx(a), HyCx(a), HxCy(a), HyCy(a),
 HxCx(b), HyCx(b), HxCy(b), HyCy(b)]

Note
----
This calculation is designed specifically to analyze data from the experiment
found in the reference and can be run with either small_protein_flag='y' or 'n'.

Lewis Kay experiment: hmqc_CH3_exchange_bigprotein_*00_lek_v2

Journal of the American Chemical Society (2004), 126, 3964-73

"""
from __future__ import absolute_import

import functools
from functools import reduce

import lmfit
import numpy as np

from chemex import constants, parameters, peaks
from chemex.bases import mq_2st, util as bases_util
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
compute_liouvillian = mq_2st.compute_liouvillian
p_180x_i = mq_2st.p_180x_i
p_180y_i = mq_2st.p_180y_i
p_180x_s = mq_2st.p_180x_s
p_180y_s = mq_2st.p_180y_s
check_par = base_profile.check_par
compute_propagator_from_time_series = bases_util.compute_propagators_from_time_series
ParameterName = parameters.ParameterName

two_pi = 2.0 * np.pi
t_zeta = 1.0 / (8.0 * 125.3)

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'time_t2': float,
    'smallflg': str,
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.ncycs = measurements['ncycs']
        self.val = measurements['intensities']
        self.err = measurements['intensities_err']

        self.h_larmor_frq = check_par(exp_details, 'h_larmor_frq', float)
        self.temperature = check_par(exp_details, 'temperature', float)
        self.time_t2 = check_par(exp_details, 'time_t2', float)
        self.small_protein_flg = (check_par(exp_details, 'small_protein_flg', str).lower() == 'y')

        self.experiment_name = check_par(exp_details, 'experiment_name')
        self.tau_cp_list = np.array([self.time_t2 / (4.0 * ncyc) if ncyc else 0.0 for ncyc in self.ncycs])

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance_s, self.resonance_i = self.peak.resonances

        self.ppm_to_rads_i = self.h_larmor_frq * two_pi * constants.xi_ratio[self.resonance_i['atom']]
        self.ppm_to_rads_s = self.h_larmor_frq * two_pi * constants.xi_ratio[self.resonance_s['atom']]

        kwargs1 = {'temperature': self.temperature}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance_i['name']}
        kwargs3 = {'temperature': self.temperature, 'nuclei': self.resonance_s['name']}
        kwargs4 = {'temperature': self.temperature, 'nuclei': self.peak.assignment, 'h_larmor_frq': self.h_larmor_frq}

        self.map_names = {
            'pb': ParameterName('pb', **kwargs1).to_full_name(),
            'kex_ab': ParameterName('kex_ab', **kwargs1).to_full_name(),
            'dw_i_ab': ParameterName('dw_ab', **kwargs2).to_full_name(),
            'dw_s_ab': ParameterName('dw_ab', **kwargs3).to_full_name(),
            'lambda_mq_a': ParameterName('lambda_mq_a', **kwargs4).to_full_name(),
            'dlambda_mq_ab': ParameterName('dlambda_mq_ab', **kwargs4).to_full_name(),
            'mu_mq_a': ParameterName('mu_mq_a', **kwargs4).to_full_name(),
        }

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
            (self.map_names['dw_s_ab'], 0.0, True, None, None, None),
            (self.map_names['lambda_mq_a'], 10.0, True, 0.0, None, None),
            (self.map_names['dlambda_mq_ab'], 0.0, False, None, None, None),
            (self.map_names['mu_mq_a'], 0.0, False, 0.0, None, None),
        )

        return parameters

    def _calculate_profile(self, pb, kex_ab, dw_i_ab, dw_s_ab,
                           lambda_mq_a, dlambda_mq_ab, mu_mq_a):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb, pc : float
            Fractional population of state B and C.
        kex_ab, kex_bc, kex_ac : float
            Exchange rates between states A, B and C in /s.
        dw_i_ab, dw_i_ac : float
            Chemical shift difference between states A and B, A and C in rad/s.
        rho_i_a : float
            Longitudinal relaxation rate of states A in /s.
        lamdda_i_a : float
            Transverse relaxation rate of state A in /s.
        dlambda_i_ab, dlambda_i_ac : float
            Transverse relaxation rate difference between states A and B, A and C in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block
        """

        omega_i_a = 0.0
        omega_i_b = dw_i_ab * self.ppm_to_rads_i
        omega_s_a = 0.0
        omega_s_b = dw_s_ab * self.ppm_to_rads_s
        lambda_mq_b = dlambda_mq_ab

        pa = 1.0 - pb
        mag_eq = np.array([[0.0, pa, 0.0, 0.0,
                            0.0, pb, 0.0, 0.0]]).T


        # Liouvillians

        kwargs = {
            'pb': pb,
            'kex_ab': kex_ab,
            'omega_i_a': omega_i_a, 'omega_s_a': omega_s_a, 'lambda_mq_a': lambda_mq_a, 'mu_mq_a': mu_mq_a,
            'omega_i_b': omega_i_b, 'omega_s_b': omega_s_b, 'lambda_mq_b': lambda_mq_b, 'mu_mq_b': mu_mq_a,
        }

        l_free = compute_liouvillian(**kwargs)

        # Propagators

        time_series = list(self.tau_cp_list)

        if self.small_protein_flg:
            time_series.append(t_zeta)

        p_free_list = compute_propagator_from_time_series(l_free, time_series)

        if self.small_protein_flg:
            p_zeta = reduce(dot, [p_free_list[t_zeta], p_180x_i, p_180x_s, p_free_list[t_zeta]])
        else:
            p_zeta = None

        # Simulate the CPMG block as function of ncyc

        profile = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):

            if ncyc == 0:

                mag = mag_eq

            else:

                p_free = p_free_list[tau_cp]
                p_cp = matrix_power(dot(dot(p_free, p_180y_s), p_free), int(ncyc))
                mag = reduce(dot, [p_cp, p_180x_i, p_cp, mag_eq])

            if self.small_protein_flg:
                mag = p_zeta.dot(mag)

            profile.append(mag[1, 0])

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
