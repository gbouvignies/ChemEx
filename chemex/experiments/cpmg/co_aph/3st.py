# -*- coding: utf-8 -*-

"""13CO - Pure Anti-phase Carbonyl 13C CPMG (3-state)

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 12x12,
2-spin exchange matrix:

[ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
  COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b),
  COx(c), COy(c), COz(c), 2COxNz(c), 2COyNz(c), 2COzNz(c)]

Because of the length of the shaped pulses used during the CPMG blocks, off-
resonance effects are taken into account only for the 90-degree pulses that
create COxNz before the CPMG and COzNz after the CPMG. The calculation is
designed explicitly for analyzing the Kay laboratory pulse sequence:

CO_CPMG_SCFilter_x00_dfh1

And can be run with or without sidechain CO inversion via the Inv_CO flag for
uniformly 13C-labeled proteins.

Journal of Biomolecular NMR (2008) 42, 35-47
"""

import functools

import lmfit
import numpy as np
from scipy import linalg

from chemex import constants, parameters, peaks
from chemex.bases import iph_aph_3st, util as bases_util
from chemex.experiments import base_profile
from chemex.experiments.cpmg import plotting

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

reduce = functools.reduce
dot = np.dot
matrix_power = np.linalg.matrix_power
compute_liouvillian = iph_aph_3st.compute_liouvillian
check_par = base_profile.check_par
compute_propagator_from_time_series = bases_util.compute_propagators_from_time_series
ParameterName = parameters.ParameterName

two_pi = 2.0 * np.pi
p_180x_i = np.diag([1.0, -1.0, 1.0, 1.0, -1.0, 1.0,
                    1.0, -1.0, 1.0, 1.0, -1.0, 1.0,
                    1.0, -1.0, 1.0, 1.0, -1.0, 1.0])

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'carrier': float,
    'pw': float,
    'time_t2': float,
    'time_equil': float,
    'taucc': float,
    'sidechain_flg': str
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
        self.sidechain_flg = check_par(exp_details, 'sidechain_flg', str).lower()
        self.taucc = check_par(exp_details, 'taucc', float, default=9.09e-3)

        self.experiment_name = check_par(exp_details, 'experiment_name')
        self.tau_cp_list = np.array([self.time_t2 / (4.0 * ncyc) if ncyc else 0.0 for ncyc in self.ncycs])

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance_i, self.resonance_s = self.peak.resonances

        self.ppm_to_rads = self.h_larmor_frq * two_pi * constants.xi_ratio[self.resonance_i['atom']]

        kwargs1 = {'temperature': self.temperature}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance_i['name']}
        kwargs3 = {'temperature': self.temperature, 'nuclei': self.resonance_i['name'],
                   'h_larmor_frq': self.h_larmor_frq}
        kwargs4 = {'temperature': self.temperature, 'nuclei': self.peak.assignment}
        kwargs5 = {'temperature': self.temperature, 'nuclei': self.peak.assignment, 'h_larmor_frq': self.h_larmor_frq}
        kwargs6 = {'temperature': self.temperature, 'nuclei': self.resonance_s['name'],
                   'h_larmor_frq': self.h_larmor_frq}

        self.map_names = {
            'pb': ParameterName('pb', **kwargs1).to_full_name(),
            'pc': ParameterName('pc', **kwargs1).to_full_name(),
            'kex_ab': ParameterName('kex_ab', **kwargs1).to_full_name(),
            'kex_bc': ParameterName('kex_bc', **kwargs1).to_full_name(),
            'kex_ac': ParameterName('kex_ac', **kwargs1).to_full_name(),
            'cs_i_a': ParameterName('cs_a', **kwargs2).to_full_name(),
            'dw_i_ab': ParameterName('dw_ab', **kwargs2).to_full_name(),
            'dw_i_ac': ParameterName('dw_ac', **kwargs2).to_full_name(),
            'lambda_i_a': ParameterName('lambda_a', **kwargs3).to_full_name(),
            'dlambda_i_ab': ParameterName('dlambda_ab', **kwargs3).to_full_name(),
            'dlambda_i_ac': ParameterName('dlambda_ac', **kwargs3).to_full_name(),
            'rho_s_a': ParameterName('rho_a', **kwargs6).to_full_name(),
            'rho_is_a': ParameterName('rho_is_a', **kwargs5).to_full_name(),
            'eta_i_a': ParameterName('eta_a', **kwargs3).to_full_name(),
            'delta_i_a': ParameterName('delta_a', **kwargs3).to_full_name(),
            'j_a': ParameterName('j_a', **kwargs4).to_full_name(),
            'dj_ab': ParameterName('dj_ab', **kwargs4).to_full_name(),
            'dj_ac': ParameterName('dj_ac', **kwargs4).to_full_name(),
        }

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['pc'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (self.map_names['kex_bc'], 200.0, True, 0.0, None, None),
            (self.map_names['kex_ac'], 0.0, False, 0.0, None, None),
            (self.map_names['cs_i_a'], 0.0, False, None, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
            (self.map_names['dw_i_ac'], 0.0, True, None, None, None),
            (self.map_names['lambda_i_a'], 10.0, True, 0.0, None, None),
            (self.map_names['dlambda_i_ab'], 0.0, False, None, None, None),
            (self.map_names['dlambda_i_ac'], 0.0, False, None, None, None),
            (self.map_names['rho_s_a'], 1.0, False, 0.0, None, None),
            (self.map_names['rho_is_a'], 3.0, False, 0.0, None, None),
            (self.map_names['eta_i_a'], 0.0, False, None, None, None),
            (self.map_names['delta_i_a'], 0.0, False, None, None, None),
            (self.map_names['j_a'], 15.0, False, None, None, None),
            (self.map_names['dj_ab'], 0.0, False, None, None, None),
            (self.map_names['dj_ac'], 0.0, False, None, None, None),
        )

        return parameters

    def _calculate_profile(self, pb, pc, kex_ab, kex_bc, kex_ac, cs_i_a, dw_i_ab, dw_i_ac, rho_s_a, lambda_i_a,
                           dlambda_i_ab, dlambda_i_ac, rho_is_a, eta_i_a, delta_i_a, j_a, dj_ab, dj_ac):
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

        domega_i_ab = dw_i_ab * self.ppm_to_rads
        domega_i_ac = dw_i_ac * self.ppm_to_rads
        omega_i_a = (cs_i_a - self.carrier) * self.ppm_to_rads
        omega_i_b = omega_i_a + domega_i_ab
        omega_i_c = omega_i_a + domega_i_ac
        lambda_i_b = lambda_i_a + dlambda_i_ab
        lambda_i_c = lambda_i_a + dlambda_i_ac

        rho_i_a = rho_is_a - rho_s_a
        rhoa_i_a = lambda_i_a - rho_s_a
        rhoa_i_b = lambda_i_b - rho_s_a
        rhoa_i_c = lambda_i_c - rho_s_a

        pa = 1.0 - pb - pc

        # 2HzNz

        mag_eq = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, pa,
                            0.0, 0.0, 0.0, 0.0, 0.0, pb,
                            0.0, 0.0, 0.0, 0.0, 0.0, pc, ]]).T

        w1 = 2.0 * np.pi / (4.0 * self.pw)

        # Liouvillians

        kwargs = {
            'pb': pb, 'pc': pc,
            'kex_ab': kex_ab, 'kex_bc': kex_bc, 'kex_ac': kex_ac,
            'rho_i_a': rho_i_a, 'lambda_i_a': lambda_i_a, 'omega_i_a': omega_i_a,
            'rho_i_b': rho_i_a, 'lambda_i_b': lambda_i_b, 'omega_i_b': omega_i_b,
            'rho_i_c': rho_i_a, 'lambda_i_c': lambda_i_c, 'omega_i_c': omega_i_c,
            'rhoa_i_a': rhoa_i_a, 'rho_is_a': rho_is_a, 'eta_i_a': eta_i_a, 'delta_i_a': delta_i_a, 'j_a': j_a,
            'rhoa_i_b': rhoa_i_b, 'rho_is_b': rho_is_a, 'eta_i_b': eta_i_a, 'delta_i_b': delta_i_a, 'j_b': j_a + dj_ab,
            'rhoa_i_c': rhoa_i_c, 'rho_is_c': rho_is_a, 'eta_i_c': eta_i_a, 'delta_i_c': delta_i_a, 'j_c': j_a + dj_ac,
        }

        l_free = compute_liouvillian(**kwargs)
        l_pw1y = compute_liouvillian(omega1y_i=+w1, **kwargs)
        l_mw1y = compute_liouvillian(omega1y_i=-w1, **kwargs)

        # Propagators

        p_90py = linalg.expm(l_pw1y * self.pw)
        p_90my = linalg.expm(l_mw1y * self.pw)
        p_180py = matrix_power(p_90py, 2)
        p_180my = matrix_power(p_90my, 2)
        p_180pmy = 0.5 * (p_180py + p_180my)

        t_neg = -2.0 * self.pw / np.pi
        time_series = [t_neg, self.taucc, self.time_equil]
        time_series.extend(self.tau_cp_list)

        p_free_list = compute_propagator_from_time_series(l_free, time_series)

        p_neg = p_free_list[t_neg]
        p_taucc = p_free_list[self.taucc]
        p_equil = p_free_list[self.time_equil]

        if self.sidechain_flg == 'y':
            p_flip = p_180pmy
        else:
            p_flip = reduce(dot, [p_90my, p_taucc, p_180pmy, p_taucc, p_90py])

        # Simulate the CPMG block as function of ncyc

        profile = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):

            if ncyc == 0:
                mag = reduce(dot, [p_equil, p_90py, p_flip, p_90py, mag_eq])

            else:
                p_free = p_free_list[tau_cp]
                p_cpx = matrix_power(p_free.dot(p_180x_i).dot(p_free), int(ncyc))
                mag = reduce(dot, [p_equil, p_90py, p_neg, p_cpx, p_neg, p_flip, p_neg, p_cpx, p_neg, p_90py, mag_eq])

            # 2HzNz (A)
            profile.append(mag[5, 0])

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

        iter_vals = zip(self.ncycs, self.val, self.err, values)

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
