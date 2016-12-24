"""13CO - Pure Anti-phase Carbonyl 13C CPMG (2-state)

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 12x12,
2-spin exchange matrix:

[ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
  COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b)]

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

import numpy as np
from scipy import linalg

from chemex.util import expmm
from chemex import parameters
from chemex.bases import util
from chemex.experiments import base_profile
from chemex.experiments.cpmg import cpmg_profile


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):
        super().__init__(profile_name, measurements, exp_details)

        self.carrier = base_profile.check_par(exp_details, 'carrier', float)
        self.time_equil = base_profile.check_par(exp_details, 'time_equil', float)
        self.sidechain_flg = base_profile.check_par(exp_details, 'sidechain_flg', str).lower()
        self.taucc = base_profile.check_par(exp_details, 'taucc', float, default=9.09e-3)

        self.tau_cp_list = np.array(
            [self.time_t2 / (4.0 * ncyc) if ncyc > 0 else 0.0 for ncyc in self.ncycs])
        self.t_neg = -2.0 * self.pw / np.pi
        self.time_series = [self.t_neg, self.taucc, self.time_equil]
        self.time_series.extend(self.tau_cp_list)

        if '3st' in self.model:
            from chemex.bases.three_state import ixyzsz
            self.base = ixyzsz
        else:
            from chemex.bases.two_state import ixyzsz
            self.base = ixyzsz

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.peak,
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total, )

        kwargs = {
            'temperature': self.temperature,
            'nuclei': self.resonance_s['name'],
            'h_larmor_frq': self.h_larmor_frq
        }

        self.map_names['r1_s_a'] = parameters.ParameterName('r1_a', **kwargs).to_full_name()
        self.map_names['r1_s_b'] = parameters.ParameterName('r1_b', **kwargs).to_full_name()

        r1_s_b = self.map_names['r1_s_a']
        r1_i_a = '{r1a_a} - {r1_s_a}'.format(**self.map_names)
        r1_i_b = '{r1a_b} - {r1_s_b}'.format(**self.map_names)
        r2a_i_a = '{r2_i_a} - {r1_s_a}'.format(**self.map_names)
        r2a_i_b = '{r2_i_b} - {r1_s_b}'.format(**self.map_names)

        self.default_params.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['r1_s_a'], 1.5, False, 0.0, None, None),
            (self.map_names['r1_s_b'], 1.5, None, 0.0, None, r1_s_b),
            (self.map_names['r1_i_a'], 0.0, None, 0.0, None, r1_i_a),
            (self.map_names['r1_i_b'], 0.0, None, 0.0, None, r1_i_b),
            (self.map_names['r2a_i_a'], 0.0, None, 0.0, None, r2a_i_a),
            (self.map_names['r2a_i_b'], 0.0, None, 0.0, None, r2a_i_b), )

        if '3st' in self.model:
            self.map_names['r1_s_c'] = parameters.ParameterName('r1_c', **kwargs).to_full_name()

            r1_s_c = self.map_names['r1_s_a']
            r1_i_c = '{r1a_c} - {r1_s_c}'.format(**self.map_names)
            r2a_i_c = '{r2_i_c} - {r1_s_c}'.format(**self.map_names)

            self.default_params.add_many(
                # Name, Value, Vary, Min, Max, Expr
                (self.map_names['r1_s_c'], 1.5, None, 0.0, None, r1_s_c),
                (self.map_names['r1_i_c'], 0.0, None, 0.0, None, r1_i_c),
                (self.map_names['r2a_i_c'], 0.0, None, 0.0, None, r2a_i_c), )

    def calculate_unscaled_profile(self, **kwargs):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb
            Fractional population of state B and C.
        kex_ab
            Exchange rates between states A and B in /s.
        dw_i_ab
            Chemical shift difference between states A and B in rad/s.
        r1_i_a
            Longitudinal relaxation rate of states A in /s.
        r2_i_a
            Transverse relaxation rate of state A in /s.
        dr2_i_ab
            Transverse relaxation rate difference between states A and B in /s.
        cs_i_a
            Resonance position of state A in ppm.

        Returns
        -------
        out
            Intensity after the CEST block

        """
        cs_i = np.array([kwargs.get(key, 0.0) for key in ('cs_i_a', 'cs_i_b', 'cs_i_c', 'cs_i_d')])
        omega_i_a, omega_i_b, omega_i_c, omega_i_d = (cs_i - self.carrier) * self.ppm_i

        # Liouvillians
        l_free = self.base.compute_liouvillian(
            omega_i_a=omega_i_a,
            omega_i_b=omega_i_b,
            omega_i_c=omega_i_c,
            omega_i_d=omega_i_d,
            **kwargs)
        l_pw1y = l_free + self.base.compute_liouvillian(omega1y_i=+self.omega1_i)
        l_mw1y = l_free + self.base.compute_liouvillian(omega1y_i=-self.omega1_i)

        # Propagators
        p_90py = expmm(l_pw1y * self.pw)
        p_90my = expmm(l_mw1y * self.pw)
        p_180py = np.linalg.matrix_power(p_90py, 2)
        p_180my = np.linalg.matrix_power(p_90my, 2)
        p_180pmy = 0.5 * (p_180py + p_180my)
        p_180x_perfect = self.base.p_180x_i_perfect

        p_free_list = util.compute_propagators_from_time_series(l_free, self.time_series)
        p_neg = p_free_list[self.t_neg]
        p_taucc = p_free_list[self.taucc]
        p_equil = p_free_list[self.time_equil]

        # 2CzNz
        mag0 = self.base.compute_equilibrium_2izsz(**kwargs)

        # Simulate the CPMG block as function of ncyc
        profile = []

        if self.sidechain_flg == 'y':
            p_flip = p_180pmy
        else:
            p_flip = functools.reduce(np.dot, [p_90my, p_taucc, p_180pmy, p_taucc, p_90py])

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):
            if ncyc == 0:
                mag = functools.reduce(np.dot, [p_equil, p_90py, p_flip, p_90py, mag0])

            else:
                p_free = p_free_list[tau_cp]
                p_cpx = np.linalg.matrix_power(p_free.dot(p_180x_perfect).dot(p_free), int(ncyc))
                mag = functools.reduce(np.dot, [
                    p_equil, p_90py, p_neg, p_cpx, p_neg, p_flip, p_neg, p_cpx, p_neg, p_90py, mag0
                ])

            profile.append(np.float64(mag[self.base.index_2izsz_a]))

        return np.asarray(profile)
