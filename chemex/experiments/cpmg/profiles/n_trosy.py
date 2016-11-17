"""15N - N-H TROSY CPMG

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the 12x12, two spin matrix:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b) ]

Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:

Proc Natl Acad Sci USA (2007) 104, 18473-7
"""

import functools

import numpy as np
from scipy import linalg

from chemex.bases import util
from chemex.experiments import base_profile
from chemex.experiments.cpmg import cpmg_profile


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):
        super().__init__(profile_name, measurements, exp_details)

        self.carrier = base_profile.check_par(exp_details, 'carrier', float)
        self.time_equil = base_profile.check_par(exp_details, 'time_equil', float)
        self.taub = base_profile.check_par(exp_details, 'taub', float)

        self.t_neg = -2.0 * self.pw / np.pi
        self.time_series = [self.t_neg, self.time_equil, self.taub]
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

        r2a_i_a = '{r2_i_a} + {r1a_a} - {r1_i_a}'.format(**self.map_names)
        r2a_i_b = '{r2_i_b} + {r1a_b} - {r1_i_b}'.format(**self.map_names)

        self.default_params.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['r2a_i_a'], 0.0, None, 0.0, None, r2a_i_a),
            (self.map_names['r2a_i_b'], 0.0, None, 0.0, None, r2a_i_b), )

        if '3st' in self.model:
            r2a_i_c = '{r2_i_c} + {r1a_c} - {r1_i_c}'.format(**self.map_names)

            self.default_params.add_many(
                # Name, Value, Vary, Min, Max, Expr
                (self.map_names['r2a_i_c'], 0.0, None, 0.0, None, r2a_i_c), )

    def calculate_unscaled_profile(self, **kwargs):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex_ab : float
            Exchange rates between states A and B in /s.
        dw_i_ab : float
            Chemical shift difference between states A and B in rad/s.
        r1_i_a : float
            Longitudinal relaxation rate of states A in /s.
        r2_i_a : float
            Transverse relaxation rate of state A in /s.
        dr2_i_ab : float
            Transverse relaxation rate difference between states A and B in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
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
        l_pw1x = l_free + self.base.compute_liouvillian(omega1x_i=+self.omega1_i)
        l_pw1y = l_free + self.base.compute_liouvillian(omega1y_i=+self.omega1_i)
        l_mw1x = l_free + self.base.compute_liouvillian(omega1x_i=-self.omega1_i)
        l_mw1y = l_free + self.base.compute_liouvillian(omega1y_i=-self.omega1_i)

        # Propagators
        p_90px = linalg.expm(l_pw1x * self.pw)
        p_90py = linalg.expm(l_pw1y * self.pw)
        p_90mx = linalg.expm(l_mw1x * self.pw)
        p_90my = linalg.expm(l_mw1y * self.pw)
        p_180px = np.linalg.matrix_power(p_90px, 2)
        p_180py = np.linalg.matrix_power(p_90py, 2)
        p_180x_s = self.base.p_180x_s

        p_free_list = util.compute_propagators_from_time_series(l_free, self.time_series)

        p_equil = p_free_list[self.time_equil]
        p_neg = p_free_list[self.t_neg]
        p_taub = p_free_list[self.taub]

        p_element = functools.reduce(
            np.dot, [p_180x_s, p_taub, p_90py, p_90px, p_180x_s, p_90px, p_90py, p_taub])
        p_element_pc = 0.5 * (p_90px.dot(p_element).dot(p_90py) + p_90mx.dot(p_element).dot(p_90my))

        # 2HzNz
        mag0 = self.base.compute_equilibrium_2izsz(**kwargs)

        # Simulate the CPMG block as function of ncyc
        profile = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):
            if ncyc == 0.0:
                mag = functools.reduce(np.dot, [p_equil, p_90py, p_element, p_90px, mag0])

            else:
                p_free = p_free_list[tau_cp]
                p_cpx = np.linalg.matrix_power(p_free.dot(p_180px).dot(p_free), int(ncyc))
                p_cpy = np.linalg.matrix_power(p_free.dot(p_180py).dot(p_free), int(ncyc))
                mag = functools.reduce(np.dot, [
                    p_equil, p_90py, p_neg, p_cpx, p_neg, p_element_pc, p_neg, p_cpy, p_neg, p_90px,
                    mag0
                ])

            # Trosy (A)
            profile.append(mag[5, 0] - mag[2, 0])

        return np.asarray(profile)
