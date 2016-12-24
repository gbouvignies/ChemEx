"""Pure in-phase CPMG.

Analyzes chemical exchange in the presence of high power 1H CW decoupling during
the CPMG block. This keeps the spin system purely in-phase throughout, and is
calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]

Notes
-----
Off resonance effects are taken into account.

The calculation is designed specifically to analyze the experiment found in the
reference:

Journal of Physical Chemistry B (2008), 112, 5898-5904

"""

import functools

import numpy as np
from scipy import linalg

from chemex.util import expmm
from chemex.bases import util
from chemex.experiments import base_profile
from chemex.experiments.cpmg import cpmg_profile


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):
        super().__init__(profile_name, measurements, exp_details)

        self.carrier = base_profile.check_par(exp_details, 'carrier', float)
        self.time_equil = base_profile.check_par(exp_details, 'time_equil', float)

        self.t_neg = -2.0 * self.pw / np.pi
        self.time_series = [self.t_neg, self.time_equil]
        self.time_series.extend(self.tau_cp_list)

        if '3st' in self.model:
            from chemex.bases.three_state import ixyz
            self.base = ixyz
        else:
            from chemex.bases.two_state import ixyz
            self.base = ixyz

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.resonance_i['name'],
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total, )

        self.default_params[self.map_names['r1_i_a']].set(vary=False)

        for name in ('r2_i_b', 'r2_i_c', 'r2_i_d'):
            if name in self.map_names:
                param = self.default_params[self.map_names[name]]
                param.set(min=param.min, max=param.max, expr=self.map_names['r2_i_a'])

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

        # Calculation of the different liouvillians
        l_free = self.base.compute_liouvillian(
            omega_i_a=omega_i_a,
            omega_i_b=omega_i_b,
            omega_i_c=omega_i_c,
            omega_i_d=omega_i_d,
            **kwargs)
        l_pw1x = l_free + self.base.compute_liouvillian(omega1x_i=+self.omega1_i)
        l_pw1y = l_free + self.base.compute_liouvillian(omega1y_i=+self.omega1_i)
        l_mw1x = l_free + self.base.compute_liouvillian(omega1x_i=-self.omega1_i)

        # Calculation of all the needed propagators
        p_90px = expmm(l_pw1x, self.pw)
        p_180px = np.linalg.matrix_power(p_90px, 2)
        p_180py = expmm(l_pw1y, 2.0 * self.pw)
        p_180mx = expmm(l_mw1x, 2.0 * self.pw)
        p_180pmx = 0.5 * (p_180px + p_180mx)  # +/- phase cycling

        p_free_list = util.compute_propagators_from_time_series(l_free, self.time_series)

        p_equil = p_free_list[self.time_equil]
        p_neg = p_free_list[self.t_neg]

        # Simulate the CPMG block as function of ncyc
        mag0 = self.base.compute_equilibrium(**kwargs)

        profile = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):
            if ncyc == 0:
                mag = functools.reduce(np.dot, [p_equil, p_90px, p_180pmx, p_90px, mag0])

            elif ncyc == -1:
                p_free = p_free_list[tau_cp]
                mag = functools.reduce(
                    np.dot,
                    [p_equil, p_90px, p_neg, p_free, p_180pmx, p_free, p_neg, p_90px, mag0])

            else:
                p_free = p_free_list[tau_cp]
                p_cp = np.linalg.matrix_power(np.dot(np.dot(p_free, p_180py), p_free), int(ncyc))
                mag = functools.reduce(
                    np.dot, [p_equil, p_90px, p_neg, p_cp, p_180pmx, p_cp, p_neg, p_90px, mag0])

            profile.append(np.float64(mag[self.base.index_iz_a]))

        return np.asarray(profile)
