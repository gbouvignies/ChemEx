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

import functools

import numpy as np

from chemex.experiments import base_profile
from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import util


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details):
        super().__init__(name, measurements, exp_details)

        self.small_protein_flg = (
            base_profile.check_par(exp_details, "small_protein_flg", str).lower() == "y"
        )
        self.t_zeta = 1.0 / (8.0 * 125.3)
        self.time_series = [self.t_zeta]
        self.time_series.extend(self.tau_cp_list)

        if "3st" in self.model:
            from chemex.spindynamics.three_state import ixysxy

            self.base = ixysxy
        else:
            from chemex.spindynamics.two_state import ixysxy

            self.base = ixysxy

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.peak,
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total,
        )

    def calculate_unscaled_profile(self, **kwargs):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb, pc : float
            Fractional population of state B and C.
        kex_ab, kex_bc, kex_ac : float
            Exchange rates between states A, B and C in /s.
        dw_i_ab, dw_i_ac : float
            Chemical shift difference between states A and B, A and C in rad/s.
        r1_i_a : float
            Longitudinal relaxation rate of states A in /s.
        r2_i_a : float
            Transverse relaxation rate of state A in /s.
        dr2_i_ab, dr2_i_ac : float
            Transverse relaxation rate difference between states A and B, A and C in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """
        omega_i_a, omega_i_b, omega_i_c, omega_i_d = (
            np.array(
                [
                    kwargs.get(key, 0.0)
                    for key in ("cs_i_a", "cs_i_b", "cs_i_c", "cs_i_d")
                ]
            )
            * self.ppm_i
        )

        omega_s_a, omega_s_b, omega_s_c, omega_s_d = (
            np.array(
                [
                    kwargs.get(key, 0.0)
                    for key in ("cs_s_a", "cs_s_b", "cs_s_c", "cs_s_d")
                ]
            )
            * self.ppm_s
        )

        # Liouvillians
        l_free = self.base.compute_liouvillian(
            omega_i_a=omega_i_a,
            omega_i_b=omega_i_b,
            omega_i_c=omega_i_c,
            omega_i_d=omega_i_d,
            omega_s_a=omega_s_a,
            omega_s_b=omega_s_b,
            omega_s_c=omega_s_c,
            omega_s_d=omega_s_d,
            **kwargs
        )

        # Propagators
        p_180x_s = self.base.p_180x_s
        p_180x_i = self.base.p_180x_i
        p_180y_i = self.base.p_180y_i

        p_free_list = util.compute_propagators_from_time_series(
            l_free, self.time_series
        )
        p_zeta = p_free_list[self.t_zeta]

        # 2HxCy
        mag0 = self.base.compute_equilibrium_iysx(**kwargs)

        # Simulate the CPMG block as function of ncyc
        profile = []

        if self.small_protein_flg:
            mag0 = functools.reduce(np.dot, [p_zeta, p_180x_s, p_180x_i, p_zeta, mag0])

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):
            if ncyc == 0:
                mag = mag0

            else:
                p_free = p_free_list[tau_cp]
                p_cp = np.linalg.matrix_power(
                    np.dot(np.dot(p_free, p_180y_i), p_free), int(ncyc)
                )
                mag = functools.reduce(np.dot, [p_cp, p_180x_s, p_cp, mag0])

            # 2HxCy[A]
            profile.append(np.float64(mag[self.base.index_iysx_a]))

        return np.asarray(profile)
