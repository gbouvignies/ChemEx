"""15N - N-H anti-TROSY CPMG

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the 12x12, two spin matrix:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b) ]

Note
----
Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:

Proc Natl Acad Sci USA (2007) 104, 18473-7
"""

import functools

import numpy as np
from scipy import linalg

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default, util


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details):
        super().__init__(name, measurements, exp_details)

        self.carrier = self.check_par(exp_details, "carrier", float)
        self.time_equil = self.check_par(exp_details, "time_equil", float)
        self.taub = self.check_par(exp_details, "taub", float)

        self.basis = basis.Basis(
            model_name=self.model,
            spin_state="ixyzsz",
            peak=self.peak,
            equilibrium=False,
        )

        self.map_names, self.default_params = default.create_params(
            basis=self.basis,
            model=self.model,
            nuclei=self.peak,
            nh_constraints=True,
            **self.conditions,
        )

        for name in ("dw_i_ab", "dw_i_ac", "dw_i_ad", "r2_i_a"):
            if name in self.map_names:
                self.default_params[self.map_names[name]].set(vary=True)

        self.detect = (
            self.basis.magnetizations["2izsz_a"] + self.basis.magnetizations["iz_a"]
        ).T

        self.t_neg = -2.0 * self.pw / np.pi
        self.time_series = [self.t_neg, self.time_equil, self.taub]
        self.time_series.extend(self.tau_cp_list)

        self.liouvillian_carrier = self.basis.compute_liouvillian(
            carrier_i=self.carrier
        )

        self.liouvillian_b1_ix = self.basis.compute_liouvillian(w1x_i=self.w1_i)
        self.liouvillian_b1_iy = self.basis.compute_liouvillian(w1y_i=self.w1_i)

    def calculate_unscaled_profile(self, **parvals):
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

        # Calculation of the different liouvillians
        l_free = self.liouvillian_carrier + self.basis.compute_liouvillian(**parvals)
        l_pw1x = l_free + self.liouvillian_b1_ix
        l_pw1y = l_free + self.liouvillian_b1_iy
        l_mw1x = l_free - self.liouvillian_b1_ix
        l_mw1y = l_free - self.liouvillian_b1_iy

        # Calculation of all the needed propagators
        p_90px = linalg.expm(l_pw1x * self.pw)
        p_90py = linalg.expm(l_pw1y * self.pw)
        p_90mx = linalg.expm(l_mw1x * self.pw)
        p_90my = linalg.expm(l_mw1y * self.pw)
        p_180px = np.linalg.matrix_power(p_90px, 2)
        p_180py = np.linalg.matrix_power(p_90py, 2)
        p_180x_s = self.basis.pp180["sx"]

        p_free_list = util.compute_propagators_from_time_series(
            l_free, self.time_series
        )

        p_equil = p_free_list[self.time_equil]
        p_neg = p_free_list[self.t_neg]
        p_taub = p_free_list[self.taub]

        # Simulate the CPMG block as function of ncyc
        mag0 = self.basis.compute_mag_eq(term="2izsz", **parvals)

        p_element = functools.reduce(
            np.dot, [p_180x_s, p_taub, p_90mx, p_90py, p_180x_s, p_90py, p_90mx, p_taub]
        )

        p_element_pc = 0.5 * (
            p_90px.dot(p_element).dot(p_90py) + p_90mx.dot(p_element).dot(p_90my)
        )

        intensities = []

        for ncyc, tau_cp in zip(self.ncycs, self.tau_cp_list):

            if ncyc == 0:
                mag = functools.reduce(
                    np.dot, [self.detect, p_equil, p_90py, p_element, p_90px, mag0]
                )

            else:
                p_free = p_free_list[tau_cp]
                p_cpx = np.linalg.matrix_power(
                    p_free.dot(p_180px).dot(p_free), int(ncyc)
                )
                p_cpy = np.linalg.matrix_power(
                    p_free.dot(p_180py).dot(p_free), int(ncyc)
                )
                mag = functools.reduce(
                    np.dot,
                    [
                        self.detect,
                        p_equil,
                        p_90py,
                        p_neg,
                        p_cpx,
                        p_neg,
                        p_element_pc,
                        p_neg,
                        p_cpy,
                        p_neg,
                        p_90px,
                        mag0,
                    ],
                )

            intensities.append(np.float64(mag))

        return np.asarray(intensities)
