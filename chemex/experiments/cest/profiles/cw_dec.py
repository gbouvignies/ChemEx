"""15N CEST with CW decoupling

Analyzes chemical exchange in the presence of 1H CW decoupling during the CEST
block. Magnetization evolution is calculated using the 30*30 two-spin matrix:

[ I{xyz}, S{xyz}, 2I{xyz}S{xyz} ]{a, b, ...}

Notes
-----
The calculation is designed specifically to analyze the experiment found in the
reference:

J Phys Chem B (2012), 116, 14311-7
"""

import numpy as np
from scipy import linalg

from chemex.experiments import base_profile
from chemex.experiments.cest import cest_profile


class Profile(cest_profile.CESTProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):

        super().__init__(profile_name, measurements, exp_details)

        self.carrier_dec = base_profile.check_par(exp_details, 'carrier_dec', float)
        omega1_dec = base_profile.check_par(exp_details, 'b1_frq_dec', float) * 2.0 * np.pi

        if '3st' in self.model:
            from chemex.bases.three_state import ixyzsxyz
            self.base = ixyzsxyz
        else:
            from chemex.bases.two_state import ixyzsxyz
            self.base = ixyzsxyz

        self.liouvillians_b1 = np.array([self.base.compute_liouvillian(omega1x_i=omega1, omega1x_s=omega1_dec)
                                         for omega1 in self.omega1s])

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.peak,
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total,
        )

        r2a_i_a = '{r2_i_a} + {r1a_a} - {r1_i_a}'.format(**self.map_names)
        r2a_i_b = '{r2_i_b} + {r1a_b} - {r1_i_b}'.format(**self.map_names)
        r2a_s_a = '{r2_s_a} - {r1_i_a}'.format(**self.map_names)
        r2a_s_b = '{r2_s_b} - {r1_i_b}'.format(**self.map_names)

        self.default_params.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['r2a_i_a'], 0.0, False, 0.0, None, r2a_i_a),
            (self.map_names['r2a_i_b'], 0.0, False, 0.0, None, r2a_i_b),
            (self.map_names['r2a_s_a'], 0.0, False, 0.0, None, r2a_s_a),
            (self.map_names['r2a_s_b'], 0.0, False, 0.0, None, r2a_s_b),
        )

        if '3st' in self.model:
            r2a_i_c = '{r2_i_c} + {r1a_c} - {r1_i_c}'.format(**self.map_names)
            r2a_s_c = '{r2_s_c} - {r1_i_c}'.format(**self.map_names)

            self.default_params.add_many(
                # Name, Value, Vary, Min, Max, Expr
                (self.map_names['r2a_i_c'], 0.0, False, 0.0, None, r2a_i_c),
                (self.map_names['r2a_s_c'], 0.0, False, 0.0, None, r2a_s_c),
            )

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
        omega_i_cars = (cs_i - self.carrier) * self.ppm_i

        cs_s = np.array([kwargs.get(key, 0.0) for key in ('cs_s_a', 'cs_s_b', 'cs_s_c', 'cs_s_d')])
        omega_s_a, omega_s_b, omega_s_c, omega_s_d = (cs_s - self.carrier_dec) * self.ppm_s

        start = self.base.index_iz
        end = self.base.index_iz_a
        mesh = np.ix_(end, start)

        mag0 = self.base.compute_equilibrium(**kwargs)

        profile = []

        for index, b1_offset in enumerate(self.b1_offsets):

            if self.reference[index]:

                magz_a = mag0[end]

            else:

                omega_i_a, omega_i_b, omega_i_c, omega_i_d = omega_i_cars - 2.0 * np.pi * b1_offset

                louvillian_nob1 = self.base.compute_liouvillian(
                    omega_i_a=omega_i_a, omega_i_b=omega_i_b, omega_i_c=omega_i_c, omega_i_d=omega_i_d,
                    omega_s_a=omega_s_a, omega_s_b=omega_s_b, omega_s_c=omega_s_c, omega_s_d=omega_s_d,
                    **kwargs
                )

                liouvillians = self.liouvillians_b1 + louvillian_nob1

                if self.b1_inh == np.inf:
                    s, vr = linalg.eig(liouvillians[0])
                    vri = linalg.inv(vr)

                    sl = [i for i, omega_i_eig in enumerate(s.imag) if abs(omega_i_eig) < 0.1]

                    propagator = (vr[np.ix_(end, sl)]
                                  .dot(np.diag(np.exp(s[sl] * self.time_t1)))
                                  .dot(vri[np.ix_(sl, start)]))

                else:
                    propagators = [linalg.expm(liouvillian * self.time_t1) for liouvillian in liouvillians]
                    propagator = np.average(propagators, weights=self.b1_weights, axis=0)[mesh]

                magz_a = propagator.dot(mag0[start]).real

            profile.append(np.float64(magz_a))

        return np.asarray(profile)
