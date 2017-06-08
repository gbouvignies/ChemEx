"""Pure in-phase CEST in uniformly C13-labeled samples.

Analyzes 15N chemical exchange in the presence of 1H composite decoupling
during the CEST block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]

Notes
-----
The calculation is designed specifically to analyze the experiment found in
the references:

J Am Chem Soc (2012), 134, 8148-61
Angew Chem (2013), 52, 4156-9
JMB (2014), 426, 763-74

"""

import numpy as np
from scipy import linalg

from chemex import constants
from chemex.experiments.cest import cest_profile, util


class Profile(cest_profile.CESTProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):
        super().__init__(profile_name, measurements, exp_details)

        self.multiplet = util.calc_multiplet(constants.j_couplings[self.resonance_i['symbol']][
            self.resonance_i['nucleus']])

        if '3st' in self.model:
            from chemex.bases.three_state import ixyz
            self.base = ixyz
        else:
            from chemex.bases.two_state import ixyz
            self.base = ixyz

        self.liouvillians_b1 = np.array(
            [self.base.compute_liouvillian(omega1x_i=omega1) for omega1 in self.omega1s])

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.peak,
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total, )

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

        start = self.base.index_iz
        end = self.base.index_iz_a
        mesh = np.ix_(end, start)

        mag0 = self.base.compute_equilibrium(**kwargs)

        profile = []

        for index, b1_offset in enumerate(self.b1_offsets):
            if self.reference[index]:
                magz_a = mag0[end]

            else:
                magz_a = 0.0

                for j, weight in self.multiplet:
                    omega_i_a, omega_i_b, omega_i_c, omega_i_d = omega_i_cars - 2.0 * np.pi * b1_offset + j

                    louvillian_nob1 = self.base.compute_liouvillian(
                        omega_i_a=omega_i_a,
                        omega_i_b=omega_i_b,
                        omega_i_c=omega_i_c,
                        omega_i_d=omega_i_d,
                        **kwargs)

                    liouvillians = self.liouvillians_b1 + louvillian_nob1

                    if self.b1_inh == np.inf:
                        s, vr = linalg.eig(liouvillians[0])
                        vri = linalg.inv(vr)

                        sl = [i for i, omega_i_eig in enumerate(s.imag) if abs(omega_i_eig) < 0.1]

                        propagator = (vr[np.ix_(end, sl)].dot(
                            np.diag(np.exp(s[sl] * self.time_t1))).dot(vri[np.ix_(sl, start)]))

                    else:
                        propagators = [
                            linalg.expm(liouvillian * self.time_t1) for liouvillian in liouvillians
                        ]
                        propagator = np.average(propagators, weights=self.b1_weights, axis=0)[mesh]

                    magz_a += weight * propagator.dot(mag0[start]).real

            profile.append(np.float64(magz_a))

        return np.asarray(profile)
