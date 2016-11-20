"""Pure anti-phase CEST.

Analyzes chemical exchange during the CEST block. This is calculated using the
12x12, single spin matrix:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b) ]

"""

import numpy as np
from scipy import linalg

from chemex import parameters
from chemex.experiments.cest import cest_profile


class Profile(cest_profile.CESTProfile):
    """TODO: class docstring."""

    def __init__(self, profile_name, measurements, exp_details):
        super().__init__(profile_name, measurements, exp_details)

        if '3st' in self.model:
            from chemex.bases.three_state import ixyzsz
            self.base = ixyzsz
        else:
            from chemex.bases.two_state import ixyzsz
            self.base = ixyzsz

        self.liouvillians_b1 = np.array(
            [self.base.compute_liouvillian(omega1x_i=omega1) for omega1 in self.omega1s])

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

        start = self.base.index_2izsz_a
        end = self.base.index_2izsz_a
        mesh = np.ix_(end, start)

        # 2HzNz
        mag0 = self.base.compute_2izsz_a(**kwargs)

        profile = []

        for index, b1_offset in enumerate(self.b1_offsets):
            if self.reference[index]:
                magz_a = mag0[end]

            else:
                omega_i_a, omega_i_b, omega_i_c, omega_i_d = omega_i_cars - 2.0 * np.pi * b1_offset

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

                    propagator = (vr[np.ix_(end, sl)].dot(np.diag(np.exp(s[sl] * self.time_t1)))
                                  .dot(vri[np.ix_(sl, start)]))

                else:
                    propagators = [
                        linalg.expm(liouvillian * self.time_t1) for liouvillian in liouvillians
                    ]
                    propagator = np.average(propagators, weights=self.b1_weights, axis=0)[mesh]

                magz_a = propagator.dot(mag0[start]).real

            profile.append(np.float64(magz_a))

        return np.asarray(profile)
