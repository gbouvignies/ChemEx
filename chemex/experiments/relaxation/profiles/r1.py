"""R1 experiment

[To be filled]

"""

import lmfit
import numpy as np

from chemex.bases import util
from chemex.experiments.relaxation import relaxation_profile


class Profile(relaxation_profile.RelaxationProfile):
    def __init__(self, profile_name, measurements, exp_details):

        super().__init__(profile_name, measurements, exp_details)

        if '3st' in self.model:
            from chemex.bases.three_state import iz
            self.base = iz
        else:
            from chemex.bases.two_state import iz
            self.base = iz

        self.map_names, self.default_params = self.base.create_default_params(
            model=self.model,
            nuclei=self.resonance_i['name'],
            temperature=self.temperature,
            h_larmor_frq=self.h_larmor_frq,
            p_total=self.p_total,
            l_total=self.l_total,
        )

    def create_default_parameters(self):

        params = lmfit.Parameters()

        params.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['kab'], 1.0e9, True, 0.0, None, None),
            (self.map_names['kd_ab'], 100.0, True, 0.0, None, None),
            (self.map_names['r1_i_a'], 1.0, True, 0.0, None, None),
            (self.map_names['r1_i_b'], 1.0, True, 0.0, None, None),
        )

        return params

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

        mag0 = self.base.compute_equilibrium_iz(**kwargs)
        liouvillian = self.base.compute_liouvillian(**kwargs)
        propagators = util.compute_propagators_from_time_series(liouvillian, self.delays)

        profile = []

        for delay in self.delays:
            mag = np.dot(propagators[delay], mag0)

            profile.append(np.float64(mag[self.base.index_iz_a]))

        return np.asarray(profile)
