"""13C(methyl) - H to C CPMG

Measures methyl carbon chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initally anti-phase and is read out as in-phase. Because of the P-element only
even ncyc should be recorded. The calculation uses a 12x12 basis set:

[Cx(a), Cy(a), Cz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
 Cx(b), Cy(b), Cz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

Off resonance effects are taken into account. The calculation is designed
explicitly for analyzing the Lewis Kay pulse sequence:

HtoC_CH3_exchange_*00_lek_ILV

Journal of Biomolecular NMR (2007) 38, 79-88
"""
import numpy as np
from numpy import linalg as la

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG2

_EXP_DETAILS = {"taub": {"default": 1.99e-3, "type": float}}


class ProfileCPMGCH3H2C(ProfileCPMG2):
    """TODO: class docstring."""

    EXP_DETAILS = dict(**ProfileCPMG2.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixyzsz"
    CONSTRAINTS = "nh"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.taub = self.exp_details["taub"]

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the delays in the experiments
        self.delays += [self.taub]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: Write docstring"""

        self.liouv.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_taub = delays[self.taub]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180_s = self.liouv.perfect180["sx"]

        # Calculate starting magnetization vector
        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")

        palmer = d_taub @ p90[0] @ p180_s @ p90[0] @ d_taub

        # Calculating the cpmg trains
        cp1 = {0: self.liouv.identity}
        cp2 = {0: self.liouv.identity}

        for ncyc in set(self.data["ncycs"][~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180[[1, 0]] @ tau_cp
            cp_trains = la.matrix_power(echo, int(ncyc))
            cp1[ncyc] = cp_trains[0] @ d_neg
            cp2[ncyc] = d_neg @ cp_trains[1]

        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ p90[1]
                @ cp2[ncyc]
                @ palmer
                @ cp1[ncyc]
                @ p90[0]
                @ mag0
            )
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)
