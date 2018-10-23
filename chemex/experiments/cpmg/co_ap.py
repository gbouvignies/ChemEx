"""13CO - Pure Anti-phase Carbonyl 13C CPMG

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 12x12,
2-spin exchange matrix:

[ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
  COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b)]

Notes
-----

Because of the length of the shaped pulses used during the CPMG blocks, off-
resonance effects are taken into account only for the 90-degree pulses that
create COxNz before the CPMG and COzNz after the CPMG. The calculation is
designed explicitly for analyzing the Kay laboratory pulse sequence:

CO_CPMG_SCFilter_x00_dfh1

And can be run with or without sidechain CO inversion via the Inv_CO flag for
uniformly 13C-labeled proteins.

Reference
---------

Journal of Biomolecular NMR (2008) 42, 35-47

"""
import numpy as np
from numpy import linalg as la

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG2

_EXP_DETAILS = {
    "sidechain": {"type": str, "default": "False"},
    "taucc": {"type": float, "default": 9.09e-3},
}


class ProfileCPMGCOAP(ProfileCPMG2):
    """TODO: class docstring."""

    EXP_DETAILS = dict(**ProfileCPMG2.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixyzsz"
    CONSTRAINTS = "hn_ap"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.taucc = self.exp_details["taucc"]
        self.sidechain = self.get_bool(self.exp_details["sidechain"])

        # Set the row vector for detection
        self.detect = self.liouv.detect["2izsz_a"]

        # Set the delays in the experiments
        self.delays += [self.taucc]

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
        d_taucc = delays[self.taucc]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling

        # Calculate starting magnetization vector
        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")

        # Calculate the flip block
        if self.sidechain:
            p_flip = p180pmy
        else:
            p_flip = p90[3] @ d_taucc @ p180pmy @ d_taucc @ p90[1]

        # Calculating the cpmg trains
        cp = {0: self.liouv.identity}

        for ncyc in set(self.data["ncycs"][~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180[[1, 0]] @ tau_cp
            cp_train = la.matrix_power(echo, int(ncyc))
            cp[ncyc] = d_neg @ cp_train @ d_neg

        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ p90[1]
                @ cp[ncyc]
                @ p_flip
                @ cp[ncyc]
                @ p90[1]
                @ mag0
            )
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)
