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

Reference
---------

Journal of the American Chemical Society (2004), 126, 3964-73

"""
import numpy as np
from numpy import linalg as la

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG1

_EXP_DETAILS = {"small_protein": {"type": str, "default": "False"}}


class ProfileCPMGCH3MQ(ProfileCPMG1):
    """TODO: class docstring."""

    EXP_DETAILS = dict(**ProfileCPMG1.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixysxy"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.small_protein = self.get_bool(self.exp_details["small_protein"])

        # Set the delays in the experiments
        self.t_zeta = 1.0 / (8.0 * 125.3)
        self.delays = [self.t_zeta] + list(self.tau_cps.values())

        # Set the row vector for detection
        self.detect = self.liouv.detect["2iysx_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_mq_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: class docstring."""

        self.liouv.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_zeta = delays[self.t_zeta]

        # Calculation of the propagators corresponding to all the pulses
        p180_sx = self.liouv.perfect180["sx"]
        p180_ix = self.liouv.perfect180["ix"]
        p180_iy = self.liouv.perfect180["iy"]

        # Calculate starting magnetization vector
        mag0 = self.liouv.compute_mag_eq(params_local, term="2iysx")

        if self.small_protein:
            mag0 = d_zeta @ p180_sx @ p180_ix @ d_zeta @ mag0

        # Calculating the cpmg trains
        cp = {0: self.liouv.identity}

        for ncyc in set(self.data["ncycs"][~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180_iy @ tau_cp
            cp_train = la.matrix_power(echo, int(ncyc))
            cp[ncyc] = cp_train @ p180_sx @ cp_train

        profile = [
            self.liouv.collapse(self.detect @ cp[ncyc] @ mag0)
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)
