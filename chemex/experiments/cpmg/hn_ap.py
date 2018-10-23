"""1H - Pure Anti-phase Proton CPMG

Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 12x12,
2-spin exchange matrix:

[ Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
  Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b)]

Note
----
Off resonance effects are taken into account. The calculation is designed
explicitly for analyzing the Lewis Kay pulse sequence:

H1_CPMG_Rex_hsqc_lek_x00

with antiphase_flg set to 'y'

Journal of Biomolecular NMR (2011) 50, 13-8
"""
import numpy as np
from numpy import linalg as la

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG2


class ProfileCPMGHNAP(ProfileCPMG2):
    """TODO: class docstring."""

    SPIN_SYSTEM = "ixyzsz"
    CONSTRAINTS = "hn_ap"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        # Set the row vector for detection
        self.detect = self.liouv.detect["2izsz_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: class docstring."""

        self.liouv.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Calculate starting magnetization vector
        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")

        # Calculating the cpmg trains
        cp1 = {0: self.liouv.identity}
        cp2 = {0: self.liouv.identity}

        for ncyc in set(self.data["ncycs"][~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180[1] @ tau_cp
            cp_train = la.matrix_power(echo, int(ncyc))
            cp1[ncyc] = cp_train @ d_neg
            cp2[ncyc] = d_neg @ cp_train

        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ p90[0]
                @ cp2[ncyc]
                @ p180pmx
                @ cp1[ncyc]
                @ p90[0]
                @ mag0
            )
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)
