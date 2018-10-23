"""15N CEST with CW decoupling.

Analyzes chemical exchange in the presence of 1H CW decoupling during the CEST
block. Magnetization evolution is calculated using the 30*30 two-spin matrix:

[ I{xyz}, S{xyz}, 2I{xyz}S{xyz} ]{a, b, ...}

Notes
-----
The calculation is designed specifically to analyze the experiment found in the
reference:

Bouvignies and Kay. J Phys Chem B (2012), 116:14311-7

"""
import numpy as np

from chemex.experiments.cest.base_cest import ProfileCEST

_EXP_DETAILS = {"carrier_dec": {"type": float}, "b1_frq_dec": {"type": float}}


class ProfileCESTNIPHCW(ProfileCEST):
    """Profile for CEST with CW decoupling."""

    EXP_DETAILS = dict(**ProfileCEST.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixyzsxyz"
    CONSTRAINTS = "nh"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.liouv.w1_s = 2 * np.pi * self.exp_details["b1_frq_dec"]
        self.liouv.carrier_s = self.exp_details["carrier_dec"]

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(
                ("dw", "r1_i_a", "r2_i_a", "r2_mq_a", "etaxy_i_a", "etaz_i_a")
            ):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, offsets=None):
        """Calculate the CEST profile in the presence of exchange.

        TODO: Parameters
        ----------

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        self.liouv.update(params_local)

        reference = self.reference
        carriers_i = self.carriers_i

        if offsets is not None:
            reference = np.zeros_like(offsets, dtype=np.bool)
            carriers_i = self.offsets_to_ppm(offsets)

        mag0 = self.liouv.compute_mag_eq(params_local, term="iz")

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            if not ref:
                cest = self.liouv.pulse_is(self.time_t1, 0.0, 0.0, self.dephasing)
            else:
                cest = self.liouv.identity
            mag = self.liouv.collapse(self.detect @ cest @ mag0)
            profile.append(mag)

        return np.asarray(profile)
