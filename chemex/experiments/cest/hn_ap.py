"""
Pure anti-phase CEST
====================

Analyzes chemical exchange during the CEST block. This is calculated using the
12x12, single spin matrix:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b) ]

Reference
=========

Sekhar, Rosenzweig, Bouvignies and Kay. PNAS (2016) 113:E2794-E2801

"""
import numpy as np

from chemex.experiments.cest.base_cest import ProfileCEST

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t1": {"type": float},
    "b1_frq": {"type": float},
    "b1_inh": {"default": np.inf, "type": float},
    "b1_inh_res": {"default": 11, "type": int},
    "filter_offsets": {"default": 0.0, "type": float},
    "filter_bandwidths": {"default": 0.0, "type": float},
}


class ProfileCESTHNAP(ProfileCEST):
    """Profile for pure in-phase CEST."""

    SPIN_SYSTEM = "ixyzsz"
    CONSTRAINTS = "hn_ap"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        # Set the row vector for detection
        self.detect = self.liouv.detect["2izsz_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r1_i_a", "r2")):
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

        # As the CEST block is after t1 evolution, the excited state
        # magnetization is set to 0.
        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")
        mag0[6:] = 0.0

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            if not ref:
                cest = self.liouv.pulse_i(self.time_t1, 0.0, self.dephasing)
            else:
                cest = self.liouv.identity
            mag = self.liouv.collapse(self.detect @ cest @ mag0)
            profile.append(mag)

        return np.asarray(profile)
