"""
Pure In-phase D-CEST
====================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]


Reference
---------

Yuwen, Kay and Bouvignies. ChemPhysChem (2018) 19:1707-1710

"""
import numpy as np

from chemex.experiments.cest.base_cest import ProfileCEST

_EXP_DETAILS = {"pw90_dante": {"type": float}, "sw_dante": {"type": float}}


class ProfileCESTXIPDANTE(ProfileCEST):
    """Profile for pure in-phase CEST."""

    EXP_DETAILS = dict(**ProfileCEST.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixyz"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.sw_dante = self.exp_details["sw_dante"]
        pw90_dante = self.exp_details["pw90_dante"]
        b1_frq = self.exp_details["b1_frq"]

        self.pw_dante = 4.0 * pw90_dante * b1_frq / self.sw_dante
        self.tau_dante = 1.0 / self.sw_dante - self.pw_dante
        self.ncyc_dante = int(self.time_t1 * self.sw_dante + 0.1)

        # Set the liouvillian
        self.liouv.w1_i = 2.0 * np.pi / (4.0 * pw90_dante)

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

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

        mag0 = self.liouv.compute_mag_eq(params_local, term="iz")

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            if not ref:
                p_delay = self.liouv.delays(self.tau_dante)
                p_pulse = self.liouv.pulse_i(self.pw_dante, 0.0)
                dcest = np.linalg.matrix_power(p_pulse @ p_delay, self.ncyc_dante)
            else:
                dcest = self.liouv.identity
            mag = self.liouv.collapse(self.detect @ dcest @ mag0)
            profile.append(mag)

        return np.asarray(profile)

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether or not the point should be
        considered in the calculation."""

        cs_a = params[self.map_names["cs_i_a"]].value

        filter_offsets = np.asarray(self.exp_details["filter_offsets"]).reshape(-1)
        filter_bandwidths = np.asarray(self.exp_details["filter_bandwidths"]).reshape(
            -1
        )

        for offset, bandwidth in zip(filter_offsets, filter_bandwidths):
            nu_offsets = (
                (cs_a - self.carrier) * self.ppms_i / (2.0 * np.pi)
                - self.data["offsets"]
                + offset
            )
            nu_offsets = (
                nu_offsets + 0.5 * self.sw_dante
            ) % self.sw_dante - 0.5 * self.sw_dante

            self.mask = np.logical_and(self.mask, abs(nu_offsets) > bandwidth * 0.5)
