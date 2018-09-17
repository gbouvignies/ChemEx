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

from chemex.experiments.cest import cest_profile
from chemex.spindynamics import basis, constants, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t1": {"type": float},
    "pw90_dante": {"type": float},
    "sw_dante": {"type": float},
    "b1_frq": {"type": float},
    "b1_inh": {"default": np.inf, "type": float},
    "b1_inh_res": {"default": 11, "type": int},
    "cn_label": {"default": False, "type": bool},
    "filter_offsets": {"default": 0.0, "type": float},
    "filter_bandwidths": {"default": 0.0, "type": float},
}


class Profile(cest_profile.CESTProfile):
    """Profile for pure in-phase CEST."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, measurements, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        pw90_dante = self.exp_details["pw90_dante"]
        self.sw_dante = self.exp_details["sw_dante"]
        b1_frq = self.exp_details["b1_frq"]
        time_t1 = self.exp_details["time_t1"]

        self.pw_dante = 4.0 * pw90_dante * b1_frq / self.sw_dante
        self.tau_dante = 1.0 / self.sw_dante - self.pw_dante
        self.ncyc_dante = int(time_t1 * self.sw_dante + 0.1)

        self.liouv = basis.Liouvillian(
            system="ixyz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.w1_i = 2.0 * np.pi / (4.0 * pw90_dante)
        self.liouv.w1_i_inh = self.exp_details["b1_inh"]
        self.liouv.w1_i_inh_res = self.exp_details["b1_inh_res"]

        self.carriers_i = self.b1_offsets_to_ppm()
        self.detect = self.liouv.detect["iz_a"]

        if self.exp_details["cn_label"]:
            symbol, nucleus = self.peak.symbols["i"], self.peak.nuclei["i"]
            j_values, j_weights = constants.get_multiplet(symbol, nucleus)
            self.liouv.j_eff_i = j_values
            self.liouv.j_eff_i_weights = j_weights

        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
        )

        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r1_i_a", "r2")):
                self.params[full_name].set(vary=True)

    def calculate_unscaled_profile(self, b1_offsets=None, **parvals):
        """Calculate the CEST profile in the presence of exchange.

        TODO: Parameters
        ----------

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        reference = self.reference
        carriers_i = self.carriers_i

        if b1_offsets is not None:
            carriers_i = (
                self.exp_details["carrier"]
                + 2.0 * np.pi * b1_offsets / self.liouv.ppms["i"]
            )
            reference = [False for _ in b1_offsets]

        mag0 = self.liouv.compute_mag_eq(term="iz", **parvals)

        self.liouv.update(**parvals)

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            mag = mag0.copy()

            if not ref:
                p_delay = self.liouv.delays(self.tau_dante)
                p_pulse = self.liouv.pulse_i(self.pw_dante, 0.0)
                mag = np.linalg.matrix_power(p_pulse @ p_delay, self.ncyc_dante) @ mag

            mag = self.detect @ mag
            profile.append(np.float64(self.liouv.collapse(mag)))

        return np.asarray(profile)

    def filter_points(self, params=None):
        """Evaluate some criteria to know whether or not the point should be
        considered in the calculation."""

        cs = params[self.map_names["cs_i_a"]].value

        filter_offsets = np.asarray(self.exp_details["filter_offsets"]).reshape(-1)
        filter_bandwidths = np.asarray(self.exp_details["filter_bandwidths"]).reshape(
            -1
        )

        for offset, bandwidth in zip(filter_offsets, filter_bandwidths):
            nu_offsets = (
                (cs - self.exp_details["carrier"])
                * self.liouv.ppms["i"]
                / (2.0 * np.pi)
                - self.b1_offsets
                + offset
            )
            nu_offsets = (
                nu_offsets + 0.5 * self.sw_dante
            ) % self.sw_dante - 0.5 * self.sw_dante

            self.mask = np.logical_and(self.mask, abs(nu_offsets) > bandwidth * 0.5)

    def b1_offsets_to_ppm(self, b1_offsets=None):
        """Convert B1 offset from Hz to ppm."""
        if b1_offsets is None:
            b1_offsets = self.b1_offsets

        return (
            2.0 * np.pi * b1_offsets / self.liouv.ppms["i"]
            + self.exp_details["carrier"]
        )
