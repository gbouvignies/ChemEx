"""
Pure In-phase CEST
==================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]


Reference
---------

Vallurupalli, Bouvignies and Kay. J Am Chem Soc (2012) 134:8148-8161

"""

import numpy as np

from chemex.experiments.cest import cest_profile
from chemex.spindynamics import basis, constants, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t1": {"type": float},
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

        self.time_t1 = self.exp_details["time_t1"]

        self.liouv = basis.Liouvillian(
            system="ixyz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.w1_i = 2 * np.pi * self.exp_details["b1_frq"]
        self.liouv.w1_i_inh = self.exp_details["b1_inh"]
        self.liouv.w1_i_inh_res = self.exp_details["b1_inh_res"]

        self.carriers_i = self.b1_offsets_to_ppm()
        self.detect = self.liouv.detect["iz_a"]
        self.dephasing = self.exp_details["b1_inh"] == np.inf

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

        self.liouv.update(**parvals)

        reference = self.reference
        carriers_i = self.carriers_i

        if b1_offsets is not None:
            reference = np.zeros_like(b1_offsets, dtype=np.bool)
            carriers_i = self.b1_offsets_to_ppm(b1_offsets)

        mag0 = self.liouv.compute_mag_eq(term="iz", **parvals)

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            mag = mag0.copy()
            if not ref:
                mag = self.liouv.pulse_i(self.time_t1, 0.0, self.dephasing) @ mag
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

            self.mask = np.logical_and(self.mask, abs(nu_offsets) > bandwidth * 0.5)

    def b1_offsets_to_ppm(self, b1_offsets=None):
        """Convert B1 offset from Hz to ppm."""

        if b1_offsets is None:
            b1_offsets = self.b1_offsets

        return (
            2.0 * np.pi * b1_offsets / self.liouv.ppms["i"]
            + self.exp_details["carrier"]
        )
