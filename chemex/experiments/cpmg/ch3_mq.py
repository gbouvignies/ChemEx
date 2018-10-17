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

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG
from chemex.spindynamics import basis
from chemex.spindynamics import default

EXP_DETAILS = {
    "time_t2": {"type": float},
    "small_protein_flg": {"type": str, "default": "n"},
}


class ProfileCPMG_CH3_MQ(ProfileCPMG):
    """TODO: class docstring."""

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        self.liouv = basis.Liouvillian(
            system="ixysxy",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.t_zeta = 1.0 / (8.0 * 125.3)
        self.t_cps = {
            ncyc: self.exp_details["time_t2"] / (4.0 * ncyc)
            for ncyc in self.data["ncycs"][~self.reference]
        }

        self.delays = [self.t_zeta]
        self.delays.extend(self.t_cps.values())

        self.detect = self.liouv.detect["2iysx_a"]

        self.small_protein_flg = (
            self.exp_details["small_protein_flg"].lower().startswith("y")
        )

        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
        )

        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_mq_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: class docstring."""

        self.liouv.update(params_local)

        # Calculation of all the needed propagators
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(params_local, term="2iysx")

        profile = []

        if self.small_protein_flg:
            mag0 = (
                delays[self.t_zeta]
                @ self.liouv.perfect180["sx"]
                @ self.liouv.perfect180["ix"]
                @ delays[self.t_zeta]
                @ mag0
            )

        for ncyc in self.data["ncycs"]:

            if ncyc == 0:

                mag = self.detect @ mag0

            else:

                p_cpy = np.linalg.matrix_power(
                    delays[self.t_cps[ncyc]]
                    @ self.liouv.perfect180["iy"]
                    @ delays[self.t_cps[ncyc]],
                    int(ncyc),
                )

                mag = self.detect @ p_cpy @ self.liouv.perfect180["sx"] @ p_cpy @ mag0

            profile.append(np.float64(self.liouv.collapse(mag)))

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""

        if ncycs is None:
            ncycs = self.data["ncycs"]

        return ncycs / self.exp_details["time_t2"]
