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

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG
from chemex.spindynamics import basis
from chemex.spindynamics import default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
    "time_equil": {"default": 0.0, "type": float},
}


class ProfileCPMG_HN_AP(ProfileCPMG):
    """TODO: class docstring."""

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        self.liouv = basis.Liouvillian(
            system="ixyzsz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.carrier_i = self.exp_details["carrier"]
        self.liouv.w1_i = 2.0 * np.pi / (4.0 * self.exp_details["pw90"])

        self.detect = self.liouv.detect["2izsz_a"]

        self.t_cps = {
            ncyc: self.exp_details["time_t2"] / (4.0 * ncyc) - self.exp_details["pw90"]
            for ncyc in self.data["ncycs"][~self.reference]
        }

        self.t_neg = -2.0 * self.exp_details["pw90"] / np.pi

        self.delays = [self.t_neg, self.exp_details["time_equil"]]
        self.delays.extend(self.t_cps.values())

        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
            hn_ap_constraints=True,
        )

        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: class docstring."""

        self.liouv.update(params_local)

        # Calculation of all the needed propagators
        pulses = self.liouv.pulses_90_180_i()
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))

        p_180pmx = 0.5 * (pulses["180px"] + pulses["180mx"])  # +/- phase cycling

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")

        profile = []

        for ncyc in self.data["ncycs"]:

            if ncyc == 0:

                mag = (
                    self.detect
                    @ delays[self.exp_details["time_equil"]]
                    @ pulses["90px"]
                    @ p_180pmx
                    @ pulses["90px"]
                    @ mag0
                )

            else:
                p_cp = np.linalg.matrix_power(
                    delays[self.t_cps[ncyc]]
                    @ pulses["180py"]
                    @ delays[self.t_cps[ncyc]],
                    int(ncyc),
                )
                mag = (
                    self.detect
                    @ delays[self.exp_details["time_equil"]]
                    @ pulses["90px"]
                    @ delays[self.t_neg]
                    @ p_cp
                    @ p_180pmx
                    @ p_cp
                    @ delays[self.t_neg]
                    @ pulses["90px"]
                    @ mag0
                )

            profile.append(np.float64(self.liouv.collapse(mag)))

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""
        if ncycs is None:
            ncycs = self.data["ncycs"]

        return ncycs / self.exp_details["time_t2"]
