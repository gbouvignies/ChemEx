"""15N - N-H anti-TROSY CPMG

Analyzes 15N constant-time anti-TROSY CPMG relaxation dispersion experiments
for measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the 12x12, two spin matrix:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b) ]

Note
----
Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:

Reference
---------
Proc Natl Acad Sci USA (2007) 104, 18473-7

"""

import numpy as np

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
    "taub": {"default": 2.68e-3, "type": float},
    "time_equil": {"default": 0.0, "type": float},
}


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, measurements, exp_details, model)

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

        self.detect = self.liouv.detect["2izsz_a"] + self.liouv.detect["iz_a"]

        self.t_cps = {
            ncyc: self.exp_details["time_t2"] / (4.0 * ncyc) - self.exp_details["pw90"]
            for ncyc in self.ncycs[self.ncycs > 0]
        }

        self.t_neg = -2.0 * self.exp_details["pw90"] / np.pi

        self.delays = [
            self.t_neg,
            self.exp_details["time_equil"],
            self.exp_details["taub"],
        ]
        self.delays.extend(self.t_cps.values())

        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
            nh_constraints=True,
        )

        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def calculate_unscaled_profile(self, **parvals):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex_ab : float
            Exchange rates between states A and B in /s.
        dw_i_ab : float
            Chemical shift difference between states A and B in rad/s.
        r1_i_a : float
            Longitudinal relaxation rate of states A in /s.
        r2_i_a : float
            Transverse relaxation rate of state A in /s.
        dr2_i_ab : float
            Transverse relaxation rate difference between states A and B in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        self.liouv.update(**parvals)

        # Calculation of all the needed propagators
        pulses = self.liouv.pulses_90_180_i()
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(term="2izsz", **parvals)

        profile = []

        p_element = (
            self.liouv.perfect180["sx"]
            @ delays[self.exp_details["taub"]]
            @ pulses["90mx"]
            @ pulses["90py"]
            @ self.liouv.perfect180["sx"]
            @ pulses["90py"]
            @ pulses["90mx"]
            @ delays[self.exp_details["taub"]]
        )

        p_element_pc = 0.5 * (
            pulses["90px"] @ p_element @ pulses["90py"]
            + pulses["90mx"] @ p_element @ pulses["90my"]
        )

        for ncyc in self.ncycs:

            if ncyc == 0:

                mag = (
                    self.detect
                    @ delays[self.exp_details["time_equil"]]
                    @ pulses["90py"]
                    @ p_element
                    @ pulses["90px"]
                    @ mag0
                )

            else:

                p_cpx = np.linalg.matrix_power(
                    delays[self.t_cps[ncyc]]
                    @ pulses["180px"]
                    @ delays[self.t_cps[ncyc]],
                    int(ncyc),
                )

                p_cpy = np.linalg.matrix_power(
                    delays[self.t_cps[ncyc]]
                    @ pulses["180py"]
                    @ delays[self.t_cps[ncyc]],
                    int(ncyc),
                )

                mag = (
                    self.detect
                    @ delays[self.exp_details["time_equil"]]
                    @ pulses["90py"]
                    @ delays[self.t_neg]
                    @ p_cpx
                    @ delays[self.t_neg]
                    @ p_element_pc
                    @ delays[self.t_neg]
                    @ p_cpy
                    @ delays[self.t_neg]
                    @ pulses["90px"]
                    @ mag0
                )

            profile.append(np.float64(self.liouv.collapse(mag)))

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""
        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.exp_details["time_t2"]
