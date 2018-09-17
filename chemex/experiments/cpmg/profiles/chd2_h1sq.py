"""1H(methyl - 13CHD2) - Pure Anti-Phase Proton CPMG

Measures methyl proton chemical exchange recorded on site-specifically
13CHD2-labeled proteins in a highly deuterated background. Magnetization is
initally anti-phase and is read out as anti-phase prior to 13C evolution. The
calculation uses a 12x12 basis set:

[ Hx(a), Hy(a), Hz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
  Hx(b), Hy(b), Hz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

Note
----

Off resonance effects are taken into account. The calculation is designed
explicitly for analyzing the Lewis Kay pulse sequence:

CHD2_H_SQ_exchange_hsqc_lek_*00_enhanced

Reference
---------

Journal of the American Chemical Society (2010) 132, 10992-5
"""

import numpy as np

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
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

        self.detect = self.liouv.detect["2izsz_a"]

        self.t_cps = {
            ncyc: self.exp_details["time_t2"] / (4.0 * ncyc) - self.exp_details["pw90"]
            for ncyc in self.ncycs[self.ncycs > 0]
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

        p_180pmx = 0.5 * (pulses["180px"] + pulses["180mx"])  # +/- phase cycling

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(term="2izsz", **parvals)

        profile = []

        for ncyc in self.ncycs:

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
            ncycs = self.ncycs

        return ncycs / self.exp_details["time_t2"]
