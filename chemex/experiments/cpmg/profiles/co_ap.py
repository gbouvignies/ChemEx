"""13CO - Pure Anti-phase Carbonyl 13C CPMG

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 12x12,
2-spin exchange matrix:

[ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
  COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b)]

Notes
-----

Because of the length of the shaped pulses used during the CPMG blocks, off-
resonance effects are taken into account only for the 90-degree pulses that
create COxNz before the CPMG and COzNz after the CPMG. The calculation is
designed explicitly for analyzing the Kay laboratory pulse sequence:

CO_CPMG_SCFilter_x00_dfh1

And can be run with or without sidechain CO inversion via the Inv_CO flag for
uniformly 13C-labeled proteins.

Reference
---------

Journal of Biomolecular NMR (2008) 42, 35-47

"""

import numpy as np

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
    "time_equil": {"type": float, "default": 0.0},
    "sidechain_flg": {"type": str, "default": "n"},
    "taucc": {"type": float, "default": 9.09e-3},
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
            ncyc: self.exp_details["time_t2"] / (4.0 * ncyc)
            for ncyc in self.ncycs[self.ncycs > 0]
        }

        self.t_neg = -2.0 * self.exp_details["pw90"] / np.pi

        self.delays = [
            self.t_neg,
            self.exp_details["time_equil"],
            self.exp_details["taucc"],
        ]
        self.delays.extend(self.t_cps.values())

        self.sidechain_flg = self.exp_details["sidechain_flg"].lower().startswith("y")

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

        p_180pmy = 0.5 * (pulses["180py"] + pulses["180my"])  # +/- phase cycling

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(term="2izsz", **parvals)

        profile = []

        if self.sidechain_flg:
            p_flip = p_180pmy
        else:
            p_flip = (
                pulses["90my"]
                @ delays[self.exp_details["taucc"]]
                @ p_180pmy
                @ delays[self.exp_details["taucc"]]
                @ pulses["90py"]
            )

        for ncyc in self.ncycs:

            if ncyc == 0:

                mag = (
                    self.detect
                    @ delays[self.exp_details["time_equil"]]
                    @ pulses["90py"]
                    @ p_flip
                    @ pulses["90py"]
                    @ mag0
                )

            else:

                p_cpx = np.linalg.matrix_power(
                    delays[self.t_cps[ncyc]]
                    @ self.liouv.perfect180["ix"]
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
                    @ p_flip
                    @ delays[self.t_neg]
                    @ p_cpx
                    @ delays[self.t_neg]
                    @ pulses["90py"]
                    @ mag0
                )

            profile.append(np.float64(self.liouv.collapse(mag)))

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""

        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.exp_details["time_t2"]
