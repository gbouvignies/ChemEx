"""
15N - N-H TROSY CPMG
====================

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)x(6n), two spin matrix, where n is the number of
states:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b),
  ... ]


Note
----
Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:


Reference
---------
Proc Natl Acad Sci USA (2007) 104, 18473-7


Experimental parameters
-----------------------
  * h_larmor_frq (1H Larmor frequency, in MHz)
  * temperature  (sample temperature, in Celsius)
  * p_total      (optional: protein concentration, in M)
  * l_total      (optional: ligand concentration, in M)
  * time_t2      (CPMG relaxation delay in seconds)
  * carrier      (position of the 15N carrier during the CPMG period, in ppm)
  * time_t2      (CPMG relaxation delay in seconds)
  * pw90         (15N 90 degree pulse width of CPMG pulses, in seconds)
  * time_equil   (equilibration delay at the end of the CPMG period, in seconds [0.0])
  * taub         (p-element delay = 1/4J, in seconds [2.68e-3])
  * antitrosy    (perform anti-trosy CPMG RD experiment [False])


Extra parameters
----------------
  * path         (directory of the profiles)
  * error        (= 'file': uncertainties are taken from the profile files
                  = 'auto': uncertainties are calculated from duplicates)


"""

import numpy as np
from numpy import linalg as la

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
    "taub": {"default": 2.68e-3, "type": float},
    "time_equil": {"default": 0.0, "type": float},
    "antitrosy": {"default": "False", "type": str},
}


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, measurements, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        self.time_t2 = self.exp_details["time_t2"]
        self.time_eq = self.exp_details["time_equil"]
        self.pw90 = self.exp_details["pw90"]
        self.taub = self.exp_details["taub"]
        self.antitrosy = self.get_bool(self.exp_details["antitrosy"])

        # Set the liouvillian
        self.liouv = basis.Liouvillian(
            system="ixyzsz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.carrier_i = self.exp_details["carrier"]
        self.liouv.w1_i = 2.0 * np.pi / (4.0 * self.pw90)

        # Set the row vector for detection
        if self.antitrosy:
            self.detect = self.liouv.detect["2izsz_a"] + self.liouv.detect["iz_a"]
        else:
            self.detect = self.liouv.detect["2izsz_a"] - self.liouv.detect["iz_a"]

        # Set the delays in the experiments
        ncycs = self.ncycs[~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs) - self.pw90))
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.delays = [self.t_neg, self.time_eq, self.taub] + list(
            self.tau_cps.values()
        )

        # Get the parameters this profile depends on
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
        """TODO: Write docstring"""

        self.liouv.update(**parvals)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_taub = delays[self.taub]
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180_sx = self.liouv.perfect180["sx"]

        # Getting the starting magnetization
        mag0 = self.liouv.compute_mag_eq(term="2izsz", **parvals)

        # Calculating the p-element
        if self.antitrosy:
            palmer_ = (
                p180_sx @ d_taub @ p90[2] @ p90[1] @ p180_sx @ p90[1] @ p90[2] @ d_taub
            )
        else:
            palmer_ = (
                p180_sx @ d_taub @ p90[1] @ p90[0] @ p180_sx @ p90[0] @ p90[1] @ d_taub
            )

        palmer = np.mean(p90[[0, 2]] @ palmer_ @ p90[[1, 3]], axis=0)

        # Calculating the cpmg trains
        cp1 = {0: np.identity(p180[1].shape[-1])}
        cp2 = {0: np.identity(p180[0].shape[-1])}

        for ncyc in set(self.ncycs[~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180[[1, 0]] @ tau_cp
            cp_trains = la.matrix_power(echo, int(ncyc))
            cp1[ncyc], cp2[ncyc] = d_neg @ cp_trains @ d_neg

        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ p90[1]
                @ cp2[ncyc]
                @ palmer
                @ cp1[ncyc]
                @ p90[0]
                @ mag0
            )
            for ncyc in self.ncycs
        ]

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""
        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.exp_details["time_t2"]
