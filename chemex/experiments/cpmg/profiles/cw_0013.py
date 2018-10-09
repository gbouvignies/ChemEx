"""
15N - Pure In-phase Nitrogen CPMG with 0013 phase cycle
=======================================================

Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)x(3n), single spin matrix, where n is the number
of states:

[ Ix(a), Iy(a), Iz(a),
  Ix(b), Iy(b), Iz(b),
   ... ]

This version is modified such that CPMG pulses are applied with [0013] phase cycle
as shown in Daiwen's paper. The cw decoupling on 1H is assumed to be
strong enough (> 15 kHz) such that perfect 1H decoupling can be achieved.

References
----------

Jiang et al. Journal of Magnetic Resonance (2015) 257
Hansen, Vallurupalli and Kay. Journal of Physical Chemistry B (2008) 112:5898-5904

Experimental parameters
-----------------------
  * h_larmor_frq (1H Larmor frequency, in MHz)
  * temperature  (sample temperature, in Celsius)
  * p_total      (optional: protein concentration, in M)
  * l_total      (optional: ligand concentration, in M)
  * time_t2      (CPMG relaxation delay in seconds)
  * carrier      (position of the 15N carrier during the CPMG period, in ppm)
  * pw90         (15N 90 degree pulse width of CPMG pulses, in seconds)
  * time_equil   (equilibration delay at the end of the CPMG period, in seconds)
  * ncyc_max     (maximum number of cycles)


Extra parameters
----------------
  * path         (directory of the profiles)
  * error        (= 'file': uncertainties are taken from the profile files
                  = 'auto': uncertainties are calculated from duplicates)

"""

from functools import reduce

import numpy as np

from chemex.experiments.cpmg import cpmg_profile
from chemex.spindynamics import basis, default

EXP_DETAILS = {
    "carrier": {"type": float},
    "time_t2": {"type": float},
    "pw90": {"type": float},
    "time_equil": {"default": 0.0, "type": float},
    "ncyc_max": {"type": int},
}

CP_PHASES = [
    [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
    [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
]


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, measurements, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        self.time_t2 = self.exp_details["time_t2"]
        self.time_eq = self.exp_details["time_equil"]
        self.pw90 = self.exp_details["pw90"]
        self.ncyc_max = self.exp_details["ncyc_max"]

        # Set the liouvillian
        self.liouv = basis.Liouvillian(
            system="ixyz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.carrier_i = self.exp_details["carrier"]
        self.liouv.w1_i = 2.0 * np.pi / (4.0 * self.pw90)

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the delays in the experiments
        ncycs = self.ncycs[~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs) - 0.75 * self.pw90))
        self.deltas = dict(zip(ncycs, self.pw90 * (self.ncyc_max - ncycs)))
        self.deltas[0] = self.pw90 * (self.ncyc_max - 1)
        self.t_pos2 = +4.0 * self.pw90 / np.pi
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.delays = (
            [self.t_neg, self.t_pos2, self.time_eq]
            + list(self.tau_cps.values())
            + list(self.deltas.values())
        )

        # Set the phase cycling of the cpmg pulses
        self.phases = {
            ncyc: np.take(CP_PHASES, np.flip(np.arange(2 * ncyc)), mode="wrap", axis=1)
            for ncyc in ncycs
        }

        # Get the parameters this profile depends on
        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
        )

        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def calculate_unscaled_profile(self, **parvals):
        """TODO: Write docstring"""

        self.liouv.update(**parvals)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_neg = delays[self.t_neg]
        d_pos2 = delays[self.t_pos2]
        d_eq = delays[self.time_eq]
        delta = {ncyc: delays[delay] for ncyc, delay in self.deltas.items()}

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])

        # Getting the starting magnetization
        mag0 = self.liouv.compute_mag_eq(term="iz", **parvals)

        # Calculating the cpmg trains
        cp = {0: p180[[0, 3]] @ d_pos2 @ p180[[0, 1]]}

        for ncyc in set(self.ncycs[~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            phase_cp = self.phases[ncyc]

            echo = tau_cp @ p180 @ tau_cp

            cp[ncyc] = [reduce(np.matmul, echo[phases]) for phases in phase_cp]
            cp[ncyc] = d_neg @ np.array(cp[ncyc]) @ d_neg

        # Make profile
        profile = [
            self.liouv.collapse(
                self.detect @ d_eq @ delta[ncyc] @ p90[3] @ cp[ncyc] @ p90[1] @ mag0
            )
            for ncyc in self.ncycs
        ]

        return np.asarray(profile)

    def ncycs_to_nu_cpmgs(self, ncycs=None):
        """Calculate the pulsing frequency, v(CPMG), from ncyc values."""
        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.exp_details["time_t2"]
