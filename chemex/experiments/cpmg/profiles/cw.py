"""
Pure in-phase CPMG
==================

Analyzes chemical exchange in the presence of high power 1H CW decoupling during
the CPMG block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)x(3n), single spin matrix, where n is the number of
states:

[ Ix(a), Iy(a), Iz(a),
  Ix(b), Iy(b), Iz(b),
   ... ]


Notes
-----
Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:


Reference
---------
Journal of Physical Chemistry B (2008), 112, 5898-5904


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
    "time_equil": {"default": 0.0, "type": float},
}


class Profile(cpmg_profile.CPMGProfile):
    """TODO: class docstring."""

    def __init__(self, name, measurements, exp_details, model):
        super().__init__(name, measurements, exp_details, model)

        self.exp_details = self.check_exp_details(exp_details, expected=EXP_DETAILS)

        self.time_t2 = self.exp_details["time_t2"]
        self.time_eq = self.exp_details["time_equil"]
        self.pw90 = self.exp_details["pw90"]

        # Set the liouvillian
        self.liouv = basis.Liouvillian(
            system="ixyz",
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=False,
        )

        self.liouv.carrier_i = self.exp_details["carrier"]
        self.liouv.w1_i = 2.0 * np.pi / (4.0 * self.exp_details["pw90"])

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the delays in the experiments
        ncycs = self.ncycs[~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs) - self.pw90))
        self.tau_cps[-1.0] = 0.5 * self.exp_details["time_t2"]
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.delays = [self.t_neg, self.time_eq] + list(self.tau_cps.values())

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
        d_eq = delays[self.time_eq]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p_180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Simulate the CPMG block as function of ncyc
        mag0 = self.liouv.compute_mag_eq(term="iz", **parvals)

        # Calculating the cpmg trains
        cp1 = {0: np.identity(p180[1].shape[-1]), -1: delays[self.tau_cps[-1]] @ d_neg}
        cp2 = {0: np.identity(p180[0].shape[-1]), -1: d_neg @ delays[self.tau_cps[-1]]}

        for ncyc in set(self.ncycs[~self.reference]):
            tau_cp = delays[self.tau_cps[ncyc]]
            echo = tau_cp @ p180[1] @ tau_cp
            cp_trains = la.matrix_power(echo, int(ncyc))
            cp1[ncyc] = cp_trains @ d_neg
            cp2[ncyc] = d_neg @ cp_trains

        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ p90[0]
                @ cp2[ncyc]
                @ p_180pmx
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
            ncycs = np.array([ncyc if ncyc >= 0 else 0.5 for ncyc in self.ncycs])

        return ncycs / self.exp_details["time_t2"]
