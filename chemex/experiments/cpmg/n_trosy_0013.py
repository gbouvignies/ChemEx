"""
15N - N-H TROSY CPMG with 0013 phase cycle
==========================================

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

This version is modified such that CPMG pulses are applied with [0013] phase
cycle.  Off resonance effects are taken into account.


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
  * ncyc_max     (maximum number of cycles)
  * antitrosy    (perform anti-trosy CPMG RD experiment [False])
  * s3e          (s3e trosy selection [True])


Extra parameters
----------------

  * path         (directory of the profiles)
  * error        (= 'file': uncertainties are taken from the profile files
                  = 'auto': uncertainties are calculated from duplicates)


"""
from functools import reduce

import numpy as np

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG2

_EXP_DETAILS = {
    "taub": {"default": 2.68e-3, "type": float},
    "ncyc_max": {"type": int},
    "antitrosy": {"default": "False", "type": str},
    "s3e": {"default": "True", "type": str},
}


class ProfileCPMGNTROSY0013(ProfileCPMG2):
    """TODO: class docstring."""

    EXP_DETAILS = dict(**ProfileCPMG2.EXP_DETAILS, **_EXP_DETAILS)
    SPIN_SYSTEM = "ixyzsz"
    CONSTRAINTS = "nh"
    CP_PHASES = [
        [1, 1, 0, 2, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 0, 2],
        [0, 2, 3, 3, 2, 0, 3, 3, 2, 0, 3, 3, 0, 2, 3, 3],
        [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
        [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
    ]

    def __init__(self, name, data, exp_details, model):

        super().__init__(name, data, exp_details, model)

        self.taub = self.exp_details["taub"]
        self.ncyc_max = self.exp_details["ncyc_max"]
        self.antitrosy = self.get_bool(self.exp_details["antitrosy"])
        self.s3e = self.get_bool(self.exp_details["s3e"])

        # Set the row vector for detection
        if self.antitrosy:
            self.detect = self.liouv.detect["2izsz_a"] + self.liouv.detect["iz_a"]
        else:
            self.detect = self.liouv.detect["2izsz_a"] - self.liouv.detect["iz_a"]

        # Set the delays in the experiments
        ncycs = self.data["ncycs"][~self.reference]
        self.tau_cps = dict(zip(ncycs, self.time_t2 / (4.0 * ncycs) - 0.75 * self.pw90))
        self.deltas = dict(zip(ncycs, 0.5 * self.pw90 * (self.ncyc_max - ncycs)))
        self.deltas[0] = 0.5 * self.pw90 * (self.ncyc_max - 1)
        self.delays = (
            [self.t_neg, self.time_eq, self.taub]
            + list(self.tau_cps.values())
            + list(self.deltas.values())
        )

        # Set the phase cycling of the cpmg pulses
        self.phases = self._get_phases()

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: Write docstring"""

        self.liouv.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_taub = delays[self.taub]
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        delta = {ncyc: delays[delay] for ncyc, delay in self.deltas.items()}

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180_sx = self.liouv.perfect180["sx"]

        # Getting the starting magnetization
        mag0 = self._get_mag0(params_local)

        # Calculating the p-element
        if self.antitrosy:
            palmer_ = (
                d_taub @ p90[1] @ p90[0] @ p180_sx @ p90[0] @ p90[1] @ d_taub @ p180_sx
            )
        else:
            palmer_ = (
                d_taub @ p90[0] @ p90[1] @ p180_sx @ p90[1] @ p90[0] @ d_taub @ p180_sx
            )

        palmer = np.mean(p90[[0, 2]] @ palmer_ @ p90[[1, 3]], axis=0)

        # Calculating the cpmg trains
        cp1 = {0: p180[[0, 1]]}
        cp2 = {0: p180[[1, 0]]}

        for ncyc in set(self.data["ncycs"][~self.reference]):

            tau_cp = delays[self.tau_cps[ncyc]]
            phase_cp = self.phases[ncyc]

            echo = tau_cp @ p180 @ tau_cp

            cp_trains = [reduce(np.matmul, echo[phases]) for phases in phase_cp]
            cp_trains = d_neg @ np.array(cp_trains) @ d_neg

            cp1[ncyc] = cp_trains[[0, 1]]
            cp2[ncyc] = cp_trains[[2, 3]]

        # Make profile
        profile = [
            self.liouv.collapse(
                self.detect
                @ d_eq
                @ delta[ncyc]
                @ p90[3]
                @ cp2[ncyc]
                @ palmer
                @ cp1[ncyc]
                @ p90[0]
                @ delta[ncyc]
                @ mag0
            )
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)

    def _get_mag0(self, params_local):

        mag0 = self.liouv.compute_mag_eq(params_local, term="2izsz")

        if self.s3e:
            if self.antitrosy:
                mag0 += self.liouv.compute_mag_eq(params_local, term="iz")
            else:
                mag0 -= self.liouv.compute_mag_eq(params_local, term="iz")
            mag0 *= 0.5

        return mag0

    def _get_phases(self):

        phases = {}

        for ncyc in self.data["ncycs"][~self.reference]:
            phases[ncyc] = np.take(self.CP_PHASES, np.arange(ncyc), mode="wrap", axis=1)
            phases[ncyc][:2] = np.flip(phases[ncyc][:2], axis=1)

        return phases
