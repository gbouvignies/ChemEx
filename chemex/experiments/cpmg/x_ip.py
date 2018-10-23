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

from chemex.experiments.cpmg.base_cpmg import ProfileCPMG2


class ProfileCPMGXIP(ProfileCPMG2):
    """TODO: class docstring."""

    SPIN_SYSTEM = "ixyz"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r2_i_a")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """TODO: Write docstring"""

        self.liouv.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        delays = dict(zip(self.delays, self.liouv.delays(self.delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]

        # Calculation of the propagators corresponding to all the pulses
        pulses = self.liouv.pulses_90_180_i()
        p90 = np.array([pulses[name] for name in ["90px", "90py", "90mx", "90my"]])
        p180 = np.array([pulses[name] for name in ["180px", "180py", "180mx", "180my"]])
        p180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Calculate starting magnetization vector
        mag0 = self.liouv.compute_mag_eq(params_local, term="iz")

        # Calculating the cpmg trains
        cp1 = {0: self.liouv.identity, -1: delays[self.tau_cps[-1]] @ d_neg}
        cp2 = {0: self.liouv.identity, -1: d_neg @ delays[self.tau_cps[-1]]}

        for ncyc in set(self.data["ncycs"][~self.reference]):
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
                @ p180pmx
                @ cp1[ncyc]
                @ p90[0]
                @ mag0
            )
            for ncyc in self.data["ncycs"]
        ]

        return np.asarray(profile)
