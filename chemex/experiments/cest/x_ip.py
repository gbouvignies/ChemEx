"""
Pure In-phase CEST
==================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)x(3n), single spin matrix, where n is the number of
states:

[ Ix(a), Iy(a), Iz(a),
  Ix(b), Iy(b), Iz(b),
   ... ]


Reference
---------
Vallurupalli, Bouvignies and Kay. J Am Chem Soc (2012) 134:8148-8161


Experimental parameters
-----------------------
  * h_larmor_frq (1H Larmor frequency, in MHz)
  * temperature  (sample temperature, in Celsius)
  * p_total      (optional: protein concentration, in M)
  * l_total      (optional: ligand concentration, in M)
  * time_t1      (CEST relaxation delay, in seconds)
  * carrier      (position of the carrier during the CEST period, in ppm)
  * b1_frq       (B1 radio-frequency field strength, in Hz)
  * b1_inh       (B1 inhomogeneity expressed as a fraction of 'b1_inh'.
                  If not set, a faster calculation takes place assuming
                  full dephasing of the magnetization components that oscillate
                  during the irradiation period.)
  * b1_inh_res   (number of points used to simulate B1 inhomogeneity, the larger
                  the longer the calculation)


Extra parameters
----------------
  * path              (directory of the profiles)
  * error             (= 'file': uncertainties are taken from the profile files
                       = 'auto': uncertainties are calculated from the baseline)
  * filter_offsets    (list of offsets relative to the main resonance position,
                       defining regions where points are excluded from the
                       calculation, in Hz)
  * filter_bandwidths (list of values defining the exclusion range around each
                       previously defined offset, in Hz)

"""
import numpy as np

from chemex.experiments.cest.base_cest import ProfileCEST


class ProfileCESTXIP(ProfileCEST):
    """Profile for pure in-phase CEST."""

    SPIN_SYSTEM = "ixyz"

    def __init__(self, name, data, exp_details, model):
        super().__init__(name, data, exp_details, model)

        # Set the row vector for detection
        self.detect = self.liouv.detect["iz_a"]

        # Set the varying parameters by default
        for name, full_name in self.map_names.items():
            if name.startswith(("dw", "r1_i_a", "r2")):
                self.params[full_name].set(vary=True)

    def _calculate_unscaled_profile(self, params_local, offsets=None):
        """Calculate the CEST profile in the presence of exchange.

        TODO: Parameters
        ----------

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        self.liouv.update(params_local)

        reference = self.reference
        carriers_i = self.carriers_i

        if offsets is not None:
            reference = np.zeros_like(offsets, dtype=np.bool)
            carriers_i = self.offsets_to_ppm(offsets)

        mag0 = self.liouv.compute_mag_eq(params_local, term="iz")

        profile = []

        for ref, carrier_i in zip(reference, carriers_i):
            self.liouv.carrier_i = carrier_i
            if not ref:
                cest = self.liouv.pulse_i(self.time_t1, 0.0, self.dephasing)
            else:
                cest = self.liouv.identity
            mag = self.liouv.collapse(self.detect @ cest @ mag0)
            profile.append(mag)

        return np.asarray(profile)
