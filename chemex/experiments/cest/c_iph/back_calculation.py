"""
Created on Aug 15, 2011

@author: guillaume
"""

# Python Modules
import scipy as sc
from scipy.linalg import expm

# Local Modules
from chemex.experiments.util import calculate_shift_ex_2st
from chemex.caching import lru_cache
from .liouvillian import (compute_cz_eq,
                          compute_base_liouvillians,
                          compute_free_liouvillian,
                          get_cz)


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, b1_inh=0.0, b1_inh_res=5,
                         carrier=0.0, ppm_to_rads=0.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CEST block.

    Parameters
    ----------
    time_t1 : float
        Duration of the CW block.
    b1_offset : float
        Frequency offset of the applied B1 field in Hz.
    b1_frq : float
        Strength of the applied B1 field in Hz.
    b1_inh : float
        B1 field inhomogeneity in Hz.
    b1_inh_res : int
        Resolution to model B1 field inhomogeneity.
    carrier : float
        Carrier position in rad/s.
    ppm_to_rads : float
        Conversion factor from ppm to rad/s
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Calculate intensity after the CEST block

    """

    base_liouvillians, weights = compute_base_liouvillians(b1_offset, b1_frq, b1_inh, b1_inh_res)

    @lru_cache(5)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_cz=1.5, r_cxy=0.0, dr_cxy=0.0, cs=0.0):
        """
        Calculate the intensity in presence of exchange after a CEST block assuming
        initial intensity of 1.0.

        Parameters
        ----------
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_cz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_cxy : float
            Transverse relaxation rate of state a in /s.
        dr_cxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs : float
            Resonance position in rad/s.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        if abs(b1_offset) >= 10000.0:

            return 1.0 - pb

        else:

            dw *= ppm_to_rads
            mag_eq = compute_cz_eq(pb)
            exchange_induced_shift, _ = calculate_shift_ex_2st(pb=pb, kex=kex,
                                                               dw=dw,
                                                               r_ixy=r_cxy, dr_ixy=dr_cxy)
            wg = (cs - carrier) * ppm_to_rads - exchange_induced_shift

            liouvillians = base_liouvillians + compute_free_liouvillian(pb=pb, kex=kex, dw=dw,
                                                                        r_cxy=r_cxy, dr_cxy=dr_cxy,
                                                                        r_cz=r_cz, cs_offset=wg)

            propagator = sc.zeros_like(liouvillians[0])
            for liouvillian, weight in zip(liouvillians, weights):
                propagator += weight * expm(liouvillian * time_t1)
            propagator /= sum(weights)

            magz_a, _ = get_cz(sc.dot(propagator, mag_eq))

            return magz_a

    def calc_observable(i0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        i0 : float
            Initial intensity.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        return i0 * _calc_observable(**kwargs)

    return calc_observable
