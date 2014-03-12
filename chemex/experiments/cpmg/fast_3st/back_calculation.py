"""
Created on Aug 15, 2011

@author: guillaume
"""

from numpy.linalg import matrix_power
from scipy.linalg import expm

from chemex.bases.three_states.fast import P_180Y
from chemex.caching import lru_cache
from .liouvillian import compute_iy_eq, compute_liouvillians, get_iy



# Local Modules
@lru_cache()
def make_calc_observable(time_t2=0.0, ppm_to_rads=1.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CEST block.

    Parameters
    ----------
    time_T2 : float
        Time of the CPMG block.
    ncyc : integer
        Number of cycles, t-180-2t-180-t.
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Calculate intensity after the CEST block

    """

    @lru_cache(100)
    def _calc_observable(pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
                         dw_ab=0.0, dw_ac=0.0, r_ixy=5.0, dr_ixy_ab=0.0, dr_ixy_ac=0.0, ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.

        Parameters
        ----------
        i0 : float
            Initial intensity.
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        pc : float
            Fractional population of state C,
            0.0 for 0%, 1.0 for 100%
        kex_ab : float
            Exchange rate between state A and B in /s.
        kex_bc : float
            Exchange rate between state B and C in /s.
        kex_ac : float
            Exchange rate between state A and C in /s.
        dw_ab : float
            Chemical shift difference between states A and B in rad/s.
        dw_ac : float
            Chemical shift difference between states A and C in rad/s.
        r_ixy : float
            Transverse relaxation rate of state a in /s.
        dr_ixy_ab : float
            Transverse relaxation rate difference between states A and B in /s.
        dr_ixy_ac : float
            Transverse relaxation rate difference between states A and C in /s.
        cs_offset : float
            Offset from the carrier in rad/s.

        Returns
        -------
        out : float
            Intensity after the CPMG block

        """

        dw_ab *= ppm_to_rads
        dw_ac *= ppm_to_rads

        mag_eq = compute_iy_eq(pb, pc)

        if ncyc == 0:

            mag = mag_eq

        else:

            l_free = compute_liouvillians(pb=pb, pc=pc, kex_ab=kex_ab, kex_bc=kex_bc, kex_ac=kex_ac,
                                          dw_ab=dw_ab, dw_ac=dw_ac, r_ixy=r_ixy, dr_ixy_ab=dr_ixy_ab,
                                          dr_ixy_ac=dr_ixy_ac)

            t_cp = time_t2 / (4.0 * ncyc)
            p_free = expm(l_free * t_cp)

            mag = matrix_power(p_free.dot(P_180Y).dot(p_free), 2 * ncyc).dot(mag_eq)

        magy_a, _, _ = get_iy(mag)

        return magy_a

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
