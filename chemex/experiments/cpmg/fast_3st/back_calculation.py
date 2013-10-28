'''
Created on Aug 15, 2011

@author: guillaume
'''

from numpy.linalg import matrix_power
from scipy.linalg import expm

from chemex.bases.two_states.fast import P_180Y
from chemex.caching import lru_cache

from .liouvillian import compute_Iy_eq, compute_liouvillians, get_Iy


# Local Modules
@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, ppm_to_rads=1.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CEST block.

    Parameters
    ----------
    pw : float
        Pulse width for a 90 degree pulse.
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
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_Ixy=5.0, dr_Ixy=0.0, ncyc=0):
        '''
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.
                _______________________________________________________________________
        1H :   |  /   /   /   /   /   /   /   /   CW   /   /   /   /   /   /   /   /   |
        15N:    Nx { tauc  2Ny  tauc }*ncyc 2Nx { tauc  2Ny  tauc }*ncyc -Nx time_equil

        Parameters
        ----------
        I0 : float
            Initial intensity.
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_Ixy : float
            Transverse relaxation rate of state a in /s.
        dr_Ixy : float
            Transverse relaxation rate difference between states a and b in /s.

        Returns
        -------
        out : float
            Intensity after the CPMG block

        '''

        dw *= ppm_to_rads

        Ieq = compute_Iy_eq(pb)

        if ncyc == 0:

            I = Ieq

        else:

            l_free = compute_liouvillians(pb=pb, kex=kex, dw=dw,
                                          r_Ixy=r_Ixy, dr_Ixy=dr_Ixy)

            t_cp = time_t2 / (4.0 * ncyc)
            p_free = expm(l_free * t_cp)

            I = matrix_power(p_free.dot(P_180Y).dot(p_free), 2 * ncyc).dot(Ieq)

        Ia, _Ib = get_Iy(I)

        return Ia

    def calc_observable(I0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        I0 : float
            Initial intensity.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        return I0 * _calc_observable(**kwargs)

    return calc_observable
