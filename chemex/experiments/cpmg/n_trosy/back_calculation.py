"""
Created on Aug 15, 2011

@author: guillaume
"""

# Python Modules
from scipy import pi, dot
from scipy.linalg import expm
from numpy.linalg import matrix_power

# Local Modules
from chemex.caching import lru_cache
from .liouvillian import (compute_2HzNz_eq,
                          compute_liouvillians,
                          get_Trz, compute_nh_etaz)
from chemex.bases.two_states.iph_aph import P180_S


@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, time_equil=0.0, ppm_to_rads=1.0, carrier=0.0,
                         taub=2.68e-3, _id=None):
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

    @lru_cache(1)
    def make_propagators(pb=0.0, kex=0.0, dw=0.0, r_nxy=5.0, dr_nxy=0.0,
                         r_nz=1.5, r_2hznz=0.0, etaxy=0.0, etaz=0.0,
                         j_hn=0.0, dj_hn=0.0, cs_offset=0.0):

        w1 = 2.0 * pi / (4.0 * pw)
        l_free, l_w1x, l_w1y = compute_liouvillians(pb=pb, kex=kex, dw=dw,
                                                    r_nxy=r_nxy, dr_nxy=dr_nxy,
                                                    r_nz=r_nz, r_2hznz=r_2hznz,
                                                    etaxy=etaxy, etaz=etaz,
                                                    j_hn=j_hn, dj_hn=dj_hn,
                                                    cs_offset=cs_offset, w1=w1)

        p_equil = expm(l_free * time_equil)
        p_neg = expm(l_free * -2.0 * pw / pi)
        p_taub = expm(l_free * (taub - 2.0 * pw - 2.0 * pw / pi))
        p_90px = expm((l_free + l_w1x) * pw)
        p_90py = expm((l_free + l_w1y) * pw)
        p_90mx = expm((l_free - l_w1x) * pw)
        p_90my = expm((l_free - l_w1y) * pw)
        p_180px = matrix_power(p_90px, 2)
        p_180py = matrix_power(p_90py, 2)

        p_element = reduce(dot, [P180_S, p_taub, p_90py, p_90px, P180_S, p_90px, p_90py, p_taub])

        ps = (p_equil, p_neg, p_90px, p_90py, p_90mx, p_90my,
              p_180px, p_180py, p_element)

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_nxy=5.0, dr_nxy=0.0, r_nz=1.5,
                         r_2hznz=0.0, etaxy=0.0, etaz=0.0, j_hn=0.0, dj_hn=0.0,
                         cs=0.0, ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.
                _______________________________________________________________________
        1H :   |  /   /   /   /   /   /   /   /   CW   /   /   /   /   /   /   /   /   |
        15N:    Nx { tauc  2Ny  tauc }*ncyc 2Nx { tauc  2Ny  tauc }*ncyc -Nx time_equil

        Parameters
        ----------
        i0 : float
            Initial intensity.
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_nz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_nxy : float
            Transverse relaxation rate of state a in /s.
        dr_nxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs_offset : float
            Offset from the carrier in rad/s.
        Returns
        -------
        out : float
            Intensity after the CPMG block

        """

        dw *= ppm_to_rads
        cs_offset = (cs - carrier) * ppm_to_rads + pi * j_hn

        etaz_calc = compute_nh_etaz(r_nz, ppm_to_rads) * 0.0

        l_free, ps = make_propagators(pb=pb, kex=kex, dw=dw, r_nxy=r_nxy, dr_nxy=dr_nxy,
                                      r_nz=r_nz, r_2hznz=r_2hznz, etaxy=etaxy, etaz=etaz_calc,
                                      j_hn=j_hn, dj_hn=dj_hn, cs_offset=cs_offset)

        (p_equil, p_neg, p_90px, p_90py, p_90mx,
         p_90my, p_180px, p_180py, p_element) = ps

        mag_eq = compute_2HzNz_eq(pb)

        if ncyc == 0:

            # The +/- phase cycling of the first 90 and the receiver is taken care
            # by setting the thermal equilibrium to 0
            I = -reduce(dot, [p_equil, p_90py, p_element, p_90px, mag_eq])

        else:

            t_cp = time_t2 / (4.0 * ncyc) - pw
            p_free = expm(l_free * t_cp)
            p_cpx = matrix_power(p_free.dot(p_180px).dot(p_free), ncyc)
            p_cpy = matrix_power(p_free.dot(p_180py).dot(p_free), ncyc)
            p_element_pc = 0.5 * (p_90px.dot(p_element).dot(p_90py) +
                                  p_90mx.dot(p_element).dot(p_90my))

            I = -reduce(dot, [p_equil, p_90py, p_neg, p_cpx, p_neg, p_element_pc,
                              p_neg, p_cpy, p_neg, p_90px, mag_eq])

        magz_a, _magz_b = get_Trz(I)

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

