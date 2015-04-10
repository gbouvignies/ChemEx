"""
Created on June 26, 2014

@author: Mike Latham
"""

from scipy import pi, dot
from scipy.linalg import expm
from numpy.linalg import matrix_power

from ....caching import lru_cache
from .liouvillian import (compute_2hznz_eq,
                          compute_liouvillians,
                          get_2hznz, )


@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, time_equil=0.0, ppm_to_rads=1.0,
                         carrier=0.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in
    presence of exchange after a CPMG block with 'antiphase_flg' set to 'y'.

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
        Calculate intensity after the CPMG block

    """

    @lru_cache(1)
    def make_propagators(pb=0.0, kex=0.0, dw=0.0, r_hxy=5.0, dr_hxy=0.0,
                         r_nz=1.5, r_2hznz=0.0, etaxy=0.0, etaz=0.0,
                         j_hn=0.0, dj_hn=0.0, cs_offset=0.0):

        w1 = 2.0 * pi / (4.0 * pw)

        l_free, l_w1x, l_w1y = compute_liouvillians(
            pb=pb,
            kex=kex,
            dw=dw,
            r_hxy=r_hxy,
            r_nz=r_nz,
            r_2hznz=r_2hznz,
            dr_hxy=dr_hxy,
            etaxy=etaxy,
            etaz=etaz,
            j_hn=j_hn,
            dj_hn=dj_hn,
            cs_offset=cs_offset,
            w1=w1
        )

        p_equil = expm(l_free * time_equil)
        p_neg = expm(l_free * -2.0 * pw / pi)
        p_90px = expm((l_free + l_w1x) * pw)
        p_90py = expm((l_free + l_w1y) * pw)
        p_90mx = expm((l_free - l_w1x) * pw)
        p_90my = expm((l_free - l_w1y) * pw)
        p_180px = matrix_power(p_90px, 2)
        p_180mx = matrix_power(p_90mx, 2)
        p_180py = matrix_power(p_90py, 2)

        ps = (p_equil, p_neg, p_90px, p_90py, p_90mx, p_90my, p_180px, p_180mx,
              p_180py)

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_hxy=5.0, dr_hxy=0.0,
                         r_nz=1.5, r_2hznz=0.0, etaxy=0.0, etaz=0.0, j_hn=0.0,
                         dj_hn=0.0, cs=0.0, ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.
               
        1H : (2HzNz)-time_eq-Hx-(2HyNz)-{ tcp-2Hy-tcp }*ncyc-2Hx-{ tcp-2Hy-tcp }*ncyc-(2HyNz)-Hx-(2HzNz) 
        15N:    

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
        r_hxy : float
            Transverse relaxation rate of state a in /s.
        dr_hxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs_offset : float
            Offset from the carrier in rad/s.
        Returns
        -------
        out : float
            Intensity after the CPMG block

        """

        dw *= ppm_to_rads
        cs_offset = (cs - carrier) * ppm_to_rads

        l_free, ps = make_propagators(
            pb=pb,
            kex=kex,
            dw=dw,
            r_hxy=r_hxy,
            dr_hxy=dr_hxy,
            r_nz=r_nz,
            r_2hznz=r_2hznz,
            etaxy=etaxy,
            etaz=etaz,
            j_hn=j_hn,
            dj_hn=dj_hn,
            cs_offset=cs_offset
        )

        (p_equil, p_neg, p_90px, p_90py, p_90mx, p_90my, p_180px, p_180mx,
         p_180py) = ps

        mag_eq = compute_2hznz_eq(pb)

        if ncyc == 0:

            # The +/- phase cycling of the first 90 and the receiver is taken
            # care by setting the thermal equilibrium to 0
            mag = reduce(
                dot,
                [
                    p_equil,
                    p_90px,
                    0.5 * (p_180px + p_180mx),
                    p_90px,
                    mag_eq
                ]
            )

        else:

            t_cp = time_t2 / (4.0 * ncyc) - pw
            p_free = expm(l_free * t_cp)
            p_cpy = matrix_power(
                p_free
                .dot(p_180py)
                .dot(p_free), ncyc)

            mag = reduce(
                dot,
                [
                    p_equil,
                    p_90px,
                    p_neg,
                    p_cpy,
                    p_neg,
                    0.5 * (p_180px + p_180mx),
                    p_neg,
                    p_cpy,
                    p_neg,
                    p_90px,
                    mag_eq
                ]
            )

        magz_a, _magz_b = get_2hznz(mag)

        return magz_a

    def calc_observable(i0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after a CPMG block.

        Parameters
        ----------
        i0 : float
            Initial intensity.

        Returns
        -------
        out : float
            Intensity after the CPMG block

        """

        return i0 * _calc_observable(**kwargs)

    return calc_observable

