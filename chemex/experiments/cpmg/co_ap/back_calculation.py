"""
Created on Nov 12, 2014

@author: Mike Latham
"""

from scipy import pi, dot, diag
from scipy.linalg import expm2 as expm
from numpy.linalg import matrix_power

from chemex.caching import lru_cache
from .liouvillian import compute_2COzNz_eq, compute_liouvillians, get_2COzNz


P180X = diag(4 * [+1.0, -1.0, +1.0])


@lru_cache()
def make_calc_observable(pwco90=0.0, time_t2=0.0, time_equil=0.0, taucc=0.0, ppm_to_rads=1.0, sidechain_flg='N', carrier=0.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CPMG block.

    Parameters
    ----------
    pwco90 : float
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
    def make_propagators(pb=0.0, kex=0.0, dw=0.0, r_coxy=5.0, dr_coxy=0.0,
                         r_nz=1.5, r_2coznz=0.0, etaxy=0.0, etaz=0.0,
                         j_nco=0.0, dj_nco=0.0, cs_offset=0.0):

        w1 = 2.0 * pi / (4.0 * pwco90)
        l_free, l_w1x, l_w1y = compute_liouvillians(pb=pb, kex=kex, dw=dw,
                                                    r_coxy=r_coxy, dr_coxy=dr_coxy,
                                                    r_nz=r_nz, r_2coznz=r_2coznz,
                                                    etaxy=etaxy, etaz=etaz,
                                                    j_nco=j_nco, dj_nco=dj_nco,
                                                    cs_offset=cs_offset, w1=w1)

        p_equil = expm(l_free * time_equil)
        p_taucc = expm(l_free * taucc)
        p_neg = expm(l_free * -2.0 * pwco90 / pi)
        p_90py = expm((l_free + l_w1y) * pwco90)
        p_90my = expm((l_free - l_w1y) * pwco90)
        p_180px = P180X  # Perfect 180 for CPMG blocks
        p_180py = matrix_power(p_90py, 2)
        p_180my = matrix_power(p_90py, 2)

        ps = (p_equil, p_taucc, p_neg, p_90py, p_90my,
              p_180px, p_180py, p_180my)

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_coxy=5.0, dr_coxy=0.0, r_nz=1.5,
                         r_2coznz=0.0, etaxy=0.0, etaz=0.0, j_nco=0.0, dj_nco=0.0,
                         cs=0.0, ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.
               
        13CO : (2COzNz)-time_eq-COy-(2COxNz)-{ tcp-(2COx)-tcp }*ncyc-2COy-{ tcp-(2COx)-tcp }*ncyc-(2COxNz)-COy-(
        2COzNz)-time_eq
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
        r_coxy : float
            Transverse relaxation rate of state a in /s.
        dr_coxy : float
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

        l_free, ps = make_propagators(pb=pb, kex=kex, dw=dw, r_coxy=r_coxy,
                                      dr_coxy=dr_coxy, r_nz=r_nz, r_2coznz=r_2coznz, etaxy=etaxy,
                                      etaz=etaz, j_nco=j_nco, dj_nco=dj_nco, cs_offset=cs_offset)

        p_equil, p_taucc, p_neg, p_90py, p_90my, p_180px, p_180py, p_180my = ps

        mag_eq = compute_2COzNz_eq(pb)

        if sidechain_flg == 'N':
            p_flip = 0.5 * (p_180py + p_180my)
        else:
            p_flip = reduce(dot, [p_90my, p_taucc, 0.5*(p_180py + p_180my),
                                  p_taucc, p_90py])

        if ncyc == 0:

            # The +/- phase cycling of the first 90 and the receiver is taken care
            # by setting the thermal equilibrium to 0
            #I = reduce(dot, [p_equil, p_90py, 0.5 * (p_180py + p_180my), p_90py, p_equil, mag_eq])
            I = reduce(dot, [p_equil, p_90py, p_flip, p_90py, p_equil, mag_eq])

        else:

            t_cp = time_t2 / (4.0 * ncyc)
            p_free = expm(l_free * t_cp)
            p_cpx = matrix_power(p_free.dot(p_180px).dot(p_free), ncyc)

            I = reduce(dot, [p_equil, p_90py, p_neg, p_cpx, p_neg,
                             p_flip, p_neg, p_cpx, p_neg, p_90py,
                             p_equil, mag_eq])

        magz_a, _magz_b = get_2COzNz(I)

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

