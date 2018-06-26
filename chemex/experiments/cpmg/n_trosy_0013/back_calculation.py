"""
Created on Jun 15, 2017

@author: tyuwen
"""

# Python Modules
from scipy import pi, dot, eye
from scipy.linalg import expm
from numpy.linalg import matrix_power

# Local Modules
from chemex.caching import lru_cache
from .liouvillian import (compute_2hznz_eq,
                          compute_alpha_eq,
                          compute_beta_eq,
                          compute_liouvillians,
                          get_alpha, get_beta)
from chemex.bases.two_states.iph_aph import P180_S


@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, time_equil=0.0, ppm_to_rads=1.0, carrier=0.0,
                         taub=2.68e-3, ncyc_max=0, s3e_flg='y', antitrosy_flg='n', _id=None):
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
                         j_hn=-93.0, dj_hn=0.0, cs_offset=0.0):

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
        p_180mx = matrix_power(p_90mx, 2)
        p_180my = matrix_power(p_90my, 2)

        if antitrosy_flg.find('y') != -1:
            p_element = reduce(dot, [p_taub, p_90py, p_90px, P180_S, p_90px, p_90py, p_taub, P180_S])
        else:
            p_element = reduce(dot, [p_taub, p_90px, p_90py, P180_S, p_90py, p_90px, p_taub, P180_S])

        ps = (p_equil, p_neg, p_90px, p_90py, p_90mx, p_90my,
              p_180px, p_180py, p_180mx, p_180my, p_element)

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_nxy=5.0, dr_nxy=0.0, r_nz=1.5,
                         r_2hznz=0.0, etaxy=0.0, etaz=0.0, j_hn=-93.0, dj_hn=0.0,
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

        l_free, ps = make_propagators(pb=pb, kex=kex, dw=dw, r_nxy=r_nxy, dr_nxy=dr_nxy,
                                      r_nz=r_nz, r_2hznz=r_2hznz, etaxy=etaxy, etaz=etaz,
                                      j_hn=j_hn, dj_hn=dj_hn, cs_offset=cs_offset)

        (p_equil, p_neg, p_90px, p_90py, p_90mx, p_90my,
         p_180px, p_180py, p_180mx, p_180my, p_element) = ps

        p_element_pc = 0.5 * (p_90px.dot(p_element).dot(p_90py) + \
                              p_90mx.dot(p_element).dot(p_90my))

        if s3e_flg.find('y') != -1:
            if antitrosy_flg.find('y') != -1:
                mag_eq = compute_alpha_eq(pb)
            else:
                mag_eq = compute_beta_eq(pb)
        else:
            mag_eq = compute_2hznz_eq(pb)

        phase_1 = [1,1,0,2,1,1,2,0,1,1,2,0,1,1,0,2]
        phase_2 = [0,2,3,3,2,0,3,3,2,0,3,3,0,2,3,3]
        phase_3 = [0,0,1,3,0,0,3,1,0,0,3,1,0,0,1,3]
        phase_4 = [1,3,2,2,3,1,2,2,3,1,2,2,1,3,2,2]
        p_180s = [p_180px, p_180py, p_180mx, p_180my]

        if ncyc == 0:
            p_free_pw = expm(l_free * pw * 0.5)

            DELTA1 = matrix_power(p_free_pw, ncyc_max).dot(expm(l_free*-0.5*pw))

            # The +/- phase cycling of the first 90 and the receiver is taken care
            # by setting the thermal equilibrium to 0
            mag1 = reduce(dot, [p_equil, DELTA1, p_90my, p_180py, p_element_pc, p_180px, p_90px, DELTA1, mag_eq])
            mag2 = reduce(dot, [p_equil, DELTA1, p_90my, p_180px, p_element_pc, p_180py, p_90px, DELTA1, mag_eq])

        else:
            t_cp = time_t2 / (4.0 * ncyc) - pw*0.75
            p_free = expm(l_free * t_cp)
            p_free_pw = expm(l_free * pw * 0.5)
            DELTA = matrix_power(p_free_pw, ncyc_max-ncyc)
            p_cp1,p_cp2,p_cp3,p_cp4 = [eye(12)]*4

            for m in range(ncyc):
                p_cp1 = p_free.dot(p_180s[phase_1[m%len(phase_1)]]).dot(p_free).dot(p_cp1)
                p_cp2 = p_free.dot(p_180s[phase_2[m%len(phase_2)]]).dot(p_free).dot(p_cp2)
                p_cp3 = p_free.dot(p_180s[phase_3[(ncyc-1-m)%len(phase_3)]]).dot(p_free).dot(p_cp3)
                p_cp4 = p_free.dot(p_180s[phase_4[(ncyc-1-m)%len(phase_4)]]).dot(p_free).dot(p_cp4)
                
            mag1 = reduce(dot, [p_equil, DELTA, p_90my, p_neg, p_cp3, p_neg, \
                                p_element_pc, p_neg, p_cp1, p_neg, p_90px, DELTA, mag_eq])
            mag2 = reduce(dot, [p_equil, DELTA, p_90my, p_neg, p_cp4, p_neg, \
                                p_element_pc, p_neg, p_cp2, p_neg, p_90px, DELTA, mag_eq])

        mag = (mag1+mag2)/2.0

        if antitrosy_flg.find('y') != -1:
            magz_a, _magz_b = get_alpha(mag)
        else:
            magz_a, _magz_b = get_beta(mag)

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

