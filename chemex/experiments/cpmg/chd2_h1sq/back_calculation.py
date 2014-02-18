"""
Created on Feb 12, 2014


@author: mike latham
"""

# Python Modules
from scipy import pi, dot
from scipy.linalg import expm
from numpy.linalg import matrix_power

# Local Modules
from chemex.caching import lru_cache
from .liouvillian import \
    compute_2hzcz_eq, \
    compute_liouvillians, \
    get_2hzcz


@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, ppm_to_rads=1.0, carrier=0.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after the CPMG-type pulse train.

    Parameters
    ----------
    pw : float
        Pulse width for a 90 degree pulse
    time_t2 : float
        Time of the CPMG block
    ppm_to_rads : float
        Conversion factor to go from ppm to rad/s
    carrier : float
        Position of the carrier in the indirect dimension in ppm
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Function for calculating intensity after the CPMG train pulse

    """

    @lru_cache(1)
    def make_propagators(pb=0.0, kex=0.0, dw=0.0, r_hxy=5.0, dr_hxy=0.0,
                         r_cz=1.5, r_2hzcz=0.0, etaxy=0.0, etaz=0.0,
                         j_hc=0.0, cs_offset=0.0):

        w1 = 2.0 * pi / (4.0 * pw)

        l_free, l_w1x, l_w1y = compute_liouvillians(
            pb=pb,
            kex=kex,
            dw=dw,
            r_hxy=r_hxy,
            dr_hxy=dr_hxy,
            r_cz=r_cz,
            r_2hzcz=r_2hzcz,
            etaxy=etaxy,
            etaz=etaz,
            j_hc=j_hc,
            cs_offset=cs_offset,
            w1=w1
        )

        p_neg = expm(l_free * -2.0 * pw / pi)
        p_90px = expm((l_free + l_w1x) * pw)
        p_90py = expm((l_free + l_w1y) * pw)
        p_90mx = expm((l_free - l_w1x) * pw)
        p_180py = matrix_power(p_90py, 2)
        p_180pmx = 0.5 * (matrix_power(p_90px, 2) +
                          matrix_power(p_90mx, 2))

        ps = p_neg, p_90px, p_90py, p_90mx, p_180py, p_180pmx

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_hxy=5.0, dr_hxy=0.0, r_cz=1.5,
                         r_2hzcz=0.0, etaxy=0.0, etaz=0.0, j_hc=0.0, cs=0.0, ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.

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
        r_cz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_hxy : float
            Transverse relaxation rate of state a in /s.
        dr_hxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs_offset : float
            Offset from the carrier in rad/s.
        ncyc : integer
            Number of cycles, t-180-2t-180-t

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
            r_cz=r_cz,
            r_2hzcz=r_2hzcz,
            etaxy=etaxy,
            etaz=etaz,
            j_hc=j_hc,
            cs_offset=cs_offset
        )

        p_neg, p_90px, p_90py, p_90mx, p_180py, p_180pmx = ps

        mag_eq = compute_2hzcz_eq(pb)

        if ncyc == 0:

            # The +/- phase cycling of the first 90 and the receiver is taken care
            # by setting the thermal equilibrium to 0
            I = -reduce(dot, [p_90px, p_180pmx, p_90px, mag_eq])

        else:

            t_cp = time_t2 / (4.0 * ncyc) - pw
            p_free = expm(l_free * t_cp)
            p_cpy = matrix_power(p_free.dot(p_180py).dot(p_free), ncyc)

            I = -reduce(dot, [p_90px, p_neg, p_cpy, p_180pmx, p_cpy, p_neg, p_90px, mag_eq])

        magz_a, _magz_b = get_2hzcz(I)

        return magz_a

    def calc_observable(i0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after the CPMG pulse train.

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

