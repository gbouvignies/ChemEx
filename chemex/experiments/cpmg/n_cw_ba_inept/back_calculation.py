import scipy as sp
from scipy import pi, dot
from scipy.linalg import expm
from numpy.linalg import matrix_power

from chemex.caching import lru_cache
from chemex.experiments import misc
from .liouvillian import (
    compute_liouvillians,
    get_nz
)


correct_intensities = misc.correct_intensities


@lru_cache()
def make_calc_observable(pw=0.0, time_t2=0.0, time_equil=0.0, ppm_to_rads=1.0,
                         ppm_to_rads_h=1.0, carrier=0.0, _id=None):
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
                         r_nz=1.5, cs_offset=0.0):

        w1 = 2.0 * pi / (4.0 * pw)
        l_free, l_w1x, l_w1y = compute_liouvillians(pb=pb, kex=kex, dw=dw,
                                                    r_nxy=r_nxy, dr_nxy=dr_nxy,
                                                    r_nz=r_nz,
                                                    cs_offset=cs_offset, w1=w1)

        p_equil = expm(l_free * time_equil)
        p_neg = expm(l_free * -2.0 * pw / pi)
        p_90px = expm((l_free + l_w1x) * pw)
        p_90py = expm((l_free + l_w1y) * pw)
        p_90mx = expm((l_free - l_w1x) * pw)
        p_180pmx = 0.5 * (matrix_power(p_90px, 2) +
                          matrix_power(p_90mx, 2))
        p_180py = matrix_power(p_90py, 2)

        ps = (p_equil, p_neg, p_90px, p_90mx, p_180pmx, p_180py)

        return l_free, ps

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, dw_h=0.0, r_nxy=5.0,
                         dr_nxy=0.0, r_nz=1.5, cs=0.0, ncyc=0):
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

        cs_offset = (cs - carrier) * ppm_to_rads
        dw_h *= ppm_to_rads_h
        dw *= ppm_to_rads

        kab = pb * kex
        kba = kex - kab

        liouv_h1 = sp.array([[-kab, +kba],
                             [+kab, -kba - dr_nxy + 1j * dw_h]])

        liouv_n1 = sp.array([[-kab, +kba],
                             [+kab, -kba - dr_nxy + 1j * dw]])

        magz = sp.array([[1.0 - pb], [pb]])

        t1 = 2.25e-3
        t2 = 2.75e-3

        magz = (
            expm(liouv_h1.conjugate() * t1)
            .dot(expm(liouv_h1 * t1))
            .dot(magz)
        ).real

        magz = (
            expm(liouv_n1.conjugate() * t2)
            .dot(expm(liouv_n1 * t2))
            .dot(magz)
        ).real

        l_free, ps = make_propagators(pb=pb, kex=kex, dw=dw, r_nxy=r_nxy,
                                      dr_nxy=dr_nxy, r_nz=r_nz,
                                      cs_offset=cs_offset)

        p_equil, p_neg, p_90px, p_90mx, p_180pmx, p_180py = ps

        mag_eq = sp.array([[0.0, 0.0, magz[0, 0], 0.0, 0.0, magz[1, 0]]]).T

        if ncyc == 0:
            # The +/- phase cycling of the first 90 and the receiver is taken care
            # by setting the thermal equilibrium to 0
            mag = reduce(dot, [p_equil, p_90px, p_180pmx, p_90px, mag_eq])

        else:

            t_cp = time_t2 / (4.0 * ncyc) - pw
            p_free = expm(l_free * t_cp)
            p_cp = matrix_power(p_free.dot(p_180py).dot(p_free), ncyc)

            mag = reduce(
                dot,
                [
                    p_equil,
                    p_90px,
                    p_neg,
                    p_cp,
                    p_180pmx,
                    p_cp,
                    p_neg,
                    p_90px,
                    mag_eq
                ]
            )

        magz_a, magz_b = get_nz(mag)

        magz_a, magz_b = correct_intensities(
            magz_a,
            magz_b,
            pb,
            kex,
            dw,
            r_nxy,
            dr_nxy
        )

        magz = sp.array([[magz_a], [magz_b]])

        magz = (
            expm(liouv_n1.conjugate() * t2)
            .dot(expm(liouv_n1 * t2))
            .dot(magz)
        ).real

        magz = (
            expm(liouv_h1.conjugate() * t1)
            .dot(expm(liouv_h1 * t1))
            .dot(magz)
        ).real

        magz_a, magz_b = magz[:, 0]

        magz_a, magz_b = correct_intensities(
            magz_a,
            magz_b,
            pb,
            kex,
            dw,
            r_nxy,
            dr_nxy
        )

        return magz_b

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
