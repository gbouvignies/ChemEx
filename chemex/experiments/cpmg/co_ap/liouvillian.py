"""
Created on Nov 12, 2014

@author: Mike Latham
"""

from scipy import zeros, asarray, pi

from chemex.bases.two_states.iph_aph import (R_IXY, R_2SZIXY, DR_XY, R_IZ, R_2SZIZ, CS, DW, J, DJ, ETAXY, ETAZ, KAB,
                                             KBA, W1X, W1Y)


def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0, r_coxy=5.0, dr_coxy=0.0, r_nz=1.5, r_2coznz=5.0, etaxy=0.0, etaz=0.0,
                         j_nco=-93.0, dj_nco=0.0, cs_offset=0.0, w1=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 12x12 cartesian basis, that is:
        {COx, COy, COz, 2COxNz, 2COyNz, 2COzNz}{a,b}.
    Here the thermal equilibrium is assumed to be 0. This is justified because of
    the +/- phase cycling of the first 90 degree pulse at the beginning of the
    cpmg block.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.
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
    r_2coznz : float
        2-spin order longitudinal relaxation rate in /s.
    etaxy : float
        Transverse cross-correlated relaxation rate in /s.
    etaz : float
        Longitudinal cross-correlated relaxation rate in /s.
    j_nco : float
        Scalar coupling between N and HN in Hz.
    dj_nco : float
        Scalar coupling difference between states a and b in Hz.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    """

    kab = kex * pb
    kba = kex - kab

    r_2coxynz = r_coxy - r_nz
    r_coz = r_2coznz - r_nz

    l_free = (
        R_IXY * r_coxy +
        R_2SZIXY * r_2coxynz +
        DR_XY * dr_coxy +
        R_IZ * r_coz +
        R_2SZIZ * r_2coznz +
        CS * cs_offset +
        DW * dw +
        J * pi * j_nco +
        DJ * pi * dj_nco +
        ETAXY * etaxy +
        ETAZ * etaz +
        KAB * kab +
        KBA * kba
    )

    l_w1x, l_w1y = w1 * asarray([W1X, W1Y])

    return l_free, l_w1x, l_w1y


def compute_2COzNz_eq(pb):
    mag_eq = zeros((12, 1))
    mag_eq[5, 0] += (1.0 - pb)
    mag_eq[11, 0] += pb

    return mag_eq


def get_2COzNz(I):
    magz_a = I[5, 0]
    magz_b = I[11, 0]

    return magz_a, magz_b


