"""
Created on Aug 26, 2011

@author: guillaume
"""

from scipy import zeros, asarray, pi

from chemex.bases.two_states.iph_aph import \
    R_IXY, R_2SZIXY, DR_XY, R_IZ, R_2SZIZ, CS, DW, \
    J, ETAXY, ETAZ, KAB, KBA, W1X, W1Y


def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0, r_cxy=5.0, dr_cxy=0.0,
                         r_cz=1.5, r_2hzcz=5.0, etaxy=0.0, etaz=0.0,
                         j_hc=-93.0, cs_offset=0.0, w1=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 12x12 cartesian basis, that is:
        {Nx, Ny, Nz, 2HzNx, 2HzNy, 2HzNz}{a,b}.
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
    r_cz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_cxy : float
        Transverse relaxation rate of state a in /s.
    dr_cxy : float
        Transverse relaxation rate difference between states a and b in /s.
    r_2hzcz : float
        2-spin order longitudinal relaxation rate in /s.
    etaxy : float
        Transverse cross-correlated relaxation rate in /s.
    etaz : float
        Longitudinal cross-correlated relaxation rate in /s.
    j_hc : float
        Scalar coupling between N and HN in Hz.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.ndarray
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    """

    kab = kex * pb
    kba = kex - kab

    r_2hzcxy = r_cxy + r_2hzcz - r_cz

    l_free = R_IXY * r_cxy
    l_free += R_2SZIXY * r_2hzcxy
    l_free += DR_XY * dr_cxy
    l_free += R_IZ * r_cz
    l_free += R_2SZIZ * r_2hzcz
    l_free += CS * cs_offset
    l_free += DW * dw
    l_free += J * pi * j_hc
    l_free += ETAXY * etaxy
    l_free += ETAZ * etaz
    l_free += KAB * kab
    l_free += KBA * kba

    l_w1x, l_w1y = w1 * asarray([W1X, W1Y])

    return l_free, l_w1x, l_w1y


def compute_2hzcz_eq(pb):
    mag_eq = zeros((12, 1))
    mag_eq[5, 0] += (1.0 - pb)
    mag_eq[11, 0] += pb

    return mag_eq


def get_cz(mag):
    magz_a = mag[2, 0]
    magz_b = mag[8, 0]

    return magz_a, magz_b
