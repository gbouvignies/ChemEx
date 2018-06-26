"""
Created on Jun 15, 2017

@author: tyuwen
"""

# Imports
from scipy import zeros, asarray

from chemex.bases.two_states.iph import \
    R_IXY, DR_IXY, R_IZ, CS, DW, KAB, KBA, W1X, W1Y


def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0,
                         r_nz=1.5, r_nxy=5.0, dr_nxy=0.0,
                         cs_offset=0.0, w1=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 6x6 cartesian basis, that is {Nx, Ny, Nz}{a,b}.
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
    r_nxy : float
        Transverse relaxation rate of state a in /s.
    dr_nxy : float
        Transverse relaxation rate difference between states a and b in /s.
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

    l_free = R_IXY * r_nxy
    l_free += DR_IXY * dr_nxy
    l_free += R_IZ * r_nz
    l_free += CS * cs_offset
    l_free += DW * dw
    l_free += KAB * kab
    l_free += KBA * kba

    l_w1x, l_w1y = w1 * asarray([W1X, W1Y])

    return l_free, l_w1x, l_w1y


def compute_nz_eq(pb):
    """
    Returns the equilibrium magnetization vector.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.

    Returns
    -------
    out: numpy.matrix
        Magnetization vector at equilibrium.

    """

    mag_eq = zeros((6, 1))
    mag_eq[2, 0] += 1.0
    mag_eq[5, 0] += pb/(1.0 - pb)

    return mag_eq


def get_nz(mag):
    """
    Returns the amount of magnetization along z.

    Parameters
    ----------
    mag : ndarray
        Magnetization vector.

    Returns
    -------
    magz_a, magz_b : float
        Amount of magnetization in state a and b along z.

    """

    magz_a = mag[2, 0]
    magz_b = mag[5, 0]

    return magz_a, magz_b

