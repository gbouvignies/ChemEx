"""
Created on Sep 1, 2011

@author: guillaume
"""

from scipy import zeros

from chemex.bases.two_states.fast import R_IXY, DR_IXY, DW, KAB, KBA


def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0,
                         r_ixy=5.0, dr_ixy=0.0):
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

    l_free = R_IXY * r_ixy
    l_free += DR_IXY * dr_ixy
    l_free += DW * dw
    l_free += KAB * kab
    l_free += KBA * kba

    return l_free


def compute_iy_eq(pb):
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

    mag_eq = zeros((4, 1))
    mag_eq[1, 0] += (1.0 - pb)
    mag_eq[3, 0] += pb

    return mag_eq


def get_iy(mag):
    """
    Returns the amount of magnetization along z.

    Parameters
    ----------
    mag : ndarray
        Magnetization vector.

    Returns
    -------
    magy_a, magy_b : float
        Amount of magnetization in state a and b along z.

    """

    magy_a = mag[1, 0]
    magy_b = mag[3, 0]

    return magy_a, magy_b

