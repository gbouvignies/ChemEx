'''
Created on Sep 1, 2011

@author: guillaume
'''

from scipy import zeros

from chemex.bases.two_states.fast import R_IXY, DR_IXY, DW, KAB, KBA


def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0,
                         r_Ixy=5.0, dr_Ixy=0.0):
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
    r_Nz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_Nxy : float
        Transverse relaxation rate of state a in /s.
    dr_Nxy : float
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

    l_free = R_IXY * r_Ixy
    l_free += DR_IXY * dr_Ixy
    l_free += DW * dw
    l_free += KAB * kab
    l_free += KBA * kba

    return l_free


def compute_Iy_eq(pb):
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

    Ieq = zeros((4, 1))
    Ieq[1, 0] += (1.0 - pb)
    Ieq[3, 0] += pb

    return Ieq


def get_Iy(I):
    """
    Returns the amount of magnetization along z.

    Parameters
    ----------
    I : ndarray
        Magnetization vector.

    Returns
    -------
    Ia, Ib : float
        Amount of magnetization in state a and b along z.

    """

    Ia = I[1, 0]
    Ib = I[3, 0]

    return Ia, Ib

