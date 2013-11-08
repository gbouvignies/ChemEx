'''
Created on Sep 1, 2011

@author: guillaume
'''

from scipy import zeros

from chemex.bases.three_states.fast import R_IXY, DR_IXY_AB, DR_IXY_AC, DW_AB, DW_AC, KAB, KBA, KBC, KCB, KAC, KCA


def compute_liouvillians(pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0, dw_ab=0.0, dw_ac=0.0,
                         r_ixy=5.0, dr_ixy_ab=0.0, dr_ixy_ac=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 6x6 cartesian basis, that is {Ix, Iy, Iz}{a,b}.
    Here the thermal equilibrium is assumed to be 0. This is justified because of
    the +/- phase cycling of the first 90 degree pulse at the beginning of the
    cpmg block.

    Parameters
    ----------
    pb : float
        Fractional population of state B,
        0.0 for 0%, 1.0 for 100%
    pc : float
        Fractional population of state C,
        0.0 for 0%, 1.0 for 100%
    kex_ab : float
        Exchange rate between state A and B in /s.
    kex_bc : float
        Exchange rate between state B and C in /s.
    kex_ac : float
        Exchange rate between state A and C in /s.
    dw_ab : float
        Chemical shift difference between states A and B in rad/s.
    dw_ac : float
        Chemical shift difference between states A and C in rad/s.
    r_ixy : float
        Transverse relaxation rate of state a in /s.
    dr_ixy_ab : float
        Transverse relaxation rate difference between states A and B in /s.
    dr_ixy_ac : float
        Transverse relaxation rate difference between states A and C in /s.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.array
        Liouvillian describing free precession of one
        isolated spin in presence of three-site exchange.

    """

    pa = 1.0 - pb - pc

    kab = kex_ab * pb / (pa + pb)
    kba = kex_ab * pa / (pa + pb)

    kbc = kex_bc * pc / (pb + pc)
    kcb = kex_bc * pb / (pb + pc)

    kac = kex_ac * pc / (pa + pc)
    kca = kex_ac * pa / (pa + pc)

    l_free = R_IXY * r_ixy
    l_free += DR_IXY_AB * dr_ixy_ab
    l_free += DR_IXY_AC * dr_ixy_ac
    l_free += DW_AB * dw_ab
    l_free += DW_AC * dw_ac
    l_free += KAB * kab
    l_free += KBA * kba
    l_free += KBC * kbc
    l_free += KCB * kcb
    l_free += KAC * kac
    l_free += KCA * kca

    return l_free


def compute_Iy_eq(pb, pc):
    """
    Returns the equilibrium magnetization vector.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.
    pc : float
        Fractional population of state C.
        0.0 for 0%, 1.0 for 100%.

    Returns
    -------
    out: numpy.matrix
        Magnetization vector at equilibrium.

    """

    mag_eq = zeros((6, 1))
    mag_eq[1, 0] += (1.0 - pb - pc)
    mag_eq[3, 0] += pb
    mag_eq[5, 0] += pc

    return mag_eq


def get_Iy(I):
    """
    Returns the amount of magnetization along y.

    Parameters
    ----------
    I : ndarray
        Magnetization vector.

    Returns
    -------
    magy_a, magy_b : float
        Amount of magnetization in state a and b along z.

    """

    magy_a = I[1, 0]
    magy_b = I[3, 0]
    magy_c = I[5, 0]

    return magy_a, magy_b, magy_c

