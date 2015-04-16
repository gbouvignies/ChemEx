"""
Created on May 1, 2013

@author: guillaume
"""

# Imports
from itertools import product

from scipy import pi, zeros, linspace, asarray
from scipy.stats import norm

from chemex.bases.two_states.single_spin import R_IXY, DR_IXY, R_IZ, W, DW, KAB, \
    KBA, W1X


def compute_base_liouvillians(b1_offset=0.0, b1_frq=0.0, b1_inh=0.0, b1_inh_res=5, multiplet=None):
    w1, w1_inh, w1_offset = 2.0 * pi * asarray([b1_frq, b1_inh, b1_offset])

    w1s = linspace(-2.0, 2.0, b1_inh_res) * w1_inh + w1
    weights1 = norm.pdf(w1s, w1, w1_inh)

    liouvillians = [
        (j - w1_offset) * W + w1 * W1X
        for w1, (j, _) in product(w1s, multiplet)
    ]

    weights = [
        weight1 * weight2
        for weight1, (_, weight2) in product(weights1, multiplet)
    ]

    return liouvillians, weights


def compute_free_liouvillian(pb=0.0, kex=0.0, dw=0.0, r_cz=1.5, r_cxy=5.0, dr_cxy=0.0, cs_offset=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 6x6 cartesian basis, that is {Cx, Cy, Cz}{a,b}.
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
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    l_free : numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    """

    kab = kex * pb
    kba = kex - kab

    l_free = R_IXY * r_cxy
    l_free += DR_IXY * dr_cxy
    l_free += R_IZ * r_cz
    l_free += W * cs_offset
    l_free += DW * dw
    l_free += KAB * kab
    l_free += KBA * kba

    return l_free


def compute_cz_eq(pb):
    """
    Returns the equilibrium magnetization vector.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.

    Returns
    -------
    mag_eq : numpy.matrix
        Magnetization vector at equilibrium.

    """

    mag_eq = zeros((6, 1))
    mag_eq[2, 0] += (1.0 - pb)
    mag_eq[5, 0] += pb

    return mag_eq


def get_cz(mag):
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

