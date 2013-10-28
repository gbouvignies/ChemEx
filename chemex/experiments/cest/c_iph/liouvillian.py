"""
Created on May 1, 2013

@author: guillaume
"""

# Imports
from scipy import (pi,
                   zeros,
                   linspace,
                   asarray)
from scipy.stats import norm

from chemex.bases.two_states.iph import (R_IXY, DR_IXY, R_IZ,
                                         CS, DW, KAB, KBA, W1X)


def compute_base_liouvillians(B1_offset=0.0, B1_frq=0.0, B1_inh=0.0, B1_inh_res=5):

    w1, w1_inh, w1_offset = 2.0 * pi * asarray([B1_frq, B1_inh, B1_offset])

    w1s = linspace(-2.0, 2.0, B1_inh_res) * w1_inh + w1
    weights = norm.pdf(w1s, w1, w1_inh)

    liouvillians = [-w1_offset * CS + w1 * W1X
                    for w1 in w1s]

    return liouvillians, weights


def compute_liouvillian_free_precession(pb=0.0, kex=0.0, dw=0.0,
                                        r_Cz=1.5, r_Cxy=5.0, dr_Cxy=0.0,
                                        cs_offset=0.0):
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
    r_Cz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_Cxy : float
        Transverse relaxation rate of state a in /s.
    dr_Cxy : float
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

    L = R_IXY * r_Cxy
    L += DR_IXY * dr_Cxy
    L += R_IZ * r_Cz
    L += CS * cs_offset
    L += DW * dw
    L += KAB * kab
    L += KBA * kba

    return L


def compute_Cz_eq(pb):
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

    Ieq = zeros((6, 1))
    Ieq[2, 0] += (1.0 - pb)
    Ieq[5, 0] += pb

    return Ieq


def get_Cz(I):
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

    Ia = I[2, 0]
    Ib = I[5, 0]

    return Ia, Ib
