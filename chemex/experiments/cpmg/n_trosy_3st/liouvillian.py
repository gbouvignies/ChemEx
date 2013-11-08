'''
Created on Aug 26, 2011

@author: guillaume
'''


# Imports
from scipy import zeros, asarray, pi

from chemex.bases.three_states.iph_aph import \
    R_IXY, R_2SZIXY, DR_IXY_AB, DR_IXY_AC, R_IZ, R_2SZIZ, CS, DW_AB, DW_AC, \
    J, DJ_AB, DJ_AC, ETAXY, ETAZ, KAB, KBA, KBC, KCB, KAC, KCA, W1X, W1Y


# Functions
def compute_liouvillians(
        pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
        dw_ab=0.0, dw_ac=0.0, r_nxy=5.0, dr_nxy_ab=0.0, dr_nxy_ac=0.0,
        r_nz=1.5, r_2hznz=5.0, etaxy=0.0, etaz=0.0, j_hn=-93.0,
        dj_hn_ab=0.0, dj_hn_ac=0.0, cs_offset=0.0, w1=0.0):
    '''
    Compute the exchange matrix (Liouvillian)

    The function assumes a 3-site (A <-> B) exchanging system.
    The matrix is written in 12x12 cartesian basis, that is:
        {Nx, Ny, Nz, 2HzNx, 2HzNy, 2HzNz}{a,b}.
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
    r_nxy : float
        Transverse relaxation rate of state a in /s.
    dr_nxy_ab : float
        Transverse relaxation rate difference between states A and B in /s.
    dr_nxy_ac : float
        Transverse relaxation rate difference between states A and C in /s.
    r_nz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_2hznz : float
        2-spin order longitudinal relaxation rate in /s.
    etaxy : float
        Transverse cross-correlated relaxation rate in /s.
    etaz : float
        Longitudinal cross-correlated relaxation rate in /s.
    j_hn : float
        Scalar coupling between N and HN in Hz.
    dj_hn_ab : float
        Scalar coupling difference between states A and B in Hz.
    dj_hn_ac : float
        Scalar coupling difference between states A and C in Hz.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    '''

    pa = 1.0 - pb - pc

    kab = kex_ab * pb / (pa + pb)
    kba = kex_ab * pa / (pa + pb)

    kbc = kex_bc * pc / (pb + pc)
    kcb = kex_bc * pb / (pb + pc)

    kac = kex_ac * pc / (pa + pc)
    kca = kex_ac * pa / (pa + pc)

    r_2HzNxy = r_nxy + r_2hznz - r_nz

    l_free = R_IXY * r_nxy
    l_free += DR_IXY_AB * dr_nxy_ab
    l_free += DR_IXY_AC * dr_nxy_ac
    l_free += R_2SZIZ * r_2hznz
    l_free += R_IZ * r_nz
    l_free += R_2SZIXY * r_2HzNxy
    l_free += CS * cs_offset
    l_free += DW_AB * dw_ab
    l_free += DW_AC * dw_ac
    l_free += J * pi * j_hn
    l_free += DJ_AB * pi * dj_hn_ab
    l_free += DJ_AC * pi * dj_hn_ac
    l_free += ETAXY * etaxy
    l_free += ETAZ * etaz
    l_free += KAB * kab
    l_free += KBA * kba
    l_free += KBC * kbc
    l_free += KCB * kcb
    l_free += KAC * kac
    l_free += KCA * kca

    l_w1x, l_w1y = w1 * asarray([W1X, W1Y])

    return l_free, l_w1x, l_w1y


def compute_2HzNz_eq(pb, pc):
    mag_eq = zeros((18, 1))
    mag_eq[5, 0] += (1.0 - pb - pc)
    mag_eq[11, 0] += pb
    mag_eq[17, 0] += pc

    return mag_eq


def get_Trz(I):
    magz_a = I[5, 0] - I[2, 0]
    magz_b = I[11, 0] - I[8, 0]
    magz_c = I[17, 0] - I[14, 0]

    return magz_a, magz_b, magz_c

