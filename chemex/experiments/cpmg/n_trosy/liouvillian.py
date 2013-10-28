'''
Created on Aug 26, 2011

@author: guillaume
'''


# Imports
from scipy import (zeros,
                   asarray,
                   pi, cos, sqrt)
from scipy.constants import hbar, mu_0

from chemex.constants import gamma
from chemex.bases.two_states.iph_aph import (R_IXY, R_2SZIXY, DR_XY,
                                             R_IZ, R_2SZIZ, CS, DW,
                                             J, DJ, ETAXY, ETAZ,
                                             KAB, KBA, W1X, W1Y)

# Functions

def compute_liouvillians(pb=0.0, kex=0.0, dw=0.0, r_Nxy=5.0, dr_Nxy=0.0,
                         r_Nz=1.5, r_2HzNz=5.0, etaxy=0.0, etaz=0.0,
                         j_HN=-93.0, dj_HN=0.0, cs_offset=0.0, w1=0.0):

    '''
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
    r_Nz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_Nxy : float
        Transverse relaxation rate of state a in /s.
    dr_Nxy : float
        Transverse relaxation rate difference between states a and b in /s.
    r_2HzNz : float
        2-spin order longitudinal relaxation rate in /s.
    etaxy : float
        Transverse cross-correlated relaxation rate in /s.
    etaz : float
        Longitudinal cross-correlated relaxation rate in /s.
    j_HN : float
        Scalar coupling between N and HN in Hz.
    dj_HN : float
        Scalar coupling difference between states a and b in Hz.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    '''

    kab = kex * pb
    kba = kex - kab

    r_2HzNxy = r_Nxy + r_2HzNz - r_Nz

    l_free = R_IXY * r_Nxy
    l_free += R_2SZIXY * r_2HzNxy
    l_free += DR_XY * dr_Nxy
    l_free += R_IZ * r_Nz
    l_free += R_2SZIZ * r_2HzNz
    l_free += CS * cs_offset
    l_free += DW * dw
    l_free += J * pi * j_HN
    l_free += DJ * pi * dj_HN
    l_free += ETAXY * etaxy
    l_free += ETAZ * etaz
    l_free += KAB * kab
    l_free += KBA * kba

    l_w1x, l_w1y = w1 * asarray([W1X, W1Y])

    return l_free, l_w1x, l_w1y


def compute_2HzNz_eq(pb):

    Ieq = zeros((12, 1))
    Ieq[5, 0] += (1.0 - pb)
    Ieq[11, 0] += pb

    return Ieq


def get_Trz(I):

    Ia = I[5, 0] - I[2, 0]
    Ib = I[11, 0] - I[8, 0]

    return Ia, Ib


def compute_nh_etaz(r_Nz, ppm_to_rads):
    # TODO: replace with appropriate code for HN(dip)/H(csa) calculation
    '''
    THIS IS NOT PRESENTLY IMPLEMENTED.
    X-CORRELATION RATES ARE ZERO UNLESS EXPLICITLY SPECIFIED

    Approximates etaxy and etaz (NH-dipolar/N-csa cross-correlated relaxation
    rates) using inphase and longitudinal rates.

    Arguments:
    r_Hxy        -- transverse relaxationrate of N15 nucleus in /s,
                    r_Hxy = 5.0 (default)
    r_Hz         -- longitudinal relaxation rate of N15 nucleus in /s,
                    r_Hz = 1.5 (default)
    ppm_to_rads  -- unit conversion factor for desired nucleus & field


    Returns:
        etaxy, etaz
        float, float
    '''

    delta_csa_nh = -166.0  # ppm
    r_nh = 1.04e-10  # meters
    sqrt3 = sqrt(3.0)
    geo_factor = 0.5 * (3.0 * cos(19.6 / 180.0 * pi) ** 2 - 1.0)

    cc = delta_csa_nh * ppm_to_rads / sqrt3
    dd = -mu_0 * hbar * gamma['H'] * gamma['N'] / (4.0 * pi * r_nh ** 3)

    cc2 = cc ** 2
    dd2 = dd ** 2
    ccdd = cc * dd

    jwn = r_Nz / (cc2 + 0.75 * dd2)

    etaz = sqrt3 * geo_factor * ccdd * jwn

    return etaz

