'''
Created on Feb 22, 2012

@author: Mike Latham
@author: Guillaume Bouvignies
'''

from scipy import mat, zeros
from chemex.bases.two_states.mq import R_MQ, DR_MQ, DWI, DWS, KAB, KBA


def compute_liouvillian(pb=0.0, kex=0.0, dwc=0.0, dwh=0.0, r_mq=10.0, dr_mq=0.0):
    '''
    Compute the exchange matrix (Liouvillian).

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 8x8 basis, that is:
        {2HxCx, 2HxCy, 2HyCx, 2HyCy}{A,B}

    Parameters
    ----------
    pb : float
        Fractional population of state B
    kex : float
        Exchange rate between state A and B in /s
        kex = kab + kba
    dwc : float
        Carbon chemical shift difference between states A and B in rad/s
    dwh : float
        Proton chemical shift difference between states A and B in rad/s
    r_mq : float
        Multiple quantum relaxation rate in /s
    dr_mq : float
        Multiple quantum relaxation rate difference between A and B in /s

    Returns
    -------
    out : numpy.ndarray
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    '''

    kab = kex * pb
    kba = kex * (1.0 - pb)

    l_free = R_MQ * r_mq
    l_free += DR_MQ * dr_mq
    l_free += DWI * dwh
    l_free += DWS * dwc
    l_free += KAB * kab
    l_free += KBA * kba

    return l_free


def compute_2hxcy_eq(pb):
    mag_eq = mat(zeros((8, 1)))
    mag_eq[1, 0] += 1.0 - pb
    mag_eq[5, 0] += pb

    return mag_eq


def get_2hxcy(mag):
    mag_a = mag[1, 0]
    mag_b = mag[5, 0]

    return mag_a, mag_b
