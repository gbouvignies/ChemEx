'''
Created on Feb 22, 2012

@author: Mike Latham
@author: Guillaume Bouvignies
'''

# Imports
from scipy import mat, zeros

from chemex.caching import lru_cache

#Functions

@lru_cache()
def compute_liouvillian(pb=0.0, kex=0.0, dwc=0.0, dwh=0.0,
                        r_mq=10.0, csc_offset=0.0, csh_offset=0.0):
    '''Compute the exchange matrix (Liouvillian).

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 8x8 basis, that is:
        {2HxCx, 2HxCy, 2HyCx, 2HyCy}{A,B}

    Keyword arguments:
    pb                -- population of state B,
                         0.0 for 0% (default),
                         1.0 for 100%
    kex               -- exchange rate between state A and B in /s,
                         kex = kab + kba = 0.0 (default)
    dwc               -- carbon chemical shift difference between
                         states A and B in rad/s,
                         dw = 0.0 rad/s (default)
    dwh               -- proton chemical shift difference between
                         states A and B in rad/s,
                         dw = 0.0 rad/s (default)
    r_mq              -- multiple quantum relaxation rate in /s,
                         r_MQ = 10.0 (default)
    csc_offset        -- carbon chemical shift
                         csc_offset = 0.0 (default)
    csh_offset        -- proton chemical shift
                         csh_offset = 0.0 (default)

    Returns: numpy.matrix
    '''

    L = mat(zeros((8, 8)))

    k_ab = kex * pb
    k_ba = kex * (1.0 - pb)

    # Site A
    # Add auto-relaxation
    L[0, 0] -= r_mq
    L[1, 1] -= r_mq
    L[2, 2] -= r_mq
    L[3, 3] -= r_mq
    R_MQ = ones
    # Add chemical shift difference
    L[0, 1] -= csc_offset
    L[1, 0] += csc_offset
    L[2, 3] -= csc_offset
    L[3, 2] += csc_offset

    L[0, 2] -= csh_offset
    L[2, 0] += csh_offset
    L[1, 3] -= csh_offset
    L[3, 1] += csh_offset

    # Site B
    # Add auto-relaxation
    L[4, 4] -= r_mq
    L[5, 5] -= r_mq
    L[6, 6] -= r_mq
    L[7, 7] -= r_mq
    # Add chemical shift difference
    L[4, 5] -= csc_offset + dwc
    L[5, 4] += csc_offset + dwc
    L[6, 7] -= csc_offset + dwc
    L[7, 6] += csc_offset + dwc

    L[4, 6] -= csh_offset + dwh
    L[6, 4] += csh_offset + dwh
    L[5, 7] -= csh_offset + dwh
    L[7, 5] += csh_offset + dwh

    # Add chemical exchange
    # HxCx
    L[0, 0] -= k_ab
    L[4, 0] += k_ab
    L[4, 4] -= k_ba
    L[0, 4] += k_ba
    # HxCy
    L[1, 1] -= k_ab
    L[5, 1] += k_ab
    L[5, 5] -= k_ba
    L[1, 5] += k_ba
    # HyCx
    L[2, 2] -= k_ab
    L[6, 2] += k_ab
    L[6, 6] -= k_ba
    L[2, 6] += k_ba
    # HyCy
    L[3, 3] -= k_ab
    L[7, 3] += k_ab
    L[7, 7] -= k_ba
    L[3, 7] += k_ba

    return L


def compute_2HxCy_eq(pb):
    mag_eq = mat(zeros((8, 1)))
    mag_eq[1, 0] += (1.0 - pb)
    mag_eq[5, 0] += pb

    return mag_eq


def get_2HxCy(mag):
    mag_a = mag[1, 0]
    mag_b = mag[5, 0]

    return mag_a, mag_b
