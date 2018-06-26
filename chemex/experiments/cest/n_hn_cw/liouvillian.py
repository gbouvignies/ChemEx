from math import pi

from scipy import zeros, cos, sin
from scipy.linalg import expm

from chemex import caching
from chemex.constants import scalar_couplings
from chemex.bases.two_states.full import (
    R_HXY, R_HZ, R_NXY_A, R_NXY_B, R_NZ, R_2HXYNZ, R_2HZNXY_A, R_2HZNXY_B,
    R_2HXYNXY, R_2HZNZ, CS_H_A, CS_H_B, CS_N_A, CS_N_B, J_HN, ETAZ, ETAXY,
    W1X_H, W1Y_H, W1X_N, W1Y_N, KAB, KBA
)


JHN = scalar_couplings['amide_HN']


@caching.lru_cache()
def compute_liouvillian_free_precession(pb=0.0, kex=0.0, dw_h=0.0, dw_n=0.0,
                                        r_nxy=5.0, dr_nxy=0.0, r_nz=1.5,
                                        r_2hznz=0.0, r_2hxynxy=0.0, r_hxy=10.0,
                                        r_hz=1.0, etaxy=0.0, etaz=0.0,
                                        cs_offset_h=0.0, cs_offset_n=0.0,
                                        j_hn=JHN):
    """Compute the exchange matrix (Liouvillian).
    
    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 13x13 basis, that is E/2, {Ix, Iy, Nz, 2IxSz, 2IySz, 2NzSz}{A,B}
    
    Keyword arguments:
    pb -- population of state B,
          0.0 for 0% (default),
          1.0 for 100%
    kex -- exchange rate between state A and B in /s,
           kex = kab + kba = 0.0 (default)
    dw -- chemical shift difference between states A and B in rad/s,
          dw = 0.0 rad/s (default)
    r_Nxy -- transverse relaxation rate in /s,
             r_Nxy = 5.0 (default)
    r_Nz -- longitudinal relaxation rate in /s,
            r_Nz = 1.5 (default)
    r_2HzNz -- 2-spin order longitudinal relaxation rate in /s,
               r_2HzNz = 5.0 (default)
    etaxy -- transverse cross-correlated relaxation rate in /s,
               eta_xy = 0.0 (default)
    etaz -- longitudinal cross-correlated relaxation rate in /s,
               eta_z = 0.0 (default)
    j_HN -- scalar coupling between I and S in Hz,
            j_HN = -92 (default)
    cs_offset -- offset from the carrier in rad/s,
                 frqr_offset = 0.0 (default)
    atom_name -- nucleus that is considered
                 atom_name = 'N' (default)

    Returns: numpy.matrix
    """

    kab = kex * pb
    kba = kex - kab

    r_2sf = r_2hznz - r_nz
    r_2hznxy = r_nxy + r_2sf
    r_2hxynz = r_hxy - r_nz

    liouvillian = (
        R_HXY * r_hxy +
        R_HZ * r_hz +
        R_NXY_A * r_nxy +
        R_NXY_B * (r_nxy + dr_nxy) +
        R_NZ * r_nz +
        R_2HXYNZ * r_2hxynz +
        R_2HZNXY_A * r_2hznxy +
        R_2HZNXY_B * (r_2hznxy + dr_nxy) +
        R_2HXYNXY * r_2hxynxy +
        R_2HZNZ * r_2hznz +
        CS_H_A * cs_offset_h +
        CS_H_B * (cs_offset_h + dw_h) +
        CS_N_A * cs_offset_n +
        CS_N_B * (cs_offset_n + dw_n) +
        J_HN * pi * j_hn +
        ETAXY * etaxy +
        ETAZ * etaz +
        KAB * kab +
        KBA * kba
    )

    return liouvillian


def set_nz(pb):
    mag_eq = zeros((30, 1))

    mag_eq[5, 0] = (1.0 - pb)
    mag_eq[5 + 15, 0] = pb

    return mag_eq


def get_nz(mag):
    ia = mag[5, 0]
    ib = mag[5 + 15, 0]

    return ia, ib


def make_pulse_nh(liouvillian, w1_h=0.0, phase_h=0.0, w1_n=0.0, phase_n=0.0,
                  pw=0.0):
    phase_rad_h = 0.5 * phase_h * pi
    phase_rad_n = 0.5 * phase_n * pi

    w1x_h = w1_h * cos(phase_rad_h)
    w1y_h = w1_h * sin(phase_rad_h)

    w1x_n = w1_n * cos(phase_rad_n)
    w1y_n = w1_n * sin(phase_rad_n)

    new_liouvillian = (
        liouvillian +
        w1x_h * W1X_H +
        w1y_h * W1Y_H +
        w1x_n * W1X_N +
        w1y_n * W1Y_N
    )

    return expm(new_liouvillian * pw)
