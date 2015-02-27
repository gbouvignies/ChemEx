from scipy import pi, zeros, linspace, asarray
from scipy.stats import norm

from chemex.bases.three_states.iph import R_IXY, DR_IXY_AB, DR_IXY_AC, R_IZ, \
    CS, DW_AB, DW_AC, KAB, KBA, KAC, KCA, KBC, KCB, W1X


def compute_base_liouvillians(b1_offset=0.0, b1_frq=0.0, b1_inh=0.0,
                              b1_inh_res=5):

    w1, w1_inh, w1_offset = 2.0 * pi * asarray([b1_frq, b1_inh, b1_offset])

    w1s = linspace(-2.0, 2.0, b1_inh_res) * w1_inh + w1
    weights = norm.pdf(w1s, w1, w1_inh)

    liouvillians = [-w1_offset * CS + w1 * W1X for w1 in w1s]

    return liouvillians, weights


def compute_free_liouvillian(pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
                             dw_ab=0.0, dw_ac=0.0, r_nz=1.5, r_nxy=5.0,
                             dr_nxy_ab=0.0, dr_nxy_ac=0.0, cs_offset=0.0):
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

    pa = 1.0 - pb - pc

    kab = kex_ab * pb / (pa + pb)
    kba = kex_ab * pa / (pa + pb)

    kbc = kex_bc * pc / (pb + pc)
    kcb = kex_bc * pb / (pb + pc)

    kac = kex_ac * pc / (pa + pc)
    kca = kex_ac * pa / (pa + pc)

    l_free = (
        R_IXY * r_nxy +
        DR_IXY_AB * dr_nxy_ab +
        DR_IXY_AC * dr_nxy_ac +
        R_IZ * r_nz +
        CS * cs_offset +
        DW_AB * dw_ab +
        DW_AC * dw_ac +
        KAB * kab +
        KBA * kba +
        KBC * kbc +
        KCB * kcb +
        KAC * kac +
        KCA * kca
    )

    return l_free


def compute_nz_eq(pb, pc):
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

    mag_eq = zeros((9, 1))
    mag_eq[2, 0] += (1.0 - pb - pc)
    mag_eq[5, 0] += pb
    mag_eq[8, 0] += pc

    return mag_eq


def get_nz(mag):
    """
    Returns the amount of magnetization along z.

    Parameters
    ----------
    mag : ndarray
        Magnetization vector.

    Returns
    -------
    magz_a, magz_b, magz_c : float
        Amount of magnetization in state A, B and C along z.

    """

    magz_a = mag[2, 0]
    magz_b = mag[5, 0]
    magz_c = mag[8, 0]

    return magz_a, magz_b, magz_c
