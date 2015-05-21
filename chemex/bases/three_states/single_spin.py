import numpy as np


# Define the basis for the liouvillian
# Axes: _XY, __Z
_XY, __Z = (
    np.diag([1.0, 1.0, 0.0]),
    np.diag([0.0, 0.0, 1.0]),
)

# States: B, C or A & B & C
___B, ___C, _ABC = (
    np.diag([0.0, 1.0, 0.0]),
    np.diag([0.0, 0.0, 1.0]),
    np.diag([1.0, 1.0, 1.0]),
)

# Auto-relaxation rates
R_IXY, R_IZ, DR_IXY_AB, DR_IXY_AC = (
    np.kron(_ABC, -_XY),
    np.kron(_ABC, -__Z),
    np.kron(___B, -_XY),
    np.kron(___C, -_XY),
)

# Chemical shifts
_W = [[+0.0, -1.0, +0.0],
      [+1.0, +0.0, +0.0],
      [+0.0, +0.0, +0.0]]

W, DW_AB, DW_AC = (
    np.kron(_ABC, _W),
    np.kron(___B, _W),
    np.kron(___C, _W),
)

# Exchange rates
KAB = np.kron([[-1.0, +0.0, +0.0],
               [+1.0, +0.0, +0.0],
               [+0.0, +0.0, +0.0]], np.eye(3))

KBA = np.kron([[+0.0, +1.0, +0.0],
               [+0.0, -1.0, +0.0],
               [+0.0, +0.0, +0.0]], np.eye(3))

KBC = np.kron([[+0.0, +0.0, +0.0],
               [+0.0, -1.0, +0.0],
               [+0.0, +1.0, +0.0]], np.eye(3))

KCB = np.kron([[+0.0, +0.0, +0.0],
               [+0.0, +0.0, +1.0],
               [+0.0, +0.0, -1.0]], np.eye(3))

KAC = np.kron([[-1.0, +0.0, +0.0],
               [+0.0, +0.0, +0.0],
               [+1.0, +0.0, +0.0]], np.eye(3))

KCA = np.kron([[+0.0, +0.0, +1.0],
               [+0.0, +0.0, +0.0],
               [+0.0, +0.0, -1.0]], np.eye(3))

# B1 field along x (in rad/s)
W1X = np.kron(_ABC, [[+0.0, +0.0, +0.0],
                     [+0.0, +0.0, -1.0],
                     [+0.0, +1.0, +0.0]])

# B1 field along y (in rad/s)
W1Y = np.kron(_ABC, [[+0.0, +0.0, -1.0],
                     [+0.0, +0.0, +0.0],
                     [+1.0, +0.0, +0.0]])


def compute_liouvillian(pb, pc, kex_ab, kex_bc, kex_ac, dw_ab, dw_ac, r_iz,
                        r_ixy, dr_ixy_ab, dr_ixy_ac, w, w1x=0.0, w1y=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 3-site (A <-> B) exchanging system.
    The matrix is written in 9x9 cartesian basis, that is {Nx, Ny, Nz}{a, b, c}.
    Here the thermal equilibrium is assumed to be 0. This is justified because of
    the +/- phase cycling of the first 90 degree pulse at the beginning of the
    cpmg block.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
    pc : float
        Fractional population of state C.
    kex_ab : float
        Exchange rate between state A and B in /s.
    kex_bc : float
        Exchange rate between state B and C in /s.
    kex_ac : float
        Exchange rate between state A and C in /s.
    w : float
        Resonance offset frequency of state A in rad/s.
    dw_ab : float
        Chemical shift difference between states A and B in rad/s.
    dw_ac : float
        Chemical shift difference between states A and C in rad/s.
    r_iz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_ixy : float
        Transverse relaxation rate of state a in /s.
    dr_ixy_ab : float
        Transverse relaxation rate difference between states A and B in /s.
    dr_ixy_ac : float
        Transverse relaxation rate difference between states A and C in /s.
    w1x : float, optional
        Angular frequency of the applied radio-frequency along x in rad/s.
    w1y: float, optional
        Angular frequency of the applied radio-frequency along y in rad/s.

    Returns
    -------
    out: numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of 3-site exchange.

    """

    pa = 1.0 - pb - pc

    kab = kex_ab * pb / (pa + pb)
    kba = kex_ab * pa / (pa + pb)

    kbc = kex_bc * pc / (pb + pc)
    kcb = kex_bc * pb / (pb + pc)

    kac = kex_ac * pc / (pa + pc)
    kca = kex_ac * pa / (pa + pc)

    liouvillian = (
        R_IXY * r_ixy +
        DR_IXY_AB * dr_ixy_ab +
        DR_IXY_AC * dr_ixy_ac +
        R_IZ * r_iz +
        W * w +
        DW_AB * dw_ab +
        DW_AC * dw_ac +
        w1x * W1X +
        w1y * W1Y +
        KAB * kab +
        KBA * kba +
        KBC * kbc +
        KCB * kcb +
        KAC * kac +
        KCA * kca
    )

    return liouvillian
