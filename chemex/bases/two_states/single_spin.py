import numpy as np

from scipy import linalg

exp = np.exp
dot = np.dot
diag = np.diag


# Define the basis for the liouvillian
# Axes: _XY, __Z
_XY, __Z = (
    np.diag([1.0, 1.0, 0.0]),
    np.diag([0.0, 0.0, 1.0]),
)

# States: B or both A & B
__B, _AB = (
    np.diag([0.0, 1.0]),
    np.diag([1.0, 1.0]),
)

# Auto-relaxation rates
R_IXY, R_IZ, DR_IXY = (
    np.kron(_AB, -_XY),
    np.kron(_AB, -__Z),
    np.kron(__B, -_XY),
)

# Chemical shifts
_W = [[+0.0, -1.0, +0.0],
      [+1.0, +0.0, +0.0],
      [+0.0, +0.0, +0.0]]

W, DW = (
    np.kron(_AB, _W),
    np.kron(__B, _W),
)

# Exchange rates
KAB = np.kron([[-1.0, 0.0],
               [+1.0, 0.0]], np.eye(3))

KBA = np.kron([[0.0, +1.0],
               [0.0, -1.0]], np.eye(3))

# B1 field along x
W1X = np.kron(np.eye(2), [[+0.0, +0.0, +0.0],
                          [+0.0, +0.0, -1.0],
                          [+0.0, +1.0, +0.0]])

# B1 field along y
W1Y = np.kron(np.eye(2), [[+0.0, +0.0, +1.0],
                          [+0.0, +0.0, +0.0],
                          [-1.0, +0.0, +0.0]])


def compute_liouvillian(pb=0.0, kex=0.0, r_iz=1.5, r_ixy=5.0, dr_ixy=0.0, w=0.0,
                        dw=0.0, w1x=0.0, w1y=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system. The matrix is
    written in 6x6 cartesian basis, that is {Nx, Ny, Nz}{a,b}. Note that the
    thermal equilibrium is assumed to be 0.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.
    kex : float
        Exchange rate between state A and B in /s.
    r_ixy : float
        Transverse relaxation rate of state A in /s.
    dr_ixy : float
        Transverse relaxation rate difference between states A and B in /s.
    r_iz : float
        Longitudinal relaxation rate of states A and B in /s.d
    dw: float
        Resonance offset frequency difference between states A and B in rad/s.
    w1x : float
        Angular frequency of the applied radio-frequency along x in rad/s.
    w1y: float
        Angular frequency of the applied radio-frequency along y in rad/s.

    Returns
    -------
    liouvillian : numpy.matrix
        Liouvillian of one isolated spin in presence of two-site exchange and a
        radio-frequency field.

    """

    kab = kex * pb
    kba = kex - kab

    liouvillian = (
        R_IXY * r_ixy +
        DR_IXY * dr_ixy +
        R_IZ * r_iz +
        W * w +
        DW * dw +
        W1X * w1x +
        W1Y * w1y +
        KAB * kab +
        KBA * kba
    )

    return liouvillian


def compute_propagators_from_time_series(liouvillian, times):
    s, vr = linalg.eig(liouvillian)
    vri = linalg.inv(vr)

    times_uniq = set(times)
    propagators = {
        time: dot(dot(vr, diag(exp(s * time))), vri).real
        for time in times_uniq}

    return propagators
