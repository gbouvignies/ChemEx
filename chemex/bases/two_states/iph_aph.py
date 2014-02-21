"""
Created on 2013-10-03

@author: guillaume
"""

# Imports
from scipy import eye, kron, diag


# Define the basis for the liouvillian
# Axes: _XY, __Z
_XY, __Z = (diag([1.0, 1.0, 0.0]),
            diag([0.0, 0.0, 1.0]))

# States: B or both A & B
__B, _AB = (diag([0.0, 1.0]),
            diag([1.0, 1.0]))

# Coherences: inphase,  antiphase or both
__I, __A, _IA = (diag([1.0, 0.0]),
                 diag([0.0, 1.0]),
                 diag([1.0, 1.0]))

# Auto-relaxation rates
R_IXY, R_2SZIXY, R_IZ, R_2SZIZ, DR_XY = (kron(_AB, -kron(__I, _XY)),
                                         kron(_AB, -kron(__A, _XY)),
                                         kron(_AB, -kron(__I, __Z)),
                                         kron(_AB, -kron(__A, __Z)),
                                         kron(__B, -kron(_IA, _XY)))

# Chemical shifts & Scalar couplings
_CS = [[+0.0, -1.0, +0.0],
       [+1.0, +0.0, +0.0],
       [+0.0, +0.0, +0.0]]

__J = [[+0.0, +1.0],
       [+1.0, +0.0]]

CS, DW, J, DJ = (kron(_AB, kron(_IA, _CS)),
                 kron(__B, kron(_IA, _CS)),
                 kron(_AB, kron(__J, _CS)),
                 kron(__B, kron(__J, _CS)))

# Cross-correlated relaxation rates
_ETA = [[+0.0, -1.0],
        [-1.0, +0.0]]

ETAXY, ETAZ = (kron(_AB, kron(_ETA, _XY)),
               kron(_AB, kron(_ETA, __Z)))

# Exchange rates
KAB = kron([[-1.0, 0.0],
            [+1.0, 0.0]], eye(6))

KBA = kron([[0.0, +1.0],
            [0.0, -1.0]], eye(6))

# B1 field along x
W1X = kron(_AB, kron(_IA, [[+0.0, +0.0, +0.0],
                           [+0.0, +0.0, -1.0],
                           [+0.0, +1.0, +0.0]]))

# B1 field along y
W1Y = kron(_AB, kron(_IA, [[+0.0, +0.0, +1.0],
                           [+0.0, +0.0, +0.0],
                           [-1.0, +0.0, +0.0]]))

# 180 degree pulse on S
P180_S = kron(eye(2), kron(diag([+1.0, -1.0]), eye(3)))


