'''
Created on 2013-10-03

@author: guillaume
'''

# Imports
from scipy import eye, kron, diag


# Define the basis for the liouvillian
# Axes: _XY, __Z
_XY, __Z = (diag([1.0, 1.0, 0.0]),
            diag([0.0, 0.0, 1.0]))

# States: B or both A & B
__B, _AB = (diag([0.0, 1.0]),
            diag([1.0, 1.0]))

# Auto-relaxation rates
R_IXY, R_IZ, DR_IXY = (kron(_AB, -_XY),
                       kron(_AB, -__Z),
                       kron(__B, -_XY))
# Chemical shifts
_CS = [[+0.0, -1.0, +0.0],
       [+1.0, +0.0, +0.0],
       [+0.0, +0.0, +0.0]]

CS, DW = (kron(_AB, _CS),
          kron(__B, _CS))

# Exchange rates
KAB = kron([[-1.0, 0.0],
            [+1.0, 0.0]], eye(3))

KBA = kron([[0.0, +1.0],
            [0.0, -1.0]], eye(3))

# B1 field along x
W1X = kron(eye(2), [[+0.0, +0.0, +0.0],
                    [+0.0, +0.0, -1.0],
                    [+0.0, +1.0, +0.0]])

# B1 field along y
W1Y = kron(eye(2), [[+0.0, +0.0, +1.0],
                    [+0.0, +0.0, +0.0],
                    [-1.0, +0.0, +0.0]])


