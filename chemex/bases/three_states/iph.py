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

# States: B, C or A & B & C
___B, ___C, _ABC = (diag([0.0, 1.0, 0.0]),
                    diag([0.0, 0.0, 1.0]),
                    diag([1.0, 1.0, 1.0]))

# Auto-relaxation rates
R_IXY, R_IZ, DR_IXY_AB, DR_IXY_AC = (kron(_ABC, -_XY),
                                     kron(_ABC, -__Z),
                                     kron(___B, -_XY),
                                     kron(___B, -_XY))
# Chemical shifts
_CS = [[+0.0, -1.0, +0.0],
       [+1.0, +0.0, +0.0],
       [+0.0, +0.0, +0.0]]

CS, DW_AB, DW_AC = (kron(_ABC, _CS),
                    kron(___B, _CS),
                    kron(___C, _CS))

# Exchange rates
KAB = kron([[-1.0, +0.0, +0.0],
            [+1.0, +0.0, +0.0],
            [+0.0, +0.0, +0.0]], eye(3))

KBA = kron([[+0.0, +1.0, +0.0],
            [+0.0, -1.0, +0.0],
            [+0.0, +0.0, +0.0]], eye(3))

KBC = kron([[+0.0, +0.0, +0.0],
            [+0.0, -1.0, +0.0],
            [+0.0, +1.0, +0.0]], eye(3))

KCB = kron([[+0.0, +0.0, +0.0],
            [+0.0, +0.0, +1.0],
            [+0.0, +0.0, -1.0]], eye(3))

KAC = kron([[-1.0, +0.0, +0.0],
            [+0.0, +0.0, +0.0],
            [+1.0, +0.0, +0.0]], eye(3))

KCA = kron([[+0.0, +0.0, +1.0],
            [+0.0, +0.0, +0.0],
            [+0.0, +0.0, -1.0]], eye(3))

# B1 field along x
W1X = kron(_ABC, [[+0.0, +0.0, +0.0],
                  [+0.0, +0.0, -1.0],
                  [+0.0, +1.0, +0.0]])
