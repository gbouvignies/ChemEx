"""
Created on 2013-10-03

@author: guillaume
"""

# Imports
from scipy import eye, kron, diag


# Define the basis for the liouvillian
# States: B, C or all states
___B, ___C, _ABC = (
    diag([0.0, 1.0, 0.0]),
    diag([0.0, 0.0, 1.0]),
    diag([1.0, 1.0, 1.0]),
)

# Auto-relaxation rates
R_IXY, DR_IXY_AB, DR_IXY_AC = (
    kron(_ABC, -eye(2)),
    kron(___B, -eye(2)),
    kron(___C, -eye(2)),
)
# Chemical shifts
_CS = [[+0.0, -1.0],
       [+1.0, +0.0]]

DW_AB = kron(___B, _CS)
DW_AC = kron(___C, _CS)

# Exchange rates
KAB = kron([[-1.0, +0.0, +0.0],
            [+1.0, +0.0, +0.0],
            [+0.0, +0.0, +0.0]], eye(2))

KBA = kron([[+0.0, +1.0, +0.0],
            [+0.0, -1.0, +0.0],
            [+0.0, +0.0, +0.0]], eye(2))

KBC = kron([[+0.0, +0.0, +0.0],
            [+0.0, -1.0, +0.0],
            [+0.0, +1.0, +0.0]], eye(2))

KCB = kron([[+0.0, +0.0, +0.0],
            [+0.0, +0.0, +1.0],
            [+0.0, +0.0, -1.0]], eye(2))

KAC = kron([[-1.0, +0.0, +0.0],
            [+0.0, +0.0, +0.0],
            [+1.0, +0.0, +0.0]], eye(2))

KCA = kron([[+0.0, +0.0, +1.0],
            [+0.0, +0.0, +0.0],
            [+0.0, +0.0, -1.0]], eye(2))


# 180 degree y pulse
P_180Y = diag([1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
