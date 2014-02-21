"""
Created on 2013-10-03

@author: guillaume
"""

# Imports
from scipy import eye, kron, diag


# Define the basis for the liouvillian
# States: B or both A & B
__B, _AB = (
    diag([0.0, 1.0]),
    diag([1.0, 1.0]),
)

# Auto-relaxation rates
R_IXY, DR_IXY = (
    kron(_AB, -eye(2)),
    kron(__B, -eye(2)),
)
# Chemical shifts
_CS = [[+0.0, -1.0],
       [+1.0, +0.0]]

DW = kron(__B, _CS)

# Exchange rates
KAB = kron([[-1.0, 0.0],
            [+1.0, 0.0]], eye(2))

KBA = kron([[0.0, +1.0],
            [0.0, -1.0]], eye(2))

# 180 degree y pulse
P_180Y = diag([1.0, -1.0, 1.0, -1.0])
