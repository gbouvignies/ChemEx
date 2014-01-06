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
R_MQ, DR_MQ = (
    kron(_AB, -eye(4)),
    kron(__B, -eye(4)),
)
# Chemical shifts
_CS = [[+0.0, -1.0],
       [+1.0, +0.0]]

DWI = kron(__B, kron(eye(2), _CS))
DWS = kron(__B, kron(_CS, eye(2)))

# Exchange rates
KAB = kron([[-1.0, 0.0],
            [+1.0, 0.0]], eye(4))

KBA = kron([[0.0, +1.0],
            [0.0, -1.0]], eye(4))

