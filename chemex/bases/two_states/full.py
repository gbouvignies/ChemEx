from scipy import array, zeros, ones, kron, eye

# Temporary matrices to help build the liouvillian basis
TMP1 = array([[+0.0, -1.0],
              [+1.0, +0.0]])

TMP2 = array([[+0.0, +0.0, -1.0],
              [+0.0, +0.0, +0.0],
              [+1.0, +0.0, +0.0]])

# Define the basis for the liouvillian
(
    R_HXY,
    R_HZ,
    R_NXY_A, R_NXY_B,
    R_NZ,
    R_2HXYNZ,
    R_2HZNXY_A,
    R_2HZNXY_B,
    R_2HXYNXY,
    R_2HZNZ,
    CS_H_A,
    CS_H_B,
    CS_N_A,
    CS_N_B,
    J_HN,
    ETAZ,
    ETAXY,
    W1X_H,
    W1Y_H,
    W1X_N,
    W1Y_N,
) = (zeros((30, 30)) for _ in range(21))

# State A
# Auto-relaxation rates
R_HXY[0:2, 0:2] -= eye(2)
R_HZ[2, 2] -= 1.0
R_NXY_A[3:5, 3:5] -= eye(2)
R_NZ[5, 5] -= 1.0
R_2HXYNZ[6:8, 6:8] -= eye(2)
R_2HZNXY_A[8:10, 8:10] -= eye(2)
R_2HXYNXY[10:14, 10:14] -= eye(4)
R_2HZNZ[14, 14] -= 1.0
#
# Chemical shifts
CS_H_A[0:2, 0:2] += TMP1  # {Hx,Hy}
CS_H_A[6:8, 6:8] += TMP1  # {2HxNz,2HyNz}
CS_H_A[10:13, 10:13] += TMP2  # {2HxNx,2HyNx}
CS_H_A[11:14, 11:14] += TMP2  # {2HxNy,2HyNy}
#
CS_N_A[3:5, 3:5] += TMP1  # {Nx,Ny}
CS_N_A[8:10, 8:10] += TMP1  # {2HzNx,2HzSy}
CS_N_A[10:12, 10:12] += TMP1  # {2HxNx,2HxNy}
CS_N_A[12:14, 12:14] += TMP1  # {2HyNx,2HyNy}
#
# Scalar coupling
J_HN[0:2, 6:8] += TMP1
J_HN[6:8, 0:2] += TMP1
J_HN[3:5, 8:10] += TMP1
J_HN[8:10, 3:5] += TMP1
#
# Transverse dipole-CSA cross-relaxation
ETAXY[3:5, 8:10] -= eye(2)
ETAXY[8:10, 3:5] -= eye(2)
#
# Longitudinal dipole-CSA cross-relaxation
ETAZ[5, 14] -= 1.0
ETAZ[14, 5] -= 1.0

#
# State B
R_HXY[15:30, 15:30] += R_HXY[0:15, 0:15]
R_HZ[15:30, 15:30] += R_HZ[0:15, 0:15]
R_NXY_B[15:30, 15:30] += R_NXY_A[0:15, 0:15]
R_NZ[15:30, 15:30] += R_NZ[0:15, 0:15]
R_2HXYNZ[15:30, 15:30] += R_2HXYNZ[0:15, 0:15]
R_2HZNXY_B[15:30, 15:30] += R_2HZNXY_A[0:15, 0:15]
R_2HXYNXY[15:30, 15:30] += R_2HXYNXY[0:15, 0:15]
R_2HZNZ[15:30, 15:30] += R_2HZNZ[0:15, 0:15]
CS_H_B[15:30, 15:30] += CS_H_A[0:15, 0:15]
CS_N_B[15:30, 15:30] += CS_N_A[0:15, 0:15]
J_HN[15:30, 15:30] += J_HN[0:15, 0:15]
ETAZ[15:30, 15:30] += ETAZ[0:15, 0:15]
ETAXY[15:30, 15:30] += ETAXY[0:15, 0:15]

# Exchange rates
KAB = kron([[-1.0, 0.0],
            [+1.0, 0.0]], eye(15))

KBA = kron([[0.0, +1.0],
            [0.0, -1.0]], eye(15))


# 1HN B1 field along x
W1X_H[1, 2], W1X_H[7, 14], W1X_H[13, 9], W1X_H[12, 8] = -ones(4)
W1X_H[2, 1], W1X_H[14, 7], W1X_H[9, 13], W1X_H[8, 12] = +ones(4)
W1X_H[15:30, 15:30] += W1X_H[0:15, 0:15]

# 1HN B1 field along y
W1Y_H[0, 2], W1Y_H[6, 14], W1Y_H[10, 8], W1Y_H[11, 9] = +ones(4)
W1Y_H[2, 0], W1Y_H[14, 6], W1Y_H[8, 10], W1Y_H[9, 11] = -ones(4)
W1Y_H[15:30, 15:30] += W1Y_H[0:15, 0:15]

# 15N B1 field along x
W1X_N[4, 5], W1X_N[11, 6], W1X_N[13, 7], W1X_N[9, 14] = -ones(4)
W1X_N[5, 4], W1X_N[6, 11], W1X_N[7, 13], W1X_N[14, 9] = +ones(4)
W1X_N[15:30, 15:30] += W1X_N[0:15, 0:15]

# 15N B1 field along y
W1Y_N[3, 5], W1Y_N[10, 6], W1Y_N[12, 7], W1Y_N[8, 14] = +ones(4)
W1Y_N[5, 3], W1Y_N[6, 10], W1Y_N[7, 12], W1Y_N[14, 8] = -ones(4)
W1Y_N[15:30, 15:30] += W1Y_N[0:15, 0:15]

#
# Some cleaning
del (TMP1)
del (TMP2)
