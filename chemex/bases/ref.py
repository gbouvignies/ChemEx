import numpy as np

pi = np.pi


# Operator basis:
# {Ix, Iy, Iz, Sx, Sy, Sz,
#  2IxSz, 2IySz, 2IzSx, 2IzSy,
#  2IxSx, 2IxSy, 2IySx, 2IySy,
#  2IzSz}

mat_a = np.diag([1.0, 0.0, 0.0, 0.0])
mat_b = np.diag([0.0, 1.0, 0.0, 0.0])
mat_c = np.diag([0.0, 0.0, 1.0, 0.0])
mat_d = np.diag([0.0, 0.0, 0.0, 1.0])
mat_all = np.eye(4)

mat_lambda_i = np.zeros((15, 15))
mat_lambda_i[[0, 1], [0, 1]] = -1.0
mat_lambda_i_a = np.kron(mat_a, mat_lambda_i)
mat_lambda_i_b = np.kron(mat_b, mat_lambda_i)
mat_lambda_i_c = np.kron(mat_c, mat_lambda_i)
mat_lambda_i_d = np.kron(mat_d, mat_lambda_i)

mat_rho_i = np.zeros((15, 15))
mat_rho_i[2, 2] = -1.0
mat_rho_i_a = np.kron(mat_a, mat_rho_i)
mat_rho_i_b = np.kron(mat_b, mat_rho_i)
mat_rho_i_c = np.kron(mat_c, mat_rho_i)
mat_rho_i_d = np.kron(mat_d, mat_rho_i)

mat_lambda_s = np.zeros((15, 15))
mat_lambda_s[[3, 4], [3, 4]] = -1.0
mat_lambda_s_a = np.kron(mat_a, mat_lambda_s)
mat_lambda_s_b = np.kron(mat_b, mat_lambda_s)
mat_lambda_s_c = np.kron(mat_c, mat_lambda_s)
mat_lambda_s_d = np.kron(mat_d, mat_lambda_s)

mat_rho_s = np.zeros((15, 15))
mat_rho_s[5, 5] = -1.0
mat_rho_s_a = np.kron(mat_a, mat_rho_s)
mat_rho_s_b = np.kron(mat_b, mat_rho_s)
mat_rho_s_c = np.kron(mat_c, mat_rho_s)
mat_rho_s_d = np.kron(mat_d, mat_rho_s)

mat_rhoa_i = np.zeros((15, 15))
mat_rhoa_i[[6, 7], [6, 7]] = -1.0
mat_rhoa_i_a = np.kron(mat_a, mat_rhoa_i)
mat_rhoa_i_b = np.kron(mat_b, mat_rhoa_i)
mat_rhoa_i_c = np.kron(mat_c, mat_rhoa_i)
mat_rhoa_i_d = np.kron(mat_d, mat_rhoa_i)

mat_rhoa_s = np.zeros((15, 15))
mat_rhoa_s[[8, 9], [8, 9]] = -1.0
mat_rhoa_s_a = np.kron(mat_a, mat_rhoa_s)
mat_rhoa_s_b = np.kron(mat_b, mat_rhoa_s)
mat_rhoa_s_c = np.kron(mat_c, mat_rhoa_s)
mat_rhoa_s_d = np.kron(mat_d, mat_rhoa_s)

mat_lambda_mq = np.zeros((15, 15))
mat_lambda_mq[[10, 11, 12, 13], [10, 11, 12, 13]] = -1.0
mat_lambda_mq_a = np.kron(mat_a, mat_lambda_mq)
mat_lambda_mq_b = np.kron(mat_b, mat_lambda_mq)
mat_lambda_mq_c = np.kron(mat_c, mat_lambda_mq)
mat_lambda_mq_d = np.kron(mat_d, mat_lambda_mq)

mat_rho_is = np.zeros((15, 15))
mat_rho_is[14, 14] = -1.0
mat_rho_is_a = np.kron(mat_a, mat_rho_is)
mat_rho_is_b = np.kron(mat_b, mat_rho_is)
mat_rho_is_c = np.kron(mat_c, mat_rho_is)
mat_rho_is_d = np.kron(mat_d, mat_rho_is)

mat_omega_i = np.zeros((15, 15))
mat_omega_i[[1, 7, 12, 13], [0, 6, 10, 11]] = +1.0
mat_omega_i[[0, 6, 10, 11], [1, 7, 12, 13]] = -1.0
mat_omega_i_a = np.kron(mat_a, mat_omega_i)
mat_omega_i_b = np.kron(mat_b, mat_omega_i)
mat_omega_i_c = np.kron(mat_c, mat_omega_i)
mat_omega_i_d = np.kron(mat_d, mat_omega_i)

mat_omega_s = np.zeros((15, 15))
mat_omega_s[[4, 9, 11, 13], [3, 8, 10, 12]] = +1.0
mat_omega_s[[3, 8, 10, 12], [4, 9, 11, 13]] = -1.0
mat_omega_s_a = np.kron(mat_a, mat_omega_s)
mat_omega_s_b = np.kron(mat_b, mat_omega_s)
mat_omega_s_c = np.kron(mat_c, mat_omega_s)
mat_omega_s_d = np.kron(mat_d, mat_omega_s)

mat_j = np.zeros((15, 15))
mat_j[[7, 9, 1, 4], [0, 3, 6, 8]] = +np.pi
mat_j[[0, 3, 6, 8], [7, 9, 1, 4]] = -np.pi
mat_j_a = np.kron(mat_a, mat_j)
mat_j_b = np.kron(mat_b, mat_j)
mat_j_c = np.kron(mat_c, mat_j)
mat_j_d = np.kron(mat_d, mat_j)

mat_eta_i = np.zeros((15, 15))
mat_eta_i[[6, 7, 0, 1], [0, 1, 6, 7]] = -1.0
mat_eta_i_a = np.kron(mat_a, mat_eta_i)
mat_eta_i_b = np.kron(mat_b, mat_eta_i)
mat_eta_i_c = np.kron(mat_c, mat_eta_i)
mat_eta_i_d = np.kron(mat_d, mat_eta_i)

mat_eta_s = np.zeros((15, 15))
mat_eta_s[[8, 9, 3, 4], [3, 4, 8, 9]] = -1.0
mat_eta_s_a = np.kron(mat_a, mat_eta_s)
mat_eta_s_b = np.kron(mat_b, mat_eta_s)
mat_eta_s_c = np.kron(mat_c, mat_eta_s)
mat_eta_s_d = np.kron(mat_d, mat_eta_s)

mat_delta_i = np.zeros((15, 15))
mat_delta_i[[14, 2], [2, 14]] = -1.0
mat_delta_i_a = np.kron(mat_a, mat_delta_i)
mat_delta_i_b = np.kron(mat_b, mat_delta_i)
mat_delta_i_c = np.kron(mat_c, mat_delta_i)
mat_delta_i_d = np.kron(mat_d, mat_delta_i)

mat_delta_s = np.zeros((15, 15))
mat_delta_s[[14, 5], [5, 14]] = -1.0
mat_delta_s_a = np.kron(mat_a, mat_delta_s)
mat_delta_s_b = np.kron(mat_b, mat_delta_s)
mat_delta_s_c = np.kron(mat_c, mat_delta_s)
mat_delta_s_d = np.kron(mat_d, mat_delta_s)

mat_sigma = np.zeros((15, 15))
mat_sigma[[5, 0], [0, 5]] = -1.0
mat_sigma_a = np.kron(mat_a, mat_sigma)
mat_sigma_b = np.kron(mat_b, mat_sigma)
mat_sigma_c = np.kron(mat_c, mat_sigma)
mat_sigma_d = np.kron(mat_d, mat_sigma)

mat_mu_mq = np.zeros((15, 15))
mat_mu_mq[[13, 10], [10, 13]] = +1.0
mat_mu_mq[[12, 11], [11, 12]] = -1.0
mat_mu_mq_a = np.kron(mat_a, mat_mu_mq)
mat_mu_mq_b = np.kron(mat_b, mat_mu_mq)
mat_mu_mq_c = np.kron(mat_c, mat_mu_mq)
mat_mu_mq_d = np.kron(mat_d, mat_mu_mq)

mat_omega1x_i = np.zeros((15, 15))
mat_omega1x_i[[2, 14, 8, 9], [1, 7, 12, 13]] = +1.0
mat_omega1x_i[[1, 7, 12, 13], [2, 14, 8, 9]] = -1.0
mat_omega1x_i = np.kron(mat_all, mat_omega1x_i)

mat_omega1y_i = np.zeros((15, 15))
mat_omega1y_i[[0, 10, 11, 6], [2, 8, 9, 14]] = +1.0
mat_omega1y_i[[2, 8, 9, 14], [0, 10, 11, 6]] = -1.0
mat_omega1y_i = np.kron(mat_all, mat_omega1y_i)

mat_omega1x_s = np.zeros((15, 15))
mat_omega1x_s[[5, 14, 6, 7], [4, 9, 11, 13]] = +1.0
mat_omega1x_s[[4, 9, 11, 13], [5, 14, 6, 7]] = -1.0
mat_omega1x_s = np.kron(mat_all, mat_omega1x_s)

mat_omega1y_s = np.zeros((15, 15))
mat_omega1y_s[[3, 10, 12, 8], [5, 6, 7, 14]] = +1.0
mat_omega1y_s[[5, 6, 7, 14], [3, 10, 12, 8]] = -1.0
mat_omega1y_s = np.kron(mat_all, mat_omega1y_s)

mat_kab = np.kron([[-1.0, +0.0, +0.0, +0.0],
                   [+1.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kba = np.kron([[+0.0, +1.0, +0.0, +0.0],
                   [+0.0, -1.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kac = np.kron([[-1.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+1.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kca = np.kron([[+0.0, +0.0, +1.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, -1.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kad = np.kron([[-1.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+1.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kda = np.kron([[+0.0, +0.0, +0.0, +1.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, -1.0]], np.eye(15))

mat_kbc = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, -1.0, +0.0, +0.0],
                   [+0.0, +1.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kcb = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +1.0, +0.0],
                   [+0.0, +0.0, -1.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kbd = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, -1.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +1.0, +0.0, +0.0]], np.eye(15))

mat_kdb = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +1.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, -1.0]], np.eye(15))

mat_kcd = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, -1.0, +0.0],
                   [+0.0, +0.0, +1.0, +0.0]], np.eye(15))

mat_kdc = np.kron([[+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +0.0],
                   [+0.0, +0.0, +0.0, +1.0],
                   [+0.0, +0.0, +0.0, -1.0]], np.eye(15))

def compute_liouvillian(
        pb, pc, pd, kex_ab, kex_ac, kex_ad, kex_bc, kex_bd, kex_cd,
        lambda_i_a=0.0, rho_i_a=0.0, rhoa_i_a=0.0, eta_i_a=0.0, delta_i_a=0.0, omega_i_a=0.0,
        lambda_i_b=0.0, rho_i_b=0.0, rhoa_i_b=0.0, eta_i_b=0.0, delta_i_b=0.0, omega_i_b=0.0,
        lambda_i_c=0.0, rho_i_c=0.0, rhoa_i_c=0.0, eta_i_c=0.0, delta_i_c=0.0, omega_i_c=0.0,
        lambda_i_d=0.0, rho_i_d=0.0, rhoa_i_d=0.0, eta_i_d=0.0, delta_i_d=0.0, omega_i_d=0.0,
        lambda_s_a=0.0, rho_s_a=0.0, rhoa_s_a=0.0, eta_s_a=0.0, delta_s_a=0.0, omega_s_a=0.0,
        lambda_s_b=0.0, rho_s_b=0.0, rhoa_s_b=0.0, eta_s_b=0.0, delta_s_b=0.0, omega_s_b=0.0,
        lambda_s_c=0.0, rho_s_c=0.0, rhoa_s_c=0.0, eta_s_c=0.0, delta_s_c=0.0, omega_s_c=0.0,
        lambda_s_d=0.0, rho_s_d=0.0, rhoa_s_d=0.0, eta_s_d=0.0, delta_s_d=0.0, omega_s_d=0.0,
        lambda_mq_a=0.0, rho_is_a=0.0, mu_mq_a=0.0, sigma_a=0.0, j_a=0.0,
        lambda_mq_b=0.0, rho_is_b=0.0, mu_mq_b=0.0, sigma_b=0.0, j_b=0.0,
        lambda_mq_c=0.0, rho_is_c=0.0, mu_mq_c=0.0, sigma_c=0.0, j_c=0.0,
        lambda_mq_d=0.0, rho_is_d=0.0, mu_mq_d=0.0, sigma_d=0.0, j_d=0.0,
        omega1x_i=0.0, omega1y_i=0.0, omega1x_s=0.0, omega1y_s=0.0):
    """TODO: make docstringx
    """

    pa = 1.0 - pb - pc - pd

    kex_ab_ = kex_ab / (pa + pb)
    kex_ac_ = kex_ac / (pa + pc)
    kex_ad_ = kex_ad / (pa + pd)
    kex_bc_ = kex_bc / (pb + pc)
    kex_bd_ = kex_bd / (pb + pd)
    kex_cd_ = kex_cd / (pc + pd)

    kab, kba = kex_ab_ * np.asarray([pb, pa])
    kac, kca = kex_ac_ * np.asarray([pc, pa])
    kad, kda = kex_ad_ * np.asarray([pd, pa])
    kbc, kcb = kex_bc_ * np.asarray([pc, pb])
    kbd, kdb = kex_bd_ * np.asarray([pd, pb])
    kcd, kdc = kex_cd_ * np.asarray([pd, pc])

    liouvillian = (
        mat_lambda_i_a * lambda_i_a +
        mat_lambda_i_b * lambda_i_b +
        mat_lambda_i_c * lambda_i_c +
        mat_lambda_i_d * lambda_i_d +

        mat_lambda_s_a * lambda_s_a +
        mat_lambda_s_b * lambda_s_b +
        mat_lambda_s_c * lambda_s_c +
        mat_lambda_s_d * lambda_s_d +

        mat_rho_i_a * rho_i_a +
        mat_rho_i_b * rho_i_b +
        mat_rho_i_c * rho_i_c +
        mat_rho_i_d * rho_i_d +

        mat_rho_s_a * rho_s_a +
        mat_rho_s_b * rho_s_b +
        mat_rho_s_c * rho_s_c +
        mat_rho_s_d * rho_s_d +

        mat_rhoa_i_a * rhoa_i_a +
        mat_rhoa_i_b * rhoa_i_b +
        mat_rhoa_i_c * rhoa_i_c +
        mat_rhoa_i_d * rhoa_i_d +

        mat_rhoa_s_a * rhoa_s_a +
        mat_rhoa_s_b * rhoa_s_b +
        mat_rhoa_s_c * rhoa_s_c +
        mat_rhoa_s_d * rhoa_s_d +

        mat_eta_i_a * eta_i_a +
        mat_eta_i_b * eta_i_b +
        mat_eta_i_c * eta_i_c +
        mat_eta_i_d * eta_i_d +

        mat_eta_s_a * eta_s_a +
        mat_eta_s_b * eta_s_b +
        mat_eta_s_c * eta_s_c +
        mat_eta_s_d * eta_s_d +

        mat_delta_i_a * delta_i_a +
        mat_delta_i_b * delta_i_b +
        mat_delta_i_c * delta_i_c +
        mat_delta_i_d * delta_i_d +

        mat_delta_s_a * delta_s_a +
        mat_delta_s_b * delta_s_b +
        mat_delta_s_c * delta_s_c +
        mat_delta_s_d * delta_s_d +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +
        mat_omega_i_d * omega_i_d +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +
        mat_omega_s_c * omega_s_c +
        mat_omega_s_d * omega_s_d +

        mat_rho_is_a * rho_is_a +
        mat_rho_is_b * rho_is_b +
        mat_rho_is_c * rho_is_c +
        mat_rho_is_d * rho_is_d +

        mat_lambda_mq_a * lambda_mq_a +
        mat_lambda_mq_b * lambda_mq_b +
        mat_lambda_mq_c * lambda_mq_c +
        mat_lambda_mq_d * lambda_mq_d +

        mat_mu_mq_a * mu_mq_a +
        mat_mu_mq_b * mu_mq_b +
        mat_mu_mq_c * mu_mq_c +
        mat_mu_mq_d * mu_mq_d +

        mat_sigma_a * sigma_a +
        mat_sigma_b * sigma_b +
        mat_sigma_c * sigma_c +
        mat_sigma_d * sigma_d +

        mat_j_a * j_a +
        mat_j_b * j_b +
        mat_j_c * j_c +
        mat_j_d * j_d +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +
        mat_omega1x_s * omega1x_s +
        mat_omega1y_s * omega1y_s +

        mat_kab * kab +
        mat_kba * kba +
        mat_kac * kac +
        mat_kca * kca +
        mat_kad * kad +
        mat_kda * kda +
        mat_kbc * kbc +
        mat_kcb * kcb +
        mat_kbd * kbd +
        mat_kdb * kdb +
        mat_kcd * kcd +
        mat_kdc * kdc
    )

    return liouvillian
