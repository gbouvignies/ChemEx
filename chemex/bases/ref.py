"""The ref module conatains the reference matrices and code for calculating the
Liouvillian."""

import numpy as np

pi = np.pi

mat_a = np.diag([1.0, 0.0, 0.0, 0.0])
mat_b = np.diag([0.0, 1.0, 0.0, 0.0])
mat_c = np.diag([0.0, 0.0, 1.0, 0.0])
mat_d = np.diag([0.0, 0.0, 0.0, 1.0])
mat_all = np.eye(4)

mat_r2_i = np.zeros((15, 15))
mat_r2_i[[0, 1], [0, 1]] = -1.0
mat_r2_i_a, mat_r2_i_b, mat_r2_i_c, mat_r2_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_r2_i_a[:-1, :-1] = np.kron(mat_a, mat_r2_i)
mat_r2_i_b[:-1, :-1] = np.kron(mat_b, mat_r2_i)
mat_r2_i_c[:-1, :-1] = np.kron(mat_c, mat_r2_i)
mat_r2_i_d[:-1, :-1] = np.kron(mat_d, mat_r2_i)

mat_r1_i = np.zeros((15, 15))
mat_r1_i[2, 2] = -1.0
mat_r1_i_a, mat_r1_i_b, mat_r1_i_c, mat_r1_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_r1_i_a[:-1, :-1] = np.kron(mat_a, mat_r1_i)
mat_r1_i_b[:-1, :-1] = np.kron(mat_b, mat_r1_i)
mat_r1_i_c[:-1, :-1] = np.kron(mat_c, mat_r1_i)
mat_r1_i_d[:-1, :-1] = np.kron(mat_d, mat_r1_i)

mat_r2_s = np.zeros((15, 15))
mat_r2_s[[3, 4], [3, 4]] = -1.0
mat_r2_s_a, mat_r2_s_b, mat_r2_s_c, mat_r2_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_r2_s_a[:-1, :-1] = np.kron(mat_a, mat_r2_s)
mat_r2_s_b[:-1, :-1] = np.kron(mat_b, mat_r2_s)
mat_r2_s_c[:-1, :-1] = np.kron(mat_c, mat_r2_s)
mat_r2_s_d[:-1, :-1] = np.kron(mat_d, mat_r2_s)

mat_r1_s = np.zeros((15, 15))
mat_r1_s[5, 5] = -1.0
mat_r1_s_a, mat_r1_s_b, mat_r1_s_c, mat_r1_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_r1_s_a[:-1, :-1] = np.kron(mat_a, mat_r1_s)
mat_r1_s_b[:-1, :-1] = np.kron(mat_b, mat_r1_s)
mat_r1_s_c[:-1, :-1] = np.kron(mat_c, mat_r1_s)
mat_r1_s_d[:-1, :-1] = np.kron(mat_d, mat_r1_s)

mat_r2a_i = np.zeros((15, 15))
mat_r2a_i[[6, 7], [6, 7]] = -1.0
mat_r2a_i_a, mat_r2a_i_b, mat_r2a_i_c, mat_r2a_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_r2a_i_a[:-1, :-1] = np.kron(mat_a, mat_r2a_i)
mat_r2a_i_b[:-1, :-1] = np.kron(mat_b, mat_r2a_i)
mat_r2a_i_c[:-1, :-1] = np.kron(mat_c, mat_r2a_i)
mat_r2a_i_d[:-1, :-1] = np.kron(mat_d, mat_r2a_i)

mat_r2a_s = np.zeros((15, 15))
mat_r2a_s[[8, 9], [8, 9]] = -1.0
mat_r2a_s_a, mat_r2a_s_b, mat_r2a_s_c, mat_r2a_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_r2a_s_a[:-1, :-1] = np.kron(mat_a, mat_r2a_s)
mat_r2a_s_b[:-1, :-1] = np.kron(mat_b, mat_r2a_s)
mat_r2a_s_c[:-1, :-1] = np.kron(mat_c, mat_r2a_s)
mat_r2a_s_d[:-1, :-1] = np.kron(mat_d, mat_r2a_s)

mat_r2_mq = np.zeros((15, 15))
mat_r2_mq[[10, 11, 12, 13], [10, 11, 12, 13]] = -1.0
mat_r2_mq_a, mat_r2_mq_b, mat_r2_mq_c, mat_r2_mq_d = (np.zeros((61, 61)) for _ in range(4))
mat_r2_mq_a[:-1, :-1] = np.kron(mat_a, mat_r2_mq)
mat_r2_mq_b[:-1, :-1] = np.kron(mat_b, mat_r2_mq)
mat_r2_mq_c[:-1, :-1] = np.kron(mat_c, mat_r2_mq)
mat_r2_mq_d[:-1, :-1] = np.kron(mat_d, mat_r2_mq)

mat_r1a = np.zeros((15, 15))
mat_r1a[14, 14] = -1.0
mat_r1a_a, mat_r1a_b, mat_r1a_c, mat_r1a_d = (np.zeros((61, 61)) for _ in range(4))
mat_r1a_a[:-1, :-1] = np.kron(mat_a, mat_r1a)
mat_r1a_b[:-1, :-1] = np.kron(mat_b, mat_r1a)
mat_r1a_c[:-1, :-1] = np.kron(mat_c, mat_r1a)
mat_r1a_d[:-1, :-1] = np.kron(mat_d, mat_r1a)

mat_omega_i = np.zeros((15, 15))
mat_omega_i[[1, 7, 12, 13], [0, 6, 10, 11]] = +1.0
mat_omega_i[[0, 6, 10, 11], [1, 7, 12, 13]] = -1.0
mat_omega_i_a, mat_omega_i_b, mat_omega_i_c, mat_omega_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_omega_i_a[:-1, :-1] = np.kron(mat_a, mat_omega_i)
mat_omega_i_b[:-1, :-1] = np.kron(mat_b, mat_omega_i)
mat_omega_i_c[:-1, :-1] = np.kron(mat_c, mat_omega_i)
mat_omega_i_d[:-1, :-1] = np.kron(mat_d, mat_omega_i)

mat_omega_s = np.zeros((15, 15))
mat_omega_s[[4, 9, 11, 13], [3, 8, 10, 12]] = +1.0
mat_omega_s[[3, 8, 10, 12], [4, 9, 11, 13]] = -1.0
mat_omega_s_a, mat_omega_s_b, mat_omega_s_c, mat_omega_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_omega_s_a[:-1, :-1] = np.kron(mat_a, mat_omega_s)
mat_omega_s_b[:-1, :-1] = np.kron(mat_b, mat_omega_s)
mat_omega_s_c[:-1, :-1] = np.kron(mat_c, mat_omega_s)
mat_omega_s_d[:-1, :-1] = np.kron(mat_d, mat_omega_s)

mat_j = np.zeros((15, 15))
mat_j[[7, 9, 1, 4], [0, 3, 6, 8]] = +np.pi
mat_j[[0, 3, 6, 8], [7, 9, 1, 4]] = -np.pi
mat_j_a, mat_j_b, mat_j_c, mat_j_d = (np.zeros((61, 61)) for _ in range(4))
mat_j_a[:-1, :-1] = np.kron(mat_a, mat_j)
mat_j_b[:-1, :-1] = np.kron(mat_b, mat_j)
mat_j_c[:-1, :-1] = np.kron(mat_c, mat_j)
mat_j_d[:-1, :-1] = np.kron(mat_d, mat_j)

mat_etaxy_i = np.zeros((15, 15))
mat_etaxy_i[[6, 7, 0, 1], [0, 1, 6, 7]] = -1.0
mat_etaxy_i_a, mat_etaxy_i_b, mat_etaxy_i_c, mat_etaxy_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_etaxy_i_a[:-1, :-1] = np.kron(mat_a, mat_etaxy_i)
mat_etaxy_i_b[:-1, :-1] = np.kron(mat_b, mat_etaxy_i)
mat_etaxy_i_c[:-1, :-1] = np.kron(mat_c, mat_etaxy_i)
mat_etaxy_i_d[:-1, :-1] = np.kron(mat_d, mat_etaxy_i)

mat_etaxy_s = np.zeros((15, 15))
mat_etaxy_s[[8, 9, 3, 4], [3, 4, 8, 9]] = -1.0
mat_etaxy_s_a, mat_etaxy_s_b, mat_etaxy_s_c, mat_etaxy_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_etaxy_s_a[:-1, :-1] = np.kron(mat_a, mat_etaxy_s)
mat_etaxy_s_b[:-1, :-1] = np.kron(mat_b, mat_etaxy_s)
mat_etaxy_s_c[:-1, :-1] = np.kron(mat_c, mat_etaxy_s)
mat_etaxy_s_d[:-1, :-1] = np.kron(mat_d, mat_etaxy_s)

mat_etaz_i = np.zeros((15, 15))
mat_etaz_i[[14, 2], [2, 14]] = -1.0
mat_etaz_i_a, mat_etaz_i_b, mat_etaz_i_c, mat_etaz_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_etaz_i_a[:-1, :-1] = np.kron(mat_a, mat_etaz_i)
mat_etaz_i_b[:-1, :-1] = np.kron(mat_b, mat_etaz_i)
mat_etaz_i_c[:-1, :-1] = np.kron(mat_c, mat_etaz_i)
mat_etaz_i_d[:-1, :-1] = np.kron(mat_d, mat_etaz_i)

mat_etaz_s = np.zeros((15, 15))
mat_etaz_s[[14, 5], [5, 14]] = -1.0
mat_etaz_s_a, mat_etaz_s_b, mat_etaz_s_c, mat_etaz_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_etaz_s_a[:-1, :-1] = np.kron(mat_a, mat_etaz_s)
mat_etaz_s_b[:-1, :-1] = np.kron(mat_b, mat_etaz_s)
mat_etaz_s_c[:-1, :-1] = np.kron(mat_c, mat_etaz_s)
mat_etaz_s_d[:-1, :-1] = np.kron(mat_d, mat_etaz_s)

mat_sigma = np.zeros((15, 15))
mat_sigma[[5, 0], [0, 5]] = -1.0
mat_sigma_a, mat_sigma_b, mat_sigma_c, mat_sigma_d = (np.zeros((61, 61)) for _ in range(4))
mat_sigma_a[:-1, :-1] = np.kron(mat_a, mat_sigma)
mat_sigma_b[:-1, :-1] = np.kron(mat_b, mat_sigma)
mat_sigma_c[:-1, :-1] = np.kron(mat_c, mat_sigma)
mat_sigma_d[:-1, :-1] = np.kron(mat_d, mat_sigma)

mat_mu_mq = np.zeros((15, 15))
mat_mu_mq[[13, 10], [10, 13]] = +1.0
mat_mu_mq[[12, 11], [11, 12]] = -1.0
mat_mu_mq_a, mat_mu_mq_b, mat_mu_mq_c, mat_mu_mq_d = (np.zeros((61, 61)) for _ in range(4))
mat_mu_mq_a[:-1, :-1] = np.kron(mat_a, mat_mu_mq)
mat_mu_mq_b[:-1, :-1] = np.kron(mat_b, mat_mu_mq)
mat_mu_mq_c[:-1, :-1] = np.kron(mat_c, mat_mu_mq)
mat_mu_mq_d[:-1, :-1] = np.kron(mat_d, mat_mu_mq)

mat_omega1x_i_ = np.zeros((15, 15))
mat_omega1x_i_[[2, 14, 8, 9], [1, 7, 12, 13]] = +1.0
mat_omega1x_i_[[1, 7, 12, 13], [2, 14, 8, 9]] = -1.0
mat_omega1x_i = np.zeros((61, 61))
mat_omega1x_i[:-1, :-1] = np.kron(mat_all, mat_omega1x_i_)

mat_omega1y_i_ = np.zeros((15, 15))
mat_omega1y_i_[[0, 10, 11, 6], [2, 8, 9, 14]] = +1.0
mat_omega1y_i_[[2, 8, 9, 14], [0, 10, 11, 6]] = -1.0
mat_omega1y_i = np.zeros((61, 61))
mat_omega1y_i[:-1, :-1] = np.kron(mat_all, mat_omega1y_i_)

mat_omega1x_s_ = np.zeros((15, 15))
mat_omega1x_s_[[5, 14, 6, 7], [4, 9, 11, 13]] = +1.0
mat_omega1x_s_[[4, 9, 11, 13], [5, 14, 6, 7]] = -1.0
mat_omega1x_s = np.zeros((61, 61))
mat_omega1x_s[:-1, :-1] = np.kron(mat_all, mat_omega1x_s_)

mat_omega1y_s_ = np.zeros((15, 15))
mat_omega1y_s_[[3, 10, 12, 8], [5, 6, 7, 14]] = +1.0
mat_omega1y_s_[[5, 6, 7, 14], [3, 10, 12, 8]] = -1.0
mat_omega1y_s = np.zeros((61, 61))
mat_omega1y_s[:-1, :-1] = np.kron(mat_all, mat_omega1y_s_)

mat_eq_i_a, mat_eq_i_b, mat_eq_i_c, mat_eq_i_d = (np.zeros((61, 61)) for _ in range(4))
mat_eq_i_a[2, -1] = mat_eq_i_b[17, -1] = mat_eq_i_c[32, -1] = mat_eq_i_d[47, -1] = +1.0

mat_eq_s_a, mat_eq_s_b, mat_eq_s_c, mat_eq_s_d = (np.zeros((61, 61)) for _ in range(4))
mat_eq_s_a[5, -1] = mat_eq_s_b[20, -1] = mat_eq_s_c[35, -1] = mat_eq_s_d[50, -1] = +1.0

mat_eq_is_a, mat_eq_is_b, mat_eq_is_c, mat_eq_is_d = (np.zeros((61, 61)) for _ in range(4))
mat_eq_is_a[14, -1] = mat_eq_is_b[29, -1] = mat_eq_is_c[44, -1] = mat_eq_is_d[59, -1] = +1.0


mat_kab, mat_kba, mat_kac, mat_kca, mat_kad, mat_kda, mat_kbc, mat_kcb, mat_kbd, mat_kdb, mat_kcd, mat_kdc = (
    np.zeros((61, 61)) for _ in range(12)
)

# yapf: disable

mat_kab[:-1,:-1] = np.kron([[-1.0, +0.0, +0.0, +0.0],
                            [+1.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kba[:-1,:-1] = np.kron([[+0.0, +1.0, +0.0, +0.0],
                            [+0.0, -1.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kac[:-1,:-1] = np.kron([[-1.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+1.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kca[:-1,:-1] = np.kron([[+0.0, +0.0, +1.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, -1.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kad[:-1,:-1] = np.kron([[-1.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+1.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kda[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +1.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, -1.0]], np.eye(15))

mat_kbc[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, -1.0, +0.0, +0.0],
                            [+0.0, +1.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kcb[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +1.0, +0.0],
                            [+0.0, +0.0, -1.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0]], np.eye(15))

mat_kbd[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, -1.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +1.0, +0.0, +0.0]], np.eye(15))

mat_kdb[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +1.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, -1.0]], np.eye(15))

mat_kcd[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, -1.0, +0.0],
                            [+0.0, +0.0, +1.0, +0.0]], np.eye(15))

mat_kdc[:-1,:-1] = np.kron([[+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +0.0],
                            [+0.0, +0.0, +0.0, +1.0],
                            [+0.0, +0.0, +0.0, -1.0]], np.eye(15))


def compute_liouvillian(
        kab=0.0, kba=0.0, kac=0.0, kca=0.0, kad=0.0, kda=0.0,
        kbc=0.0, kcb=0.0, kbd=0.0, kdb=0.0, kcd=0.0, kdc=0.0,
        r2_i_a=0.0, r1_i_a=0.0, r2a_i_a=0.0, etaxy_i_a=0.0, etaz_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, r2a_i_b=0.0, etaxy_i_b=0.0, etaz_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, r1_i_c=0.0, r2a_i_c=0.0, etaxy_i_c=0.0, etaz_i_c=0.0, omega_i_c=0.0,
        r2_i_d=0.0, r1_i_d=0.0, r2a_i_d=0.0, etaxy_i_d=0.0, etaz_i_d=0.0, omega_i_d=0.0,
        r2_s_a=0.0, r1_s_a=0.0, r2a_s_a=0.0, etaxy_s_a=0.0, etaz_s_a=0.0, omega_s_a=0.0,
        r2_s_b=0.0, r1_s_b=0.0, r2a_s_b=0.0, etaxy_s_b=0.0, etaz_s_b=0.0, omega_s_b=0.0,
        r2_s_c=0.0, r1_s_c=0.0, r2a_s_c=0.0, etaxy_s_c=0.0, etaz_s_c=0.0, omega_s_c=0.0,
        r2_s_d=0.0, r1_s_d=0.0, r2a_s_d=0.0, etaxy_s_d=0.0, etaz_s_d=0.0, omega_s_d=0.0,
        r2_mq_a=0.0, r1a_a=0.0, mu_mq_a=0.0, sigma_a=0.0, j_a=0.0,
        r2_mq_b=0.0, r1a_b=0.0, mu_mq_b=0.0, sigma_b=0.0, j_b=0.0,
        r2_mq_c=0.0, r1a_c=0.0, mu_mq_c=0.0, sigma_c=0.0, j_c=0.0,
        r2_mq_d=0.0, r1a_d=0.0, mu_mq_d=0.0, sigma_d=0.0, j_d=0.0,
        omega1x_i=0.0, omega1y_i=0.0, omega1x_s=0.0, omega1y_s=0.0,
        eq_i=0.0, eq_s=0.0):
    """Compute the liouvillian."""

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +
        mat_r2_i_c * r2_i_c +
        mat_r2_i_d * r2_i_d +

        mat_r2_s_a * r2_s_a +
        mat_r2_s_b * r2_s_b +
        mat_r2_s_c * r2_s_c +
        mat_r2_s_d * r2_s_d +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +
        mat_r1_i_c * r1_i_c +
        mat_r1_i_d * r1_i_d +

        mat_r1_s_a * r1_s_a +
        mat_r1_s_b * r1_s_b +
        mat_r1_s_c * r1_s_c +
        mat_r1_s_d * r1_s_d +

        mat_r2a_i_a * r2a_i_a +
        mat_r2a_i_b * r2a_i_b +
        mat_r2a_i_c * r2a_i_c +
        mat_r2a_i_d * r2a_i_d +

        mat_r2a_s_a * r2a_s_a +
        mat_r2a_s_b * r2a_s_b +
        mat_r2a_s_c * r2a_s_c +
        mat_r2a_s_d * r2a_s_d +

        mat_etaxy_i_a * etaxy_i_a +
        mat_etaxy_i_b * etaxy_i_b +
        mat_etaxy_i_c * etaxy_i_c +
        mat_etaxy_i_d * etaxy_i_d +

        mat_etaxy_s_a * etaxy_s_a +
        mat_etaxy_s_b * etaxy_s_b +
        mat_etaxy_s_c * etaxy_s_c +
        mat_etaxy_s_d * etaxy_s_d +

        mat_etaz_i_a * etaz_i_a +
        mat_etaz_i_b * etaz_i_b +
        mat_etaz_i_c * etaz_i_c +
        mat_etaz_i_d * etaz_i_d +

        mat_etaz_s_a * etaz_s_a +
        mat_etaz_s_b * etaz_s_b +
        mat_etaz_s_c * etaz_s_c +
        mat_etaz_s_d * etaz_s_d +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +
        mat_omega_i_d * omega_i_d +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +
        mat_omega_s_c * omega_s_c +
        mat_omega_s_d * omega_s_d +

        mat_r1a_a * r1a_a +
        mat_r1a_b * r1a_b +
        mat_r1a_c * r1a_c +
        mat_r1a_d * r1a_d +

        mat_r2_mq_a * r2_mq_a +
        mat_r2_mq_b * r2_mq_b +
        mat_r2_mq_c * r2_mq_c +
        mat_r2_mq_d * r2_mq_d +

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

        mat_eq_i_a * (r1_i_a * eq_i + sigma_a * eq_s) +
        mat_eq_i_b * (r1_i_b * eq_i + sigma_b * eq_s) +
        mat_eq_i_c * (r1_i_c * eq_i + sigma_c * eq_s) +
        mat_eq_i_d * (r1_i_d * eq_i + sigma_d * eq_s) +

        mat_eq_s_a * (r1_s_a * eq_s + sigma_a * eq_i) +
        mat_eq_s_b * (r1_s_b * eq_s + sigma_b * eq_i) +
        mat_eq_s_c * (r1_s_c * eq_s + sigma_c * eq_i) +
        mat_eq_s_d * (r1_s_d * eq_s + sigma_d * eq_i) +

        mat_eq_is_a * (etaz_i_a * eq_i + etaz_s_a * eq_s) +
        mat_eq_is_b * (etaz_i_b * eq_i + etaz_s_b * eq_s) +
        mat_eq_is_c * (etaz_i_c * eq_i + etaz_s_c * eq_s) +
        mat_eq_is_d * (etaz_i_d * eq_i + etaz_s_d * eq_s) +

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
