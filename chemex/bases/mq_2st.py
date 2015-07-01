import numpy as np

from chemex.bases import ref

pi = np.pi


# Operator basis:
# {2IxSx, 2IxSy, 2IySx, 2IySy}{a, b}

indexes = [10, 11, 12, 13,
           25, 26, 27, 28, ]

mat_omega_i_a = ref.mat_omega_i_a[np.ix_(indexes, indexes)]
mat_omega_i_b = ref.mat_omega_i_b[np.ix_(indexes, indexes)]
mat_omega_s_a = ref.mat_omega_s_a[np.ix_(indexes, indexes)]
mat_omega_s_b = ref.mat_omega_s_b[np.ix_(indexes, indexes)]
mat_lambda_mq_a = ref.mat_lambda_mq_a[np.ix_(indexes, indexes)]
mat_lambda_mq_b = ref.mat_lambda_mq_b[np.ix_(indexes, indexes)]
mat_mu_mq_a = ref.mat_mu_mq_a[np.ix_(indexes, indexes)]
mat_mu_mq_b = ref.mat_mu_mq_b[np.ix_(indexes, indexes)]
mat_kab = ref.mat_kab[np.ix_(indexes, indexes)]
mat_kba = ref.mat_kba[np.ix_(indexes, indexes)]
mat_kbc = ref.mat_kbc[np.ix_(indexes, indexes)]
mat_kcb = ref.mat_kcb[np.ix_(indexes, indexes)]
mat_kac = ref.mat_kac[np.ix_(indexes, indexes)]
mat_kca = ref.mat_kca[np.ix_(indexes, indexes)]

p_180x_i = np.diag([1.0, 1.0, -1.0, -1.0,
                    1.0, 1.0, -1.0, -1.0, ])
p_180y_i = np.diag([-1.0, -1.0, 1.0, 1.0,
                    -1.0, -1.0, 1.0, 1.0, ])
p_180x_s = np.diag([1.0, -1.0, 1.0, -1.0,
                    1.0, -1.0, 1.0, -1.0, ])
p_180y_s = np.diag([-1.0, 1.0, -1.0, 1.0,
                    -1.0, 1.0, -1.0, 1.0, ])


def compute_liouvillian(
        pb, kex_ab,
        omega_i_a=0.0, omega_s_a=0.0,
        omega_i_b=0.0, omega_s_b=0.0,
        lambda_mq_a=0.0, mu_mq_a=0.0,
        lambda_mq_b=0.0, mu_mq_b=0.0, ):
    """TODO: make docstringx
    """

    pa = 1.0 - pb
    kab, kba = kex_ab * np.asarray([pb, pa])

    liouvillian = (
        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +

        mat_lambda_mq_a * lambda_mq_a +
        mat_lambda_mq_b * lambda_mq_b +

        mat_mu_mq_a * mu_mq_a +
        mat_mu_mq_b * mu_mq_b +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
