import numpy as np

from chemex.bases import ref

indexes = [0, 1, 2, 15, 16, 17]

mat_lambda_i_a = ref.mat_lambda_i_a[np.ix_(indexes, indexes)]
mat_lambda_i_b = ref.mat_lambda_i_b[np.ix_(indexes, indexes)]
mat_rho_i_a = ref.mat_rho_i_a[np.ix_(indexes, indexes)]
mat_rho_i_b = ref.mat_rho_i_b[np.ix_(indexes, indexes)]
mat_omega_i_a = ref.mat_omega_i_a[np.ix_(indexes, indexes)]
mat_omega_i_b = ref.mat_omega_i_b[np.ix_(indexes, indexes)]
mat_omega1x_i = ref.mat_omega1x_i[np.ix_(indexes, indexes)]
mat_omega1y_i = ref.mat_omega1y_i[np.ix_(indexes, indexes)]
mat_kab = ref.mat_kab[np.ix_(indexes, indexes)]
mat_kba = ref.mat_kba[np.ix_(indexes, indexes)]


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        lambda_i_a=0.0, rho_i_a=0.0, omega_i_a=0.0,
        lambda_i_b=0.0, rho_i_b=0.0, omega_i_b=0.0,
        omega1x_i=0.0, omega1y_i=0.0):
    """TODO: make docstring
    """

    pa = 1.0 - pb
    kab, kba = kex_ab * np.asarray([pb, pa])

    liouvillian = (
        mat_lambda_i_a * lambda_i_a +
        mat_lambda_i_b * lambda_i_b +

        mat_rho_i_a * rho_i_a +
        mat_rho_i_b * rho_i_b +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
