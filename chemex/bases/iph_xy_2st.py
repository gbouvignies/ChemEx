import numpy as np

from chemex.bases import ref

indexes = [0, 1, 15, 16]

mat_lambda_i_a = ref.mat_lambda_i_a[np.ix_(indexes, indexes)]
mat_lambda_i_b = ref.mat_lambda_i_b[np.ix_(indexes, indexes)]
mat_omega_i_a = ref.mat_omega_i_a[np.ix_(indexes, indexes)]
mat_omega_i_b = ref.mat_omega_i_b[np.ix_(indexes, indexes)]
mat_kab = ref.mat_kab[np.ix_(indexes, indexes)]
mat_kba = ref.mat_kba[np.ix_(indexes, indexes)]

p_180x_i = np.diag([1.0, -1.0,
                    1.0, -1.0])

p_180y_i = np.diag([-1.0, 1.0,
                    -1.0, 1.0])


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        lambda_i_a=0.0, omega_i_a=0.0,
        lambda_i_b=0.0, omega_i_b=0.0):
    """TODO: make docstring
    """

    pa = 1.0 - pb
    kab, kba = kex_ab * np.asarray([pb, pa])

    liouvillian = (
        mat_lambda_i_a * lambda_i_a +
        mat_lambda_i_b * lambda_i_b +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
