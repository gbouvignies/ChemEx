import numpy as np

from chemex.bases import ref

indexes = [0, 1, 15, 16, 30, 31]

mat_lambda_i_a = ref.mat_lambda_i_a[np.ix_(indexes, indexes)]
mat_lambda_i_b = ref.mat_lambda_i_b[np.ix_(indexes, indexes)]
mat_lambda_i_c = ref.mat_lambda_i_c[np.ix_(indexes, indexes)]
mat_omega_i_a = ref.mat_omega_i_a[np.ix_(indexes, indexes)]
mat_omega_i_b = ref.mat_omega_i_b[np.ix_(indexes, indexes)]
mat_omega_i_c = ref.mat_omega_i_c[np.ix_(indexes, indexes)]
mat_kab = ref.mat_kab[np.ix_(indexes, indexes)]
mat_kba = ref.mat_kba[np.ix_(indexes, indexes)]
mat_kbc = ref.mat_kbc[np.ix_(indexes, indexes)]
mat_kcb = ref.mat_kcb[np.ix_(indexes, indexes)]
mat_kac = ref.mat_kac[np.ix_(indexes, indexes)]
mat_kca = ref.mat_kca[np.ix_(indexes, indexes)]

p_180x_i = np.diag([1.0, -1.0,
                    1.0, -1.0,
                    1.0, -1.0])

p_180y_i = np.diag([-1.0, 1.0,
                    -1.0, 1.0,
                    -1.0, 1.0])


def compute_liouvillian(
        pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
        lambda_i_a=0.0, omega_i_a=0.0,
        lambda_i_b=0.0, omega_i_b=0.0,
        lambda_i_c=0.0, omega_i_c=0.0):
    """TODO: make docstring
    """

    pa = 1.0 - pb - pc

    kex_ab_ = kex_ab / (pa + pb)
    kex_bc_ = kex_bc / (pb + pc)
    kex_ac_ = kex_ac / (pa + pc)

    kab, kba = kex_ab_ * np.asarray([pb, pa])
    kbc, kcb = kex_bc_ * np.asarray([pc, pb])
    kac, kca = kex_ac_ * np.asarray([pc, pa])

    liouvillian = (
        mat_lambda_i_a * lambda_i_a +
        mat_lambda_i_b * lambda_i_b +
        mat_lambda_i_c * lambda_i_c +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +

        mat_kab * kab +
        mat_kba * kba +
        mat_kbc * kbc +
        mat_kcb * kcb +
        mat_kac * kac +
        mat_kca * kca
    )

    return liouvillian
