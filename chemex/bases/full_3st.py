from __future__ import absolute_import
# Operator basis:
# {Ix, Iy, Iz, Sx, Sy, Sz,
#  2IxSz, 2IySz, 2IzSx, 2IzSy,
#  2IxSx, 2IxSy, 2IySx, 2IySy,
#  2IzSz}{A,B,C}

import numpy as np

from chemex.bases import ref

pi = np.pi

mat_lambda_i_a = ref.mat_lambda_i_a[0:45, 0:45]
mat_lambda_i_b = ref.mat_lambda_i_b[0:45, 0:45]
mat_lambda_i_c = ref.mat_lambda_i_c[0:45, 0:45]
mat_rho_i_a = ref.mat_rho_i_a[0:45, 0:45]
mat_rho_i_b = ref.mat_rho_i_b[0:45, 0:45]
mat_rho_i_c = ref.mat_rho_i_c[0:45, 0:45]
mat_lambda_s_a = ref.mat_lambda_s_a[0:45, 0:45]
mat_lambda_s_b = ref.mat_lambda_s_b[0:45, 0:45]
mat_lambda_s_c = ref.mat_lambda_s_c[0:45, 0:45]
mat_rho_s_a = ref.mat_rho_s_a[0:45, 0:45]
mat_rho_s_b = ref.mat_rho_s_b[0:45, 0:45]
mat_rho_s_c = ref.mat_rho_s_c[0:45, 0:45]
mat_rhoa_i_a = ref.mat_rhoa_i_a[0:45, 0:45]
mat_rhoa_i_b = ref.mat_rhoa_i_b[0:45, 0:45]
mat_rhoa_i_c = ref.mat_rhoa_i_c[0:45, 0:45]
mat_rhoa_s_a = ref.mat_rhoa_s_a[0:45, 0:45]
mat_rhoa_s_b = ref.mat_rhoa_s_b[0:45, 0:45]
mat_rhoa_s_c = ref.mat_rhoa_s_c[0:45, 0:45]
mat_lambda_mq_a = ref.mat_lambda_mq_a[0:45, 0:45]
mat_lambda_mq_b = ref.mat_lambda_mq_b[0:45, 0:45]
mat_lambda_mq_c = ref.mat_lambda_mq_c[0:45, 0:45]
mat_rho_is_a = ref.mat_rho_is_a[0:45, 0:45]
mat_rho_is_b = ref.mat_rho_is_b[0:45, 0:45]
mat_rho_is_c = ref.mat_rho_is_c[0:45, 0:45]
mat_omega_i_a = ref.mat_omega_i_a[0:45, 0:45]
mat_omega_i_b = ref.mat_omega_i_b[0:45, 0:45]
mat_omega_i_c = ref.mat_omega_i_c[0:45, 0:45]
mat_omega_s_a = ref.mat_omega_s_a[0:45, 0:45]
mat_omega_s_b = ref.mat_omega_s_b[0:45, 0:45]
mat_omega_s_c = ref.mat_omega_s_c[0:45, 0:45]
mat_j_a = ref.mat_j_a[0:45, 0:45]
mat_j_b = ref.mat_j_b[0:45, 0:45]
mat_j_c = ref.mat_j_c[0:45, 0:45]
mat_eta_i_a = ref.mat_eta_i_a[0:45, 0:45]
mat_eta_i_b = ref.mat_eta_i_b[0:45, 0:45]
mat_eta_i_c = ref.mat_eta_i_c[0:45, 0:45]
mat_eta_s_a = ref.mat_eta_s_a[0:45, 0:45]
mat_eta_s_b = ref.mat_eta_s_b[0:45, 0:45]
mat_eta_s_c = ref.mat_eta_s_c[0:45, 0:45]
mat_delta_i_a = ref.mat_delta_i_a[0:45, 0:45]
mat_delta_i_b = ref.mat_delta_i_b[0:45, 0:45]
mat_delta_i_c = ref.mat_delta_i_c[0:45, 0:45]
mat_delta_s_a = ref.mat_delta_s_a[0:45, 0:45]
mat_delta_s_b = ref.mat_delta_s_b[0:45, 0:45]
mat_delta_s_c = ref.mat_delta_s_c[0:45, 0:45]
mat_sigma_a = ref.mat_sigma_a[0:45, 0:45]
mat_sigma_b = ref.mat_sigma_b[0:45, 0:45]
mat_sigma_c = ref.mat_sigma_c[0:45, 0:45]
mat_mu_mq_a = ref.mat_mu_mq_a[0:45, 0:45]
mat_mu_mq_b = ref.mat_mu_mq_b[0:45, 0:45]
mat_mu_mq_c = ref.mat_mu_mq_c[0:45, 0:45]
mat_omega1x_i = ref.mat_omega1x_i[0:45, 0:45]
mat_omega1y_i = ref.mat_omega1y_i[0:45, 0:45]
mat_omega1x_s = ref.mat_omega1x_s[0:45, 0:45]
mat_omega1y_s = ref.mat_omega1y_s[0:45, 0:45]
mat_kab = ref.mat_kab[0:45, 0:45]
mat_kba = ref.mat_kba[0:45, 0:45]
mat_kac = ref.mat_kac[0:45, 0:45]
mat_kca = ref.mat_kca[0:45, 0:45]
mat_kbc = ref.mat_kbc[0:45, 0:45]
mat_kcb = ref.mat_kcb[0:45, 0:45]


def compute_liouvillian(
        pb, pc, kex_ab, kex_ac, kex_bc,
        lambda_i_a=0.0, rho_i_a=0.0, rhoa_i_a=0.0, eta_i_a=0.0, delta_i_a=0.0, omega_i_a=0.0,
        lambda_i_b=0.0, rho_i_b=0.0, rhoa_i_b=0.0, eta_i_b=0.0, delta_i_b=0.0, omega_i_b=0.0,
        lambda_i_c=0.0, rho_i_c=0.0, rhoa_i_c=0.0, eta_i_c=0.0, delta_i_c=0.0, omega_i_c=0.0,
        lambda_s_a=0.0, rho_s_a=0.0, rhoa_s_a=0.0, eta_s_a=0.0, delta_s_a=0.0, omega_s_a=0.0,
        lambda_s_b=0.0, rho_s_b=0.0, rhoa_s_b=0.0, eta_s_b=0.0, delta_s_b=0.0, omega_s_b=0.0,
        lambda_s_c=0.0, rho_s_c=0.0, rhoa_s_c=0.0, eta_s_c=0.0, delta_s_c=0.0, omega_s_c=0.0,
        lambda_mq_a=0.0, rho_is_a=0.0, mu_mq_a=0.0, sigma_a=0.0, j_a=0.0,
        lambda_mq_b=0.0, rho_is_b=0.0, mu_mq_b=0.0, sigma_b=0.0, j_b=0.0,
        lambda_mq_c=0.0, rho_is_c=0.0, mu_mq_c=0.0, sigma_c=0.0, j_c=0.0,
        omega1x_i=0.0, omega1y_i=0.0, omega1x_s=0.0, omega1y_s=0.0):
    """TODO: make docstring
    """

    pa = 1.0 - pb - pc
    kab, kba = kex_ab * np.asarray([pb, pa])
    kac, kca = kex_ac * np.asarray([pc, pa])
    kbc, kcb = kex_bc * np.asarray([pc, pb])

    liouvillian = (
        mat_lambda_i_a * lambda_i_a +
        mat_lambda_i_b * lambda_i_b +
        mat_lambda_i_c * lambda_i_c +

        mat_lambda_s_a * lambda_s_a +
        mat_lambda_s_b * lambda_s_b +
        mat_lambda_s_c * lambda_s_c +

        mat_rho_i_a * rho_i_a +
        mat_rho_i_b * rho_i_b +
        mat_rho_i_c * rho_i_c +

        mat_rho_s_a * rho_s_a +
        mat_rho_s_b * rho_s_b +
        mat_rho_s_c * rho_s_c +

        mat_rhoa_i_a * rhoa_i_a +
        mat_rhoa_i_b * rhoa_i_b +
        mat_rhoa_i_c * rhoa_i_c +

        mat_rhoa_s_a * rhoa_s_a +
        mat_rhoa_s_b * rhoa_s_b +
        mat_rhoa_s_c * rhoa_s_c +

        mat_eta_i_a * eta_i_a +
        mat_eta_i_b * eta_i_b +
        mat_eta_i_c * eta_i_c +

        mat_eta_s_a * eta_s_a +
        mat_eta_s_b * eta_s_b +
        mat_eta_s_c * eta_s_c +

        mat_delta_i_a * delta_i_a +
        mat_delta_i_b * delta_i_b +
        mat_delta_i_c * delta_i_c +

        mat_delta_s_a * delta_s_a +
        mat_delta_s_b * delta_s_b +
        mat_delta_s_c * delta_s_c +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +
        mat_omega_s_c * omega_s_c +

        mat_rho_is_a * rho_is_a +
        mat_rho_is_b * rho_is_b +
        mat_rho_is_c * rho_is_c +

        mat_lambda_mq_a * lambda_mq_a +
        mat_lambda_mq_b * lambda_mq_b +
        mat_lambda_mq_c * lambda_mq_c +

        mat_mu_mq_a * mu_mq_a +
        mat_mu_mq_b * mu_mq_b +
        mat_mu_mq_c * mu_mq_c +

        mat_sigma_a * sigma_a +
        mat_sigma_b * sigma_b +
        mat_sigma_c * sigma_c +

        mat_j_a * j_a +
        mat_j_b * j_b +
        mat_j_c * j_c +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +
        mat_omega1x_s * omega1x_s +
        mat_omega1y_s * omega1y_s +

        mat_kab * kab +
        mat_kba * kba +
        mat_kac * kac +
        mat_kca * kca +
        mat_kbc * kbc +
        mat_kcb * kcb
    )

    return liouvillian
