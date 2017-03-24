"""TODO: module docstring."""

import lmfit
import numpy as np
from scipy import linalg

from chemex import parameters
from chemex.bases import ref_eq as ref
from chemex.bases.four_state import exchange_model

# yapf: disable
indexes = [0, 1, 2,
           15, 16, 17,
           30, 31, 32,
           45, 46, 47,
           60]
mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]
mat_eq_i_a = ref.mat_eq_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]
mat_eq_i_b = ref.mat_eq_i_b[mesh]

mat_r2_i_c = ref.mat_r2_i_c[mesh]
mat_r1_i_c = ref.mat_r1_i_c[mesh]
mat_omega_i_c = ref.mat_omega_i_c[mesh]
mat_eq_i_c = ref.mat_eq_i_c[mesh]

mat_r2_i_d = ref.mat_r2_i_d[mesh]
mat_r1_i_d = ref.mat_r1_i_d[mesh]
mat_omega_i_d = ref.mat_omega_i_d[mesh]
mat_eq_i_d = ref.mat_eq_i_d[mesh]

mat_omega1x_i = ref.mat_omega1x_i[mesh]
mat_omega1y_i = ref.mat_omega1y_i[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]
mat_kac = ref.mat_kac[mesh]
mat_kca = ref.mat_kca[mesh]
mat_kad = ref.mat_kad[mesh]
mat_kda = ref.mat_kda[mesh]
mat_kbc = ref.mat_kbc[mesh]
mat_kcb = ref.mat_kcb[mesh]
mat_kbd = ref.mat_kbd[mesh]
mat_kdb = ref.mat_kdb[mesh]
mat_kcd = ref.mat_kcd[mesh]
mat_kdc = ref.mat_kdc[mesh]

index_iz = [2, 5, 8, 11]
index_iz_eq = [2, 5, 8, 11, 12]
index_iz_a = [2]
index_iz_b = [5]
index_iz_c = [8]
index_iz_d = [11]


def compute_liouvillian(
        pb=0.0, pc=0.0, pd=0.0, kex_ab=0.0, kex_ac=0.0, kex_ad=0.0, kex_bc=0.0, kex_bd=0.0, kex_cd=0.0,
        r2_i_a=0.0, r1_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, r1_i_c=0.0, omega_i_c=0.0,
        r2_i_d=0.0, r1_i_d=0.0, omega_i_d=0.0,
        omega1x_i=0.0, omega1y_i=0.0,
        **kwargs):
    """Compute the Liouvillian."""
    pa = 1.0 - pb - pc - pd

    kab = kba = kac = kca = kad = kda = kbc = kcb = kbd = kdb = kcd = kdc = 0.0

    if pa + pb:
        kab, kba = kex_ab / (pa + pb) * np.asarray([pb, pa])
    if pa + pc:
        kac, kca = kex_ac / (pa + pc) * np.asarray([pc, pa])
    if pa + pd:
        kad, kda = kex_ad / (pa + pd) * np.asarray([pd, pa])
    if pb + pc:
        kbc, kcb = kex_bc / (pb + pc) * np.asarray([pc, pb])
    if pb + pd:
        kbd, kdb = kex_bd / (pb + pd) * np.asarray([pd, pb])
    if pc + pd:
        kcd, kdc = kex_cd / (pc + pd) * np.asarray([pd, pc])

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +
        mat_r2_i_c * r2_i_c +
        mat_r2_i_d * r2_i_d +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +
        mat_r1_i_c * r1_i_c +
        mat_r1_i_d * r1_i_d +

        mat_eq_i_a * r1_i_a * pa +
        mat_eq_i_b * r1_i_b * pb +
        mat_eq_i_c * r1_i_c * pc +
        mat_eq_i_d * r1_i_d * pd +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +
        mat_omega_i_d * omega_i_d +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +

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


# yapf: enable


def compute_equilibrium_after_d1(
        time_d1=0.0,
        pb=0.0, pc=0.0, pd=0.0,
        kex_ab=0.0, kex_ac=0.0, kex_ad=0.0, kex_bc=0.0, kex_bd=0.0, kex_cd=0.0,
        r1_i_a=0.0, r1_i_b=0.0, r1_i_c=0.0, r1_i_d=0.0, **kwargs):
    """Compute the equilibrium magnetization."""

    mag0 = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]).T

    liouvillian_free = compute_liouvillian(
        pb=pb, pc=pc, pd=pd,
        kex_ab=kex_ab, kex_ac=kex_ac, kex_ad=kex_ad,
        kex_bc=kex_bc, kex_bd=kex_bd, kex_cd=kex_cd,
        r1_i_a=r1_i_a, r1_i_b=r1_i_b, r1_i_c=r1_i_c, r1_i_d=r1_i_d)

    return linalg.expm(liouvillian_free * time_d1).dot(mag0)


def create_default_params(model=None,
                          temperature=None,
                          nuclei=None,
                          h_larmor_frq=None,
                          p_total=None,
                          l_total=None):
    """Create the default experimental and fitting parameters."""
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs2 = {'temperature': temperature, 'nuclei': nuclei}
    kwargs3 = {'temperature': temperature, 'nuclei': nuclei, 'h_larmor_frq': h_larmor_frq}

    map_names = {
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'pc': parameters.ParameterName('pc', **kwargs1).to_full_name(),
        'pd': parameters.ParameterName('pd', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'kex_ac': parameters.ParameterName('kex_ac', **kwargs1).to_full_name(),
        'kex_ad': parameters.ParameterName('kex_ad', **kwargs1).to_full_name(),
        'kex_bc': parameters.ParameterName('kex_bc', **kwargs1).to_full_name(),
        'kex_bd': parameters.ParameterName('kex_bd', **kwargs1).to_full_name(),
        'kex_cd': parameters.ParameterName('kex_cd', **kwargs1).to_full_name(),
        'cs_i_a': parameters.ParameterName('cs_a', **kwargs2).to_full_name(),
        'r2_i_a': parameters.ParameterName('r2_a', **kwargs3).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kwargs2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
        'cs_i_c': parameters.ParameterName('cs_c', **kwargs2).to_full_name(),
        'r2_i_c': parameters.ParameterName('r2_c', **kwargs3).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kwargs3).to_full_name(),
        'cs_i_d': parameters.ParameterName('cs_d', **kwargs2).to_full_name(),
        'r2_i_d': parameters.ParameterName('r2_d', **kwargs3).to_full_name(),
        'r1_i_d': parameters.ParameterName('r1_d', **kwargs3).to_full_name(),
    }

    r1_i_b = r1_i_c = r1_i_d = map_names['r1_i_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.025, True, 0.0, 1.0, None),
        (map_names['pc'], 0.0, False, 0.0, 1.0, None),
        (map_names['pd'], 0.0, False, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['kex_ac'], 0.0, False, 0.0, None, None),
        (map_names['kex_ad'], 0.0, False, 0.0, None, None),
        (map_names['kex_bc'], 0.0, False, 0.0, None, None),
        (map_names['kex_bd'], 0.0, False, 0.0, None, None),
        (map_names['kex_cd'], 0.0, False, 0.0, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['cs_i_b'], 0.0, True, None, None, None),
        (map_names['r2_i_b'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_b'], 1.0, True, 0.0, None, r1_i_b),
        (map_names['cs_i_c'], 0.0, True, None, None, None),
        (map_names['r2_i_c'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_c'], 1.0, True, 0.0, None, r1_i_c),
        (map_names['cs_i_d'], 0.0, True, None, None, None),
        (map_names['r2_i_d'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_d'], 1.0, True, 0.0, None, r1_i_d), )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total, l_total)

    return map_names, params
