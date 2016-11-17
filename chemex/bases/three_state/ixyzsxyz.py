"""TODO: module docstring.

# Operator basis:
# {Ix, Iy, Iz, Sx, Sy, Sz,
#  2IxSz, 2IySz, 2IzSx, 2IzSy,
#  2IxSx, 2IxSy, 2IySx, 2IySy,
#  2IzSz}{A,B,C}
"""

import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.three_state import exchange_model

mat_r2_i_a = ref.mat_r2_i_a[0:45, 0:45]
mat_r1_i_a = ref.mat_r1_i_a[0:45, 0:45]
mat_r2_s_a = ref.mat_r2_s_a[0:45, 0:45]
mat_r1_s_a = ref.mat_r1_s_a[0:45, 0:45]
mat_r2a_i_a = ref.mat_r2a_i_a[0:45, 0:45]
mat_r2a_s_a = ref.mat_r2a_s_a[0:45, 0:45]
mat_r2_mq_a = ref.mat_r2_mq_a[0:45, 0:45]
mat_r1a_a = ref.mat_r1a_a[0:45, 0:45]
mat_omega_i_a = ref.mat_omega_i_a[0:45, 0:45]
mat_omega_s_a = ref.mat_omega_s_a[0:45, 0:45]
mat_j_a = ref.mat_j_a[0:45, 0:45]
mat_etaxy_i_a = ref.mat_etaxy_i_a[0:45, 0:45]
mat_etaxy_s_a = ref.mat_etaxy_s_a[0:45, 0:45]
mat_etaz_i_a = ref.mat_etaz_i_a[0:45, 0:45]
mat_etaz_s_a = ref.mat_etaz_s_a[0:45, 0:45]
mat_sigma_a = ref.mat_sigma_a[0:45, 0:45]
mat_mu_mq_a = ref.mat_mu_mq_a[0:45, 0:45]

mat_r2_i_b = ref.mat_r2_i_b[0:45, 0:45]
mat_r1_i_b = ref.mat_r1_i_b[0:45, 0:45]
mat_r2_s_b = ref.mat_r2_s_b[0:45, 0:45]
mat_r1_s_b = ref.mat_r1_s_b[0:45, 0:45]
mat_r2a_i_b = ref.mat_r2a_i_b[0:45, 0:45]
mat_r2a_s_b = ref.mat_r2a_s_b[0:45, 0:45]
mat_r2_mq_b = ref.mat_r2_mq_b[0:45, 0:45]
mat_r1a_b = ref.mat_r1a_b[0:45, 0:45]
mat_omega_i_b = ref.mat_omega_i_b[0:45, 0:45]
mat_omega_s_b = ref.mat_omega_s_b[0:45, 0:45]
mat_j_b = ref.mat_j_b[0:45, 0:45]
mat_etaxy_i_b = ref.mat_etaxy_i_b[0:45, 0:45]
mat_etaxy_s_b = ref.mat_etaxy_s_b[0:45, 0:45]
mat_etaz_i_b = ref.mat_etaz_i_b[0:45, 0:45]
mat_etaz_s_b = ref.mat_etaz_s_b[0:45, 0:45]
mat_sigma_b = ref.mat_sigma_b[0:45, 0:45]
mat_mu_mq_b = ref.mat_mu_mq_b[0:45, 0:45]

mat_r2_i_c = ref.mat_r2_i_c[0:45, 0:45]
mat_r1_i_c = ref.mat_r1_i_c[0:45, 0:45]
mat_r2_s_c = ref.mat_r2_s_c[0:45, 0:45]
mat_r1_s_c = ref.mat_r1_s_c[0:45, 0:45]
mat_r2a_i_c = ref.mat_r2a_i_c[0:45, 0:45]
mat_r2a_s_c = ref.mat_r2a_s_c[0:45, 0:45]
mat_r2_mq_c = ref.mat_r2_mq_c[0:45, 0:45]
mat_r1a_c = ref.mat_r1a_c[0:45, 0:45]
mat_omega_i_c = ref.mat_omega_i_c[0:45, 0:45]
mat_omega_s_c = ref.mat_omega_s_c[0:45, 0:45]
mat_j_c = ref.mat_j_c[0:45, 0:45]
mat_etaxy_i_c = ref.mat_etaxy_i_c[0:45, 0:45]
mat_etaxy_s_c = ref.mat_etaxy_s_c[0:45, 0:45]
mat_etaz_i_c = ref.mat_etaz_i_c[0:45, 0:45]
mat_etaz_s_c = ref.mat_etaz_s_c[0:45, 0:45]
mat_sigma_c = ref.mat_sigma_c[0:45, 0:45]
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

index_iz = [2, 17, 32]
index_iz_a = [2]


# yapf: disable
def compute_liouvillian(
        pb=0.0, pc=0.0, kex_ab=0.0, kex_ac=0.0, kex_bc=0.0,
        r2_i_a=0.0, r1_i_a=0.0, r2a_i_a=0.0, etaxy_i_a=0.0, etaz_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, r2a_i_b=0.0, etaxy_i_b=0.0, etaz_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, r1_i_c=0.0, r2a_i_c=0.0, etaxy_i_c=0.0, etaz_i_c=0.0, omega_i_c=0.0,
        r2_s_a=0.0, r1_s_a=0.0, r2a_s_a=0.0, etaxy_s_a=0.0, etaz_s_a=0.0, omega_s_a=0.0,
        r2_s_b=0.0, r1_s_b=0.0, r2a_s_b=0.0, etaxy_s_b=0.0, etaz_s_b=0.0, omega_s_b=0.0,
        r2_s_c=0.0, r1_s_c=0.0, r2a_s_c=0.0, etaxy_s_c=0.0, etaz_s_c=0.0, omega_s_c=0.0,
        r2_mq_a=0.0, r1a_a=0.0, mu_mq_a=0.0, sigma_a=0.0, j_a=0.0,
        r2_mq_b=0.0, r1a_b=0.0, mu_mq_b=0.0, sigma_b=0.0, j_b=0.0,
        r2_mq_c=0.0, r1a_c=0.0, mu_mq_c=0.0, sigma_c=0.0, j_c=0.0,
        omega1x_i=0.0, omega1y_i=0.0, omega1x_s=0.0, omega1y_s=0.0,
        **kwargs):
    """Compute the Liouvillian."""
    pa = 1.0 - pb - pc

    kab = kba = kbc = kcb = kac = kca = 0.0

    if pa + pb:
        kab, kba = kex_ab / (pa + pb) * np.array([pb, pa])
    if pb + pc:
        kbc, kcb = kex_bc / (pb + pc) * np.array([pc, pb])
    if pa + pc:
        kac, kca = kex_ac / (pa + pc) * np.array([pc, pa])

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +
        mat_r2_i_c * r2_i_c +

        mat_r2_s_a * r2_s_a +
        mat_r2_s_b * r2_s_b +
        mat_r2_s_c * r2_s_c +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +
        mat_r1_i_c * r1_i_c +

        mat_r1_s_a * r1_s_a +
        mat_r1_s_b * r1_s_b +
        mat_r1_s_c * r1_s_c +

        mat_r2a_i_a * r2a_i_a +
        mat_r2a_i_b * r2a_i_b +
        mat_r2a_i_c * r2a_i_c +

        mat_r2a_s_a * r2a_s_a +
        mat_r2a_s_b * r2a_s_b +
        mat_r2a_s_c * r2a_s_c +

        mat_etaxy_i_a * etaxy_i_a +
        mat_etaxy_i_b * etaxy_i_b +
        mat_etaxy_i_c * etaxy_i_c +

        mat_etaxy_s_a * etaxy_s_a +
        mat_etaxy_s_b * etaxy_s_b +
        mat_etaxy_s_c * etaxy_s_c +

        mat_etaz_i_a * etaz_i_a +
        mat_etaz_i_b * etaz_i_b +
        mat_etaz_i_c * etaz_i_c +

        mat_etaz_s_a * etaz_s_a +
        mat_etaz_s_b * etaz_s_b +
        mat_etaz_s_c * etaz_s_c +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +
        mat_omega_s_c * omega_s_c +

        mat_r1a_a * r1a_a +
        mat_r1a_b * r1a_b +
        mat_r1a_c * r1a_c +

        mat_r2_mq_a * r2_mq_a +
        mat_r2_mq_b * r2_mq_b +
        mat_r2_mq_c * r2_mq_c +

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
# yapf: enable


def compute_equilibrium(pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((45, 1))
    mag0[index_iz] = [[1.0 - pb - pc], [pb], [pc]]

    return mag0


def create_default_params(model=None,
                          temperature=None,
                          nuclei=None,
                          h_larmor_frq=None,
                          p_total=None,
                          l_total=None):
    """Create the default experimental and fitting parameters."""
    resonance_i, resonance_s = nuclei.resonances
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs2 = {'temperature': temperature, 'nuclei': resonance_i['name']}
    kwargs3 = {
        'temperature': temperature,
        'nuclei': resonance_i['name'],
        'h_larmor_frq': h_larmor_frq
    }
    kwargs4 = {'temperature': temperature, 'nuclei': resonance_s['name']}
    kwargs5 = {
        'temperature': temperature,
        'nuclei': resonance_s['name'],
        'h_larmor_frq': h_larmor_frq
    }
    kwargs6 = {'temperature': temperature, 'nuclei': nuclei.assignment}
    kwargs7 = {
        'temperature': temperature,
        'nuclei': nuclei.assignment,
        'h_larmor_frq': h_larmor_frq
    }

    map_names = {
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'pc': parameters.ParameterName('pc', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'kex_bc': parameters.ParameterName('kex_bc', **kwargs1).to_full_name(),
        'kex_ac': parameters.ParameterName('kex_ac', **kwargs1).to_full_name(),
        'dw_i_ab': parameters.ParameterName('dw_ab', **kwargs2).to_full_name(),
        'dw_i_ac': parameters.ParameterName('dw_ac', **kwargs2).to_full_name(),
        'dw_s_ab': parameters.ParameterName('dw_ab', **kwargs4).to_full_name(),
        'dw_s_ac': parameters.ParameterName('dw_ac', **kwargs4).to_full_name(),
        'cs_i_a': parameters.ParameterName('cs_a', **kwargs2).to_full_name(),
        'r2_i_a': parameters.ParameterName('r2_a', **kwargs3).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'r2a_i_a': parameters.ParameterName('r2a_a', **kwargs3).to_full_name(),
        'etaxy_i_a': parameters.ParameterName('etaxy_a', **kwargs3).to_full_name(),
        'etaz_i_a': parameters.ParameterName('etaz_a', **kwargs3).to_full_name(),
        'cs_s_a': parameters.ParameterName('cs_a', **kwargs4).to_full_name(),
        'r2_s_a': parameters.ParameterName('r2_a', **kwargs5).to_full_name(),
        'r1_s_a': parameters.ParameterName('r1_a', **kwargs5).to_full_name(),
        'r2a_s_a': parameters.ParameterName('r2a_a', **kwargs5).to_full_name(),
        'etaxy_s_a': parameters.ParameterName('etaxy_a', **kwargs5).to_full_name(),
        'etaz_s_a': parameters.ParameterName('etaz_a', **kwargs5).to_full_name(),
        'r1a_a': parameters.ParameterName('r1a_a', **kwargs7).to_full_name(),
        'r2_mq_a': parameters.ParameterName('r2_mq_a', **kwargs7).to_full_name(),
        'mu_mq_a': parameters.ParameterName('mu_mq_a', **kwargs7).to_full_name(),
        'sigma_a': parameters.ParameterName('sigma_a', **kwargs7).to_full_name(),
        'j_a': parameters.ParameterName('j_a', **kwargs6).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kwargs2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
        'r2a_i_b': parameters.ParameterName('r2a_b', **kwargs3).to_full_name(),
        'etaxy_i_b': parameters.ParameterName('etaxy_b', **kwargs3).to_full_name(),
        'etaz_i_b': parameters.ParameterName('etaz_b', **kwargs3).to_full_name(),
        'cs_s_b': parameters.ParameterName('cs_b', **kwargs4).to_full_name(),
        'r2_s_b': parameters.ParameterName('r2_b', **kwargs5).to_full_name(),
        'r1_s_b': parameters.ParameterName('r1_b', **kwargs5).to_full_name(),
        'r2a_s_b': parameters.ParameterName('r2a_b', **kwargs5).to_full_name(),
        'etaxy_s_b': parameters.ParameterName('etaxy_b', **kwargs5).to_full_name(),
        'etaz_s_b': parameters.ParameterName('etaz_b', **kwargs5).to_full_name(),
        'r1a_b': parameters.ParameterName('r1a_b', **kwargs7).to_full_name(),
        'r2_mq_b': parameters.ParameterName('r2_mq_b', **kwargs7).to_full_name(),
        'mu_mq_b': parameters.ParameterName('mu_mq_b', **kwargs7).to_full_name(),
        'sigma_b': parameters.ParameterName('sigma_b', **kwargs7).to_full_name(),
        'j_b': parameters.ParameterName('j_b', **kwargs6).to_full_name(),
        'cs_i_c': parameters.ParameterName('cs_c', **kwargs2).to_full_name(),
        'r2_i_c': parameters.ParameterName('r2_c', **kwargs3).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kwargs3).to_full_name(),
        'r2a_i_c': parameters.ParameterName('r2a_c', **kwargs3).to_full_name(),
        'etaxy_i_c': parameters.ParameterName('etaxy_c', **kwargs3).to_full_name(),
        'etaz_i_c': parameters.ParameterName('etaz_c', **kwargs3).to_full_name(),
        'cs_s_c': parameters.ParameterName('cs_c', **kwargs4).to_full_name(),
        'r2_s_c': parameters.ParameterName('r2_c', **kwargs5).to_full_name(),
        'r1_s_c': parameters.ParameterName('r1_c', **kwargs5).to_full_name(),
        'r2a_s_c': parameters.ParameterName('r2a_c', **kwargs5).to_full_name(),
        'etaxy_s_c': parameters.ParameterName('etaxy_c', **kwargs5).to_full_name(),
        'etaz_s_c': parameters.ParameterName('etaz_c', **kwargs5).to_full_name(),
        'r1a_c': parameters.ParameterName('r1a_c', **kwargs7).to_full_name(),
        'r2_mq_c': parameters.ParameterName('r2_mq_c', **kwargs7).to_full_name(),
        'mu_mq_c': parameters.ParameterName('mu_mq_c', **kwargs7).to_full_name(),
        'sigma_c': parameters.ParameterName('sigma_c', **kwargs7).to_full_name(),
        'j_c': parameters.ParameterName('j_c', **kwargs6).to_full_name(),
    }

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    cs_i_c = '{cs_i_a} + {dw_i_ac}'.format(**map_names)
    cs_s_b = '{cs_s_a} + {dw_s_ab}'.format(**map_names)
    cs_s_c = '{cs_s_a} + {dw_s_ac}'.format(**map_names)
    r1_i_c = r1_i_b = map_names['r1_i_a']
    r2_i_c = r2_i_b = map_names['r2_i_a']
    r1_s_c = r1_s_b = map_names['r1_s_a']
    r2_s_c = r2_s_b = map_names['r2_s_a']
    r1a_c = r1a_b = map_names['r1a_a']
    r2a_i_c = r2a_i_b = map_names['r2a_i_a']
    r2a_s_c = r2a_s_b = map_names['r2a_s_a']
    etaxy_i_c = etaxy_i_b = map_names['etaxy_i_a']
    etaz_i_c = etaz_i_b = map_names['etaz_i_a']
    etaxy_s_c = etaxy_s_b = map_names['etaxy_s_a']
    etaz_s_c = etaz_s_b = map_names['etaz_s_a']
    sigma_c = sigma_b = map_names['sigma_a']
    mu_mq_c = mu_mq_b = map_names['mu_mq_a']
    r2_mq_c = r2_mq_b = map_names['r2_mq_a']
    j_c = j_b = map_names['j_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.025, True, 0.0, 1.0, None),
        (map_names['pc'], 0.025, True, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['kex_bc'], 200.0, True, 0.0, None, None),
        (map_names['kex_ac'], 0.0, False, 0.0, None, None),
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_i_ac'], 0.0, True, None, None, None),
        (map_names['dw_s_ab'], 0.0, False, None, None, None),
        (map_names['dw_s_ac'], 0.0, False, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['r2a_i_a'], 15.0, False, 0.0, None, None),
        (map_names['etaxy_i_a'], 0.0, True, None, None, None),
        (map_names['etaz_i_a'], 0.0, True, None, None, None),
        (map_names['cs_s_a'], 0.0, False, None, None, None),
        (map_names['r2_s_a'], 10.0, False, 0.0, None, None),
        (map_names['r1_s_a'], 1.0, False, 0.0, None, None),
        (map_names['r2a_s_a'], 15.0, False, 0.0, None, None),
        (map_names['etaxy_s_a'], 0.0, False, None, None, None),
        (map_names['etaz_s_a'], 0.0, False, None, None, None),
        (map_names['mu_mq_a'], 0.0, False, None, None, None),
        (map_names['sigma_a'], 0.0, False, None, None, None),
        (map_names['r2_mq_a'], 10.0, False, 0.0, None, None),
        (map_names['r1a_a'], 2.0, False, 0.0, None, None),
        (map_names['j_a'], -93.0, False, None, None, None),
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
        (map_names['r1_i_b'], 1.0, None, 0.0, None, r1_i_b),
        (map_names['r2a_i_b'], 15.0, None, 0.0, None, r2a_i_b),
        (map_names['etaxy_i_b'], 0.0, None, None, None, etaxy_i_b),
        (map_names['etaz_i_b'], 0.0, None, None, None, etaz_i_b),
        (map_names['cs_s_b'], 0.0, None, None, None, cs_s_b),
        (map_names['r2_s_b'], 10.0, None, 0.0, None, r2_s_b),
        (map_names['r1_s_b'], 1.0, None, 0.0, None, r1_s_b),
        (map_names['r2a_s_b'], 15.0, None, 0.0, None, r2a_s_b),
        (map_names['etaxy_s_b'], 0.0, None, None, None, etaxy_s_b),
        (map_names['etaz_s_b'], 0.0, None, None, None, etaz_s_b),
        (map_names['mu_mq_b'], 0.0, None, None, None, mu_mq_b),
        (map_names['sigma_b'], 0.0, None, None, None, sigma_b),
        (map_names['r2_mq_b'], 10.0, None, 0.0, None, r2_mq_b),
        (map_names['r1a_b'], 2.0, None, 0.0, None, r1a_b),
        (map_names['j_b'], -93.0, None, None, None, j_b),
        (map_names['cs_i_c'], 0.0, None, None, None, cs_i_c),
        (map_names['r2_i_c'], 10.0, None, 0.0, None, r2_i_c),
        (map_names['r1_i_c'], 1.0, None, 0.0, None, r1_i_c),
        (map_names['r2a_i_c'], 15.0, None, 0.0, None, r2a_i_c),
        (map_names['etaxy_i_c'], 0.0, None, None, None, etaxy_i_c),
        (map_names['etaz_i_c'], 0.0, None, None, None, etaz_i_c),
        (map_names['cs_s_c'], 0.0, None, None, None, cs_s_c),
        (map_names['r2_s_c'], 10.0, None, 0.0, None, r2_s_c),
        (map_names['r1_s_c'], 1.0, None, 0.0, None, r1_s_c),
        (map_names['r2a_s_c'], 15.0, None, 0.0, None, r2a_s_c),
        (map_names['etaxy_s_c'], 0.0, None, None, None, etaxy_s_c),
        (map_names['etaz_s_c'], 0.0, None, None, None, etaz_s_c),
        (map_names['mu_mq_c'], 0.0, None, None, None, mu_mq_c),
        (map_names['sigma_c'], 0.0, None, None, None, sigma_c),
        (map_names['r2_mq_c'], 10.0, None, 0.0, None, r2_mq_c),
        (map_names['r1a_c'], 2.0, None, 0.0, None, r1a_c),
        (map_names['j_c'], -93.0, None, None, None, j_c), )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total,
                                                     l_total)

    return map_names, params
