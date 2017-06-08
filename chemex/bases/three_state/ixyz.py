"""TODO: module docstring."""

import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.three_state import exchange_model

# yapf: disable
indexes = [0, 1, 2,
           15, 16, 17,
           30, 31, 32]
mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]

mat_r2_i_c = ref.mat_r2_i_c[mesh]
mat_r1_i_c = ref.mat_r1_i_c[mesh]
mat_omega_i_c = ref.mat_omega_i_c[mesh]

mat_omega1x_i = ref.mat_omega1x_i[mesh]
mat_omega1y_i = ref.mat_omega1y_i[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]
mat_kbc = ref.mat_kbc[mesh]
mat_kcb = ref.mat_kcb[mesh]
mat_kac = ref.mat_kac[mesh]
mat_kca = ref.mat_kca[mesh]

index_iz = [2, 5, 8]
index_iz_a = [2]


def compute_liouvillian(
        kab=0.0, kba=0.0, kac=0.0, kca=0.0, kbc=0.0, kcb=0.0,
        r2_i_a=0.0, r1_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, r1_i_c=0.0, omega_i_c=0.0,
        omega1x_i=0.0, omega1y_i=0.0,
        **kwargs):
    """Compute the Liouvillian."""

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +
        mat_r2_i_c * r2_i_c +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +
        mat_r1_i_c * r1_i_c +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +

        mat_kab * kab +
        mat_kba * kba +
        mat_kbc * kbc +
        mat_kcb * kcb +
        mat_kac * kac +
        mat_kca * kca
    )

    return liouvillian
# yapf: enable


def compute_equilibrium(pa=0.0, pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    return np.array([[0.0, 0.0, pa, 0.0, 0.0, pb, 0.0, 0.0, pc]]).T


def create_default_params(model=None,
                          temperature=None,
                          nuclei=None,
                          h_larmor_frq=None,
                          p_total=None,
                          l_total=None):
    """Create the default experimental and fitting parameters."""

    map_names, params = exchange_model.create_exchange_params(model, temperature, p_total, l_total)

    resonance_i, resonance_s = nuclei.resonances

    kw1 = {'temperature': temperature, 'nuclei': resonance_i['name']}

    map_names.update({
        'dw_i_ab': parameters.ParameterName('dw_ab', **kw1).to_full_name(),
        'dw_i_ac': parameters.ParameterName('dw_ac', **kw1).to_full_name(),
        'cs_i_a': parameters.ParameterName('cs_a', **kw1).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kw1).to_full_name(),
        'cs_i_c': parameters.ParameterName('cs_c', **kw1).to_full_name(),
    })

    kw2 = {'temperature': temperature, 'nuclei': resonance_i['name'], 'h_larmor_frq': h_larmor_frq}

    map_names.update({
        'r2_i_a': parameters.ParameterName('r2_a', **kw2).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kw2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kw2).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kw2).to_full_name(),
        'r2_i_c': parameters.ParameterName('r2_c', **kw2).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kw2).to_full_name(),
    })

    params.add_many(
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_i_ac'], 0.0, False, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None), )

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    cs_i_c = '{cs_i_a} + {dw_i_ac}'.format(**map_names)
    r1_i_c = r1_i_b = map_names['r1_i_a']
    r2_i_c = r2_i_b = map_names['r2_i_a']

    params.add_many(
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
        (map_names['r1_i_b'], 1.0, None, 0.0, None, r1_i_b),
        (map_names['cs_i_c'], 0.0, None, None, None, cs_i_c),
        (map_names['r2_i_c'], 10.0, None, 0.0, None, r2_i_c),
        (map_names['r1_i_c'], 1.0, None, 0.0, None, r1_i_c), )

    return map_names, params
