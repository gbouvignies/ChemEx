"""TODO: module docstring.

# Operator basis: # {Ix, Iy, Iz, 2IxSz, 2IySz, 2IzSz}{a, b, c}

"""

import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.three_state import exchange_model

# yapf: disable
indexes = [0, 1, 2, 6, 7, 14,
           15, 16, 17, 21, 22, 29,
           30, 31, 32, 36, 37, 44]
mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_r2a_i_a = ref.mat_r2a_i_a[mesh]
mat_r1a_a = ref.mat_r1a_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]
mat_j_a = ref.mat_j_a[mesh]
mat_etaxy_i_a = ref.mat_etaxy_i_a[mesh]
mat_etaz_i_a = ref.mat_etaz_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_r2a_i_b = ref.mat_r2a_i_b[mesh]
mat_r1a_b = ref.mat_r1a_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]
mat_j_b = ref.mat_j_b[mesh]
mat_etaxy_i_b = ref.mat_etaxy_i_b[mesh]
mat_etaz_i_b = ref.mat_etaz_i_b[mesh]

mat_r2_i_c = ref.mat_r2_i_c[mesh]
mat_r1_i_c = ref.mat_r1_i_c[mesh]
mat_r2a_i_c = ref.mat_r2a_i_c[mesh]
mat_r1a_c = ref.mat_r1a_c[mesh]
mat_omega_i_c = ref.mat_omega_i_c[mesh]
mat_j_c = ref.mat_j_c[mesh]
mat_etaxy_i_c = ref.mat_etaxy_i_c[mesh]
mat_etaz_i_c = ref.mat_etaz_i_c[mesh]
mat_omega1x_i = ref.mat_omega1x_i[mesh]
mat_omega1y_i = ref.mat_omega1y_i[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]
mat_kbc = ref.mat_kbc[mesh]
mat_kcb = ref.mat_kcb[mesh]
mat_kac = ref.mat_kac[mesh]
mat_kca = ref.mat_kca[mesh]

index_iz = [2, 8, 14]
index_iz_a = [2]

index_2izsz = [5, 11, 17]
index_2izsz_a = [5]

p_180x_s = np.diag([1.0, 1.0, 1.0, -1.0, -1.0, -1.0,
                    1.0, 1.0, 1.0, -1.0, -1.0, -1.0,
                    1.0, 1.0, 1.0, -1.0, -1.0, -1.0])

p_180y_s = p_180x_s[:]

p_180x_i_perfect = np.diag([1.0, -1.0, -1.0, 1.0, -1.0, -1.0,
                            1.0, -1.0, -1.0, 1.0, -1.0, -1.0,
                            1.0, -1.0, -1.0, 1.0, -1.0, -1.0])


def compute_liouvillian(
        kab=0.0, kba=0.0, kac=0.0, kca=0.0, kbc=0.0, kcb=0.0,
        r2_i_a=0.0, r1_i_a=0.0, r2a_i_a=0.0, etaxy_i_a=0.0, etaz_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, r2a_i_b=0.0, etaxy_i_b=0.0, etaz_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, r1_i_c=0.0, r2a_i_c=0.0, etaxy_i_c=0.0, etaz_i_c=0.0, omega_i_c=0.0,
        r1a_a=0.0, j_a=0.0,
        r1a_b=0.0, j_b=0.0,
        r1a_c=0.0, j_c=0.0,
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

        mat_r2a_i_a * r2a_i_a +
        mat_r2a_i_b * r2a_i_b +
        mat_r2a_i_c * r2a_i_c +

        mat_etaxy_i_a * etaxy_i_a +
        mat_etaxy_i_b * etaxy_i_b +
        mat_etaxy_i_c * etaxy_i_c +

        mat_etaz_i_a * etaz_i_a +
        mat_etaz_i_b * etaz_i_b +
        mat_etaz_i_c * etaz_i_c +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +
        mat_omega_i_c * omega_i_c +

        mat_r1a_a * r1a_a +
        mat_r1a_b * r1a_b +
        mat_r1a_c * r1a_c +

        mat_j_a * j_a +
        mat_j_b * j_b +
        mat_j_c * j_c +

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


def compute_equilibrium_2izsz(pa=0.0, pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((18, 1))
    mag0[index_2izsz] = [[pa], [pb], [pc]]

    return mag0


def compute_2izsz_a(pa=0.0, **kwargs):
    """Compute the equilibrium magnetization for anti-TROSY."""
    mag0 = np.zeros((18, 1))
    mag0[index_2izsz] = [[pa], [0.0], [0.0]]

    return mag0


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
        'r2a_i_a': parameters.ParameterName('r2a_a', **kw2).to_full_name(),
        'etaxy_i_a': parameters.ParameterName('etaxy_a', **kw2).to_full_name(),
        'etaz_i_a': parameters.ParameterName('etaz_a', **kw2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kw2).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kw2).to_full_name(),
        'r2a_i_b': parameters.ParameterName('r2a_b', **kw2).to_full_name(),
        'etaxy_i_b': parameters.ParameterName('etaxy_b', **kw2).to_full_name(),
        'etaz_i_b': parameters.ParameterName('etaz_b', **kw2).to_full_name(),
        'r2_i_c': parameters.ParameterName('r2_c', **kw2).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kw2).to_full_name(),
        'r2a_i_c': parameters.ParameterName('r2a_c', **kw2).to_full_name(),
        'etaxy_i_c': parameters.ParameterName('etaxy_c', **kw2).to_full_name(),
        'etaz_i_c': parameters.ParameterName('etaz_c', **kw2).to_full_name(),
    })

    kw5 = {'temperature': temperature, 'nuclei': nuclei.assignment}

    map_names.update({
        'j_a': parameters.ParameterName('j_a', **kw5).to_full_name(),
        'j_b': parameters.ParameterName('j_b', **kw5).to_full_name(),
        'j_c': parameters.ParameterName('j_c', **kw5).to_full_name(),
    })

    kw6 = {'temperature': temperature, 'nuclei': nuclei.assignment, 'h_larmor_frq': h_larmor_frq}

    map_names.update({
        'r1a_a': parameters.ParameterName('r1a_a', **kw6).to_full_name(),
        'r1a_b': parameters.ParameterName('r1a_b', **kw6).to_full_name(),
        'r1a_c': parameters.ParameterName('r1a_c', **kw6).to_full_name(),
    })

    params.add_many(
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_i_ac'], 0.0, False, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['r2a_i_a'], 15.0, False, 0.0, None, None),
        (map_names['etaxy_i_a'], 0.0, True, None, None, None),
        (map_names['etaz_i_a'], 0.0, True, None, None, None),
        (map_names['r1a_a'], 2.0, False, 0.0, None, None),
        (map_names['j_a'], -93.0, False, None, None, None), )

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    cs_i_c = '{cs_i_a} + {dw_i_ac}'.format(**map_names)
    r1_i_c = r1_i_b = map_names['r1_i_a']
    r2_i_c = r2_i_b = map_names['r2_i_a']
    r1a_c = r1a_b = map_names['r1a_a']
    r2a_i_c = r2a_i_b = map_names['r2a_i_a']
    etaxy_i_c = etaxy_i_b = map_names['etaxy_i_a']
    etaz_i_c = etaz_i_b = map_names['etaz_i_a']
    j_c = j_b = map_names['j_a']

    params.add_many(
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
        (map_names['r1_i_b'], 1.0, None, 0.0, None, r1_i_b),
        (map_names['r2a_i_b'], 15.0, None, 0.0, None, r2a_i_b),
        (map_names['etaxy_i_b'], 0.0, None, None, None, etaxy_i_b),
        (map_names['etaz_i_b'], 0.0, None, None, None, etaz_i_b),
        (map_names['r1a_b'], 2.0, None, 0.0, None, r1a_b),
        (map_names['j_b'], -93.0, None, None, None, j_b),
        (map_names['cs_i_c'], 0.0, None, None, None, cs_i_c),
        (map_names['r2_i_c'], 10.0, None, 0.0, None, r2_i_c),
        (map_names['r1_i_c'], 1.0, None, 0.0, None, r1_i_c),
        (map_names['r2a_i_c'], 15.0, None, 0.0, None, r2a_i_c),
        (map_names['etaxy_i_c'], 0.0, None, None, None, etaxy_i_c),
        (map_names['etaz_i_c'], 0.0, None, None, None, etaz_i_c),
        (map_names['r1a_c'], 2.0, None, 0.0, None, r1a_c),
        (map_names['j_c'], -93.0, None, None, None, j_c), )

    return map_names, params
