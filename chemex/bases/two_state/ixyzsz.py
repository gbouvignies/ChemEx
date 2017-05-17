"""TODO: module docstring.

# Operator basis: # {Ix, Iy, Iz, 2IxSz, 2IySz, 2IzSz}{a, b,}

"""

import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.two_state import exchange_model

# yapf: disable
indexes = [0, 1, 2, 6, 7, 14,
           15, 16, 17, 21, 22, 29]
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

mat_omega1x_i = ref.mat_omega1x_i[mesh]
mat_omega1y_i = ref.mat_omega1y_i[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]

index_iz = [2, 8]
index_iz_a = [2]

index_2izsz = [5, 11]
index_2izsz_a = [5]

p_180x_s = np.diag([1.0, 1.0, 1.0, -1.0, -1.0, -1.0,
                    1.0, 1.0, 1.0, -1.0, -1.0, -1.0])

p_180y_s = p_180x_s[:]

p_180x_i_perfect = np.diag([1.0, -1.0, -1.0, 1.0, -1.0, -1.0,
                            1.0, -1.0, -1.0, 1.0, -1.0, -1.0])


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        r2_i_a=0.0, r1_i_a=0.0, r2a_i_a=0.0, etaxy_i_a=0.0, etaz_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, r2a_i_b=0.0, etaxy_i_b=0.0, etaz_i_b=0.0, omega_i_b=0.0,
        r1a_a=0.0, j_a=0.0,
        r1a_b=0.0, j_b=0.0,
        omega1x_i=0.0, omega1y_i=0.0,
        **kwargs):
    """Compute the Liouvillian."""
    pa = 1.0 - pb

    kab, kba = kex_ab * np.array([pb, pa])

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +

        mat_r2a_i_a * r2a_i_a +
        mat_r2a_i_b * r2a_i_b +

        mat_etaxy_i_a * etaxy_i_a +
        mat_etaxy_i_b * etaxy_i_b +

        mat_etaz_i_a * etaz_i_a +
        mat_etaz_i_b * etaz_i_b +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_r1a_a * r1a_a +
        mat_r1a_b * r1a_b +

        mat_j_a * j_a +
        mat_j_b * j_b +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
# yapf: enable


def compute_equilibrium_2izsz(pb=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((12, 1))
    mag0[index_2izsz] = [[1.0 - pb], [pb]]

    return mag0


def compute_2izsz_a(pb=0.0, **kwargs):
    """Compute the equilibrium magnetization for anti-TROSY."""
    mag0 = np.zeros((12, 1))
    mag0[index_2izsz] = [[1.0 - pb], [0.0]]

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
    kwargs4 = {'temperature': temperature, 'nuclei': nuclei.assignment}
    kwargs5 = {
        'temperature': temperature,
        'nuclei': nuclei.assignment,
        'h_larmor_frq': h_larmor_frq
    }

    map_names = {
        'pa': parameters.ParameterName('pa', **kwargs1).to_full_name(),
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'dw_i_ab': parameters.ParameterName('dw_ab', **kwargs2).to_full_name(),
        'cs_i_a': parameters.ParameterName('cs_a', **kwargs2).to_full_name(),
        'r2_i_a': parameters.ParameterName('r2_a', **kwargs3).to_full_name(),
        'r2a_i_a': parameters.ParameterName('r2a_a', **kwargs3).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'r1a_a': parameters.ParameterName('r1a_a', **kwargs5).to_full_name(),
        'etaxy_i_a': parameters.ParameterName('etaxy_a', **kwargs3).to_full_name(),
        'etaz_i_a': parameters.ParameterName('etaz_a', **kwargs3).to_full_name(),
        'j_a': parameters.ParameterName('j_a', **kwargs4).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kwargs2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kwargs3).to_full_name(),
        'r2a_i_b': parameters.ParameterName('r2a_b', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
        'r1a_b': parameters.ParameterName('r1a_b', **kwargs5).to_full_name(),
        'etaxy_i_b': parameters.ParameterName('etaxy_b', **kwargs3).to_full_name(),
        'etaz_i_b': parameters.ParameterName('etaz_b', **kwargs3).to_full_name(),
        'j_b': parameters.ParameterName('j_b', **kwargs4).to_full_name(),
    }

    pa = '1.0 - {pb}'.format(**map_names)
    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    r1_i_b = map_names['r1_i_a']
    r2_i_b = map_names['r2_i_a']
    r2a_i_b = map_names['r2a_i_a']
    r1a_b = map_names['r1a_a']
    etaxy_i_b = map_names['etaxy_i_a']
    etaz_i_b = map_names['etaz_i_a']
    j_b = map_names['j_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.05, True, 0.0, 1.0, None),
        (map_names['pa'], 0.00, True, 0.0, 1.0, pa),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1a_a'], 3.0, False, 0.0, None, None),
        (map_names['etaxy_i_a'], 0.0, False, None, None, None),
        (map_names['etaz_i_a'], 0.0, False, None, None, None),
        (map_names['j_a'], 135.0, False, None, None, None),
        (map_names['r2a_i_a'], 3.0, False, 0.0, None, None),
        (map_names['r1_i_a'], 3.0, False, 0.0, None, None),
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
        (map_names['r1a_b'], 3.0, None, 0.0, None, r1a_b),
        (map_names['etaxy_i_b'], 0.0, None, None, None, etaxy_i_b),
        (map_names['etaz_i_b'], 0.0, None, None, None, etaz_i_b),
        (map_names['j_b'], 135.0, None, None, None, j_b),
        (map_names['r2a_i_b'], 3.0, None, 0.0, None, r2a_i_b),
        (map_names['r1_i_b'], 3.0, None, 0.0, None, r1_i_b), )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total,
                                                     l_total)

    return map_names, params
