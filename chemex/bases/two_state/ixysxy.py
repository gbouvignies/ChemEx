"""TODO: module docstring.

# Operator basis: # {2IxSx, 2IxSy, 2IySx, 2IySy}{a, b}

"""

import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.two_state import exchange_model

# yapf: disable
indexes = [10, 11, 12, 13,
           25, 26, 27, 28, ]
mesh = np.ix_(indexes, indexes)

mat_omega_i_a = ref.mat_omega_i_a[mesh]
mat_omega_s_a = ref.mat_omega_s_a[mesh]
mat_r2_mq_a = ref.mat_r2_mq_a[mesh]
mat_mu_mq_a = ref.mat_mu_mq_a[mesh]

mat_omega_i_b = ref.mat_omega_i_b[mesh]
mat_omega_s_b = ref.mat_omega_s_b[mesh]
mat_r2_mq_b = ref.mat_r2_mq_b[mesh]
mat_mu_mq_b = ref.mat_mu_mq_b[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]

index_iysx = [2, 6]
index_iysx_a = [2]

p_180x_i = np.diag([1.0, 1.0, -1.0, -1.0,
                    1.0, 1.0, -1.0, -1.0])

p_180y_i = np.diag([-1.0, -1.0, 1.0, 1.0,
                    -1.0, -1.0, 1.0, 1.0])

p_180x_s = np.diag([1.0, -1.0, 1.0, -1.0,
                    1.0, -1.0, 1.0, -1.0])

p_180y_s = np.diag([-1.0, 1.0, -1.0, 1.0,
                    -1.0, 1.0, -1.0, 1.0])


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        omega_i_a=0.0, omega_s_a=0.0,
        omega_i_b=0.0, omega_s_b=0.0,
        r2_mq_a=0.0, mu_mq_a=0.0,
        r2_mq_b=0.0, mu_mq_b=0.0,
        **kwargs):
    """Compute the Liouvillian."""
    pa = 1.0 - pb

    kab, kba = kex_ab * np.array([pb, pa])

    liouvillian = (
        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_omega_s_a * omega_s_a +
        mat_omega_s_b * omega_s_b +

        mat_r2_mq_a * r2_mq_a +
        mat_r2_mq_b * r2_mq_b +

        mat_mu_mq_a * mu_mq_a +
        mat_mu_mq_b * mu_mq_b +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
# yapf: enable


def compute_equilibrium_iysx(pb=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((8, 1))
    mag0[index_iysx] = [[1.0 - pb], [pb]]

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
        'cs_i_a': parameters.ParameterName('cs_a', **kw1).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kw1).to_full_name(),
    })

    kw3 = {'temperature': temperature, 'nuclei': resonance_s['name']}

    map_names.update({
        'dw_s_ab': parameters.ParameterName('dw_ab', **kw3).to_full_name(),
        'cs_s_a': parameters.ParameterName('cs_a', **kw3).to_full_name(),
        'cs_s_b': parameters.ParameterName('cs_b', **kw3).to_full_name(),
    })

    kw6 = {'temperature': temperature, 'nuclei': nuclei.assignment, 'h_larmor_frq': h_larmor_frq}

    map_names.update({
        'r2_mq_a': parameters.ParameterName('r2_mq_a', **kw6).to_full_name(),
        'mu_mq_a': parameters.ParameterName('mu_mq_a', **kw6).to_full_name(),
        'r2_mq_b': parameters.ParameterName('r2_mq_b', **kw6).to_full_name(),
        'mu_mq_b': parameters.ParameterName('mu_mq_b', **kw6).to_full_name(),
    })

    params.add_many(
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_s_ab'], 0.0, False, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['cs_s_a'], 0.0, False, None, None, None),
        (map_names['mu_mq_a'], 0.0, False, None, None, None),
        (map_names['r2_mq_a'], 10.0, False, 0.0, None, None), )

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    cs_s_b = '{cs_s_a} + {dw_s_ab}'.format(**map_names)
    mu_mq_b = map_names['mu_mq_a']
    r2_mq_b = map_names['r2_mq_a']

    params.add_many(
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['cs_s_b'], 0.0, None, None, None, cs_s_b),
        (map_names['mu_mq_b'], 0.0, None, None, None, mu_mq_b),
        (map_names['r2_mq_b'], 10.0, None, 0.0, None, r2_mq_b), )

    return map_names, params
