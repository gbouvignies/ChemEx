"""TODO: module docstring."""

import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.three_state import exchange_model

# yapf: disable
indexes = [0, 1,
           15, 16,
           30, 31]
mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]

mat_r2_i_c = ref.mat_r2_i_c[mesh]
mat_omega_i_c = ref.mat_omega_i_c[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]
mat_kbc = ref.mat_kbc[mesh]
mat_kcb = ref.mat_kcb[mesh]
mat_kac = ref.mat_kac[mesh]
mat_kca = ref.mat_kca[mesh]

index_ix = [0, 2, 4]
index_ix_a = [0]

p_180x_i = np.diag([1.0, -1.0,
                    1.0, -1.0,
                    1.0, -1.0])

p_180y_i = np.diag([-1.0, 1.0,
                    -1.0, 1.0,
                    -1.0, 1.0])


def compute_liouvillian(
        kab=0.0, kba=0.0, kac=0.0, kca=0.0, kbc=0.0, kcb=0.0,
        r2_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, omega_i_b=0.0,
        r2_i_c=0.0, omega_i_c=0.0):
    """Compute the Liouvillian."""

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +
        mat_r2_i_c * r2_i_c +

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
# yapf: enable


def compute_equilibrium_x(pa=0.0, pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    return np.array([[pa, 0.0, pb, 0.0, pc, 0.0]]).T


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
        'r2_i_b': parameters.ParameterName('r2_b', **kw2).to_full_name(),
        'r2_i_c': parameters.ParameterName('r2_c', **kw2).to_full_name(),
    })

    params.add_many(
        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_i_ac'], 0.0, False, None, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None), )

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    cs_i_c = '{cs_i_a} + {dw_i_ac}'.format(**map_names)
    r2_i_c = r2_i_b = map_names['r2_i_a']

    params.add_many(
        (map_names['cs_i_b'], 0.0, None, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
        (map_names['cs_i_c'], 0.0, None, None, None, cs_i_c),
        (map_names['r2_i_c'], 10.0, None, 0.0, None, r2_i_c), )

    return map_names, params
