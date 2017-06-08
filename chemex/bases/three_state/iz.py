"""TODO: module docstring."""

import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.three_state import exchange_model

indexes = [2, 17, 32]
mesh = np.ix_(indexes, indexes)

mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_r1_i_c = ref.mat_r1_i_c[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]
mat_kbc = ref.mat_kbc[mesh]
mat_kcb = ref.mat_kcb[mesh]
mat_kac = ref.mat_kac[mesh]
mat_kca = ref.mat_kca[mesh]

index_iz = [0, 1, 2]
index_iz_a = [0]


# yapf: disable
def compute_liouvillian(
        kab=0.0, kba=0.0, kac=0.0, kca=0.0, kbc=0.0, kcb=0.0,
        r1_i_a=0.0, r1_i_b=0.0, r1_i_c=0.0):
    """Compute the Liouvillian."""

    liouvillian = (

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +
        mat_r1_i_c * r1_i_c +

        mat_kab * kab +
        mat_kba * kba +
        mat_kbc * kbc +
        mat_kcb * kcb +
        mat_kac * kac +
        mat_kca * kca
    )

    return liouvillian
# yapf: enable


def compute_equilibrium_iz(pa=0.0, pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((3, 1))
    mag0[index_iz] = [[pa], [pb], [pc]]

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

    kw2 = {'temperature': temperature, 'nuclei': resonance_i['name'], 'h_larmor_frq': h_larmor_frq}

    map_names.update({
        'r1_i_a': parameters.ParameterName('r1_a', **kw2).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kw2).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kw2).to_full_name(),
    })

    params.add_many((map_names['r1_i_a'], 1.0, True, 0.0, None, None), )

    r1_i_c = r1_i_b = map_names['r1_i_a']

    params.add_many(
        (map_names['r1_i_b'], 1.0, None, 0.0, None, r1_i_b),
        (map_names['r1_i_c'], 1.0, None, 0.0, None, r1_i_c), )

    return map_names, params
