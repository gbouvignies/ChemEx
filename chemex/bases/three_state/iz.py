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
        pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
        r1_i_a=0.0, r1_i_b=0.0, r1_i_c=0.0):
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


def compute_equilibrium_iz(pb=0.0, pc=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.zeros((3, 1))
    mag0[index_iz] = [[1.0 - pb - pc], [pb], [pc]]

    return mag0


def create_default_params(model=None,
                          temperature=None,
                          nuclei=None,
                          h_larmor_frq=None,
                          p_total=None,
                          l_total=None):
    """Create the default experimental and fitting parameters."""
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs3 = {'temperature': temperature, 'nuclei': nuclei, 'h_larmor_frq': h_larmor_frq}

    map_names = {
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'pc': parameters.ParameterName('pc', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'kex_ac': parameters.ParameterName('kex_ac', **kwargs1).to_full_name(),
        'kex_bc': parameters.ParameterName('kex_bc', **kwargs1).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
        'r1_i_c': parameters.ParameterName('r1_c', **kwargs3).to_full_name(),
    }

    r1_i_c = r1_i_b = map_names['r1_i_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.025, True, 0.0, 1.0, None),
        (map_names['pc'], 0.025, True, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['kex_ac'], 0.0, False, 0.0, None, None),
        (map_names['kex_bc'], 200.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['r1_i_b'], 1.0, None, 0.0, None, r1_i_b),
        (map_names['r1_i_c'], 1.0, None, 0.0, None, r1_i_c), )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total,
                                                     l_total)

    return map_names, params
