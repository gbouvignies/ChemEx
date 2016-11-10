import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.two_state import exchange_model

indexes = [2, 17]
mesh = np.ix_(indexes, indexes)

mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]

index_iz = [0, 1]
index_iz_a = [0]


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        r1_i_a=0.0,
        r1_i_b=0.0, ):
    """TODO: make docstring
    """

    pa = 1.0 - pb
    kab, kba = kex_ab * np.asarray([pb, pa])

    liouvillian = (

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian


def compute_equilibrium_iz(pb=0.0, **kwargs):
    mag0 = np.zeros((2, 1))
    mag0[index_iz] = [[1.0 - pb], [pb]]
    return mag0


def create_default_params(model=None, temperature=None, nuclei=None, h_larmor_frq=None, p_total=None, l_total=None):
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs3 = {'temperature': temperature, 'nuclei': nuclei, 'h_larmor_frq': h_larmor_frq}

    map_names = {
        'pb'    : parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
    }

    r1_i_b = map_names['r1_i_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.05, True, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['r1_i_b'], 1.0, True, 0.0, None, r1_i_b),
    )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total, l_total)

    return map_names, params
