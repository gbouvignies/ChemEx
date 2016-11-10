import lmfit
import numpy as np

from chemex import parameters
from chemex.bases import ref
from chemex.bases.two_state import exchange_model

indexes = [0, 1,
           15, 16]

mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]

index_ix = [0, 2]
index_ix_a = [0]

p_180x_i = np.diag([1.0, -1.0,
                    1.0, -1.0])

p_180y_i = np.diag([-1.0, 1.0,
                    -1.0, 1.0])


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        r2_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, omega_i_b=0.0,
        **kwargs
):
    """TODO: make docstring
    """

    pa = 1.0 - pb

    kab, kba = kex_ab * np.array([pb, pa])

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian


def compute_equilibrium_x(pb=0.0, **kwargs):
    return np.array([[1.0 - pb, 0.0, pb, 0.0]]).T


def create_default_params(model=None, temperature=None, nuclei=None, h_larmor_frq=None, p_total=None, l_total=None):
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs2 = {'temperature': temperature, 'nuclei': nuclei}
    kwargs3 = {'temperature': temperature, 'nuclei': nuclei, 'h_larmor_frq': h_larmor_frq}

    map_names = {
        'pb'     : parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'kex_ab' : parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),

        'dw_i_ab': parameters.ParameterName('dw_ab', **kwargs2).to_full_name(),

        'cs_i_a' : parameters.ParameterName('cs_a', **kwargs2).to_full_name(),
        'r2_i_a' : parameters.ParameterName('r2_a', **kwargs3).to_full_name(),

        'cs_i_b' : parameters.ParameterName('cs_b', **kwargs2).to_full_name(),
        'r2_i_b' : parameters.ParameterName('r2_b', **kwargs3).to_full_name(),
    }

    cs_i_b = '{cs_i_a} + {dw_i_ab}'.format(**map_names)
    r2_i_b = map_names['r2_i_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.025, True, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),

        (map_names['dw_i_ab'], 0.0, True, None, None, None),
        (map_names['dw_i_ac'], 0.0, True, None, None, None),

        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),

        (map_names['cs_i_b'], 0.0, True, None, None, cs_i_b),
        (map_names['r2_i_b'], 10.0, None, 0.0, None, r2_i_b),
    )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total, l_total)

    return map_names, params
