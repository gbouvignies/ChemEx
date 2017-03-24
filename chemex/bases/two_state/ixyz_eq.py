"""TODO: module docstring."""

import lmfit
import numpy as np
from scipy import linalg

from chemex import parameters
from chemex.bases import ref_eq as ref
from chemex.bases.two_state import exchange_model

# yapf: disable
indexes = [0, 1, 2,
           15, 16, 17,
           60]
mesh = np.ix_(indexes, indexes)

mat_r2_i_a = ref.mat_r2_i_a[mesh]
mat_r1_i_a = ref.mat_r1_i_a[mesh]
mat_omega_i_a = ref.mat_omega_i_a[mesh]
mat_eq_i_a = ref.mat_eq_i_a[mesh]

mat_r2_i_b = ref.mat_r2_i_b[mesh]
mat_r1_i_b = ref.mat_r1_i_b[mesh]
mat_omega_i_b = ref.mat_omega_i_b[mesh]
mat_eq_i_b = ref.mat_eq_i_b[mesh]

mat_omega1x_i = ref.mat_omega1x_i[mesh]
mat_omega1y_i = ref.mat_omega1y_i[mesh]

mat_kab = ref.mat_kab[mesh]
mat_kba = ref.mat_kba[mesh]

index_iz = [2, 5]
index_iz_eq = [2, 5, 6]
index_iz_a = [2]
index_iz_b = [5]


def compute_liouvillian(
        pb=0.0, kex_ab=0.0,
        r2_i_a=0.0, r1_i_a=0.0, omega_i_a=0.0,
        r2_i_b=0.0, r1_i_b=0.0, omega_i_b=0.0,
        omega1x_i=0.0, omega1y_i=0.0,
        **kwargs):
    """Compute the Liouvillian."""
    kab, kba = kex_ab * np.array([pb, 1.0 - pb])

    liouvillian = (
        mat_r2_i_a * r2_i_a +
        mat_r2_i_b * r2_i_b +

        mat_r1_i_a * r1_i_a +
        mat_r1_i_b * r1_i_b +

        mat_eq_i_a * r1_i_a * (1.0 - pb) +
        mat_eq_i_b * r1_i_b * pb +

        mat_omega_i_a * omega_i_a +
        mat_omega_i_b * omega_i_b +

        mat_omega1x_i * omega1x_i +
        mat_omega1y_i * omega1y_i +

        mat_kab * kab +
        mat_kba * kba
    )

    return liouvillian
# yapf: enable


def compute_equilibrium_after_d1(pb=0.0, kex_ab=0.0, r1_i_a=0.0, r1_i_b=0.0, time_d1=0.0, **kwargs):
    """Compute the equilibrium magnetization."""
    mag0 = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]).T
    liouvillian_free = compute_liouvillian(pb=pb, kex_ab=kex_ab, r1_i_a=r1_i_a, r1_i_b=r1_i_b)
    return linalg.expm(liouvillian_free * time_d1).dot(mag0)


def create_default_params(model=None,
                          temperature=None,
                          nuclei=None,
                          h_larmor_frq=None,
                          p_total=None,
                          l_total=None):
    """Create the default experimental and fitting parameters."""
    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}
    kwargs2 = {'temperature': temperature, 'nuclei': nuclei}
    kwargs3 = {'temperature': temperature, 'nuclei': nuclei, 'h_larmor_frq': h_larmor_frq}

    map_names = {
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
        'cs_i_a': parameters.ParameterName('cs_a', **kwargs2).to_full_name(),
        'r2_i_a': parameters.ParameterName('r2_a', **kwargs3).to_full_name(),
        'r1_i_a': parameters.ParameterName('r1_a', **kwargs3).to_full_name(),
        'cs_i_b': parameters.ParameterName('cs_b', **kwargs2).to_full_name(),
        'r2_i_b': parameters.ParameterName('r2_b', **kwargs3).to_full_name(),
        'r1_i_b': parameters.ParameterName('r1_b', **kwargs3).to_full_name(),
    }

    r1_i_b = map_names['r1_i_a']

    params = lmfit.Parameters()

    params.add_many(
        # Name, Value, Vary, Min, Max, Expr
        (map_names['pb'], 0.05, True, 0.0, 1.0, None),
        (map_names['kex_ab'], 200.0, True, 0.0, None, None),
        (map_names['cs_i_a'], 0.0, False, None, None, None),
        (map_names['r2_i_a'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_a'], 1.0, True, 0.0, None, None),
        (map_names['cs_i_b'], 0.0, True, None, None, None),
        (map_names['r2_i_b'], 10.0, True, 0.0, None, None),
        (map_names['r1_i_b'], 1.0, True, 0.0, None, r1_i_b), )

    map_names, params = exchange_model.update_params(params, map_names, model, temperature, p_total,
                                                     l_total)

    return map_names, params
