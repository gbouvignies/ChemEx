
from chemex.caching import lru_cache
from chemex.experiments.util import calculate_shift_ex_2st


@lru_cache()
def calc_observable(pb=0.0, kex=0.0, dw_n=0.0, ppm_to_rads_n_1=1.0,
                    ppm_to_rads_n_2=1.0):
    """Returns: float"""

    dw_n_1 = dw_n * ppm_to_rads_n_1
    dw_n_2 = dw_n * ppm_to_rads_n_2

    shift_sq_1 = calculate_shift_ex_2st(pb, kex, dw_n_1)[0] / ppm_to_rads_n_1
    shift_sq_2 = calculate_shift_ex_2st(pb, kex, dw_n_2)[0] / ppm_to_rads_n_2

    return (shift_sq_1 - shift_sq_2) * 1e3



