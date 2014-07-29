from chemex.caching import lru_cache
from chemex.experiments.misc import correct_chemical_shift


@lru_cache()
def calc_observable(pb=0.0, kex=0.0, dw_h=0.0, dw_n=0.0, ppm_to_rads_h=1.0, ppm_to_rads_n=1.0):
    ''' Returns: float '''

    dw_h *= ppm_to_rads_h
    dw_n *= ppm_to_rads_n

    shift_sq = correct_chemical_shift(pb, kex, dw_n)[0]
    shift_mq = 0.5 * (correct_chemical_shift(pb, kex, dw_n + dw_h)[0] +
                      correct_chemical_shift(pb, kex, dw_n - dw_h)[0])

    return (shift_sq - shift_mq) / ppm_to_rads_n * 1e3



