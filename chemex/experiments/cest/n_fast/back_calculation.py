"""
Created on March 25, 2014

@author: guillaume
"""

# Python Modules
import scipy as sc

# Local Modules
from chemex.caching import lru_cache


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0,
                         carrier=0.0, ppm_to_rads=0.0, _id=None):
    @lru_cache(20)
    def calc_observable(i0=1.0, pb=0.0, kex=0.0, dw=0.0, r_nz=1.5, r_nxy=0.0, dr_nxy=0.0, cs=0.0):

        if abs(b1_offset) >= 10000.0:

            return (1.0 - pb) * i0

        else:

            kab = kex * pb
            kba = kex - kab

            wg = (cs - carrier) * ppm_to_rads - b1_offset * 2.0 * sc.pi
            w1 = b1_frq * 2.0 * sc.pi

            dw *= ppm_to_rads

            a1 = (dw ** 2 + dr_nxy ** 2) * kba
            a2 = (wg ** 2 + w1 ** 2 + kba ** 2) * dr_nxy
            a3 = ((wg + dw) ** 2 + w1 ** 2 + kba ** 2 + dr_nxy ** 2) * kba
            a4 = dr_nxy * w1 ** 2
            rex = kab * (a1 + a2) / (a3 + a4)
            r1rho = (r_nz * wg ** 2 + (r_nxy + rex) * w1 ** 2) / (wg ** 2 + w1 ** 2)

            return wg ** 2 / (wg ** 2 + w1 ** 2) * sc.exp(-time_t1 * r1rho) * i0

    return calc_observable

