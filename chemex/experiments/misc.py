"""
Created on Sep 1, 2011

@author: guillaume
"""

# Imports
import collections

import numpy as np
from scipy import array, pi

from chemex.caching import lru_cache

# Constants
SIGN = array([1.0, -1.0])

@lru_cache()
def correct_chemical_shift(pb=0.0, kex=0.0, dw=0.0, r_Ixy=0.0, dr_Ixy=0.0):
    """Corrects major and minor peak positions in presence of exchange."""

    kab = kex * pb
    kba = kex - kab

    k2ab = r_Ixy + kab
    k2ba = r_Ixy + dr_Ixy - 1j * dw + kba

    k2ex = k2ab + k2ba
    fac = ((k2ab - k2ba) ** 2 + 4.0 * kab * kba) ** 0.5

    nu1, nu2 = (0.5 * (-k2ex + SIGN * fac)).imag

    if abs(nu1) > abs(nu2):
        nu1, nu2 = nu2, nu1

    return nu1, nu2


def correct_intensities(Ia=1.0, Ib=0.0, pb=0.0, kex=0.0, dw=0.0, r_Ixy=0.0, dr_Ixy=0.0):
    """Corrects major and minor peak intensities in presence of exchange."""

    kab = kex * pb
    kba = kex - kab

    k2ab = r_Ixy + kab
    k2ba = r_Ixy + dr_Ixy - 1j * dw + kba

    k2ex = k2ab + k2ba
    fac = ((k2ab - k2ba) ** 2 + 4.0 * kab * kba) ** 0.5

    nu1, nu2 = 0.5 * (-k2ex + SIGN * fac)

    Ia_c = +((kab - nu2 - k2ab) * Ia + (kba + nu1 + k2ab) * Ib) / (nu1 - nu2)
    Ib_c = -((kab - nu1 - k2ab) * Ia + (kba + nu2 + k2ab) * Ib) / (nu1 - nu2)

    Ia_c, Ib_c = (Ib_c, Ia_c) if abs(nu1.imag) > abs(nu2.imag) else (Ia_c, Ib_c)

    signa = 1.0 if abs(np.angle(Ia_c)) <= 0.5 * pi else -1.0
    signb = 1.0 if abs(np.angle(Ib_c)) <= 0.5 * pi else -1.0

    return signa * abs(Ia_c), signb * abs(Ib_c)


def calc_peak_intensity(pb=0.0, kex=0.0, dw=0.0, intensities=list()):
    """
    Calculates the intensity of the volume of the peak in presence
    of chemical exchange.
    """

    Ia, Ib = intensities
    Ia, Ib = correct_intensities(Ia=Ia, Ib=Ib, pb=pb, kex=kex, dw=dw, r_Ixy=0.0)

    return Ia


def get_par(par_name, par, par_indexes, par_fixed=list()):

    if par_name in par_indexes:
        return par[par_indexes[par_name]]
    else:
        return par_fixed[par_name]


def calc_multiplet(couplings=None, mult=None):

    if mult is None:
        mult = [0.0]

    if couplings:
        couplings = list(couplings)
        j = couplings.pop()
        mult = [frq + sign * j * pi
                for frq in mult
                for sign in (1.0, -1.0)]

        return calc_multiplet(couplings, mult)

    else:
        counter = collections.Counter(mult)
        nb_component = sum(counter.values())

        return tuple((val, count / float(nb_component)) for val, count in sorted(counter.items()))


