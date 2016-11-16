"""Utility functions for CEST experiments."""

import collections

import numpy as np
from scipy import interpolate, linalg, signal, stats


def estimate_noise(x):
    """Estimate the uncertainty in the CEST profile.
    # TODO: add reference to this method...
    """

    n = len(x)

    fda = [
        [1, -1],
        [1, -2, 1],
        [1, -3, 3, -1],
        [1, -4, 6, -4, 1],
        [1, -5, 10, -10, 5, -1],
        [1, -6, 15, -20, 15, -6, 1]
    ]

    fda = [np.array(a_fda) / linalg.norm(a_fda) for a_fda in fda]

    perc = np.array([0.05] + list(np.arange(0.1, 0.40, 0.025)))
    z = stats.norm.ppf(1.0 - perc)

    sigma_est = []

    for fdai in fda:

        noisedata = signal.convolve(x, fdai, mode='valid')
        ntrim = len(noisedata)

        if ntrim >= 2:

            noisedata.sort()

            p = 0.5 + np.arange(1, ntrim + 1)
            p /= ntrim + 0.5

            q = []

            f = interpolate.interp1d(p, noisedata, 'linear')

            for a_perc, a_z in zip(perc, z):
                try:
                    val = (f(1.0 - a_perc) - f(a_perc)) / (2.0 * a_z)
                    q.append(val)
                except ValueError:
                    pass

            sigma_est.append(np.median(q))

    noisevar = np.median(sigma_est) ** 2
    noisevar /= (1.0 + 15.0 * (n + 1.225) ** -1.245)

    return np.sqrt(noisevar)


def calc_multiplet(couplings=None, multiplet=None):
    """Calculate the multiplet pattern."""
    couplings = list(couplings)

    if multiplet is None:
        multiplet = [0.0]

    if couplings:
        j = couplings.pop()
        multiplet_updated = [frq + sign * j * np.pi for frq in multiplet for sign in (1.0, -1.0)]

        return calc_multiplet(couplings, multiplet_updated)

    else:
        counter = collections.Counter(multiplet)
        nb_component = sum(counter.values())
        multiplet = tuple((val, count / float(nb_component)) for val, count in sorted(counter.items()))

        return multiplet
