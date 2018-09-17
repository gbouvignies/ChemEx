"""Utility functions for bases."""
import collections
import re

import numpy as np
from scipy import linalg

SIGN = np.array([1.0, -1.0])


def compute_propagator_standard(liouvillian, time):
    """TODO: function docstring."""

    return linalg.expm(liouvillian * time)


def compute_propagator_beff(liouvillian, time):
    """TODO: function docstring."""

    s, vr = linalg.eig(liouvillian)
    vri = linalg.inv(vr)

    sl = np.where(abs(s.imag) < 1e-6)[0]

    vri = vri[sl, :]
    diag = np.diag(np.exp(s[sl] * time))
    vr = vr[:, sl]

    return vr.dot(diag).dot(vri).real


def compute_propagators_from_time_series(liouvillian, times):
    """TODO: function docstring."""

    s, vr = linalg.eig(liouvillian)
    vri = linalg.inv(vr)

    propagators = {
        t: vr.dot(np.diag(np.exp(s * t))).dot(vri).real
        for t in set(times)
        if abs(t) != np.inf
    }

    return propagators


def calculate_shift_ex_2st(pb=0.0, kex_ab=0.0, domega_i_ab=0.0, r2_i_a=0.0, r2_i_b=0.0):
    """Correct major and minor peak positions in presence of exchange."""
    pa = 1.0 - pb
    kab, kba = kex_ab * np.asarray([pb, pa])

    k2ab = r2_i_a + kab
    k2ba = r2_i_b + kba - 1j * domega_i_ab

    k2ex = k2ab + k2ba
    fac = ((k2ab - k2ba) ** 2 + 4.0 * kab * kba) ** 0.5

    nu1 = (0.5 * (-k2ex + fac)).imag
    nu2 = (0.5 * (-k2ex - fac)).imag

    condition = abs(nu1) < abs(nu2)

    nu1_ = np.where(condition, nu1, nu2)
    nu2_ = np.where(condition, nu2, nu1)

    return nu1_, nu2_


def correct_intensities(
    magz_a=1.0, magz_b=0.0, pb=0.0, kex=0.0, dw=0.0, r2_i=0.0, dr2_i=0.0
):
    """Correct major and minor peak intensities in presence of exchange."""
    kab = kex * pb
    kba = kex - kab

    k2ab = r2_i + kab
    k2ba = r2_i + dr2_i - 1j * dw + kba

    k2ex = k2ab + k2ba
    fac = ((k2ab - k2ba) ** 2 + 4.0 * kab * kba) ** 0.5

    nu1, nu2 = 0.5 * (-k2ex + SIGN * fac)

    magz_a_c = +((kab - nu2 - k2ab) * magz_a + (kba + nu1 + k2ab) * magz_b) / (
        nu1 - nu2
    )

    magz_b_c = -((kab - nu1 - k2ab) * magz_a + (kba + nu2 + k2ab) * magz_b) / (
        nu1 - nu2
    )

    if abs(nu1.imag) > abs(nu2.imag):
        magz_a_c, magz_b_c = magz_b_c, magz_a_c

    if abs(np.angle(magz_a_c)) <= 0.5 * np.pi:
        magz_a_c = +abs(magz_a_c)
    else:
        magz_a_c = -abs(magz_a_c)

    if abs(np.angle(magz_b_c)) <= 0.5 * np.pi:
        magz_b_c = +abs(magz_b_c)
    else:
        magz_b_c = -abs(magz_b_c)

    return magz_a_c, magz_b_c


def parse_model(name):
    Model = collections.namedtuple("Model", ["name", "state_nb", "kind"])
    match = re.match("(\d)st\.(\w+)", name, re.IGNORECASE)
    if match:
        state_nb = int(match.group(1))
        kind = match.group(2)
        return Model(name, state_nb, kind)
    else:
        raise NameError("The model name '{}' is not correct.".format(name))
