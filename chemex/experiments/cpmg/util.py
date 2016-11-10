import numpy as np
from scipy import math


def correction(n):
    k = n // 2

    if n == 2 * k:
        factor = (
            math.sqrt(math.pi * (2 * k - 1) / 2.0) *
            math.factorial(2 * k - 2) / ((2 ** (2 * k - 2)) * (math.factorial(k - 1) ** 2))
        )
    else:
        factor = (
            math.sqrt(k / math.pi) *
            (2 ** (2 * k - 1)) * (math.factorial(k - 1) ** 2) / math.factorial(2 * k - 1)
        )

    return factor


def estimate_noise(profile):
    intensity_dict = {}

    for ncyc, intensity in zip(profile.ncycs, profile.val):
        intensity_dict.setdefault(ncyc, []).append(intensity)

    std_list = []
    for duplicates in intensity_dict.values():
        n_duplicates = len(duplicates)
        if n_duplicates > 1:
            std_list.append(np.std(duplicates, ddof=1) * correction(n_duplicates))

    if std_list:
        error = np.mean(std_list)
    else:
        error = np.mean(profile.err)

    return error
