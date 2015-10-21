from __future__ import absolute_import

import numpy as np

from six.moves import zip

factor = {
    2: np.sqrt(np.pi / 2.0),
    3: 2.0 / np.sqrt(np.pi),
    4: 0.5 * np.sqrt(3.0 * np.pi / 2.0),
}

def estimate_noise(profile):

    intensity_dict = {}

    for ncyc, intensity in zip(profile.ncycs, profile.val):
        intensity_dict.setdefault(ncyc, []).append(intensity)

    std_list = []
    for duplicates in intensity_dict.values():
        n_duplicates = len(duplicates)
        if n_duplicates > 1:
            std_list.append(np.std(duplicates, ddof=1) * factor[n_duplicates])

    if std_list:
        error = np.mean(std_list)
    else:
        error = np.mean(profile.err)

    return error
