import numpy as np


def estimate_noise(profile):
    intensity_dict = {}

    for ncyc, intensity in zip(profile.ncycs, profile.val):
        intensity_dict.setdefault(ncyc, []).append(intensity)

    std_list = [
        np.std(duplicates, ddof=1)
        for duplicates in intensity_dict.values()
        if len(duplicates) > 1]

    if std_list:
        error = np.mean(std_list)
    else:
        error = np.mean(profile.err)

    return error
