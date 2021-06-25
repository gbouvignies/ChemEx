import numpy as np


def get_scale(values_ref, errors_ref, values_to_scale):
    values_err2 = values_to_scale / errors_ref ** 2
    try:
        scale = sum(values_ref * values_err2) / sum(values_to_scale * values_err2)
    except ZeroDivisionError:
        scale = 0.0
    return scale


def get_grid(values, size=400, extension=0.0):
    value_min = np.min(values)
    value_max = np.max(values)
    extra = (value_max - value_min) * extension
    return np.sort(
        np.concatenate(
            (np.linspace(value_min - extra, value_max + extra, size), values)
        )
    )
