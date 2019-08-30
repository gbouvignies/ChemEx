def get_scale(values_ref, errors_ref, values_to_scale):
    values_err2 = values_to_scale / errors_ref ** 2
    try:
        scale = sum(values_ref * values_err2) / sum(values_to_scale * values_err2)
    except ZeroDivisionError:
        scale = 0.0
    return scale
