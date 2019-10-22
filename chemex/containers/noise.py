import numpy as np
from scipy import interpolate
from scipy import linalg as la
from scipy import signal
from scipy import stats


def _variance_from_duplicates(data):
    """Estimate the variance of duplicate points

    Estimate the uncertainty using the pooled standard deviation
    Ref: http://goldbook.iupac.org/html/P/P04758.html
    """
    groups = {}
    x_name, y_name, e_name = data.dtype.names
    for x, y in data[[x_name, y_name]]:
        groups.setdefault(x, []).append(y)
    variances, weights = [], []
    for group in groups.values():
        group_size = len(group)
        if group_size > 1:
            variances.append(np.var(group, ddof=1))
            weights.append(group_size - 1)
    if variances:
        variance = np.average(variances, weights=weights)
    else:
        variance = np.mean(data[e_name])
    return variance


def _variance_from_scatter(data):
    """Estimate the uncertainty in the CEST profile.

    Adapted from:
    https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise
    """
    x_name, y_name, *_ = data.dtype.names
    data_sorted = np.sort(data, order=x_name)
    values = data_sorted[y_name]
    size = values.size
    fda = [
        [1, -1],
        [1, -2, 1],
        [1, -3, 3, -1],
        [1, -4, 6, -4, 1],
        [1, -5, 10, -10, 5, -1],
        [1, -6, 15, -20, 15, -6, 1],
    ]
    fda = [np.array(a_fda) / la.norm(a_fda) for a_fda in fda]
    percents = np.array([0.05] + list(np.arange(0.1, 0.40, 0.025)))
    percent_points = stats.norm.ppf(1.0 - percents)
    sigma_est = []
    for fdai in fda:
        noisedata = sorted(signal.convolve(values, fdai, mode="valid"))
        ntrim = len(noisedata)
        if ntrim >= 2:
            xaxis = (0.5 + np.arange(1, ntrim + 1)) / (ntrim + 0.5)
            sigmas = []
            function = interpolate.interp1d(xaxis, noisedata, "linear")
            for a_perc, a_z in zip(percents, percent_points):
                try:
                    val = (function(1.0 - a_perc) - function(a_perc)) / (2.0 * a_z)
                    sigmas.append(val)
                except ValueError:
                    pass
            sigma_est.append(np.median(sigmas))
    variance = np.median(sigma_est) ** 2 / (1.0 + 15.0 * (size + 1.225) ** -1.245)
    return max(variance, 1e-8)


estimate_noise_variance = {
    "scatter": _variance_from_scatter,
    "duplicates": _variance_from_duplicates,
}
