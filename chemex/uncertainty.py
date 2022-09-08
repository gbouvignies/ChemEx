from __future__ import annotations

import contextlib
from collections import defaultdict
from collections.abc import Callable
from typing import Any

import numpy as np
from scipy import interpolate
from scipy import signal
from scipy import stats
from scipy.linalg import norm

from chemex.containers.data import Data


def _variance_from_duplicates(data: Data) -> float:
    """Estimate the variance of duplicate points.

    Estimate the uncertainty using the pooled standard deviation.
    Reference: https://goldbook.iupac.org/html/P/P04758.html

    If no duplicates are found, the average of the experimental error
    values is returned

    """
    groups: defaultdict[Any, list[float]] = defaultdict(list)
    for x, y in zip(data.metadata, data.exp):
        groups[x].append(y)
    variances, weights = [], []
    for group in groups.values():
        group_size = len(group)
        if group_size > 1:
            variances.append(np.var(group, ddof=1))
            weights.append(group_size - 1)
    if not variances:
        return float(np.mean(data.err))
    return float(np.average(variances, weights=weights))


def _variance_from_scatter(data: Data) -> float:
    """Estimate the uncertainty in the CEST profile.

    Adapted from:
    https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise

    """
    _, exp, _ = zip(*sorted(zip(data.metadata, data.exp, data.err)))
    size = len(exp)
    fda = [
        [1, -1],
        [1, -2, 1],
        [1, -3, 3, -1],
        [1, -4, 6, -4, 1],
        [1, -5, 10, -10, 5, -1],
        [1, -6, 15, -20, 15, -6, 1],
    ]
    fda = [np.array(a_fda) / norm(a_fda) for a_fda in fda]
    percents = np.array([0.05] + list(np.arange(0.1, 0.4, 0.025)))
    percent_points = stats.norm.ppf(1.0 - percents)
    sigma_est = []
    for fdai in fda:
        noisedata = sorted(signal.convolve(exp, fdai, mode="valid"))
        ntrim = len(noisedata)
        if ntrim >= 2:
            xaxis = (0.5 + np.arange(1, ntrim + 1)) / (ntrim + 0.5)
            sigmas = []
            function = interpolate.interp1d(xaxis, noisedata, "linear")
            for a_perc, a_z in zip(percents, percent_points):
                with contextlib.suppress(ValueError):
                    val = (function(1.0 - a_perc) - function(a_perc)) / (2.0 * a_z)
                    sigmas.append(val)
            sigma_est.append(np.median(sigmas))
    variance = np.median(sigma_est) ** 2 / (1.0 + 15.0 * (size + 1.225) ** -1.245)
    return np.max([variance, 1e-8])


estimate_noise_variance: dict[str, Callable[[Data], float]] = {
    "scatter": _variance_from_scatter,
    "duplicates": _variance_from_duplicates,
}
