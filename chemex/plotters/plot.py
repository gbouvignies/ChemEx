from __future__ import annotations

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from chemex.containers.data import Data

_RED100 = "#FFCDD2"
_RED300 = "#E57373"
_RED700 = "#D32F2F"


def _create_fig(name: str) -> tuple[Figure, Axes, Axes]:
    fig = Figure()
    ax1, ax2 = fig.subplots(2, 1, sharex="all", gridspec_kw={"height_ratios": [1, 4]})
    fig.align_labels()
    fig.suptitle(f"{str(name).upper()}")
    ax1.set_ylabel("Residuals")
    ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y", useMathText=True)

    return fig, ax1, ax2


def get_grid(values: np.ndarray, size: int = 400, extension: float = 0.0) -> np.ndarray:
    value_min = np.min(values)
    value_max = np.max(values)
    extra = (value_max - value_min) * extension
    return np.linspace(value_min - extra, value_max + extra, size)


def _plot_fit(data_calc: Data, ax2: Axes):
    fit_x = data_calc.metadata
    fit_y = data_calc.calc
    range_x = get_grid(fit_x, 2, 0.02)
    ax2.set_xlim(*range_x)
    ax2.plot(fit_x, fit_y, linestyle="-", color=_RED300)


def _plot_exp(data_exp: Data, ax1: Axes, ax2: Axes):
    exp_x = data_exp.metadata
    exp_y = data_exp.exp
    exp_e = abs(data_exp.err)
    res_y = data_exp.exp - data_exp.calc

    m_sel = data_exp.mask
    m_inf = exp_e.sum(axis=1) > 1e16

    s1 = m_sel & ~m_inf
    ax1.errorbar(exp_x[s1], res_y[s1], exp_e[s1].T, fmt=".", color=_RED700, zorder=3)
    ax2.errorbar(exp_x[s1], exp_y[s1], exp_e[s1].T, fmt=".", color=_RED700, zorder=3)

    s3 = ~m_sel
    ax1.errorbar(exp_x[s3], res_y[s3], exp_e[s3].T, fmt=".", color=_RED100, zorder=3)
    ax2.errorbar(exp_x[s3], exp_y[s3], exp_e[s3].T, fmt=".", color=_RED100, zorder=3)

    range1_y = ax1.get_ylim()
    range2_y = ax2.get_ylim()

    exp_e[exp_e == np.inf] = 500.0
    s2 = m_sel & m_inf
    ax1.errorbar(exp_x[s2], res_y[s2], exp_e[s2].T, fmt=".", color=_RED700, zorder=3)
    ax2.errorbar(exp_x[s2], exp_y[s2], exp_e[s2].T, fmt=".", color=_RED700, zorder=3)

    ax1.set_ylim(range1_y)
    ax2.set_ylim(range2_y)


def plot_profile(name: str, data_exp: Data, data_calc: Data) -> Figure:
    fig, ax1, ax2 = _create_fig(name)

    _plot_fit(data_calc, ax2)

    if data_exp.size:
        _plot_exp(data_exp, ax1, ax2)

    for axis in (ax1, ax2):
        axis.axhline(0, color="k", linewidth=0.5, zorder=1)

    return fig
