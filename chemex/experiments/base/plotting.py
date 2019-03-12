"""The plotting module contains experiment-independent settings and
functions."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from chemex import peaks

PALETTE = {
    "Red": {
        "50": "#FFEBEE",
        "100": "#FFCDD2",
        "200": "#EF9A9A",
        "300": "#E57373",
        "400": "#EF5350",
        "500": "#F44336",
        "600": "#E53935",
        "700": "#D32F2F",
        "800": "#C62828",
        "900": "#B71C1C",
    },
    "Grey": {
        "50": "#FAFAFA",
        "100": "#F5F5F5",
        "200": "#EEEEEE",
        "300": "#E0E0E0",
        "400": "#BDBDBD",
        "500": "#9E9E9E",
        "600": "#757575",
        "700": "#616161",
        "800": "#424242",
        "900": "#212121",
    },
}


def set_lim(values, scale):
    """Provide a range that contains all the value and adds a margin."""

    values_ = values[np.isfinite(values)]
    v_min, v_max = min(values_), max(values_)
    margin = (v_max - v_min) * scale
    v_min, v_max = v_min - margin, v_max + margin

    return v_min, v_max


def group_data(dataset):
    """Group the data resonance specifically."""
    data_grouped = dict()

    for profile in dataset:
        resonance_id = profile.name
        peak = peaks.Peak(resonance_id)
        data_grouped[peak] = profile

    return data_grouped


def plot_data(data, params, path):
    """Plot all data types."""
    subsets = dict()

    for profile in data:
        subsets.setdefault(profile.plot_data, []).append(profile)

    for plot, dataset in subsets.items():
        plot(dataset, params, path)

    return


def plot_results_brute(result, best_vals=True, varlabels=None, output="results_brute"):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square value per parameter and contour
    plots for all combination of two parameters.

    Inspired by the `corner` package (https://github.com/dfm/corner.py).

    """
    npars = len(result.var_names)
    _, axes = plt.subplots(npars, npars)

    if not varlabels:
        varlabels = result.var_names
    if best_vals and isinstance(best_vals, bool):
        best_vals = result.params

    for i, par1 in enumerate(result.var_names):
        for j, par2 in enumerate(result.var_names):

            # parameter vs chi2 in case of only one parameter
            if npars == 1:
                axes.plot(result.brute_grid, result.brute_Jout, "o", ms=3)
                axes.set_ylabel(r"$\chi^{2}$")
                axes.set_xlabel(varlabels[i])
                if best_vals:
                    axes.axvline(best_vals[par1].value, ls="dashed", color="r")

            # parameter vs chi2 profile on top
            elif i == j and j < npars - 1:
                if i == 0:
                    axes[0, 0].axis("off")
                axis = axes[i, j + 1]
                red_axis = tuple([a for a in range(npars) if a != i])
                axis.plot(
                    np.unique(result.brute_grid[i]),
                    np.minimum.reduce(result.brute_Jout, axis=red_axis),
                    "o",
                    ms=3,
                )
                axis.set_ylabel(r"$\chi^{2}$")
                axis.yaxis.set_label_position("right")
                axis.yaxis.set_ticks_position("right")
                axis.set_xticks([])
                if best_vals:
                    axis.axvline(best_vals[par1].value, ls="dashed", color="r")

            # parameter vs chi2 profile on the left
            elif j == 0 and i > 0:
                axis = axes[i, j]
                red_axis = tuple([a for a in range(npars) if a != i])
                axis.plot(
                    np.minimum.reduce(result.brute_Jout, axis=red_axis),
                    np.unique(result.brute_grid[i]),
                    "o",
                    ms=3,
                )
                axis.invert_xaxis()
                axis.set_ylabel(varlabels[i])
                if i != npars - 1:
                    axis.set_xticks([])
                elif i == npars - 1:
                    axis.set_xlabel(r"$\chi^{2}$")
                if best_vals:
                    axis.axhline(best_vals[par1].value, ls="dashed", color="r")

            # contour plots for all combinations of two parameters
            elif j > i:
                axis = axes[j, i + 1]
                red_axis = tuple([a for a in range(npars) if a != i and a != j])
                X, Y = np.meshgrid(
                    np.unique(result.brute_grid[i]), np.unique(result.brute_grid[j])
                )
                lvls1 = np.linspace(
                    result.brute_Jout.min(),
                    np.median(result.brute_Jout) / 2.0,
                    7,
                    dtype="int",
                )
                lvls2 = np.linspace(
                    np.median(result.brute_Jout) / 2.0,
                    np.median(result.brute_Jout),
                    3,
                    dtype="int",
                )
                lvls = np.unique(np.concatenate((lvls1, lvls2)))
                axis.contourf(
                    X.T,
                    Y.T,
                    np.minimum.reduce(result.brute_Jout, axis=red_axis),
                    lvls,
                    norm=LogNorm(),
                )
                axis.set_yticks([])
                if best_vals:
                    axis.axvline(best_vals[par1].value, ls="dashed", color="r")
                    axis.axhline(best_vals[par2].value, ls="dashed", color="r")
                    axis.plot(best_vals[par1].value, best_vals[par2].value, "rs", ms=3)
                if j != npars - 1:
                    axis.set_xticks([])
                elif j == npars - 1:
                    axis.set_xlabel(varlabels[i])
                if j - i >= 2:
                    axes[i, j].axis("off")

    plt.savefig(f"{output}")
