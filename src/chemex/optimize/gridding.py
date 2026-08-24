from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.pyplot import get_cmap

from chemex.parameters.name import ParamName
from chemex.typing import Array


@dataclass
class GridResult:
    grid: dict[str, Array]
    chisqr: Array


def _reshape_chisqr(
    grid_ref: dict[str, Array],
    grid_result: GridResult,
) -> Array:
    keys = list(grid_result.grid)
    order = [keys.index(key) for key in grid_ref if key in keys]
    axes_to_reduce = tuple(sorted(set(range(len(keys))) - set(order)))
    order.extend(axes_to_reduce)

    # After transpose, axes are shuffled
    axes_to_reduce = tuple(order.index(index) for index in axes_to_reduce)

    chisqr_final = grid_result.chisqr.transpose(order)
    chisqr_final = np.minimum.reduce(chisqr_final, axis=axes_to_reduce)

    shape = tuple(
        len(grid_ref[key]) if key in grid_result.grid else 1 for key in grid_ref
    )
    return chisqr_final.reshape(shape)


def _get_grids(
    grid: dict[str, Array],
    grid_results: list[GridResult],
) -> list[dict[str, Array]]:
    grid_params = {tuple(sorted(grid_result.grid)) for grid_result in grid_results}
    grid_params_tmp = grid_params.copy()
    for params1, params2 in permutations(grid_params, 2):
        if set(params1) <= set(params2):
            grid_params_tmp.remove(params1)
    return [{key: grid[key] for key in params} for params in grid_params_tmp]


def combine_grids(
    grid: dict[str, Array],
    grid_results: list[GridResult],
) -> list[GridResult]:
    grids = _get_grids(grid, grid_results)

    results: list[GridResult] = []

    for grid_ref in grids:
        shape = tuple(len(values) for values in grid_ref.values())

        chisqr_sum = np.zeros(shape)

        for grid_result in grid_results:
            chisqr_final = _reshape_chisqr(grid_ref, grid_result)
            chisqr_sum += chisqr_final

        result = GridResult(grid_ref, chisqr_sum)
        results.append(result)

    return results


def make_grids_nd(
    grid: dict[str, Array],
    grids_combined: list[GridResult],
    ndim: int,
    *,
    parameter_names: Mapping[str, ParamName],
) -> list[GridResult]:
    grids: list[GridResult] = []
    ids = sorted(grid, key=parameter_names.__getitem__)
    for selection in combinations(ids, ndim):
        for grid_result in grids_combined:
            if set(selection) <= set(grid_result.grid):
                grid_ref = {fname: grid[fname] for fname in selection}
                chisqr_final = _reshape_chisqr(grid_ref, grid_result)
                grids.append(GridResult(grid_ref, chisqr_final))
                break
    return grids


def plot_grid_1d(
    grids_1d: list[GridResult],
    path: Path,
    *,
    parameter_names: Mapping[str, ParamName],
    accepted_values: Mapping[str, float],
) -> None:
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square values per parameter.

    """
    if not grids_1d:
        return

    with PdfPages(str(path / "grid_1d.pdf")) as pdf:
        for grid_result in grids_1d:
            ((id_, values),) = list(grid_result.grid.items())
            _fig, ax = plt.subplots(1, 1)
            ax.plot(values, grid_result.chisqr, "o", ms=3)
            ax.axvline(
                accepted_values[id_],
                ls="dashed",
                color=(0.5, 0.5, 0.5),
            )
            ax.set_xlabel(str(parameter_names[id_]))
            ax.set_ylabel(r"$\chi^{2}$")
            pdf.savefig()
            plt.close()


def plot_grid_2d(
    grids_2d: list[GridResult],
    path: Path,
    *,
    parameter_names: Mapping[str, ParamName],
    accepted_values: Mapping[str, float],
) -> None:
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square contour
    plots for all combination of two parameters.

    """
    if not grids_2d:
        return

    with PdfPages(str(path / "grid_2d.pdf")) as pdf:
        for grid_result in grids_2d:
            (id_x, values_x), (id_y, values_y) = grid_result.grid.items()
            fig, ax = plt.subplots(1, 1)
            grid_x, grid_y = np.meshgrid(values_x, values_y)
            ax.axvline(
                accepted_values[id_x],
                ls="dashed",
                color=(0.5, 0.5, 0.5),
                zorder=-1,
            )
            ax.axhline(
                accepted_values[id_y],
                ls="dashed",
                color=(0.5, 0.5, 0.5),
                zorder=-1,
            )
            cs = ax.scatter(
                grid_x,
                grid_y,
                c=grid_result.chisqr.T,
                norm=LogNorm(),
                cmap=get_cmap("viridis_r"),
            )
            cbar = fig.colorbar(cs)
            cbar.set_label(r"$\chi^{2}$")
            ax.set_xlabel(str(parameter_names[id_x]))
            ax.set_ylabel(str(parameter_names[id_y]))
            pdf.savefig()
            plt.close()
