from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import combinations, permutations, product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from lmfit.parameter import Parameters
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm
from rich.progress import track

from chemex.containers.experiments import Experiments
from chemex.messages import print_group_name, print_running_grid
from chemex.optimize.grouping import Group, create_groups
from chemex.optimize.helper import (
    calculate_statistics,
    execute_post_fit,
    execute_post_fit_groups,
    print_header,
    print_values,
)
from chemex.optimize.minimizer import minimize
from chemex.parameters import database
from chemex.typing import ArrayFloat


@dataclass
class GridResult:
    grid: dict[str, ArrayFloat]
    chisqr: ArrayFloat


def _set_param_values(
    params: Parameters,
    fnames: Iterable[str],
    values: tuple[float, ...],
):
    for fname, value in zip(fnames, values, strict=True):
        params[fname].value = value


def run_group_grid(
    group: Group,
    grid: dict[str, ArrayFloat],
    path: Path,
    fitmethod: str,
) -> GridResult:
    group_ids = group.experiments.param_ids
    group_params = database.build_lmfit_params(group_ids)
    group_grid = {
        param_id: values for param_id, values in grid.items() if param_id in group_ids
    }

    grid_ids = tuple(group_grid)
    shape = tuple(len(values) for values in group_grid.values())
    grid_size = np.prod(shape)

    basename = group.path if group.path != Path() else Path("grid")
    filename = path / "Grid" / f"{basename}.out"
    filename.parent.mkdir(parents=True, exist_ok=True)

    best_chisqr = np.inf
    best_params = group_params

    with filename.open("w", encoding="utf-8") as fileout:
        fileout.write(print_header(group_grid))

        chisqr_list: list[float] = []

        grid_values = product(*group_grid.values())

        for values in track(grid_values, total=float(grid_size), description="   "):
            _set_param_values(group_params, grid_ids, values)
            optimized_params = minimize(group.experiments, group_params, fitmethod)
            stats = calculate_statistics(group.experiments, optimized_params)
            chisqr: float = stats.get("chisqr", np.inf)
            chisqr_list.append(chisqr)
            fileout.write(print_values(values, chisqr))
            fileout.flush()

            if chisqr < best_chisqr:
                best_chisqr = chisqr
                best_params = optimized_params

    chisqr_array = np.array(chisqr_list).reshape(shape)
    database.update_from_parameters(best_params)

    return GridResult(group_grid, chisqr_array)


def _reshape_chisqr(
    grid_ref: dict[str, ArrayFloat],
    grid_result: GridResult,
) -> ArrayFloat:
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
    grid: dict[str, ArrayFloat],
    grid_results: list[GridResult],
) -> list[dict[str, ArrayFloat]]:
    grid_params = {tuple(sorted(grid_result.grid)) for grid_result in grid_results}
    grid_params_tmp = grid_params.copy()
    for params1, params2 in permutations(grid_params, 2):
        if set(params1) <= set(params2):
            grid_params_tmp.remove(params1)
    return [{key: grid[key] for key in params} for params in grid_params_tmp]


def combine_grids(
    grid: dict[str, ArrayFloat],
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


def set_params_from_grid(grids_1d: Iterable[GridResult]):
    par_values = {}
    for grid_result in grids_1d:
        id_, values = next(iter(grid_result.grid.items()))
        par_values[id_] = values[grid_result.chisqr.argmin()]
    database.set_param_values(par_values)


def make_grids_nd(
    grid: dict[str, ArrayFloat],
    grids_combined: list[GridResult],
    ndim: int,
) -> list[GridResult]:
    grids: list[GridResult] = []
    parameters = database.get_parameters(grid)
    ids = sorted(grid, key=lambda x: parameters[x].param_name)
    for selection in combinations(ids, ndim):
        for grid_result in grids_combined:
            if set(selection) <= set(grid_result.grid):
                grid_ref = {fname: grid[fname] for fname in selection}
                chisqr_final = _reshape_chisqr(grid_ref, grid_result)
                grids.append(GridResult(grid_ref, chisqr_final))
                break
    return grids


def plot_grid_1d(grids_1d: list[GridResult], path: Path) -> None:
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square values per parameter.

    """
    if not grids_1d:
        return

    with PdfPages(str(path / "grid_1d.pdf")) as pdf:
        for grid_result in grids_1d:
            parameters = database.get_parameters(grid_result.grid)
            ((id_, values),) = list(grid_result.grid.items())
            _fig, ax = plt.subplots(1, 1)
            ax.plot(values, grid_result.chisqr, "o", ms=3)
            best_value = parameters[id_].value
            if best_value is not None:
                ax.axvline(best_value, ls="dashed", color=(0.5, 0.5, 0.5))
            ax.set_xlabel(str(parameters[id_].param_name))
            ax.set_ylabel(r"$\chi^{2}$")
            pdf.savefig()
            plt.close()


def plot_grid_2d(grids_2d: list[GridResult], path: Path) -> None:
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square contour
    plots for all combination of two parameters.

    """
    if not grids_2d:
        return

    with PdfPages(str(path / "grid_2d.pdf")) as pdf:
        for grid_result in grids_2d:
            parameters = database.get_parameters(grid_result.grid)
            (id_x, values_x), (id_y, values_y) = grid_result.grid.items()
            fig, ax = plt.subplots(1, 1)
            grid_x, grid_y = np.meshgrid(values_x, values_y)
            best_value_x = parameters[id_x].value
            best_value_y = parameters[id_y].value
            if best_value_x is not None:
                ax.axvline(best_value_x, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1)
            if best_value_y is not None:
                ax.axhline(best_value_y, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1)
            cs = ax.scatter(
                grid_x,
                grid_y,
                c=grid_result.chisqr.T,
                norm=LogNorm(),
                cmap=get_cmap("viridis_r"),
            )
            cbar = fig.colorbar(cs)
            cbar.set_label(r"$\chi^{2}$")
            ax.set_xlabel(str(parameters[id_x].param_name))
            ax.set_ylabel(str(parameters[id_y].param_name))
            pdf.savefig()
            plt.close()


def run_grid(
    experiments: Experiments,
    grid_raw: list[str],
    path: Path,
    plot: str,
    fitmethod: str,
) -> None:
    print_running_grid()

    grid = database.parse_grid(grid_raw)

    groups = create_groups(experiments)

    grid_results: list[GridResult] = []
    for group in groups:
        if message := group.message:
            print_group_name(message)
        grid_result = run_group_grid(group, grid, path, fitmethod)
        grid_results.append(grid_result)

    grids_combined = combine_grids(grid, grid_results)
    grids_1d = make_grids_nd(grid, grids_combined, 1)
    grids_2d = make_grids_nd(grid, grids_combined, 2)
    set_params_from_grid(grids_1d)
    plot_grid_1d(grids_1d, path / "Grid")
    plot_grid_2d(grids_2d, path / "Grid")

    if len(groups) > 1:
        execute_post_fit_groups(experiments, path, plot)
    else:
        execute_post_fit(experiments, path, plot != "nothing")
