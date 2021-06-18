import collections as co
import itertools as it

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from tqdm.contrib import itertools as tqdmit

import chemex.optimize.fitting as cof
import chemex.optimize.grouping as coc
import chemex.optimize.helper as coh


GridResult = co.namedtuple("GridResult", "grid chisqr")


def run_grid(group, grid, path, fitmethod):

    g_exp = group["experiments"]
    g_par = group["params"]

    g_grid = {k: v for k, v in grid.items() if k in g_par}
    fnames = set(g_grid)
    shape = tuple(len(values) for values in g_grid.values())
    grid_n = np.prod(shape)

    filename = path / "Grid" / f"{str(group['path'])}.out"
    filename.parent.mkdir(parents=True, exist_ok=True)

    with filename.open("w") as fileout:

        fileout.write(_print_header(g_par, g_grid))

        chisqr = []

        for values in tqdmit.product(*g_grid.values(), total=grid_n):
            _set_param_values(g_par, fnames, values)
            params = cof.minimize(g_exp, g_par, fitmethod, verbose=False)
            chisqr.append(coh.calculate_statistics(g_exp, params).get("chisqr"))
            fileout.write(_print_values(values, chisqr))
            fileout.flush()

    chisqr = np.array(chisqr).reshape(shape)

    return GridResult(g_grid, chisqr)


def _set_param_values(params, fnames, values):
    for fname, value in zip(fnames, values):
        params[fname].value = value


def _print_header(params, grid):
    pnames = (str(params[fname].user_data["pname"]) for fname in grid)
    header_pnames = " ".join(f"{pname:>25s}" for pname in pnames)
    return f"# {header_pnames} {'chisqr':>25s}\n"


def _print_values(values, chisqr):
    body_values = " ".join(f"{value:>25.8e}" for value in values)
    return f"  {body_values} {chisqr[-1]:>25.8e}\n"


def combine_grids(experiments, params, grid_results, path):

    grid = {}
    for grid_result in grid_results:
        grid.update(grid_result.grid)

    for fname in grid:
        params[fname].vary = True

    reduced_grids = []

    while grid_results:

        for group in coc.group_pnames(experiments, params):

            g_names = get_common_names(grid_results, group)

            g_grid = {name: grid.pop(name) for name in g_names}

            chisqr_list = []
            for grid_result in grid_results:
                axis = tuple(
                    index
                    for index, name in enumerate(grid_result.grid)
                    if name not in g_grid
                )
                chisqr_red = np.minimum.reduce(grid_result.chisqr, axis=axis)
                chisqr_list.append(chisqr_red)

            reduced_grid = GridResult(g_grid, sum(chisqr_list))

            reduced_grids.append(reduced_grid)

            argmin = np.unravel_index(
                reduced_grid.chisqr.argmin(), reduced_grid.chisqr.shape
            )

            for index, fname in zip(argmin, g_grid):
                params[fname].value = g_grid[fname][index]
                params[fname].vary = False

            grid_results = trim_grid(grid_results, g_grid, argmin)

    _plot_grid(reduced_grids, params, path)

    return params


def get_common_names(grid_results, group):
    name_sets = (
        set(grid) for grid, _ in grid_results if set(grid) and set(grid) <= group
    )
    return set.intersection(*name_sets)


def trim_grid(grid_results, grid, argmin):
    grids_trimmed = []
    for grid_, chisqr in grid_results:
        axis = tuple(index for index, name in enumerate(grid_) if name in grid)
        indexes = {ax: index for ax, index in zip(axis, argmin)}
        ix = tuple(indexes.get(dim, slice(None)) for dim in range(chisqr.ndim))
        grid_trimmed = {
            name: values for name, values in grid_.items() if name not in grid
        }
        if grid_trimmed:
            chisqr_trimmed = chisqr[ix]
            grids_trimmed.append(GridResult(grid_trimmed, chisqr_trimmed))
    return grids_trimmed


def _plot_grid(grids, params, path):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square value per parameter and contour
    plots for all combination of two parameters.
    """

    path.mkdir(parents=True, exist_ok=True)

    _plot_grid_1d(grids, params, path)
    _plot_grid_2d(grids, params, path)


def _plot_grid_1d(grids, params, path):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path / "grid_1d.pdf") as pdf:
        for grid_values, chisqr in grids:
            for i, (fname, values) in enumerate(grid_values.items()):
                red_axis = tuple(a for a in range(len(grid_values)) if a != i)
                chi2_vector = np.minimum.reduce(chisqr, axis=red_axis)
                _fig, ax = plt.subplots(1, 1)
                ax.plot(values, chi2_vector, "o", ms=3)
                min_index = np.argmin(chi2_vector)
                ax.axvline(values[min_index], ls="dashed", color=(0.5, 0.5, 0.5))
                ax.set_xlabel(str(params[fname].user_data["pname"]))
                ax.set_ylabel(r"$\chi^{2}$")
                pdf.savefig()
                plt.close()


def _plot_grid_2d(grids, params, path):

    with PdfPages(path / "grid_2d.pdf") as pdf:
        for grid_values, chisqr in grids:
            combinations = it.combinations(enumerate(grid_values.items()), 2)
            for (i, (n1, v1)), (j, (n2, v2)) in combinations:
                fig, ax = plt.subplots(1, 1)
                grid_x, grid_y = np.meshgrid(v1, v2)
                red_axis = tuple(a for a in range(len(grid_values)) if a not in (i, j))
                chi2_matrix = np.minimum.reduce(chisqr, axis=red_axis).T
                min_indexes = np.unravel_index(
                    np.argmin(chi2_matrix), chi2_matrix.shape
                )
                min_x = grid_x[min_indexes]
                min_y = grid_y[min_indexes]
                ax.axvline(min_x, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1)
                ax.axhline(min_y, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1)
                cs = ax.scatter(
                    grid_x,
                    grid_y,
                    c=chi2_matrix,
                    norm=LogNorm(),
                    cmap=cm.viridis_r,
                )
                cbar = fig.colorbar(cs)
                cbar.set_label(r"$\chi^{2}$")
                ax.set_xlabel(str(params[n1].user_data["pname"]))
                ax.set_ylabel(str(params[n2].user_data["pname"]))
                pdf.savefig()
                plt.close()
