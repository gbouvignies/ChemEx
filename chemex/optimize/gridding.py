import collections as co
import itertools as it

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from tqdm.contrib import itertools as tqdmit

import chemex.optimize.fitting as cof
import chemex.optimize.helper as coh


GridResult = co.namedtuple("GridResult", "grid chisqr")


def run_grid(group, grid, path, fitmethod):

    g_exp = group["experiments"]
    g_par = group["params"]

    g_grid = {fname: values for fname, values in grid.items() if fname in g_par}

    fnames = tuple(g_grid)
    shape = tuple(len(values) for values in g_grid.values())
    grid_n = np.prod(shape)

    filename = path / "Grid" / f"{str(group['path'])}.out"
    filename.parent.mkdir(parents=True, exist_ok=True)

    with filename.open("w") as fileout:

        fileout.write(coh.print_header(g_par, g_grid))

        chisqr = []

        for values in tqdmit.product(*g_grid.values(), total=grid_n):
            _set_param_values(g_par, fnames, values)
            params = cof.minimize(g_exp, g_par, fitmethod, verbose=False)
            chisqr.append(coh.calculate_statistics(g_exp, params).get("chisqr"))
            fileout.write(coh.print_values(values, chisqr[-1]))
            fileout.flush()
        print()

    chisqr = np.array(chisqr).reshape(shape)

    return GridResult(g_grid, chisqr)


def _set_param_values(params, fnames, values):
    for fname, value in zip(fnames, values):
        params[fname].value = value


def combine_grids(grid, grid_results):

    grids = _get_grids(grid, grid_results)

    results = []

    for grid_ref in grids:

        shape = tuple(len(values) for values in grid_ref.values())

        chisqr_sum = np.zeros(shape)

        for a_grid, chisqr in grid_results:
            chisqr_final = _reshape_chisqr(grid_ref, a_grid, chisqr)
            chisqr_sum += chisqr_final

        result = GridResult(grid_ref, chisqr_sum)
        results.append(result)

    return results


def _reshape_chisqr(grid_ref, grid, chisqr):
    keys = list(grid)
    order = [keys.index(key) for key in grid_ref if key in keys]
    axes_to_reduce = tuple(sorted(set(range(len(keys))) - set(order)))
    order.extend(axes_to_reduce)

    # After transpose, axes are shuffled
    axes_to_reduce = tuple(order.index(index) for index in axes_to_reduce)

    chisqr_final = chisqr.transpose(order)
    chisqr_final = np.minimum.reduce(chisqr_final, axis=axes_to_reduce)

    shape = tuple(len(grid_ref[key]) if key in grid else 1 for key in grid_ref)
    chisqr_final = chisqr_final.reshape(shape)
    return chisqr_final


def make_grids_1d(params, grid, grids_combined):
    grids_2d = []
    fnames = sorted(grid, key=lambda x: params[x].user_data["pname"])
    for pair in it.combinations(fnames, 2):
        for a_grid, chisqr in grids_combined:
            if set(pair) <= set(a_grid):
                grid_ref = {fname: grid[fname] for fname in pair}
                chisqr_final = _reshape_chisqr(grid_ref, a_grid, chisqr)
                grids_2d.append(GridResult(grid_ref, chisqr_final))
    return grids_2d


def _get_grids(grid, grid_results):
    grid_params = {tuple(sorted(grid_result.grid)) for grid_result in grid_results}
    grid_params_tmp = grid_params.copy()
    for params1, params2 in it.permutations(grid_params, 2):
        if set(params1) <= set(params2):
            grid_params_tmp.remove(params1)
    return [{key: grid[key] for key in params} for params in grid_params_tmp]


def set_params_from_grid(grids_1d, params):
    for grid, chisqr in grids_1d:
        fname, values = list(grid.items())[0]
        params[fname].value = values[chisqr.argmin()]
        params[fname].vary = False


def make_grids_nd(grid, grids_combined, params, ndim):
    grids = []
    fnames = sorted(grid, key=lambda x: params[x].user_data["pname"])
    for selection in it.combinations(fnames, ndim):
        for a_grid, chisqr in grids_combined:
            if set(selection) <= set(a_grid):
                grid_ref = {fname: grid[fname] for fname in selection}
                chisqr_final = _reshape_chisqr(grid_ref, a_grid, chisqr)
                grids.append(GridResult(grid_ref, chisqr_final))
                break
    return grids


def plot_grid_1d(grids_1d, params, path):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square values per parameter.

    """

    if not grids_1d:
        return

    with PdfPages(path / "grid_1d.pdf") as pdf:
        for grid, chisqr in grids_1d:
            fname, values = list(grid.items())[0]
            _fig, ax = plt.subplots(1, 1)
            ax.plot(values, chisqr, "o", ms=3)
            ax.axvline(params[fname].value, ls="dashed", color=(0.5, 0.5, 0.5))
            ax.set_xlabel(str(params[fname].user_data["pname"]))
            ax.set_ylabel(r"$\chi^{2}$")
            pdf.savefig()
            plt.close()


def plot_grid_2d(grids_2d, params, path):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square contour
    plots for all combination of two parameters.

    """

    if not grids_2d:
        return

    with PdfPages(path / "grid_2d.pdf") as pdf:
        for grid, chisqr in grids_2d:
            (fname_x, values_x), (fname_y, values_y) = grid.items()
            fig, ax = plt.subplots(1, 1)
            grid_x, grid_y = np.meshgrid(values_x, values_y)
            ax.axvline(
                params[fname_x].value, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1
            )
            ax.axhline(
                params[fname_y].value, ls="dashed", color=(0.5, 0.5, 0.5), zorder=-1
            )
            cs = ax.scatter(
                grid_x, grid_y, c=chisqr.T, norm=LogNorm(), cmap=cm.viridis_r
            )
            cbar = fig.colorbar(cs)
            cbar.set_label(r"$\chi^{2}$")
            ax.set_xlabel(str(params[fname_x].user_data["pname"]))
            ax.set_ylabel(str(params[fname_y].user_data["pname"]))
            pdf.savefig()
            plt.close()
