"""Plots for exact profiled chi-square GRID surfaces."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
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
    """One plot-ready exact profile surface."""

    grid: dict[str, Array]
    chisqr: Array


def plot_grid_1d(
    grids_1d: list[GridResult],
    path: Path,
    *,
    parameter_names: Mapping[str, ParamName],
    accepted_values: Mapping[str, float],
) -> None:
    """Visualize exact one-dimensional profiled chi-square surfaces."""
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
    """Visualize exact two-dimensional profiled chi-square surfaces."""
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
