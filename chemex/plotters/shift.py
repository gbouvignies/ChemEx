from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.figure import Figure

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename

_RED700 = "#D32F2F"


def plot_shift(name_pdf: Path, data: Data):
    fig = Figure()
    ax = fig.subplots(1, 1)
    fig.align_labels()
    ax.errorbar(data.calc, data.exp, yerr=data.err, fmt=".", color=_RED700)
    val_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    val_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_aspect("equal", "box")
    ax.plot([val_min, val_max], [val_min, val_max], color="k", linewidth=0.5, zorder=1)
    ax.set_xlabel(r"$Δδ_\mathregular{fit}$ (ppb)")
    ax.set_ylabel(r"$Δδ_\mathregular{exp}$ (ppb)")
    fig.savefig(name_pdf)


def create_plot_data(profiles: list[Profile]) -> Data:

    exp = np.array([profile.data.exp for profile in profiles]).flatten()
    err = np.array([profile.data.err for profile in profiles]).flatten()
    calc = np.array([profile.data.calc for profile in profiles]).flatten()

    data_plot = Data(exp=exp, err=err)
    data_plot.calc = calc

    return data_plot


class ShiftPlotter:
    def __init__(self, filename: Path, **_extra: Any):
        self.filename = filename

    def plot(self, path: Path, profiles: list[Profile]) -> None:

        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")

        print_plot_filename(name_pdf, extra=False)

        data_plot = create_plot_data(profiles)
        plot_shift(name_pdf, data_plot)
