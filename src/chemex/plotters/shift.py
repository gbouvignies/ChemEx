from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.figure import Figure

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename
from chemex.printers.data import ShiftUnit

_RED700 = "#D32F2F"
_AXIS_UNITS: dict[ShiftUnit, str] = {
    "PPB": "ppb",
    "PPM": "ppm",
}


def _axis_unit_label(unit: ShiftUnit) -> str:
    if unit not in _AXIS_UNITS:
        msg = f"Unsupported shift unit {unit!r}; expected one of {tuple(_AXIS_UNITS)}."
        raise ValueError(msg)
    return _AXIS_UNITS[unit]


def plot_shift(
    name_pdf: Path,
    data: Data,
    *,
    unit: ShiftUnit = "PPB",
) -> None:
    unit_label = _axis_unit_label(unit)
    fig = Figure()
    ax = fig.subplots(1, 1)
    fig.align_labels()
    ax.errorbar(data.calc, data.exp, yerr=data.err, fmt=".", color=_RED700)
    val_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    val_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_aspect("equal", "box")
    ax.plot([val_min, val_max], [val_min, val_max], color="k", linewidth=0.5, zorder=1)
    ax.set_xlabel(rf"$Δδ_\mathregular{{fit}}$ ({unit_label})")
    ax.set_ylabel(rf"$Δδ_\mathregular{{exp}}$ ({unit_label})")
    fig.savefig(str(name_pdf))


def create_plot_data(profiles: list[Profile]) -> Data:
    exp = np.array([profile.data.exp for profile in profiles]).flatten()
    err = np.array([profile.data.err for profile in profiles]).flatten()
    calc = np.array([profile.data.calc for profile in profiles]).flatten()

    data_plot = Data(exp=exp, err=err)
    data_plot.calc = calc

    return data_plot


class ShiftPlotter:
    def __init__(
        self,
        filename: Path,
        *,
        unit: ShiftUnit = "PPB",
        **_extra: Any,
    ) -> None:
        _axis_unit_label(unit)
        self.filename = filename
        self.unit: ShiftUnit = unit

    def plot(self, output_stem: Path, profiles: list[Profile]) -> None:
        name_pdf = output_stem.with_suffix(".pdf")

        print_plot_filename(name_pdf, extra=False)

        data_plot = create_plot_data(profiles)
        plot_shift(name_pdf, data_plot, unit=self.unit)

    def plot_simulation(self, output_stem: Path, profiles: list[Profile]) -> None:
        self.plot(output_stem, profiles)
