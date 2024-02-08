from __future__ import annotations

from collections import defaultdict
from contextlib import ExitStack
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename
from chemex.plotters.plot import get_grid, plot_profile
from chemex.printers.plot import PlotPrinter, data_plot_printers


def plot_exsy(file_pdf: PdfPages, name: str, data_exp: Data, data_calc: Data):
    fig = plot_profile(name, data_exp, data_calc)
    ax2 = fig.axes[1]
    ax2.set_xlabel(r"Time (s)")
    ax2.set_ylabel(r"Intensity")
    ax2.set_xlim(0.0, max(data_exp.metadata) + min(data_calc.metadata))
    file_pdf.savefig(fig)


def create_plot_data_exp(profile: Profile) -> dict[tuple[str, str], Data]:
    data = profile.data

    intensity0 = data.exp[np.argmax(abs(data.exp))]

    data_values = defaultdict(list)
    for metadata, exp, calc, err, mask in zip(
        data.metadata,
        data.exp,
        data.calc,
        data.err,
        data.mask,
        strict=True,
    ):
        states = tuple(metadata[["states1", "states2"]])
        data_values[states].append((metadata["times"], exp, calc, err, mask))

    data_exp_dict = {}
    dtype = [
        ("times", "f8"),
        ("exp", "f8"),
        ("calc", "f8"),
        ("err", "f8"),
        ("mask", "?"),
    ]
    for states, values in data_values.items():
        values_array = np.array(values, dtype=dtype)
        values_array["exp"] /= intensity0
        values_array["calc"] /= intensity0
        errors = values_array["err"] / abs(intensity0)
        errorbars = np.array([-errors, errors]).transpose()

        data_exp = Data(
            exp=values_array["exp"],
            err=errorbars,
            metadata=values_array["times"],
        )
        data_exp.calc = values_array["calc"]
        data_exp.mask = values_array["mask"]
        data_exp.sort()

        data_exp_dict[states] = data_exp

    return data_exp_dict


def create_plot_data_calc(profile: Profile) -> dict[tuple[str, str], Data]:
    spectrometer = profile.spectrometer
    data = profile.data

    times = get_grid(data.metadata["times"], 100, 0.0)
    states = np.unique(data.metadata[["states1", "states2"]])
    metadata_list = [
        (time, state1, state2) for time, (state1, state2) in product(times, states)
    ]
    dtype = [("times", "f8"), ("states1", "U1"), ("states2", "U1")]
    metadata = np.array(metadata_list, dtype=dtype)

    filler = np.zeros_like(metadata["times"])
    data_fit = Data(exp=filler, err=filler, metadata=metadata)

    intensity0 = data.exp[np.argmax(abs(data.exp))]
    scale = profile.data.scale / intensity0
    data_fit.calc = scale * profile.pulse_sequence.calculate(spectrometer, data_fit)

    data_values = defaultdict(list)
    for metadata, calc in zip(data_fit.metadata, data_fit.calc, strict=False):
        states = tuple(metadata[["states1", "states2"]])
        data_values[states].append((metadata["times"], calc))

    data_calc_dict = {}
    dtype = [("times", "f8"), ("calc", "f8")]
    for states, values in data_values.items():
        values_array = np.array(values, dtype=dtype)
        filler = np.zeros_like(values_array["calc"])

        data_calc = Data(exp=filler, err=filler, metadata=values_array["times"])
        data_calc.calc = values_array["calc"]
        data_calc.sort()

        data_calc_dict[states] = data_calc

    return data_calc_dict


class EXSYPlotter:
    def __init__(self, filename: Path, **_extra: Any) -> None:
        self.filename = filename
        self.printer: PlotPrinter = data_plot_printers["relaxation"]

    def plot(self, path: Path, profiles: list[Profile]) -> None:
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        name_exp = basename.with_suffix(".exp")
        name_fit = basename.with_suffix(".fit")

        print_plot_filename(name_pdf)

        with ExitStack() as stack:
            file_pdf = stack.enter_context(PdfPages(str(name_pdf)))
            file_calc = stack.enter_context(name_fit.open("w", encoding="utf-8"))
            file_exp = stack.enter_context(name_exp.open("w", encoding="utf-8"))
            for profile in sorted(profiles):
                data_exp_dict = create_plot_data_exp(profile)
                data_calc_dict = create_plot_data_calc(profile)
                for states in data_exp_dict:
                    data_exp = data_exp_dict[states]
                    data_calc = data_calc_dict[states]
                    name = f"{profile.spin_system!s}_{''.join(states)}"
                    plot_exsy(file_pdf, name, data_exp, data_calc)
                    file_exp.write(self.printer.print_exp(name, data_exp))
                    file_calc.write(self.printer.print_calc(name, data_calc))
