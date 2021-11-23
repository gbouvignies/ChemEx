from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename
from chemex.plotters.plot import get_grid
from chemex.plotters.plot import plot_profile
from chemex.printers.plot import data_plot_printers
from chemex.printers.plot import PlotPrinter


def plot_relaxation(file_pdf: PdfPages, name: str, data_exp: Data, data_calc: Data):
    fig = plot_profile(name, data_exp, data_calc)
    ax2 = fig.axes[1]
    ax2.set_xlabel(r"Time (s)")
    ax2.set_ylabel(r"Intensity")
    ax2.set_xlim(0.0, max(data_exp.metadata) + min(data_calc.metadata))
    file_pdf.savefig(fig)


def create_plot_data_exp(profile: Profile) -> Data:
    data = profile.data

    times = data.metadata

    intensity0 = data.exp[np.argmax(abs(data.exp))]
    intensities = data.exp / intensity0
    intensities_calc = data.calc / intensity0

    errors = data.err / abs(intensity0)
    errorbars = np.array([-errors, errors]).transpose()

    data_exp = Data(exp=intensities, err=errorbars, metadata=times)
    data_exp.calc = intensities_calc
    data_exp.mask = data.mask
    data_exp.sort()

    return data_exp


def create_plot_data_calc(profile: Profile) -> Data:
    spectrometer = profile.spectrometer
    data = profile.data

    times = get_grid(data.metadata, 100, 0.02)

    filler = np.zeros_like(times)
    data_fit = Data(exp=filler, err=filler, metadata=times)

    intensity0 = data.exp[np.argmax(abs(data.exp))]
    scale = profile.data.scale / intensity0
    data_fit.calc = scale * profile.pulse_sequence.calculate(spectrometer, data_fit)

    return data_fit


class RelaxationPlotter:
    def __init__(self, filename: Path, **_extra: Any):
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
            file_calc = stack.enter_context(name_fit.open("w"))
            file_exp = stack.enter_context(name_exp.open("w"))
            for profile in sorted(profiles):
                data_exp = create_plot_data_exp(profile)
                data_calc = create_plot_data_calc(profile)
                plot_relaxation(file_pdf, str(profile.name), data_exp, data_calc)
                file_exp.write(self.printer.print_exp(str(profile.name), data_exp))
                file_calc.write(self.printer.print_calc(str(profile.name), data_calc))
