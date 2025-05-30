from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Generic, Protocol, TypeVar

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename
from chemex.nmr.spectrometer import Spectrometer
from chemex.plotters.plot import get_grid, plot_profile
from chemex.printers.plot import PlotPrinter, data_plot_printers
from chemex.typing import ArrayFloat

_GREY400 = "#BDBDBD"
_LINESTYLES = [
    (0, ()),  # solid
    (0, (5, 1)),  # densely dashed
    (0, (5, 5)),  # dashed
    (0, (5, 10)),  # loosely dashed
    (0, (1, 1)),  # densely dotted
    (0, (1, 5)),  # dotted
    (0, (1, 10)),  # loosely dotted
]


class CestExperimentSettings(Protocol):
    sw: float


class CestExperimentConfig(Protocol):
    experiment: CestExperimentSettings


T = TypeVar("T", bound=CestExperimentConfig)


def estimate_sigma(values: ArrayFloat) -> float:
    """Estimates standard deviation using median to exclude outliers.

    Up to 50% can be bad.

    Reference:
    Rousseeuw, Peter & Croux, Christophe. (1993). Alternatives to Median Absolute
    Deviation. Journal of the American Statistical Association. 88. 1273 - 1283.
    10.1080/01621459.1993.10476408.

    """
    if not all(values):
        return 0.0
    _values = values.reshape(1, -1)
    return float(1.1926 * np.median(np.median(abs(_values - _values.T), axis=0)))


@dataclass
class CircularShift:
    spectrometer: Spectrometer
    sw: float
    sw_ppm: float = field(init=False)

    def __post_init__(self) -> None:
        sw_ppm, sw_ref = self.spectrometer.offsets_to_ppms(np.array([self.sw, 0.0]))
        self.sw_ppm = sw_ppm - sw_ref

    def centre(self, array: ArrayFloat, position: float) -> ArrayFloat:
        if self.sw == np.inf:
            return array
        cs_min = position - self.sw_ppm / 2.0
        return (array - cs_min) % self.sw_ppm + cs_min


def add_resonance_positions(
    ax1: Axes,
    ax2: Axes,
    cs_values: ArrayFloat,
    centre: float,
    circular_shift: CircularShift,
) -> tuple[Axes, Axes]:
    kwargs2 = {"color": _GREY400, "linewidth": 0.75, "zorder": -1}
    cs_shifted = circular_shift.centre(cs_values, centre)
    for a_cs, a_cs_shifted, lstyle in zip(
        cs_values,
        cs_shifted,
        _LINESTYLES,
        strict=False,
    ):
        ax1.axvline(a_cs_shifted, linestyle=lstyle, **kwargs2)
        ax2.axvline(a_cs_shifted, linestyle=lstyle, **kwargs2)
        if a_cs != a_cs_shifted:
            x, _ = ax2.transLimits.transform((a_cs_shifted, 0))
            ax2.text(x - 0.02, 0.95, "*", transform=ax2.transAxes)
    return ax1, ax2


def plot_dcest(
    file_pdf: PdfPages,
    name: str,
    data_exp: Data,
    data_calc: Data,
    cs_values: ArrayFloat,
    circular_shift: CircularShift,
) -> None:
    if data_exp.size > 0:
        residuals = data_exp.exp - data_exp.calc
        sigma = estimate_sigma(residuals)
    else:
        sigma = 0.0
    centre = float(np.mean(cs_values))
    data_exp.metadata = circular_shift.centre(data_exp.metadata, centre)
    data_calc.metadata = circular_shift.centre(data_calc.metadata, centre)
    data_exp.sort()
    data_calc.sort()
    fig = plot_profile(name, data_exp, data_calc)
    ax1, ax2 = fig.axes
    xmin, xmax = sorted(ax2.get_xlim())
    ax1.set_xlim(xmax, xmin)
    ax2.set_xlim(xmax, xmin)
    ax1.xaxis.set_major_locator(MaxNLocator(6))
    ax2.xaxis.set_major_locator(MaxNLocator(6))
    ax2.set_xlabel(r"$B_1$ position (ppm)")
    ax2.set_ylabel(r"$I/I_0$")
    kwargs1 = {"facecolor": (0, 0, 0, 0.1), "edgecolor": "none"}
    ax1.fill_between(ax1.get_xlim(), -1.0 * sigma, 1.0 * sigma, **kwargs1)
    ax1.fill_between(ax1.get_xlim(), -2.0 * sigma, 2.0 * sigma, **kwargs1)

    ax1, ax2 = add_resonance_positions(ax1, ax2, cs_values, centre, circular_shift)

    if data_exp.size == 0:
        fig.delaxes(ax1)

    file_pdf.savefig(fig)


def create_plot_data_exp(profile: Profile) -> Data:
    spectrometer = profile.spectrometer
    data = profile.data
    refs = data.refs

    ppms = spectrometer.offsets_to_ppms(data.metadata[~refs])

    intensity0 = np.mean(data.exp[refs])
    intensities = data.exp[~refs] / intensity0
    intensities_calc = data.calc[~refs] / intensity0

    errors = data.err[~refs] / abs(intensity0)
    errorbars = np.array([errors, errors]).transpose()

    data_exp = Data(exp=intensities, err=errorbars, metadata=ppms)
    data_exp.calc = intensities_calc
    data_exp.mask = data.mask[~refs]
    data_exp.sort()

    return data_exp


def create_plot_data_calc(profile: Profile) -> Data:
    spectrometer = profile.spectrometer
    data = profile.data
    refs = data.refs

    offsets_plot = get_grid(data.metadata[~refs], 400, 0.02)

    filler = np.full_like(offsets_plot, 0.0)

    data_for_calculation = Data(exp=filler, err=filler, metadata=offsets_plot)
    data_for_calculation.calc = data.scale * profile.pulse_sequence.calculate(
        spectrometer,
        data_for_calculation,
    )

    ppms = spectrometer.offsets_to_ppms(offsets_plot)
    intensity0 = np.mean(data.exp[refs])
    data_fit = Data(exp=filler, err=filler, metadata=ppms)
    data_fit.calc = data_for_calculation.calc / intensity0

    return data_fit


def get_state_positions(spectrometer: Spectrometer) -> ArrayFloat:
    names = (f"cs_i_{state}" for state in "abcdef")
    return np.array(
        [
            spectrometer.par_values[name]
            for name in names
            if name in spectrometer.par_values
        ],
    )


class CestPlotter(Generic[T]):
    def __init__(self, filename: Path, config: T, **_extra: Any) -> None:
        self.filename = filename
        self.config = config
        self.printer: PlotPrinter = data_plot_printers["cest"]

    def _plot_profile(
        self,
        pdf: PdfPages,
        profile: Profile,
        data_exp: Data,
        data_calc: Data,
    ) -> None:
        spectrometer = profile.spectrometer
        sw = self.config.experiment.sw
        dip_positions = get_state_positions(spectrometer)
        circular_shift = CircularShift(spectrometer, sw)
        plot_dcest(
            pdf,
            str(profile.spin_system),
            data_exp,
            data_calc,
            dip_positions,
            circular_shift,
        )

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
                data_exp = create_plot_data_exp(profile)
                data_calc = create_plot_data_calc(profile)
                self._plot_profile(file_pdf, profile, data_exp, data_calc)
                file_exp.write(
                    self.printer.print_exp(str(profile.spin_system), data_exp),
                )
                file_calc.write(
                    self.printer.print_calc(str(profile.spin_system), data_calc),
                )

    def plot_simulation(self, path: Path, profiles: list[Profile]) -> None:
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        name_sim = basename.with_suffix(".sim")

        print_plot_filename(name_pdf, extra=False)

        with ExitStack() as stack:
            file_pdf = stack.enter_context(PdfPages(str(name_pdf)))
            file_sim = stack.enter_context(name_sim.open("w", encoding="utf-8"))
            for profile in sorted(profiles):
                data_exp = Data(np.array([]), np.array([]), np.array([]))
                data_calc = create_plot_data_calc(profile)
                self._plot_profile(file_pdf, profile, data_exp, data_calc)
                file_sim.write(
                    self.printer.print_calc(str(profile.spin_system), data_calc),
                )
