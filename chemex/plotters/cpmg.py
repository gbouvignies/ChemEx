from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path
from typing import Any, Generic, Protocol, TypeVar

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from chemex.containers.data import Data
from chemex.containers.profile import Profile
from chemex.messages import print_plot_filename
from chemex.plotters.plot import plot_profile
from chemex.printers.plot import PlotPrinter, data_plot_printers
from chemex.typing import ArrayFloat

LARGE_ERROR = 1e15


rng = np.random.default_rng()


class CpmgExperimentSettings(Protocol):
    time_t2: float
    even_ncycs: bool = False


class CpmgExperimentConfig(Protocol):
    experiment: CpmgExperimentSettings


T = TypeVar("T", bound=CpmgExperimentConfig)


def plot_cpmg(file_pdf: PdfPages, name: str, data_exp: Data, data_calc: Data):
    fig = plot_profile(name, data_exp, data_calc)
    ax2 = fig.axes[1]
    ax2.set_xlabel(r"$\nu_\mathregular{CPMG}$ (Hz)")
    ax2.set_ylabel(r"$R_{2,\mathregular{eff}}$ (s$^{-1}$)")
    ax2.set_xlim(0.0, max(data_calc.metadata) + min(data_calc.metadata))
    if data_exp.size == 0:
        fig.delaxes(fig.axes[0])
    file_pdf.savefig(fig)


def intensities_to_rates(
    intensities: ArrayFloat,
    intensities0: ArrayFloat,
    time_t2: float,
) -> ArrayFloat:
    normalized_intensities = intensities / np.mean(intensities0, axis=-1, keepdims=True)

    # If the normalized intensity is negative, no rate can be estimated.
    # The rate is then set to infinity by default.
    rates = np.full_like(intensities, np.inf)
    neg = normalized_intensities <= 0
    rates[~neg] = -np.log(normalized_intensities[~neg]) / time_t2

    return rates


def calculate_exp_rates(data: Data, time_t2: float) -> ArrayFloat:
    intensities = data.exp[~data.refs]
    intensities0 = data.exp[data.refs]
    return intensities_to_rates(intensities, intensities0, time_t2)


def calculate_calc_rates(data: Data, time_t2: float) -> ArrayFloat:
    intensities = data.calc[~data.refs]
    intensities0 = data.exp[data.refs]
    return intensities_to_rates(intensities, intensities0, time_t2)


def calculate_errorbars(data: Data, rates: ArrayFloat, time_t2: float) -> ArrayFloat:
    randn = rng.standard_normal(size=(10000, 1))
    randn0 = rng.standard_normal(size=(10000, 1))

    intensities = data.exp[~data.refs]
    intensities0 = data.exp[data.refs]

    intensities_err = data.err[~data.refs]
    intensities0_err = data.err[data.refs]

    intensities_ensemble = intensities + intensities_err * randn
    intensities0_ensemble = intensities0 + intensities0_err * randn0

    rates_ensemble: ArrayFloat = intensities_to_rates(
        intensities_ensemble,
        intensities0_ensemble,
        time_t2,
    )

    # Set infinity rates to high values as `np.percentile` does not take infinity
    rates_ensemble[rates_ensemble == np.inf] = 1e16

    errors = np.percentile(rates_ensemble - rates, [15.9, 84.1], axis=0).transpose()

    # Set back large errors to infinity
    errors[errors > LARGE_ERROR] = np.inf

    return np.abs(errors)


def ncycs_to_nu_cpmgs(ncycs: ArrayFloat, time_t2: float) -> ArrayFloat:
    modified_ncycs = ncycs.copy()
    modified_ncycs[modified_ncycs == -1] = 0.5
    return modified_ncycs[modified_ncycs != 0] / time_t2


def create_plot_data_exp(profile: Profile, config: CpmgExperimentConfig) -> Data:
    time_t2 = config.experiment.time_t2
    data = profile.data
    refs = data.refs

    nu_cpmgs = ncycs_to_nu_cpmgs(data.metadata, time_t2)
    rates = calculate_exp_rates(data, time_t2)
    rates_calc = calculate_calc_rates(data, time_t2)
    errorbars = calculate_errorbars(data, rates, time_t2)

    data_exp = Data(exp=rates, err=errorbars, metadata=nu_cpmgs)
    data_exp.calc = rates_calc
    data_exp.mask = data.mask[~refs]
    data_exp.sort()

    return data_exp


def create_plot_data_calc(profile: Profile, config: CpmgExperimentConfig) -> Data:
    spectrometer = profile.spectrometer
    time_t2 = config.experiment.time_t2
    data = profile.data
    refs = data.refs

    step = 2 if config.experiment.even_ncycs else 1
    ncycs = np.arange(2, max(data.metadata) + 1, step, dtype=np.float64)
    ncycs = np.asarray(sorted(set(ncycs) | set(data.metadata[~refs])))

    filler = np.zeros_like(ncycs)
    data_for_calculation = Data(exp=filler, err=filler, metadata=ncycs)

    intensities = profile.data.scale * profile.pulse_sequence.calculate(
        spectrometer,
        data_for_calculation,
    )

    nu_cpmgs = ncycs_to_nu_cpmgs(ncycs, time_t2)
    data_fit = Data(exp=filler, err=filler, metadata=nu_cpmgs)
    data_fit.calc = intensities_to_rates(intensities, data.exp[refs], time_t2)

    return data_fit


class CpmgPlotter(Generic[T]):
    def __init__(self, filename: Path, config: T, **_extra: Any) -> None:
        self.filename = filename
        self.config = config
        self.printer: PlotPrinter = data_plot_printers["cpmg"]

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
                data_exp = create_plot_data_exp(profile, self.config)
                data_calc = create_plot_data_calc(profile, self.config)
                plot_cpmg(file_pdf, str(profile.spin_system), data_exp, data_calc)
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
                data_calc = create_plot_data_calc(profile, self.config)
                plot_cpmg(file_pdf, str(profile.spin_system), data_exp, data_calc)
                file_sim.write(
                    self.printer.print_calc(str(profile.spin_system), data_calc),
                )
