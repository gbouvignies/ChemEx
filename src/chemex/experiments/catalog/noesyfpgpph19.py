from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import ExperimentSettings
from chemex.containers.data import Data
from chemex.containers.dataset import Dataset, load_exsy_dataset
from chemex.experiments.experiment_types import (
    ExperimentSupport,
    ProfileCalculation,
    experiment_type,
    register_experiment_type,
)
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.exsy import EXSYPlotter
from chemex.printers.data import EXSYPrinter
from chemex.typing import Array

EXPERIMENT_NAME = "noesyfpgpph19"


class Noesyfpgpph19Settings(ExperimentSettings):
    name: Literal["noesyfpgpph19"]


class Noesyfpgpph19Config(
    ExperimentConfiguration[
        Noesyfpgpph19Settings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        return ToBeFitted(rates=["r1_i"], model_free=["tauc", "s2"])


def build_spectrometer(
    config: Noesyfpgpph19Config,
    spin_system: SpinSystem,
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="iz", spin_system="hn", model=config.model)
    return Spectrometer.from_spin_system(spin_system, basis, conditions)


class Noesyfpgpph19Sequence:
    """Sequence for NOESY-FPGP-PH19 experiment."""

    def __init__(self, settings: Noesyfpgpph19Settings) -> None:
        self.settings = settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        times = data.metadata["times"]
        states1 = data.metadata["states1"]
        states2 = data.metadata["states2"]

        # Getting the starting magnetization
        equilibrium = spectrometer.get_equilibrium()

        # Calculate delay propagators
        delays = spectrometer.delays(times)

        intensities: list[float] = []
        for state1, state2, delay in zip(states1, states2, delays, strict=False):
            start = spectrometer.keep(equilibrium, [f"iz_{state1}"])
            spectrometer.detection = f"[iz_{state2}]"
            intensities.append(spectrometer.detect(delay @ start))

        return np.array(intensities)

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def create_profile_calculation(
    config: Noesyfpgpph19Config,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Noesyfpgpph19Sequence(config.experiment),
    )


def load_dataset(base_path: Path, settings: Noesyfpgpph19Config) -> Dataset:
    return load_exsy_dataset(base_path, settings)


def create_filterer(
    *,
    config: Noesyfpgpph19Config,
    spectrometer: Spectrometer,
) -> PlanesFilterer:
    return PlanesFilterer(config=config, spectrometer=spectrometer)


def create_printer() -> EXSYPrinter:
    return EXSYPrinter()


def create_plotter(
    *,
    filename: Path,
    config: Noesyfpgpph19Config,
) -> EXSYPlotter:
    return EXSYPlotter(filename=filename, config=config)


EXSY_SUPPORT = ExperimentSupport[Noesyfpgpph19Config](
    load_dataset=load_dataset,
    create_filterer=create_filterer,
    create_printer=create_printer,
    create_plotter=create_plotter,
)

EXPERIMENT_TYPE = experiment_type(
    name=EXPERIMENT_NAME,
    config_type=Noesyfpgpph19Config,
    support=EXSY_SUPPORT,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
