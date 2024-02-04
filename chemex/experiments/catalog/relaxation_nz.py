from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import ExperimentSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.relaxation import RelaxationPlotter
from chemex.printers.data import RelaxationPrinter
from chemex.typing import ArrayBool, ArrayFloat

EXPERIMENT_NAME = "relaxation_nz"


class RelaxationNzSettings(ExperimentSettings):
    name: Literal["relaxation_nz"]

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class RelaxationNzConfig(
    ExperimentConfiguration[
        RelaxationNzSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: RelaxationNzConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="iz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class RelaxationNzSequence:
    settings: RelaxationNzSettings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        times = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_equilibrium()

        # Return profile
        delays = spectrometer.delays(times)
        return np.array([spectrometer.detect(delay @ start) for delay in delays])

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=RelaxationNzConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=RelaxationNzSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=RelaxationPrinter,
        plotter_creator=RelaxationPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
