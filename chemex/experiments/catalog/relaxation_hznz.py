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

EXPERIMENT_NAME = "relaxation_hznz"


class RelaxationHzNzSettings(ExperimentSettings):
    name: Literal["relaxation_hznz"]

    @property
    def detection(self) -> str:
        return f"[2izsz_{self.observed_state}]"


class RelaxationHzNzConfig(
    ExperimentConfiguration[
        RelaxationHzNzSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r1a_is_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: RelaxationHzNzConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="izsz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class RelaxationHzNzSequence:
    settings: RelaxationHzNzSettings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        times = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(["2izsz"])

        # Return profile
        delays = spectrometer.delays(0.25 * np.array(times))
        p180_i = spectrometer.perfect180_i[0]
        p180_s = spectrometer.perfect180_s[0]
        return np.array(
            [
                spectrometer.detect(
                    delay @ p180_s @ delay @ p180_i @ delay @ p180_s @ delay @ start,
                )
                for delay in delays
            ],
        )

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=RelaxationHzNzConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=RelaxationHzNzSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=RelaxationPrinter,
        plotter_creator=RelaxationPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
