from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import ExperimentConfig
from chemex.configuration.experiment import RelaxationSettings
from chemex.configuration.experiment import ToBeFitted
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators
from chemex.experiments.factories import factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.liouvillian import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters import RelaxationPlotter
from chemex.printers.data import RelaxationPrinter

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "relaxation_nz"


class RelaxationNzSettings(RelaxationSettings):
    name: Literal["relaxation_nz"]
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class RelaxationNzConfig(
    ExperimentConfig[RelaxationNzSettings, RelaxationDataSettings]
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r1_i_{state}"], model_free=[f"tauc_{state}", f"s2_{state}"]
        )


def build_spectrometer(
    config: RelaxationNzConfig, spin_system: SpinSystem
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

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        times = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_equilibrium()

        # Return profile
        delays = spectrometer.delays(times)
        return np.array([spectrometer.detect(delay @ start) for delay in delays])

    def is_reference(self, metadata: NDArrayFloat) -> NDArrayBool:
        return np.full_like(metadata, False, dtype=np.bool_)


def register():
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
