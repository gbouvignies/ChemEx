from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import (
    ExperimentConfig,
    RelaxationSettings,
    ToBeFitted,
)
from chemex.containers.dataset import load_exsy_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.plotters.exsy import EXSYPlotter
from chemex.printers.data import EXSYPrinter

if TYPE_CHECKING:
    from chemex.containers.data import Data
    from chemex.parameters.spin_system import SpinSystem

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "noesyfpgpph19"


class Noesyfpgpph19Settings(RelaxationSettings):
    name: Literal["noesyfpgpph19"]


class Noesyfpgpph19Config(
    ExperimentConfig[Noesyfpgpph19Settings, RelaxationDataSettings],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        return ToBeFitted(rates=["r1_i"], model_free=["tauc", "s2"])


def build_spectrometer(
    config: Noesyfpgpph19Config,
    spin_system: SpinSystem,
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="iz", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)

    return Spectrometer(liouvillian)


@dataclass
class Noesyfpgpph19Sequence:
    settings: Noesyfpgpph19Settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:
        times = data.metadata["times"]
        states1 = data.metadata["states1"]
        states2 = data.metadata["states2"]

        # Getting the starting magnetization
        equilibrium = spectrometer.get_equilibrium()

        # Calculate delay propagators
        delays = spectrometer.delays(times)

        intensities = []
        for state1, state2, delay in zip(states1, states2, delays):
            start = spectrometer.keep(equilibrium, [f"iz_{state1}"])
            spectrometer.detection = f"[iz_{state2}]"
            intensities.append(spectrometer.detect(delay @ start))

        return np.array(intensities)

    def is_reference(self, metadata: NDArrayFloat) -> NDArrayBool:
        return np.full_like(metadata, False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=Noesyfpgpph19Config,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Noesyfpgpph19Sequence,
        dataset_creator=load_exsy_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=EXSYPrinter,
        plotter_creator=EXSYPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
