from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import ShiftDataSettings
from chemex.configuration.experiment import ExperimentSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_shift_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import NoFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.shift import ShiftPlotter
from chemex.printers.data import ShiftPrinter
from chemex.typing import ArrayBool, ArrayFloat

EXPERIMENT_NAME = "shift_15n_sq"


class Shift15NSqSettings(ExperimentSettings):
    name: Literal["shift_15n_sq"]

    @property
    def cs_i_name(self) -> str:
        return f"cs_i_{self.observed_state}"


class Shift15NSqConfig(
    ExperimentConfiguration[
        Shift15NSqSettings, ConditionsWithValidations, ShiftDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"cs_i_{state}"],
            model_free=[f"cs_i_{state}"],
        )


def build_spectrometer(
    config: Shift15NSqConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="ixy", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    return Spectrometer(liouvillian)


def _find_nearest(array: ArrayFloat, value: float) -> float:
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


@dataclass
class Shift15NSqSequence:
    settings: Shift15NSqSettings

    def calculate(self, spectrometer: Spectrometer, _data: Data) -> ArrayFloat:
        ppm_i = spectrometer.liouvillian.ppm_i
        ref_shift_i = spectrometer.par_values[self.settings.cs_i_name] * ppm_i
        shifts = spectrometer.calculate_shifts()
        shift_sq = _find_nearest(shifts, ref_shift_i)
        return np.array([shift_sq / ppm_i])

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=Shift15NSqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Shift15NSqSequence,
        dataset_creator=load_shift_dataset,
        filterer_creator=NoFilterer,
        printer_creator=ShiftPrinter,
        plotter_creator=ShiftPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
