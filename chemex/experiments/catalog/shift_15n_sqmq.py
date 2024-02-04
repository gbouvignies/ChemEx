from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration
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

EXPERIMENT_NAME = "shift_15n_sqmq"


class Shift15NSqMqSettings(ExperimentSettings):
    name: Literal["shift_15n_sqmq"]

    @property
    def cs_i_name(self) -> str:
        return f"cs_i_{self.observed_state}"

    @property
    def cs_s_name(self) -> str:
        return f"cs_s_{self.observed_state}"


class Shift15NSqMqConfig(
    ExperimentConfiguration[
        Shift15NSqMqSettings,
        ConditionsWithValidations,
        ShiftDataSettings,
    ],
):
    ...


def build_spectrometer(
    config: Shift15NSqMqConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="ixy_ixysxy", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    return Spectrometer(liouvillian)


def _find_nearest(array: ArrayFloat, value: float) -> float:
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


@dataclass
class Shift15NSqMqSequence:
    settings: Shift15NSqMqSettings

    def calculate(self, spectrometer: Spectrometer, _data: Data) -> ArrayFloat:
        ref_shift_i = (
            spectrometer.par_values[self.settings.cs_i_name]
            * spectrometer.liouvillian.ppm_i
        )
        ref_shift_s = (
            spectrometer.par_values[self.settings.cs_s_name]
            * spectrometer.liouvillian.ppm_s
        )
        ref_shift_dq = ref_shift_i + ref_shift_s
        ref_shift_zq = ref_shift_i - ref_shift_s
        shifts = spectrometer.calculate_shifts()
        shift_sq = _find_nearest(shifts, ref_shift_i)
        shift_dq = _find_nearest(shifts, ref_shift_dq)
        shift_zq = _find_nearest(shifts, ref_shift_zq)
        return np.array(
            [
                1e3
                * (shift_sq - 0.5 * (shift_dq + shift_zq))
                / spectrometer.liouvillian.ppm_i,
            ],
        )

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=Shift15NSqMqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Shift15NSqMqSequence,
        dataset_creator=load_shift_dataset,
        filterer_creator=NoFilterer,
        printer_creator=ShiftPrinter,
        plotter_creator=ShiftPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
