from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from chemex.configuration.data import ShiftDataSettings
from chemex.configuration.experiment import ExperimentConfig, ShiftSettings
from chemex.containers.dataset import load_shift_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import NoFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.plotters.shift import ShiftPlotter
from chemex.printers.data import ShiftPrinter

if TYPE_CHECKING:
    from chemex.containers.data import Data
    from chemex.parameters.spin_system import SpinSystem

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "shift_15n_sqmq"


class Shift15NSqMqSettings(ShiftSettings):
    name: Literal["shift_15n_sqmq"]
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def cs_i_name(self) -> str:
        return f"cs_i_{self.observed_state}"

    @property
    def cs_s_name(self) -> str:
        return f"cs_s_{self.observed_state}"


class Shift15NSqMqConfig(ExperimentConfig[Shift15NSqMqSettings, ShiftDataSettings]):
    ...


def build_spectrometer(
    config: Shift15NSqMqConfig, spin_system: SpinSystem
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="ixy_ixysxy", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    return Spectrometer(liouvillian)


def _find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


@dataclass
class Shift15NSqMqSequence:
    settings: Shift15NSqMqSettings

    def calculate(self, spectrometer: Spectrometer, _data: Data) -> np.ndarray:
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
                / spectrometer.liouvillian.ppm_i
            ]
        )

    def is_reference(self, metadata: NDArrayFloat) -> NDArrayBool:
        return np.full_like(metadata, False, dtype=np.bool_)


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
