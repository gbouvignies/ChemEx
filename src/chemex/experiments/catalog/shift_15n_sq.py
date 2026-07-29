from __future__ import annotations

from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import ShiftDataSettings
from chemex.configuration.experiment import ExperimentSettings
from chemex.containers.data import Data
from chemex.experiments.experiment_types import (
    ProfileCalculation,
    register_experiment_type,
    shift_type,
)
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

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

    basis = Basis(type="ixy", spin_system="nh", model=config.model)
    return Spectrometer.from_spin_system(spin_system, basis, conditions)


def _find_nearest(array: Array, value: float) -> float:
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


class Shift15NSqSequence:
    """Sequence for 15N single-quantum chemical shift measurement."""

    def __init__(self, settings: Shift15NSqSettings) -> None:
        self.settings = settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        ppm_i = spectrometer.ppm_i
        ref_shift_i = spectrometer.par_values[self.settings.cs_i_name] * ppm_i
        shifts = spectrometer.analysis.calculate_shifts()
        shift_sq = _find_nearest(shifts, ref_shift_i)
        return np.array([shift_sq / ppm_i])

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def create_profile_calculation(
    config: Shift15NSqConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Shift15NSqSequence(config.experiment),
    )


EXPERIMENT_TYPE = shift_type(
    name=EXPERIMENT_NAME,
    config_type=Shift15NSqConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
