from __future__ import annotations

from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration
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
): ...


def build_spectrometer(
    config: Shift15NSqMqConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    conditions = config.conditions

    basis = Basis(type="ixy_ixysxy", spin_system="nh", model=config.model)
    return Spectrometer.from_spin_system(spin_system, basis, conditions)


class Shift15NSqMqSequence:
    """Sequence for 15N single-quantum and multiple-quantum chemical shift measurement."""

    def __init__(self, settings: Shift15NSqMqSettings) -> None:
        self.settings = settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        del data
        observed_state = self.settings.observed_state
        shift_sq = spectrometer.analysis.calculate_shift(
            observed_state,
            coherence="sq_i",
            relative=True,
        )
        shift_dq = spectrometer.analysis.calculate_shift(
            observed_state,
            coherence="dq_is",
            relative=True,
        )
        shift_zq = spectrometer.analysis.calculate_shift(
            observed_state,
            coherence="zq_is",
            relative=True,
        )
        return np.array(
            [
                1e3 * (shift_sq - 0.5 * (shift_dq + shift_zq)) / spectrometer.ppm_i,
            ],
        )

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def create_profile_calculation(
    config: Shift15NSqMqConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Shift15NSqMqSequence(config.experiment),
    )


EXPERIMENT_TYPE = shift_type(
    name=EXPERIMENT_NAME,
    config_type=Shift15NSqMqConfig,
    unit="PPB",
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
