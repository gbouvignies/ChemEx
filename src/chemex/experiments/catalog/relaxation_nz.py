from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import RelaxationSettings
from chemex.containers.data import Data
from chemex.experiments.experiment_types import (
    ProfileCalculation,
    register_experiment_type,
    relaxation_type,
)
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

EXPERIMENT_NAME = "relaxation_nz"


class RelaxationNzSettings(RelaxationSettings):
    """Settings for longitudinal 15N relaxation (Nz) experiment."""

    name: Literal["relaxation_nz"]

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Initial magnetization terms for the experiment.

        Returns:
            List of initial state terms for the Liouvillian calculation.

        """
        return self.get_start_terms("iz")

    @computed_field
    @property
    def detection(self) -> str:
        """Detection mode for the observable magnetization.

        Returns:
            Detection term for the Liouvillian calculation.

        """
        return self.get_detection_expression("[iz]")


class RelaxationNzConfig(
    ExperimentConfiguration[
        RelaxationNzSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
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

    basis = Basis(type="iz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.detection = settings.detection

    return spectrometer


class RelaxationNzSequence:
    """Sequence for longitudinal 15N relaxation (Nz) experiment."""

    def __init__(self, settings: RelaxationNzSettings) -> None:
        self.settings = settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        times = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        # Return profile
        delays = spectrometer.delays(times)
        return np.array([spectrometer.detect(delay @ start) for delay in delays])

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def create_profile_calculation(
    config: RelaxationNzConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=RelaxationNzSequence(config.experiment),
    )


EXPERIMENT_TYPE = relaxation_type(
    name=EXPERIMENT_NAME,
    config_type=RelaxationNzConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
