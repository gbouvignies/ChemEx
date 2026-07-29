from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, CestSettings
from chemex.configuration.types import ChemicalShift, Delay
from chemex.containers.data import Data
from chemex.experiments.experiment_types import (
    ProfileCalculation,
    cest_type,
    register_experiment_type,
)
from chemex.nmr.basis import Basis
from chemex.nmr.constants import get_multiplet
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

EXPERIMENT_NAME = "cest_13c"
OFFSET_REF = 1e4


class Cest13CSettings(CestSettings, B1InhomogeneityMixin):
    """Pure in-phase 13C CEST experiment settings."""

    name: Literal["cest_13c"]
    time_t1: Delay = Field(description="CEST relaxation delay (seconds)")
    carrier: ChemicalShift = Field(description="13C carrier position during CEST (ppm)")

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms."""
        return self.get_start_terms("iz")

    @computed_field
    @property
    def detection(self) -> str:
        """Detection operator."""
        return self.get_detection_expression("[iz]")


class Cest13CConfig(
    ExperimentConfiguration[
        Cest13CSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(config: Cest13CConfig, spin_system: SpinSystem) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="ch", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        symbol = spin_system.symbols["i"]
        atom = spin_system.atoms["i"]
        spectrometer.jeff_i = get_multiplet(symbol, atom.name)

    return spectrometer


class Cest13CSequence:
    """Sequence for CEST 13C experiment."""

    def __init__(self, settings: Cest13CSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        intensities: dict[float, Array] = {}

        for offset in set(offsets):
            intensities[offset] = start

            if self.is_reference(offset):
                continue

            spectrometer.offset_i = offset

            intensities[offset] = (
                spectrometer.pulse_i(self.settings.time_t1, 0.0) @ intensities[offset]
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets],
        )


def create_profile_calculation(
    config: Cest13CConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Cest13CSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=Cest13CConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
