from __future__ import annotations

from typing import ClassVar, Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, CestSettings
from chemex.configuration.types import Frequency
from chemex.containers.data import Data
from chemex.experiments.experiment_types import (
    ProfileCalculation,
    cest_type,
    register_experiment_type,
)
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

EXPERIMENT_NAME = "cest_1hn_ap"

OFFSET_REF = 1e4


class Cest1HnApSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for anti-phase 1H-15N CEST experiment."""

    name: Literal["cest_1hn_ap"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="1H carrier position in Hz")
    start_from_observed_by_default: ClassVar[bool] = True

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization term (anti-phase)."""
        return self.get_start_terms("2izsz")

    @computed_field
    @property
    def detection(self) -> str:
        """Detection operator (anti-phase)."""
        return self.get_detection_expression("[2izsz]")


class Cest1HnApConfig(
    ExperimentConfiguration[
        Cest1HnApSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}", f"r1_s_{state}", f"etaxy_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: Cest1HnApConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="hn", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.detection = settings.detection

    return spectrometer


class Cest1HnApSequence:
    """Sequence for anti-phase 1H-15N CEST experiment."""

    def __init__(self, settings: Cest1HnApSettings) -> None:
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
            [spectrometer.detect(intensities[offset]) for offset in offsets]
        )


def create_profile_calculation(
    config: Cest1HnApConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Cest1HnApSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=Cest1HnApConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
