from __future__ import annotations

from typing import Literal

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
from chemex.nmr.constants import get_multiplet
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

EXPERIMENT_NAME = "cest_15n_tr"

OFFSET_REF = 1e4


class Cest15NTrSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for TROSY-based 15N CEST experiment."""

    name: Literal["cest_15n_tr"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="15N carrier position in Hz")
    antitrosy: bool = False

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Start from the TROSY or ANTI-TROSY component.

        TROSY: (2IzSz + Iz) / 2
        ANTITROSY: (2IzSz - Iz) / 2.
        """
        if not self.antitrosy:
            return self.get_start_terms("2izsz", "-iz")
        return self.get_start_terms("2izsz", "iz")

    @computed_field
    @property
    def detection(self) -> str:
        """Detection operator for TROSY or ANTI-TROSY component."""
        if self.antitrosy:
            return self.get_detection_expression("[2izsz] + [iz]")
        return self.get_detection_expression("[2izsz] - [iz]")


class Cest15NTrConfig(
    ExperimentConfiguration[
        Cest15NTrSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state

        to_be_fitted = ToBeFitted(
            rates=["r2_i", f"r1_i_{state}", "r1a_is"],
            model_free=[f"tauc_{state}", f"s2_{state}", "khh"],
        )

        if self.experiment.antitrosy:
            to_be_fitted.rates.extend([f"etaxy_i_{state}", f"etaz_i_{state}"])

        return to_be_fitted


def build_spectrometer(
    config: Cest15NTrConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        spectrometer.jeff_i = get_multiplet("", "n")

    return spectrometer


class Cest15NTrSequence:
    """Sequence for TROSY-based 15N CEST experiment."""

    def __init__(self, settings: Cest15NTrSettings) -> None:
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
    config: Cest15NTrConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Cest15NTrSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=Cest15NTrConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
