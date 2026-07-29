from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, CestSettings
from chemex.configuration.types import B1Field, Delay, Frequency
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
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import Array

EXPERIMENT_NAME = "cest_15n_cw"

OFFSET_REF = 1e4


class Cest15NCwSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for constant-wave (CW) 15N CEST experiment."""

    name: Literal["cest_15n_cw"]
    time_t1: Delay = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="15N carrier position in Hz")
    carrier_dec: Frequency = Field(description="1H decoupling carrier position in Hz")
    b1_frq_dec: B1Field = Field(description="1H decoupling B1 field strength in Hz")

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms for the experiment."""
        return self.get_start_terms("iz")

    @computed_field
    @property
    def detection(self) -> str:
        """Detection operator for the experiment."""
        return self.get_detection_expression("[iz]")


class Cest15NCwConfig(
    ExperimentConfiguration[
        Cest15NCwSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=[
                "r2_i",
                f"r2_s_{state}",
                f"r1_i_{state}",
                f"r1_s_{state}",
                f"r2mq_is_{state}",
                f"etaxy_i_{state}",
                f"etaz_i_{state}",
            ],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: Cest15NCwConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsxyz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.b1_s = settings.b1_frq_dec
    spectrometer.carrier_s = settings.carrier_dec
    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        spectrometer.jeff_i = get_multiplet("", "n")

    return spectrometer


class Cest15NCwSequence:
    """Sequence for constant-wave 15N CEST experiment."""

    def __init__(self, settings: Cest15NCwSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(
            terms=self.settings.start_terms, atom=Nucleus.N15
        )

        intensities: dict[float, Array] = {}

        for offset in set(offsets):
            intensities[offset] = start

            if self.is_reference(offset):
                continue

            spectrometer.offset_i = offset

            intensities[offset] = (
                spectrometer.pulse_is(self.settings.time_t1, 0.0, 0.0)
                @ intensities[offset]
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets],
        )


def create_profile_calculation(
    config: Cest15NCwConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Cest15NCwSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=Cest15NCwConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
