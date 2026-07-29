from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettingsNoRef
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

EXPERIMENT_NAME = "cest_ch3_1h_ip_ap"

OFFSET_REF = 1e4


class CestCh31HIpApSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for in-phase/anti-phase CH3 1H CEST experiment."""

    name: Literal["cest_ch3_1h_ip_ap"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    d1: float = Field(description="Relaxation delay in seconds")
    carrier: Frequency = Field(description="1H carrier position in Hz")
    taua: float = 2.00e-3

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting longitudinal identity components."""
        return self.get_start_terms("ie")

    @property
    def detection(self) -> str:
        return self.get_detection_expression("[2izsz]")


class CestCh31HIpApConfig(
    ExperimentConfiguration[
        CestCh31HIpApSettings,
        ConditionsWithValidations,
        CestDataSettingsNoRef,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=[
                "r2_i",
                f"r1_i_{state}",
                f"r1_s_{state}",
                f"etaxy_i_{state}",
                f"etaz_i_{state}",
            ],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: CestCh31HIpApConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz_eq", spin_system="hc", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.detection = settings.detection

    return spectrometer


class CestCh31HIpApSequence:
    """Sequence for in-phase/anti-phase CH3 1H CEST experiment."""

    def __init__(self, settings: CestCh31HIpApSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata
        spectrometer.offset_i = 0.0

        d_d1, d_taua = spectrometer.delays([self.settings.d1, self.settings.taua])

        pp90_i = spectrometer.perfect90_i
        pp180_isx = spectrometer.perfect180_i[0] @ spectrometer.perfect180_s[0]

        start = d_d1 @ spectrometer.get_start_magnetization(
            terms=self.settings.start_terms,
        )
        start = spectrometer.keep(start, components=["ie", "iz"])

        intensities: dict[float, float] = {}

        for offset in set(offsets):
            spectrometer.offset_i = offset
            mag = spectrometer.pulse_i(self.settings.time_t1, 0.0) @ start
            if self.is_reference(offset):
                inept = pp90_i[3] @ d_taua @ pp180_isx @ d_taua @ pp90_i[0]
                mag = inept @ mag
            intensities[offset] = spectrometer.detect(mag)

        return np.array([intensities[offset] for offset in offsets])


def create_profile_calculation(
    config: CestCh31HIpApConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=CestCh31HIpApSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=CestCh31HIpApConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
