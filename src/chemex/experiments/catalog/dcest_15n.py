from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import DCestSettings
from chemex.configuration.types import ChemicalShift
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

Dataset = list[tuple[SpinSystem, Data]]

EXPERIMENT_NAME = "dcest_15n"

OFFSET_REF = 1e4


class DCest15NSettings(DCestSettings):
    """D-CEST 15N experiment settings with DANTE pulse train."""

    name: Literal["dcest_15n"]
    carrier: ChemicalShift = Field(description="15N carrier position during CEST (ppm)")

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms based on kinetic model."""
        if self.start_state is not None:
            return self.get_start_terms("iz")
        starts = {"2st_hd": ["iz_a"], "4st_hd": ["iz_a", "iz_b"]}
        return starts.get(self.model_name, ["iz"])

    @computed_field
    @property
    def detection(self) -> str:
        """Detection operator for the experiment."""
        return self.get_detection_expression("[iz]")


class DCest15NConfig(
    ExperimentConfiguration[
        DCest15NSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(config: DCest15NConfig, spin_system: SpinSystem) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier

    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.get_b1_model(),
    )

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        symbol = spin_system.symbols["i"]
        atom = spin_system.atoms["i"]
        spectrometer.jeff_i = get_multiplet(symbol, atom.name)

    return spectrometer


class DCest15NSequence:
    """Sequence for D-CEST 15N experiment."""

    def __init__(self, settings: DCest15NSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata
        spectrometer.offset_i = 0.0

        start = spectrometer.get_start_magnetization(
            self.settings.start_terms, Nucleus.N15
        )

        d_eq = (
            spectrometer.delays(self.settings.time_equil)
            if self.settings.time_equil > 0
            else spectrometer.identity
        )

        intensities: dict[float, Array] = {}

        for offset in set(offsets):
            if self.is_reference(offset):
                intensities[offset] = d_eq @ start
                continue

            spectrometer.offset_i = offset

            p_delay = spectrometer.delays(self.settings.tau_dante)
            p_pulse = spectrometer.pulse_i(self.settings.pw_dante, 0.0)

            intensities[offset] = (
                d_eq @ matrix_power(p_delay @ p_pulse, self.settings.ncyc_dante) @ start
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets],
        )


def create_profile_calculation(
    config: DCest15NConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=DCest15NSequence(config.experiment),
    )


EXPERIMENT_TYPE = cest_type(
    name=EXPERIMENT_NAME,
    config_type=DCest15NConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
