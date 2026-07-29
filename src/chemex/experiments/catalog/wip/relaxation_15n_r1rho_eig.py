from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, DetectionSettings
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

EXPERIMENT_NAME = "wip.relaxation_15n_r1rho_eig"


class Relaxation15NR1RhoSettings(DetectionSettings, B1InhomogeneityMixin):
    name: Literal["wip.relaxation_15n_r1rho_eig"]

    carrier: float
    b1_frq: float
    legacy_b1_inh_scale_default: ClassVar[float | None] = np.inf
    legacy_b1_inh_res_default: ClassVar[int | None] = 11

    @property
    def detection(self) -> str:
        return self.get_detection_expression("[iz]")


class Relaxation15NR1RhoConfig(
    ExperimentConfiguration[
        Relaxation15NR1RhoSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state
        return ToBeFitted(
            rates=[f"r2_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: Relaxation15NR1RhoConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier
    spectrometer.set_b1_i_inhomogeneity(
        settings.get_b1_nominal(),
        settings.b1_distribution,
    )

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Relaxation15NR1RhoSequence:
    settings: Relaxation15NR1RhoSettings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        times = data.metadata

        r1rho = spectrometer.analysis.calculate_r1rho()

        # Return profile
        return np.exp(-r1rho * times)

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def create_profile_calculation(
    config: Relaxation15NR1RhoConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Relaxation15NR1RhoSequence(config.experiment),
    )


EXPERIMENT_TYPE = relaxation_type(
    name=EXPERIMENT_NAME,
    config_type=Relaxation15NR1RhoConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
