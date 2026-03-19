from __future__ import annotations

# pyright: reportGeneralTypeIssues=false
from dataclasses import dataclass
from functools import cached_property
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import CestSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import CestFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.constants import get_multiplet
from chemex.nmr.distributions.gaussian import GaussianDistributionConfig
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.printers.data import CestPrinter
from chemex.typing import Array

EXPERIMENT_NAME = "cest_15n_test"

OFFSET_REF = 1e4


class Cest15NSettings(CestSettings):
    name: Literal["cest_15n_test"]
    time_t1: float
    carrier: float
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11

    @cached_property
    def start_terms(self) -> list[str]:
        return [f"iz{self.suffix}"]

    @cached_property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cest15NTestConfig(
    ExperimentConfiguration[
        Cest15NSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: Cest15NTestConfig, spin_system: SpinSystem
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh", model=config.model)
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.set_b1_i_inhomogeneity(
        settings.b1_frq,
        GaussianDistributionConfig(
            scale=settings.b1_inh_scale,
            res=settings.b1_inh_res,
        ),
    )

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        liouvillian.jeff_i = get_multiplet("", "n")

    return spectrometer


@dataclass
class Cest15NTestSequence:
    settings: Cest15NSettings

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
            [spectrometer.detect_ls_2st(intensities[offset]) for offset in offsets],
        )


def register() -> None:
    creators = Creators(
        config_creator=Cest15NTestConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest15NTestSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
