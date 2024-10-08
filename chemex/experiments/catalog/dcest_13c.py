from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Literal

import numpy as np
from numpy.linalg import matrix_power

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import MFCestSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import CestFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.constants import get_multiplet
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.printers.data import CestPrinter
from chemex.typing import ArrayBool, ArrayFloat

EXPERIMENT_NAME = "dcest_13c"

OFFSET_REF = 1e4


class DCest13CSettings(MFCestSettings):
    name: Literal["dcest_13c"]
    time_t1: float
    time_equil: float = 0.0
    carrier: float
    pw90: float
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11

    @cached_property
    def pw_dante(self) -> float:
        return 4.0 * self.pw90 * self.b1_frq / self.sw

    @cached_property
    def tau_dante(self) -> float:
        return 1.0 / self.sw - self.pw_dante

    @cached_property
    def ncyc_dante(self) -> int:
        return int(self.time_t1 * self.sw + 0.1)

    @cached_property
    def start_terms(self) -> list[str]:
        return [f"iz{self.suffix}"]

    @cached_property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class DCest13CConfig(
    ExperimentConfiguration[
        DCest13CSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(config: DCest13CConfig, spin_system: SpinSystem) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1.0 / (4.0 * settings.pw90)
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res
    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        symbol = spin_system.symbols["i"]
        atom = spin_system.atoms["i"]
        liouvillian.jeff_i = get_multiplet(symbol, atom.name)

    return spectrometer


@dataclass
class DCest13CSequence:
    settings: DCest13CSettings

    @staticmethod
    def is_reference(metadata: ArrayFloat) -> ArrayBool:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(terms=self.settings.start_terms)

        d_eq = (
            spectrometer.delays(self.settings.time_equil)
            if self.settings.time_equil > 0
            else spectrometer.identity
        )

        intensities: dict[float, ArrayFloat] = {}

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


def register() -> None:
    creators = Creators(
        config_creator=DCest13CConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=DCest13CSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
