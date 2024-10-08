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

EXPERIMENT_NAME = "coscest_13c"

OFFSET_REF = 1e4


class CosCest13CSettings(MFCestSettings):
    name: Literal["coscest_13c"]
    time_t1: float
    time_equil: float = 0.0
    carrier: float
    cos_n: int
    cos_res: int = 10
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11

    @cached_property
    def start_terms(self) -> list[str]:
        return [f"iz{self.suffix}"]

    @cached_property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class CosCest13CConfig(
    ExperimentConfiguration[
        CosCest13CSettings, ConditionsWithValidations, CestDataSettings
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
    config: CosCest13CConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res
    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        symbol = spin_system.symbols["i"]
        atom = spin_system.atoms["i"]
        liouvillian.jeff_i = get_multiplet(symbol, atom.name)

    return spectrometer


@dataclass
class CosCest13CSequence:
    settings: CosCest13CSettings

    @staticmethod
    def is_reference(metadata: ArrayFloat) -> ArrayBool:
        return np.abs(metadata) > OFFSET_REF

    def _calc_cosine_shape(self, spectrometer: Spectrometer) -> ArrayFloat:
        time_t1 = self.settings.time_t1
        sw = self.settings.sw
        cos_n = self.settings.cos_n
        cos_res = self.settings.cos_res

        dt = 1.0 / (cos_res * sw)
        n_periods = int(time_t1 * sw)
        n_left = int((time_t1 * sw - n_periods) * cos_res)
        double_periods = n_periods // 2
        extra_period = n_periods % 2
        phase1 = 2 if cos_n % 2 == 0 else 0
        phase_left = 0 if extra_period else phase1

        grid = np.linspace(-np.pi, np.pi, cos_res, endpoint=False)

        n_values = (np.arange(cos_n) - 0.5 * (cos_n - 1)).reshape(-1, 1)
        amplitudes = np.cos(n_values * grid).sum(axis=0)
        phases = np.zeros(cos_res)

        base_pulse = spectrometer.shaped_pulse_i(cos_res * dt, amplitudes, phases)

        pulse = matrix_power(base_pulse[0] @ base_pulse[phase1], double_periods)

        if extra_period:
            pulse = base_pulse[phase1] @ pulse

        if n_left:
            pulse_left = spectrometer.shaped_pulse_i(
                n_left * dt,
                amplitudes[:n_left],
                phases[:n_left],
            )
            pulse = pulse_left[phase_left] @ pulse

        return pulse

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        offsets = data.metadata

        start = spectrometer.get_equilibrium()

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

            intensities[offset] = d_eq @ self._calc_cosine_shape(spectrometer) @ start

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets],
        )


def register() -> None:
    creators = Creators(
        config_creator=CosCest13CConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CosCest13CSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
