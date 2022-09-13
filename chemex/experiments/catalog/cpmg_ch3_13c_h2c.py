from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from numpy.typing import NDArray

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettingsEvenNcycs
from chemex.configuration.experiment import ExperimentConfig
from chemex.configuration.experiment import ToBeFitted
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators
from chemex.experiments.factories import factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.liouvillian import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters import CpmgPlotter
from chemex.printers.data import CpmgPrinter

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "cpmg_ch3_13c_h2c"


class CpmgCh313CH2cSettings(CpmgSettingsEvenNcycs):
    name: Literal["cpmg_ch3_13c_h2c"]
    time_t2: float
    carrier: float
    pw90: float
    taub: float = 2.0e-3
    time_equil: float = 0.0
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def start(self) -> list[str]:
        return [f"2izsz_{self.observed_state}"]

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class CpmgCh313CH2cConfig(
    ExperimentConfig[CpmgCh313CH2cSettings, RelaxationDataSettings]
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: CpmgCh313CH2cConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CpmgCh313CH2cSequence:
    settings: CpmgCh313CH2cSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        delays = [
            self.settings.t_neg,
            self.settings.taub,
            self.settings.time_equil,
            *tau_cps.values(),
        ]
        return tau_cps, delays

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_neg = delays[self.settings.t_neg]
        d_eq = delays[self.settings.time_equil]
        d_taub = delays[self.settings.taub]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180_sx = spectrometer.perfect180_s[0]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(terms=self.settings.start)

        # Calculating the p-element
        palmer = d_taub @ p90[0] @ p180_sx @ p90[0] @ d_taub

        # Calculating the instensities as a function of ncyc
        intensities = {
            0.0: spectrometer.detect(d_eq @ p90[1] @ palmer @ p90[0] @ start)
        }
        part1 = d_neg @ p90[0] @ start
        part2 = d_eq @ p90[1] @ d_neg
        for ncyc in set(ncycs) - {0.0}:
            echo = d_cp[ncyc] @ p180[[1, 0]] @ d_cp[ncyc]
            cpmg1, cpmg2 = matrix_power(echo, int(ncyc))
            end = part2 @ cpmg2 @ palmer @ cpmg1 @ part1
            intensities[ncyc] = spectrometer.detect(end)

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=CpmgCh313CH2cConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgCh313CH2cSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
