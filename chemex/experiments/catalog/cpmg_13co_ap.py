from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from numpy.typing import NDArray

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettings
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


EXPERIMENT_NAME = "cpmg_13co_ap"


class Cpmg13CoApSettings(CpmgSettings):
    name: Literal["cpmg_13co_ap"]
    time_t2: float
    carrier: float
    pw90: float
    time_equil: float = 0.0
    refocusing: bool = False
    taucc: float = 9.09e-3
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def detection(self) -> str:
        return f"[2izsz_{self.observed_state}]"

    @property
    def even_ncycs(self) -> bool:
        return self.refocusing


class Cpmg13CoApConfig(ExperimentConfig[Cpmg13CoApSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg13CoApConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Cpmg13CoApSequence:
    settings: Cpmg13CoApSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        delays = [
            self.settings.t_neg,
            self.settings.time_equil,
            self.settings.taucc,
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
        d_taucc = delays[self.settings.taucc]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        perfect180x = spectrometer.perfect180_i[0]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(["2izsz"])

        # Calculate the flip block
        if self.settings.refocusing:
            p_flip0 = p90[3] @ d_taucc @ p180pmy @ d_taucc @ p90[1]
            p_flip = d_neg @ p_flip0 @ d_neg
        else:
            p_flip = p_flip0 = p180pmy

        # Calculating the instensities as a function of ncyc
        intensities = {
            0.0: spectrometer.detect(d_eq @ p90[1] @ p_flip0 @ p90[1] @ start)
        }

        part1 = d_neg @ p90[1] @ start
        part2 = d_eq @ p90[1] @ d_neg

        for ncyc in set(ncycs) - {0.0}:
            echo = d_cp[ncyc] @ perfect180x @ d_cp[ncyc]
            cpmg = matrix_power(echo, int(ncyc))
            intensities[ncyc] = spectrometer.detect(
                part2 @ cpmg @ p_flip @ cpmg @ part1
            )

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=Cpmg13CoApConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg13CoApSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
