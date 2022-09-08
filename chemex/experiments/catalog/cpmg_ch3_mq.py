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


EXPERIMENT_NAME = "cpmg_ch3_mq"


class CpmgCh3MqSettings(CpmgSettings):
    name: Literal["cpmg_ch3_mq"]
    time_t2: float
    t_zeta: float = 1.0 / (8.0 * 125.3)
    small_protein: bool = False
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[2iysx_{self.observed_state}]"


class CpmgCh3MqConfig(ExperimentConfig[CpmgCh3MqSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2mq_is_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: CpmgCh3MqConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixysxy", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CpmgCh3MqSequence:
    settings: CpmgCh3MqSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        ncycs_no_refs = ncycs[ncycs > 0.0]
        tau_cps = {ncyc: self.settings.time_t2 / (4.0 * ncyc) for ncyc in ncycs_no_refs}
        delays = [self.settings.t_zeta, *tau_cps.values()]
        return tau_cps, delays

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_zeta = delays[self.settings.t_zeta]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p180_sx = spectrometer.perfect180_s[0]
        p180_ix = spectrometer.perfect180_i[0]
        p180_iy = spectrometer.perfect180_i[1]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(terms=["2iysx"])

        # Calculating the instensities as a function of ncyc
        part1 = start
        if self.settings.small_protein:
            part1 = d_zeta @ p180_sx @ p180_ix @ d_zeta @ part1
        intensities = {0.0: spectrometer.detect(p180_sx @ part1)}
        for ncyc in set(ncycs) - {0.0}:
            echo = d_cp[ncyc] @ p180_iy @ d_cp[ncyc]
            cpmg = matrix_power(echo, int(ncyc))
            intensities[ncyc] = spectrometer.detect(cpmg @ p180_sx @ cpmg @ part1)

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=CpmgCh3MqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgCh3MqSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
