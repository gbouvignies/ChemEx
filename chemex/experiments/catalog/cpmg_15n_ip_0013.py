from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import Literal

import numpy as np
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
Delays = tuple[dict[float, float], dict[float, float], list[float]]


EXPERIMENT_NAME = "cpmg_15n_ip_0013"


class Cpmg15N0013IpSettings(CpmgSettings):
    name: Literal["cpmg_15n_ip_0013"]
    time_t2: float
    carrier: float
    pw90: float
    time_equil: float = 0.0
    ncyc_max: int
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def t_pos(self) -> float:
        return 4.0 * self.pw90 / np.pi

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cpmg15N0013IpConfig(
    ExperimentConfig[Cpmg15N0013IpSettings, RelaxationDataSettings]
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg15N0013IpConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Cpmg15N0013IpSequence:
    settings: Cpmg15N0013IpSettings

    def _get_delays(self, ncycs: np.ndarray) -> Delays:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - 0.75 * self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        deltas = {
            ncyc: self.settings.pw90 * (self.settings.ncyc_max - ncyc)
            + self.settings.time_equil
            for ncyc in ncycs_no_ref
        }
        deltas[0.0] = (
            self.settings.pw90 * (self.settings.ncyc_max - 1) + self.settings.time_equil
        )
        delays = [
            self.settings.t_neg,
            self.settings.t_pos,
            self.settings.time_equil,
            *deltas.values(),
            *tau_cps.values(),
        ]
        return tau_cps, deltas, delays

    def _get_phases(self, ncyc):
        cp_phases = np.array(
            [
                [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
                [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
            ]
        )
        indexes = np.flip(np.arange(2 * int(ncyc)))
        return np.take(cp_phases, indexes, mode="wrap", axis=1)

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_neg = delays[self.settings.t_neg]
        d_pos = delays[self.settings.t_pos]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i

        # Getting the starting magnetization
        start = spectrometer.get_equilibrium()

        # Calculating the instensities as a function of ncyc
        intensities = {
            0.0: spectrometer.detect(
                d_delta[0]
                @ p90[3]
                @ p180[[0, 3]]
                @ d_pos
                @ p180[[0, 1]]
                @ p90[1]
                @ start
            )
        }
        for ncyc in set(ncycs) - {0.0}:
            phases = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = reduce(np.matmul, echo[phases.T])
            intensities[ncyc] = spectrometer.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start
            )

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=Cpmg15N0013IpConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg15N0013IpSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
