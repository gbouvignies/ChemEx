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


EXPERIMENT_NAME = "cpmg_15n_tr"


class Cpmg15NTrSettings(CpmgSettings):
    name: Literal["cpmg_15n_tr"]
    time_t2: float
    carrier: float
    pw90: float
    time_equil: float = 0.0
    taub: float = 2.68e-3
    antitrosy: bool = False
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def taub_eff(self) -> float:
        return self.taub - 2.0 * self.pw90 - 2.0 * self.pw90 / np.pi

    @property
    def detection(self) -> str:
        if self.antitrosy:
            return f"[2izsz_{self.observed_state}] + [iz_{self.observed_state}]"
        else:
            return f"[2izsz_{self.observed_state}] - [iz_{self.observed_state}]"


class Cpmg15NTrConfig(ExperimentConfig[Cpmg15NTrSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state

        to_be_fitted = ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])

        if self.experiment.antitrosy:
            to_be_fitted.rates.append(f"etaxy_i_{state}")
            to_be_fitted.model_free.append(f"s2_{state}")

        return to_be_fitted


def build_spectrometer(
    config: Cpmg15NTrConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Cpmg15NTrSequence:
    settings: Cpmg15NTrSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        delays = [
            self.settings.t_neg,
            self.settings.taub_eff,
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
        d_taub = delays[self.settings.taub_eff]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180_sx = spectrometer.perfect180_s[0]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(["2izsz"])

        # Calculating the p-element
        p1, p2 = (2, 1) if self.settings.antitrosy else (1, 0)
        palmer0 = (
            p180_sx @ d_taub @ p90[p1] @ p90[p2] @ p180_sx @ p90[p2] @ p90[p1] @ d_taub
        )
        palmer = np.mean(p90[[0, 2]] @ palmer0 @ p90[[1, 3]], axis=0)

        # Calculating the instensities as a function of ncyc
        part1 = p90[0] @ start
        part2 = d_eq @ p90[1]

        intst = {0.0: spectrometer.detect(part2 @ palmer0 @ part1)}

        for ncyc in set(ncycs) - {0.0}:
            echo = d_cp[ncyc] @ p180[[1, 0]] @ d_cp[ncyc]
            cpmg1, cpmg2 = d_neg @ matrix_power(echo, int(ncyc)) @ d_neg
            intst[ncyc] = spectrometer.detect(part2 @ cpmg2 @ palmer @ cpmg1 @ part1)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=Cpmg15NTrConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg15NTrSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
