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


EXPERIMENT_NAME = "cpmg_ch3_1h_tq"


class CpmgCh31HTqSettings(CpmgSettings):
    name: Literal["cpmg_ch3_1h_tq"]
    time_t2: float
    carrier: float
    pw90: float
    tauc: float = 0.67e-3  # ~ 1/(12*J[HC])
    comp180_flg: bool = True
    ipap_flg: bool = False
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[2ixsz_{self.observed_state}]"


class CpmgCh31HTqConfig(ExperimentConfig[CpmgCh31HTqSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state

        to_be_fitted = ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])

        if self.experiment.ipap_flg:
            to_be_fitted.rates.append(f"r2a_i_{state}")
            to_be_fitted.model_free.append(f"s2_{state}")

        return to_be_fitted


def build_spectrometer(
    config: CpmgCh31HTqConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", extension="tq", spin_system="hc")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CpmgCh31HTqSequence:
    settings: CpmgCh31HTqSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        ncyc_no_ref = ncycs[ncycs > 0.0]
        factor = 2.0 if self.settings.comp180_flg else 1.0
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - factor * self.settings.pw90
            for ncyc in ncyc_no_ref
        }
        delays = [self.settings.tauc, *tau_cps.values()]
        return tau_cps, delays

    def _get_phases(self, ncyc: float) -> tuple[np.ndarray, np.ndarray]:
        cp_phases1 = [0, 1, 0, 1, 1, 0, 1, 0, 2, 3, 2, 3, 3, 2, 3, 2]
        cp_phases2 = [0, 3, 0, 3, 3, 0, 3, 0, 2, 1, 2, 1, 1, 2, 1, 2]
        indexes = np.arange(int(ncyc))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap")
        phases2 = np.take(cp_phases2, indexes, mode="wrap")
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(["2ixsz"], atom="h")

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_tauc = delays[self.settings.tauc]

        # Calculation of the spectrometers corresponding to all the pulses
        p180 = spectrometer.p180_i
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        p180_sx = spectrometer.perfect180_s[0]
        if self.settings.comp180_flg:
            p180_cp1 = spectrometer.p9018090_i_1
            p180_cp2 = spectrometer.p9018090_i_2
        else:
            p180_cp1 = p180_cp2 = p180

        # Calculating the intensities as a function of ncyc
        if self.settings.ipap_flg:
            part1 = d_tauc @ start
            part2 = d_tauc
            centre0 = 0.5 * (p180pmy + p180_sx @ p180pmy @ p180_sx)
        else:
            part1 = start
            part2 = spectrometer.identity
            centre0 = p180pmy

        intensities = {0.0: spectrometer.detect(part2 @ centre0 @ part1)}

        for ncyc in set(ncycs) - {0.0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo1 = d_cp[ncyc] @ p180_cp1 @ d_cp[ncyc]
            echo2 = d_cp[ncyc] @ p180_cp2 @ d_cp[ncyc]
            cpmg1 = reduce(np.matmul, echo1[phases1])
            cpmg2 = reduce(np.matmul, echo2[phases2])
            centre = cpmg2 @ p180pmy @ cpmg1
            if self.settings.ipap_flg:
                centre = 0.5 * (centre + p180_sx @ centre @ p180_sx)
            intensities[ncyc] = spectrometer.detect(part2 @ centre @ part1)

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=CpmgCh31HTqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgCh31HTqSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
