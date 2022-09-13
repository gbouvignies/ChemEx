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


EXPERIMENT_NAME = "cpmg_ch3_1h_sq"


class CpmgCh31HSqSettings(CpmgSettings):
    name: Literal["cpmg_ch3_1h_sq"]
    time_t2: float
    carrier: float
    pw90: float
    ncyc_max: int
    taua: float = 2.0e-3
    comp180_flg: bool = True
    ipap_flg: bool = False
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[iy_{self.observed_state}]"


class CpmgCh31HSqConfig(ExperimentConfig[CpmgCh31HSqSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state

        to_be_fitted = ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])

        if self.experiment.ipap_flg:
            to_be_fitted.rates.append(f"r2a_i_{state}")
            to_be_fitted.model_free.append(f"s2_{state}")

        return to_be_fitted


def build_spectrometer(
    config: CpmgCh31HSqConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="hc")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CpmgCh31HSqSequence:
    settings: CpmgCh31HSqSettings

    def _get_delays(self, ncycs: np.ndarray) -> tuple[dict[float, float], list[float]]:
        frac = 7.0 / 3.0 if self.settings.comp180_flg else 1.0
        ncyc_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: (
                self.settings.time_t2
                - frac * self.settings.pw90 * 4.0 * self.settings.ncyc_max
            )
            / (4.0 * ncyc)
            for ncyc in ncyc_no_ref
        }
        delays = [self.settings.taua, *tau_cps.values()]
        return tau_cps, delays

    def _get_phases(self) -> tuple[np.ndarray, np.ndarray]:
        cp_phases1 = np.array([[0, 1], [1, 0]])
        cp_phases2 = np.array([[0, 3], [3, 0]])
        indexes = np.arange(int(self.settings.ncyc_max))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap", axis=1)
        phases2 = np.take(cp_phases2, indexes, mode="wrap", axis=1)
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_taua = delays[self.settings.taua]

        # Calculation of the spectrometers corresponding to all the pulses
        perfect180y = spectrometer.perfect180_i[1]
        p180 = spectrometer.p180_i
        p180c_py = spectrometer.p9018090_i_1[1]
        p180c_my = spectrometer.p9018090_i_2[3]
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        if self.settings.comp180_flg:
            p180_cp1 = spectrometer.p9024090_i_1
            p180_cp2 = spectrometer.p9024090_i_2
        else:
            p180_cp1 = p180_cp2 = p180

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(terms=["iy"])
        start = perfect180y @ d_taua @ d_taua @ start
        start = spectrometer.keep(start, ["2ixsz_a", "2iysz_a"])

        # Calculating the intensities as a function of ncyc
        if self.settings.ipap_flg:
            intensities = {
                0.0: spectrometer.detect(
                    d_taua @ (p180pmy @ p180c_py + p180c_my @ p180pmy) @ d_taua @ start
                )
            }
        else:
            intensities = {
                0.0: spectrometer.detect(d_taua @ d_taua @ p180pmy @ p180c_py @ start)
            }

        phases1, phases2 = self._get_phases()
        for ncyc in set(ncycs) - {0.0}:
            phases1_1 = phases1.T[-int(ncyc) :]
            phases1_2 = phases1.T[: -int(ncyc)]
            phases2_1 = phases2.T[: int(ncyc)]
            phases2_2 = phases2.T[int(ncyc) :]
            echo1 = d_cp[ncyc] @ p180_cp1 @ d_cp[ncyc]
            echo2 = d_cp[ncyc] @ p180_cp2 @ d_cp[ncyc]
            cpmg1 = reduce(np.matmul, echo1[phases1_1])
            cpmg2 = reduce(np.matmul, echo2[phases2_1])
            if ncyc < self.settings.ncyc_max:
                cpmg1 = reduce(np.matmul, p180_cp1[phases1_2]) @ cpmg1
                cpmg2 = cpmg2 @ reduce(np.matmul, p180_cp2[phases2_2])
            centre = cpmg2 @ p180pmy @ cpmg1
            if self.settings.ipap_flg:
                intensities[ncyc] = spectrometer.detect(
                    d_taua @ (centre @ p180c_py + p180c_my @ centre) @ d_taua @ start
                )
            else:
                intensities[ncyc] = spectrometer.detect(
                    d_taua @ d_taua @ centre @ p180c_py @ start
                )

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=CpmgCh31HSqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgCh31HSqSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
