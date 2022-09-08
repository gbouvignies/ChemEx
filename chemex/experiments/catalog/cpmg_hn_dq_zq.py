from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import Literal

import numpy as np
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


EXPERIMENT_NAME = "cpmg_hn_dq_zq"


class CpmgHNDqZqSettings(CpmgSettingsEvenNcycs):
    name: Literal["cpmg_hn_dq_zq"]
    time_t2: float
    carrier_h: float
    carrier_n: float
    pw90_h: float
    pw90_n: float
    dq_flg: bool
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        if self.dq_flg:
            return f"[2ixsx_{self.observed_state}] - [2iysy_{self.observed_state}]"
        else:
            return f"[2ixsx_{self.observed_state}] + [2iysy_{self.observed_state}]"


class CpmgHNDqZqConfig(ExperimentConfig[CpmgHNDqZqSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r2mq_is_{state}", f"mu_is_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: CpmgHNDqZqConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsxyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier_n
    spectrometer.carrier_s = settings.carrier_h
    spectrometer.b1_i = 1 / (4.0 * settings.pw90_n)
    spectrometer.b1_s = 1 / (4.0 * settings.pw90_h)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CpmgHNDqZqSequence:
    settings: CpmgHNDqZqSettings

    def _get_tau_cps(self, ncycs: np.ndarray) -> dict[float, float]:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        return dict(
            zip(
                ncycs_no_ref,
                self.settings.time_t2 / (4.0 * ncycs_no_ref)
                - 7.0 / 3.0 * self.settings.pw90_n,
            )
        )

    def _get_phases(self, ncyc: float) -> tuple[np.ndarray, np.ndarray]:
        nu_cpmg = ncyc / self.settings.time_t2
        if nu_cpmg < 51.0:
            cp_phases1 = [0, 1, 0, 1]
            cp_phases2 = [1, 0, 1, 0]
        elif nu_cpmg < 255.0:
            cp_phases1 = [0]
            cp_phases2 = [1]
        else:
            cp_phases1 = [0, 1, 0, 1, 1, 0, 1, 0]
            cp_phases2 = [1, 0, 1, 0, 0, 1, 0, 1]
        indexes = np.arange(2 * int(ncyc))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap")
        phases2 = np.take(cp_phases2, np.flip(indexes), mode="wrap")
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps = self._get_tau_cps(ncycs)
        d_cp = dict(zip(tau_cps.keys(), spectrometer.delays(list(tau_cps.values()))))

        # Calculation of the spectrometers corresponding to all the pulses
        p9024090_1 = spectrometer.p9024090_nh_1[[0, 1], [0, 1]]
        p9024090_2 = spectrometer.p9024090_nh_2[[0, 1], [0, 1]]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(["2ixsx"])

        # Calculating the cpmg trains
        intensities = {0.0: spectrometer.detect(start)}
        for ncyc in set(ncycs) - {0.0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo1 = d_cp[ncyc] @ p9024090_1 @ d_cp[ncyc]
            echo2 = d_cp[ncyc] @ p9024090_2 @ d_cp[ncyc]
            cpmg1 = reduce(np.matmul, echo1[phases1])
            cpmg2 = reduce(np.matmul, echo2[phases2])
            intensities[ncyc] = spectrometer.detect(0.5 * (cpmg1 + cpmg2) @ start)

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=CpmgHNDqZqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgHNDqZqSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
