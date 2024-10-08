from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, reduce
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import (
    CpmgSettingsEvenNcycs,
)
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cpmg import CpmgPlotter
from chemex.printers.data import CpmgPrinter
from chemex.typing import ArrayBool, ArrayFloat, ArrayInt

EXPERIMENT_NAME = "cpmg_hn_dq_zq"

NU_CPMG_LIMIT_1 = 51.0
NU_CPMG_LIMIT_2 = 255.0


class CpmgHNDqZqSettings(CpmgSettingsEvenNcycs):
    name: Literal["cpmg_hn_dq_zq"]
    time_t2: float
    carrier_h: float
    carrier_n: float
    pw90_h: float
    pw90_n: float
    dq_flg: bool

    @cached_property
    def start_terms(self) -> list[str]:
        return [f"2ixsx{self.suffix}"]

    @cached_property
    def detection(self) -> str:
        if self.dq_flg:
            return f"[2ixsx_{self.observed_state}] - [2iysy_{self.observed_state}]"
        return f"[2ixsx_{self.observed_state}] + [2iysy_{self.observed_state}]"


class CpmgHNDqZqConfig(
    ExperimentConfiguration[
        CpmgHNDqZqSettings, ConditionsWithValidations, RelaxationDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r2mq_is_{state}", f"mu_is_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: CpmgHNDqZqConfig,
    spin_system: SpinSystem,
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

    def _get_tau_cps(self, ncycs: ArrayFloat) -> dict[float, float]:
        ncycs_no_ref = ncycs[ncycs > 0]
        return dict(
            zip(
                ncycs_no_ref,
                self.settings.time_t2 / (4.0 * ncycs_no_ref)
                - 7.0 / 3.0 * self.settings.pw90_n,
                strict=True,
            ),
        )

    def _get_phases(self, ncyc: float) -> tuple[ArrayInt, ArrayInt]:
        nu_cpmg = ncyc / self.settings.time_t2
        if nu_cpmg < NU_CPMG_LIMIT_1:
            cp_phases1 = [0, 1, 0, 1]
            cp_phases2 = [1, 0, 1, 0]
        elif nu_cpmg < NU_CPMG_LIMIT_2:
            cp_phases1 = [0]
            cp_phases2 = [1]
        else:
            cp_phases1 = [0, 1, 0, 1, 1, 0, 1, 0]
            cp_phases2 = [1, 0, 1, 0, 0, 1, 0, 1]
        indexes = np.arange(2 * int(ncyc))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap")
        phases2 = np.take(cp_phases2, np.flip(indexes), mode="wrap")
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps = self._get_tau_cps(ncycs)
        ncyc_list, delay_list = zip(*tau_cps.items(), strict=True)
        d_cp = dict(zip(ncyc_list, spectrometer.delays(delay_list), strict=True))

        # Calculation of the spectrometers corresponding to all the pulses
        p9024090_1 = spectrometer.p9024090_nh_1[[0, 1], [0, 1]]
        p9024090_2 = spectrometer.p9024090_nh_2[[0, 1], [0, 1]]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(self.settings.start_terms)

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
    def is_reference(metadata: ArrayFloat) -> ArrayBool:
        return metadata == 0


def register() -> None:
    creators = Creators(
        config_creator=CpmgHNDqZqConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgHNDqZqSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
