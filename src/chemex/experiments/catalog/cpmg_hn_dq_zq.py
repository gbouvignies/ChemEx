from __future__ import annotations

from functools import reduce
from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettingsEvenNcycs
from chemex.configuration.types import Delay, Frequency, PulseWidth
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
from chemex.typing import Array

EXPERIMENT_NAME = "cpmg_hn_dq_zq"

NU_CPMG_LIMIT_1 = 51.0
NU_CPMG_LIMIT_2 = 255.0


class CpmgHNDqZqSettings(CpmgSettingsEvenNcycs):
    """Settings for HN DQ/ZQ CPMG relaxation dispersion experiment."""

    name: Literal["cpmg_hn_dq_zq"]
    time_t2: Delay = Field(description="Total CPMG relaxation delay (seconds)")
    carrier_h: Frequency = Field(description="1H carrier frequency (Hz)")
    carrier_n: Frequency = Field(description="15N carrier frequency (Hz)")
    pw90_h: PulseWidth = Field(description="1H 90-degree pulse width (seconds)")
    pw90_n: PulseWidth = Field(description="15N 90-degree pulse width (seconds)")
    dq_flg: bool = Field(
        description="Flag for double-quantum (True) vs zero-quantum (False)"
    )

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Initial magnetization terms for the experiment.

        Returns:
            List of initial state terms for the Liouvillian calculation.

        """
        return [f"2ixsx{self.suffix_start}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection mode for the observable magnetization.

        Returns:
            Detection term for the Liouvillian calculation.

        """
        suffix = self.suffix_detect
        if self.dq_flg:
            return f"[2ixsx{suffix}] - [2iysy{suffix}]"
        return f"[2ixsx{suffix}] + [2iysy{suffix}]"


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


class CpmgHNDqZqSequence:
    """Sequence for HN DQ/ZQ CPMG relaxation dispersion experiment."""

    def __init__(self, settings: CpmgHNDqZqSettings) -> None:
        self.settings = settings
        self._phase_cache: dict[float, tuple[Array, Array]] = {}

    def _get_tau_cps(self, ncycs: Array) -> dict[float, float]:
        ncycs_no_ref = ncycs[ncycs > 0]
        return dict(
            zip(
                ncycs_no_ref,
                self.settings.time_t2 / (4.0 * ncycs_no_ref)
                - 7.0 / 3.0 * self.settings.pw90_n,
                strict=True,
            ),
        )

    def _get_phases(self, ncyc: float) -> tuple[Array, Array]:
        ncyc_key = float(ncyc)
        if ncyc_key not in self._phase_cache:
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
            self._phase_cache[ncyc_key] = (phases1, phases2)
        return self._phase_cache[ncyc_key]

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
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
    def is_reference(metadata: Array) -> Array:
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
