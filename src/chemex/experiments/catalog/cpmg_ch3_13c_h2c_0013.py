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

EXPERIMENT_NAME = "cpmg_ch3_13c_h2c_0013"


class CpmgCh313CH2c0013Settings(CpmgSettingsEvenNcycs):
    """Settings for CH3 13C H2C CPMG experiment with 0013 variant."""

    name: Literal["cpmg_ch3_13c_h2c_0013"]
    time_t2: Delay = Field(description="Total CPMG relaxation delay (seconds)")
    carrier: Frequency = Field(description="13C carrier frequency (Hz)")
    pw90: PulseWidth = Field(description="13C 90-degree pulse width (seconds)")
    taub: float = Field(
        default=2.0e-3,
        gt=0.0,
        description="Delay for coherence evolution (seconds)",
    )
    time_equil: float = Field(
        default=5.0e-3,
        ge=0.0,
        description="Equilibration delay (seconds)",
    )
    time_grad: float = Field(
        default=1.2e-3,
        ge=0.0,
        description="Gradient duration (seconds)",
    )

    @computed_field  # type: ignore[misc]
    @property
    def t_neg(self) -> float:
        """Calculate negative delay compensation for pulse imperfections.

        Returns:
            Negative delay time in seconds.

        """
        return -2.0 * self.pw90 / np.pi

    @computed_field  # type: ignore[misc]
    @property
    def t_pos(self) -> float:
        """Calculate positive delay compensation for pulse imperfections.

        Returns:
            Positive delay time in seconds.

        """
        return 4.0 * self.pw90 / np.pi

    @computed_field  # type: ignore[misc]
    @property
    def taub_eff(self) -> float:
        """Calculate effective taub delay accounting for pulse width.

        Returns:
            Effective delay time in seconds.

        """
        return self.taub - 2.0 * self.pw90 / np.pi

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Initial magnetization terms for the experiment.

        Returns:
            List of initial state terms for the Liouvillian calculation.

        """
        return [f"2izsz{self.suffix_start}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection mode for the observable magnetization.

        Returns:
            Detection term for the Liouvillian calculation.

        """
        return f"[iz{self.suffix_detect}]"


class CpmgCh313CH2c0013Config(
    ExperimentConfiguration[
        CpmgCh313CH2c0013Settings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: CpmgCh313CH2c0013Config,
    spin_system: SpinSystem,
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


class CpmgCh313CH2c0013Sequence:
    """Sequence for CH3 13C H2C CPMG experiment with 0013 variant."""

    def __init__(self, settings: CpmgCh313CH2c0013Settings) -> None:
        self.settings = settings
        self._phase_cache: dict[float, tuple[Array, Array]] = {}

    def _get_delays(
        self, ncycs: Array
    ) -> tuple[dict[float, float], dict[float, float], list[float]]:
        ncycs_above_one = ncycs[ncycs > 1]
        ncyc_max = max(ncycs)
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - 0.75 * self.settings.pw90
            for ncyc in ncycs_above_one
        }
        tau_cps[1.0] = (
            self.settings.time_t2 / 4.0 - (2.0 / np.pi + 1.5) * self.settings.pw90
        )

        ncycs_1 = np.floor(np.floor(ncycs * 0.5) * 2 + 0.1)
        ncycs_2 = np.floor(np.floor((ncycs + 1) * 0.5) * 2 + 0.1)
        ncycs_3 = np.unique(np.concatenate((ncycs_1, ncycs_2)))
        ncycs_4 = ncycs_3[ncycs_3 > 0]
        deltas = {
            ncyc: 0.5 * self.settings.pw90 * (ncyc_max + 1 - ncyc) for ncyc in ncycs_4
        }
        deltas[0.0] = 0.5 * self.settings.pw90 * (ncyc_max - 1)

        delays = [
            self.settings.t_neg,
            self.settings.t_pos,
            self.settings.taub_eff,
            self.settings.time_equil,
            self.settings.time_grad,
            *deltas.values(),
            *tau_cps.values(),
        ]
        return tau_cps, deltas, delays

    # Define [0013] phase cycle for CPMG pulses
    def _get_phases(self, ncyc: float) -> tuple[Array, Array]:
        ncyc_key = float(ncyc)
        if ncyc_key not in self._phase_cache:
            cp_phases1 = np.array(
                [
                    [1, 1, 2, 0],
                    [2, 0, 3, 3],
                ],
            )
            cp_phases2 = np.array(
                [
                    [2, 0, 1, 1],
                    [3, 3, 2, 0],
                ],
            )
            indexes = np.arange(int(ncyc))
            phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap", axis=1)
            phases2 = np.take(cp_phases2, np.flip(indexes), mode="wrap", axis=1)
            self._phase_cache[ncyc_key] = (phases1, phases2)
        return self._phase_cache[ncyc_key]

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
        d_neg = delays[self.settings.t_neg]
        d_pos = delays[self.settings.t_pos]
        d_eq = delays[self.settings.time_equil]
        d_taub = delays[self.settings.taub_eff]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_grad = delays[self.settings.time_grad]

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180_sx = spectrometer.perfect180_s[0]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(terms=self.settings.start_terms)

        # Calculating the p-element with additional purge gradients
        zfilter = spectrometer.zfilter
        p_element = (
            d_grad
            @ zfilter
            @ p90[1]
            @ d_taub
            @ p90[0]
            @ p180_sx
            @ p90[0]
            @ d_taub
            @ p90[0]
            @ p180_sx
            @ d_grad
            @ zfilter
        )

        # Calculating the instensities as a function of ncyc
        intensities = {
            0.0: spectrometer.detect(
                d_eq
                @ d_delta[0]
                @ p90[0]
                @ p180[[0, 3]]
                @ d_pos
                @ p180[[2, 3]]
                @ (p90[0] - p90[2])
                * 0.5
                @ p_element
                @ d_delta[0]
                @ p90[0]
                @ p180[[1, 0]]
                @ d_pos
                @ p180[[1, 2]]
                @ (p90[0] - p90[2])
                * 0.5
                @ d_eq
                @ start,
            ),
        }

        for ncyc in set(ncycs) - {0.0}:
            ncyc_1 = np.floor(np.floor(ncyc * 0.5) * 2 + 0.1)
            ncyc_2 = np.floor(np.floor((ncyc + 1) * 0.5) * 2 + 0.1)
            phases1, phases2 = self._get_phases(ncyc_1)
            phases3, phases4 = self._get_phases(ncyc_2)

            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            if ncyc_1 > 0:
                cpmg1 = d_neg @ reduce(np.matmul, echo[phases1.T]) @ d_neg
                cpmg2 = d_neg @ reduce(np.matmul, echo[phases2.T]) @ d_neg
            else:
                cpmg1 = p180[[1, 0]] @ d_pos @ p180[[1, 2]]
                cpmg2 = p180[[0, 3]] @ d_pos @ p180[[2, 3]]
            cpmg3 = d_neg @ reduce(np.matmul, echo[phases3.T]) @ d_neg
            cpmg4 = d_neg @ reduce(np.matmul, echo[phases4.T]) @ d_neg

            if int(ncyc) % 2 == 0:
                if int(ncyc_1) % 4 == 0:
                    intst1 = spectrometer.detect(
                        d_eq
                        @ d_delta[ncyc_1]
                        @ p90[0]
                        @ cpmg1
                        @ (p90[0] - p90[2])
                        * 0.5
                        @ p_element
                        @ d_delta[ncyc_1]
                        @ p90[0]
                        @ cpmg1
                        @ (p90[0] - p90[2])
                        * 0.5
                        @ d_eq
                        @ start
                    )
                else:
                    intst1 = spectrometer.detect(
                        d_eq
                        @ d_delta[ncyc_1]
                        @ p90[0]
                        @ cpmg2
                        @ (p90[0] - p90[2])
                        * 0.5
                        @ p_element
                        @ d_delta[ncyc_1]
                        @ p90[0]
                        @ cpmg1
                        @ (p90[0] - p90[2])
                        * 0.5
                        @ d_eq
                        @ start
                    )

                intst2 = intst1
            elif int(ncyc_1) % 4 == 0:
                intst1 = spectrometer.detect(
                    d_eq
                    @ d_delta[ncyc_2]
                    @ p90[0]
                    @ cpmg3
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ p_element
                    @ d_delta[ncyc_1]
                    @ p90[0]
                    @ cpmg1
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ d_eq
                    @ start
                )

                intst2 = spectrometer.detect(
                    d_eq
                    @ d_delta[ncyc_1]
                    @ p90[0]
                    @ cpmg2
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ p_element
                    @ d_delta[ncyc_2]
                    @ p90[0]
                    @ cpmg3
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ d_eq
                    @ start
                )
            else:
                intst1 = spectrometer.detect(
                    d_eq
                    @ d_delta[ncyc_2]
                    @ p90[0]
                    @ cpmg4
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ p_element
                    @ d_delta[ncyc_1]
                    @ p90[0]
                    @ cpmg1
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ d_eq
                    @ start
                )

                intst2 = spectrometer.detect(
                    d_eq
                    @ d_delta[ncyc_1]
                    @ p90[0]
                    @ cpmg1
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ p_element
                    @ d_delta[ncyc_2]
                    @ p90[0]
                    @ cpmg3
                    @ (p90[0] - p90[2])
                    * 0.5
                    @ d_eq
                    @ start
                )

            intensities[ncyc] = (intst1 + intst2) * 0.5

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return metadata == 0


def register() -> None:
    creators = Creators(
        config_creator=CpmgCh313CH2c0013Config,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CpmgCh313CH2c0013Sequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
