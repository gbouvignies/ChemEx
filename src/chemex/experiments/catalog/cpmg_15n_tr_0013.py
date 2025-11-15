from __future__ import annotations

from functools import reduce
from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettings
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

Delays = tuple[dict[float, float], dict[float, float], list[float]]


EXPERIMENT_NAME = "cpmg_15n_tr_0013"


class Cpmg15NTr0013Settings(CpmgSettings):
    """Settings for 15N TROSY-based CPMG experiment with 0013 variant."""

    name: Literal["cpmg_15n_tr_0013"]
    time_t2: Delay = Field(description="Total CPMG relaxation delay (seconds)")
    carrier: Frequency = Field(description="15N carrier frequency (Hz)")
    pw90: PulseWidth = Field(description="15N 90-degree pulse width (seconds)")
    ncyc_max: int = Field(gt=0, description="Maximum number of CPMG cycles")
    time_equil: float = Field(
        default=0.0,
        ge=0.0,
        description="Equilibration delay (seconds)",
    )
    taub: float = Field(
        default=2.68e-3,
        gt=0.0,
        description="Delay for TROSY evolution (seconds)",
    )
    antitrosy: bool = Field(
        default=False,
        description="Flag for antitrosy component selection",
    )
    s3e: bool = Field(
        default=True,
        description="Flag for sensitivity-enhanced version",
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
    def taub_eff(self) -> float:
        """Calculate effective taub delay accounting for pulse widths.

        Returns:
            Effective delay time in seconds.

        """
        return self.taub - 2.0 * self.pw90 - 2.0 * self.pw90 / np.pi

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Initial magnetization terms for the experiment.

        Returns:
            List of initial state terms for the Liouvillian calculation.

        """
        suffix = self.suffix_start
        if self.s3e:
            if self.antitrosy:
                return [f"2izsz{suffix}", f"iz{suffix}"]
            return [f"2izsz{suffix}", f"-iz{suffix}"]
        return [f"2izsz{suffix}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection mode for the observable magnetization.

        Returns:
            Detection term for the Liouvillian calculation.

        """
        suffix = self.suffix_detect
        if self.antitrosy:
            return f"[2izsz{suffix}] + [iz{suffix}]"
        return f"[2izsz{suffix}] - [iz{suffix}]"


class Cpmg15NTr0013Config(
    ExperimentConfiguration[
        Cpmg15NTr0013Settings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state

        to_be_fitted = ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])

        if self.experiment.antitrosy:
            to_be_fitted.rates.append(f"etaxy_i_{state}")
            to_be_fitted.model_free.append(f"s2_{state}")

        return to_be_fitted


def build_spectrometer(
    config: Cpmg15NTr0013Config,
    spin_system: SpinSystem,
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


class Cpmg15NTr0013Sequence:
    """Sequence for 15N TROSY-based CPMG experiment with 0013 variant."""

    def __init__(self, settings: Cpmg15NTr0013Settings) -> None:
        self.settings = settings

    def _get_delays(self, ncycs: Array) -> Delays:
        ncycs_no_ref = ncycs[ncycs > 0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - 0.75 * self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        deltas = {
            ncyc: 0.5 * self.settings.pw90 * (self.settings.ncyc_max - ncyc)
            for ncyc in ncycs_no_ref
        }
        deltas[0.0] = 0.5 * self.settings.pw90 * (self.settings.ncyc_max - 1)
        delays = [
            self.settings.t_neg,
            self.settings.taub_eff,
            self.settings.time_equil,
            *deltas.values(),
            *tau_cps.values(),
        ]
        return tau_cps, deltas, delays

    def _get_phases(self, ncyc: float) -> tuple[Array, Array]:
        cp_phases1 = np.array(
            [
                [1, 1, 0, 2, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 0, 2],
                [0, 2, 3, 3, 2, 0, 3, 3, 2, 0, 3, 3, 0, 2, 3, 3],
            ],
        )
        cp_phases2 = np.array(
            [
                [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
                [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
            ],
        )
        indexes = np.arange(int(ncyc))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap", axis=1)
        phases2 = np.take(cp_phases2, indexes, mode="wrap", axis=1)
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
        d_neg = delays[self.settings.t_neg]
        d_eq = delays[self.settings.time_equil]
        d_taub = delays[self.settings.taub_eff]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180_sx = spectrometer.perfect180_s[0]

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        # Calculating the p-element
        if self.settings.antitrosy:
            p180_is = p90[1] @ p90[0] @ p180_sx @ p90[0] @ p90[1]
        else:
            p180_is = p90[0] @ p90[1] @ p180_sx @ p90[1] @ p90[0]
        palmer = np.mean(
            p90[[0, 2]] @ d_taub @ p180_is @ d_taub @ p180_sx @ p90[[1, 3]],
            axis=0,
        )

        # Calculating the instensities as a function of ncyc
        intst = {
            0.0: spectrometer.detect(
                d_eq
                @ d_delta[0]
                @ p90[3]
                @ p180[[1, 0]]
                @ palmer
                @ p180[[0, 1]]
                @ p90[0]
                @ d_delta[0]
                @ start,
            ),
        }
        for ncyc in set(ncycs) - {0.0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg1 = d_neg @ reduce(np.matmul, echo[phases1.T]) @ d_neg
            cpmg2 = d_neg @ reduce(np.matmul, echo[phases2.T]) @ d_neg
            intst[ncyc] = spectrometer.detect(
                d_eq
                @ d_delta[ncyc]
                @ p90[3]
                @ cpmg2
                @ palmer
                @ cpmg1
                @ p90[0]
                @ d_delta[ncyc]
                @ start,
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return metadata == 0


def register() -> None:
    creators = Creators(
        config_creator=Cpmg15NTr0013Config,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg15NTr0013Sequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
