from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
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

EXPERIMENT_NAME = "cpmg_13c_ip"


class Cpmg13CIpSettings(CpmgSettings):
    """Settings for in-phase 13C CPMG relaxation dispersion experiment."""

    name: Literal["cpmg_13c_ip"]
    time_t2: Delay = Field(description="Total CPMG relaxation delay in seconds")
    carrier: Frequency = Field(description="13C carrier position in Hz")
    pw90: PulseWidth = Field(description="90-degree pulse width in seconds")
    time_equil: Delay = 0.0

    @computed_field  # type: ignore[misc]
    @property
    def t_neg(self) -> float:
        """Negative time delay for CPMG element."""
        return -2.0 * self.pw90 / np.pi

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms for the experiment."""
        return [f"iz{self.suffix_start}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection operator for the experiment."""
        return f"[iz{self.suffix_detect}]"


class Cpmg13CIpConfig(
    ExperimentConfiguration[
        Cpmg13CIpSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg13CIpConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


class Cpmg13CIpSequence:
    """Sequence for in-phase 13C CPMG relaxation dispersion experiment."""

    def __init__(self, settings: Cpmg13CIpSettings) -> None:
        self.settings = settings

    def _get_delays(self, ncycs: Array) -> tuple[dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        # An ncyc of -1 corresponds to the experiment without the CPMG
        # refocusisng pulses
        tau_cps[-1.0] = 0.5 * self.settings.time_t2
        delays = [self.settings.t_neg, self.settings.time_equil, *tau_cps.values()]
        return tau_cps, delays

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
        d_neg = delays[self.settings.t_neg]
        d_eq = delays[self.settings.time_equil]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        # Calculating the instensities as a function of ncyc
        part1 = d_neg @ p90[0] @ start
        part2 = d_eq @ p90[0] @ d_neg
        intst = {
            0.0: spectrometer.detect(d_eq @ p90[0] @ p180pmx @ p90[0] @ start),
            -1.0: spectrometer.detect(part2 @ d_cp[-1] @ p180pmx @ d_cp[-1] @ part1),
        }
        for ncyc in set(ncycs) - {0.0, -1.0}:
            echo = d_cp[ncyc] @ p180[1] @ d_cp[ncyc]
            cpmg = matrix_power(echo, int(ncyc))
            intst[ncyc] = spectrometer.detect(part2 @ cpmg @ p180pmx @ cpmg @ part1)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return metadata == 0


def register() -> None:
    creators = Creators(
        config_creator=Cpmg13CIpConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg13CIpSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
