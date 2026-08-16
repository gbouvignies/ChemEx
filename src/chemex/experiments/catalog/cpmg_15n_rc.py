from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettingsEvenNcycs
from chemex.configuration.types import Delay, Frequency, PulseWidth
from chemex.containers.data import Data
from chemex.experiments.experiment_types import (
    ProfileCalculation,
    cpmg_type,
    register_experiment_type,
)
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.typing import Array

EXPERIMENT_NAME = "cpmg_15n_rc"

# Phases of the refocusing pulses within one repetition of the CPMG loops, in
# the order they are applied. 0/1/2/3 denote x/y/-x/-y. The anti-phase block
# uses the 1102 cycle and the in-phase block the same cycle rotated by 90
# degrees, so that in both blocks half of the pulses are perpendicular to the
# magnetization.
_ANTIPHASE_CYCLE = (1, 1, 0, 2)
_ANTIPHASE_REMAINDER = (1, 1)
_INPHASE_CYCLE = (0, 0, 1, 3)
_INPHASE_REMAINDER = (0, 0)


class Cpmg15NRcSettings(CpmgSettingsEvenNcycs):
    """Settings for the relaxation-compensated 15N CPMG experiment.

    Two CPMG blocks of equal duration are separated by a P-element that
    interconverts anti-phase and in-phase 15N magnetization. The magnetization
    is anti-phase during the first block and in-phase during the second one, so
    the measured rate is the average of the two irrespective of the CPMG
    frequency.

    The compensation delay equalizes the time spent along z during the
    refocusing pulses across the whole series, which is why the largest ncyc is
    needed.

    References:
        D. Long, M. Liu and D. Yang. J. Am. Chem. Soc. 130, 2432-2433 (2008).
        T. Yuwen and L. E. Kay. J. Biomol. NMR 73, 641-650 (2019).

    """

    name: Literal["cpmg_15n_rc"]
    time_t2: Delay = Field(description="Total CPMG relaxation delay in seconds")
    carrier: Frequency = Field(description="15N carrier position in Hz")
    pw90: PulseWidth = Field(
        description="15N 90-degree pulse width in seconds, at the CPMG power level",
    )
    ncyc_max: int = Field(
        gt=0,
        description="Largest ncyc of the series, used for the R1 compensation delay",
    )
    taub: float = 2.68e-3

    @property
    def pw180(self) -> float:
        """Width of the CPMG refocusing pulse."""
        return 2.0 * self.pw90

    @computed_field
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms (anti-phase, from the refocused INEPT)."""
        return self.get_start_terms("2izsz")

    @computed_field
    @property
    def detection(self) -> str:
        """The last 15N 90-degree pulse stores in-phase magnetization along z."""
        return self.get_detection_expression("[iz]")


class Cpmg15NRcConfig(
    ExperimentConfiguration[
        Cpmg15NRcSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.primary_state

        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg15NRcConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    # The anti-phase block and the J evolution during the CPMG trains require
    # the two-spin longitudinal basis.
    basis = Basis(type="ixyzsz", spin_system="nh", model=config.model)
    spectrometer = Spectrometer.from_spin_system(spin_system, basis, conditions)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


def _echo_train(echoes: list[Array], phases: tuple[int, ...]) -> Array:
    """Chain spin echoes applied in the given order of refocusing phases."""
    propagator = echoes[phases[0]]
    for phase in phases[1:]:
        propagator = echoes[phase] @ propagator
    return propagator


class Cpmg15NRcSequence:
    """Sequence for the relaxation-compensated 15N CPMG experiment."""

    def __init__(self, settings: Cpmg15NRcSettings) -> None:
        self.settings = settings

    @staticmethod
    def even_ncyc(ncyc: float) -> int:
        """Round ncyc down to the even value used by the pulse sequence.

        The number of refocusing pulses per block is set by an integer loop
        counter equal to ncyc / 2, and the CPMG frequency is recomputed from it
        before the inter-pulse delay is derived. Valid CPMG frequencies are
        multiples of 2 / time_t2, for which this rounding is a no-op.
        """
        return 2 * (int(ncyc) // 2)

    def _get_delays(
        self,
        ncycs: Array,
    ) -> tuple[dict[float, float], dict[float, float], list[float]]:
        settings = self.settings

        tau_cps = {
            ncyc: settings.time_t2 / (4.0 * self.even_ncyc(ncyc))
            - 0.375 * settings.pw180
            for ncyc in ncycs
            if self.even_ncyc(ncyc) > 0
        }
        deltas = {
            ncyc: 0.25 * (settings.ncyc_max - self.even_ncyc(ncyc)) * settings.pw180
            for ncyc in ncycs
        }
        delays = [settings.taub, *tau_cps.values(), *deltas.values()]

        return tau_cps, deltas, delays

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        ncycs = data.metadata
        settings = self.settings

        # Calculation of the propagators corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
        d_taub = delays[settings.taub]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_comp = {ncyc: delays[delay] for ncyc, delay in deltas.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180_ix = spectrometer.perfect180_i[0]
        p180_sx = spectrometer.perfect180_s[0]

        # Calculating the P-element, which lets J evolve for 1 / (2 J) and
        # converts the anti-phase magnetization of the first block into the
        # in-phase magnetization of the second one
        p_element = d_taub @ p180_sx @ p180_ix @ d_taub

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(settings.start_terms)

        # Calculating the intensities as a function of ncyc
        intst: dict[float, float] = {}

        for ncyc in set(ncycs):
            d_r1 = d_comp[ncyc]
            magnetization = p90[0] @ d_r1 @ start

            if self.even_ncyc(ncyc) > 0:
                echoes = [d_cp[ncyc] @ p180[phase] @ d_cp[ncyc] for phase in range(4)]
                counter1, counter2 = divmod(self.even_ncyc(ncyc) // 2, 2)
                cpmg1 = matrix_power(
                    _echo_train(echoes, _ANTIPHASE_REMAINDER),
                    counter2,
                ) @ matrix_power(_echo_train(echoes, _ANTIPHASE_CYCLE), counter1)
                cpmg2 = matrix_power(
                    _echo_train(echoes, _INPHASE_REMAINDER),
                    counter2,
                ) @ matrix_power(_echo_train(echoes, _INPHASE_CYCLE), counter1)
                magnetization = cpmg2 @ p_element @ cpmg1 @ magnetization
            else:
                magnetization = p_element @ magnetization

            intst[ncyc] = spectrometer.detect(d_r1 @ p90[1] @ magnetization)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return metadata == 0


def create_profile_calculation(
    config: Cpmg15NRcConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    return ProfileCalculation(
        spectrometer=build_spectrometer(config, spin_system),
        pulse_sequence=Cpmg15NRcSequence(config.experiment),
    )


EXPERIMENT_TYPE = cpmg_type(
    name=EXPERIMENT_NAME,
    config_type=Cpmg15NRcConfig,
)(create_profile_calculation)


def register() -> None:
    register_experiment_type(EXPERIMENT_TYPE)
