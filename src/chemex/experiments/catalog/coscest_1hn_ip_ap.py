from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettingsNoRef
from chemex.configuration.experiment import B1InhomogeneityMixin, MFCestSettings
from chemex.configuration.types import Frequency
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import CestFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.printers.data import CestPrinter
from chemex.typing import Array

EXPERIMENT_NAME = "coscest_1hn_ip_ap"

OFFSET_REF = 1e4


class CosCest1HnIpApSettings(MFCestSettings, B1InhomogeneityMixin):
    """Settings for cosine-modulated 1H-15N in-phase/anti-phase CEST experiment."""

    name: Literal["coscest_1hn_ip_ap"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="1H carrier position in Hz")
    cos_n: int = Field(description="Number of cosine cycles")
    cos_res: int = 10
    d1: float = Field(description="Relaxation delay in seconds")
    taua: float = 2.38e-3

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection operator (anti-phase)."""
        return f"[2izsz{self.suffix_detect}]"


class CosCest1HnIpApConfig(
    ExperimentConfiguration[
        CosCest1HnIpApSettings,
        ConditionsWithValidations,
        CestDataSettingsNoRef,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[
                "r2_i",
                f"r1_i_{state}",
                f"r1_s_{state}",
                f"etaxy_i_{state}",
                f"etaz_i_{state}",
            ],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: CosCest1HnIpApConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz_eq", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier

    # Set the B1 inhomogeneity distribution
    distribution = settings.get_b1_distribution()
    liouvillian.set_b1_i_distribution(distribution)

    spectrometer.detection = settings.detection

    return spectrometer


class CosCest1HnIpApSequence:
    """Sequence for cosine-modulated 1H-15N in-phase/anti-phase CEST experiment."""

    def __init__(self, settings: CosCest1HnIpApSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def _calc_cosine_shape(self, spectrometer: Spectrometer) -> Array:
        time_t1 = self.settings.time_t1
        sw = self.settings.sw
        cos_n = self.settings.cos_n
        cos_res = self.settings.cos_res

        dt = 1.0 / (cos_res * sw)
        n_periods = int(time_t1 * sw)
        n_left = int((time_t1 * sw - n_periods) * cos_res)
        double_periods, extra_period = divmod(n_periods, 2)
        phase1 = 2 if cos_n % 2 == 0 else 0
        phase_left = 0 if extra_period else phase1

        grid = np.linspace(-np.pi, np.pi, cos_res, endpoint=False)

        n_values = (np.arange(cos_n) - 0.5 * (cos_n - 1)).reshape(-1, 1)
        amplitudes = np.cos(n_values * grid).sum(axis=0)
        phases = np.zeros(cos_res)

        base_pulse = spectrometer.shaped_pulse_i(cos_res * dt, amplitudes, phases)

        pulse = matrix_power(base_pulse[0] @ base_pulse[phase1], double_periods)

        if extra_period:
            pulse = base_pulse[phase1] @ pulse

        if n_left:
            pulse_left = spectrometer.shaped_pulse_i(
                n_left * dt,
                amplitudes[:n_left],
                phases[:n_left],
            )
            pulse = pulse_left[phase_left] @ pulse

        return pulse

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata
        spectrometer.offset_i = 0.0

        d_d1, d_taua, d_t1 = spectrometer.delays(
            [self.settings.d1, self.settings.taua, self.settings.time_t1],
        )

        pp90_i = spectrometer.perfect90_i
        pp180_isx = spectrometer.perfect180_i[0] @ spectrometer.perfect180_s[0]

        start = d_d1 @ spectrometer.get_start_magnetization(terms=["ie"])
        start = spectrometer.keep(start, components=["ie", "iz"])

        intensities: dict[float, float] = {}

        for offset in set(offsets):
            spectrometer.offset_i = offset
            if self.is_reference(offset):
                inept = pp90_i[3] @ d_taua @ pp180_isx @ d_taua @ pp90_i[0]
                cest = inept @ d_t1
            else:
                cest = self._calc_cosine_shape(spectrometer)
            intensities[offset] = spectrometer.detect(cest @ start)

        return np.array([intensities[offset] for offset in offsets])


def register() -> None:
    creators = Creators(
        config_creator=CosCest1HnIpApConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CosCest1HnIpApSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
