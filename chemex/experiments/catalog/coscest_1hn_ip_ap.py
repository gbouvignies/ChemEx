from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.linalg import matrix_power
from numpy.typing import NDArray

from chemex.configuration.data import CestDataSettingsNoRef
from chemex.configuration.experiment import CestSettings
from chemex.configuration.experiment import ExperimentConfig
from chemex.configuration.experiment import ToBeFitted
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators
from chemex.experiments.factories import factories
from chemex.filterers import CestFilterer
from chemex.nmr.liouvillian import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters import CestPlotter
from chemex.printers.data import CestPrinter

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "coscest_1hn_ip_ap"


class CosCest1HnIpApSettings(CestSettings):
    name: Literal["coscest_1hn_ip_ap"]
    time_t1: float
    carrier: float
    sw: float
    cos_n: int
    cos_res: int = 10
    d1: float
    taua: float = 2.38e-3
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[2izsz_{self.observed_state}]"


class CosCest1HnIpApConfig(
    ExperimentConfig[CosCest1HnIpApSettings, CestDataSettingsNoRef]
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
    config: CosCest1HnIpApConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz_eq", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class CosCest1HnIpApSequence:
    settings: CosCest1HnIpApSettings

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return np.abs(metadata) > 1e4

    def _calc_cosine_shape(self, spectrometer: Spectrometer) -> np.ndarray:
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
                n_left * dt, amplitudes[:n_left], phases[:n_left]
            )
            pulse = pulse_left[phase_left] @ pulse

        return pulse

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        offsets = data.metadata
        spectrometer.offset_i = 0.0

        d_d1, d_taua, d_t1 = spectrometer.delays(
            [self.settings.d1, self.settings.taua, self.settings.time_t1]
        )

        pp90_i = spectrometer.perfect90_i
        pp180_isx = spectrometer.perfect180_i[0] @ spectrometer.perfect180_s[0]

        start = d_d1 @ spectrometer.get_start_magnetization(terms=["ie"])
        start = spectrometer.keep(start, components=["ie", "iz"])

        intensities = {}
        for offset in set(offsets):
            spectrometer.offset_i = offset
            if self.is_reference(offset):
                inept = pp90_i[3] @ d_taua @ pp180_isx @ d_taua @ pp90_i[0]
                cest = inept @ d_t1
            else:
                cest = self._calc_cosine_shape(spectrometer)
            intensities[offset] = spectrometer.detect(cest @ start)
        return np.array([intensities[offset] for offset in offsets])


def register():
    creators = Creators(
        config_creator=CosCest1HnIpApConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=CosCest1HnIpApSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
