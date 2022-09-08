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


EXPERIMENT_NAME = "cest_1hn_ip_ap"


class Cest1HnIpApSettings(CestSettings):
    name: Literal["cest_1hn_ip_ap"]
    time_t1: float
    carrier: float
    d1: float
    taua: float = 2.38e-3
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11
    eta_block: int = 0
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[2izsz_{self.observed_state}]"

    @property
    def taud(self) -> float:
        return max(self.d1 - self.time_t1, 0.0) if self.eta_block > 0 else self.d1

    @property
    def taue(self) -> float:
        return 0.5 * self.time_t1 / self.eta_block


class Cest1HnIpApConfig(ExperimentConfig[Cest1HnIpApSettings, CestDataSettingsNoRef]):
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
    config: Cest1HnIpApConfig, spin_system: SpinSystem
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
class Cest1HnIpApSequence:
    settings: Cest1HnIpApSettings

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return np.abs(metadata) > 1e4

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        offsets = data.metadata
        spectrometer.offset_i = 0.0

        d_taud, d_taua = spectrometer.delays([self.settings.taud, self.settings.taua])

        pp90_i = spectrometer.perfect90_i
        pp180_sx = spectrometer.perfect180_s[0]
        pp180_isx = spectrometer.perfect180_i[0] @ spectrometer.perfect180_s[0]

        start = d_taud @ spectrometer.get_start_magnetization(terms=["ie"])
        start = spectrometer.keep(start, components=["ie", "iz"])

        intensities = {}
        for offset in set(offsets):
            spectrometer.offset_i = offset
            if self.settings.eta_block > 0:
                d_2taue = spectrometer.delays(2.0 * self.settings.taue)
                p_taue = spectrometer.pulse_i(self.settings.taue, 0.0)
                cest_block = p_taue @ pp180_sx @ d_2taue @ pp180_sx @ p_taue
                cest = matrix_power(cest_block, self.settings.eta_block)
            else:
                cest = spectrometer.pulse_i(self.settings.time_t1, 0.0)
            if self.is_reference(offset):
                inept = pp90_i[3] @ d_taua @ pp180_isx @ d_taua @ pp90_i[0]
                cest = inept @ cest
            intensities[offset] = spectrometer.detect(cest @ start)
        return np.array([intensities[offset] for offset in offsets])


def register():
    creators = Creators(
        config_creator=Cest1HnIpApConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest1HnIpApSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
