from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import CestSettings
from chemex.configuration.experiment import ExperimentConfig
from chemex.configuration.experiment import ToBeFitted
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators
from chemex.experiments.factories import factories
from chemex.filterers import CestFilterer
from chemex.nmr.constants import get_multiplet
from chemex.nmr.liouvillian import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters import CestPlotter
from chemex.printers.data import CestPrinter

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "cest_15n_cw"


class Cest15NCwSettings(CestSettings):
    name: Literal["cest_15n_cw"]
    time_t1: float
    carrier: float
    carrier_dec: float
    b1_frq_dec: float
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cest15NCwConfig(ExperimentConfig[Cest15NCwSettings, CestDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[
                "r2_i",
                f"r2_s_{state}",
                f"r1_i_{state}",
                f"r1_s_{state}",
                f"r2mq_is_{state}",
                f"etaxy_i_{state}",
                f"etaz_i_{state}",
            ],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: Cest15NCwConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsxyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res

    spectrometer.b1_s = settings.b1_frq_dec
    spectrometer.carrier_s = settings.carrier_dec
    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        liouvillian.jeff_i = get_multiplet("", "n")

    return spectrometer


@dataclass
class Cest15NCwSequence:
    settings: Cest15NCwSettings

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return np.abs(metadata) > 1e4

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        offsets = data.metadata

        start = spectrometer.get_start_magnetization(terms=["iz"], atom="n")

        intensities = {}

        for offset in set(offsets):

            intensities[offset] = start

            if self.is_reference(offset):
                continue

            spectrometer.offset_i = offset

            intensities[offset] = (
                spectrometer.pulse_is(self.settings.time_t1, 0.0, 0.0)
                @ intensities[offset]
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets]
        )


def register():
    creators = Creators(
        config_creator=Cest15NCwConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest15NCwSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
