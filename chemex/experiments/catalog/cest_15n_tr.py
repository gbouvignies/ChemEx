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


EXPERIMENT_NAME = "cest_15n_tr"


class Cest15NTrSettings(CestSettings):
    name: Literal["cest_15n_tr"]
    time_t1: float
    carrier: float
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11
    antitrosy: bool = False
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        if self.antitrosy:
            return f"[2izsz_{self.observed_state}] + [iz_{self.observed_state}]"
        else:
            return f"[2izsz_{self.observed_state}] - [iz_{self.observed_state}]"


class Cest15NTrConfig(ExperimentConfig[Cest15NTrSettings, CestDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state

        to_be_fitted = ToBeFitted(
            rates=["r2_i", f"r1_i_{state}", "r1a_is"],
            model_free=[f"tauc_{state}", f"s2_{state}", "khh"],
        )

        if self.experiment.antitrosy:
            to_be_fitted.rates.extend([f"etaxy_i_{state}", f"etaz_i_{state}"])

        return to_be_fitted


def build_spectrometer(
    config: Cest15NTrConfig, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        liouvillian.jeff_i = get_multiplet("", "n")

    return spectrometer


@dataclass
class Cest15NTrSequence:
    settings: Cest15NTrSettings

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return np.abs(metadata) > 1e4

    def _get_start(self, spectrometer: Spectrometer) -> np.ndarray:
        """
        TROSY: (2IzSz - Iz) / 2
        ANTITROSY: (2IzSz + Iz) / 2
        """
        start = 0.5 * spectrometer.get_start_magnetization(["2izsz", "iz"])
        if not self.settings.antitrosy:
            start -= spectrometer.get_start_magnetization(["iz"])
        return start

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        offsets = data.metadata

        start = self._get_start(spectrometer)

        intensities = {}

        for offset in set(offsets):

            intensities[offset] = start

            if self.is_reference(offset):
                continue

            spectrometer.offset_i = offset

            intensities[offset] = (
                spectrometer.pulse_i(self.settings.time_t1, 0.0) @ intensities[offset]
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets]
        )


def register():
    creators = Creators(
        config_creator=Cest15NTrConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest15NTrSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)
