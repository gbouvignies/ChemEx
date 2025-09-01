from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import ExperimentSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.relaxation import RelaxationPlotter
from chemex.printers.data import RelaxationPrinter
from chemex.typing import ArrayBool, ArrayFloat

EXPERIMENT_NAME = "wip.relaxation_15n_r1rho_eig"


class Relaxation15NR1RhoSettings(ExperimentSettings):
    name: Literal["wip.relaxation_15n_r1rho_eig"]

    carrier: float
    b1_frq: float
    b1_inh_scale: float = np.inf
    b1_inh_res: int = 11
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Relaxation15NR1RhoConfig(
    ExperimentConfiguration[
        Relaxation15NR1RhoSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r2_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: Relaxation15NR1RhoConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res

    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Relaxation15NR1RhoSequence:
    settings: Relaxation15NR1RhoSettings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        times = data.metadata

        r1rho = spectrometer.liouvillian.calculate_r1rho()

        # Return profile
        return np.exp(-r1rho * times)

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=Relaxation15NR1RhoConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Relaxation15NR1RhoSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=RelaxationPrinter,
        plotter_creator=RelaxationPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
