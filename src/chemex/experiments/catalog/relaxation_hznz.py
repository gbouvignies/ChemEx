from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import RelaxationSettings
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
from chemex.typing import Array

EXPERIMENT_NAME = "relaxation_hznz"


class RelaxationHzNzSettings(RelaxationSettings):
    """Settings for longitudinal 1H-15N two-spin order (HzNz) relaxation experiment."""

    name: Literal["relaxation_hznz"]

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
        return f"[2izsz{self.suffix_detect}]"


class RelaxationHzNzConfig(
    ExperimentConfiguration[
        RelaxationHzNzSettings,
        ConditionsWithValidations,
        RelaxationDataSettings,
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=[f"r1a_is_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: RelaxationHzNzConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="izsz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.detection = settings.detection

    return spectrometer


class RelaxationHzNzSequence:
    """Sequence for longitudinal 1H-15N two-spin order (HzNz) relaxation experiment."""

    def __init__(self, settings: RelaxationHzNzSettings) -> None:
        self.settings = settings

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        times = data.metadata

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        # Return profile
        delays = spectrometer.delays(0.25 * np.array(times))
        p180_i = spectrometer.perfect180_i[0]
        p180_s = spectrometer.perfect180_s[0]
        return np.array(
            [
                spectrometer.detect(
                    delay @ p180_s @ delay @ p180_i @ delay @ p180_s @ delay @ start,
                )
                for delay in delays
            ],
        )

    def is_reference(self, metadata: Array) -> Array:
        return np.full_like(metadata, fill_value=False, dtype=np.bool_)


def register() -> None:
    creators = Creators(
        config_creator=RelaxationHzNzConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=RelaxationHzNzSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=RelaxationPrinter,
        plotter_creator=RelaxationPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
