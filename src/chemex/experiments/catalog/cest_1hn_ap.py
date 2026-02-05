from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field, computed_field

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, CestSettings
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

EXPERIMENT_NAME = "cest_1hn_ap"

OFFSET_REF = 1e4


class Cest1HnApSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for anti-phase 1H-15N CEST experiment."""

    name: Literal["cest_1hn_ap"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="1H carrier position in Hz")
    cs_evolution_prior: bool = True

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization term (anti-phase)."""
        return [f"2izsz{self.suffix_start}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection operator (anti-phase)."""
        return f"[2izsz{self.suffix_detect}]"


class Cest1HnApConfig(
    ExperimentConfiguration[
        Cest1HnApSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}", f"r1_s_{state}", f"etaxy_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}", f"khh_{state}"],
        )


def build_spectrometer(
    config: Cest1HnApConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier

    # Set the B1 inhomogeneity distribution
    distribution = settings.get_b1_distribution()
    liouvillian.set_b1_i_distribution(distribution)

    spectrometer.detection = settings.detection

    return spectrometer


class Cest1HnApSequence:
    """Sequence for anti-phase 1H-15N CEST experiment."""

    def __init__(self, settings: Cest1HnApSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(self.settings.start_terms)
        intensities: dict[float, Array] = {}

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


def register() -> None:
    creators = Creators(
        config_creator=Cest1HnApConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest1HnApSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
