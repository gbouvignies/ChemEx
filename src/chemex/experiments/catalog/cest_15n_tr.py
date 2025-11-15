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
from chemex.nmr.constants import get_multiplet
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.printers.data import CestPrinter
from chemex.typing import Array

EXPERIMENT_NAME = "cest_15n_tr"

OFFSET_REF = 1e4


class Cest15NTrSettings(CestSettings, B1InhomogeneityMixin):
    """Settings for TROSY-based 15N CEST experiment."""

    name: Literal["cest_15n_tr"]
    time_t1: float = Field(description="Length of the CEST block in seconds")
    carrier: Frequency = Field(description="15N carrier position in Hz")
    antitrosy: bool = False

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Start from the TROSY or ANTI-TROSY component.

        TROSY: (2IzSz + Iz) / 2
        ANTITROSY: (2IzSz - Iz) / 2.
        """
        suffix = self.suffix_start
        if not self.antitrosy:
            return [f"2izsz{suffix}", f"-iz{suffix}"]
        return [f"2izsz{suffix}", f"iz{suffix}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection operator for TROSY or ANTI-TROSY component."""
        suffix = self.suffix_detect
        if self.antitrosy:
            return f"[2izsz{suffix}] + [iz{suffix}]"
        return f"[2izsz{suffix}] - [iz{suffix}]"


class Cest15NTrConfig(
    ExperimentConfiguration[
        Cest15NTrSettings, ConditionsWithValidations, CestDataSettings
    ],
):
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
    config: Cest15NTrConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier

    # Set the B1 inhomogeneity distribution
    distribution = settings.get_b1_distribution()
    liouvillian.set_b1_i_distribution(distribution)

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        liouvillian.jeff_i = get_multiplet("", "n")

    return spectrometer


class Cest15NTrSequence:
    """Sequence for TROSY-based 15N CEST experiment."""

    def __init__(self, settings: Cest15NTrSettings) -> None:
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
            ).real.astype(float)

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets],
        )


def register() -> None:
    creators = Creators(
        config_creator=Cest15NTrConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest15NTrSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
