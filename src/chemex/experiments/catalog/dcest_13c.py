from __future__ import annotations

from typing import Literal, Self

import numpy as np
from numpy.linalg import matrix_power
from pydantic import AliasChoices, Field, computed_field, model_validator

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import B1InhomogeneityMixin, MFCestSettings
from chemex.configuration.types import B1Field, ChemicalShift, Delay
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

EXPERIMENT_NAME = "dcest_13c"

OFFSET_REF = 1e4


class DCest13CSettings(MFCestSettings, B1InhomogeneityMixin):
    """D-CEST 13C experiment settings with DANTE pulse train."""

    name: Literal["dcest_13c"]
    time_t1: Delay = Field(description="CEST relaxation delay (seconds)")
    time_equil: Delay = Field(
        default=0.0,
        description="Equilibration delay at end of CEST period (seconds)",
    )
    carrier: ChemicalShift = Field(description="13C carrier position during CEST (ppm)")
    # Effective B1 for the equivalent continuous-wave CEST irradiation.
    # Accept both 'b1_eff' (preferred) and 'b1_frq' (alias) in configs.
    b1_eff: B1Field = Field(
        ...,
        validation_alias=AliasChoices("b1_eff", "b1_frq"),
        description="Effective B1 field for DANTE excitation (Hz)",
    )

    @model_validator(mode="after")
    def _validate_dcest_fields(self) -> Self:
        """Validate D-CEST requires pw90 for hardware B1."""
        if self.pw90 is None:
            msg = (
                "D-CEST experiments require 'pw90' (hardware B1 pulse width). "
                "Both 'pw90' and 'b1_eff' (alias: 'b1_frq') must be specified."
            )
            raise ValueError(msg)
        return self

    # D-CEST convention: the B1 distribution should be centered on pw90.
    def get_b1_nominal(self) -> float:  # type: ignore[override]
        """Get nominal B1 from hardware pw90 for distribution centering."""
        if self.pw90 is None:
            msg = (
                "For D-CEST, 'pw90' must be specified to define the nominal B1 "
                "used for the B1 distribution."
            )
            raise ValueError(msg)
        return 1.0 / (4.0 * float(self.pw90))

    @computed_field  # type: ignore[misc]
    @property
    def pw_dante(self) -> float:
        """Duration of each DANTE pulse (seconds)."""
        if self.pw90 is None:
            msg = "pw90 must be specified for DANTE pulse calculation"
            raise ValueError(msg)
        return 4.0 * float(self.pw90) * float(self.b1_eff) / float(self.sw)

    @computed_field  # type: ignore[misc]
    @property
    def tau_dante(self) -> float:
        """Inter-pulse delay in DANTE train (seconds)."""
        return 1.0 / self.sw - self.pw_dante

    @computed_field  # type: ignore[misc]
    @property
    def ncyc_dante(self) -> int:
        """Number of DANTE pulses in the train."""
        return int(self.time_t1 * self.sw + 0.1)

    @computed_field  # type: ignore[misc]
    @property
    def start_terms(self) -> list[str]:
        """Starting magnetization terms."""
        return [f"iz{self.suffix_start}"]

    @computed_field  # type: ignore[misc]
    @property
    def detection(self) -> str:
        """Detection operator for the experiment."""
        return f"[iz{self.suffix_detect}]"


class DCest13CConfig(
    ExperimentConfiguration[
        DCest13CSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(config: DCest13CConfig, spin_system: SpinSystem) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="ch")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    # Set the B1 inhomogeneity distribution
    distribution = settings.get_b1_distribution()
    liouvillian.set_b1_i_distribution(distribution)

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        symbol = spin_system.symbols["i"]
        atom = spin_system.atoms["i"]
        liouvillian.jeff_i = get_multiplet(symbol, atom.name)

    return spectrometer


class DCest13CSequence:
    """Sequence for D-CEST 13C experiment."""

    def __init__(self, settings: DCest13CSettings) -> None:
        self.settings = settings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(terms=self.settings.start_terms)

        d_eq = (
            spectrometer.delays(self.settings.time_equil)
            if self.settings.time_equil > 0
            else spectrometer.identity
        )

        intensities: dict[float, Array] = {}

        for offset in set(offsets):
            if self.is_reference(offset):
                intensities[offset] = d_eq @ start
                continue

            spectrometer.offset_i = offset

            p_delay = spectrometer.delays(self.settings.tau_dante)
            p_pulse = spectrometer.pulse_i(self.settings.pw_dante, 0.0)

            intensities[offset] = (
                d_eq @ matrix_power(p_delay @ p_pulse, self.settings.ncyc_dante) @ start
            )

        return np.array(
            [spectrometer.detect(intensities[offset]) for offset in offsets]
        )


def register() -> None:
    creators = Creators(
        config_creator=DCest13CConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=DCest13CSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)
