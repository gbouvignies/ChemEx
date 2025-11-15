from __future__ import annotations

import math
from typing import Any, Literal, TypeVar

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, computed_field, model_validator

from chemex.configuration.types import OptionalB1Field, PulseWidth
from chemex.configuration.utils import key_to_lower
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import (
    get_b1_distribution as get_distribution_from_registry,
)

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"

    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class B1InhomogeneityMixin(BaseModel):
    """Unified B1 inhomogeneity configuration for all experiments.

    This mixin should be inherited by experiment settings that support
    B1 field distributions for modeling inhomogeneity effects.

    Usage:
      [experiment]
      pw90 = <hardware pulse width, seconds>
      b1_frq = <effective B1, Hz or kHz>

    [experiment.b1_distribution]
    type = "gaussian"  # or "skewed", "custom", etc.
    scale = 0.1        # Fractional std. dev. (e.g., 0.1 for ~10%)
    res = 11           # Number of grid points
    # ...other plugin-specific parameters...
    """

    model_config = ConfigDict(
        extra="allow",  # Mixins allow additional fields
        validate_assignment=True,  # Validate on attribute assignment
    )

    pw90: PulseWidth | None = Field(
        default=None,
        description="90-degree pulse width (seconds) for hardware B1 calibration",
        examples=[45e-6, 96e-6],
    )
    b1_frq: OptionalB1Field = Field(
        default=None,
        description="B1 radio-frequency field strength (Hz or kHz)",
        examples=[25.0, 700.0],
    )
    b1_distribution: dict[str, Any] | None = Field(
        default=None,
        description="B1 distribution configuration (type, scale, res, etc.)",
    )

    def get_b1_nominal(self) -> float:
        """Nominal B1 value for distribution generation."""
        if self.b1_frq is not None:
            return float(self.b1_frq)
        if self.pw90 is not None:
            return 1.0 / (4.0 * float(self.pw90))
        msg = "Either 'b1_frq' or 'pw90' must be provided"
        raise ValueError(msg)

    def get_b1_distribution(self) -> Distribution:
        """Generate B1 distribution using the plugin registry."""
        nominal = self.get_b1_nominal()

        if not self.b1_distribution:
            # Single-point distribution (no inhomogeneity)
            return Distribution(
                values=np.array([nominal]),
                weights=np.array([1.0]),
            )

        dist_type = self.b1_distribution.get("type", "gaussian")
        params = dict(self.b1_distribution)
        params.pop("type", None)

        return get_distribution_from_registry(dist_type, nominal, **params)


class RelaxationSettings(ExperimentSettings):
    """Base settings for relaxation-based experiments (CEST, DCEST, CPMG, etc).

    Attributes:
        cs_evolution_prior (bool):
            Controls the initial condition (equilibrium vs non-equilibrium).
            If False (default), the initial magnetization assumes thermal
            equilibrium with all states populated according to their equilibrium
            populations. If True, only the observed state is initially populated
            (non-equilibrium condition), which is appropriate when chemical shift
            evolution occurs before the CEST/CPMG element in the pulse sequence,
            breaking the equilibrium assumption. This can affect the accuracy of
            extracted kinetic parameters, especially for slow exchange rates
            (see Yuwen et al., J. Biomol. NMR 2016, 65:143-156).
        detect_all_states (bool):
            If True, detection will use all states (e.g., '[iz]'),
            otherwise only the observed state (e.g., '[iz_a]').

    Properties:
        suffix_start: Suffix for starting terms, respects cs_evolution_prior.
        suffix_detect: Suffix for detection, respects detect_all_states.
            Returns '' if detect_all_states is True, else '_{observed_state}'.

    """

    cs_evolution_prior: bool = False
    detect_all_states: bool = False

    @computed_field  # type: ignore[misc]
    @property
    def suffix_start(self) -> str:
        """Suffix for starting terms.

        Returns '_{observed_state}' if cs_evolution_prior is True (non-equilibrium
        initial condition with only the observed state populated), else ''
        (equilibrium initial condition with all states populated).
        """
        return f"_{self.observed_state}" if self.cs_evolution_prior else ""

    @computed_field  # type: ignore[misc]
    @property
    def suffix_detect(self) -> str:
        """Suffix for detection.

        Returns '' if detect_all_states is True, else '_{observed_state}'.
        """
        return "" if self.detect_all_states else f"_{self.observed_state}"


class CpmgSettings(RelaxationSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(RelaxationSettings):
    even_ncycs: bool = True


class CestSettings(RelaxationSettings):
    sw: float = math.inf


class MFCestSettings(RelaxationSettings):
    sw: float
