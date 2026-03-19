from __future__ import annotations

import math
from typing import Any, Literal, TypeVar

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

from chemex.configuration.b1_config import parse_distribution_config
from chemex.configuration.types import OptionalB1Field, PulseWidth
from chemex.configuration.utils import key_to_lower
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import DistributionConfig

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"
    model_name: str = Field(default="", exclude=True)

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
    b1_distribution: DistributionConfig | None = Field(
        default=None,
        description="B1 distribution configuration (type, scale, res, etc.)",
    )
    b1_inh_scale: float | None = Field(
        default=None,
        ge=0.0,
        description="Legacy Gaussian B1 inhomogeneity scale",
    )
    b1_inh_res: int | None = Field(
        default=None,
        ge=1,
        description="Legacy Gaussian B1 inhomogeneity resolution",
    )

    @field_validator("b1_distribution", mode="before")
    @classmethod
    def _parse_b1_distribution(cls, value: object) -> DistributionConfig | None:
        if value is None:
            return None
        if isinstance(value, dict) and "type" not in value:
            value = {"type": "gaussian", **value}
        return parse_distribution_config(value)

    @model_validator(mode="before")
    @classmethod
    def _normalize_legacy_b1_fields(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data

        normalized = dict(data)
        has_legacy = (
            normalized.get("b1_inh_scale") is not None
            or normalized.get("b1_inh_res") is not None
        )
        if not has_legacy:
            return normalized
        if normalized.get("b1_distribution") is not None:
            msg = (
                "Use either 'b1_distribution' or the legacy "
                "'b1_inh_scale'/'b1_inh_res' fields, not both"
            )
            raise ValueError(msg)

        scale = (
            normalized["b1_inh_scale"]
            if normalized.get("b1_inh_scale") is not None
            else 0.1
        )
        if math.isinf(scale):
            normalized["b1_distribution"] = {"type": "dephasing"}
            return normalized

        res = (
            normalized["b1_inh_res"]
            if normalized.get("b1_inh_res") is not None
            else 11
        )
        normalized["b1_distribution"] = {
            "type": "gaussian",
            "scale": scale,
            "res": res,
        }
        return normalized

    def get_b1_nominal(self) -> float:
        """Nominal B1 value for distribution generation."""
        if self.b1_frq is not None:
            return float(self.b1_frq)
        if self.pw90 is not None:
            return 1.0 / (4.0 * float(self.pw90))
        msg = "Either 'b1_frq' or 'pw90' must be provided"
        raise ValueError(msg)

    def get_b1_distribution(self) -> Distribution:
        """Generate B1 distribution using typed plugin configuration."""
        nominal = self.get_b1_nominal()

        if self.b1_distribution is None:
            # Single-point distribution (no inhomogeneity)
            return Distribution(
                values=np.array([nominal]),
                weights=np.array([1.0]),
            )

        return self.b1_distribution.get_distribution(nominal)


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

    @computed_field
    @property
    def suffix_start(self) -> str:
        """Suffix for starting terms.

        Returns '_{observed_state}' if cs_evolution_prior is True (non-equilibrium
        initial condition with only the observed state populated), else ''
        (equilibrium initial condition with all states populated).
        """
        return f"_{self.observed_state}" if self.cs_evolution_prior else ""

    @computed_field
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
