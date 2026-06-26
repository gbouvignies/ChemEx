from __future__ import annotations

import math
from string import ascii_lowercase
from typing import Annotated, Any, ClassVar, TypeVar

import numpy as np
from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)

from chemex.configuration.b1_config import parse_distribution_config
from chemex.configuration.types import OptionalB1Field, PulseWidth
from chemex.configuration.utils import key_to_lower
from chemex.models.model import ModelSpec
from chemex.nmr._engine.detection import build_state_detection_expression
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import DistributionConfig

T = TypeVar("T")
_SUPPORTED_STATES = tuple(ascii_lowercase[:6])


def _model_states_from_context(info: ValidationInfo) -> tuple[str, ...]:
    context = info.context
    model = context.get("model") if isinstance(context, dict) else None
    if isinstance(model, ModelSpec):
        return tuple(model.states)
    return _SUPPORTED_STATES


def _validate_state(value: object, info: ValidationInfo) -> str:
    field = info.field_name or "state"
    if not isinstance(value, str):
        msg = f"'{field}' must be a state string"
        raise ValueError(msg)  # noqa: TRY004

    state = value.lower()
    available_states = _model_states_from_context(info)
    if state not in available_states:
        available_list = ", ".join(repr(item) for item in available_states)
        msg = (
            f"Unknown {field.replace('_', ' ')}: {state!r}. "
            f"Available states: {available_list}"
        )
        raise ValueError(msg)
    return state


def _validate_state_selection(
    value: object,
    info: ValidationInfo,
) -> str | tuple[str, ...]:
    field = info.field_name or "state"
    if isinstance(value, str):
        return _validate_state(value, info)
    if not isinstance(value, (list, tuple)):
        msg = (
            f"'{field}' must be a state string or a non-empty list "
            "of state strings"
        )
        raise ValueError(msg)  # noqa: TRY004
    if not value:
        msg = f"'{field}' must contain at least one state"
        raise ValueError(msg)
    if not all(isinstance(state, str) for state in value):
        msg = f"Every entry in '{field}' must be a string"
        raise ValueError(msg)

    states = tuple(_validate_state(state, info) for state in value)
    duplicates = sorted(state for state in set(states) if states.count(state) > 1)
    if duplicates:
        duplicate_list = ", ".join(repr(state) for state in duplicates)
        msg = f"'{field}' contains duplicate states: {duplicate_list}"
        raise ValueError(msg)
    return states


def _validate_start_state_selection(
    value: object,
    info: ValidationInfo,
) -> str | tuple[str, ...]:
    if isinstance(value, (list, tuple)) and not value:
        return ()
    return _validate_state_selection(value, info)


_State = Annotated[str, BeforeValidator(_validate_state)]
_StateSelection = Annotated[
    str | tuple[str, ...],
    BeforeValidator(_validate_state_selection),
]
_StartStateSelection = Annotated[
    str | tuple[str, ...],
    BeforeValidator(_validate_start_state_selection),
]


def _legacy_b1_distribution_input(
    scale: float | None,
    res: int | None,
) -> dict[str, object]:
    current_scale = 0.1 if scale is None else scale
    if math.isinf(current_scale):
        return {"type": "dephasing"}
    return {
        "type": "gaussian",
        "scale": current_scale,
        "res": 11 if res is None else res,
    }


def normalize_b1_eff_alias(data: Any) -> Any:
    """Map legacy D-CEST b1_frq input to b1_eff without populating both fields."""
    if not isinstance(data, dict) or "b1_frq" not in data:
        return data

    b1_frq = data["b1_frq"]
    if b1_frq is None:
        return data

    normalized = dict(data)
    normalized.pop("b1_frq")
    if "b1_eff" in normalized and normalized["b1_eff"] != b1_frq:
        msg = "Use either 'b1_eff' or its legacy alias 'b1_frq', not both"
        raise ValueError(msg)
    normalized.setdefault("b1_eff", b1_frq)
    return normalized


class ExperimentSettings(BaseModel):
    observed_state: _State = "a"
    model_name: str = Field(default="", exclude=True)

    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class DetectionSettings(ExperimentSettings):
    """Settings for experiments that detect a final magnetization vector."""

    observed_state: _StateSelection = "a"

    @property
    def observed_states(self) -> tuple[str, ...]:
        """Normalized states contributing to final-magnetization detection."""
        if isinstance(self.observed_state, str):
            return (self.observed_state,)
        return self.observed_state

    @property
    def primary_state(self) -> str:
        """First observed state, used where a single reference state is needed."""
        return self.observed_states[0]

    def get_detection_expression(self, expression: str) -> str:
        """Apply the selected observed states to a detection expression."""
        return build_state_detection_expression(expression, self.observed_states)


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
        extra="forbid",
        validate_assignment=True,  # Validate on attribute assignment
    )

    legacy_b1_inh_scale_default: ClassVar[float | None] = None
    legacy_b1_inh_res_default: ClassVar[int | None] = None

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
        distribution = normalized.get("b1_distribution")
        has_legacy = (
            normalized.get("b1_inh_scale") is not None
            or normalized.get("b1_inh_res") is not None
        )
        if has_legacy and distribution is not None and not isinstance(
            distribution, DistributionConfig
        ):
            msg = (
                "Use either 'b1_distribution' or the legacy "
                "'b1_inh_scale'/'b1_inh_res' fields, not both"
            )
            raise ValueError(msg)

        if distribution is not None:
            return normalized

        scale = normalized.get("b1_inh_scale", cls.legacy_b1_inh_scale_default)
        res = normalized.get("b1_inh_res", cls.legacy_b1_inh_res_default)
        if scale is None and res is None:
            return normalized

        if "b1_inh_scale" not in normalized and scale is not None:
            normalized["b1_inh_scale"] = scale
        if "b1_inh_res" not in normalized and res is not None:
            normalized["b1_inh_res"] = res
        normalized["b1_distribution"] = _legacy_b1_distribution_input(
            scale,
            res,
        )
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


class RelaxationSettings(DetectionSettings):
    """Base settings for relaxation-based experiments (CEST, DCEST, CPMG, etc).

    Attributes:
        start_state:
            Optional state or states used for non-equilibrium starting
            magnetization. If omitted, the experiment-specific default is used;
            this is thermal equilibrium for most experiments. Set it when
            chemical shift evolution or another preparation step occurs before
            the measured relaxation/exchange element and the equilibrium
            assumption no longer matches the pulse sequence.

    """

    start_from_observed_by_default: ClassVar[bool] = False
    start_state: _StartStateSelection | None = None

    @model_validator(mode="before")
    @classmethod
    def _reject_removed_state_flags(
        cls,
        data: Any,
        info: ValidationInfo,
    ) -> Any:
        if not isinstance(data, dict):
            return data

        values = {str(key).lower(): value for key, value in data.items()}
        detect_all_states = values.get("detect_all_states")
        if detect_all_states is True:
            states = ", ".join(
                f'"{state}"' for state in _model_states_from_context(info)
            )
            msg = (
                "'detect_all_states' has been removed; use "
                f"'observed_state = [{states}]'"
            )
            raise ValueError(msg)
        if detect_all_states is False:
            msg = (
                "'detect_all_states' has been removed; remove it to retain "
                "single-state detection through 'observed_state'"
            )
            raise ValueError(msg)
        if "detect_all_states" in values:
            msg = (
                "'detect_all_states' has been removed; configure detection with "
                "'observed_state' instead"
            )
            raise ValueError(msg)

        cs_evolution_prior = values.get("cs_evolution_prior")
        if cs_evolution_prior is True:
            msg = (
                "'cs_evolution_prior' has been removed; set 'start_state' to "
                "the observed state or states"
            )
            raise ValueError(msg)
        if cs_evolution_prior is False:
            replacement = (
                "use 'start_state = []' to retain equilibrium preparation"
                if cls.start_from_observed_by_default
                else "remove it to retain equilibrium preparation"
            )
            msg = f"'cs_evolution_prior' has been removed; {replacement}"
            raise ValueError(msg)
        if "cs_evolution_prior" in values:
            msg = (
                "'cs_evolution_prior' has been removed; configure preparation "
                "with 'start_state' instead"
            )
            raise ValueError(msg)
        return data

    @property
    def start_states(self) -> tuple[str, ...]:
        """Normalized non-equilibrium starting states, or empty for equilibrium."""
        if self.start_state is None:
            if self.start_from_observed_by_default:
                return self.observed_states
            return ()
        if isinstance(self.start_state, str):
            return (self.start_state,)
        return self.start_state

    def get_start_terms(self, *components: str) -> list[str]:
        """Apply starting-state selection to basis components."""
        states = self.start_states
        if not states:
            return list(components)
        return [
            f"{component}_{state}"
            for state in states
            for component in components
        ]


class CpmgSettings(RelaxationSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(RelaxationSettings):
    even_ncycs: bool = True


class CestSettings(RelaxationSettings):
    sw: float = math.inf


class MFCestSettings(RelaxationSettings):
    sw: float
