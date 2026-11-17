from __future__ import annotations

import math
from string import ascii_lowercase
from typing import Annotated, Any, ClassVar, Self, TypeVar

from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    PrivateAttr,
    ValidationInfo,
    computed_field,
    field_validator,
    model_validator,
)

from chemex.configuration.b1_config import (
    B1DistributionConfig,
    parse_distribution_config,
)
from chemex.configuration.types import (
    B1Field,
    Delay,
    OptionalB1Field,
    PositiveFloat,
    PulseWidth,
)
from chemex.configuration.utils import key_to_lower
from chemex.models.model import ModelSpec
from chemex.nmr._engine.detection import build_state_detection_expression
from chemex.nmr.b1 import B1RuntimeModel, HomogeneousB1Model
from chemex.nmr.constants import Distribution

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
        raise ValueError(msg)  # noqa: TRY004 - Pydantic field-validation error

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
        msg = f"'{field}' must be a state string or a non-empty list of state strings"
        raise ValueError(msg)  # noqa: TRY004 - Pydantic field-validation error
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
    scale: object,
    res: object,
) -> tuple[dict[str, object], str]:
    if isinstance(scale, (int, float)) and not isinstance(scale, bool):
        if math.isinf(scale):
            return {"type": "dephasing"}, '{ type = "dephasing" }'
        if scale == 0.0:
            return {"type": "homogeneous"}, '{ type = "homogeneous" }'
    if res == 1 and not isinstance(res, bool):
        return {"type": "homogeneous"}, '{ type = "homogeneous" }'

    current_scale: object = 0.1 if scale is None else scale
    scale_text = (
        f"{current_scale:g}"
        if isinstance(current_scale, (int, float))
        and not isinstance(current_scale, bool)
        else repr(current_scale)
    )
    if res is None:
        return (
            {"type": "gaussian", "scale": current_scale},
            f'{{ type = "gaussian", scale = {scale_text} }}',
        )
    return (
        {"type": "gaussian", "scale": current_scale, "res": res},
        f'{{ type = "gaussian", scale = {scale_text}, res = {res} }}',
    )


def normalize_b1_eff_alias(data: Any) -> Any:
    """Map legacy D-CEST b1_frq input to b1_eff without populating both fields."""
    if not isinstance(data, dict):
        return data

    normalized = key_to_lower(data)
    if "b1_frq" not in normalized:
        return normalized

    b1_frq = normalized["b1_frq"]
    if b1_frq is None:
        return normalized

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
    _legacy_b1_deprecation: tuple[tuple[str, ...], str] | None = PrivateAttr(
        default=None
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
    b1_distribution: B1DistributionConfig | None = Field(
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
    def _parse_b1_distribution(cls, value: object) -> B1DistributionConfig | None:
        if value is None:
            return None
        return parse_distribution_config(value)

    @model_validator(mode="before")
    @classmethod
    def _normalize_legacy_b1_fields(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data

        normalized = key_to_lower(data)
        distribution = normalized.get("b1_distribution")
        has_legacy = (
            normalized.get("b1_inh_scale") is not None
            or normalized.get("b1_inh_res") is not None
        )
        if has_legacy and distribution is not None:
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
        legacy_distribution, _ = _legacy_b1_distribution_input(scale, res)
        normalized["b1_distribution"] = legacy_distribution
        return normalized

    @model_validator(mode="after")
    def _finalize_legacy_b1_fields(self) -> Self:
        legacy_fields = tuple(
            field_name
            for field_name in ("b1_inh_scale", "b1_inh_res")
            if getattr(self, field_name) is not None
        )
        if not legacy_fields:
            return self

        cls = type(self)
        if (
            cls.legacy_b1_inh_scale_default is None
            and cls.legacy_b1_inh_res_default is None
        ):
            _, replacement = _legacy_b1_distribution_input(
                self.b1_inh_scale,
                self.b1_inh_res,
            )
            self._legacy_b1_deprecation = (legacy_fields, replacement)

        object.__setattr__(self, "b1_inh_scale", None)
        object.__setattr__(self, "b1_inh_res", None)
        return self

    @property
    def legacy_b1_deprecation(self) -> tuple[tuple[str, ...], str] | None:
        """Return legacy fields and suggested replacement when explicitly used."""
        return self._legacy_b1_deprecation

    def get_b1_nominal(self) -> float:
        """Nominal B1 value for distribution generation."""
        if self.b1_frq is not None:
            return float(self.b1_frq)
        if self.pw90 is not None:
            return 1.0 / (4.0 * float(self.pw90))
        msg = "Either 'b1_frq' or 'pw90' must be provided"
        raise ValueError(msg)

    def get_b1_oscillation_inputs(self) -> tuple[float, float]:
        """Return (effective_B1_Hz, saturation_duration_s) for auto-res."""
        t_sat = float(getattr(self, "time_t1", 0.5))
        return self.get_b1_nominal(), t_sat

    def get_b1_model(self) -> B1RuntimeModel:
        """Lower b1_distribution to an immutable runtime B1 model."""
        if self.b1_distribution is None:
            return HomogeneousB1Model()
        b1_eff_hz, t_sat_s = self.get_b1_oscillation_inputs()
        return self.b1_distribution.to_b1_model(
            b1_eff_hz=b1_eff_hz,
            t_sat_s=t_sat_s,
        )

    def get_b1_distribution(self) -> Distribution:
        """Generate concrete B1 nodes and weights for this experiment."""
        nominal = self.get_b1_nominal()
        model = self.get_b1_model()
        if model is None:
            return HomogeneousB1Model().get_distribution(nominal)
        return model.get_distribution(nominal)


class CestB1Settings(B1InhomogeneityMixin):
    """Mandatory explicit B1 configuration for CEST-family experiments."""

    allows_dephasing: ClassVar[bool] = True
    requires_cw_b1_xor: ClassVar[bool] = True

    @model_validator(mode="after")
    def _validate_cest_b1_settings(self) -> Self:
        if self.requires_cw_b1_xor:
            if self.b1_frq is None and self.pw90 is None:
                msg = (
                    "Either 'b1_frq' (Hz) or 'pw90' (s) must be provided in "
                    "[experiment] to define the B1 field."
                )
                raise ValueError(msg)
            if self.b1_frq is not None and self.pw90 is not None:
                msg = "Use either 'b1_frq' or 'pw90' in [experiment], not both."
                raise ValueError(msg)

        if self.b1_distribution is None:
            options = [
                'b1_distribution = { type = "gaussian", scale = 0.1 }',
                (
                    'b1_distribution = { type = "skewed", scale = 0.1, '
                    'skewness = -0.5, center = "mode" }'
                ),
            ]
            if self.allows_dephasing:
                options.append('b1_distribution = { type = "dephasing" }')
            options.extend(
                (
                    (
                        'b1_distribution = { type = "custom", scales = [...], '
                        "weights = [...] }"
                    ),
                    'b1_distribution = { type = "homogeneous" }',
                )
            )
            msg = (
                "CEST-family experiments require an explicit 'b1_distribution' "
                "in [experiment]. Choose one of: " + "; ".join(options)
            )
            raise ValueError(msg)

        if self.b1_distribution.type == "dephasing" and not self.allows_dephasing:
            exp_name = getattr(self, "name", "this experiment")
            msg = (
                f'b1_distribution = {{ type = "dephasing" }} is not physically valid '
                f"for {exp_name!r} because pulse-train or shaped-pulse saturation "
                "cannot be approximated by continuous-wave transverse dephasing. "
                'Use type = "gaussian", type = "skewed", type = "custom", or '
                'type = "homogeneous".'
            )
            raise ValueError(msg)

        return self


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
        return [f"{component}_{state}" for state in states for component in components]


class CpmgSettings(RelaxationSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(RelaxationSettings):
    even_ncycs: bool = True


class CestSettings(RelaxationSettings):
    sw: float = math.inf


class MFCestSettings(RelaxationSettings):
    sw: PositiveFloat = Field(description="Spectral width (Hz)")


class DCestSettings(MFCestSettings, CestB1Settings):
    """Shared B1 and DANTE pulse-train settings for D-CEST experiments."""

    allows_dephasing: ClassVar[bool] = False
    requires_cw_b1_xor: ClassVar[bool] = False

    time_t1: Delay = Field(description="CEST relaxation delay (seconds)")
    time_equil: Delay = Field(
        default=0.0,
        description="Equilibration delay at end of CEST period (seconds)",
    )
    b1_eff: B1Field = Field(
        ...,
        description="Effective B1 field for DANTE excitation (Hz)",
    )

    _normalize_b1_eff_alias = model_validator(mode="before")(normalize_b1_eff_alias)

    @model_validator(mode="after")
    def _validate_dcest_fields(self) -> Self:
        """Validate D-CEST hardware pw90 and physical DANTE duty cycle."""
        if self.pw90 is None:
            msg = (
                "D-CEST experiments require 'pw90' (hardware B1 pulse width). "
                "Both 'pw90' and 'b1_eff' (alias: 'b1_frq') must be specified."
            )
            raise ValueError(msg)
        if self.tau_dante < 0.0:
            duty = 4.0 * float(self.pw90) * float(self.b1_eff)
            msg = (
                f"D-CEST DANTE pulse width ({self.pw_dante:.3e} s) exceeds "
                f"1 / sw ({1.0 / self.sw:.3e} s); duty cycle 4 * pw90 * b1_eff = "
                f"{duty:.3f} must be <= 1.0."
            )
            raise ValueError(msg)
        return self

    # D-CEST convention: the B1 inhomogeneity distribution is centered on the
    # hardware B1 set by pw90, while DANTE nutation cycles track b1_eff.
    def get_b1_nominal(self) -> float:
        """Get nominal B1 from hardware pw90 for distribution centering."""
        if self.pw90 is None:
            msg = (
                "For D-CEST, 'pw90' must be specified to define the nominal B1 "
                "used for the B1 distribution."
            )
            raise ValueError(msg)
        return 1.0 / (4.0 * float(self.pw90))

    def get_b1_oscillation_inputs(self) -> tuple[float, float]:
        """Return effective DANTE CW field and train duration for auto-res."""
        return float(self.b1_eff), self.ncyc_dante / float(self.sw)

    @computed_field
    @property
    def pw_dante(self) -> float:
        """Duration of each DANTE pulse (seconds)."""
        if self.pw90 is None:
            msg = "pw90 must be specified for DANTE pulse calculation"
            raise ValueError(msg)
        return 4.0 * float(self.pw90) * float(self.b1_eff) / float(self.sw)

    @computed_field
    @property
    def tau_dante(self) -> float:
        """Inter-pulse delay in DANTE train (seconds)."""
        return 1.0 / self.sw - self.pw_dante

    @computed_field
    @property
    def ncyc_dante(self) -> int:
        """Number of DANTE pulses in the train."""
        return int(self.time_t1 * self.sw + 0.1)
