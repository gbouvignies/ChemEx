"""B1 field inhomogeneity configuration schemas.

Defines the static discriminated union of supported B1 distribution models and
lowers validated configurations into pure runtime B1 models with fixed
quadrature resolution.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Annotated, Any, Literal, cast

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    field_validator,
    model_validator,
)

from chemex.nmr.b1 import (
    B1RuntimeModel,
    CustomB1Model,
    DephasingB1Model,
    GaussianB1Model,
    HomogeneousB1Model,
    SkewedB1Model,
)
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions._quadrature import (
    CUSTOM_MAX_POINTS,
    MAX_SKEWNESS,
    RES_MAX,
    resolve_gaussian_res,
    resolve_skewed_res,
    validate_skewness,
)

VALID_B1_DISTRIBUTION_TYPES: tuple[str, ...] = (
    "homogeneous",
    "dephasing",
    "gaussian",
    "skewed",
    "custom",
)


def _normalize_mapping(data: Mapping[Any, object]) -> dict[str, object]:
    normalized = {str(k).lower(): v for k, v in data.items()}
    for key in ("type", "center"):
        if isinstance(normalized.get(key), str):
            normalized[key] = str(normalized[key]).lower()
    return normalized


def _normalize_str_keys(data: object) -> object:
    if not isinstance(data, Mapping):
        return data
    return _normalize_mapping(data)


class DistributionConfig(BaseModel):
    """Base configuration model for B1 field inhomogeneity specifications."""

    model_config = ConfigDict(extra="forbid", frozen=True, str_to_lower=True)

    type: str

    @model_validator(mode="before")
    @classmethod
    def _normalize_inputs(cls, data: object) -> object:
        return _normalize_str_keys(data)

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float,
        t_sat_s: float = 0.5,
    ) -> B1RuntimeModel:
        """Lower configuration to an immutable runtime B1 model with locked res."""
        msg = f"{type(self).__name__} must implement 'to_b1_model()'"
        raise NotImplementedError(msg)

    def get_distribution(
        self,
        value: float,
        *,
        b1_eff_hz: float | None = None,
        t_sat_s: float = 0.5,
    ) -> Distribution:
        """Build the concrete B1 distribution for a nominal field value."""
        effective_hz = float(value) if b1_eff_hz is None else float(b1_eff_hz)
        runtime_model = self.to_b1_model(b1_eff_hz=effective_hz, t_sat_s=t_sat_s)
        if runtime_model is None:
            return HomogeneousB1Model().get_distribution(value)
        return runtime_model.get_distribution(value)


class HomogeneousDistributionConfig(DistributionConfig):
    """Configuration for an ideal single-point homogeneous B1 field."""

    type: Literal["homogeneous"] = "homogeneous"

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float = 0.0,
        t_sat_s: float = 0.5,
    ) -> HomogeneousB1Model:
        return HomogeneousB1Model()


class DephasingConfig(DistributionConfig):
    """Configuration for complete transverse dephasing in pure-CW CEST."""

    type: Literal["dephasing"] = "dephasing"

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float = 0.0,
        t_sat_s: float = 0.5,
    ) -> DephasingB1Model:
        return DephasingB1Model()


class GaussianDistributionConfig(DistributionConfig):
    """Configuration for positive-conditioned Gaussian B1 distribution."""

    type: Literal["gaussian"] = "gaussian"
    scale: float = Field(default=0.1, ge=0.0, le=1.0, allow_inf_nan=False)
    res: int | None = Field(default=None, ge=1, le=RES_MAX)

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float,
        t_sat_s: float = 0.5,
    ) -> GaussianB1Model:
        resolved_res = resolve_gaussian_res(
            self.scale,
            b1_eff_hz=b1_eff_hz,
            t_sat_s=t_sat_s,
            res=self.res,
        )
        return GaussianB1Model(scale=self.scale, res=resolved_res)


class SkewedDistributionConfig(DistributionConfig):
    """Configuration for positive-conditioned Azzalini skew-normal B1."""

    type: Literal["skewed"] = "skewed"
    scale: float = Field(default=0.1, ge=0.0, le=1.0, allow_inf_nan=False)
    skewness: float = Field(
        ...,
        allow_inf_nan=False,
        description=(
            "Theoretical third standardized moment of the parent Azzalini "
            f"skew-normal distribution (|skewness| < {MAX_SKEWNESS:.12g})"
        ),
    )
    center: Literal["mode", "mean"] = Field(
        ...,
        description=(
            "Which summary of the parent skew-normal law equals nominal B1 "
            "('mode' or 'mean')"
        ),
    )
    res: int | None = Field(default=None, ge=1, le=RES_MAX)

    @model_validator(mode="before")
    @classmethod
    def _require_explicit_center(cls, data: object) -> object:
        normalized = _normalize_str_keys(data)
        if isinstance(normalized, Mapping):
            center = normalized.get("center")
            if center is None:
                msg = (
                    "'skewed' B1 distribution requires an explicit 'center' "
                    '("mode" or "mean"), e.g.: '
                    'b1_distribution = { type = "skewed", scale = 0.1, '
                    'skewness = -0.5, center = "mode" }'
                )
                raise ValueError(msg)
            if center not in ("mode", "mean"):
                msg = (
                    f"Invalid 'center' {center!r} for 'skewed' B1 distribution; "
                    'must be "mode" or "mean".'
                )
                raise ValueError(msg)
        return normalized

    @field_validator("skewness")
    @classmethod
    def _check_skewness_bounds(cls, value: float) -> float:
        return validate_skewness(value)

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float,
        t_sat_s: float = 0.5,
    ) -> SkewedB1Model:
        resolved_res = resolve_skewed_res(
            self.scale,
            self.skewness,
            self.center,
            b1_eff_hz=b1_eff_hz,
            t_sat_s=t_sat_s,
            res=self.res,
        )
        return SkewedB1Model(
            scale=self.scale,
            skewness=self.skewness,
            center=self.center,
            res=resolved_res,
        )


class CustomDistributionConfig(DistributionConfig):
    """Configuration for custom B1 distribution with explicit scales and weights."""

    type: Literal["custom"] = "custom"
    scales: list[float] = Field(min_length=1, max_length=CUSTOM_MAX_POINTS)
    weights: list[float] = Field(min_length=1, max_length=CUSTOM_MAX_POINTS)

    @field_validator("scales")
    @classmethod
    def _validate_scales_positive(cls, values: list[float]) -> list[float]:
        if any(not _is_finite_positive(scale) for scale in values):
            msg = "All scales must be finite and positive"
            raise ValueError(msg)
        return values

    @field_validator("weights")
    @classmethod
    def _validate_weights_positive(cls, values: list[float]) -> list[float]:
        if any(not _is_finite_positive(weight) for weight in values):
            msg = "All weights must be finite and positive"
            raise ValueError(msg)
        return values

    @model_validator(mode="after")
    def _validate_same_length(self) -> CustomDistributionConfig:
        if len(self.scales) != len(self.weights):
            msg = "scales and weights must have the same length"
            raise ValueError(msg)
        return self

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float = 0.0,
        t_sat_s: float = 0.5,
    ) -> CustomB1Model:
        return CustomB1Model(
            scales=tuple(self.scales),
            weights=tuple(self.weights),
        )


def _is_finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


type B1DistributionModelConfig = (
    HomogeneousDistributionConfig
    | DephasingConfig
    | GaussianDistributionConfig
    | SkewedDistributionConfig
    | CustomDistributionConfig
)

B1DistributionConfig = Annotated[
    B1DistributionModelConfig,
    Field(discriminator="type"),
]

_DISTRIBUTION_ADAPTER: TypeAdapter[B1DistributionModelConfig] = TypeAdapter(
    B1DistributionConfig
)


def _validate_distribution_type(dist_type: object) -> None:
    lowered = dist_type.lower() if isinstance(dist_type, str) else None
    if lowered is None:
        valid = ", ".join(repr(name) for name in VALID_B1_DISTRIBUTION_TYPES)
        msg = f"B1 distribution table must specify 'type' (one of: {valid})"
        raise ValueError(msg)

    if lowered == "hermite":
        msg = (
            "B1 distribution type 'hermite' has been consolidated into 'gaussian', "
            "which uses smooth positive-support Gauss-Legendre quadrature with "
            "automatic resolution. Replace with: "
            'b1_distribution = { type = "gaussian", scale = ... }'
        )
        raise ValueError(msg)
    if lowered in ("beta", "truncated_skewed"):
        msg = (
            f"B1 distribution type {lowered!r} has been removed because hard upper "
            "cutoffs at nominal B1 are unphysical. Use "
            'type = "skewed" (with explicit center = "mode" or "mean") for '
            'asymmetric profiles, type = "gaussian" for symmetric profiles, or '
            'type = "custom" for empirical nutation profiles.'
        )
        raise ValueError(msg)
    if lowered not in VALID_B1_DISTRIBUTION_TYPES:
        valid = ", ".join(VALID_B1_DISTRIBUTION_TYPES)
        msg = f"Unknown B1 distribution type: {dist_type}. Valid types: {valid}"
        raise ValueError(msg)


def _malformed_distribution_message(value: object) -> str:
    if isinstance(value, str):
        return (
            f"'b1_distribution' must be a table (got string {value!r}), e.g.: "
            'b1_distribution = { type = "gaussian", scale = 0.1 }'
        )
    valid = ", ".join(repr(name) for name in VALID_B1_DISTRIBUTION_TYPES)
    return (
        f"'b1_distribution' must be a table specifying 'type' "
        f"(one of: {valid}), got {type(value).__name__}"
    )


def parse_distribution_config(value: object) -> B1DistributionModelConfig:
    """Parse a B1 distribution config against the curated discriminated union."""
    if isinstance(value, DistributionConfig):
        return cast("B1DistributionModelConfig", value)
    if isinstance(value, Mapping):
        normalized = _normalize_mapping(value)
        _validate_distribution_type(normalized.get("type"))
        return _DISTRIBUTION_ADAPTER.validate_python(normalized)
    msg = _malformed_distribution_message(value)
    raise ValueError(msg)
