"""B1 field inhomogeneity configuration.

This module wires ChemEx configuration to the B1 distribution plugin system.
It supports both simple (float) and advanced (flat TOML table) configs,
building the discriminated union of distribution configs dynamically from the
registered plugins. Adding a new distribution requires no changes here.
"""

from __future__ import annotations

import operator
from functools import reduce
from typing import TYPE_CHECKING, Annotated, Self

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    field_validator,
    model_validator,
)

if TYPE_CHECKING:
    from chemex.nmr.constants import Distribution

from chemex.nmr.distributions import registry  # Triggers plugin auto-registration
from chemex.nmr.distributions.registry import DistributionConfig


def _build_distribution_union() -> object:
    """Construct a Union of all registered distribution config classes."""
    classes = registry.get_all_config_classes()
    if not classes:
        msg = "No B1 distribution configs registered. Check plugin loading."
        raise RuntimeError(msg)
    return reduce(operator.or_, classes)


def _default_distribution() -> DistributionConfig:
    config = registry.get_config_class("gaussian").model_validate({"type": "gaussian"})
    if isinstance(config, DistributionConfig):
        return config
    msg = "Invalid default B1 distribution configuration"
    raise TypeError(msg)


def _build_distribution_adapter() -> TypeAdapter:
    """Construct a discriminated TypeAdapter for the registered distribution configs."""
    distribution_union = _build_distribution_union()
    discriminated_union = Annotated.__getitem__(
        (distribution_union, Field(discriminator="type"))
    )
    return TypeAdapter(discriminated_union)


# Discriminated union of all dynamically-registered distribution configs
_DISTRIBUTION_ADAPTER = _build_distribution_adapter()


class B1FieldConfig(BaseModel):
    """Complete B1 field configuration including nominal value and distribution.

    This model is used when b1_frq is specified as a table in TOML.

    Attributes
    ----------
    value : float | None
        Nominal B1 field value (in kHz). Either this or pw90 must be specified.
    pw90 : float | None
        90-degree pulse width (in seconds). If specified, value is calculated
        as 1/(4*pw90).
    distribution : DistributionConfig
        Distribution configuration (dynamically typed based on registered plugins)

    """

    model_config = ConfigDict(
        extra="forbid",  # Reject unknown fields at top level
        validate_assignment=True,  # Validate on attribute assignment
        str_strip_whitespace=True,  # Auto-strip strings
    )

    value: float | None = Field(
        default=None,
        gt=0.0,
        description="Nominal B1 field strength (Hz or kHz)",
        examples=[25.0, 700.0],
    )
    pw90: float | None = Field(
        default=None,
        gt=0.0,
        description="90-degree pulse width (seconds)",
        examples=[45e-6, 96e-6],
    )
    # Default to Gaussian distribution if available
    distribution: DistributionConfig = Field(
        default_factory=_default_distribution,
        description="B1 distribution configuration",
    )

    # Accept only flat TOML schema: all distribution params at top level
    @model_validator(mode="before")
    @classmethod
    def _normalize_flat_schema(cls, data: object) -> object:
        if not isinstance(data, dict):
            return data

        flat_schema: dict[str, object] = {}
        for key, value in data.items():
            if not isinstance(key, str):
                return data
            flat_schema[key] = value

        # Reject nested schema (where distribution is a sub-dict)
        if "distribution" in flat_schema and isinstance(
            flat_schema["distribution"], dict
        ):
            msg = (
                "Nested distribution schema is not supported. "
                "Use flat schema: put all keys (value, type, scale, etc.) "
                "at the same level under [experiment.b1_frq]"
            )
            raise ValueError(msg)

        # If top-level contains a distribution type alongside value or pw90,
        # fold keys into a nested 'distribution' object internally
        if ("value" in flat_schema or "pw90" in flat_schema) and "type" in flat_schema:
            # Move all keys except 'value' and 'pw90' into the 'distribution' sub-dict
            distribution = {
                key: value
                for key, value in flat_schema.items()
                if key not in {"value", "pw90"}
            }
            result: dict[str, object] = {"distribution": distribution}
            if "value" in flat_schema:
                result["value"] = flat_schema["value"]
            if "pw90" in flat_schema:
                result["pw90"] = flat_schema["pw90"]
            return result

        return flat_schema

    @model_validator(mode="after")
    def _validate_value_or_pw90(self) -> Self:
        """Ensure exactly one of value or pw90 is specified."""
        has_value = self.value is not None
        has_pw90 = self.pw90 is not None

        if not (has_value ^ has_pw90):  # XOR: exactly one must be True
            msg = (
                "Exactly one of 'value' or 'pw90' must be specified in B1 field "
                "config, not both or neither"
            )
            raise ValueError(msg)
        return self

    # Parse the distribution using the dynamic union adapter
    @field_validator("distribution", mode="before")
    @classmethod
    def _parse_distribution(cls, value: object) -> DistributionConfig:
        """Parse distribution config using dynamic union adapter."""
        if isinstance(value, DistributionConfig):
            return value
        parsed = _DISTRIBUTION_ADAPTER.validate_python(value)
        if isinstance(parsed, DistributionConfig):
            return parsed
        msg = "Invalid B1 distribution configuration"
        raise TypeError(msg)

    @property
    def b1_nominal(self) -> float:
        """Get the nominal B1 field value.

        Returns
        -------
        float
            B1 field value in kHz (calculated from pw90 if necessary)

        """
        if self.value is not None:
            return self.value
        # At this point, pw90 is guaranteed to be not None
        return 1.0 / (4.0 * self.pw90)  # type: ignore[operator]

    def get_distribution(self) -> Distribution:
        """Generate the B1 distribution based on the configuration.

        Returns
        -------
        Distribution
            B1 values and weights for numerical integration

        """
        return self.distribution.get_distribution(self.b1_nominal)


# Type alias for the b1_frq field that accepts either float or table
# When float: use default Gaussian with separate b1_inh_scale/b1_inh_res
# When dict/table: use B1FieldConfig
B1Field = float | B1FieldConfig
