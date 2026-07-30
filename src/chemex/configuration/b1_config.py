"""B1 field inhomogeneity configuration.

This module wires ChemEx configuration to the B1 distribution plugin system,
building the discriminated union of distribution configs dynamically from the
registered plugins and parsing TOML input against it. Adding a new
distribution requires no changes here.
"""

from __future__ import annotations

import operator
from functools import reduce
from typing import Annotated, cast

from pydantic import Field, TypeAdapter

from chemex.nmr.distributions import registry  # Triggers plugin auto-registration
from chemex.nmr.distributions.registry import DistributionConfig


def _build_distribution_union() -> object:
    """Construct a Union of all registered distribution config classes."""
    classes = registry.get_all_config_classes()
    if not classes:
        msg = "No B1 distribution configs registered. Check plugin loading."
        raise RuntimeError(msg)
    return reduce(operator.or_, classes)


def _build_distribution_adapter() -> TypeAdapter:
    """Construct a discriminated TypeAdapter for the registered distribution configs."""
    distribution_union = _build_distribution_union()
    discriminated_union = Annotated.__getitem__(
        (distribution_union, Field(discriminator="type"))
    )
    return TypeAdapter(discriminated_union)


# Discriminated union of all dynamically-registered distribution configs
_DISTRIBUTION_ADAPTER = _build_distribution_adapter()


def parse_distribution_config(value: object) -> DistributionConfig:
    """Parse a distribution config using the registered plugin union."""
    if isinstance(value, DistributionConfig):
        return value
    if isinstance(value, dict):
        raw_config = cast("dict[str, object]", value)
        distribution_type = raw_config.get("type")
        if isinstance(distribution_type, str):
            registry.get_generator(distribution_type)
    parsed = _DISTRIBUTION_ADAPTER.validate_python(value)
    if isinstance(parsed, DistributionConfig):
        return parsed
    msg = "Invalid B1 distribution configuration"
    raise TypeError(msg)
