"""Homogeneous single-point B1 field specification."""

from __future__ import annotations

from typing import Any

import numpy as np

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions._quadrature import validate_nominal_and_scale


def generate(value: float) -> Distribution:
    """Generate a single-point ideal homogeneous B1 distribution."""
    validate_nominal_and_scale(value, 0.0)
    return Distribution(
        values=np.array([float(value)], dtype=float),
        weights=np.array([1.0], dtype=float),
    )


def __getattr__(name: str) -> Any:
    if name == "HomogeneousDistributionConfig":
        from chemex.configuration.b1_config import HomogeneousDistributionConfig

        return HomogeneousDistributionConfig
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
