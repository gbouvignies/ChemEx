"""Custom B1 distribution with user-defined relative scales and weights."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from chemex.nmr.constants import Distribution


def generate(
    value: float,
    scales: Sequence[float],
    weights: Sequence[float],
) -> Distribution:
    """Generate custom B1 distribution from user-specified scales and weights."""
    if not math.isfinite(value) or value < 0.0:
        msg = f"Nominal B1 value must be finite and non-negative, got {value}"
        raise ValueError(msg)

    scales_array = np.asarray(scales, dtype=float)
    weights_array = np.asarray(weights, dtype=float)

    if scales_array.ndim != 1 or scales_array.size == 0:
        msg = "Custom B1 scales must be a non-empty 1D sequence"
        raise ValueError(msg)
    if weights_array.ndim != 1 or weights_array.size != scales_array.size:
        msg = "scales and weights must have the same non-empty length"
        raise ValueError(msg)
    if not np.all(np.isfinite(scales_array)) or np.any(scales_array <= 0.0):
        msg = "All scales must be finite and positive"
        raise ValueError(msg)
    if not np.all(np.isfinite(weights_array)) or np.any(weights_array <= 0.0):
        msg = "All weights must be finite and positive"
        raise ValueError(msg)

    normalized_weights = weights_array / np.sum(weights_array)
    return Distribution(float(value) * scales_array, normalized_weights)


def __getattr__(name: str) -> Any:
    if name == "CustomDistributionConfig":
        from chemex.configuration.b1_config import CustomDistributionConfig

        return CustomDistributionConfig
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
