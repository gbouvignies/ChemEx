"""Custom B1 distribution with user-defined scales and weights.

Allows users to specify explicit B1 scaling factors and their weights
directly in the TOML configuration.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, Field, field_validator
from pydantic_core.core_schema import ValidationInfo

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import registry


def generate(
    value: float,
    scales: list[float],
    weights: list[float],
) -> Distribution:
    """Generate custom B1 distribution from user-specified scales and weights.

    Parameters
    ----------
    value : float
        Nominal B1 field value (kHz or Hz)
    scales : list[float]
        Scaling factors relative to nominal value (e.g., [0.95, 1.0, 1.05])
    weights : list[float]
        Weights for each scale (will be normalized to sum to 1.0)

    Returns
    -------
    Distribution
        B1 values and normalized weights

    Examples
    --------
    >>> # Three-point distribution around nominal
    >>> dist = generate(
    ...     value=15.95,
    ...     scales=[0.95, 1.0, 1.05],
    ...     weights=[0.25, 0.5, 0.25]
    ... )
    >>> print(f"B1 values: {dist.values}")
    >>> print(f"Weights sum: {dist.weights.sum()}")

    Notes
    -----
    - Weights are automatically normalized to sum to 1.0
    - Scales are multiplicative factors applied to the nominal value
    - Use scales=[1.0] and weights=[1.0] for no inhomogeneity

    """
    scales_array = np.array(scales, dtype=float)
    weights_array = np.array(weights, dtype=float)

    # Normalize weights
    weights_array = weights_array / weights_array.sum()

    # Compute B1 values
    b1_values = value * scales_array

    return Distribution(b1_values, weights_array)


class CustomDistributionConfig(BaseModel):
    """Configuration for custom B1 distribution with user-defined points.

    Attributes
    ----------
    type : Literal["custom"]
        Distribution type identifier
    scales : list[float]
        Scaling factors relative to nominal value (e.g., [0.95, 1.0, 1.05])
        Each scale is multiplied by the nominal B1 value
    weights : list[float]
        Weights for each scale point (automatically normalized to sum to 1.0)
        Must have the same length as scales

    """

    type: Literal["custom"] = "custom"
    scales: list[float] = Field(min_length=1)
    weights: list[float] = Field(min_length=1)

    @field_validator("scales")
    @classmethod
    def _validate_scales_positive(cls, v: list[float]) -> list[float]:
        """Ensure all scales are positive."""
        if any(s <= 0.0 for s in v):
            msg = "All scales must be positive"
            raise ValueError(msg)
        return v

    @field_validator("weights")
    @classmethod
    def _validate_weights_positive(cls, v: list[float]) -> list[float]:
        """Ensure all weights are positive."""
        if any(w <= 0.0 for w in v):
            msg = "All weights must be positive"
            raise ValueError(msg)
        return v

    @field_validator("weights")
    @classmethod
    def _validate_same_length(cls, v: list[float], info: ValidationInfo) -> list[float]:
        """Ensure scales and weights have the same length."""
        if "scales" in info.data and len(v) != len(info.data["scales"]):
            msg = "scales and weights must have the same length"
            raise ValueError(msg)
        return v

    def get_distribution(self, value: float) -> Distribution:
        """Generate the custom distribution."""
        return generate(value, self.scales, self.weights)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("custom", generate, CustomDistributionConfig)
