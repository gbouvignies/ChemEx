"""Hermite B1 distribution using Gauss-Hermite quadrature.

Advanced method using Hermite polynomial roots for optimal integration
of Gaussian-weighted functions.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, Field
from scipy.special import roots_hermite

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import registry


def generate(
    value: float,
    scale: float,
    res: int,
) -> Distribution:
    """Generate Gaussian distribution using Gauss-Hermite quadrature.

    Uses Hermite polynomial roots for optimal integration of Gaussian-weighted
    functions.

    Parameters
    ----------
    value : float
        Nominal B1 field value (kHz or Hz)
    scale : float
        Standard deviation as fraction of value (e.g., 0.10 = 10%)
    res : int
        Number of integration points (Hermite nodes)

    Returns
    -------
    Distribution
        B1 values and weights for Gauss-Hermite integration

    Notes
    -----
    - Uses Gauss-Hermite quadrature nodes and weights
    - Optimal for integrating functions with Gaussian weighting
    - Weights are pre-normalized for the Gaussian distribution

    """
    if scale in (0.0, np.inf) or res <= 1:
        return Distribution(np.array([value]), np.array([1.0]))

    nodes, weights = roots_hermite(res)

    distribution = nodes * np.sqrt(2) * scale + 1.0

    # Normalize weights to sum to 1.0
    weights = weights / weights.sum()

    return Distribution(distribution * value, weights)


class HermiteDistributionConfig(BaseModel):
    """Configuration for B1 distribution using Gauss-Hermite quadrature.

    Advanced method using Hermite polynomial roots for optimal integration
    of Gaussian-weighted functions.

    Attributes
    ----------
    type : Literal["hermite"]
        Distribution type identifier
    scale : float
        Standard deviation as fraction of nominal value (e.g., 0.10 = 10%)
    res : int
        Number of integration points (Hermite nodes)

    """

    type: Literal["hermite"] = "hermite"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    res: int = Field(default=11, ge=1, le=51)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the Hermite distribution."""
        return generate(value, self.scale, self.res)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("hermite", generate, HermiteDistributionConfig)
