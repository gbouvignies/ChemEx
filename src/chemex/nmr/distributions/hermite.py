"""Hermite B1 distribution using Gauss-Hermite quadrature.

Advanced method using Hermite polynomial roots for optimal integration
of Gaussian-weighted functions, switching to positive-truncated Gaussian
Gauss-Legendre quadrature when wide scales would place Hermite nodes at B1 < 0.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.polynomial.legendre import leggauss
from pydantic import Field
from scipy import stats
from scipy.special import roots_hermite

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import DistributionConfig, registry


def generate(
    value: float,
    scale: float,
    res: int,
) -> Distribution:
    """Generate Gaussian B1 distribution using Gauss-Hermite quadrature.

    Uses exact Gauss-Hermite nodes and weights whenever all nodes lie at
    non-negative B1 values. When ``scale`` and ``res`` are large enough that
    unbounded Gauss-Hermite nodes would become negative, evaluates the
    Gaussian density restricted to ``B1 >= 0`` via ``res``-point Gauss-Legendre
    quadrature so all quadrature points remain non-negative.
    """
    if value < 0.0:
        msg = f"Nominal B1 value must be non-negative, got {value}"
        raise ValueError(msg)
    if scale < 0.0:
        msg = f"Scale must be non-negative, got {scale}"
        raise ValueError(msg)

    if value == 0.0 or scale in (0.0, np.inf) or res <= 1:
        return Distribution(np.array([float(value)]), np.array([1.0]))

    nodes, weights = roots_hermite(res)
    distribution = nodes * np.sqrt(2.0) * scale + 1.0

    if distribution.min() >= 0.0:
        weights = weights / weights.sum()
        return Distribution(distribution * value, weights)

    # Positive-truncated Gaussian quadrature on [0, rel_max] so B1 >= 0 always.
    tail_sigma = max(5.0, float(-nodes[0] * np.sqrt(2.0)))
    rel_max = 1.0 + tail_sigma * scale
    gl_nodes, gl_weights = leggauss(res)
    rel_values = 0.5 * rel_max * (gl_nodes + 1.0)
    pdf_values = stats.norm.pdf(rel_values, loc=1.0, scale=scale)
    pos_weights = 0.5 * rel_max * gl_weights * pdf_values
    pos_weights = pos_weights / pos_weights.sum()
    return Distribution(rel_values * value, pos_weights)


class HermiteDistributionConfig(DistributionConfig):
    """Configuration for B1 distribution using Gauss-Hermite quadrature.

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
