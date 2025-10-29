"""Gaussian B1 distribution with linear grid.

This is the original ChemEx B1 inhomogeneity distribution method.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, Field
from scipy import stats

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import registry


def generate(
    value: float,
    scale: float,
    res: int,
) -> Distribution:
    """Generate simple Gaussian distribution with linear grid.

    Uses a linear grid from -2σ to +2σ with Gaussian PDF weights.
    This is the original ChemEx implementation for backward compatibility.

    Parameters
    ----------
    value : float
        Nominal B1 field value (kHz or Hz)
    scale : float
        Standard deviation as fraction of value (e.g., 0.10 = 10%)
    res : int
        Number of grid points

    Returns
    -------
    Distribution
        B1 values and Gaussian PDF weights

    Notes
    -----
    - Grid spans [-2, 2] standard deviations
    - Weights are normalized Gaussian PDF values
    - Original ChemEx implementation for backward compatibility

    """
    if scale not in (0.0, np.inf) and res > 1:
        grid = np.linspace(-2.0, 2.0, res)
        distribution = grid * scale + 1.0
    else:
        grid = np.array([0.0])
        distribution = np.array([1.0])

    weights = stats.norm.pdf(grid)
    weights /= weights.sum()

    return Distribution(distribution * value, weights)


class GaussianDistributionConfig(BaseModel):
    """Configuration for simple Gaussian B1 distribution with linear grid.

    This is the original ChemEx B1 inhomogeneity method using a linear grid
    from -2 to +2 standard deviations with Gaussian PDF weights.

    Attributes
    ----------
    type : Literal["gaussian"]
        Distribution type identifier
    scale : float
        Standard deviation as fraction of nominal value (e.g., 0.10 = 10%)
    res : int
        Number of grid points

    """

    type: Literal["gaussian"] = "gaussian"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    res: int = Field(default=11, ge=1, le=51)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the Gaussian distribution."""
        return generate(value, self.scale, self.res)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("gaussian", generate, GaussianDistributionConfig)
