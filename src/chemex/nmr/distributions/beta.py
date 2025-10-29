"""Beta B1 distribution for modeling upper-bounded B1 inhomogeneity.

The beta distribution is bounded on [0, value], making it ideal for modeling
B1 field inhomogeneity where the nominal value represents an upper limit
(you cannot exceed the applied RF power) but the field can degrade below that.

This is physically appropriate for:
- Coil inhomogeneity (field is weaker away from center, never stronger)
- Sample loading effects (RF penetration limited by nominal power)
- Nutation experiments showing asymmetric profiles with sharp upper cutoff

For unbounded or symmetric distributions, use gaussian, hermite, or skewed instead.
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
    res: int = 11,
    mean_frac: float = 0.95,
) -> Distribution:
    """Generate upper-bounded B1 distribution using Beta distribution.

    Creates a distribution on [0, value] where 'value' is the upper bound
    (nominal B1 field). The distribution is parameterized by its mean position
    and scale (coefficient of variation).

    Parameters
    ----------
    value : float
        Nominal B1 field value (upper bound of distribution, in kHz or Hz)
    scale : float
        Coefficient of variation (CV = std/mean) controlling width.
        Typical range: 0.01 to 0.5
        - scale = 0.05: tight around mean (CV=5%)
        - scale = 0.10: moderate spread (CV=10%)
        - scale = 0.20: wide spread (CV=20%)
    res : int
        Number of integration points (quantile nodes with equal weights)
    mean_frac : float
        Mean position as fraction of nominal value.
        Default 0.95 means distribution peaks near 95% of nominal.
        Range: 0.5 to 0.99 (must be < 1.0 for proper beta shape)

    Returns
    -------
    Distribution
        B1 values on [0, value] with equal weights (quantile discretization)

    Notes
    -----
    The beta distribution parameters α and β are derived from the mean and CV:
    - mean_rel = mean_frac (on [0,1] scale)
    - variance = (mean_rel * scale)² 
    - α, β computed via method of moments

    Physical interpretation:
    - Upper bound = value: You cannot exceed the nominal RF power
    - Mean < value: Average field is reduced due to inhomogeneity
    - Asymmetric: Field degrades below nominal more than it varies above

    Examples
    --------
    >>> # Tight distribution peaked at 95% of 20 kHz nominal
    >>> dist = generate(value=20.0, scale=0.05, res=11, mean_frac=0.95)
    >>> # Wider distribution peaked at 90% of 20 kHz
    >>> dist = generate(value=20.0, scale=0.15, res=11, mean_frac=0.90)

    """
    # Handle edge cases
    if scale in (0.0, np.inf) or res <= 1:
        return Distribution(np.array([value], dtype=float), np.array([1.0], dtype=float))

    # Clamp inputs to valid ranges
    scale = np.clip(scale, 0.01, 0.5)
    mean_frac = np.clip(mean_frac, 0.5, 0.99)

    # Target mean and std on [0, 1] scale
    mu = mean_frac
    sigma = mu * scale  # CV definition: sigma/mu = scale

    # Convert to Beta parameters via method of moments
    # For Beta(α, β) on [0,1]: mean = α/(α+β), var = αβ/[(α+β)²(α+β+1)]
    t = mu * (1 - mu) / (sigma**2) - 1

    if t <= 0:
        # Variance too large for this mean; use fallback
        # Choose shape that gives desired mean with moderate spread
        alpha = 2 * mu
        beta_param = 2 * (1 - mu)
    else:
        alpha = mu * t
        beta_param = (1 - mu) * t

    # Generate quantile nodes on [0, 1]
    # Mid-quantiles: p_i = (i - 0.5) / res for i = 1, ..., res
    p = (np.arange(1, res + 1) - 0.5) / res
    rel_vals = stats.beta.ppf(p, alpha, beta_param)

    # Scale to [0, value]
    b1_values = value * rel_vals

    # Equal weights (quantile discretization)
    weights = np.ones(res, dtype=float) / res

    return Distribution(b1_values, weights)


class BetaDistributionConfig(BaseModel):
    """Configuration for upper-bounded B1 distribution using Beta distribution.

    The beta distribution is bounded on [0, nominal], making it ideal for
    modeling B1 field inhomogeneity where the nominal value is an upper limit.

    Attributes
    ----------
    type : Literal["beta"]
        Distribution type identifier
    scale : float
        Coefficient of variation (CV = std/mean), controlling width.
        Recommended range: 0.01 to 0.5
    res : int
        Number of integration points (1 to 51)
    mean_frac : float
        Mean position as fraction of nominal value (0.5 to 0.99).
        Default 0.95 places peak near 95% of nominal.

    Examples
    --------
    TOML configuration:

    # Moderate spread, peaked at 95% of nominal
    b1_distribution = { type = "beta", scale = 0.10, res = 11, mean_frac = 0.95 }

    # Wider spread, peaked at 90% of nominal
    b1_distribution = { type = "beta", scale = 0.15, res = 11, mean_frac = 0.90 }

    """

    type: Literal["beta"] = "beta"
    scale: float = Field(default=0.1, ge=0.01, le=0.5)
    res: int = Field(default=11, ge=1, le=51)
    mean_frac: float = Field(default=0.95, ge=0.5, lt=1.0)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the Beta distribution."""
        return generate(value, self.scale, self.res, self.mean_frac)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("beta", generate, BetaDistributionConfig)
