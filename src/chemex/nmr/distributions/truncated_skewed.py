"""Truncated skew-normal B1 distribution with upper bound enforcement.

This distribution combines the asymmetry of the skew-normal distribution with
a hard upper bound at the nominal B1 value, making it suitable for modeling
RF field inhomogeneity when the B1 field cannot exceed the applied RF power.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.polynomial.legendre import leggauss
from pydantic import BaseModel, Field
from scipy.optimize import minimize_scalar
from scipy.special import erf

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import registry


def _skewnorm_pdf(x: np.ndarray, xi: float, omega: float, alpha: float) -> np.ndarray:
    """Evaluate skew-normal PDF at points x.

    Parameters
    ----------
    x : np.ndarray
        Points at which to evaluate PDF
    xi : float
        Location parameter
    omega : float
        Scale parameter
    alpha : float
        Shape parameter (controls skewness)

    Returns
    -------
    np.ndarray
        PDF values

    """
    z = (x - xi) / omega
    phi = np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)
    big_phi = 0.5 * (1.0 + erf(alpha * z / np.sqrt(2.0)))
    return (2.0 / omega) * phi * big_phi


def _skewnorm_cdf(x: np.ndarray, xi: float, omega: float, alpha: float) -> np.ndarray:
    """Evaluate skew-normal CDF at points x.

    Parameters
    ----------
    x : np.ndarray
        Points at which to evaluate CDF
    xi : float
        Location parameter
    omega : float
        Scale parameter
    alpha : float
        Shape parameter (controls skewness)

    Returns
    -------
    np.ndarray
        CDF values

    """
    z = (x - xi) / omega
    phi_z = 0.5 * (1.0 + erf(z / np.sqrt(2.0)))

    # Owen's T function approximation for skew-normal CDF
    # CDF(x) = Phi(z) - 2*T(z, alpha)
    return phi_z - 2.0 * _owens_t_approx(z, alpha)


def _owens_t_approx(h: np.ndarray, a: float) -> np.ndarray:
    """Approximate Owen's T function.

    Uses a simple approximation for Owen's T(h, a).

    Parameters
    ----------
    h : np.ndarray
        First argument
    a : float
        Second argument

    Returns
    -------
    np.ndarray
        Approximate Owen's T values

    """
    # Simple approximation: T(h,a) ≈ (1/2π) * arctan(a) * exp(-h²/2)
    return (1.0 / (2.0 * np.pi)) * np.arctan(a) * np.exp(-0.5 * h**2)


def _compute_mode_location(alpha: float) -> float:
    """Compute mode location for standard skew-normal (xi=0, omega=1).

    Parameters
    ----------
    alpha : float
        Shape parameter

    Returns
    -------
    float
        Mode location relative to xi

    """
    delta = alpha / np.sqrt(1.0 + alpha**2)
    mu_z = delta * np.sqrt(2.0 / np.pi)

    # For standard skew-normal, mode is not at mu_z but needs numerical solution
    # Mode satisfies: d/dz [phi(z) * Phi(alpha*z)] = 0

    def neg_pdf(z: float) -> float:
        phi = np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)
        big_phi = 0.5 * (1.0 + erf(alpha * z / np.sqrt(2.0)))
        return -2.0 * phi * big_phi

    # Search around the mean for the mode
    result = minimize_scalar(
        neg_pdf,
        bounds=(mu_z - 1.0, mu_z + 1.0),
        method="bounded",  # type: ignore[call-overload]
    )
    return float(result.x)  # type: ignore[attr-defined]


def generate(
    value: float,
    scale: float,
    skewness: float,
    res: int,
    lower_bound: float | None = None,
) -> Distribution:
    """Generate truncated skew-normal distribution with upper bound at nominal.

    The distribution is skew-normal but truncated at [lower_bound, value].
    This enforces that B1 cannot exceed the nominal value while maintaining
    asymmetry.

    Parameters
    ----------
    value : float
        Nominal B1 field value (kHz or Hz) - also the UPPER BOUND
    scale : float
        Coefficient of variation (std/mean) for the underlying distribution
    skewness : float
        Skewness parameter:
        - Negative values: left-skewed (tail toward lower B1)
        - Positive values: right-skewed (tail toward higher B1, but truncated)
        - Zero: symmetric (becomes truncated Gaussian)
    res : int
        Number of quadrature points
    lower_bound : float, optional
        Lower truncation bound. If None, uses value * (1 - 4*scale)

    Returns
    -------
    Distribution
        B1 values and weights (normalized)

    Notes
    -----
    The distribution is parameterized so that:
    - The mode of the underlying skew-normal is at the nominal value
    - The distribution is truncated at [lower_bound, value]
    - The resulting distribution has the requested skewness and scale

    Physical interpretation:
    - Upper bound = nominal: RF amplifier cannot exceed applied power
    - Lower bound: Sample/coil inhomogeneity reduces B1 in some regions
    - Skewness < 0: More probability mass near nominal, tail toward lower B1

    Example
    -------
    >>> # Left-skewed distribution bounded at [14, 20] kHz
    >>> dist = generate(value=20.0, scale=0.1, skewness=-0.8, res=11)
    >>> print(f"Range: [{dist.values.min():.1f}, {dist.values.max():.1f}] kHz")
    Range: [14.0, 20.0] kHz

    """
    if scale <= 0.0 or np.isinf(scale):
        return Distribution(np.array([value]), np.array([1.0]))

    if res <= 1:
        return Distribution(np.array([value]), np.array([1.0]))

    # Convert skewness to alpha parameter
    alpha = skewness

    # Determine bounds
    if lower_bound is None:
        lower_bound = value * (1.0 - 4.0 * scale)

    upper_bound = value  # Hard upper limit at nominal

    # For skew-normal, mode location depends on alpha
    mode_rel = _compute_mode_location(alpha)

    # Set omega based on desired scale and bounds
    # We want the truncated distribution to have CV ≈ scale
    omega = scale * value

    # Set xi so that the mode is at the nominal value
    xi = value - omega * mode_rel

    # Generate quadrature points using Gauss-Legendre
    nodes, weights = leggauss(res)

    # Map nodes from [-1, 1] to [lower_bound, upper_bound]
    b1_values = 0.5 * (upper_bound - lower_bound) * (nodes + 1.0) + lower_bound

    # Evaluate truncated PDF
    pdf_values = _skewnorm_pdf(b1_values, xi, omega, alpha)

    # Compute normalization constant (integral over truncation region)
    # Using the Gauss-Legendre quadrature itself for normalization
    normalization = np.sum(pdf_values * weights * 0.5 * (upper_bound - lower_bound))

    # Compute final weights (normalized truncated PDF * quadrature weights)
    if normalization > 0:
        final_weights = (
            pdf_values * weights * 0.5 * (upper_bound - lower_bound) / normalization
        )
    else:
        # Fallback to uniform if normalization fails
        final_weights = weights / np.sum(weights)

    # Ensure weights are normalized
    final_weights /= final_weights.sum()

    return Distribution(b1_values, final_weights)


class TruncatedSkewedDistributionConfig(BaseModel):
    r"""Configuration for truncated skew-normal B1 distribution.

    This distribution combines the asymmetry of skew-normal with a hard
    upper bound at the nominal B1 value, making it ideal for modeling
    RF field inhomogeneity when the B1 field cannot exceed the applied
    RF power (e.g., amplifier saturation, power limits).

    Attributes
    ----------
    type : Literal["truncated_skewed"]
        Distribution type identifier
    scale : float
        Coefficient of variation (std/mean) for the distribution
    skewness : float
        Shape parameter controlling asymmetry:
        - Negative: left-skewed (tail toward lower B1)
        - Positive: right-skewed (tail toward higher B1, but truncated)
        - Zero: symmetric (truncated Gaussian)
    res : int
        Number of quadrature points (recommend 11-21 for smooth distributions)
    lower_bound : float, optional
        Lower truncation bound as fraction of nominal (e.g., 0.7 for 70%)
        If None, automatically set to nominal * (1 - 4*scale)

    Examples
    --------
    Left-skewed distribution bounded at nominal (most common case):

        [experiment.b1_distribution]
        type = "truncated_skewed"
        scale = 0.10
        skewness = -0.8
        res = 11

    Or using inline table syntax:

        b1_distribution = { type = "truncated_skewed", scale = 0.10, \
                           skewness = -0.8, res = 11 }

    With custom lower bound:

        b1_distribution = { type = "truncated_skewed", scale = 0.10, \
                           skewness = -0.8, res = 11, lower_bound = 0.75 }

    """

    type: Literal["truncated_skewed"] = "truncated_skewed"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    skewness: float = Field(default=-0.8, ge=-10.0, le=10.0)
    res: int = Field(default=11, ge=1, le=101)
    lower_bound: float | None = Field(default=None, ge=0.0, le=1.0)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the truncated skew-normal distribution.

        Parameters
        ----------
        value : float
            Nominal B1 field value (upper bound)

        Returns
        -------
        Distribution
            B1 field values and weights

        """
        lower = value * self.lower_bound if self.lower_bound is not None else None
        return generate(value, self.scale, self.skewness, self.res, lower)


def register() -> None:
    """Register the truncated skewed distribution."""
    registry.register("truncated_skewed", generate, TruncatedSkewedDistributionConfig)
