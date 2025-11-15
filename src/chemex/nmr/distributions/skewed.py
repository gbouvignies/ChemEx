"""Skewed B1 distribution using skew-normal distribution.

Uses skew-normal distribution with configurable skewness parameter.
The mode of the distribution is positioned at the nominal B1 value.
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
    """Probability density function of skew-normal distribution.

    Parameters
    ----------
    x : np.ndarray
        Values at which to evaluate PDF
    xi : float
        Location parameter
    omega : float
        Scale parameter (must be > 0)
    alpha : float
        Shape parameter (controls skewness)

    Returns
    -------
    np.ndarray
        PDF values at x

    Notes
    -----
    PDF formula: f(x) = (2/omega) * phi((x-xi)/omega) * Phi(alpha*(x-xi)/omega)
    where phi is standard normal PDF and Phi is standard normal CDF

    """
    z = (x - xi) / omega

    # Standard normal PDF
    phi = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * z * z)

    # Standard normal CDF using error function
    phi_cdf = 0.5 * (1 + erf(alpha * z / np.sqrt(2)))

    # Skew-normal PDF
    pdf = (2 / omega) * phi * phi_cdf

    return pdf


def _compute_mode_location(alpha: float) -> float:
    """Compute the mode location of standard skew-normal distribution.

    For a skew-normal with xi=0 and omega=1, find the location of the mode
    (maximum of the PDF). This must be solved numerically.

    Parameters
    ----------
    alpha : float
        Shape parameter (controls skewness)

    Returns
    -------
    float
        Mode location relative to xi (in units of omega)

    """
    # For symmetric case (alpha=0), mode is at 0
    if abs(alpha) < 0.01:
        return 0.0

    # For skewed case, find the maximum of the PDF numerically
    # Search in a reasonable range
    search_range = (-3, 3) if alpha > 0 else (-3, 3)

    result = minimize_scalar(
        lambda x: -_skewnorm_pdf(np.array([x]), 0.0, 1.0, alpha)[0],
        bounds=search_range,
        method="bounded",
    )

    return result.x


def generate(
    value: float,
    scale: float,
    skewness: float,
    res: int,
) -> Distribution:
    """Generate skewed B1 field distribution for ChemEx integration.

    Uses skew-normal distribution with **mode at the nominal value**, specified
    width, and configurable skewness. Falls back to Hermite when |skewness| < 0.01.

    Parameters
    ----------
    value : float
        Nominal B1 field value (mode of the distribution, in kHz or Hz)
    scale : float
        Width parameter as fraction of value (e.g., 0.10 = 10% of nominal)
        Approximately corresponds to std deviation for small skewness
    skewness : float
        Skewness parameter (dimensionless)
        - skewness < 0: left tail extended (negative skew)
        - skewness ≈ 0: symmetric (Gaussian)
        - skewness > 0: right tail extended (positive skew)
        Typical range: -1.5 to +1.5
    res : int
        Number of integration points (recommended: 11-21)

    Returns
    -------
    Distribution
        B1 values and weights for numerical integration

    Notes
    -----
    The distribution is parameterized so that the **mode** (peak of the PDF)
    is positioned exactly at the nominal B1 value. This ensures that the most
    probable B1 value matches the target.

    """
    # Handle edge cases
    if scale in (0.0, np.inf) or res <= 1:
        return Distribution(np.array([value]), np.array([1.0]))

    # For nearly symmetric distributions, use Hermite (more efficient)
    if abs(skewness) < 0.01:
        from chemex.nmr.distributions.hermite import generate as hermite_generate

        return hermite_generate(value, scale, res)

    # Convert skewness to delta parameter
    if abs(skewness) < 2.0:
        delta = np.sign(skewness) * np.power(
            abs(skewness) / (0.995 * ((4 - np.pi) / 2)), 1 / 3
        )
    else:
        # Asymptotic approximation for extreme skewness
        delta = np.sign(skewness) * (1.0 - 0.5 / abs(skewness))

    # Ensure delta is in valid range
    delta = np.clip(delta, -0.995, 0.995)

    # Convert delta to alpha (shape parameter)
    alpha = delta / np.sqrt(1 - delta**2)

    # Estimate omega (scale parameter) from width
    skew_factor = 1.0 + 0.15 * abs(skewness)
    omega = scale * skew_factor

    # Compute mode location for standard skew-normal (xi=0, omega=1)
    mode_rel = _compute_mode_location(alpha)

    # Set xi so that mode is at relative position 1.0 (nominal value)
    # mode_actual = xi + omega * mode_rel = 1.0
    # Therefore: xi = 1.0 - omega * mode_rel
    xi = 1.0 - omega * mode_rel

    # Get Gauss-Legendre quadrature nodes and weights
    nodes, weights_gl = leggauss(res)

    # Define integration bounds (capture 99.9% of probability)
    rel_b1_min = xi - 4.0 * omega
    rel_b1_max = xi + 4.0 * omega

    # Ensure bounds are positive and reasonable
    rel_b1_min = max(0.3, rel_b1_min)  # Don't go below 30% of nominal
    rel_b1_max = min(2.0, rel_b1_max)  # Don't exceed 200% of nominal

    # Transform Gauss-Legendre nodes from [-1, 1] to [rel_b1_min, rel_b1_max]
    rel_b1_values = 0.5 * (rel_b1_max - rel_b1_min) * nodes + 0.5 * (
        rel_b1_max + rel_b1_min
    )

    # Calculate skew-normal PDF at each relative B1 value
    pdf_values = _skewnorm_pdf(rel_b1_values, xi, omega, alpha)

    # Weight by PDF and Jacobian of transformation
    weights = 0.5 * (rel_b1_max - rel_b1_min) * weights_gl * pdf_values

    # Normalize weights to sum to 1.0
    weights = weights / np.sum(weights)

    # Convert relative values to absolute B1 values
    b1_values = rel_b1_values * value

    return Distribution(b1_values, weights)


class SkewedDistributionConfig(BaseModel):
    """Configuration for skewed B1 distribution.

    Attributes
    ----------
    type : Literal["skewed"]
        Distribution type identifier
    scale : float
        Width parameter as fraction of nominal value
    res : int
        Number of integration points
    skewness : float
        Skewness parameter

    """

    type: Literal["skewed"] = "skewed"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    res: int = Field(default=11, ge=1, le=51)
    skewness: float = Field(default=0.0, ge=-2.0, le=2.0)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the skewed distribution."""
        return generate(value, self.scale, self.skewness, self.res)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("skewed", generate, SkewedDistributionConfig)
