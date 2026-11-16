"""Skewed B1 distribution using the Azzalini skew-normal distribution.

Uses an Azzalini skew-normal distribution parameterized by its theoretical
standardized third cumulant (skewness, |gamma_1| < 0.99527...) with the mode
positioned at the nominal B1 value.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.polynomial.legendre import leggauss
from pydantic import Field, field_validator
from scipy.optimize import minimize_scalar
from scipy.special import ndtr

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import DistributionConfig, registry

# Theoretical maximum skewness magnitude for an Azzalini skew-normal random
# variable: ((4 - pi) / 2) * (2 / (pi - 2))**(3 / 2) ≈ 0.9952717464311563.
MAX_SKEWNESS: float = float(0.5 * (4.0 - np.pi) * (2.0 / (np.pi - 2.0)) ** 1.5)


def _validate_skewness(skewness: float) -> float:
    """Validate that skewness lies strictly inside the skew-normal domain."""
    if not np.isfinite(skewness) or abs(skewness) >= MAX_SKEWNESS:
        msg = (
            f"Skewness {skewness} is outside the valid Azzalini skew-normal range "
            f"(-{MAX_SKEWNESS:.6f}, {MAX_SKEWNESS:.6f})"
        )
        raise ValueError(msg)
    return float(skewness)


def _skewness_to_alpha_and_omega(skewness: float, scale: float) -> tuple[float, float]:
    """Invert Azzalini's skewness formula analytically to (alpha, omega).

    For X ~ SN(xi, omega, alpha) with delta = alpha / sqrt(1 + alpha^2):
        gamma_1 = ((4 - pi) / 2) * (delta * sqrt(2 / pi))^3 / (1 - 2*delta^2/pi)^(3/2)
        Var(X)  = omega^2 * (1 - 2*delta^2/pi) = scale^2
    Letting c = (2 * |gamma_1| / (4 - pi))^(2 / 3) gives:
        2*delta^2/pi = c / (1 + c),   1 - 2*delta^2/pi = 1 / (1 + c),
        1 - delta^2  = (1 - ((pi - 2) / 2) * c) / (1 + c) > 0,
        alpha = sign(gamma_1) * sqrt((pi / 2) * c / (1 - ((pi - 2) / 2) * c)),
        omega = scale * sqrt(1 + c).
    """
    _validate_skewness(skewness)
    if skewness == 0.0:
        return 0.0, float(scale)
    c = (2.0 * abs(skewness) / (4.0 - np.pi)) ** (2.0 / 3.0)
    denom = 1.0 - 0.5 * (np.pi - 2.0) * c
    alpha = float(np.copysign(np.sqrt(0.5 * np.pi * c / denom), skewness))
    omega = float(scale * np.sqrt(1.0 + c))
    return alpha, omega


def _skewnorm_pdf(x: np.ndarray, xi: float, omega: float, alpha: float) -> np.ndarray:
    """Probability density function of the skew-normal distribution."""
    z = (x - xi) / omega
    phi = (1.0 / np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * z * z)
    return (2.0 / omega) * phi * ndtr(alpha * z)


def _compute_mode_location(alpha: float) -> float:
    """Compute the mode of SN(0, 1, alpha), bracketed between 0 and E[Z]."""
    if alpha == 0.0:
        return 0.0

    delta = alpha / np.sqrt(1.0 + alpha * alpha)
    mu_z = float(delta * np.sqrt(2.0 / np.pi))
    low, high = (0.0, mu_z) if mu_z >= 0.0 else (mu_z, 0.0)

    result = minimize_scalar(
        lambda x: -float(_skewnorm_pdf(np.array([x]), 0.0, 1.0, alpha)[0]),
        bounds=(low, high),
        method="bounded",
    )
    return float(result.x)


def generate(
    value: float,
    scale: float,
    skewness: float,
    res: int,
) -> Distribution:
    """Generate a skewed B1 field distribution using Azzalini's skew-normal law.

    Uses a skew-normal distribution with its mode positioned at ``value``,
    relative standard deviation ``scale``, and third standardized moment
    ``skewness`` strictly inside ``(-MAX_SKEWNESS, MAX_SKEWNESS)``
    (``|skewness| < 0.99527...``). Exact ``skewness == 0.0`` delegates to
    Gauss-Hermite quadrature; nonzero skewness uses smooth Gauss-Legendre
    quadrature truncated below ``B1 = 0``.
    """
    if not np.isfinite(value) or value < 0.0:
        msg = f"Nominal B1 value must be finite and non-negative, got {value}"
        raise ValueError(msg)
    if np.isnan(scale) or scale < 0.0:
        msg = f"Scale must be non-negative, got {scale}"
        raise ValueError(msg)
    _validate_skewness(skewness)

    if value == 0.0 or scale in (0.0, np.inf) or res <= 1:
        return Distribution(np.array([float(value)]), np.array([1.0]))

    if skewness == 0.0:
        from chemex.nmr.distributions.hermite import generate as hermite_generate

        return hermite_generate(value, scale, res)

    alpha, omega = _skewness_to_alpha_and_omega(skewness, scale)
    mode_rel = _compute_mode_location(alpha)

    # Place the skew-normal mode at relative B1 = 1.0 (the nominal value).
    xi = 1.0 - omega * mode_rel

    # Adapt quadrature support to the skew-normal tail decay while enforcing B1 >= 0.
    tail_factor = float(np.sqrt(1.0 + alpha * alpha))
    z_min = -5.0 / tail_factor if alpha > 0.0 else -5.0
    z_max = 5.0 / tail_factor if alpha < 0.0 else 5.0
    rel_b1_min = max(0.0, xi + z_min * omega)
    rel_b1_max = xi + z_max * omega

    nodes, weights_gl = leggauss(res)
    rel_b1_values = 0.5 * (rel_b1_max - rel_b1_min) * nodes + 0.5 * (
        rel_b1_max + rel_b1_min
    )

    pdf_values = _skewnorm_pdf(rel_b1_values, xi, omega, alpha)
    weights = 0.5 * (rel_b1_max - rel_b1_min) * weights_gl * pdf_values
    weights_sum = float(np.sum(weights))
    if weights_sum > 0.0:
        weights = weights / weights_sum
    else:
        weights = weights_gl / np.sum(weights_gl)

    return Distribution(rel_b1_values * value, weights)


class SkewedDistributionConfig(DistributionConfig):
    """Configuration for Azzalini skew-normal B1 distribution.

    Attributes
    ----------
    type : Literal["skewed"]
        Distribution type identifier
    scale : float
        Standard deviation as fraction of nominal value
    res : int
        Number of integration points
    skewness : float
        Theoretical skewness (third standardized moment), strictly inside
        ``(-0.995271..., +0.995271...)``

    """

    type: Literal["skewed"] = "skewed"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    res: int = Field(default=11, ge=1, le=51)
    skewness: float = Field(
        default=0.0,
        description=(
            "Theoretical third standardized moment of the Azzalini skew-normal "
            "distribution (|skewness| < 0.995272)"
        ),
    )

    @field_validator("skewness")
    @classmethod
    def _check_skewness_bounds(cls, value: float) -> float:
        return _validate_skewness(value)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the skewed distribution."""
        return generate(value, self.scale, self.skewness, self.res)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("skewed", generate, SkewedDistributionConfig)
