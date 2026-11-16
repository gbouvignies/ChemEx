"""Truncated skew-normal B1 distribution with upper bound enforcement.

This distribution combines the asymmetry of the skew-normal distribution with
a hard upper bound at the nominal B1 value, making it suitable for modeling
RF field inhomogeneity when the B1 field cannot exceed the applied RF power.

Note: For backward compatibility in TOML configurations, the ``skewness``
field specifies the Azzalini skew-normal slant/shape parameter ``alpha`` in
``[-10, 10]`` of the pre-truncation density (not the third standardized
moment ``gamma_1`` used by ``type = "skewed"``).
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.polynomial.legendre import leggauss
from pydantic import Field
from scipy.optimize import minimize_scalar
from scipy.special import ndtr

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import DistributionConfig, registry


def _skewnorm_pdf(x: np.ndarray, xi: float, omega: float, alpha: float) -> np.ndarray:
    """Evaluate skew-normal PDF at points x."""
    z = (x - xi) / omega
    phi = np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)
    return (2.0 / omega) * phi * ndtr(alpha * z)


def _compute_mode_location(alpha: float) -> float:
    """Compute mode location for standard skew-normal SN(0, 1, alpha)."""
    if alpha == 0.0:
        return 0.0

    delta = alpha / np.sqrt(1.0 + alpha**2)
    mu_z = float(delta * np.sqrt(2.0 / np.pi))
    low, high = (0.0, mu_z) if mu_z >= 0.0 else (mu_z, 0.0)

    result = minimize_scalar(
        lambda z: -float(_skewnorm_pdf(np.array([z]), 0.0, 1.0, alpha)[0]),
        bounds=(low, high),
        method="bounded",
    )
    return float(result.x)


def generate(
    value: float,
    scale: float,
    skewness: float,
    res: int,
    lower_bound: float | None = None,
) -> Distribution:
    """Generate truncated skew-normal distribution with upper bound at nominal.

    The distribution is skew-normal with slant/shape parameter
    ``alpha = skewness`` and mode at ``value``, truncated to
    ``[lower_bound, value]`` with ``lower_bound >= 0`` so B1 quadrature
    points are always non-negative and never exceed ``value``.

    Parameters
    ----------
    value : float
        Nominal B1 field value (kHz or Hz) - also the hard upper bound
    scale : float
        Relative scale parameter (omega / value) of the pre-truncation density
    skewness : float
        Skew-normal slant/shape parameter ``alpha`` (not third-moment skewness):
        - Negative values: left-skewed (tail toward lower B1)
        - Positive values: right-skewed (tail toward higher B1, truncated at value)
        - Zero: symmetric (half/truncated Gaussian)
    res : int
        Number of quadrature points
    lower_bound : float, optional
        Lower truncation bound in the same units as ``value``. If ``None``,
        uses ``max(0.0, value * (1 - 4 * scale))``.

    Returns
    -------
    Distribution
        Non-negative B1 values and normalized quadrature weights.

    """
    if not np.isfinite(value) or value < 0.0:
        msg = f"Nominal B1 value must be finite and non-negative, got {value}"
        raise ValueError(msg)
    if np.isnan(scale) or scale < 0.0:
        msg = f"Scale must be non-negative, got {scale}"
        raise ValueError(msg)
    if not np.isfinite(skewness):
        msg = f"Skewness slant alpha must be finite, got {skewness}"
        raise ValueError(msg)

    if lower_bound is not None and (
        not np.isfinite(lower_bound) or lower_bound < 0.0 or lower_bound > value
    ):
        msg = f"lower_bound ({lower_bound}) must satisfy 0 <= lower_bound <= value ({value})"
        raise ValueError(msg)

    if value == 0.0 or scale == 0.0 or np.isinf(scale) or res <= 1:
        return Distribution(np.array([float(value)]), np.array([1.0]))

    alpha = float(skewness)

    if lower_bound is None:
        lower_bound = max(0.0, value * (1.0 - 4.0 * scale))

    upper_bound = float(value)
    if lower_bound == upper_bound:
        return Distribution(np.array([upper_bound]), np.array([1.0]))

    mode_rel = _compute_mode_location(alpha)
    omega = scale * value
    xi = value - omega * mode_rel

    nodes, weights = leggauss(res)
    b1_values = 0.5 * (upper_bound - lower_bound) * (nodes + 1.0) + lower_bound
    pdf_values = _skewnorm_pdf(b1_values, xi, omega, alpha)

    normalization = float(
        np.sum(pdf_values * weights * 0.5 * (upper_bound - lower_bound))
    )
    if normalization > 0.0:
        final_weights = (
            pdf_values * weights * 0.5 * (upper_bound - lower_bound) / normalization
        )
    else:
        final_weights = weights / np.sum(weights)

    final_weights /= final_weights.sum()
    return Distribution(b1_values, final_weights)


class TruncatedSkewedDistributionConfig(DistributionConfig):
    r"""Configuration for truncated skew-normal B1 distribution.

    Combines a skew-normal density with a hard upper bound at the nominal B1
    value and a non-negative lower bound. Note that ``skewness`` here is the
    skew-normal slant/shape parameter ``alpha`` in ``[-10, 10]`` of the
    pre-truncation density (retained under the ``skewness`` key for TOML
    compatibility), rather than the third standardized moment ``gamma_1``.

    Attributes
    ----------
    type : Literal["truncated_skewed"]
        Distribution type identifier
    scale : float
        Relative scale parameter (omega / nominal) of the underlying skew-normal
    skewness : float
        Skew-normal slant/shape parameter ``alpha`` in ``[-10, 10]`` controlling
        asymmetry of the pre-truncation density
    res : int
        Number of quadrature points (recommend 11-21 for smooth distributions)
    lower_bound : float, optional
        Lower truncation bound as fraction of nominal in ``[0, 1]`` (e.g., 0.7).
        If ``None``, automatically set to ``max(0, 1 - 4 * scale) * nominal``.

    """

    type: Literal["truncated_skewed"] = "truncated_skewed"
    scale: float = Field(default=0.1, ge=0.0, le=1.0)
    skewness: float = Field(
        default=-0.8,
        ge=-10.0,
        le=10.0,
        description=(
            "Skew-normal slant/shape parameter alpha in [-10, 10] of the "
            "pre-truncation density (distinct from third-moment skewness)"
        ),
    )
    res: int = Field(default=11, ge=1, le=101)
    lower_bound: float | None = Field(default=None, ge=0.0, le=1.0)

    def get_distribution(self, value: float) -> Distribution:
        """Generate the truncated skew-normal distribution."""
        lower = value * self.lower_bound if self.lower_bound is not None else None
        return generate(value, self.scale, self.skewness, self.res, lower)


def register() -> None:
    """Register the truncated skewed distribution."""
    registry.register("truncated_skewed", generate, TruncatedSkewedDistributionConfig)
