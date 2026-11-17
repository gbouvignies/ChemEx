"""Pure numerical quadrature for Gaussian and skew-normal B1 distributions.

Both distributions represent a parent Normal or Azzalini Skew-Normal law
conditioned on physical non-negative support ``B1 >= 0``:

    p_+(b) = p(b) / (1 - F(0)),   b >= 0,

and integrated over the positive-conditioned quantile bracket
``[F_+^{-1}(eps), F_+^{-1}(1 - eps)]`` with ``eps = 1e-5`` via standardized
Gauss-Legendre quadrature.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import log_ndtr, ndtr, ndtri
from scipy.stats import skewnorm

from chemex.nmr.constants import Distribution

QUAD_TAIL_EPS: float = 1e-5
AUTO_RES_MIN: int = 13
AUTO_RES_MAX: int = 651

# Theoretical maximum skewness magnitude for an Azzalini skew-normal variable:
# ((4 - pi) / 2) * (2 / (pi - 2))**(3 / 2) ~= 0.9952717464311563.
MAX_SKEWNESS: float = float(0.5 * (4.0 - np.pi) * (2.0 / (np.pi - 2.0)) ** 1.5)

_HALF_LOG_2PI = 0.5 * math.log(2.0 * math.pi)
_SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)


def validate_nominal_and_scale(value: float, scale: float) -> None:
    """Validate that nominal B1 and relative scale are finite and non-negative."""
    if not math.isfinite(value) or value < 0.0:
        msg = f"Nominal B1 value must be finite and non-negative, got {value}"
        raise ValueError(msg)
    if not math.isfinite(scale) or scale < 0.0:
        msg = f"Scale must be finite and non-negative, got {scale}"
        raise ValueError(msg)


def validate_skewness(skewness: float) -> float:
    """Validate that skewness lies strictly inside the Azzalini domain."""
    if not math.isfinite(skewness) or abs(skewness) >= MAX_SKEWNESS:
        msg = (
            f"Skewness {skewness} is outside the valid Azzalini skew-normal range "
            f"(-{MAX_SKEWNESS:.6f}, {MAX_SKEWNESS:.6f})"
        )
        raise ValueError(msg)
    return float(skewness)


def validate_center(center: str) -> Literal["mode", "mean"]:
    """Validate the centering convention for a skewed distribution."""
    if center == "mode":
        return "mode"
    if center == "mean":
        return "mean"
    msg = f"Center must be 'mode' or 'mean', got {center!r}"
    raise ValueError(msg)


def validate_res(res: int) -> int:
    """Validate an explicit quadrature node count."""
    if isinstance(res, bool) or not isinstance(res, int) or res < 1:
        msg = f"Resolution 'res' must be a positive integer, got {res}"
        raise ValueError(msg)
    return res


def skewness_to_alpha_and_sigma_z(skewness: float) -> tuple[float, float]:
    """Invert Azzalini's skewness gamma_1 to (alpha, sigma_Z)."""
    validate_skewness(skewness)
    if skewness == 0.0:
        return 0.0, 1.0
    c = (2.0 * abs(skewness) / (4.0 - math.pi)) ** (2.0 / 3.0)
    sigma_z = 1.0 / math.sqrt(1.0 + c)
    denom = 1.0 - 0.5 * (math.pi - 2.0) * c
    alpha = math.copysign(math.sqrt(0.5 * math.pi * c / denom), skewness)
    return alpha, sigma_z


def skewness_to_alpha_and_omega(skewness: float, scale: float) -> tuple[float, float]:
    """Invert Azzalini's skewness gamma_1 and relative scale to (alpha, omega)."""
    alpha, sigma_z = skewness_to_alpha_and_sigma_z(skewness)
    return alpha, float(scale) / sigma_z


def compute_skewnorm_mode_z(alpha: float) -> float:
    """Compute the mode of standard skew-normal SN(0, 1, alpha) to ~1e-16."""
    if alpha == 0.0:
        return 0.0

    a = abs(alpha)
    delta = a / math.sqrt(1.0 + a * a)
    mu_z = delta * _SQRT_2_OVER_PI

    if a < 1.0e-3:
        z = mu_z * (1.0 - (2.0 / math.pi) * a * a)
    else:
        sigma_z = math.sqrt(max(1.0 - mu_z * mu_z, 1.0e-30))
        gamma1 = 0.5 * (4.0 - math.pi) * (mu_z**3) / (sigma_z**3)
        z = mu_z - 0.5 * gamma1 * sigma_z - 0.5 * math.exp(-2.0 * math.pi / a)

    low, high = 0.0, mu_z
    if not (low < z < high):
        z = 0.5 * mu_z

    # Solve f(z) = -z + a * phi(a*z) / Phi(a*z) = 0 with safeguarded Newton steps.
    for _ in range(25):
        az = a * z
        ratio = math.exp(
            -0.5 * az * az - _HALF_LOG_2PI - float(log_ndtr(np.float64(az)))
        )
        f_val = -z + a * ratio
        if f_val == 0.0:
            break
        if f_val > 0.0:
            low = z
        else:
            high = z
        deriv = -1.0 - (a * a) * ratio * (az + ratio)
        step = f_val / deriv
        z_next = z - step
        if abs(step) <= 1.0e-15 or abs(high - low) <= 1.0e-15:
            if low <= z_next <= high:
                z = z_next
            break
        if not (low < z_next < high):
            z_next = 0.5 * (low + high)
        z = z_next

    return math.copysign(z, alpha)


def _standardized_bracket(
    scale: float,
    skewness: float,
    center: Literal["mode", "mean"],
) -> tuple[float, float, float, float, float]:
    """Return (alpha, sigma_z, m_anchor, z_lo, z_hi) for B1 >= 0 support."""
    if skewness == 0.0:
        z0 = -1.0 / scale
        f0 = float(ndtr(np.float64(z0)))
        s0 = float(ndtr(np.float64(-z0)))
        z_lo = float(ndtri(np.float64(f0 + QUAD_TAIL_EPS * s0)))
        z_hi = -float(ndtri(np.float64(QUAD_TAIL_EPS * s0)))
        return 0.0, 1.0, 0.0, z_lo, z_hi

    alpha, sigma_z = skewness_to_alpha_and_sigma_z(skewness)
    if center == "mean":
        delta = alpha / math.sqrt(1.0 + alpha * alpha)
        m_anchor = delta * _SQRT_2_OVER_PI
    else:
        m_anchor = compute_skewnorm_mode_z(alpha)

    z0 = m_anchor - sigma_z / scale
    f0 = float(skewnorm.cdf(z0, alpha))
    s0 = float(skewnorm.sf(z0, alpha))
    z_lo = float(skewnorm.ppf(f0 + QUAD_TAIL_EPS * s0, alpha))
    z_hi = float(skewnorm.isf(QUAD_TAIL_EPS * s0, alpha))
    return alpha, sigma_z, m_anchor, z_lo, z_hi


def auto_quadrature_res(
    span_s: float,
    b1_eff_hz: float,
    t_sat_s: float,
    alpha: float = 0.0,
) -> int:
    """Choose odd Gauss-Legendre resolution for <= 1e-4 CEST integration error."""
    cyc_q = abs(span_s * b1_eff_hz * t_sat_s)
    abs_alpha = abs(alpha)
    d_rabi = 0.5 * math.pi * cyc_q + 3.5 * (cyc_q ** (1.0 / 3.0))
    d_skew = (
        4.0 * math.sqrt(abs_alpha) + 7.5 * abs_alpha / (1.0 + abs_alpha)
        if abs_alpha > 0.0
        else 0.0
    )
    min_nodes = math.ceil(max(float(AUTO_RES_MIN), 5.5 + math.hypot(d_rabi, d_skew)))
    n_odd = min_nodes if min_nodes % 2 == 1 else min_nodes + 1
    return min(AUTO_RES_MAX, n_odd)


def resolve_gaussian_res(
    scale: float,
    *,
    b1_eff_hz: float,
    t_sat_s: float,
    res: int | None = None,
) -> int:
    """Resolve the fixed integer resolution for a Gaussian B1 distribution."""
    validate_nominal_and_scale(max(0.0, float(b1_eff_hz)), scale)
    if scale == 0.0 or res == 1:
        return 1
    if res is not None:
        return validate_res(res)
    _, sigma_z, _, z_lo, z_hi = _standardized_bracket(scale, 0.0, "mean")
    span_s = (scale / sigma_z) * (z_hi - z_lo)
    return auto_quadrature_res(span_s, float(b1_eff_hz), float(t_sat_s), 0.0)


def resolve_skewed_res(
    scale: float,
    skewness: float,
    center: Literal["mode", "mean"],
    *,
    b1_eff_hz: float,
    t_sat_s: float,
    res: int | None = None,
) -> int:
    """Resolve the fixed integer resolution for a skew-normal B1 distribution."""
    validate_nominal_and_scale(max(0.0, float(b1_eff_hz)), scale)
    validate_skewness(skewness)
    validate_center(center)
    if scale == 0.0 or res == 1:
        return 1
    if res is not None:
        return validate_res(res)
    alpha, sigma_z, _, z_lo, z_hi = _standardized_bracket(scale, skewness, center)
    span_s = (scale / sigma_z) * (z_hi - z_lo)
    return auto_quadrature_res(span_s, float(b1_eff_hz), float(t_sat_s), alpha)


def build_skewed_or_gaussian(
    value: float,
    scale: float,
    *,
    skewness: float = 0.0,
    center: Literal["mode", "mean"] = "mean",
    res: int,
) -> Distribution:
    """Build non-negative Gauss-Legendre quadrature for Gaussian/skew-normal B1."""
    validate_nominal_and_scale(value, scale)
    validate_skewness(skewness)
    checked_center = validate_center(center)
    checked_res = validate_res(res)

    if value == 0.0 or scale == 0.0 or checked_res == 1:
        return Distribution(np.array([float(value)]), np.array([1.0]))

    alpha, sigma_z, m_anchor, z_lo, z_hi = _standardized_bracket(
        scale,
        skewness,
        checked_center,
    )
    x_gl, v_gl = leggauss(checked_res)
    z_nodes = 0.5 * (z_hi - z_lo) * x_gl + 0.5 * (z_hi + z_lo)
    log_w = np.log(v_gl) - 0.5 * z_nodes * z_nodes
    if alpha != 0.0:
        log_w = log_w + log_ndtr(alpha * z_nodes)
    log_w = log_w - np.max(log_w)
    raw_weights = np.exp(log_w)
    weights = raw_weights / np.sum(raw_weights)
    u_nodes = (z_nodes - m_anchor) / sigma_z
    scales = np.maximum(0.0, 1.0 + scale * u_nodes)
    return Distribution(float(value) * scales, weights)
