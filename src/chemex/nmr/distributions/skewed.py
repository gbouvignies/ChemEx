"""Skewed B1 distribution using the Azzalini skew-normal law."""

from __future__ import annotations

from typing import Any, Literal

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions._quadrature import (
    MAX_SKEWNESS,
    build_skewed_or_gaussian,
    skewness_to_alpha_and_omega,
)

_skewness_to_alpha_and_omega = skewness_to_alpha_and_omega


def generate(
    value: float,
    scale: float,
    skewness: float,
    res: int = 21,
    *,
    center: Literal["mode", "mean"] = "mode",
) -> Distribution:
    """Generate a positive-conditioned Azzalini skew-normal B1 distribution.

    Models an Azzalini skew-normal distribution conditioned on ``B1 >= 0``,
    parameterized by parent standard deviation ``scale * value``, third
    standardized moment ``skewness`` in ``(-MAX_SKEWNESS, MAX_SKEWNESS)``,
    and ``center`` placing either the parent ``"mode"`` or parent ``"mean"`` at
    ``value``. Exact ``skewness == 0.0`` reduces identically to ``gaussian``.
    """
    return build_skewed_or_gaussian(
        value,
        scale,
        skewness=skewness,
        center=center,
        res=res,
    )


def __getattr__(name: str) -> Any:
    if name == "SkewedDistributionConfig":
        from chemex.configuration.b1_config import SkewedDistributionConfig

        return SkewedDistributionConfig
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)


__all__ = ["MAX_SKEWNESS", "_skewness_to_alpha_and_omega", "generate"]
