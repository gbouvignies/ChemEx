"""Gaussian B1 distribution with positive-support Gauss-Legendre quadrature."""

from __future__ import annotations

from typing import Any

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions._quadrature import build_skewed_or_gaussian


def generate(
    value: float,
    scale: float,
    res: int = 15,
) -> Distribution:
    """Generate a positive-conditioned Gaussian B1 distribution.

    Models a Normal distribution with parent mean ``value`` and parent standard
    deviation ``scale * value`` conditioned on ``B1 >= 0``, integrated by
    Gauss-Legendre quadrature over its ``[1e-5, 1 - 1e-5]`` quantile support.
    """
    return build_skewed_or_gaussian(
        value,
        scale,
        skewness=0.0,
        center="mean",
        res=res,
    )


def __getattr__(name: str) -> Any:
    if name == "GaussianDistributionConfig":
        from chemex.configuration.b1_config import GaussianDistributionConfig

        return GaussianDistributionConfig
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
