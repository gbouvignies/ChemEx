"""Dephasing mode for extreme B1 field inhomogeneity in continuous-wave CEST."""

from __future__ import annotations

from typing import Any

import numpy as np

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions._quadrature import validate_nominal_and_scale


def generate(value: float) -> Distribution:
    """Generate a single-point distribution with transverse dephasing enabled.

    For positive pulse durations, the propagator zeroes eigenmodes with a
    non-negligible imaginary eigenvalue while retaining exponential evolution
    for modes aligned along the effective field.
    """
    validate_nominal_and_scale(value, 0.0)
    return Distribution(
        values=np.array([float(value)], dtype=float),
        weights=np.array([1.0], dtype=float),
        dephasing=True,
    )


def __getattr__(name: str) -> Any:
    if name == "DephasingConfig":
        from chemex.configuration.b1_config import DephasingConfig

        return DephasingConfig
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
