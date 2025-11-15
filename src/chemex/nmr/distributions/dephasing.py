"""Dephasing mode for B1 field inhomogeneity.

This simulates extreme B1 inhomogeneity by dephasing all magnetization
that is not along the effective field vector. This is not a true distribution
but rather a special calculation mode used in some CEST experiments.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.registry import registry


def generate(
    value: float,
) -> Distribution:
    """Generate a 'distribution' that triggers dephasing mode.

    In dephasing mode, the propagator calculation adjusts eigenvalues
    to effectively dephase all coherences except those along the
    effective field (eigenvalues with non-zero imaginary parts are
    scaled to cause rapid dephasing).

    Parameters
    ----------
    value : float
        Nominal B1 field value (not used, but required for API consistency)

    Returns
    -------
    Distribution
        Single-point distribution with dephasing flag set to True

    Notes
    -----
    This is equivalent to the old `b1_inh_scale = inf` hack, but
    implemented as an explicit distribution type for clarity.

    """
    return Distribution(
        values=np.array([value]),
        weights=np.array([1.0]),
        dephasing=True,
    )


class DephasingConfig(BaseModel):
    """Configuration for dephasing mode (no additional parameters needed)."""

    type: Literal["dephasing"]

    def get_distribution(self, value: float) -> Distribution:
        """Generate the dephasing distribution."""
        return generate(value)


def register() -> None:
    """Register this distribution with the registry."""
    registry.register("dephasing", generate, DephasingConfig)
