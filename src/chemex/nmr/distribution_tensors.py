"""Helpers for broadcasting weighted NMR distributions over Liouvillian tensors."""

from __future__ import annotations

from dataclasses import dataclass

from chemex.nmr.constants import Distribution
from chemex.typing import Array


@dataclass(frozen=True, slots=True)
class DistributionAxis:
    """Distribution values and weights reshaped for tensor broadcasting."""

    distribution: Distribution
    values: Array
    weights: Array

    @classmethod
    def create(
        cls,
        distribution: Distribution,
        *,
        ndim: int,
    ) -> DistributionAxis:
        if ndim < 1:
            msg = f"'ndim' must be positive, got {ndim}"
            raise ValueError(msg)
        shape = (-1,) + (1,) * (ndim - 1)
        return cls(
            distribution=distribution,
            values=distribution.values.reshape(shape),
            weights=distribution.weights.reshape(shape),
        )


@dataclass(frozen=True, slots=True)
class B1DistributionState:
    """Broadcasted B1 values, weights, and scaled RF Liouvillian terms."""

    axis: DistributionAxis
    l_x: Array
    l_y: Array

    @classmethod
    def build(
        cls,
        distribution: Distribution,
        *,
        matrix_x: Array | float = 0.0,
        matrix_y: Array | float = 0.0,
    ) -> B1DistributionState:
        axis = DistributionAxis.create(distribution, ndim=3)
        return cls(
            axis=axis,
            l_x=matrix_x * axis.values,
            l_y=matrix_y * axis.values,
        )


@dataclass(frozen=True, slots=True)
class JeffDistributionState:
    """Broadcasted effective-coupling values, weights, and scaled Liouvillian term."""

    axis: DistributionAxis
    liouvillian: Array

    @classmethod
    def build(
        cls,
        distribution: Distribution,
        *,
        matrix: Array | float = 0.0,
    ) -> JeffDistributionState:
        axis = DistributionAxis.create(distribution, ndim=4)
        return cls(
            axis=axis,
            liouvillian=matrix * axis.values,
        )
