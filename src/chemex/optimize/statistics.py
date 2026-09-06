"""Shared counting semantics for fit statistics and covariance scaling."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class FitStatisticsCounts:
    """Counts that distinguish controlled and analytically fitted quantities."""

    residual_count: int
    controlled_coordinate_count: int
    profiled_normalization_count: int

    def __post_init__(self) -> None:
        if (
            min(
                self.residual_count,
                self.controlled_coordinate_count,
                self.profiled_normalization_count,
            )
            < 0
        ):
            raise ValueError("Fit-statistics counts cannot be negative")

    @property
    def estimated_parameter_count(self) -> int:
        """Return all quantities estimated from the observations."""
        return self.controlled_coordinate_count + self.profiled_normalization_count

    @property
    def effective_observation_count(self) -> int:
        """Return observations remaining after profiled normalizations."""
        return self.residual_count - self.profiled_normalization_count

    @property
    def residual_degrees_of_freedom(self) -> int:
        """Return nominal residual degrees of freedom, ``N - P - G``."""
        return self.residual_count - self.estimated_parameter_count

    @property
    def positive_residual_degrees_of_freedom(self) -> int | None:
        """Return residual degrees of freedom only when summary use is valid."""
        degrees_of_freedom = self.residual_degrees_of_freedom
        return degrees_of_freedom if degrees_of_freedom > 0 else None

    def reduced_chi_square(self, chi_square: float) -> float:
        """Return reduced chi-square, or NaN without positive residual DOF."""
        degrees_of_freedom = self.positive_residual_degrees_of_freedom
        return (
            math.nan if degrees_of_freedom is None else chi_square / degrees_of_freedom
        )
