"""Pure runtime models for B1 fields and inhomogeneity profiles."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal, Protocol, runtime_checkable

import numpy as np

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions import custom as custom_dist
from chemex.nmr.distributions import dephasing as dephasing_dist
from chemex.nmr.distributions import gaussian as gaussian_dist
from chemex.nmr.distributions import homogeneous as homogeneous_dist
from chemex.nmr.distributions import skewed as skewed_dist


@dataclass(frozen=True)
class HomogeneousB1Model:
    """Single-point homogeneous B1 field."""

    def get_distribution(self, value: float) -> Distribution:
        return homogeneous_dist.generate(value)


@dataclass(frozen=True)
class DephasingB1Model:
    """Single-point B1 model with complete transverse dephasing."""

    def get_distribution(self, value: float) -> Distribution:
        return dephasing_dist.generate(value)


@dataclass(frozen=True)
class GaussianB1Model:
    """Positive-support Gaussian B1 model with fixed quadrature node count."""

    scale: float = 0.0
    res: int = 1

    def get_distribution(self, value: float) -> Distribution:
        return gaussian_dist.generate(value, self.scale, self.res)


@dataclass(frozen=True)
class SkewedB1Model:
    """Positive-support Azzalini skew-normal B1 model with fixed node count."""

    scale: float
    skewness: float
    center: Literal["mode", "mean"]
    res: int

    def get_distribution(self, value: float) -> Distribution:
        return skewed_dist.generate(
            value,
            self.scale,
            self.skewness,
            self.res,
            center=self.center,
        )


@dataclass(frozen=True)
class CustomB1Model:
    """Discrete empirical B1 distribution scaled relative to nominal B1."""

    scales: tuple[float, ...]
    weights: tuple[float, ...]

    def get_distribution(self, value: float) -> Distribution:
        return custom_dist.generate(value, self.scales, self.weights)


@dataclass(frozen=True)
class FixedDistributionModel:
    """Rescale a sampled distribution when the nominal B1 field changes."""

    scales: tuple[float, ...]
    weights: tuple[float, ...]
    dephasing: bool = False

    @classmethod
    def from_distribution(
        cls,
        distribution: Distribution,
        nominal: float,
    ) -> FixedDistributionModel:
        values = np.asarray(distribution.values, dtype=float)
        weights = np.asarray(distribution.weights, dtype=float)
        scales = np.ones_like(values) if nominal == 0.0 else values / nominal
        return cls(
            tuple(scales.tolist()), tuple(weights.tolist()), distribution.dephasing
        )

    def get_distribution(self, value: float) -> Distribution:
        values = value * np.asarray(self.scales, dtype=float)
        weights = np.asarray(self.weights, dtype=float)
        return Distribution(values, weights, dephasing=self.dephasing)


B1RuntimeModel = (
    HomogeneousB1Model
    | DephasingB1Model
    | GaussianB1Model
    | SkewedB1Model
    | CustomB1Model
    | FixedDistributionModel
    | None
)


@runtime_checkable
class SupportsToB1Model(Protocol):
    """Configuration objects that lower to a pure runtime B1 model."""

    def to_b1_model(
        self,
        *,
        b1_eff_hz: float,
        t_sat_s: float = 0.5,
    ) -> B1RuntimeModel: ...


B1DistributionModel = B1RuntimeModel | SupportsToB1Model


def resolve_runtime_b1_model(
    distribution: B1DistributionModel,
    *,
    nominal: float,
    t_sat_s: float = 0.5,
) -> B1RuntimeModel:
    """Normalize a runtime B1 model or configuration into a pure frozen model."""
    if isinstance(distribution, SupportsToB1Model):
        return distribution.to_b1_model(b1_eff_hz=nominal, t_sat_s=t_sat_s)
    return distribution


@dataclass(frozen=True)
class B1Profile:
    """Nominal B1 field together with the immutable model used to build it."""

    nominal: float
    distribution: B1RuntimeModel = None

    @classmethod
    def gaussian(
        cls,
        nominal: float,
        *,
        scale: float = 0.0,
        res: int = 11,
    ) -> B1Profile:
        resolved_res = 1 if scale == 0.0 else res
        return cls(
            nominal=nominal,
            distribution=GaussianB1Model(scale=scale, res=resolved_res),
        )

    def build_distribution(self) -> Distribution:
        if self.distribution is None:
            return Distribution(np.array([self.nominal]), np.array([1.0]))
        return self.distribution.get_distribution(self.nominal)

    def with_distribution(self, distribution: B1DistributionModel) -> B1Profile:
        runtime_model = resolve_runtime_b1_model(distribution, nominal=self.nominal)
        return replace(self, distribution=runtime_model)

    def with_nominal(self, nominal: float) -> B1Profile:
        return replace(self, nominal=nominal)
