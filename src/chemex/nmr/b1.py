"""Runtime helpers for B1 fields and inhomogeneity models."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.gaussian import GaussianDistributionConfig
from chemex.nmr.distributions.registry import DistributionConfig


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
        return cls(tuple(scales.tolist()), tuple(weights.tolist()), distribution.dephasing)

    def get_distribution(self, value: float) -> Distribution:
        values = value * np.asarray(self.scales, dtype=float)
        weights = np.asarray(self.weights, dtype=float)
        return Distribution(values, weights, dephasing=self.dephasing)


B1DistributionModel = DistributionConfig | FixedDistributionModel | None


def clone_distribution_model(
    distribution: B1DistributionModel,
) -> B1DistributionModel:
    """Copy mutable distribution configs before storing them at runtime."""
    if distribution is None or isinstance(distribution, FixedDistributionModel):
        return distribution
    return distribution.model_copy(deep=True)


@dataclass(frozen=True)
class B1Profile:
    """Nominal B1 field together with the model used to build its distribution."""

    nominal: float
    distribution: B1DistributionModel = None

    @classmethod
    def gaussian(
        cls,
        nominal: float,
        *,
        scale: float = 0.0,
        res: int = 11,
    ) -> B1Profile:
        return cls(
            nominal=nominal,
            distribution=GaussianDistributionConfig(scale=scale, res=res),
        )

    def build_distribution(self) -> Distribution:
        if self.distribution is None:
            return Distribution(np.array([self.nominal]), np.array([1.0]))
        return self.distribution.get_distribution(self.nominal)

    def with_distribution(self, distribution: B1DistributionModel) -> B1Profile:
        return replace(self, distribution=clone_distribution_model(distribution))

    def with_nominal(self, nominal: float) -> B1Profile:
        return replace(self, nominal=nominal)

    def with_gaussian(
        self,
        *,
        scale: float | None = None,
        res: int | None = None,
    ) -> B1Profile:
        current_scale = float(getattr(self.distribution, "scale", 0.0))
        current_res = self.res if self.res > 1 else 11
        return type(self).gaussian(
            self.nominal,
            scale=current_scale if scale is None else scale,
            res=current_res if res is None else res,
        )

    @property
    def scale(self) -> float:
        return float(getattr(self.distribution, "scale", 0.0))

    @property
    def res(self) -> int:
        distribution = self.distribution
        if distribution is None:
            return 1
        res = getattr(distribution, "res", None)
        if isinstance(res, int):
            return res
        scales = getattr(distribution, "scales", None)
        if isinstance(scales, tuple):
            return len(scales)
        return 1
