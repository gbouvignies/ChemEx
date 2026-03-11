from __future__ import annotations

from collections.abc import Callable, Hashable, Iterable
from dataclasses import dataclass, field

from cachetools import LRUCache
from lmfit import Parameters as ParametersLF

from chemex.typing import Array


@dataclass
class ProfileEvaluator:
    maxsize: int = 5
    cache: LRUCache = field(init=False)

    def __post_init__(self) -> None:
        self.cache = LRUCache(maxsize=self.maxsize)

    def clear(self) -> None:
        self.cache.clear()

    def residuals(
        self,
        params: ParametersLF,
        *,
        param_ids: Iterable[str],
        data_revision: int,
        calculate: Callable[[ParametersLF], Array],
        exp: Array,
        err: Array,
        mask: Array,
    ) -> Array:
        key: tuple[Hashable, ...] = (
            *(params[param_id].value for param_id in param_ids),
            data_revision,
        )
        cached_residuals = self.cache.get(key)
        if cached_residuals is not None:
            return cached_residuals

        residuals = (calculate(params) - exp) / err
        masked_residuals = residuals[mask]
        self.cache[key] = masked_residuals
        return masked_residuals
