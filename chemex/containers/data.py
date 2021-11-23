from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from dataclasses import field
from random import choices

import numpy as np
from numpy.typing import NDArray

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayInt = NDArray[np.int_]
NDArrayBool = NDArray[np.bool_]


def get_scale(exp: np.ndarray, err: np.ndarray, calc: np.ndarray) -> float:
    try:
        calc_err2 = calc / err**2
        return sum(exp * calc_err2) / sum(calc * calc_err2)
    except ZeroDivisionError:
        return 1.0


@dataclass
class Data:
    exp: NDArrayFloat
    err: NDArrayFloat
    metadata: np.ndarray = field(default_factory=lambda: np.array([]))
    scale: float = 1.0
    size: int = field(init=False)
    calc: NDArrayFloat = field(init=False)
    mask: NDArrayBool = field(init=False)
    refs: NDArrayBool = field(init=False)

    def __post_init__(self):
        self.size = self.exp.size
        self.calc = np.full_like(self.exp, 1e32)
        self.mask = np.full_like(self.exp, True, dtype=np.bool_)
        self.refs = np.full_like(self.exp, False, dtype=np.bool_)

    def scale_calc(self) -> None:
        mask = self.mask
        scale = get_scale(self.exp[mask], self.err[mask], self.calc[mask])
        self.calc *= scale
        self.scale = scale

    def monte_carlo(self) -> Data:
        """Generate a data set to run Monte-Carlo simulation"""
        data = deepcopy(self)
        data.exp = np.random.normal(self.calc, self.err)
        return data

    def bootstrap(self: Data) -> Data:
        """Generate a data set to run Bootstrap simulation"""
        indexes = np.arange(self.metadata.size)
        pool1 = indexes[self.refs & self.mask]
        pool2 = indexes[~self.refs & self.mask]
        bs_indexes: list[int] = sorted(
            [*choices(pool1, k=pool1.size), *choices(pool2, k=pool2.size)]
        )
        data = deepcopy(self)
        data.metadata[self.mask] = self.metadata[bs_indexes]
        data.exp[self.mask] = self.exp[bs_indexes]
        data.err[self.mask] = self.err[bs_indexes]
        return data

    def sort(self) -> None:
        sorted_indexes = self.metadata.argsort()
        self.metadata = self.metadata[sorted_indexes]
        self.exp = self.exp[sorted_indexes]
        self.err = self.err[sorted_indexes]
        self.calc = self.calc[sorted_indexes]
        self.refs = self.refs[sorted_indexes]
        self.mask = self.mask[sorted_indexes]

    def any_duplicate(self):
        return np.unique(self.metadata).size != self.metadata.size

    def __add__(self: Data, other: Data) -> Data:
        # Make sure that both datasets have the same number of points
        assert self.size == other.size
        result = deepcopy(self)
        result.exp = self.exp + other.exp
        result.err = np.sqrt(self.err**2 + other.err**2)
        return result
