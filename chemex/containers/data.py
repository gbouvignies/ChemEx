from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from random import choices
from typing import Any, Self

import numpy as np
from numpy.typing import NDArray

from chemex.typing import ArrayBool, ArrayFloat

rng = np.random.default_rng()


@dataclass
class Data:
    """Dataset with experimental, calculated data and metadata.

    Attributes:
        exp (ArrayFloat): Experimental data array.
        err (ArrayFloat): Error array for experimental data.
        metadata (NDArray[Any]): Metadata for data points.
        calc (ArrayFloat): Calculated data, set after instantiation.
        calc (ArrayFloat): Unscaled calculated data, set after instantiation.
        mask (ArrayBool): Mask array for data selection, set after instantiation.
        refs (ArrayBool): Array of reference points, set after instantiation.
    """

    exp: ArrayFloat
    err: ArrayFloat
    metadata: NDArray[Any] = field(default_factory=lambda: np.array([]))
    calc: ArrayFloat = field(init=False)
    calc_unscaled: ArrayFloat = field(init=False)
    mask: ArrayBool = field(init=False)
    refs: ArrayBool = field(init=False)

    def __post_init__(self) -> None:
        """Initialize computed attributes of the Data class."""
        self.calc = np.full_like(self.exp, fill_value=1e32, dtype=np.float64)
        self.calc_unscaled = np.full_like(self.exp, fill_value=1e32, dtype=np.float64)
        self.mask = np.full_like(self.exp, fill_value=True, dtype=np.bool_)
        self.refs = np.full_like(self.exp, fill_value=False, dtype=np.bool_)

    @cached_property
    def size(self) -> int:
        """Return the number of data points."""
        return self.exp.size

    @property
    def scale(self) -> float:
        """Calculates the scale factor between experimental and calculated data.

        This method computes the scale factor to align calculated data with experimental
        data, considering associated errors. Includes safeguards for numerical
        stability.

        Returns:
            float: Scale factor for aligning calculated and experimental data.
        """
        expe = self.exp[self.mask]
        calc = self.calc_unscaled[self.mask]
        error = self.err[self.mask]

        # Adding a small constant to avoid division by zero or very small numbers
        error_safe = error + np.finfo(float).eps

        expe_norm = expe / error_safe
        calc_norm = calc / error_safe

        dot_product_num = np.dot(expe_norm, calc_norm)
        dot_product_den = np.dot(calc_norm, calc_norm)

        return dot_product_num / dot_product_den

    def monte_carlo(self) -> Self:
        """Generate Data instance for Monte-Carlo simulation.

        Returns:
            Self: New Data instance for Monte-Carlo simulation.
        """
        data = deepcopy(self)
        data.exp = rng.normal(self.calc, self.err)
        return data

    def bootstrap(self) -> Self:
        """Generate Data instance using Bootstrap resampling.

        Returns:
            Self: New Data instance for Bootstrap simulation.
        """
        indexes = np.arange(self.metadata.size)
        pool1 = indexes[self.refs & self.mask]
        pool2 = indexes[~self.refs & self.mask]
        bs_indexes: list[int] = sorted(
            [*choices(pool1, k=pool1.size), *choices(pool2, k=pool2.size)],
        )
        data = deepcopy(self)
        data.metadata[self.mask] = self.metadata[bs_indexes]
        data.exp[self.mask] = self.exp[bs_indexes]
        data.err[self.mask] = self.err[bs_indexes]
        return data

    def sort(self) -> None:
        """Sort the data based on metadata values."""
        sorted_indexes = self.metadata.argsort()
        self.metadata = self.metadata[sorted_indexes]
        self.exp = self.exp[sorted_indexes]
        self.err = self.err[sorted_indexes]
        self.calc = self.calc[sorted_indexes]
        self.refs = self.refs[sorted_indexes]
        self.mask = self.mask[sorted_indexes]

    def any_duplicate(self) -> bool:
        """Check for duplicate entries in metadata.

        Returns:
            bool: True if duplicates exist, False otherwise.
        """
        return np.unique(self.metadata).size != self.metadata.size

    def __add__(self, other: Self) -> Self:
        """Add two Data instances, combining their experimental data and errors.

        Both Data instances must have the same number of points.

        Args:
            other (Self): Another Data instance.

        Returns:
            Self: New Data instance representing the sum of the two datasets.

        Raises:
            ValueError: If the sizes of the datasets do not match.
        """
        # Check if the sizes of the datasets match
        if self.size != other.size:
            msg = "Sizes of the datasets must be identical to add them."
            raise ValueError(msg)

        # Perform addition if sizes match
        result = deepcopy(self)
        result.exp = self.exp + other.exp
        result.err = np.sqrt(self.err**2 + other.err**2)

        return result
