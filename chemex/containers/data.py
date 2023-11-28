from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from random import choices
from typing import TYPE_CHECKING, Any, TypeVar

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from chemex.typing import ArrayBool, ArrayFloat

Self = TypeVar("Self", bound="Data")

rng = np.random.default_rng()


def get_scale(exp: ArrayFloat, err: ArrayFloat, calc: ArrayFloat) -> float:
    """Calculate and return the scale factor for the calculated data.

    Args:
        exp (ArrayFloat): The experimental data array.
        err (ArrayFloat): The error associated with the experimental data.
        calc (ArrayFloat): The calculated data array.

    Returns:
        float: The calculated scale factor.

    Raises:
        ZeroDivisionError: If division by zero occurs in the scale calculation.
    """
    if (calc == 0).all():
        return 1.0
    try:
        calc_err2 = calc / err**2
        return sum(exp * calc_err2) / sum(calc * calc_err2)
    except ZeroDivisionError:
        return 1.0


@dataclass
class Data:
    """Dataset with experimental, calculated data and metadata.

    Attributes:
        exp (ArrayFloat): Experimental data.
        err (ArrayFloat): Error in experimental data.
        metadata (NDArray[Any]): Data metadata.
        scale (float): Scale factor for calculated data.
        size (int): Data points count, set post-instantiation.
        calc (ArrayFloat): Calculated data, set post-instantiation.
        mask (ArrayBool): Mask array, set post-instantiation.
        refs (ArrayBool): Reference points array, set post-instantiation.
    """

    exp: ArrayFloat
    err: ArrayFloat
    metadata: NDArray[Any] = field(default_factory=lambda: np.array([]))
    scale: float = 1.0
    size: int = field(init=False)
    calc: ArrayFloat = field(init=False)
    mask: ArrayBool = field(init=False)
    refs: ArrayBool = field(init=False)

    def __post_init__(self):
        """Initialize computed attributes of the Data class."""
        self.size = self.exp.size
        self.calc = np.full_like(self.exp, 1e32)
        self.mask = np.full_like(self.exp, True, dtype=np.bool_)
        self.refs = np.full_like(self.exp, False, dtype=np.bool_)

    def scale_calc(self) -> None:
        """Scale the calculated data based on the experimental data and mask."""
        mask = self.mask
        scale = get_scale(self.exp[mask], self.err[mask], self.calc[mask])
        self.calc *= scale
        self.scale = scale

    def prepare_for_simulation(self) -> None:
        """Prepare the dataset for simulation.

        Set experimental data to calculated values and errors to zero.
        """
        self.exp = self.calc.copy()
        self.err *= 0.0
        if any(self.refs):
            scale = float(100.0 / np.mean(self.exp[self.refs]))
            self.exp *= scale
            self.calc *= scale
            self.scale *= scale

    def monte_carlo(self: Self) -> Self:
        """Generate Data instance for Monte-Carlo simulation.

        Returns:
            Self: New Data instance for Monte-Carlo simulation.
        """
        data = deepcopy(self)
        data.exp = rng.normal(self.calc, self.err)
        return data

    def bootstrap(self: Self) -> Self:
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

    def any_duplicate(self):
        """Check for duplicate entries in metadata.

        Returns:
            bool: True if duplicates exist, False otherwise.
        """
        return np.unique(self.metadata).size != self.metadata.size

    def __add__(self: Self, other: Self) -> Self:
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
