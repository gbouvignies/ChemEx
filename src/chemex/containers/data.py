from copy import deepcopy
from random import choices
from typing import Any, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, computed_field, field_serializer

from chemex.typing import Array

rng = np.random.default_rng()


class Data(BaseModel):
    """Dataset with experimental, calculated data and metadata.

    Attributes:
        exp (Array): Experimental data array.
        err (Array): Error array for experimental data.
        metadata (Array): Metadata for data points.
        calc (Array): Calculated data, set after instantiation.
        calc_unscaled (Array): Unscaled calculated data, set after instantiation.
        mask (Array): Mask array for data selection, set after instantiation.
        refs (Array): Array of reference points, set after instantiation.

    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True,  # Allow numpy arrays
        validate_assignment=True,  # Validate on assignment
    )

    exp: Array
    err: Array
    metadata: Array = Field(default_factory=lambda: np.array([]))
    calc: Array = Field(init=False, default=None)  # type: ignore[call-arg]
    calc_unscaled: Array = Field(init=False, default=None)  # type: ignore[call-arg]
    mask: Array = Field(init=False, default=None)  # type: ignore[call-arg]
    refs: Array = Field(init=False, default=None)  # type: ignore[call-arg]

    def __init__(self, **data: Array) -> None:
        """Initialize the Data instance and set computed arrays."""
        super().__init__(**data)
        self.calc = np.full_like(self.exp, fill_value=1e32, dtype=np.float64)
        self.calc_unscaled = np.full_like(self.exp, fill_value=1e32, dtype=np.float64)
        self.mask = np.full_like(self.exp, fill_value=True, dtype=np.bool_)
        self.refs = np.full_like(self.exp, fill_value=False, dtype=np.bool_)

    # Serialize numpy arrays to lists for JSON output
    @field_serializer("exp", "err", "metadata", "calc", "calc_unscaled", "mask", "refs")
    def serialize_array(self, value: Array, _info: Any) -> list:  # noqa: ANN401
        """Convert numpy arrays to lists for JSON serialization."""
        if value is None:
            return []
        return value.tolist()

    @computed_field  # type: ignore[misc]
    @property
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
