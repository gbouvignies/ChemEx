from __future__ import annotations

from collections.abc import Hashable
from functools import total_ordering
from typing import Annotated, Literal, Self, TypeVar

from pydantic import (
    BaseModel,
    BeforeValidator,
    Field,
    NonNegativeFloat,
    ValidationError,
    field_validator,
    model_validator,
)

from chemex.configuration.utils import key_to_lower, to_lower

T = TypeVar("T")
LabelType = Annotated[Literal["1h", "2h", "13c", "15n"], BeforeValidator(to_lower)]


@total_ordering
class Conditions(BaseModel, frozen=True):
    """Represents experimental conditions for NMR measurements.

    Attributes:
        h_larmor_frq (NonNegativeFloat | None): Larmor frequency of Hydrogen in MHz.
        temperature (float | None): Experimental temperature in degrees Celsius.
        p_total (NonNegativeFloat | None): Total concentration of protein in the sample.
        l_total (NonNegativeFloat | None): Total ligand concentration in the sample.
        d2o (float | None): Fraction of D2O in the solvent, between 0 and 1.
        label (tuple[LabelType, ...]): Tuple of NMR active isotopes used in the
                                       experiment.
    """

    h_larmor_frq: NonNegativeFloat | None = None
    temperature: float | None = None
    p_total: NonNegativeFloat | None = None
    l_total: NonNegativeFloat | None = None
    d2o: float | None = Field(gt=0.0, lt=1.0, default=None)
    label: tuple[LabelType, ...] = ()

    def rounded(self) -> Self:
        """Returns a new instance with rounded h_larmor_frq and temperature."""
        h_larmor_frq = round(self.h_larmor_frq, 1) if self.h_larmor_frq else None
        temperature = round(self.temperature, 1) if self.temperature else None
        return self.model_copy(
            update={"h_larmor_frq": h_larmor_frq, "temperature": temperature},
        )

    def match(self, other: Self) -> bool:
        """Checks if the current instance is equivalent to another."""
        return self == self & other

    @property
    def search_keys(self) -> set[Hashable]:
        """Creates a set of hashable search keys based on conditions."""
        return {
            self.h_larmor_frq,
            self.temperature,
            self.p_total,
            self.l_total,
            self.d2o,
        }

    def _conditions_list(self) -> list[tuple[str, str]]:
        """Helper method to create a list of condition keys and values."""
        conditions: list[tuple[str, str]] = []
        if self.temperature is not None:
            conditions.append(("T", f"{self.temperature:.1f}C"))
        if self.h_larmor_frq is not None:
            conditions.append(("B0", f"{self.h_larmor_frq:.1f}MHz"))
        if self.p_total is not None:
            conditions.append(("[P]", f"{self.p_total:e}M"))
        if self.l_total is not None:
            conditions.append(("[L]", f"{self.l_total:e}M"))
        if self.d2o is not None:
            conditions.append(("D2O", f"{self.d2o:.4f}"))

        return conditions

    @property
    def section(self) -> str:
        """Generates a string representation of the conditions."""
        return ", ".join(f"{key}->{value}" for key, value in self._conditions_list())

    @property
    def folder(self) -> str:
        """Generates a folder name representation of the conditions."""
        return "_".join(value for _, value in self._conditions_list())

    @property
    def is_deuterated(self) -> bool:
        """Checks if the conditions include deuterium."""
        return "2h" in self.label

    def select_conditions(self, conditions_selection: tuple[str, ...]) -> Self:
        """Selects specific conditions based on given keys.

        Args:
            conditions_selection (tuple[str, ...]): Keys to select conditions by.

        Returns:
            A new instance of Conditions with selected conditions.
        """
        return type(self).model_construct(
            **{
                key: value
                for key, value in self.model_dump().items()
                if key in conditions_selection
            },
        )

    def __and__(self, other: object) -> Self:
        """Defines bitwise AND operation for Conditions instances."""
        if not isinstance(other, type(self)):
            return NotImplemented
        self_dict = self.model_dump()
        other_dict = other.model_dump()
        intersection = {
            key: value for key, value in self_dict.items() if other_dict[key] == value
        }
        return type(self).model_construct(**intersection)

    def __lt__(self, other: object) -> bool:
        """Defines less than operation for Conditions instances."""
        if not isinstance(other, type(self)):
            return NotImplemented
        tuple_self = tuple(
            value if value is not None else -1e16
            for value in self.model_dump().values()
        )
        tuple_other = tuple(
            value if value is not None else -1e16
            for value in other.model_dump().values()
        )
        return tuple_self < tuple_other


class ConditionsWithValidations(Conditions, frozen=True):
    _key_to_lower = model_validator(mode="before")(key_to_lower)

    @field_validator("d2o")
    @classmethod
    def validate_d2o(cls, d2o: float | None) -> float | None:
        """Validates the d2o field for specific model requirements."""
        from chemex.models.model import model

        if "hd" in model.name and d2o is None:
            msg = 'To use the "hd" model, d2o must be provided'
            raise ValidationError(msg)
        return d2o

    @field_validator("temperature")
    @classmethod
    def validate_temperature(cls, temperature: float | None) -> float | None:
        """Validates temperature field for specific model requirements."""
        from chemex.models.model import model

        if "eyring" in model.name and temperature is None:
            msg = 'To use the "eyring" model, "temperature" must be provided'
            raise ValidationError(msg)
        return temperature

    @model_validator(mode="after")
    def validate_p_total_l_total(self) -> Self:
        """Validates both p_total and l_total for specific model requirements."""
        from chemex.models.model import model

        are_not_both_set = self.p_total is None or self.l_total is None
        if "binding" in model.name and are_not_both_set:
            msg = 'To use the "binding" model, "p_total" and "l_total" must be provided'
            raise ValueError(msg)
        return self
