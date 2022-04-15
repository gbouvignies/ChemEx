from __future__ import annotations

from collections.abc import Hashable
from collections.abc import MutableMapping
from functools import total_ordering
from typing import Any
from typing import Literal
from typing import Optional

from pydantic import BaseModel
from pydantic import Field
from pydantic import MissingError
from pydantic import PositiveFloat
from pydantic import root_validator
from pydantic import validator

from chemex.model import model


@total_ordering
class Conditions(BaseModel, frozen=True):
    h_larmor_frq: Optional[PositiveFloat]
    temperature: Optional[float] = None
    p_total: Optional[PositiveFloat] = Field(default=None)
    l_total: Optional[PositiveFloat] = Field(default=None)
    d2o: Optional[float] = Field(gt=0.0, lt=1.0, default=None)
    label: tuple[Literal["1h", "2h", "13c", "15n"], ...] = tuple()

    def rounded(self) -> Conditions:
        h_larmor_frq = round(self.h_larmor_frq, 1) if self.h_larmor_frq else None
        temperature = round(self.temperature, 1) if self.temperature else None
        return self.copy(
            update={"h_larmor_frq": h_larmor_frq, "temperature": temperature}
        )

    def match(self, other: Conditions) -> bool:
        return self == self & other

    @property
    def search_keys(self) -> set[Hashable]:
        return {
            self.h_larmor_frq,
            self.temperature,
            self.p_total,
            self.l_total,
            self.d2o,
        }

    @property
    def section(self) -> str:
        parts = []
        if self.temperature is not None:
            parts.append(f"T->{self.temperature:.1f}C")
        if self.h_larmor_frq is not None:
            parts.append(f"B0->{self.h_larmor_frq:.1f}MHz")
        if self.p_total is not None:
            parts.append(f"[P]->{self.p_total:e}M")
        if self.l_total is not None:
            parts.append(f"[L]->{self.l_total:e}M")
        if self.d2o is not None:
            parts.append(f"D2O->{self.d2o:.4f}")
        return ", ".join(parts)

    @property
    def folder(self):
        parts = []
        if self.temperature is not None:
            parts.append(f"{self.temperature:.1f}C")
        if self.h_larmor_frq is not None:
            parts.append(f"{self.h_larmor_frq:.1f}MHz")
        if self.p_total is not None:
            parts.append(f"P{self.p_total:e}M")
        if self.l_total is not None:
            parts.append(f"L{self.l_total:e}M")
        if self.d2o is not None:
            parts.append(f"D{self.d2o:.4f}")
        return "_".join(parts)

    @property
    def is_deuterated(self) -> bool:
        return "2h" in self.label

    def select_conditions(self, conditions_selection: tuple[str, ...]) -> Conditions:
        return Conditions.construct(
            **{
                key: value
                for key, value in self.dict().items()
                if key in conditions_selection
            }
        )

    def __and__(self, other: Conditions) -> Conditions:
        self_dict = self.dict()
        other_dict = other.dict()
        intersection = {
            key: value for key, value in self_dict.items() if other_dict[key] == value
        }
        return Conditions.construct(**intersection)

    def __lt__(self, other: Conditions) -> bool:
        tuple_self = tuple(
            value if value is not None else -1e16 for value in self.dict().items()
        )
        tuple_other = tuple(
            value if value is not None else -1e16 for value in other.dict().items()
        )
        return tuple_self < tuple_other


@total_ordering
class ConditionsFromFile(Conditions, frozen=True):
    @validator("d2o")
    def validate_d2o(cls, d2o):
        if "hd" in model.name and d2o is None:
            raise MissingError
        return d2o

    @validator("temperature")
    def validate_temperature(cls, temperature):
        if "eyring" in model.name and temperature is None:
            raise MissingError
        return temperature

    @validator("label", pre=True)
    def set_to_lower_case(cls, label):
        return tuple(value.lower() for value in label)

    @root_validator()
    def validate_p_total_l_total(cls, values):
        are_not_both_set = values["p_total"] is None or values["l_total"] is None
        if "binding" in model.name and are_not_both_set:
            raise ValueError("Either p_total or l_total must be provided")
        return values

    @root_validator(pre=True)
    def set_keys_to_lower_case(
        cls, values: MutableMapping[str, Any]
    ) -> MutableMapping[str, Any]:
        return {k.lower(): v for k, v in values.items()}
