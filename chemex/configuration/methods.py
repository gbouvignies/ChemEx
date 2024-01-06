from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal, TypeVar

from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    PlainValidator,
    ValidationError,
    model_validator,
)
from pydantic.types import PositiveInt

from chemex.configuration.utils import to_lower
from chemex.messages import print_method_error
from chemex.parameters.spin_system import SpinSystem
from chemex.toml import read_toml

# Type definitions
AllType = Annotated[Literal["*", "all"], BeforeValidator(str.lower)]
PydanticSpinSystem = Annotated[SpinSystem, PlainValidator(lambda x: SpinSystem(x))]
SelectionType = list[PydanticSpinSystem | AllType] | AllType | None
T = TypeVar("T")


class Statistics(BaseModel):
    mc: PositiveInt | None = None
    bs: PositiveInt | None = None
    bsn: PositiveInt | None = None

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, T]) -> dict[str, T]:
        """Model validator to convert all dictionary keys to lowercase."""
        return {to_lower(k): v for k, v in model.items()}


@dataclass
class Selection:
    include: SelectionType
    exclude: SelectionType


class Method(BaseModel):
    model_config = ConfigDict(str_to_lower=True, extra="forbid")
    fitmethod: str = "leastsq"
    include: SelectionType = None
    exclude: SelectionType = None
    fit: list[str] = Field(default_factory=list)
    fix: list[str] = Field(default_factory=list)
    constraints: list[str] = Field(default_factory=list)
    grid: list[str] = Field(default_factory=list)
    statistics: Statistics | None = None

    @property
    def selection(self) -> Selection:
        return Selection(include=self.include, exclude=self.exclude)

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, T]) -> dict[str, T]:
        """Model validator to convert all dictionary keys to lowercase."""
        return {to_lower(k): v for k, v in model.items()}


Methods = dict[str, Method]


def read_methods(filenames: Iterable[Path]) -> Methods:
    methods: Methods = {}

    for filename in filenames:
        methods_dict = read_toml(filename)
        for section, settings in methods_dict.items():
            try:
                method = Method(**settings)
            except ValidationError as error:
                options = {option for err in error.errors() for option in err["loc"]}
                print_method_error(filename, section, options)
                sys.exit()
            methods[section] = method
    return methods
