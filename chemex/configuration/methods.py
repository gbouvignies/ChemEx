from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING, Annotated, Literal

from pydantic import BeforeValidator, ConfigDict, Field, ValidationError
from pydantic.types import PositiveInt

from chemex.configuration.base import BaseModelLowerCase
from chemex.messages import print_method_error
from chemex.parameters.spin_system import PydanticSpinSystem
from chemex.toml import read_toml

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


# Type definitions
AllType = Annotated[Literal["*", "all"], BeforeValidator(str.lower)]
SelectionType = list[PydanticSpinSystem] | AllType | None


class Statistics(BaseModelLowerCase):
    mc: PositiveInt | None = None
    bs: PositiveInt | None = None
    bsn: PositiveInt | None = None


@dataclass
class Selection:
    include: SelectionType
    exclude: SelectionType


class Method(BaseModelLowerCase):
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
