from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pydantic import Field, ValidationError
from pydantic.types import PositiveInt

from chemex.configuration.base import BaseModelLowerCase
from chemex.messages import print_method_error
from chemex.parameters.name import SpinSystem
from chemex.toml import read_toml

# Type definitions
AllType = Literal["*", "all", "ALL", "All", "ALl", "AlL"]
SelectionType = list[SpinSystem] | AllType | None


class Statistics(BaseModelLowerCase):
    mc: PositiveInt | None = None
    bs: PositiveInt | None = None
    bsn: PositiveInt | None = None


@dataclass
class Selection:
    include: SelectionType
    exclude: SelectionType


class Method(BaseModelLowerCase):
    class Config:
        anystr_lower = True

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
    methods = {}

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
