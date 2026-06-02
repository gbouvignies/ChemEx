from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal, Self

from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    PlainValidator,
    ValidationError,
    field_validator,
    model_validator,
)
from pydantic.types import NonNegativeInt, PositiveInt

from chemex.configuration.utils import key_to_lower
from chemex.messages import print_method_error
from chemex.parameters.spin_system import SpinSystem
from chemex.toml import read_toml

# Type definitions
AllType = Annotated[Literal["*", "all"], BeforeValidator(str.lower)]
CoercedSpinSystem = Annotated[SpinSystem, PlainValidator(SpinSystem.from_name)]
SelectionType = list[SpinSystem] | Literal["*"] | None
McmcBurnSetting = NonNegativeInt | Literal["auto"]


class McmcSettings(BaseModel):
    steps: PositiveInt
    burn: McmcBurnSetting = "auto"
    thin: PositiveInt = 1
    walkers: PositiveInt | None = None
    seed: int | None = None
    workers: PositiveInt = 1
    update_parameters: bool = False

    _key_to_lower = model_validator(mode="before")(key_to_lower)

    @field_validator("burn", mode="before")
    @classmethod
    def parse_burn(cls, value: object) -> object:
        if isinstance(value, str) and value.lower() == "auto":
            return "auto"
        return value

    @model_validator(mode="after")
    def validate_sample_window(self) -> Self:
        burn = 0 if self.burn == "auto" else self.burn
        if burn >= self.steps:
            msg = "MCMC burn must be smaller than steps"
            raise ValueError(msg)
        retained = (self.steps - burn) // self.thin
        if retained < 1:
            msg = "MCMC settings must retain at least one sample"
            raise ValueError(msg)
        return self


class Statistics(BaseModel):
    mc: PositiveInt | None = None
    bs: PositiveInt | None = None
    bsn: PositiveInt | None = None
    mcmc: McmcSettings | None = None

    _key_to_lower = model_validator(mode="before")(key_to_lower)

    @field_validator("mcmc", mode="before")
    @classmethod
    def parse_mcmc_settings(cls, value: int | dict | None) -> dict | None:
        if isinstance(value, int):
            return {"steps": value}
        return value


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

    _key_to_lower = model_validator(mode="before")(key_to_lower)

    @field_validator("include", "exclude", mode="before")
    @classmethod
    def parse_residue_list(
        cls, value: list[str | int] | str | None
    ) -> list[SpinSystem] | Literal["*"] | None:
        if isinstance(value, list):
            for residue in value:
                if isinstance(residue, str) and residue.lower() in ("*", "all"):
                    return "*"
            return [SpinSystem.from_name(residue) for residue in value]
        if isinstance(value, str) and value.lower() in ("*", "all"):
            return "*"
        msg = f"Invalid residue list: {value}"
        raise ValueError(msg)

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
