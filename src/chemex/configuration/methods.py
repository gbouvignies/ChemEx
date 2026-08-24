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

from chemex.configuration.method_plan import MethodFormatError, MethodPlan, SourceRef
from chemex.configuration.method_v1 import adapt_v1
from chemex.configuration.method_v2 import adapt_v2
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
    workers: PositiveInt | None = None
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
        if self.update_parameters:
            msg = "MCMC posterior sampling cannot mutate the committed central fit"
            raise ValueError(msg)
        if self.seed is not None and not 0 <= self.seed < 1 << 64:
            msg = "MCMC seed must be an unsigned 64-bit integer"
            raise ValueError(msg)
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
    fitmethod: Literal["trf"] = "trf"
    include: SelectionType = None
    exclude: SelectionType = None
    fit: list[str] = Field(default_factory=list)
    fix: list[str] = Field(default_factory=list)
    constraints: list[str] = Field(default_factory=list)
    grid: list[str] = Field(default_factory=list)
    statistics: Statistics | None = None

    _key_to_lower = model_validator(mode="before")(key_to_lower)

    @field_validator("fitmethod", mode="before")
    @classmethod
    def canonicalize_fitmethod(cls, value: object) -> object:
        if isinstance(value, str):
            normalized = value.lower()
            if normalized == "least_squares":
                return "trf"
            if normalized == "trf":
                return normalized
        msg = (
            "FITMETHOD supports only 'trf'; 'least_squares' is accepted as a "
            "temporary alias"
        )
        raise ValueError(msg)

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
type _MethodSource = tuple[Path, dict[str, object]]
type _RawStep = tuple[Path, str, dict[str, object]]


def _method_format(sources: list[_MethodSource]) -> bool:
    formats: set[int] = set()
    for filename, data in sources:
        version_items = [
            value for key, value in data.items() if str(key).lower() == "format_version"
        ]
        if len(version_items) > 1:
            raise MethodFormatError(
                "Duplicate FORMAT_VERSION declaration",
                SourceRef(filename, "<file>", "FORMAT_VERSION"),
            )
        if not version_items:
            formats.add(1)
        elif (
            isinstance(version_items[0], int)
            and not isinstance(version_items[0], bool)
            and version_items[0] == 2
        ):
            formats.add(2)
        else:
            raise MethodFormatError(
                "FORMAT_VERSION supports only version 2; omit it for legacy v1",
                SourceRef(filename, "<file>", "FORMAT_VERSION"),
            )
    if len(formats) > 1:
        filename = sources[-1][0]
        raise MethodFormatError(
            "A mixed v1/v2 method invocation is not supported",
            SourceRef(filename, "<file>", "FORMAT_VERSION"),
        )
    return formats == {2}


def _explicit_model_values(model: BaseModel) -> dict[str, object]:
    return {
        field: _explicit_model_values(value) if isinstance(value, BaseModel) else value
        for field in model.model_fields_set
        if (value := getattr(model, field)) is not None
    }


def _validate_v1_step(
    filename: Path, name: str, settings: dict[str, object]
) -> dict[str, object]:
    try:
        return _explicit_model_values(Method.model_validate(settings))
    except ValidationError as error:
        first = error.errors()[0]
        field = ".".join(str(item).upper() for item in first["loc"])
        raise MethodFormatError(
            str(first["msg"]),
            SourceRef(filename, name, field or "<step>"),
        ) from error


def _raw_steps(sources: list[_MethodSource], *, is_v2: bool) -> list[_RawStep]:
    raw_steps: list[_RawStep] = []
    step_positions: dict[str, int] = {}
    for filename, data in sources:
        for name, settings in data.items():
            if str(name).lower() == "format_version":
                continue
            if not isinstance(settings, dict):
                raise MethodFormatError(
                    "Method steps must be TOML tables",
                    SourceRef(filename, str(name), "<step>"),
                )
            normalized_settings = {str(key): value for key, value in settings.items()}
            if not is_v2:
                normalized_settings = _validate_v1_step(
                    filename, name, normalized_settings
                )
            if name in step_positions:
                if is_v2:
                    raise MethodFormatError(
                        f"Duplicate v2 step {name!r}",
                        SourceRef(filename, str(name), "<step>"),
                    )
                raw_steps[step_positions[name]] = (filename, name, normalized_settings)
            else:
                step_positions[name] = len(raw_steps)
                raw_steps.append((filename, name, normalized_settings))
    return raw_steps


def read_method_plan(filenames: Iterable[Path]) -> MethodPlan:
    sources: list[_MethodSource] = [
        (filename, read_toml(filename)) for filename in filenames
    ]
    is_v2 = _method_format(sources)
    raw_steps = _raw_steps(sources, is_v2=is_v2)
    if is_v2:
        return adapt_v2(raw_steps)
    return adapt_v1(raw_steps)


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
                sys.exit(1)
            methods[section] = method
    return methods
