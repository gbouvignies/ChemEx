from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

from annotated_types import Len
from pydantic import AfterValidator, BeforeValidator, RootModel, ValidationError

from chemex.configuration.utils import ensure_list
from chemex.exceptions import ChemExError
from chemex.parameters.name import ParamName
from chemex.toml import read_toml_multi


def rename_section(section_name: str) -> str:
    if section_name == "global":
        return ""
    return f"{section_name},nuc->"


ValuesType = Annotated[
    list[float],
    Len(min_length=1, max_length=4),
    BeforeValidator(ensure_list),
]
LowerCaseString = Annotated[str, BeforeValidator(str.lower)]
ValuesDictType = dict[LowerCaseString, ValuesType]
SectionType = Annotated[LowerCaseString, AfterValidator(rename_section)]
ParamsConfigType = dict[SectionType, ValuesDictType]
ParamsConfigModel = RootModel[ParamsConfigType]


class ParameterConfigurationError(ChemExError, ValueError):
    """User parameter defaults do not satisfy the parameter-file schema."""

    def __init__(self, filenames: tuple[Path, ...], explanation: str) -> None:
        super().__init__(explanation)
        self.filenames = filenames
        self.explanation = explanation


@dataclass(frozen=True)
class DefaultSetting:
    value: float
    min: float | None = None
    max: float | None = None
    brute_step: float | None = None


DefaultType = tuple[ParamName, DefaultSetting]
DefaultListType = list[DefaultType]


def build_default_list(params_config: ParamsConfigModel) -> DefaultListType:
    defaults: DefaultListType = []
    for section, params in params_config.root.items():
        for key, values in params.items():
            pname = ParamName.from_section(f"{section}{key}")
            default_values = DefaultSetting(*values)
            defaults.append((pname, default_values))
    return defaults


def read_defaults(filenames: Iterable[Path]) -> DefaultListType:
    sources = tuple(filenames)
    config = read_toml_multi(sources)
    try:
        param_config = ParamsConfigModel.model_validate(config)
    except ValidationError as error:
        first = error.errors()[0]
        location = " -> ".join(str(item) for item in first["loc"])
        explanation = str(first["msg"])
        if location:
            explanation = f"{location}: {explanation}"
        raise ParameterConfigurationError(sources, explanation) from error
    return build_default_list(param_config)
