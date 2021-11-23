from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from pydantic import BaseModel
from pydantic import conlist
from pydantic import validator

from chemex.parameters.name import ParamName
from chemex.toml import read_toml_multi

# Type definitions
ValueListType = conlist(float, min_items=1, max_items=4)
DefaultType = tuple[ParamName, "DefaultSetting"]
DefaultListType = list[DefaultType]


@dataclass(frozen=True)
class DefaultSetting:
    value: float
    min: float | None = None
    max: float | None = None
    brute_step: float | None = None


class ParamsConfig(BaseModel):
    __root__: dict[str, dict[str, ValueListType]]

    @validator("__root__", pre=True)
    def to_lower(cls, values):
        return {
            k1.lower(): {k2.lower(): v2 for k2, v2 in v1.items()}
            for k1, v1 in values.items()
        }

    @validator("__root__", pre=True)
    def to_list(cls, values):
        for values1 in values.values():
            for key2, values2 in values1.items():
                if not isinstance(values2, Iterable):
                    values1[key2] = [values2]
        return values

    @validator("__root__")
    def reorder(cls, values):
        return {"global": values.pop("global", {}), **values}

    def to_defaults_list(self) -> DefaultListType:
        defaults_list: DefaultListType = []
        for section, settings in self.__root__.items():
            prefix = f"{section}, NUC->" if section != "global" else ""
            for key, values in settings.items():
                pname = ParamName.from_section(f"{prefix}{key}")
                default_values = DefaultSetting(*values)
                defaults_list.append((pname, default_values))
        return defaults_list


def read_defaults(filenames: Iterable[Path]) -> DefaultListType:
    config = read_toml_multi(filenames)
    return ParamsConfig.parse_obj(config).to_defaults_list()
