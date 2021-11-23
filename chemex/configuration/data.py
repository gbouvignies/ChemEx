from __future__ import annotations

from pathlib import Path
from typing import Literal
from typing import Union

from pydantic import validator

from chemex.configuration.base import BaseModelLowerCase
from chemex.parameters.spin_system import SpinSystem


class DataSettings(BaseModelLowerCase):
    path: Path = Path("./")
    scaled: bool = True


class RelaxationDataSettings(DataSettings):
    error: Literal["file", "duplicates"] = "file"
    filter_planes: list[int] = []
    profiles: dict[SpinSystem, Union[Path, list[Path]]] = {}

    @validator("profiles", pre=True)
    def make_list(cls, v):
        if isinstance(v, dict):
            for key, value in v.items():
                if isinstance(value, str):
                    v[key] = [value]
        return v


class CestDataSettings(DataSettings):
    error: Literal["file", "duplicates", "scatter"] = "file"
    filter_planes: list[int] = []
    filter_offsets: list[tuple[float, float]] = [(0.0, 0.0)]
    filter_ref_planes: bool = False
    profiles: dict[SpinSystem, list[Path]] = {}

    @validator("profiles", pre=True)
    def make_list(cls, v):
        if isinstance(v, dict):
            for key, value in v.items():
                if isinstance(value, str):
                    v[key] = [value]
        return v


class CestDataSettingsNoRef(CestDataSettings):
    filter_ref_planes: bool = True


class ShiftDataSettings(DataSettings):
    error: Literal["file", "duplicates"] = "file"
    scaled: bool = False
