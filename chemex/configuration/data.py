from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import field_validator

from chemex.configuration.base import BaseModelLowerCase
from chemex.parameters.spin_system import PydanticSpinSystem


class DataSettings(BaseModelLowerCase):
    path: Path = Path("./")
    scaled: bool = True


class RelaxationDataSettings(DataSettings):
    error: Literal["file", "duplicates"] = "file"
    filter_planes: list[int] = []
    profiles: dict[PydanticSpinSystem, Path | list[Path]] = {}

    @field_validator("profiles", mode="before")
    @classmethod
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
    profiles: dict[PydanticSpinSystem, list[Path]] = {}

    @field_validator("profiles", mode="before")
    @classmethod
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
