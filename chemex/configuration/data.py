from __future__ import annotations

from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BeforeValidator, field_validator

from chemex.configuration.base import BaseModelLowerCase, ensure_list
from chemex.parameters.spin_system import PydanticSpinSystem


class DataSettings(BaseModelLowerCase):
    path: Path = Path("./")
    scaled: bool = True


PathList = Annotated[list[Path], BeforeValidator(ensure_list)]


class RelaxationDataSettings(DataSettings):
    error: Literal["file", "duplicates"] = "file"
    filter_planes: list[int] = []
    profiles: dict[PydanticSpinSystem, PathList] = {}

    @field_validator("error", mode="before")
    def to_lower(cls, error: Any) -> str:
        if isinstance(error, str):
            return error.lower()
        return error


class CestDataSettings(DataSettings):
    error: Literal["file", "duplicates", "scatter"] = "file"
    filter_planes: list[int] = []
    filter_offsets: list[tuple[float, float]] = [(0.0, 0.0)]
    filter_ref_planes: bool = False
    profiles: dict[PydanticSpinSystem, PathList] = {}

    @field_validator("error", mode="before")
    def to_lower(cls, error: Any) -> str:
        if isinstance(error, str):
            return error.lower()
        return error


class CestDataSettingsNoRef(CestDataSettings):
    filter_ref_planes: bool = True


class ShiftDataSettings(DataSettings):
    error: Literal["file", "duplicates"] = "file"
    scaled: bool = False

    @field_validator("error", mode="before")
    def to_lower(cls, error: Any) -> str:
        if isinstance(error, str):
            return error.lower()
        return error
