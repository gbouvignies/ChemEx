from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, TypeVar

from pydantic import BaseModel, BeforeValidator, Field, PlainValidator, model_validator

from chemex.configuration.utils import ensure_list, key_to_lower, to_lower
from chemex.parameters.spin_system import SpinSystem

T = TypeVar("T")
ErrorType = Annotated[
    Literal["file", "duplicates", "scatter"],
    BeforeValidator(to_lower),
]
PathList = Annotated[list[Path], BeforeValidator(ensure_list)]
ProfilesType = dict[
    Annotated[SpinSystem, PlainValidator(SpinSystem.from_name)], PathList
]


class DataSettings(BaseModel):
    error: ErrorType = "file"
    global_error: bool = True
    path: Path = Path("./")
    scaled: bool = True

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class RelaxationDataSettings(DataSettings):
    filter_planes: list[int] = Field(default_factory=list)
    profiles: ProfilesType = Field(default_factory=dict)


class CestDataSettings(RelaxationDataSettings):
    filter_offsets: list[tuple[float, float]] = Field(default_factory=list)
    filter_ref_planes: bool = False


class CestDataSettingsNoRef(CestDataSettings):
    filter_ref_planes: bool = True


class ShiftDataSettings(DataSettings):
    scaled: bool = False
