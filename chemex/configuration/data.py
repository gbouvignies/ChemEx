from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, TypeVar

from pydantic import BaseModel, BeforeValidator, Field, PlainValidator, model_validator

from chemex.configuration.utils import ensure_list, to_lower
from chemex.parameters.spin_system import SpinSystem

T = TypeVar("T")
ErrorType = Annotated[
    Literal["file", "duplicates", "scatter"],
    BeforeValidator(to_lower),
]
PathList = Annotated[list[Path], BeforeValidator(ensure_list)]
ProfilesType = dict[
    Annotated[SpinSystem, PlainValidator(lambda x: SpinSystem(x))],
    PathList,
]


class DataSettings(BaseModel):
    error: ErrorType = "file"
    global_error: bool = True
    path: Path = Path("./")
    scaled: bool = True

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, T]) -> dict[str, T]:
        """Model validator to convert all dictionary keys to lowercase."""
        return {to_lower(k): v for k, v in model.items()}


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
