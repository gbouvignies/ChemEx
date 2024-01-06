from __future__ import annotations

import math
from typing import Literal, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from chemex.configuration.utils import to_lower

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"

    model_config = ConfigDict(str_to_lower=True)

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, T]) -> dict[str, T]:
        """Model validator to convert all dictionary keys to lowercase."""
        return {to_lower(k): v for k, v in model.items()}


class CpmgSettings(ExperimentSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(ExperimentSettings):
    even_ncycs: bool = True


class CestSettings(ExperimentSettings):
    sw: float = math.inf
