from __future__ import annotations

import math
from typing import Literal, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from chemex.configuration.utils import key_to_lower

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"

    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class CpmgSettings(ExperimentSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(ExperimentSettings):
    even_ncycs: bool = True


class CestSettings(ExperimentSettings):
    sw: float = math.inf
