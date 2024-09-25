from __future__ import annotations

import math
from functools import cached_property
from typing import Literal, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from chemex.configuration.utils import key_to_lower

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"

    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class RelaxationSettings(ExperimentSettings):
    cs_evolution_prior: bool = False

    @cached_property
    def suffix(self) -> str:
        return f"_{self.observed_state}" if self.cs_evolution_prior else ""


class CpmgSettings(RelaxationSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(RelaxationSettings):
    even_ncycs: bool = True


class CestSettings(RelaxationSettings):
    sw: float = math.inf


class MFCestSettings(RelaxationSettings):
    sw: float
