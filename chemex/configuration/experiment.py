from __future__ import annotations

from dataclasses import dataclass, field
from typing import Generic, TypeVar

import numpy as np
from pydantic import BaseModel, ConfigDict

from chemex.configuration.base import BaseModelLowerCase
from chemex.configuration.conditions import ConditionsFromFile


@dataclass
class ToBeFitted:
    rates: list[str] = field(default_factory=list)
    model_free: list[str] = field(default_factory=list)


class ExperimentNameSettings(BaseModelLowerCase):
    model_config = ConfigDict(str_to_lower=True)
    name: str


class RelaxationSettings(BaseModelLowerCase):
    model_config = ConfigDict(str_to_lower=True)
    name: str


class CpmgSettings(BaseModelLowerCase):
    model_config = ConfigDict(str_to_lower=True)
    name: str
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(CpmgSettings):
    even_ncycs: bool = True


class CestSettings(BaseModelLowerCase):
    model_config = ConfigDict(str_to_lower=True)
    name: str
    sw: float = np.inf


class ShiftSettings(BaseModelLowerCase):
    model_config = ConfigDict(str_to_lower=True)
    name: str


class ExperimentNameConfig(BaseModel):
    experiment: ExperimentNameSettings


ExperimentSettingsType = TypeVar("ExperimentSettingsType", bound=BaseModelLowerCase)
DataSettingsType = TypeVar("DataSettingsType", bound=BaseModelLowerCase)


class ExperimentConfig(BaseModel, Generic[ExperimentSettingsType, DataSettingsType]):
    experiment: ExperimentSettingsType
    conditions: ConditionsFromFile
    data: DataSettingsType

    @property
    def to_be_fitted(self) -> ToBeFitted:
        return ToBeFitted()
