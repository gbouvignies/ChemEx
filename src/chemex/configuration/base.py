"""Configuration models for experiment settings.

Imports dataclasses, typing generics, pydantic BaseModel for validation,
and utility functions.

Defines:

- ToBeFitted: Dataclass for tracking parameters to fit.
- BaseSettings: BaseModel with validator to lowercase dict keys.
- ExperimentConfiguration: Main config class composing experiment, conditions, data.
"""
from dataclasses import dataclass, field
from typing import Generic, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from chemex.configuration.utils import key_to_lower

T = TypeVar("T")
ExperimentSettings = TypeVar("ExperimentSettings", bound=BaseModel)
ConditionsSettings = TypeVar("ConditionsSettings", bound=BaseModel)
DataSettings = TypeVar("DataSettings", bound=BaseModel)


@dataclass
class ToBeFitted:
    rates: list[str] = field(default_factory=list)
    model_free: list[str] = field(default_factory=list)


class BaseSettings(BaseModel):
    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class ExperimentConfiguration(
    BaseModel,
    Generic[ExperimentSettings, ConditionsSettings, DataSettings],
):
    experiment: ExperimentSettings
    conditions: ConditionsSettings
    data: DataSettings

    @property
    def to_be_fitted(self) -> ToBeFitted:
        return ToBeFitted()
