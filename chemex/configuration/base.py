from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, model_validator


def ensure_list(variable: Any | list[Any] | None) -> list[Any]:
    if isinstance(variable, list):
        return variable
    if variable is None:
        return []
    return [variable]


def to_lower(string: Any) -> Any:
    if isinstance(string, str):
        return string.lower()
    return string


class BaseModelLowerCase(BaseModel):
    model_config = ConfigDict(str_to_lower=True)

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, Any]) -> dict[str, Any]:
        return {to_lower(k): v for k, v in model.items()}
