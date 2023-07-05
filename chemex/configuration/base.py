from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pydantic import BaseModel, ConfigDict, model_validator

if TYPE_CHECKING:
    from collections.abc import MutableMapping


class BaseModelLowerCase(BaseModel):
    model_config = ConfigDict(str_to_lower=True)

    @model_validator(mode="before")
    @classmethod
    def to_lower_case(
        cls, values: MutableMapping[str, Any]
    ) -> MutableMapping[str, Any]:
        return {k.lower(): v for k, v in values.items()}
