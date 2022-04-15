from __future__ import annotations

from collections.abc import MutableMapping
from typing import Any

from pydantic import BaseModel
from pydantic import root_validator


class BaseModelLowerCase(BaseModel):
    class Config:
        anystr_lower = True

    @root_validator(pre=True)
    def to_lower_case(
        cls, values: MutableMapping[str, Any]
    ) -> MutableMapping[str, Any]:
        return {k.lower(): v for k, v in values.items()}
