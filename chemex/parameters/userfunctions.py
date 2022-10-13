"""Pace where user functions used by lmfit are stored."""
from __future__ import annotations

from typing import Any
from typing import Callable

from chemex.configuration.conditions import Conditions
from chemex.parameters.setting import ParamLocalSetting

SettingsType = dict[str, ParamLocalSetting]
SettingMakerType = Callable[[Conditions], SettingsType]


class Registry:
    """Registry for storing all the user functions for lmfit constraints."""

    user_function_registry: dict[str, Any] = {}

    def register(self, name: str, user_functions: dict[str, Callable]):
        """Register a new set of user functions."""
        self.user_function_registry[name] = user_functions

    def get(self, name: str) -> dict[str, Callable]:
        if name not in self.user_function_registry:
            return {}
        return self.user_function_registry[name]


user_function_registry = Registry()
