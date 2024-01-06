"""Pace where user functions used by lmfit are stored."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, ClassVar

from chemex.configuration.conditions import Conditions
from chemex.parameters.setting import ParamLocalSetting

SettingsType = dict[str, ParamLocalSetting]
SettingMakerType = Callable[[Conditions], SettingsType]


class Registry:
    """Registry for storing all the user functions for lmfit constraints."""

    user_function_registry: ClassVar[dict[str, Any]] = {}

    def register(
        self, name: str, user_functions: dict[str, Callable[[Any], Any]]
    ) -> None:
        """Register a new set of user functions."""
        self.user_function_registry[name] = user_functions

    def get(self, name: str) -> dict[str, Callable[[Any], Any]]:
        if name not in self.user_function_registry:
            return {}
        return self.user_function_registry[name]


user_function_registry = Registry()
