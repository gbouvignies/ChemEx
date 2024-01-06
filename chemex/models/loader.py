"""A simple plugin loader."""
from __future__ import annotations

import importlib
from pkgutil import iter_modules
from typing import Protocol, cast

from chemex.models import kinetic


class SettingModuleInterface(Protocol):
    """Represents a setting module interface.

    A setting module plugin has a single register function.
    """

    @staticmethod
    def register() -> None:
        """Register the necessary items related to the settings."""
        ...


def import_module(name: str) -> SettingModuleInterface:
    """Imports a module given a name."""
    return cast(SettingModuleInterface, importlib.import_module(name))


def register_kinetic_settings() -> None:
    """Loads the plugins defined in the plugins list."""
    for module in iter_modules(kinetic.__path__):
        module_name = f"{kinetic.__name__}.{module.name}"
        setting_module = import_module(module_name)
        setting_module.register()
