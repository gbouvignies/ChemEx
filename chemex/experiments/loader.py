"""A simple plugin loader."""
from __future__ import annotations

import importlib
from pkgutil import iter_modules
from typing import Protocol, cast

from chemex.experiments import catalog


class ExperimentModuleInterface(Protocol):
    """Represents an experiment module interface.
    An experiment module plugin has a single register function.
    """

    @staticmethod
    def register() -> None:
        """Register the necessary items related to the experiment."""
        ...


def import_module(name: str) -> ExperimentModuleInterface:
    """Imports a module given a name."""
    return cast(ExperimentModuleInterface, importlib.import_module(name))


def register_experiments() -> None:
    """Loads the plugins defined in the plugins list."""
    for module in iter_modules(catalog.__path__):
        module_name = f"{catalog.__name__}.{module.name}"
        experiment_module = import_module(module_name)
        experiment_module.register()
