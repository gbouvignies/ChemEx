"""A simple plugin loader for experiment modules."""

from __future__ import annotations

import importlib
from collections.abc import Iterator
from importlib.util import find_spec
from pkgutil import iter_modules
from types import ModuleType
from typing import Protocol

from chemex.experiments import catalog
from chemex.experiments.catalog import wip


class ExperimentModule(Protocol):
    """Represents an experiment module interface."""

    @staticmethod
    def register() -> None:
        """Register the necessary items related to the experiment."""


def import_module(name: str) -> ExperimentModule:
    """Imports a module given a name."""
    return importlib.import_module(name)  # type: ignore


def iter_experiment_modules(package: ModuleType) -> Iterator[str]:
    """Yields experiment module names from a given package."""
    return (
        f"{package.__name__}.{module.name}"
        for module in iter_modules(package.__path__)
        if not module.ispkg
    )


def register_experiments() -> None:
    """Loads and registers all experiment plugins."""
    modules_to_register = [
        *iter_experiment_modules(catalog),
        *iter_experiment_modules(wip),
    ]

    for module_name in modules_to_register:
        if find_spec(module_name):
            experiment_module = import_module(module_name)
            experiment_module.register()
        else:
            print(f"Warning: Module {module_name} not found")


if __name__ == "__main__":
    register_experiments()
