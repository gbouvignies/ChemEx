"""Plugin loader for B1 distribution modules."""

import importlib
from collections.abc import Iterator
from importlib.util import find_spec
from pkgutil import iter_modules
from types import ModuleType
from typing import Protocol

from chemex.nmr import distributions


class DistributionModule(Protocol):
    """Represents a distribution module interface."""

    @staticmethod
    def register() -> None:
        """Register the distribution generator and config class."""


def import_module(name: str) -> DistributionModule:
    """Import a module given a name."""
    return importlib.import_module(name)  # type: ignore


def iter_distribution_modules(package: ModuleType) -> Iterator[str]:
    """Yield distribution module names from the package."""
    return (
        f"{package.__name__}.{module.name}"
        for module in iter_modules(package.__path__)
        if not module.ispkg and module.name not in ("__init__", "loader", "registry")
    )


def register_distributions() -> None:
    """Load and register all distribution plugins."""
    for module_name in iter_distribution_modules(distributions):
        if find_spec(module_name):
            distribution_module = import_module(module_name)
            distribution_module.register()


if __name__ == "__main__":
    register_distributions()
