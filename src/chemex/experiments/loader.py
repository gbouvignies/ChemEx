"""A simple plugin loader for experiment modules."""

from __future__ import annotations

import importlib
import os
from collections.abc import Iterator
from importlib.util import find_spec
from pkgutil import iter_modules
from types import ModuleType
from typing import Protocol, cast

from chemex.experiments import catalog, experiment_types
from chemex.experiments.catalog import wip

WIP_ENV_VAR = "CHEMEX_ENABLE_WIP"
LEGACY_WIP_EXPERIMENT_NAMES: frozenset[str] = frozenset({"cest_15n_test"})
_TRUTHY_ENV_VALUES = frozenset({"1", "true", "yes", "on"})


class ExperimentModule(Protocol):
    """Represents an experiment module interface."""

    @staticmethod
    def register() -> None:
        """Register the necessary items related to the experiment."""


def import_module(name: str) -> ExperimentModule:
    """Imports a module given a name."""
    module = importlib.import_module(name)
    if not callable(getattr(module, "register", None)):
        msg = f"Experiment module {name!r} does not define a callable register()"
        raise ImportError(msg)  # noqa: TRY004 - plugin contract makes this an import failure
    return cast("ExperimentModule", module)


def iter_experiment_modules(package: ModuleType) -> Iterator[str]:
    """Yields experiment module names from a given package."""
    return (
        f"{package.__name__}.{module.name}"
        for module in iter_modules(package.__path__)
        if not module.ispkg
    )


def wip_experiment_names() -> frozenset[str]:
    """Return canonical Experiment Type names under ``catalog.wip`` without importing."""
    return frozenset(
        f"wip.{module.name}"
        for module in iter_modules(wip.__path__)
        if not module.ispkg
    )


def wip_enabled_from_env() -> bool:
    """Return whether experimental WIP plugins are enabled via environment."""
    raw = os.environ.get(WIP_ENV_VAR)
    if raw is None:
        return False
    return raw.strip().lower() in _TRUTHY_ENV_VALUES


def register_experiments(*, include_wip: bool | None = None) -> None:
    """Load and register supported catalog plugins (and WIP plugins on explicit opt-in)."""
    enable_wip = wip_enabled_from_env() if include_wip is None else include_wip

    for module_name in iter_experiment_modules(catalog):
        if find_spec(module_name):
            experiment_module = import_module(module_name)
            experiment_module.register()

    if enable_wip:
        for module_name in iter_experiment_modules(wip):
            if find_spec(module_name):
                experiment_module = import_module(module_name)
                experiment_module.register()
    else:
        experiment_types.unregister_wip_experiment_types()


if __name__ == "__main__":
    register_experiments()
