"""A simple plugin loader."""

import importlib
from collections.abc import Callable
from pkgutil import iter_modules
from types import ModuleType

from chemex.models import kinetic


def _load_register(name: str) -> Callable[[], None]:
    """Load the registration callable for a kinetic settings module."""
    module: ModuleType = importlib.import_module(name)
    register = getattr(module, "register", None)
    if callable(register):
        return register
    msg = f"Module {name!r} does not define a callable 'register' function"
    raise TypeError(msg)


def register_kinetic_settings() -> None:
    """Loads the plugins defined in the plugins list."""
    for module in iter_modules(kinetic.__path__):
        module_name = f"{kinetic.__name__}.{module.name}"
        register = _load_register(module_name)
        register()
