from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Protocol

from chemex.experiments.loader import register_experiments
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import model
from chemex.parameters import database
from chemex.parameters.factory import clear_settings_cache
from chemex.typing import Array

_PLUGIN_STATE = {"registered": False}


class ModelController(Protocol):
    def reset(self) -> None: ...

    def set_model(self, name: str) -> None: ...


class ParameterStore(Protocol):
    def reset_parameters(self) -> None: ...

    def build_lmfit_params(self, param_ids: Iterable[str]) -> object: ...

    def update_from_parameters(self, parameters: object) -> None: ...

    def set_param_values(self, par_values: dict[str, float]) -> None: ...

    def set_param_defaults(self, defaults: object) -> None: ...

    def fix_all_parameters(self) -> None: ...

    def sort_parameters(self) -> None: ...

    def set_parameter_status(self, method: object) -> None: ...

    def parse_grid(self, grid_entries: list[str]) -> dict[str, Array]: ...

    def get_parameters(self, param_ids: Iterable[str]) -> dict[str, object]: ...


@dataclass
class AnalysisSession:
    """Compatibility shell around the current global runtime state."""

    model: ModelController = field(default_factory=lambda: model)
    parameters: ParameterStore = field(default_factory=lambda: database)

    @classmethod
    def create(cls) -> AnalysisSession:
        """Create a fresh analysis session bound to the global backends."""
        ensure_plugins_registered()
        session = cls()
        session.reset()
        return session

    def reset(self) -> None:
        """Clear cached runtime state before starting a new analysis."""
        clear_settings_cache()
        self.parameters.reset_parameters()
        self.model.reset()

    def set_model(self, name: str) -> None:
        """Set the active kinetics model for the current session."""
        ensure_plugins_registered()
        clear_settings_cache()
        self.model.set_model(name)


def ensure_plugins_registered() -> None:
    """Register model and experiment plugins once per process."""
    if _PLUGIN_STATE["registered"]:
        return

    register_kinetic_settings()
    register_experiments()
    _PLUGIN_STATE["registered"] = True


def reset_plugin_registration() -> None:
    """Mark plugin registration as pending for a fresh session bootstrap."""
    _PLUGIN_STATE["registered"] = False
