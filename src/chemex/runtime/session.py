from __future__ import annotations

from typing import Protocol

from chemex.configuration.parameters import DefaultListType
from chemex.experiments.loader import register_experiments
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec, ModelState
from chemex.parameters.configuration import (
    ParameterConfiguration,
    build_parameter_configuration,
)
from chemex.parameters.database import (
    ModelReader,
    ParameterStore,
    create_parameter_store,
)
from chemex.parameters.definitions import (
    ParameterDefinitions,
    ParameterDefinitionsCollector,
)
from chemex.parameters.factory import ParameterFactory
from chemex.runtime.execution import ExecutionSettings

_PLUGIN_STATE = {"registered": False}


class ModelController(ModelReader, Protocol):
    @property
    def spec(self) -> ModelSpec: ...

    def reset(self) -> None: ...

    def set_model(self, name: str) -> None: ...


class AnalysisSession:
    """Session-scoped runtime state for one analysis run."""

    def __init__(
        self,
        *,
        model: ModelController | None = None,
        parameters: ParameterStore | None = None,
        parameter_factory: ParameterFactory | None = None,
        execution: ExecutionSettings | None = None,
    ) -> None:
        self.model = ModelState() if model is None else model
        self.parameters = (
            create_parameter_store(self.model) if parameters is None else parameters
        )
        self.parameter_factory = (
            ParameterFactory(self.parameters)
            if parameter_factory is None
            else parameter_factory
        )
        if self.parameter_factory.definitions is None:
            self.parameter_factory.definitions = ParameterDefinitionsCollector()

        self.execution = ExecutionSettings() if execution is None else execution
        self.ordinary_definitions: ParameterDefinitions | None = None
        self.model_free_definitions: ParameterDefinitions | None = None
        self.ordinary_configuration: ParameterConfiguration | None = None
        self.model_free_configuration: ParameterConfiguration | None = None

    @classmethod
    def create(cls) -> AnalysisSession:
        """Create a fresh analysis session with isolated runtime state."""
        ensure_plugins_registered()
        session = cls()
        session.reset()
        return session

    def reset(self) -> None:
        """Clear cached runtime state before starting a new analysis."""
        self.parameter_factory.clear_cache()
        self.parameters.reset()
        self.model.reset()

    def set_model(self, name: str) -> None:
        """Set the active kinetics model for the current session."""
        ensure_plugins_registered()
        self.parameter_factory.clear_cache()
        self.model.set_model(name)

    def seal_configuration(self, defaults: DefaultListType) -> None:
        """Seal definitions and configurations after loading experiments and defaults."""
        definitions = self.parameter_factory.definitions
        if definitions is None:
            msg = "Cannot seal configuration: ParameterDefinitionsCollector is missing"
            raise RuntimeError(msg)
        ordinary_def, mf_def = definitions.seal()
        self.ordinary_definitions = ordinary_def
        self.model_free_definitions = mf_def

        self.ordinary_configuration = build_parameter_configuration(
            self.ordinary_definitions,
            self.parameters.get_catalog(model_free=False),
            defaults,
        )
        self.model_free_configuration = build_parameter_configuration(
            self.model_free_definitions,
            self.parameters.get_catalog(model_free=True),
            defaults,
        )


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
