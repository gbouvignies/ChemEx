from __future__ import annotations

from typing import Protocol

from chemex.configuration.methods import Method
from chemex.experiments.loader import register_experiments
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec, ModelState
from chemex.parameters.database import (
    ModelReader,
    ParameterStore,
    create_parameter_store,
)
from chemex.parameters.factory import ParameterFactory
from chemex.parameters.legacy_adapter import LegacyValuesAdapter
from chemex.parameters.parameterization import (
    ActiveParameterization,
    build_initial_analysis_values,
    compile_active_parameterization,
)
from chemex.parameters.values import AnalysisValues
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
        analysis_values: AnalysisValues | None = None,
        execution: ExecutionSettings | None = None,
    ) -> None:
        self.model = ModelState() if model is None else model
        self.analysis_values = (
            AnalysisValues() if analysis_values is None else analysis_values
        )
        self.legacy_values_adapter = LegacyValuesAdapter(self.analysis_values)
        self.parameters = (
            create_parameter_store(self.model) if parameters is None else parameters
        )
        if isinstance(self.parameters, ParameterStore):
            self.parameters.attach_legacy_values_adapter(self.legacy_values_adapter)
        self.parameter_factory = (
            ParameterFactory(self.parameters)
            if parameter_factory is None
            else parameter_factory
        )
        self.execution = ExecutionSettings() if execution is None else execution

    @classmethod
    def create(cls) -> AnalysisSession:
        """Create a fresh analysis session with isolated runtime state."""
        ensure_plugins_registered()
        session = cls()
        session.reset()
        return session

    def reset(self) -> None:
        """Clear cached runtime state before starting a new analysis."""
        self.parameter_factory.reset()
        self.parameters.reset()
        self.analysis_values.reset()
        self.legacy_values_adapter.reset()
        self.model.reset()

    def try_build_analysis_values(self) -> bool:
        """Seal configuration and initialize the non-authoritative native values."""
        if not self.parameter_factory.try_seal_configuration():
            return False
        configuration = self.parameter_factory.sealed_configuration
        parameter_model = self.parameter_factory.sealed_parameter_model
        if configuration is None or parameter_model is None:
            error = RuntimeError("Sealed native parameter model is unavailable")
            self.parameter_factory.disable_native_candidate(error)
            self.legacy_values_adapter.disable(error)
            return False
        model_identity = self.model.spec.identity
        try:
            initial_values = (
                build_initial_analysis_values(parameter_model)
                if any(item.effective_value is None for item in configuration)
                else None
            )
            self.analysis_values.initialize(
                model_identity,
                configuration,
                _native_initial_values=initial_values,
            )
        except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
            self.parameter_factory.disable_native_candidate(error)
            self.legacy_values_adapter.disable(error)
            return False
        return True

    def set_model(self, name: str) -> None:
        """Set the active kinetics model for the current session."""
        ensure_plugins_registered()
        self.parameter_factory.clear_cache()
        self.model.set_model(name)

    def compile_parameterization(
        self,
        method: Method,
        required_ids: set[str],
    ) -> ActiveParameterization:
        """Compile one native preview without changing authoritative state."""
        parameter_model = self.parameter_factory.sealed_parameter_model
        if parameter_model is None:
            msg = "The sealed native parameter model is unavailable"
            raise RuntimeError(msg)
        return compile_active_parameterization(
            parameter_model,
            self.analysis_values.snapshot(),
            method,
            required_ids,
        )

    def try_compile_parameterization(
        self,
        method: Method,
        required_ids: set[str],
    ) -> ActiveParameterization | None:
        """Best-effort native preview that cannot veto the legacy occurrence."""
        try:
            candidate = self.compile_parameterization(method, required_ids)
            frame = candidate.frame_from_snapshot(self.analysis_values.snapshot())
            candidate.resolve(frame)
        except Exception:  # noqa: BLE001 - checkpoint-1 isolation boundary
            return None
        return candidate


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
