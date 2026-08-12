"""Typed product-facing provenance for native ChemEx occurrences."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import tomllib
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Literal

import numpy as np

from chemex import __version__
from chemex.configuration.conditions import Conditions
from chemex.evaluation.native import EvaluationPlan
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ActiveParameterization,
    SealedParameterModel,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.runtime import ExecutionSettings


class NativeProvenanceError(ValueError):
    """Raised when product-facing provenance is incomplete or ambiguous."""


type StepLifecycle = Literal[
    "committed",
    "successful_no_state_change",
    "failed",
    "accepted_uncommitted",
    "cancelled",
    "interrupted",
]
type StepAuthority = Literal["committed_fit", "evaluation_only", "diagnostic_only"]
type ArtifactRole = Literal[
    "committed_restart_state",
    "report_only_fitted_values",
    "diagnostic_provenance",
    "partial_evidence",
    "product_output",
]
_BASELINE_REFERENCE_KEY = object()
_WORKFLOW_PROVENANCE_KEY = object()


def _require_text(value: str, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise NativeProvenanceError(f"{field_name} cannot be empty")
    return text


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "schema_version": 1, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _toml_key(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def _parameter_float(value: float) -> str:
    scalar = float(value)
    if math.isnan(scalar):
        raise NativeProvenanceError("Independent parameter state cannot contain NaN")
    return repr(0.0 if scalar == 0.0 else scalar)


def _normalized_executed_method(
    parameterization: ActiveParameterization,
    plan: EvaluationPlan,
) -> str:
    """Serialize the realized method semantics, not caller-authored selectors."""
    lines = [
        "FORMAT_VERSION = 2",
        "",
        "[RESOLVED]",
        f"MODEL = {_toml_key(parameterization.binder.model_name)}",
        f"PARAMETERIZATION_IDENTITY = {_toml_key(parameterization.identity)}",
        f"CONSTRAINT_PROGRAM_IDENTITY = {_toml_key(parameterization.program.fingerprint)}",
        f"EVALUATION_PLAN_IDENTITY = {_toml_key(plan.identity)}",
        "",
    ]
    for param_id in parameterization.scope_ids:
        lines.extend(
            (
                "[[RESOLVED.PARAMETERS]]",
                f"ID = {_toml_key(param_id)}",
                f"ROLE = {_toml_key(parameterization.role(param_id).value)}",
                "",
            )
        )
    for constraint in parameterization.program.constraints:
        lines.extend(
            (
                "[[RESOLVED.CONSTRAINTS]]",
                f"TARGET = {_toml_key(constraint.target_id)}",
                f"EXPRESSION = {_toml_key(constraint.expression_text)}",
                f"SOURCE = {_toml_key(constraint.source)}",
                f"DEPENDENCIES = {json.dumps(list(constraint.dependencies))}",
                "",
            )
        )
    return "\n".join(lines)


def serialize_independent_parameters(
    parameter_model: SealedParameterModel,
    independent_ids: tuple[str, ...],
    snapshot: AnalysisValuesSnapshot,
    *,
    state_kind: Literal["starting", "committed"],
) -> str:
    """Render complete independent values and bounds in ChemEx input form."""
    if (
        snapshot.model_identity != parameter_model.model_identity
        or snapshot.definitions_identity != parameter_model.definitions.identity
        or snapshot.configuration_identity != parameter_model.configuration.identity
    ):
        raise NativeProvenanceError(
            "Independent parameter snapshot differs from the sealed parameter model"
        )
    if len(set(independent_ids)) != len(independent_ids) or any(
        param_id not in parameter_model.configuration for param_id in independent_ids
    ):
        raise NativeProvenanceError(
            "Independent parameter scope must be unique and completely configured"
        )
    sections: dict[str, list[tuple[str, float, float, float]]] = {}
    for param_id in independent_ids:
        definition = parameter_model.definitions[param_id]
        name = ParamName(
            definition.name,
            SpinSystem.from_name(definition.spin_system_name),
            Conditions.model_construct(None, **dict(definition.condition_entries)),
        )
        config = parameter_model.configuration[param_id]
        if name.spin_system:
            section = name.section
            key = str(name.spin_system)
        else:
            section = "GLOBAL"
            key = name.section_res
        sections.setdefault(section, []).append(
            (key, snapshot[param_id], config.lower_bound, config.upper_bound)
        )
    if state_kind == "starting":
        heading = (
            "# Complete starting independent parameter state used by this invocation.",
            "# This is provenance for the original start, not a committed restart.",
        )
    else:
        heading = (
            "# Complete committed independent parameter state.",
            "# This file is intended to start a new ChemEx invocation.",
            "# Report-only fitted values are published separately in fitted.toml.",
        )
    lines = [*heading, "# Each array is [value, minimum, maximum].", ""]
    for section, items in sections.items():
        lines.append(f"[{_toml_key(section)}]")
        lines.extend(
            f"{_toml_key(key)} = [{_parameter_float(value)}, "
            f"{_parameter_float(lower)}, {_parameter_float(upper)}]"
            for key, value, lower, upper in items
        )
        lines.append("")
    return "\n".join(lines)


@dataclass(frozen=True, slots=True)
class ProfileSelectionRecord:
    """Model and profile selection applied to one native workflow."""

    model_name: str
    include: tuple[str, ...]
    exclude: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "model_name", _require_text(self.model_name, "model"))
        for field_name in ("include", "exclude"):
            values = tuple(
                _require_text(item, field_name) for item in getattr(self, field_name)
            )
            if len(set(values)) != len(values):
                raise NativeProvenanceError(
                    f"Profile {field_name} selection contains duplicates"
                )
            object.__setattr__(self, field_name, values)

    def to_record(self) -> dict[str, object]:
        return {
            "model": self.model_name,
            "include": list(self.include),
            "exclude": list(self.exclude),
        }


@dataclass(frozen=True, slots=True)
class PolicyRecord:
    """Named immutable policy identity used by a workflow."""

    name: str
    identity: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _require_text(self.name, "policy name"))
        object.__setattr__(
            self, "identity", _require_text(self.identity, "policy identity")
        )


@dataclass(frozen=True, slots=True)
class BudgetRecord:
    """Resolved non-fungible budget and its authoritative usage."""

    name: str
    limit: int
    used: int | None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _require_text(self.name, "budget name"))
        if (
            isinstance(self.limit, bool)
            or self.limit < 1
            or isinstance(self.used, bool)
            or (self.used is not None and (self.used < 0 or self.used > self.limit))
        ):
            raise NativeProvenanceError(
                "Budget limits must be positive and known usage within the limit"
            )


@dataclass(frozen=True, slots=True)
class SeedRecord:
    """Resolved unsigned seed and the derivation policy that owns it."""

    name: str
    value: int
    policy_identity: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _require_text(self.name, "seed name"))
        object.__setattr__(
            self,
            "policy_identity",
            _require_text(self.policy_identity, "seed policy identity"),
        )
        if isinstance(self.value, bool) or not 0 <= self.value < 2**64:
            raise NativeProvenanceError("Seeds must be unsigned 64-bit integers")


@dataclass(frozen=True, slots=True)
class BaselineReference:
    """Reference to baseline provenance without transferring its authority."""

    _construction_key: object = field(repr=False, compare=False)
    kind: Literal["occurrence", "bundle", "manifest"]
    identity: str

    def __post_init__(self) -> None:
        if self._construction_key is not _BASELINE_REFERENCE_KEY:
            raise NativeProvenanceError(
                "Baseline references must come from validated baseline records"
            )
        if self.kind not in ("occurrence", "bundle", "manifest"):
            raise NativeProvenanceError("Unsupported baseline reference kind")
        object.__setattr__(
            self,
            "identity",
            _require_text(self.identity, "baseline reference identity"),
        )

    @classmethod
    def from_occurrence(cls, occurrence: object) -> BaselineReference:
        """Reference an occurrence only after canonical record revalidation."""
        from chemex.baselines import Occurrence

        if not isinstance(occurrence, Occurrence):
            raise TypeError("Baseline occurrence reference requires Occurrence")
        canonical = Occurrence.from_record(occurrence.to_record())
        return cls(_BASELINE_REFERENCE_KEY, "occurrence", canonical.identity)

    @classmethod
    def from_result_bundle(cls, bundle: object) -> BaselineReference:
        """Reference a result bundle only after canonical record revalidation."""
        from chemex.baselines import ResultBundle

        if not isinstance(bundle, ResultBundle):
            raise TypeError("Baseline bundle reference requires ResultBundle")
        canonical = ResultBundle.from_record(bundle.to_record())
        return cls(_BASELINE_REFERENCE_KEY, "bundle", canonical.identity)

    @classmethod
    def from_result_manifest(cls, bundle: object) -> BaselineReference:
        """Reference a result manifest only after canonical bundle revalidation."""
        from chemex.baselines import ResultBundle

        if not isinstance(bundle, ResultBundle):
            raise TypeError("Baseline manifest reference requires ResultBundle")
        canonical = ResultBundle.from_record(bundle.to_record())
        return cls(_BASELINE_REFERENCE_KEY, "manifest", canonical.manifest_identity)


@dataclass(frozen=True, slots=True)
class ProvenanceEnvironment:
    """Resolved software and numerical-library environment."""

    chemex_version: str
    python_version: str
    python_implementation: str
    platform: str
    numpy_version: str
    scipy_version: str
    emcee_version: str
    numerical_libraries: tuple[tuple[str, str, str], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        names = (
            "chemex_version",
            "python_version",
            "python_implementation",
            "platform",
            "numpy_version",
            "scipy_version",
            "emcee_version",
        )
        for name in names:
            object.__setattr__(self, name, _require_text(getattr(self, name), name))
        libraries = tuple(
            tuple(_require_text(value, "numerical library field") for value in item)
            for item in self.numerical_libraries
        )
        if any(len(item) != 3 for item in libraries):
            raise NativeProvenanceError(
                "Numerical libraries require kind, name, and version"
            )
        if len({item[0] for item in libraries}) != len(libraries):
            raise NativeProvenanceError("Numerical library kinds must be unique")
        object.__setattr__(self, "numerical_libraries", libraries)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-provenance-environment",
                tuple(getattr(self, name) for name in names) + (libraries,),
            ),
        )

    @classmethod
    def from_current_process(cls) -> ProvenanceEnvironment:
        """Capture the local product and numerical runtime without lmfit fields."""

        def package_version(name: str) -> str:
            try:
                return version(name)
            except PackageNotFoundError:
                return "unavailable"

        build_dependencies = getattr(np.__config__, "CONFIG", {}).get(
            "Build Dependencies", {}
        )
        libraries = tuple(
            (
                kind,
                str(details.get("name", "unknown")),
                str(details.get("version", "unknown")),
            )
            for kind, details in sorted(build_dependencies.items())
            if isinstance(details, dict) and details.get("found", False)
        )
        return cls(
            chemex_version=__version__,
            python_version=platform.python_version(),
            python_implementation=platform.python_implementation(),
            platform=platform.platform(),
            numpy_version=package_version("numpy"),
            scipy_version=package_version("scipy"),
            emcee_version=package_version("emcee"),
            numerical_libraries=libraries,
        )

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "chemex_version": self.chemex_version,
            "python_version": self.python_version,
            "python_implementation": self.python_implementation,
            "platform": self.platform,
            "numpy_version": self.numpy_version,
            "scipy_version": self.scipy_version,
            "emcee_version": self.emcee_version,
            "numerical_libraries": [
                {"kind": kind, "name": name, "version": library_version}
                for kind, name, library_version in self.numerical_libraries
            ],
        }


@dataclass(frozen=True, slots=True)
class WorkflowProvenance:
    """Complete static and resolved provenance for one native method workflow."""

    _construction_key: object = field(repr=False, compare=False)
    normalized_method_text: str
    selection: ProfileSelectionRecord
    parameterization_identity: str
    evaluation_plan_identity: str
    policies: tuple[PolicyRecord, ...]
    budgets: tuple[BudgetRecord, ...]
    seeds: tuple[SeedRecord, ...]
    execution: ExecutionSettings
    environment: ProvenanceEnvironment
    baseline_references: tuple[BaselineReference, ...] = ()
    workflow_identity: str = field(init=False)
    method_identity: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self._construction_key is not _WORKFLOW_PROVENANCE_KEY:
            raise NativeProvenanceError(
                "Workflow provenance must be canonically constructed"
            )
        method_text = self.normalized_method_text
        if not method_text.strip():
            raise NativeProvenanceError("Normalized method text cannot be empty")
        try:
            method = tomllib.loads(method_text)
        except tomllib.TOMLDecodeError as error:
            raise NativeProvenanceError(
                "Normalized method text must be TOML"
            ) from error
        format_version = next(
            (value for key, value in method.items() if key.lower() == "format_version"),
            None,
        )
        if format_version != 2:
            raise NativeProvenanceError(
                "Normalized native method must declare FORMAT_VERSION = 2"
            )
        for name, records in (
            ("policy", self.policies),
            ("budget", self.budgets),
            ("seed", self.seeds),
        ):
            record_names = tuple(item.name for item in records)
            if len(set(record_names)) != len(record_names):
                raise NativeProvenanceError(f"{name.title()} names must be unique")
        baseline_keys = tuple(
            (item.kind, item.identity) for item in self.baseline_references
        )
        if len(set(baseline_keys)) != len(baseline_keys):
            raise NativeProvenanceError("Baseline references must be unique")
        method_identity = _identity(
            "normalized-method-v2",
            self.normalized_method_sha256,
        )
        workflow_identity = _identity(
            "native-method-step-workflow",
            {
                "method_identity": method_identity,
                "parameterization_identity": self.parameterization_identity,
                "evaluation_plan_identity": self.evaluation_plan_identity,
                "selection": self.selection.to_record(),
                "policies": [(item.name, item.identity) for item in self.policies],
                "budgets": [
                    (item.name, item.limit, item.used) for item in self.budgets
                ],
                "seeds": [
                    (item.name, item.value, item.policy_identity) for item in self.seeds
                ],
                "execution": (
                    self.execution.workers,
                    self.execution.native_threads,
                ),
            },
        )
        object.__setattr__(self, "method_identity", method_identity)
        object.__setattr__(self, "workflow_identity", workflow_identity)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-workflow-provenance",
                {
                    "workflow_identity": self.workflow_identity,
                    "method_identity": self.method_identity,
                    "normalized_method_sha256": self.normalized_method_sha256,
                    "parameterization_identity": self.parameterization_identity,
                    "evaluation_plan_identity": self.evaluation_plan_identity,
                    "selection": self.selection.to_record(),
                    "policies": [(item.name, item.identity) for item in self.policies],
                    "budgets": [
                        (item.name, item.limit, item.used) for item in self.budgets
                    ],
                    "seeds": [
                        (item.name, item.value, item.policy_identity)
                        for item in self.seeds
                    ],
                    "execution": (
                        self.execution.workers,
                        self.execution.native_threads,
                    ),
                    "environment_identity": self.environment.identity,
                    "baseline_references": baseline_keys,
                },
            ),
        )

    @classmethod
    def create(
        cls,
        *,
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
        policies: tuple[PolicyRecord, ...],
        budgets: tuple[BudgetRecord, ...],
        seeds: tuple[SeedRecord, ...],
        execution: ExecutionSettings,
        environment: ProvenanceEnvironment,
        baseline_references: tuple[BaselineReference, ...] = (),
    ) -> WorkflowProvenance:
        """Derive canonical method, selection, and workflow from executed inputs."""
        if (
            plan.parameterization_identity != parameterization.evaluator_identity
            or plan.constraint_program_identity != parameterization.program.fingerprint
        ):
            raise NativeProvenanceError(
                "Workflow provenance requires the executed plan and parameterization"
            )
        return cls(
            _WORKFLOW_PROVENANCE_KEY,
            _normalized_executed_method(parameterization, plan),
            ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            ),
            parameterization.identity,
            plan.identity,
            policies,
            budgets,
            seeds,
            execution,
            environment,
            baseline_references,
        )

    @property
    def normalized_method_sha256(self) -> str:
        return hashlib.sha256(self.normalized_method_text.encode()).hexdigest()

    def validate_execution_context(
        self,
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
    ) -> None:
        """Re-derive and validate every execution-owned method/selection field."""
        if (
            self.parameterization_identity != parameterization.identity
            or self.evaluation_plan_identity != plan.identity
            or self.normalized_method_text
            != _normalized_executed_method(parameterization, plan)
            or self.selection
            != ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            )
        ):
            raise NativeProvenanceError(
                "Workflow method and selection differ from executed artifacts"
            )

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.workflow_identity,
            "provenance_identity": self.identity,
            "method": {
                "identity": self.method_identity,
                "normalized_sha256": self.normalized_method_sha256,
                "parameterization_identity": self.parameterization_identity,
                "evaluation_plan_identity": self.evaluation_plan_identity,
            },
            "selection": self.selection.to_record(),
            "policies": [
                {"name": item.name, "identity": item.identity} for item in self.policies
            ],
            "budgets": [
                {"name": item.name, "limit": item.limit, "used": item.used}
                for item in self.budgets
            ],
            "seeds": [
                {
                    "name": item.name,
                    "value": item.value,
                    "policy_identity": item.policy_identity,
                }
                for item in self.seeds
            ],
            "execution": {
                "workers": self.execution.workers,
                "native_threads": self.execution.native_threads,
            },
            "environment": self.environment.to_record(),
            "baseline_references": [
                {"kind": item.kind, "identity": item.identity}
                for item in self.baseline_references
            ],
        }


@dataclass(frozen=True, slots=True)
class ArtifactReference:
    """Content-addressed product artifact linked by a step manifest."""

    path: str
    role: ArtifactRole
    sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", _require_text(self.path, "artifact path"))
        object.__setattr__(self, "role", _require_text(self.role, "artifact role"))
        digest = self.sha256.lower()
        if len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        ):
            raise NativeProvenanceError(
                "Artifact SHA-256 must be lowercase hexadecimal"
            )
        object.__setattr__(self, "sha256", digest)


@dataclass(frozen=True, slots=True)
class PublishedStepReference:
    """Verified link from invocation provenance to one atomic step publication."""

    path: Path
    lifecycle: StepLifecycle
    authority: StepAuthority
    manifest_identity: str
    manifest_sha256: str
    provenance: WorkflowProvenance
    independent_ids: tuple[str, ...]
    artifacts: tuple[ArtifactReference, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path).resolve())
        for field_name in (
            "lifecycle",
            "authority",
            "manifest_identity",
            "manifest_sha256",
        ):
            object.__setattr__(
                self, field_name, _require_text(getattr(self, field_name), field_name)
            )
        if len(self.manifest_sha256) != 64:
            raise NativeProvenanceError("Step manifest SHA-256 is malformed")
        if len(set(self.independent_ids)) != len(self.independent_ids):
            raise NativeProvenanceError("Published independent scope has duplicates")
