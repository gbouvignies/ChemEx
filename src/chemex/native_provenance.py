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
from threading import RLock
from typing import Literal, cast
from weakref import ReferenceType, WeakKeyDictionary, ref

import numpy as np

from chemex import __version__
from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.evaluation.native import EvaluationPlan, EvaluationResult
from chemex.optimize.direct_trf import DirectTrfExecution, DirectTrfInvocation
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
    "no_objective_data",
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
_PUBLISHED_STEP_KEY = object()


def _require_text(value: str, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise NativeProvenanceError(f"{field_name} cannot be empty")
    return text


def _record_int(value: object, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise NativeProvenanceError(f"{field_name} must be an integer")
    return value


def _record_strings(value: object, field_name: str) -> tuple[str, ...]:
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise NativeProvenanceError(f"{field_name} must be a string list")
    return tuple(cast(list[str], value))


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "schema_version": 1, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def native_step_manifest_identity(record: dict[str, object]) -> str:
    """Compute the durable step identity over every canonical child field."""
    payload = dict(record)
    payload.pop("manifest_identity", None)
    return _identity("native-step-manifest-v2", payload)


def validate_native_step_manifest_bytes(  # noqa: C901 - closed durable schema
    content: bytes,
) -> dict[str, object]:
    """Recursively validate one durable lifecycle-aware native step manifest."""
    try:
        record = tomllib.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise NativeProvenanceError(
            "Native step manifest is not canonical TOML"
        ) from error
    allowed = {
        "schema_version",
        "manifest_identity",
        "lifecycle",
        "authority",
        "publication_occurrence_identity",
        "starting_state_kind",
        "starting_occurrence_identity",
        "starting_revision",
        "accepted_result_identity",
        "accepted_occurrence_identity",
        "commit_receipt_identity",
        "commit_operation_identity",
        "commit_occurrence_identity",
        "committed_occurrence_identity",
        "committed_revision",
        "problem_identity",
        "parameterization_identity",
        "restart_record_identity",
        "operation_identity",
        "workflow",
        "method",
        "selection",
        "execution",
        "environment",
        "policies",
        "budgets",
        "seeds",
        "baseline_references",
        "components",
        "evidence",
        "artifacts",
    }
    if set(record) - allowed or record.get("schema_version") != 1:
        raise NativeProvenanceError("Native step manifest has unsupported fields")
    lifecycle = record.get("lifecycle")
    authority = record.get("authority")
    common = {
        "schema_version",
        "manifest_identity",
        "lifecycle",
        "authority",
        "publication_occurrence_identity",
        "workflow",
        "method",
        "selection",
        "execution",
        "environment",
        "artifacts",
    }
    optional = {
        "policies",
        "budgets",
        "seeds",
        "baseline_references",
        "components",
        "evidence",
    }
    committed = {
        "starting_state_kind",
        "starting_occurrence_identity",
        "starting_revision",
        "accepted_result_identity",
        "accepted_occurrence_identity",
        "commit_receipt_identity",
        "commit_operation_identity",
        "commit_occurrence_identity",
        "committed_occurrence_identity",
        "committed_revision",
        "problem_identity",
        "parameterization_identity",
        "restart_record_identity",
    }
    if lifecycle == "committed":
        if authority != "committed_fit" or not common | committed <= set(record):
            raise NativeProvenanceError("Committed manifest lacks lifecycle authority")
        forbidden = {"operation_identity"}
    elif lifecycle in {"successful_no_state_change", "no_objective_data"}:
        if authority != "evaluation_only" or not common <= set(record):
            raise NativeProvenanceError("Evaluation manifest lacks lifecycle authority")
        forbidden = committed | {"operation_identity"}
    elif lifecycle in {"failed", "accepted_uncommitted", "cancelled", "interrupted"}:
        if authority != "diagnostic_only" or not common | {"operation_identity"} <= set(
            record
        ):
            raise NativeProvenanceError("Suppressed manifest lacks lifecycle authority")
        forbidden = committed - {
            "accepted_result_identity",
            "accepted_occurrence_identity",
        }
    else:
        raise NativeProvenanceError("Native step manifest lifecycle is invalid")
    if forbidden & set(record) or set(record) - (
        common | optional | committed | {"operation_identity"}
    ):
        raise NativeProvenanceError("Native step manifest mixes lifecycle fields")
    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise NativeProvenanceError("Native step manifest requires artifacts")
    roles: list[str] = []
    for path, artifact in artifacts.items():
        if (
            not isinstance(path, str)
            or Path(path).is_absolute()
            or ".." in Path(path).parts
            or not isinstance(artifact, dict)
            or set(artifact) != {"role", "sha256"}
        ):
            raise NativeProvenanceError("Native step manifest has malformed artifacts")
        ArtifactReference(path, artifact.get("role"), artifact.get("sha256"))  # type: ignore[arg-type]
        roles.append(str(artifact.get("role")))
    restart_count = roles.count("committed_restart_state")
    fitted_count = roles.count("report_only_fitted_values")
    if lifecycle == "committed":
        restart_paths = {
            path
            for path, artifact in artifacts.items()
            if isinstance(artifact, dict)
            and artifact.get("role") == "committed_restart_state"
        }
        if (
            restart_paths
            != {
                "Parameters/restart.toml",
                "Parameters/restart-provenance.json",
            }
            or fitted_count > 1
        ):
            raise NativeProvenanceError(
                "Committed manifest lacks exact restart linkage"
            )
    elif restart_count or fitted_count:
        raise NativeProvenanceError(
            "Non-committed manifest claims fitted/restart authority"
        )
    components = record.get("components", [])
    if not isinstance(components, list) or any(
        not isinstance(item, dict)
        or set(item)
        != {"identity", "disposition", "controlled_ids", "location", "authority"}
        or item.get("disposition")
        not in {
            "succeeded",
            "failed",
            "execution_failure",
            "cancelled",
            "interrupted",
            "not_started",
        }
        or item.get("location") != "Components/index.json"
        or item.get("authority") != "diagnostic_only"
        or item.get("location") not in artifacts
        or not isinstance(artifacts.get(item.get("location")), dict)
        or cast(dict[str, object], artifacts[item["location"]]).get("role")
        != "diagnostic_provenance"
        for item in components
    ):
        raise NativeProvenanceError("Native step components are malformed")
    component_records = cast(list[dict[str, object]], components)
    component_controls = tuple(
        _record_strings(item["controlled_ids"], "component controls")
        for item in component_records
    )
    if len({str(item["identity"]) for item in component_records}) != len(
        component_records
    ) or any(
        not controls or len(set(controls)) != len(controls)
        for controls in component_controls
    ):
        raise NativeProvenanceError("Native step component identities are invalid")
    evidence = record.get("evidence", [])
    if not isinstance(evidence, list) or any(
        not isinstance(item, dict)
        or set(item) != {"family", "identity", "validity", "location"}
        or item.get("location") not in artifacts
        for item in evidence
    ):
        raise NativeProvenanceError("Native step evidence is malformed")
    evidence_records = cast(list[dict[str, object]], evidence)
    if len({str(item["family"]) for item in evidence_records}) != len(evidence_records):
        raise NativeProvenanceError("Native step evidence families must be unique")
    WorkflowProvenance.from_manifest_record(record)
    if record.get("manifest_identity") != native_step_manifest_identity(record):
        raise NativeProvenanceError(
            "Native step manifest identity does not match children"
        )
    return record


def _toml_key(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def _parameter_float(value: float) -> str:
    scalar = float(value)
    if math.isnan(scalar):
        raise NativeProvenanceError("Independent parameter state cannot contain NaN")
    return repr(0.0 if scalar == 0.0 else scalar)


def _normalized_method_semantics(method: Method) -> str:
    """Serialize only canonical method semantics; execution facts live elsewhere."""
    payload = method.model_dump(
        mode="json",
        exclude={"include", "exclude"},
        exclude_none=True,
    )
    canonical = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    return "\n".join(
        (
            "FORMAT_VERSION = 2",
            "",
            "[METHOD]",
            f"CANONICAL_JSON = {_toml_key(canonical)}",
            "",
        )
    )


def _validate_normalized_method_semantics(text: str) -> Method:
    """Require the exact canonical schema-v2 normalized method representation."""
    try:
        record = tomllib.loads(text)
        if set(record) != {"FORMAT_VERSION", "METHOD"} or record["FORMAT_VERSION"] != 2:
            raise TypeError
        method_record = record["METHOD"]
        if not isinstance(method_record, dict) or set(method_record) != {
            "CANONICAL_JSON"
        }:
            raise TypeError
        canonical_json = method_record["CANONICAL_JSON"]
        if not isinstance(canonical_json, str):
            raise TypeError
        payload = json.loads(canonical_json)
        if not isinstance(payload, dict):
            raise TypeError
        method = Method.model_validate(payload)
    except (
        json.JSONDecodeError,
        KeyError,
        TypeError,
        ValueError,
        tomllib.TOMLDecodeError,
    ) as error:
        raise NativeProvenanceError(
            "Normalized method semantics are malformed"
        ) from error
    if _normalized_method_semantics(method) != text:
        raise NativeProvenanceError("Normalized method semantics are not canonical")
    return method


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
    occurrence_identity: str
    result_bundle_identity: str

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
        object.__setattr__(
            self,
            "occurrence_identity",
            _require_text(self.occurrence_identity, "baseline occurrence identity"),
        )
        object.__setattr__(
            self,
            "result_bundle_identity",
            _require_text(
                self.result_bundle_identity, "baseline result bundle identity"
            ),
        )
        for value in (
            self.identity,
            self.occurrence_identity,
            self.result_bundle_identity,
        ):
            if len(value) != 64 or any(
                character not in "0123456789abcdef" for character in value
            ):
                raise NativeProvenanceError(
                    "Baseline provenance identities must be SHA-256 values"
                )

    @classmethod
    def from_occurrence(cls, occurrence: object) -> BaselineReference:
        """Reference an occurrence only after canonical record revalidation."""
        from chemex.baselines import Occurrence

        if not isinstance(occurrence, Occurrence):
            raise TypeError("Baseline occurrence reference requires Occurrence")
        canonical = Occurrence.from_record(occurrence.to_record())
        if (
            canonical.lifecycle != "SUCCEEDED"
            or canonical.result_bundle_identity is None
        ):
            raise NativeProvenanceError(
                "Baseline provenance requires a successful occurrence"
            )
        requested_identity = Occurrence(
            canonical.execution_specification_identity,
            canonical.case_identity,
            canonical.actual_implementation_identity,
            canonical.lane_reference,
            canonical.lane_attestation_identity,
            canonical.input_member_identities,
            canonical.attempt_token,
        ).identity
        return cls(
            _BASELINE_REFERENCE_KEY,
            "occurrence",
            canonical.identity,
            requested_identity,
            canonical.result_bundle_identity,
        )

    @classmethod
    def from_result_bundle(cls, bundle: object) -> BaselineReference:
        """Reference a result bundle only after canonical record revalidation."""
        from chemex.baselines import ResultBundle

        if not isinstance(bundle, ResultBundle):
            raise TypeError("Baseline bundle reference requires ResultBundle")
        canonical = ResultBundle.from_record(bundle.to_record())
        return cls(
            _BASELINE_REFERENCE_KEY,
            "bundle",
            canonical.identity,
            canonical.occurrence_identity,
            canonical.identity,
        )

    @classmethod
    def from_result_manifest(cls, bundle: object) -> BaselineReference:
        """Reference a result manifest only after canonical bundle revalidation."""
        from chemex.baselines import ResultBundle

        if not isinstance(bundle, ResultBundle):
            raise TypeError("Baseline manifest reference requires ResultBundle")
        canonical = ResultBundle.from_record(bundle.to_record())
        return cls(
            _BASELINE_REFERENCE_KEY,
            "manifest",
            canonical.manifest_identity,
            canonical.occurrence_identity,
            canonical.identity,
        )


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
        """Capture the local product and native numerical runtime."""

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
    workflow_type: str
    grouping_topology: tuple[tuple[str, tuple[str, ...]], ...]
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
    execution_identity: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self._construction_key is not _WORKFLOW_PROVENANCE_KEY:
            raise NativeProvenanceError(
                "Workflow provenance must be canonically constructed"
            )
        _validate_normalized_method_semantics(self.normalized_method_text)
        object.__setattr__(
            self,
            "parameterization_identity",
            _require_text(self.parameterization_identity, "parameterization identity"),
        )
        object.__setattr__(
            self,
            "evaluation_plan_identity",
            _require_text(self.evaluation_plan_identity, "evaluation plan identity"),
        )
        object.__setattr__(
            self,
            "workflow_type",
            _require_text(self.workflow_type, "workflow type"),
        )
        grouping = tuple((name, tuple(ids)) for name, ids in self.grouping_topology)
        if len({name for name, _ids in grouping}) != len(grouping) or any(
            not name or not ids or len(set(ids)) != len(ids) for name, ids in grouping
        ):
            raise NativeProvenanceError("Workflow grouping topology is invalid")
        object.__setattr__(self, "grouping_topology", grouping)
        for name, records in (
            ("policy", self.policies),
            ("budget", self.budgets),
            ("seed", self.seeds),
        ):
            record_names = tuple(item.name for item in records)
            if len(set(record_names)) != len(record_names):
                raise NativeProvenanceError(f"{name.title()} names must be unique")
        baseline_keys = tuple(
            (
                item.kind,
                item.identity,
                item.occurrence_identity,
                item.result_bundle_identity,
            )
            for item in self.baseline_references
        )
        if (
            not {"occurrence", "bundle"}
            <= {item.kind for item in self.baseline_references}
            or len(self.baseline_references) not in {2, 3}
            or (
                len(self.baseline_references) == 3
                and {item.kind for item in self.baseline_references}
                != {"occurrence", "bundle", "manifest"}
            )
            or len({item.occurrence_identity for item in self.baseline_references}) != 1
            or len({item.result_bundle_identity for item in self.baseline_references})
            != 1
            or next(
                item.identity
                for item in self.baseline_references
                if item.kind == "occurrence"
            )
            != self.baseline_references[0].occurrence_identity
            or next(
                item.identity
                for item in self.baseline_references
                if item.kind == "bundle"
            )
            != self.baseline_references[0].result_bundle_identity
        ):
            raise NativeProvenanceError(
                "Baseline provenance requires one linked occurrence/bundle pair"
            )
        method_identity = _identity(
            "normalized-method-v2",
            self.normalized_method_sha256,
        )
        workflow_identity = _identity(
            "native-method-step-workflow",
            {
                "method_identity": method_identity,
                "workflow_type": self.workflow_type,
                "grouping_topology": self.grouping_topology,
            },
        )
        execution_identity = _identity(
            "native-execution-configuration-v2",
            {
                "workflow_identity": workflow_identity,
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
                "environment_identity": self.environment.identity,
            },
        )
        object.__setattr__(self, "method_identity", method_identity)
        object.__setattr__(self, "workflow_identity", workflow_identity)
        object.__setattr__(self, "execution_identity", execution_identity)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-workflow-provenance",
                {
                    "workflow_identity": self.workflow_identity,
                    "execution_identity": self.execution_identity,
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
        method: Method,
        invocation: DirectTrfInvocation | None,
        native_execution: DirectTrfExecution | EvaluationResult,
        policies: tuple[PolicyRecord, ...],
        budgets: tuple[BudgetRecord, ...],
        seeds: tuple[SeedRecord, ...],
        environment: ProvenanceEnvironment,
        baseline_references: tuple[BaselineReference, ...],
    ) -> WorkflowProvenance:
        """Derive canonical method, selection, and workflow from executed inputs."""
        if (
            plan.parameterization_identity != parameterization.evaluator_identity
            or plan.constraint_program_identity != parameterization.program.fingerprint
        ):
            raise NativeProvenanceError(
                "Workflow provenance requires the executed plan and parameterization"
            )
        workflow_type, grouping_topology, execution_settings = (
            cls._derive_live_execution_facts(
                parameterization,
                plan,
                invocation,
                native_execution,
            )
        )
        return cls(
            _WORKFLOW_PROVENANCE_KEY,
            _normalized_method_semantics(method),
            ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            ),
            workflow_type,
            grouping_topology,
            parameterization.identity,
            plan.identity,
            policies,
            budgets,
            seeds,
            execution_settings,
            environment,
            baseline_references,
        )

    @classmethod
    def create_method_step(
        cls,
        *,
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
        method: Method,
        semantic_workflow_identity: str,
        grouping_topology: tuple[tuple[str, tuple[str, ...]], ...],
        policies: tuple[PolicyRecord, ...],
        budgets: tuple[BudgetRecord, ...],
        seeds: tuple[SeedRecord, ...],
        execution: ExecutionSettings,
        environment: ProvenanceEnvironment,
        baseline_references: tuple[BaselineReference, ...],
    ) -> WorkflowProvenance:
        """Capture one already-compiled closed method-step composition."""
        if (
            not semantic_workflow_identity
            or plan.parameterization_identity != parameterization.evaluator_identity
            or plan.constraint_program_identity != parameterization.program.fingerprint
        ):
            raise NativeProvenanceError(
                "Method-step provenance requires its exact compiled native context"
            )
        return cls(
            _WORKFLOW_PROVENANCE_KEY,
            _normalized_method_semantics(method),
            ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            ),
            f"method_step:{semantic_workflow_identity}",
            grouping_topology,
            parameterization.identity,
            plan.identity,
            policies,
            budgets,
            seeds,
            execution,
            environment,
            baseline_references,
        )

    def validate_method_step_context(
        self,
        *,
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
        method: Method,
        semantic_workflow_identity: str,
        grouping_topology: tuple[tuple[str, tuple[str, ...]], ...],
        execution: ExecutionSettings,
    ) -> None:
        """Revalidate an exact closed method-step context for publication."""
        if (
            self.parameterization_identity != parameterization.identity
            or self.evaluation_plan_identity != plan.identity
            or self.normalized_method_text != _normalized_method_semantics(method)
            or self.selection
            != ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            )
            or self.workflow_type != f"method_step:{semantic_workflow_identity}"
            or self.grouping_topology != grouping_topology
            or self.execution != execution
        ):
            raise NativeProvenanceError(
                "Method-step provenance differs from its compiled execution context"
            )

    @staticmethod
    def _derive_live_execution_facts(
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
        invocation: DirectTrfInvocation | None,
        native_execution: DirectTrfExecution | EvaluationResult,
    ) -> tuple[
        str,
        tuple[tuple[str, tuple[str, ...]], ...],
        ExecutionSettings,
    ]:
        """Extract publication facts from the exact native occurrence artifacts."""
        grouping = (("aggregate", parameterization.independent_ids),)
        if isinstance(native_execution, DirectTrfExecution):
            if (
                invocation is None
                or native_execution.problem_identity != invocation.problem_identity
                or native_execution.invocation_identity != invocation.identity
            ):
                raise NativeProvenanceError(
                    "Direct TRF provenance requires its exact invocation and execution"
                )
            return "direct_trf", grouping, invocation.execution_settings
        if (
            invocation is not None
            or native_execution.plan_identity != plan.identity
            or native_execution.parameterization_identity
            != parameterization.evaluator_identity
        ):
            raise NativeProvenanceError(
                "Evaluation provenance requires its exact evaluation result"
            )
        # Evaluation-only has no worker dispatcher or optimizer invocation. These are
        # therefore implementation facts, not publication-time configuration claims.
        return "evaluation_only", grouping, ExecutionSettings()

    @property
    def normalized_method_sha256(self) -> str:
        return hashlib.sha256(self.normalized_method_text.encode()).hexdigest()

    def validate_execution_context(
        self,
        parameterization: ActiveParameterization,
        plan: EvaluationPlan,
        method: Method,
        invocation: DirectTrfInvocation | None,
        native_execution: DirectTrfExecution | EvaluationResult,
    ) -> None:
        """Re-derive and validate every execution-owned method/selection field."""
        workflow_type, grouping_topology, execution_settings = (
            self._derive_live_execution_facts(
                parameterization,
                plan,
                invocation,
                native_execution,
            )
        )
        if (
            self.parameterization_identity != parameterization.identity
            or self.evaluation_plan_identity != plan.identity
            or self.normalized_method_text != _normalized_method_semantics(method)
            or self.selection
            != ProfileSelectionRecord(
                parameterization.binder.model_name,
                tuple(profile.identity for profile in plan.profiles),
                (),
            )
            or self.workflow_type != workflow_type
            or self.grouping_topology != grouping_topology
            or self.execution != execution_settings
        ):
            raise NativeProvenanceError(
                "Workflow method, selection, or execution provenance differs from "
                "executed artifacts"
            )

    @classmethod
    def from_manifest_record(cls, manifest: dict[str, object]) -> WorkflowProvenance:
        """Rebuild durable provenance without recreating live execution authority."""
        try:
            raw_children = tuple(
                manifest[name]
                for name in (
                    "workflow",
                    "method",
                    "selection",
                    "execution",
                    "environment",
                )
            )
            if not all(isinstance(item, dict) for item in raw_children):
                raise TypeError
            workflow, method, selection, execution, environment = (
                cast(dict[str, object], item) for item in raw_children
            )
            if (
                set(workflow)
                != {
                    "identity",
                    "execution_identity",
                    "provenance_identity",
                    "type",
                    "groups",
                }
                or set(method)
                != {
                    "identity",
                    "normalized",
                    "normalized_sha256",
                    "parameterization_identity",
                    "evaluation_plan_identity",
                }
                or set(selection) != {"model", "include", "exclude"}
                or set(execution) != {"workers", "native_threads"}
                or set(environment)
                != {
                    "identity",
                    "chemex_version",
                    "python_version",
                    "python_implementation",
                    "platform",
                    "numpy_version",
                    "scipy_version",
                    "emcee_version",
                    "numerical_libraries",
                }
            ):
                raise TypeError
            groups = workflow["groups"]
            if not isinstance(groups, list) or any(
                not isinstance(item, dict)
                or set(item) != {"identity", "controlled_ids"}
                for item in groups
            ):
                raise TypeError
            group_records = cast(list[dict[str, object]], groups)
            for name, fields in (
                ("policies", {"name", "identity"}),
                ("seeds", {"name", "value", "policy_identity"}),
                (
                    "baseline_references",
                    {
                        "kind",
                        "identity",
                        "occurrence_identity",
                        "result_bundle_identity",
                    },
                ),
            ):
                children = manifest.get(name, [])
                if not isinstance(children, list) or any(
                    not isinstance(item, dict) or set(item) != fields
                    for item in children
                ):
                    raise TypeError
            budgets = manifest.get("budgets", [])
            if not isinstance(budgets, list) or any(
                not isinstance(item, dict)
                or not {"name", "limit"} <= set(item)
                or set(item) - {"name", "limit", "used"}
                for item in budgets
            ):
                raise TypeError
            libraries = environment.get("numerical_libraries", [])
            if not isinstance(libraries, list) or any(
                not isinstance(item, dict) or set(item) != {"kind", "name", "version"}
                for item in libraries
            ):
                raise TypeError
            library_records = cast(list[dict[str, object]], libraries)
            policy_records = cast(list[dict[str, object]], manifest.get("policies", []))
            budget_records = cast(list[dict[str, object]], manifest.get("budgets", []))
            seed_records = cast(list[dict[str, object]], manifest.get("seeds", []))
            baseline_records = cast(
                list[dict[str, object]], manifest.get("baseline_references", [])
            )
            provenance = cls(
                _WORKFLOW_PROVENANCE_KEY,
                str(method["normalized"]),
                ProfileSelectionRecord(
                    str(selection["model"]),
                    _record_strings(selection["include"], "selection include"),
                    _record_strings(selection["exclude"], "selection exclude"),
                ),
                str(workflow["type"]),
                tuple(
                    (
                        str(item["identity"]),
                        _record_strings(item["controlled_ids"], "controlled IDs"),
                    )
                    for item in group_records
                ),
                str(method["parameterization_identity"]),
                str(method["evaluation_plan_identity"]),
                tuple(
                    PolicyRecord(str(item["name"]), str(item["identity"]))
                    for item in policy_records
                ),
                tuple(
                    BudgetRecord(
                        str(item["name"]),
                        _record_int(item["limit"], "budget limit"),
                        (
                            None
                            if "used" not in item
                            else _record_int(item["used"], "budget used")
                        ),
                    )
                    for item in budget_records
                ),
                tuple(
                    SeedRecord(
                        str(item["name"]),
                        _record_int(item["value"], "seed value"),
                        str(item["policy_identity"]),
                    )
                    for item in seed_records
                ),
                ExecutionSettings(
                    workers=_record_int(execution["workers"], "workers"),
                    native_threads=(
                        None
                        if execution["native_threads"] == "auto"
                        else _record_int(execution["native_threads"], "native threads")
                    ),
                ),
                ProvenanceEnvironment(
                    str(environment["chemex_version"]),
                    str(environment["python_version"]),
                    str(environment["python_implementation"]),
                    str(environment["platform"]),
                    str(environment["numpy_version"]),
                    str(environment["scipy_version"]),
                    str(environment["emcee_version"]),
                    tuple(
                        (str(item["kind"]), str(item["name"]), str(item["version"]))
                        for item in library_records
                    ),
                ),
                tuple(
                    BaselineReference(
                        _BASELINE_REFERENCE_KEY,
                        cast(
                            Literal["occurrence", "bundle", "manifest"],
                            str(item["kind"]),
                        ),
                        str(item["identity"]),
                        str(item["occurrence_identity"]),
                        str(item["result_bundle_identity"]),
                    )
                    for item in baseline_records
                ),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise NativeProvenanceError(
                "Native step workflow provenance is malformed"
            ) from error
        if (
            workflow.get("identity") != provenance.workflow_identity
            or workflow.get("execution_identity") != provenance.execution_identity
            or workflow.get("provenance_identity") != provenance.identity
            or method.get("identity") != provenance.method_identity
            or method.get("normalized_sha256") != provenance.normalized_method_sha256
            or environment.get("identity") != provenance.environment.identity
        ):
            raise NativeProvenanceError(
                "Native step workflow identities do not match their children"
            )
        return provenance

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.workflow_identity,
            "execution_identity": self.execution_identity,
            "workflow_type": self.workflow_type,
            "grouping_topology": [
                {"identity": name, "controlled_ids": list(ids)}
                for name, ids in self.grouping_topology
            ],
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
                {
                    "kind": item.kind,
                    "identity": item.identity,
                    "occurrence_identity": item.occurrence_identity,
                    "result_bundle_identity": item.result_bundle_identity,
                }
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
        if (
            Path(self.path).is_absolute()
            or ".." in Path(self.path).parts
            or Path(self.path).parts in {(), (".",)}
            or self.role
            not in {
                "committed_restart_state",
                "report_only_fitted_values",
                "diagnostic_provenance",
                "partial_evidence",
                "product_output",
            }
        ):
            raise NativeProvenanceError("Artifact link or role is invalid")
        digest = self.sha256.lower()
        if len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        ):
            raise NativeProvenanceError(
                "Artifact SHA-256 must be lowercase hexadecimal"
            )
        object.__setattr__(self, "sha256", digest)


@dataclass(frozen=True, slots=True)
class CommittedRestartRecord:
    """Durable exact binding from a restart file to one witnessed commit."""

    starting_occurrence_identity: str
    starting_revision: int
    accepted_result_identity: str
    accepted_occurrence_identity: str
    commit_operation_identity: str
    commit_occurrence_identity: str
    committed_occurrence_identity: str
    committed_revision: int
    workflow_identity: str
    problem_identity: str
    parameterization_identity: str
    parameter_model_identity: str
    configuration_identity: str
    independent_items: tuple[tuple[str, str, str, str], ...]
    restart_sha256: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        identities = (
            self.starting_occurrence_identity,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
            self.commit_operation_identity,
            self.commit_occurrence_identity,
            self.committed_occurrence_identity,
            self.workflow_identity,
            self.problem_identity,
            self.parameterization_identity,
            self.parameter_model_identity,
            self.configuration_identity,
        )
        if any(not value for value in identities):
            raise NativeProvenanceError("Committed restart lineage is incomplete")
        if (
            self.starting_revision < 0
            or self.committed_revision != self.starting_revision + 1
            or len(self.restart_sha256) != 64
            or any(
                character not in "0123456789abcdef" for character in self.restart_sha256
            )
        ):
            raise NativeProvenanceError("Committed restart revisions are inconsistent")
        items = tuple(self.independent_items)
        ids = tuple(item[0] for item in items)
        if not items or len(set(ids)) != len(ids) or tuple(sorted(ids)) != ids:
            raise NativeProvenanceError(
                "Committed restart independent frame is invalid"
            )
        for param_id, value, lower, upper in items:
            try:
                scalars = tuple(float.fromhex(item) for item in (value, lower, upper))
            except ValueError as error:
                raise NativeProvenanceError(
                    "Committed restart contains malformed binary64 values"
                ) from error
            if (
                not param_id
                or not math.isfinite(scalars[0])
                or math.isnan(scalars[1])
                or math.isnan(scalars[2])
                or tuple(item.hex() for item in scalars) != (value, lower, upper)
            ):
                raise NativeProvenanceError(
                    "Committed restart contains non-finite independent values"
                )
            if not scalars[1] <= scalars[0] <= scalars[2]:
                raise NativeProvenanceError(
                    "Committed restart value is outside configured bounds"
                )
        object.__setattr__(self, "independent_items", items)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-committed-restart-v2",
                (
                    *identities,
                    self.starting_revision,
                    self.committed_revision,
                    items,
                    self.restart_sha256,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": 2,
            "identity": self.identity,
            "starting_occurrence_identity": self.starting_occurrence_identity,
            "starting_revision": self.starting_revision,
            "accepted_result_identity": self.accepted_result_identity,
            "accepted_occurrence_identity": self.accepted_occurrence_identity,
            "commit_operation_identity": self.commit_operation_identity,
            "commit_occurrence_identity": self.commit_occurrence_identity,
            "committed_occurrence_identity": self.committed_occurrence_identity,
            "committed_revision": self.committed_revision,
            "workflow_identity": self.workflow_identity,
            "problem_identity": self.problem_identity,
            "parameterization_identity": self.parameterization_identity,
            "parameter_model_identity": self.parameter_model_identity,
            "configuration_identity": self.configuration_identity,
            "independent_items": [
                {"id": param_id, "value": value, "minimum": lower, "maximum": upper}
                for param_id, value, lower, upper in self.independent_items
            ],
            "restart_sha256": self.restart_sha256,
        }

    @classmethod
    def from_record(cls, record: object) -> CommittedRestartRecord:
        if not isinstance(record, dict):
            raise NativeProvenanceError("Committed restart provenance must be a record")
        record = {str(key): value for key, value in record.items()}
        keys = {
            "schema_version",
            "identity",
            "starting_occurrence_identity",
            "starting_revision",
            "accepted_result_identity",
            "accepted_occurrence_identity",
            "commit_operation_identity",
            "commit_occurrence_identity",
            "committed_occurrence_identity",
            "committed_revision",
            "workflow_identity",
            "problem_identity",
            "parameterization_identity",
            "parameter_model_identity",
            "configuration_identity",
            "independent_items",
            "restart_sha256",
        }
        if set(record) != keys or record.get("schema_version") != 2:
            raise NativeProvenanceError(
                "Committed restart provenance fields are invalid"
            )
        raw_items = record.get("independent_items")
        if not isinstance(raw_items, list) or any(
            not isinstance(item, dict)
            or set(item) != {"id", "value", "minimum", "maximum"}
            for item in raw_items
        ):
            raise NativeProvenanceError(
                "Committed restart independent frame is malformed"
            )
        item_records = cast(list[dict[str, object]], raw_items)
        rebuilt = cls(
            str(record["starting_occurrence_identity"]),
            _record_int(record["starting_revision"], "starting revision"),
            str(record["accepted_result_identity"]),
            str(record["accepted_occurrence_identity"]),
            str(record["commit_operation_identity"]),
            str(record["commit_occurrence_identity"]),
            str(record["committed_occurrence_identity"]),
            _record_int(record["committed_revision"], "committed revision"),
            str(record["workflow_identity"]),
            str(record["problem_identity"]),
            str(record["parameterization_identity"]),
            str(record["parameter_model_identity"]),
            str(record["configuration_identity"]),
            tuple(
                (
                    str(item["id"]),
                    str(item["value"]),
                    str(item["minimum"]),
                    str(item["maximum"]),
                )
                for item in item_records
            ),
            str(record["restart_sha256"]),
        )
        if record.get("identity") != rebuilt.identity:
            raise NativeProvenanceError(
                "Committed restart identity does not match children"
            )
        return rebuilt


def published_step_reference_identity(
    *,
    publication_occurrence_identity: str,
    lifecycle: str,
    authority: str,
    manifest_identity: str,
    manifest_sha256: str,
    provenance: WorkflowProvenance,
    independent_ids: tuple[str, ...],
    artifacts: tuple[ArtifactReference, ...],
) -> str:
    """Recompute a durable step-reference identity without minting live authority."""
    return _identity(
        "native-published-step-reference-v2",
        {
            "publication_occurrence_identity": publication_occurrence_identity,
            "lifecycle": lifecycle,
            "authority": authority,
            "manifest_identity": manifest_identity,
            "manifest_sha256": manifest_sha256,
            "workflow_provenance_identity": provenance.identity,
            "execution_identity": provenance.execution_identity,
            "independent_ids": independent_ids,
            "artifacts": tuple(
                (item.path, item.role, item.sha256) for item in artifacts
            ),
        },
    )


class _PublishedStepWitness:
    """Opaque process-local witness for one exact #608 publication occurrence."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> _PublishedStepWitness:
        raise TypeError("Published-step witnesses are minted only by publication")

    def __reduce__(self) -> tuple[object, ...]:
        raise TypeError("Published-step witnesses cannot be serialized")


@dataclass(frozen=True, slots=True)
class _PublishedStepBinding:
    occurrence_identity: str
    step_reference: ReferenceType[object] | None
    publication: ReferenceType[object]


_PUBLISHED_STEP_BINDINGS: WeakKeyDictionary[
    _PublishedStepWitness, _PublishedStepBinding
] = WeakKeyDictionary()
_PUBLISHED_STEP_LOCK = RLock()


def _mint_published_step_witness(
    occurrence_identity: str,
    publication: object,
) -> _PublishedStepWitness:
    witness = object.__new__(_PublishedStepWitness)
    with _PUBLISHED_STEP_LOCK:
        _PUBLISHED_STEP_BINDINGS[witness] = _PublishedStepBinding(
            occurrence_identity,
            None,
            ref(publication),
        )
    return witness


def _bind_published_step_reference(
    witness: _PublishedStepWitness,
    occurrence_identity: str,
    step: object,
) -> bool:
    with _PUBLISHED_STEP_LOCK:
        binding = _PUBLISHED_STEP_BINDINGS.get(witness)
        if binding is None or binding.occurrence_identity != occurrence_identity:
            return False
        current = None if binding.step_reference is None else binding.step_reference()
        if current is not None:
            return current is step
        _PUBLISHED_STEP_BINDINGS[witness] = _PublishedStepBinding(
            occurrence_identity,
            ref(step),
            binding.publication,
        )
        return True


def _is_exact_live_published_step(
    witness: _PublishedStepWitness | None,
    occurrence_identity: str,
    step: object,
) -> bool:
    if witness is None:
        return False
    with _PUBLISHED_STEP_LOCK:
        binding = _PUBLISHED_STEP_BINDINGS.get(witness)
        return (
            binding is not None
            and binding.occurrence_identity == occurrence_identity
            and binding.step_reference is not None
            and binding.step_reference() is step
            and binding.publication() is not None
        )


@dataclass(frozen=True, slots=True, weakref_slot=True)
class PublishedStepReference:
    """Verified link from invocation provenance to one atomic step publication."""

    _construction_key: object = field(repr=False, compare=False)
    publication_occurrence_identity: str
    path: Path
    lifecycle: StepLifecycle
    authority: StepAuthority
    manifest_identity: str
    manifest_sha256: str
    provenance: WorkflowProvenance
    independent_ids: tuple[str, ...]
    artifacts: tuple[ArtifactReference, ...]
    identity: str = field(init=False)
    _witness: _PublishedStepWitness | None = field(
        repr=False,
        compare=False,
        kw_only=True,
    )

    def __post_init__(self) -> None:
        if self._construction_key is not _PUBLISHED_STEP_KEY:
            raise NativeProvenanceError(
                "Published step references require an exact live publication"
            )
        _require_text(
            self.publication_occurrence_identity,
            "publication occurrence identity",
        )
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
        object.__setattr__(
            self,
            "identity",
            published_step_reference_identity(
                publication_occurrence_identity=self.publication_occurrence_identity,
                lifecycle=self.lifecycle,
                authority=self.authority,
                manifest_identity=self.manifest_identity,
                manifest_sha256=self.manifest_sha256,
                provenance=self.provenance,
                independent_ids=self.independent_ids,
                artifacts=self.artifacts,
            ),
        )
        if self._witness is None or not _bind_published_step_reference(
            self._witness,
            self.publication_occurrence_identity,
            self,
        ):
            raise NativeProvenanceError(
                "Published step reference lacks its exact live publication"
            )

    def require_exact_live_publication(self) -> None:
        """Reject reconstruction, replacement, and expired live publications."""
        if not _is_exact_live_published_step(
            self._witness,
            self.publication_occurrence_identity,
            self,
        ):
            raise NativeProvenanceError(
                "Schema-v2 writing requires the exact live publication"
            )


def _published_step_from_successful_native_publication(
    publication: object,
    publication_occurrence_identity: str,
    path: Path,
    lifecycle: StepLifecycle,
    authority: StepAuthority,
    manifest_identity: str,
    manifest_sha256: str,
    provenance: WorkflowProvenance,
    independent_ids: tuple[str, ...],
    artifacts: tuple[ArtifactReference, ...],
) -> PublishedStepReference:
    """Internal post-rename adapter from one successful #608 publication."""
    witness = _mint_published_step_witness(
        publication_occurrence_identity,
        publication,
    )
    return PublishedStepReference(
        _PUBLISHED_STEP_KEY,
        publication_occurrence_identity,
        path,
        lifecycle,
        authority,
        manifest_identity,
        manifest_sha256,
        provenance,
        independent_ids,
        artifacts,
        _witness=witness,
    )
