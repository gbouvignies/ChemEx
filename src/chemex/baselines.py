"""Small content-addressed baseline evidence harness for migration checkpoint 2.

The types in this module describe one approved ChemEx analysis without changing
the product runner.  They deliberately do not provide numerical lanes, gates,
calibration, optimizer execution, or a general artifact store.  The one runner
below invokes the existing lmfit-backed product path and freezes its output as
an observation, never as scientific truth.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import shutil
import tempfile
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from errno import EEXIST, ENOTEMPTY
from numbers import Real
from pathlib import Path
from typing import Literal, cast

from chemex.chemex import run
from chemex.cli import build_parser

_SCHEMA_VERSION = 1
_SEMANTIC_VERSION = "chemex-baseline-v1"
_CPMG_SOURCE_COMMIT = "d5ed0c87e8ce7a7f17745feea346af4dfbae7ecf"
_CPMG_LOCKFILE_HASH = "f0fb2ffc7b1a5ecd1bf7ac43956fc4861b96c058d158948b68b4e97027a6086a"
_HASH_LENGTH = 64
type OccurrenceLifecycle = Literal["REQUESTED", "SUCCEEDED", "FAILED"]
# #581 permits only these closed future comparison families.  #587 needs the
# request identity now, but this module deliberately executes no comparison.
type ComparisonFamily = Literal[
    "STRUCTURAL_EXACTNESS",
    "DETERMINISTIC_REPLAY",
    "NUMERICAL_ARTIFACT_COMPARISON",
    "STATISTICAL_COMPARISON",
    "CALIBRATION_SENSITIVITY",
    "OPERATIONAL_PERFORMANCE",
]


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {
            "kind": kind,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "record": record,
        },
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return _sha256(encoded)


def _canonical_float(value: float) -> str:
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError("Baseline semantic values must be finite binary64 values")
    return (0.0 if scalar == 0.0 else scalar).hex()


def _semantic_value(value: object) -> object:
    """Return a path-free, deterministic JSON-compatible semantic value."""
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, Real):
        return {"binary64": _canonical_float(float(value))}
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError("Baseline semantic mapping keys must be strings")
        if set(value) == {"binary64"}:
            raise ValueError("'binary64' is reserved for canonical float encoding")
        return {
            key: _semantic_value(item)
            for key, item in sorted(value.items(), key=lambda item: item[0])
        }
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray, str)):
        return [_semantic_value(item) for item in value]
    raise TypeError(f"Unsupported baseline semantic value: {type(value).__name__}")


def _canonical_json(value: object) -> str:
    return json.dumps(
        _semantic_value(value), ensure_ascii=True, separators=(",", ":"), sort_keys=True
    )


def _canonical_record_bytes(record: Mapping[str, object]) -> bytes:
    """Encode one already-validated baseline record exactly once."""
    return json.dumps(
        record,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("ascii")


def _reject_duplicate_object_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Baseline JSON records may not contain duplicate keys")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> object:
    raise ValueError(f"Baseline JSON does not permit {value}")


def _canonical_record_from_bytes(content: bytes, name: str) -> Mapping[str, object]:
    """Decode one canonical record without silently accepting equivalent bytes."""
    try:
        decoded = json.loads(
            content.decode("ascii"),
            object_pairs_hook=_reject_duplicate_object_keys,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"Malformed {name}") from error
    if not isinstance(decoded, Mapping):
        raise TypeError(f"{name} must be a record")
    record = cast("Mapping[str, object]", decoded)
    if _canonical_record_bytes(record) != content:
        raise ValueError(f"{name} is not canonically serialized")
    return record


def _canonical_strings(values: Sequence[str], name: str) -> tuple[str, ...]:
    result = tuple(values)
    if not result or any(not isinstance(value, str) or not value for value in result):
        raise ValueError(f"Baseline {name} must be non-empty strings")
    if tuple(sorted(result)) != result or len(set(result)) != len(result):
        raise ValueError(f"Baseline {name} must be unique and canonically ordered")
    return result


def _record_semantic_json(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"Baseline {field_name} must be canonical JSON")
    try:
        decoded = json.loads(value)
    except json.JSONDecodeError as error:
        raise ValueError(f"Baseline {field_name} is not JSON") from error
    canonical = json.dumps(
        _record_semantic_value(decoded),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    if canonical != value:
        raise ValueError(f"Baseline {field_name} is not canonical JSON")
    return value


def _record_semantic_value(value: object) -> object:
    """Validate a semantic value already encoded for a JSON record."""
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        raise TypeError("Baseline float values must use canonical binary64 encoding")
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError("Baseline semantic mapping keys must be strings")
        mapping = cast("Mapping[str, object]", value)
        if set(value) == {"binary64"}:
            encoded = mapping["binary64"]
            if not isinstance(encoded, str):
                raise TypeError("Baseline binary64 encoding must be a string")
            try:
                scalar = float.fromhex(encoded)
            except ValueError as error:
                raise ValueError("Baseline binary64 encoding is invalid") from error
            if _canonical_float(scalar) != encoded:
                raise ValueError("Baseline binary64 encoding is not canonical")
            return {"binary64": encoded}
        return {
            key: _record_semantic_value(item)
            for key, item in sorted(mapping.items(), key=lambda item: item[0])
        }
    if isinstance(value, list):
        return [_record_semantic_value(item) for item in value]
    raise TypeError("Baseline semantic value has an unsupported record type")


def _semantic_json_from_record(value: object, field_name: str) -> str:
    try:
        return json.dumps(
            _record_semantic_value(value),
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise type(error)(f"Malformed baseline {field_name}") from error


@dataclass(frozen=True, slots=True)
class CanonicalBaselineValue:
    """One immutable structured baseline value stored in canonical JSON form."""

    encoded: str

    @classmethod
    def from_value(cls, value: object) -> CanonicalBaselineValue:
        return cls(_canonical_json(value))

    @classmethod
    def from_record(cls, value: object, field_name: str) -> CanonicalBaselineValue:
        return cls(_semantic_json_from_record(value, field_name))

    def __post_init__(self) -> None:
        _record_semantic_json(self.encoded, "semantic value")

    def to_record_value(self) -> object:
        return json.loads(self.encoded)


def _record_exact_keys(
    record: Mapping[str, object], expected: set[str], name: str
) -> None:
    if set(record) != expected:
        raise ValueError(f"{name} record has unknown or missing fields")


def _record_hash(value: object, field_name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _HASH_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"Baseline {field_name} must be a SHA-256 digest")
    return value


def _record_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value:
        raise TypeError(f"Baseline {field_name} must be a non-empty string")
    return value


def _record_int(value: object, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise TypeError(f"Baseline {field_name} must be a non-negative integer")
    return value


def _member_identity(kind: str, role: str, content_hash: str, size: int) -> str:
    return _identity(kind, (role, content_hash, size))


def _member_record(
    role: str, content_hash: str, size: int, identity: str
) -> dict[str, object]:
    return {
        "role": role,
        "content_hash": content_hash,
        "size": size,
        "identity": identity,
    }


def _member_record_values(
    record: Mapping[str, object], name: str
) -> tuple[str, str, int]:
    _record_exact_keys(record, {"role", "content_hash", "size", "identity"}, name)
    return (
        _record_string(record.get("role"), "member role"),
        _record_hash(record.get("content_hash"), "member hash"),
        _record_int(record.get("size"), "member size"),
    )


@dataclass(frozen=True, slots=True)
class CaseSourceAuthority:
    """Frozen provenance of shipped scientific intent, not result authority."""

    source_commit: str
    lockfile_hash: str
    authority_version: str = "case-source-authority-v1"
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _record_string(self.source_commit, "case-source commit")
        _record_hash(self.lockfile_hash, "case-source lockfile hash")
        _record_string(self.authority_version, "case-source authority version")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "case-source-authority",
                (self.source_commit, self.lockfile_hash, self.authority_version),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "source_commit": self.source_commit,
            "lockfile_hash": self.lockfile_hash,
            "authority_version": self.authority_version,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> CaseSourceAuthority:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "source_commit",
                "lockfile_hash",
                "authority_version",
                "identity",
            },
            "Case-source authority",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported baseline case-source authority schema")
        authority = cls(
            _record_string(record.get("source_commit"), "case-source commit"),
            _record_hash(record.get("lockfile_hash"), "case-source lockfile hash"),
            _record_string(
                record.get("authority_version"), "case-source authority version"
            ),
        )
        if record.get("identity") != authority.identity:
            raise ValueError("Case-source authority identity does not match payload")
        return authority


@dataclass(frozen=True, slots=True, order=True)
class InputMember:
    """One logical scientific input role and immutable byte digest."""

    role: str
    content_hash: str
    size: int

    def __post_init__(self) -> None:
        _record_string(self.role, "input role")
        _record_hash(self.content_hash, "input member hash")
        _record_int(self.size, "input member size")

    @property
    def identity(self) -> str:
        return _member_identity(
            "case-input-member", self.role, self.content_hash, self.size
        )

    def to_record(self) -> dict[str, object]:
        return _member_record(self.role, self.content_hash, self.size, self.identity)

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> InputMember:
        member = cls(*_member_record_values(record, "Case input member"))
        if record.get("identity") != member.identity:
            raise ValueError("Case input-member identity does not match payload")
        return member


@dataclass(frozen=True, slots=True)
class CaseDefinition:
    """Immutable scientific intent independent of its source-tree location."""

    name: str
    source_authority: CaseSourceAuthority
    scientific_intent: str
    inputs: tuple[InputMember, ...]
    identity: str = field(init=False)

    @classmethod
    def create(
        cls,
        name: str,
        source_authority: CaseSourceAuthority,
        scientific_intent: Mapping[str, object],
        inputs: Sequence[InputMember],
    ) -> CaseDefinition:
        return cls(
            name,
            source_authority,
            _canonical_json(scientific_intent),
            tuple(sorted(inputs)),
        )

    def __post_init__(self) -> None:
        _record_string(self.name, "case name")
        _record_semantic_json(self.scientific_intent, "scientific intent")
        inputs = tuple(self.inputs)
        if (
            not inputs
            or tuple(sorted(inputs)) != inputs
            or len({member.role for member in inputs}) != len(inputs)
        ):
            raise ValueError(
                "Case definition inputs must be unique and canonically role-ordered"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "case-definition",
                (
                    self.name,
                    self.source_authority.identity,
                    self.scientific_intent,
                    tuple(member.identity for member in inputs),
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "name": self.name,
            "source_authority": self.source_authority.to_record(),
            "scientific_intent": json.loads(self.scientific_intent),
            "inputs": [member.to_record() for member in self.inputs],
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> CaseDefinition:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "name",
                "source_authority",
                "scientific_intent",
                "inputs",
                "identity",
            },
            "Case definition",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported baseline case-definition semantics")
        source = record.get("source_authority")
        raw_inputs = record.get("inputs")
        if not isinstance(source, Mapping) or not isinstance(raw_inputs, list):
            raise TypeError("Malformed baseline case-definition record")
        case = cls(
            _record_string(record.get("name"), "case name"),
            CaseSourceAuthority.from_record(cast("Mapping[str, object]", source)),
            _semantic_json_from_record(
                record.get("scientific_intent"), "scientific intent"
            ),
            tuple(
                InputMember.from_record(cast("Mapping[str, object]", item))
                for item in raw_inputs
                if isinstance(item, Mapping)
            ),
        )
        if len(case.inputs) != len(raw_inputs):
            raise TypeError("Case definition inputs must be records")
        if record.get("identity") != case.identity:
            raise ValueError("Case definition identity does not match payload")
        return case


@dataclass(frozen=True, slots=True)
class LegacyObservationImplementation:
    """The only implementation authority used by this checkpoint's anchor."""

    implementation_version: str = "legacy-lmfit-product-v1"
    package_version: str = "test-package"
    lmfit_version: str = "test-lmfit"
    source_manifest_hash: str = "0" * _HASH_LENGTH
    authority_role: str = "LegacyObservationImplementation"
    identity: str = field(init=False)

    @classmethod
    def from_current_package(cls) -> LegacyObservationImplementation:
        package_root = Path(__file__).parent
        members = tuple(
            sorted(
                (
                    path.relative_to(package_root).as_posix(),
                    _sha256(path.read_bytes()),
                )
                for path in package_root.rglob("*.py")
            )
        )
        return cls(
            package_version=importlib.metadata.version("chemex"),
            lmfit_version=importlib.metadata.version("lmfit"),
            source_manifest_hash=_identity("legacy-package-source", members),
        )

    def __post_init__(self) -> None:
        if self.authority_role != "LegacyObservationImplementation":
            raise ValueError("Baseline legacy implementation authority role is fixed")
        _record_string(self.implementation_version, "legacy implementation version")
        _record_string(self.package_version, "legacy package version")
        _record_string(self.lmfit_version, "legacy lmfit version")
        _record_hash(self.source_manifest_hash, "legacy source manifest hash")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "implementation",
                (
                    self.authority_role,
                    self.implementation_version,
                    self.package_version,
                    self.lmfit_version,
                    self.source_manifest_hash,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "authority_role": self.authority_role,
            "implementation_version": self.implementation_version,
            "package_version": self.package_version,
            "lmfit_version": self.lmfit_version,
            "source_manifest_hash": self.source_manifest_hash,
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls, record: Mapping[str, object]
    ) -> LegacyObservationImplementation:
        _record_exact_keys(
            record,
            {
                "authority_role",
                "implementation_version",
                "package_version",
                "lmfit_version",
                "source_manifest_hash",
                "identity",
            },
            "Legacy observation implementation",
        )
        implementation = cls(
            _record_string(
                record.get("implementation_version"), "legacy implementation version"
            ),
            _record_string(record.get("package_version"), "legacy package version"),
            _record_string(record.get("lmfit_version"), "legacy lmfit version"),
            _record_hash(
                record.get("source_manifest_hash"), "legacy source manifest hash"
            ),
            _record_string(record.get("authority_role"), "legacy authority role"),
        )
        if record.get("identity") != implementation.identity:
            raise ValueError("Legacy implementation identity does not match payload")
        return implementation


@dataclass(frozen=True, slots=True)
class ExecutionSpecification:
    """Frozen execution request, intentionally separate from scientific intent."""

    case_identity: str
    implementation: LegacyObservationImplementation
    workflow: CanonicalBaselineValue
    lane_reference: str
    policy: CanonicalBaselineValue
    budget: CanonicalBaselineValue
    seed: CanonicalBaselineValue
    execution_settings: CanonicalBaselineValue
    artifact_inventory: CanonicalBaselineValue
    roles: tuple[str, ...]
    claims: tuple[str, ...]
    identity: str = field(init=False)

    @classmethod
    def create(
        cls,
        case: CaseDefinition,
        implementation: LegacyObservationImplementation,
        *,
        workflow: Mapping[str, object],
        lane_reference: str,
        policy: Mapping[str, object],
        budget: Mapping[str, object],
        seed: object,
        execution_settings: Mapping[str, object],
        artifact_inventory: Mapping[str, object],
        roles: Sequence[str],
        claims: Sequence[str],
    ) -> ExecutionSpecification:
        return cls(
            case.identity,
            implementation,
            CanonicalBaselineValue.from_value(workflow),
            lane_reference,
            CanonicalBaselineValue.from_value(policy),
            CanonicalBaselineValue.from_value(budget),
            CanonicalBaselineValue.from_value(seed),
            CanonicalBaselineValue.from_value(execution_settings),
            CanonicalBaselineValue.from_value(artifact_inventory),
            _canonical_strings(roles, "execution roles"),
            _canonical_strings(claims, "execution claims"),
        )

    def __post_init__(self) -> None:
        _record_hash(self.case_identity, "execution case identity")
        _record_string(self.lane_reference, "lane reference")
        for value, name in (
            (self.workflow, "workflow"),
            (self.policy, "policy"),
            (self.budget, "budget"),
            (self.seed, "seed"),
            (self.execution_settings, "execution settings"),
            (self.artifact_inventory, "artifact inventory"),
        ):
            if not isinstance(value, CanonicalBaselineValue):
                raise TypeError(f"Execution specification {name} must be canonical")
        _canonical_strings(self.roles, "execution roles")
        _canonical_strings(self.claims, "execution claims")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "execution-specification",
                (
                    self.case_identity,
                    self.implementation.identity,
                    self.workflow.encoded,
                    self.lane_reference,
                    self.policy.encoded,
                    self.budget.encoded,
                    self.seed.encoded,
                    self.execution_settings.encoded,
                    self.artifact_inventory.encoded,
                    self.roles,
                    self.claims,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "case_identity": self.case_identity,
            "implementation": self.implementation.to_record(),
            "workflow": self.workflow.to_record_value(),
            "lane_reference": self.lane_reference,
            "policy": self.policy.to_record_value(),
            "budget": self.budget.to_record_value(),
            "seed": self.seed.to_record_value(),
            "execution_settings": self.execution_settings.to_record_value(),
            "artifact_inventory": self.artifact_inventory.to_record_value(),
            "roles": list(self.roles),
            "claims": list(self.claims),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ExecutionSpecification:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "case_identity",
                "implementation",
                "workflow",
                "lane_reference",
                "policy",
                "budget",
                "seed",
                "execution_settings",
                "artifact_inventory",
                "roles",
                "claims",
                "identity",
            },
            "Execution specification",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported baseline execution-specification semantics")
        raw_implementation = record.get("implementation")
        raw_roles = record.get("roles")
        raw_claims = record.get("claims")
        if (
            not isinstance(raw_implementation, Mapping)
            or not isinstance(raw_roles, list)
            or not isinstance(raw_claims, list)
        ):
            raise TypeError("Execution specification implementation must be a record")
        specification = cls(
            _record_hash(record.get("case_identity"), "execution case identity"),
            LegacyObservationImplementation.from_record(
                cast("Mapping[str, object]", raw_implementation)
            ),
            CanonicalBaselineValue.from_record(record.get("workflow"), "workflow"),
            _record_string(record.get("lane_reference"), "lane reference"),
            CanonicalBaselineValue.from_record(record.get("policy"), "policy"),
            CanonicalBaselineValue.from_record(record.get("budget"), "budget"),
            CanonicalBaselineValue.from_record(record.get("seed"), "seed"),
            CanonicalBaselineValue.from_record(
                record.get("execution_settings"), "execution settings"
            ),
            CanonicalBaselineValue.from_record(
                record.get("artifact_inventory"), "artifact inventory"
            ),
            _canonical_strings(cast("Sequence[str]", raw_roles), "execution roles"),
            _canonical_strings(cast("Sequence[str]", raw_claims), "execution claims"),
        )
        if record.get("identity") != specification.identity:
            raise ValueError("Execution specification identity does not match payload")
        return specification


@dataclass(frozen=True, slots=True)
class Occurrence:
    """One immutable attempt identity and its closed lifecycle record."""

    execution_specification_identity: str
    case_identity: str
    actual_implementation_identity: str
    lane_reference: str
    input_member_identities: tuple[str, ...]
    attempt_token: str
    lifecycle: OccurrenceLifecycle = "REQUESTED"
    result_bundle_identity: str | None = None
    failure_code: str | None = None

    @classmethod
    def requested(
        cls,
        specification: ExecutionSpecification,
        case: CaseDefinition,
        attempt_token: str,
    ) -> Occurrence:
        if specification.case_identity != case.identity:
            raise ValueError("Occurrence case does not match execution specification")
        return cls(
            specification.identity,
            case.identity,
            specification.implementation.identity,
            specification.lane_reference,
            tuple(sorted(member.identity for member in case.inputs)),
            attempt_token,
        )

    def __post_init__(self) -> None:
        _record_hash(self.execution_specification_identity, "occurrence specification")
        _record_hash(self.case_identity, "occurrence case")
        _record_hash(self.actual_implementation_identity, "occurrence implementation")
        _record_string(self.lane_reference, "occurrence lane reference")
        if (
            not self.input_member_identities
            or tuple(sorted(self.input_member_identities))
            != self.input_member_identities
        ):
            raise ValueError(
                "Occurrence input identities must be non-empty and ordered"
            )
        for identity in self.input_member_identities:
            _record_hash(identity, "occurrence input identity")
        _record_string(self.attempt_token, "occurrence attempt token")
        if self.lifecycle == "REQUESTED":
            if self.result_bundle_identity is not None or self.failure_code is not None:
                raise ValueError("Requested occurrence cannot have a terminal result")
        elif self.lifecycle == "SUCCEEDED":
            if self.result_bundle_identity is None or self.failure_code is not None:
                raise ValueError(
                    "Successful occurrence requires exactly one result bundle"
                )
            _record_hash(self.result_bundle_identity, "occurrence result bundle")
        elif self.lifecycle == "FAILED":
            if self.result_bundle_identity is not None or self.failure_code is None:
                raise ValueError("Failed occurrence cannot have a result bundle")
            _record_string(self.failure_code, "occurrence failure code")
        else:
            raise ValueError("Unknown baseline occurrence lifecycle")

    @property
    def identity(self) -> str:
        return _identity(
            "occurrence",
            (
                self.execution_specification_identity,
                self.case_identity,
                self.actual_implementation_identity,
                self.lane_reference,
                self.input_member_identities,
                self.attempt_token,
            ),
        )

    @property
    def lifecycle_identity(self) -> str:
        return _identity(
            "occurrence-lifecycle",
            (
                self.identity,
                self.lifecycle,
                self.result_bundle_identity,
                self.failure_code,
            ),
        )

    def succeeded(self, bundle: ResultBundle) -> Occurrence:
        if self.lifecycle != "REQUESTED" or bundle.occurrence_identity != self.identity:
            raise ValueError("Only a matching requested occurrence may succeed")
        return Occurrence(
            self.execution_specification_identity,
            self.case_identity,
            self.actual_implementation_identity,
            self.lane_reference,
            self.input_member_identities,
            self.attempt_token,
            "SUCCEEDED",
            bundle.identity,
        )

    def failed(self, failure_code: str) -> Occurrence:
        if self.lifecycle != "REQUESTED":
            raise ValueError("Only a requested occurrence may fail")
        return Occurrence(
            self.execution_specification_identity,
            self.case_identity,
            self.actual_implementation_identity,
            self.lane_reference,
            self.input_member_identities,
            self.attempt_token,
            "FAILED",
            failure_code=failure_code,
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "execution_specification_identity": self.execution_specification_identity,
            "case_identity": self.case_identity,
            "actual_implementation_identity": self.actual_implementation_identity,
            "lane_reference": self.lane_reference,
            "input_member_identities": list(self.input_member_identities),
            "attempt_token": self.attempt_token,
            "lifecycle": self.lifecycle,
            "result_bundle_identity": self.result_bundle_identity,
            "failure_code": self.failure_code,
            "identity": self.identity,
            "lifecycle_identity": self.lifecycle_identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        specification: ExecutionSpecification | None = None,
    ) -> Occurrence:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "execution_specification_identity",
                "case_identity",
                "actual_implementation_identity",
                "lane_reference",
                "input_member_identities",
                "attempt_token",
                "lifecycle",
                "result_bundle_identity",
                "failure_code",
                "identity",
                "lifecycle_identity",
            },
            "Occurrence",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported baseline occurrence schema")
        lifecycle = record.get("lifecycle")
        if lifecycle not in {"REQUESTED", "SUCCEEDED", "FAILED"}:
            raise ValueError("Unknown baseline occurrence lifecycle")
        raw_inputs = record.get("input_member_identities")
        if not isinstance(raw_inputs, list):
            raise TypeError("Occurrence input identities must be a list")
        occurrence = cls(
            _record_hash(
                record.get("execution_specification_identity"),
                "occurrence specification",
            ),
            _record_hash(record.get("case_identity"), "occurrence case"),
            _record_hash(
                record.get("actual_implementation_identity"),
                "occurrence implementation",
            ),
            _record_string(record.get("lane_reference"), "occurrence lane reference"),
            tuple(
                _record_hash(value, "occurrence input identity") for value in raw_inputs
            ),
            _record_string(record.get("attempt_token"), "occurrence attempt token"),
            cast("OccurrenceLifecycle", lifecycle),
            cast("str | None", record.get("result_bundle_identity")),
            cast("str | None", record.get("failure_code")),
        )
        if specification is not None and (
            occurrence.execution_specification_identity != specification.identity
            or occurrence.actual_implementation_identity
            != specification.implementation.identity
            or occurrence.lane_reference != specification.lane_reference
        ):
            raise ValueError(
                "Occurrence does not belong to the execution specification"
            )
        if (
            record.get("identity") != occurrence.identity
            or record.get("lifecycle_identity") != occurrence.lifecycle_identity
        ):
            raise ValueError("Occurrence identities do not match payload")
        return occurrence


@dataclass(frozen=True, slots=True, order=True)
class ResultMember:
    """One logical result artifact with a content-derived member identity."""

    role: str
    content_hash: str
    size: int

    def __post_init__(self) -> None:
        _record_string(self.role, "result member role")
        _record_hash(self.content_hash, "result member hash")
        _record_int(self.size, "result member size")

    @property
    def identity(self) -> str:
        return _member_identity(
            "result-member", self.role, self.content_hash, self.size
        )

    def to_record(self) -> dict[str, object]:
        return _member_record(self.role, self.content_hash, self.size, self.identity)

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ResultMember:
        member = cls(*_member_record_values(record, "Result member"))
        if record.get("identity") != member.identity:
            raise ValueError("Result member identity does not match payload")
        return member


@dataclass(frozen=True, slots=True)
class ResultBundle:
    """Immutable observed evidence with an explicit manifest identity."""

    occurrence_identity: str
    execution_specification_identity: str
    implementation: LegacyObservationImplementation
    members: tuple[ResultMember, ...]
    ordering_version: str = "result-member-role-lexicographic-v1"
    manifest_identity: str = field(init=False)
    identity: str = field(init=False)

    @classmethod
    def create(
        cls,
        occurrence_identity: str,
        execution_specification_identity: str,
        implementation: LegacyObservationImplementation,
        members: Sequence[ResultMember],
    ) -> ResultBundle:
        return cls(
            occurrence_identity,
            execution_specification_identity,
            implementation,
            tuple(sorted(members)),
        )

    def __post_init__(self) -> None:
        _record_hash(self.occurrence_identity, "bundle occurrence")
        _record_hash(self.execution_specification_identity, "bundle specification")
        _record_string(self.ordering_version, "bundle ordering version")
        members = tuple(self.members)
        if (
            not members
            or tuple(sorted(members)) != members
            or len({member.role for member in members}) != len(members)
        ):
            raise ValueError(
                "Result bundle members must be unique and canonically role-ordered"
            )
        manifest_identity = _identity(
            "result-manifest",
            (
                self.occurrence_identity,
                self.execution_specification_identity,
                self.implementation.identity,
                self.implementation.authority_role,
                self.ordering_version,
                tuple(member.identity for member in members),
            ),
        )
        object.__setattr__(self, "manifest_identity", manifest_identity)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "result-bundle",
                (manifest_identity, tuple(member.identity for member in members)),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "occurrence_identity": self.occurrence_identity,
            "execution_specification_identity": self.execution_specification_identity,
            "implementation": self.implementation.to_record(),
            "ordering_version": self.ordering_version,
            "members": [member.to_record() for member in self.members],
            "manifest_identity": self.manifest_identity,
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        occurrence: Occurrence | None = None,
    ) -> ResultBundle:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "occurrence_identity",
                "execution_specification_identity",
                "implementation",
                "ordering_version",
                "members",
                "manifest_identity",
                "identity",
            },
            "Result bundle",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported baseline result-bundle semantics")
        raw_implementation = record.get("implementation")
        raw_members = record.get("members")
        if not isinstance(raw_implementation, Mapping) or not isinstance(
            raw_members, list
        ):
            raise TypeError("Malformed baseline result-bundle record")
        bundle = cls(
            _record_hash(record.get("occurrence_identity"), "bundle occurrence"),
            _record_hash(
                record.get("execution_specification_identity"), "bundle specification"
            ),
            LegacyObservationImplementation.from_record(
                cast("Mapping[str, object]", raw_implementation)
            ),
            tuple(
                ResultMember.from_record(cast("Mapping[str, object]", item))
                for item in raw_members
                if isinstance(item, Mapping)
            ),
            _record_string(record.get("ordering_version"), "bundle ordering version"),
        )
        if len(bundle.members) != len(raw_members):
            raise TypeError("Result bundle members must be records")
        if occurrence is not None:
            if occurrence.lifecycle != "SUCCEEDED":
                raise ValueError("Result bundle requires a successful occurrence")
            if (
                bundle.occurrence_identity != occurrence.identity
                or bundle.execution_specification_identity
                != occurrence.execution_specification_identity
                or bundle.implementation.identity
                != occurrence.actual_implementation_identity
                or occurrence.result_bundle_identity != bundle.identity
            ):
                raise ValueError("Result bundle does not belong to occurrence")
        if (
            record.get("manifest_identity") != bundle.manifest_identity
            or record.get("identity") != bundle.identity
        ):
            raise ValueError("Result bundle identities do not match payload")
        return bundle


@dataclass(frozen=True, slots=True)
class BaselineComparisonRequest:
    """#587's identity-only request for one #581-permitted future comparison.

    It is neither a comparison execution nor a gate, and cannot select lanes,
    thresholds, calibration, or a result authority.
    """

    observation_bundle_identity: str
    subject_bundle_identity: str
    comparison_family: ComparisonFamily
    policy: str
    identity: str = field(init=False)

    @classmethod
    def create(
        cls,
        observation: ResultBundle,
        subject: ResultBundle,
        comparison_family: ComparisonFamily,
        policy: Mapping[str, object],
    ) -> BaselineComparisonRequest:
        return cls(
            observation.identity,
            subject.identity,
            comparison_family,
            _canonical_json(policy),
        )

    def __post_init__(self) -> None:
        _record_hash(self.observation_bundle_identity, "comparison observation bundle")
        _record_hash(self.subject_bundle_identity, "comparison subject bundle")
        if self.observation_bundle_identity == self.subject_bundle_identity:
            raise ValueError("Baseline comparison needs two distinct bundles")
        if self.comparison_family not in {
            "STRUCTURAL_EXACTNESS",
            "DETERMINISTIC_REPLAY",
            "NUMERICAL_ARTIFACT_COMPARISON",
            "STATISTICAL_COMPARISON",
            "CALIBRATION_SENSITIVITY",
            "OPERATIONAL_PERFORMANCE",
        }:
            raise ValueError("Unknown baseline comparison family")
        _record_semantic_json(self.policy, "comparison policy")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "baseline-comparison-request",
                (
                    self.observation_bundle_identity,
                    self.subject_bundle_identity,
                    self.comparison_family,
                    self.policy,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "observation_bundle_identity": self.observation_bundle_identity,
            "subject_bundle_identity": self.subject_bundle_identity,
            "comparison_family": self.comparison_family,
            "policy": json.loads(self.policy),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> BaselineComparisonRequest:
        _record_exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "observation_bundle_identity",
                "subject_bundle_identity",
                "comparison_family",
                "policy",
                "identity",
            },
            "Baseline comparison request",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported baseline comparison-request semantics")
        family = _record_string(record.get("comparison_family"), "comparison family")
        request = cls(
            _record_hash(
                record.get("observation_bundle_identity"),
                "comparison observation bundle",
            ),
            _record_hash(
                record.get("subject_bundle_identity"), "comparison subject bundle"
            ),
            cast("ComparisonFamily", family),
            _semantic_json_from_record(record.get("policy"), "comparison policy"),
        )
        if record.get("identity") != request.identity:
            raise ValueError(
                "Baseline comparison request identity does not match payload"
            )
        return request


@dataclass(frozen=True, slots=True)
class ArtifactContent:
    """Bytes collected before publication; never an identity based on a path."""

    role: str
    content: bytes

    @property
    def member(self) -> ResultMember:
        return ResultMember(self.role, _sha256(self.content), len(self.content))


@dataclass(frozen=True, slots=True)
class PublishedEvidence:
    case: CaseDefinition
    specification: ExecutionSpecification
    occurrence: Occurrence
    bundle: ResultBundle
    manifest_identity: str
    location: Path


class _TerminalInstallationCollision(FileExistsError):
    """The atomic terminal destination already exists for this occurrence."""


class BaselinePublicationIntegrityError(RuntimeError):
    """A terminal collision did not leave readable authoritative evidence."""


class CpmgBaselinePublisher:
    """Reservation-backed, append-only publisher for this one CPMG anchor."""

    def __init__(self, root: Path) -> None:
        self.root = root

    @property
    def _attempts(self) -> Path:
        return self.root / "attempts"

    @staticmethod
    def _required_roles(specification: ExecutionSpecification) -> tuple[str, ...]:
        inventory = specification.artifact_inventory.to_record_value()
        if not isinstance(inventory, Mapping):
            raise TypeError("CPMG artifact inventory must be a record")
        if set(inventory) != {
            "closed",
            "excluded_path_components",
            "required_roles",
            "structured_suffixes",
            "version",
        }:
            raise ValueError("CPMG artifact inventory has unknown or missing fields")
        if inventory.get("version") != "cpmg-structured-artifact-contract-v1":
            raise ValueError("CPMG artifact inventory has an unsupported version")
        if inventory.get("closed") is not True:
            raise ValueError("CPMG artifact inventory must be closed")
        raw_roles = inventory.get("required_roles")
        raw_suffixes = inventory.get("structured_suffixes")
        raw_excluded = inventory.get("excluded_path_components")
        if (
            not isinstance(raw_roles, list)
            or not isinstance(raw_suffixes, list)
            or not isinstance(raw_excluded, list)
            or not all(isinstance(value, str) for value in raw_suffixes)
            or not all(isinstance(value, str) for value in raw_excluded)
        ):
            raise TypeError("CPMG artifact inventory has malformed fields")
        return _canonical_strings(cast("Sequence[str]", raw_roles), "artifact roles")

    @staticmethod
    def _validate_lineage(
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
    ) -> None:
        if (
            specification.case_identity != case.identity
            or occurrence.execution_specification_identity != specification.identity
            or occurrence.case_identity != case.identity
            or occurrence.actual_implementation_identity
            != specification.implementation.identity
            or occurrence.lane_reference != specification.lane_reference
            or occurrence.input_member_identities
            != tuple(sorted(member.identity for member in case.inputs))
        ):
            raise ValueError("Published occurrence has incompatible case lineage")

    def _attempt_directory(self, occurrence: Occurrence) -> Path:
        return self._attempts / occurrence.identity

    def _validate_reservation(
        self,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
    ) -> Path:
        self._validate_lineage(case, specification, occurrence)
        if occurrence.lifecycle != "REQUESTED":
            raise ValueError("Only a requested occurrence can be reserved")
        attempt_directory = self._attempt_directory(occurrence)
        try:
            raw = (attempt_directory / "attempt.json").read_bytes()
        except OSError as error:
            raise ValueError("Baseline occurrence was not reserved") from error
        persisted = Occurrence.from_record(
            _canonical_record_from_bytes(raw, "Baseline attempt record"), specification
        )
        if persisted != occurrence:
            raise ValueError("Baseline attempt reservation does not match occurrence")
        return attempt_directory

    @staticmethod
    def _requested_occurrence(occurrence: Occurrence) -> Occurrence:
        return Occurrence(
            occurrence.execution_specification_identity,
            occurrence.case_identity,
            occurrence.actual_implementation_identity,
            occurrence.lane_reference,
            occurrence.input_member_identities,
            occurrence.attempt_token,
        )

    def terminal_occurrence(
        self,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        requested: Occurrence,
    ) -> Occurrence:
        """Read the single persisted lifecycle outcome for one reserved attempt."""
        attempt_directory = self._validate_reservation(case, specification, requested)
        terminal_directory = attempt_directory / "terminal"
        failure_path = terminal_directory / "failure.json"
        success_path = terminal_directory / "success" / "manifest.json"
        try:
            if failure_path.is_file() and not success_path.exists():
                terminal = Occurrence.from_record(
                    _canonical_record_from_bytes(
                        failure_path.read_bytes(), "Baseline failure record"
                    ),
                    specification,
                )
            elif success_path.is_file() and not failure_path.exists():
                manifest = _canonical_record_from_bytes(
                    success_path.read_bytes(), "Published baseline manifest"
                )
                raw_occurrence = manifest.get("occurrence")
                if not isinstance(raw_occurrence, Mapping):
                    raise TypeError("Published baseline manifest has no occurrence")
                terminal = Occurrence.from_record(
                    cast("Mapping[str, object]", raw_occurrence), specification
                )
            else:
                raise ValueError(
                    "Baseline occurrence has no unambiguous terminal record"
                )
        except OSError as error:
            raise ValueError("Baseline occurrence has no terminal record") from error
        if terminal.lifecycle not in {"SUCCEEDED", "FAILED"}:
            raise ValueError("Baseline terminal record is not terminal")
        if self._requested_occurrence(terminal) != requested:
            raise ValueError("Baseline terminal record does not match its reservation")
        return terminal

    def reserve(
        self,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
    ) -> Occurrence:
        """Persist the immutable requested-attempt record before execution."""
        self._validate_lineage(case, specification, occurrence)
        if occurrence.lifecycle != "REQUESTED":
            raise ValueError("Only requested occurrences may be reserved")
        self._attempts.mkdir(parents=True, exist_ok=True)
        destination = self._attempt_directory(occurrence)
        staging = Path(tempfile.mkdtemp(prefix=".attempt-", dir=self._attempts))
        try:
            (staging / "attempt.json").write_bytes(
                _canonical_record_bytes(occurrence.to_record())
            )
            try:
                os.replace(staging, destination)
            except FileExistsError as error:
                raise FileExistsError(
                    "Baseline occurrence is already reserved"
                ) from error
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        return occurrence

    @staticmethod
    def _publish_terminal(staging: Path, attempt_directory: Path) -> None:
        """Atomically make exactly one complete terminal state visible."""
        try:
            os.replace(staging, attempt_directory / "terminal")
        except OSError as error:
            if error.errno not in {EEXIST, ENOTEMPTY}:
                raise
            raise _TerminalInstallationCollision(
                "Baseline occurrence already has a terminal outcome"
            ) from error

    def _install_terminal(
        self,
        staging: Path,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        requested: Occurrence,
        attempt_directory: Path,
    ) -> None:
        """Install one terminal state or expose its persisted authoritative winner."""
        try:
            self._publish_terminal(staging, attempt_directory)
        except _TerminalInstallationCollision as collision:
            try:
                winner = self.terminal_occurrence(case, specification, requested)
            except (OSError, TypeError, ValueError) as error:
                raise BaselinePublicationIntegrityError(
                    "Terminal installation conflicted without valid persisted evidence"
                ) from error
            raise BaselineLifecycleConflictError(winner) from collision

    def record_failure(
        self,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        requested: Occurrence,
        failure_code: str,
    ) -> Occurrence:
        """Close one reserved occurrence as failed without publishing evidence."""
        attempt_directory = self._validate_reservation(case, specification, requested)
        failed = requested.failed(failure_code)
        staging = Path(tempfile.mkdtemp(prefix=".failure-", dir=attempt_directory))
        try:
            (staging / "failure.json").write_bytes(
                _canonical_record_bytes(failed.to_record())
            )
            self._install_terminal(
                staging, case, specification, requested, attempt_directory
            )
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        finally:
            shutil.rmtree(staging, ignore_errors=True)
        return failed

    def publish(
        self,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
        bundle: ResultBundle,
        artifacts: Sequence[ArtifactContent],
    ) -> PublishedEvidence:
        requested = self._requested_occurrence(occurrence)
        attempt_directory = self._validate_reservation(case, specification, requested)
        if (
            occurrence.lifecycle != "SUCCEEDED"
            or occurrence.result_bundle_identity != bundle.identity
        ):
            raise ValueError("Only a completed occurrence can publish a result bundle")
        if (
            bundle.occurrence_identity != occurrence.identity
            or bundle.execution_specification_identity != specification.identity
            or bundle.implementation.identity != specification.implementation.identity
        ):
            raise ValueError("Published bundle does not belong to occurrence")
        ordered_artifacts = tuple(sorted(artifacts, key=lambda item: item.role))
        if tuple(item.member for item in ordered_artifacts) != bundle.members:
            raise ValueError("Published artifacts do not match the result bundle")
        if tuple(item.role for item in ordered_artifacts) != self._required_roles(
            specification
        ):
            raise ValueError("CPMG result artifacts do not match the frozen contract")
        destination = attempt_directory / "terminal"
        staging = Path(tempfile.mkdtemp(prefix=".success-", dir=attempt_directory))
        try:
            success = staging / "success"
            success.mkdir()
            members_path = success / "members"
            members_path.mkdir()
            for artifact in ordered_artifacts:
                (members_path / artifact.member.content_hash).write_bytes(
                    artifact.content
                )
            manifest_without_identity: dict[str, object] = {
                "schema_version": _SCHEMA_VERSION,
                "semantic_version": _SEMANTIC_VERSION,
                "case": case.to_record(),
                "specification": specification.to_record(),
                "occurrence": occurrence.to_record(),
                "bundle": bundle.to_record(),
                "member_manifest_identity": bundle.manifest_identity,
            }
            manifest_identity = self._publication_manifest_identity(
                case, specification, occurrence, bundle
            )
            manifest = {
                **manifest_without_identity,
                "manifest_identity": manifest_identity,
            }
            (success / "manifest.json").write_bytes(_canonical_record_bytes(manifest))
            self._install_terminal(
                staging, case, specification, requested, attempt_directory
            )
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        return PublishedEvidence(
            case,
            specification,
            occurrence,
            bundle,
            manifest_identity,
            destination / "success",
        )

    def _find_success(self, identity: str) -> tuple[Path, Mapping[str, object]]:
        locations = (
            tuple(self._attempts.glob("*/terminal/success"))
            if self._attempts.exists()
            else ()
        )
        for candidate in locations:
            try:
                candidate_manifest = _canonical_record_from_bytes(
                    (candidate / "manifest.json").read_bytes(),
                    "Published baseline manifest",
                )
            except OSError:
                continue
            raw_bundle = candidate_manifest.get("bundle")
            if (
                isinstance(raw_bundle, Mapping)
                and raw_bundle.get("identity") == identity
            ):
                return candidate, candidate_manifest
        raise ValueError("Published baseline evidence is missing")

    @staticmethod
    def _manifest_records(
        manifest: Mapping[str, object],
    ) -> tuple[CaseDefinition, ExecutionSpecification, Occurrence, ResultBundle]:
        _record_exact_keys(
            manifest,
            {
                "schema_version",
                "semantic_version",
                "case",
                "specification",
                "occurrence",
                "bundle",
                "member_manifest_identity",
                "manifest_identity",
            },
            "Published baseline",
        )
        if (
            manifest.get("schema_version") != _SCHEMA_VERSION
            or manifest.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported published baseline manifest semantics")
        raw_case = manifest.get("case")
        raw_specification = manifest.get("specification")
        raw_occurrence = manifest.get("occurrence")
        raw_bundle = manifest.get("bundle")
        if (
            not isinstance(raw_case, Mapping)
            or not isinstance(raw_specification, Mapping)
            or not isinstance(raw_occurrence, Mapping)
            or not isinstance(raw_bundle, Mapping)
        ):
            raise TypeError("Published baseline manifest has malformed records")
        case = CaseDefinition.from_record(cast("Mapping[str, object]", raw_case))
        specification = ExecutionSpecification.from_record(
            cast("Mapping[str, object]", raw_specification)
        )
        occurrence = Occurrence.from_record(
            cast("Mapping[str, object]", raw_occurrence), specification
        )
        bundle = ResultBundle.from_record(
            cast("Mapping[str, object]", raw_bundle), occurrence
        )
        return case, specification, occurrence, bundle

    @staticmethod
    def _publication_manifest_identity(
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
        bundle: ResultBundle,
    ) -> str:
        return _identity(
            "cpmg-publication-manifest",
            (
                _SCHEMA_VERSION,
                _SEMANTIC_VERSION,
                case.identity,
                specification.identity,
                occurrence.lifecycle_identity,
                bundle.identity,
                bundle.manifest_identity,
            ),
        )

    @staticmethod
    def _validate_published_relationships(
        manifest: Mapping[str, object],
        identity: str,
        location: Path,
        case: CaseDefinition,
        specification: ExecutionSpecification,
        occurrence: Occurrence,
        bundle: ResultBundle,
    ) -> str:
        expected_manifest_identity = (
            CpmgBaselinePublisher._publication_manifest_identity(
                case, specification, occurrence, bundle
            )
        )
        if (
            specification.case_identity != case.identity
            or occurrence.case_identity != case.identity
            or occurrence.input_member_identities
            != tuple(sorted(member.identity for member in case.inputs))
            or bundle.identity != identity
            or occurrence.lifecycle != "SUCCEEDED"
            or occurrence.result_bundle_identity != bundle.identity
            or manifest.get("member_manifest_identity") != bundle.manifest_identity
            or manifest.get("manifest_identity") != expected_manifest_identity
        ):
            raise ValueError("Published baseline evidence has incompatible identity")
        try:
            raw_attempt = (location.parent.parent / "attempt.json").read_bytes()
        except OSError as error:
            raise ValueError(
                "Published baseline evidence lacks its terminal record"
            ) from error
        requested = Occurrence.from_record(
            _canonical_record_from_bytes(raw_attempt, "Baseline attempt record"),
            specification,
        )
        if requested != Occurrence(
            occurrence.execution_specification_identity,
            occurrence.case_identity,
            occurrence.actual_implementation_identity,
            occurrence.lane_reference,
            occurrence.input_member_identities,
            occurrence.attempt_token,
        ):
            raise ValueError(
                "Published baseline terminal state has incompatible identity"
            )
        return expected_manifest_identity

    @staticmethod
    def _validate_published_members(location: Path, bundle: ResultBundle) -> None:
        for member in bundle.members:
            try:
                content = (location / "members" / member.content_hash).read_bytes()
            except OSError as error:
                raise ValueError(
                    "Published baseline evidence has a missing member"
                ) from error
            if _sha256(content) != member.content_hash or len(content) != member.size:
                raise ValueError("Published baseline evidence member hash mismatch")

    def read(self, identity: str) -> PublishedEvidence:
        _record_hash(identity, "published bundle identity")
        location, manifest = self._find_success(identity)
        case, specification, occurrence, bundle = self._manifest_records(manifest)
        manifest_identity = self._validate_published_relationships(
            manifest, identity, location, case, specification, occurrence, bundle
        )
        self._validate_published_members(location, bundle)
        return PublishedEvidence(
            case, specification, occurrence, bundle, manifest_identity, location
        )


class LegacyAnchorExecutionError(RuntimeError):
    """Legacy product failure with a persisted failed occurrence record."""

    def __init__(self, occurrence: Occurrence) -> None:
        self.occurrence = occurrence
        super().__init__(f"Legacy anchor execution failed: {occurrence.failure_code}")


class BaselineLifecycleConflictError(RuntimeError):
    """An operation lost a race; its occurrence keeps the persisted terminal truth."""

    def __init__(self, occurrence: Occurrence) -> None:
        self.occurrence = occurrence
        super().__init__(
            "Baseline terminal publication conflicted with persisted "
            f"{occurrence.lifecycle} occurrence"
        )


@dataclass(frozen=True, slots=True)
class _CapturedInput:
    """One authoritative input captured as bytes before a CPMG attempt starts."""

    member: InputMember
    snapshot_relative_path: Path
    content: bytes


def _input_snapshot_path(anchor_directory: Path, source: Path) -> Path:
    try:
        return source.relative_to(anchor_directory)
    except ValueError as error:
        raise ValueError(
            "CPMG anchor inputs must remain inside the anchor directory"
        ) from error


def _read_cpmg_input(anchor_directory: Path, role: str, source: Path) -> _CapturedInput:
    content = source.read_bytes()
    return _CapturedInput(
        InputMember(role, _sha256(content), len(content)),
        _input_snapshot_path(anchor_directory, source),
        content,
    )


def _cpmg_experiment(content: bytes) -> Mapping[str, object]:
    try:
        decoded = tomllib.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise ValueError("Malformed frozen CPMG experiment input") from error
    if not isinstance(decoded, Mapping):
        raise TypeError("CPMG experiment input must be a record")
    return cast("Mapping[str, object]", decoded)


def _capture_cpmg_inputs(anchor_directory: Path) -> tuple[_CapturedInput, ...]:
    """Capture every authoritative input exactly once before creating a snapshot."""
    experiments = {
        "500mhz": anchor_directory / "Experiments" / "500mhz.toml",
        "800mhz": anchor_directory / "Experiments" / "800mhz.toml",
    }
    captured = [
        _read_cpmg_input(
            anchor_directory,
            "parameters",
            anchor_directory / "Parameters" / "parameters.toml",
        ),
        _read_cpmg_input(
            anchor_directory, "method", anchor_directory / "Methods" / "method.toml"
        ),
    ]
    experiment_inputs = {
        name: _read_cpmg_input(anchor_directory, f"experiment:{name}", path)
        for name, path in experiments.items()
    }
    captured.extend(experiment_inputs.values())
    for name, experiment_input in experiment_inputs.items():
        experiment = _cpmg_experiment(experiment_input.content)
        data = experiment.get("data")
        if not isinstance(data, Mapping):
            raise TypeError("CPMG experiment input has no data record")
        profiles = data.get("profiles")
        if not isinstance(profiles, Mapping):
            raise TypeError("CPMG experiment input has no profiles record")
        data_path = _record_string(data.get("path"), "CPMG data path")
        data_directory = experiments[name].parent / data_path
        for profile, filename in profiles.items():
            profile_name = _record_string(profile, "CPMG profile")
            captured.append(
                _read_cpmg_input(
                    anchor_directory,
                    f"data:{name}:{profile_name}",
                    data_directory / _record_string(filename, "CPMG data filename"),
                )
            )
    result = tuple(sorted(captured, key=lambda item: item.member))
    if len({item.member.role for item in result}) != len(result):
        raise ValueError("CPMG anchor has duplicate frozen input roles")
    return result


def _materialize_cpmg_snapshot(
    work_directory: Path, inputs: Sequence[_CapturedInput]
) -> Path:
    """Write the already-captured bytes into a private, read-only snapshot."""
    snapshot_directory = work_directory / "anchor"
    for captured in inputs:
        destination = snapshot_directory / captured.snapshot_relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(captured.content)
        destination.chmod(0o444)
    for directory in sorted(
        {path.parent for path in snapshot_directory.rglob("*") if path.is_file()},
        key=lambda path: len(path.parts),
        reverse=True,
    ):
        directory.chmod(0o555)
    return snapshot_directory


def _verify_cpmg_snapshot(
    snapshot_directory: Path, inputs: Sequence[_CapturedInput]
) -> None:
    """Reject a snapshot changed after its input-member records were frozen."""
    for captured in inputs:
        content = (snapshot_directory / captured.snapshot_relative_path).read_bytes()
        member = InputMember(captured.member.role, _sha256(content), len(content))
        if member != captured.member:
            raise ValueError("Frozen CPMG input snapshot drifted during execution")


def _case_from_cpmg_inputs(inputs: Sequence[_CapturedInput]) -> CaseDefinition:
    return CaseDefinition.create(
        "cpmg-15n-ip",
        CaseSourceAuthority(_CPMG_SOURCE_COMMIT, _CPMG_LOCKFILE_HASH),
        {
            "analysis_kind": "fit",
            "anchor_role": "scientific-anchor",
            "experiment": "cpmg_15n_ip",
            "model": "2st",
            "workflow": "shipped-two-step-method",
        },
        tuple(item.member for item in inputs),
    )


def cpmg_15n_ip_case(anchor_directory: Path) -> CaseDefinition:
    """Freeze #581's first shipped anchor without retaining source paths."""
    return _case_from_cpmg_inputs(_capture_cpmg_inputs(anchor_directory))


def _cpmg_required_artifact_roles(inputs: Sequence[_CapturedInput]) -> tuple[str, ...]:
    """Derive the closed output set from the two frozen shipped experiment TOMLs."""
    by_role = {item.member.role: item for item in inputs}
    profiles_by_field: list[tuple[str, ...]] = []
    for field_name in ("500mhz", "800mhz"):
        experiment = _cpmg_experiment(by_role[f"experiment:{field_name}"].content)
        data = experiment.get("data")
        if not isinstance(data, Mapping):
            raise TypeError("CPMG experiment input has no data record")
        profiles = data.get("profiles")
        if not isinstance(profiles, Mapping):
            raise TypeError("CPMG experiment input has no profiles record")
        profiles_by_field.append(
            tuple(_record_string(profile, "CPMG profile") for profile in profiles)
        )
    if profiles_by_field[0] != profiles_by_field[1]:
        raise ValueError("CPMG fields must declare the same frozen profile sequence")
    shared_outputs = (
        "Data/500mhz.dat",
        "Data/800mhz.dat",
        "Parameters/constrained.toml",
        "Parameters/fitted.toml",
        "Parameters/fixed.toml",
        "Plots/500mhz.fit",
        "Plots/800mhz.fit",
        "statistics.toml",
    )
    group_outputs = tuple(
        item for item in shared_outputs if not item.startswith("Plots/")
    )
    roles = [
        *(f"legacy-output:STEP1/{item}" for item in shared_outputs),
        *(f"legacy-output:STEP2/All/{item}" for item in shared_outputs),
    ]
    roles.extend(
        f"legacy-output:STEP2/Groups/{number}_{profile}/{item}"
        for number, profile in enumerate(profiles_by_field[0], start=1)
        for item in group_outputs
    )
    return _canonical_strings(sorted(roles), "CPMG required artifact roles")


def _cpmg_artifact_inventory(inputs: Sequence[_CapturedInput]) -> dict[str, object]:
    return {
        "version": "cpmg-structured-artifact-contract-v1",
        "closed": True,
        "excluded_path_components": ["run_info"],
        "structured_suffixes": [".dat", ".fit", ".toml", ".txt"],
        "required_roles": list(_cpmg_required_artifact_roles(inputs)),
    }


def cpmg_15n_ip_legacy_specification(
    case: CaseDefinition,
    *,
    artifact_inventory: Mapping[str, object],
    implementation: LegacyObservationImplementation | None = None,
    lane_reference: str = "unqualified-local-lane-v1",
) -> ExecutionSpecification:
    """The deliberately minimal pre-#588 execution request for the anchor."""
    return ExecutionSpecification.create(
        case,
        implementation or LegacyObservationImplementation.from_current_package(),
        workflow={
            "command": "fit",
            "method": "shipped-two-step-method",
            "model": "2st",
        },
        lane_reference=lane_reference,
        policy={"concurrency": "legacy-default", "implementation": "legacy-product"},
        budget={"policy": "legacy-product-default"},
        seed=None,
        execution_settings={
            "native_threads": "auto",
            "plot": "normal",
            "workers": "auto",
        },
        artifact_inventory=artifact_inventory,
        roles=["Scientific anchor"],
        claims=["legacy-observation-continuity"],
    )


def _assert_cpmg_request(
    case: CaseDefinition,
    specification: ExecutionSpecification,
) -> None:
    if specification.case_identity != case.identity:
        raise ValueError("Execution specification belongs to another case")
    if specification.implementation.authority_role != "LegacyObservationImplementation":
        raise ValueError("CPMG baseline requires the legacy observation implementation")


def _legacy_result_artifacts(
    output_directory: Path, specification: ExecutionSpecification
) -> tuple[ArtifactContent, ...]:
    inventory = specification.artifact_inventory.to_record_value()
    required_roles = CpmgBaselinePublisher._required_roles(specification)
    if not isinstance(inventory, Mapping):
        raise TypeError("CPMG artifact inventory must be a record")
    raw_suffixes = inventory.get("structured_suffixes")
    raw_excluded = inventory.get("excluded_path_components")
    if not isinstance(raw_suffixes, list) or not all(
        isinstance(suffix, str) for suffix in raw_suffixes
    ):
        raise ValueError("CPMG artifact inventory has malformed suffixes")
    if not isinstance(raw_excluded, list) or not all(
        isinstance(component, str) for component in raw_excluded
    ):
        raise ValueError("CPMG artifact inventory has malformed exclusions")
    allowed_suffixes = frozenset(raw_suffixes)
    excluded_components = frozenset(raw_excluded)
    artifacts = tuple(
        ArtifactContent(
            f"legacy-output:{path.relative_to(output_directory).as_posix()}",
            path.read_bytes(),
        )
        for path in sorted(output_directory.rglob("*"))
        if path.is_file()
        and not (excluded_components & set(path.relative_to(output_directory).parts))
        and path.suffix in allowed_suffixes
    )
    if tuple(item.role for item in artifacts) != required_roles:
        raise ValueError("Legacy anchor did not produce the closed required artifacts")
    return artifacts


def capture_cpmg_15n_ip_legacy_observation(
    *,
    publisher: CpmgBaselinePublisher,
    anchor_directory: Path,
    lane_reference: str = "unqualified-local-lane-v1",
    attempt_token: str,
) -> PublishedEvidence:
    """Run the unchanged CPMG analysis through the legacy product runner.

    The temporary output is sampled only after the legacy run ends.  No baseline
    artifact is supplied to the analysis, and a failure is never published.
    """
    work_directory = Path(tempfile.mkdtemp(prefix="chemex-cpmg-baseline-"))
    requested: Occurrence | None = None
    reserved = False
    try:
        captured_inputs = _capture_cpmg_inputs(anchor_directory)
        case = _case_from_cpmg_inputs(captured_inputs)
        specification = cpmg_15n_ip_legacy_specification(
            case,
            artifact_inventory=_cpmg_artifact_inventory(captured_inputs),
            lane_reference=lane_reference,
        )
        _assert_cpmg_request(case, specification)
        requested = Occurrence.requested(specification, case, attempt_token)
        publisher.reserve(case, specification, requested)
        reserved = True
        snapshot_directory = _materialize_cpmg_snapshot(work_directory, captured_inputs)
        output_directory = work_directory / "Output"
        paths = {
            item.member.role: snapshot_directory / item.snapshot_relative_path
            for item in captured_inputs
        }
        parser = build_parser()
        args = parser.parse_args(
            [
                "fit",
                "-e",
                str(paths["experiment:500mhz"]),
                str(paths["experiment:800mhz"]),
                "-p",
                str(paths["parameters"]),
                "-m",
                str(paths["method"]),
                "-o",
                str(output_directory),
            ]
        )
        run(
            args,
            argv=(
                "chemex",
                "fit",
                "-e",
                "Experiments/*.toml",
                "-p",
                "Parameters/parameters.toml",
                "-m",
                "Methods/method.toml",
                "-o",
                "Output",
            ),
        )
        artifacts = _legacy_result_artifacts(output_directory, specification)
        _verify_cpmg_snapshot(snapshot_directory, captured_inputs)
        bundle = ResultBundle.create(
            requested.identity,
            specification.identity,
            specification.implementation,
            tuple(artifact.member for artifact in artifacts),
        )
        return publisher.publish(
            case, specification, requested.succeeded(bundle), bundle, artifacts
        )
    except BaselineLifecycleConflictError:
        raise
    except Exception as error:
        if requested is None or not reserved:
            raise
        failed = publisher.record_failure(
            case, specification, requested, type(error).__name__
        )
        raise LegacyAnchorExecutionError(failed) from error
    finally:
        shutil.rmtree(work_directory, ignore_errors=True)
