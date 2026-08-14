"""Compact lifecycle-result records and fixed semantic checks for #592."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field, fields
from typing import Any, Literal, cast

# fmt: off
type _Stage = Literal["acceptance", "commit", "construction", "materialization", "primary_execution", "publication", "requested_evidence", "worker"]
type _Terminal = Literal["accepted_uncommitted", "committed", "failed", "publication_failed"]
type _Acceptance = Literal["accepted", "failed", "not_started"]
type _Commit = Literal["committed", "not_started", "rejected"]
type _Publication = Literal["diagnostics_only", "failed", "partial_only", "suppressed"]
type _Successor = Literal["available", "denied"]
type _PartialEvidence = Literal["cancelled_population_retained", "interrupted_population_retained", "not_applicable", "summary_failure_retained", "worker_failure_retained"]
# fmt: on

FAILURE_REQUIREMENTS = (
    "migration-core.failure.construction",
    "migration-core.failure.execution",
    "migration-core.failure.materialization",
    "migration-core.failure.commit",
    "migration-core.failure.partial-stochastic-evidence",
    "migration-core.failure.workflow-stop",
    "migration-core.failure.publication",
)
CAPTURE_VERSION = "migration-core-lifecycle-probes-v1"


@dataclass(frozen=True, slots=True)
class _ProbeDefinition:
    identifier: str
    requirement_indexes: tuple[int, ...]
    stage: _Stage
    terminal: _Terminal
    revisions: tuple[int, int]
    counts: tuple[int, int, int, int, int, int]
    dispositions: tuple[
        _Acceptance, _Commit, _Publication, _Successor, _PartialEvidence
    ]


# Fixed claim semantics are deliberately separate from capture execution.
# fmt: off
_PROBES = tuple(_ProbeDefinition(*row) for row in (
    ("construction-dependent-stop", (0, 5), "construction", "failed", (0, 0), (1, 0, 1, 0, 0, 3), ("not_started", "not_started", "suppressed", "denied", "not_applicable")),
    ("primary-execution-failure", (1, 6), "primary_execution", "failed", (0, 0), (1, 0, 1, 0, 0, 3), ("not_started", "not_started", "diagnostics_only", "denied", "not_applicable")),
    ("aggregate-materialization-failure", (2, 5), "materialization", "failed", (0, 0), (2, 1, 1, 0, 0, 2), ("not_started", "not_started", "suppressed", "denied", "not_applicable")),
    ("aggregate-acceptance-stop", (3, 5, 6), "acceptance", "failed", (0, 0), (2, 1, 1, 0, 0, 2), ("failed", "not_started", "suppressed", "denied", "not_applicable")),
    ("accepted-commit-rejected", (3, 5, 6), "commit", "accepted_uncommitted", (0, 0), (3, 2, 1, 0, 0, 1), ("accepted", "rejected", "suppressed", "denied", "not_applicable")),
    ("cancelled-resampling-partial-publication", (4, 6), "requested_evidence", "committed", (0, 1), (3, 2, 0, 1, 0, 2), ("accepted", "committed", "partial_only", "available", "cancelled_population_retained")),
    ("worker-resampling-failure", (1, 4), "worker", "committed", (0, 1), (4, 3, 1, 0, 0, 0), ("accepted", "committed", "partial_only", "available", "worker_failure_retained")),
    ("requested-resampling-summary-failure", (4,), "requested_evidence", "committed", (0, 1), (5, 4, 1, 0, 0, 0), ("accepted", "committed", "partial_only", "available", "summary_failure_retained")),
    ("interrupted-resampling-partial-publication", (4, 5), "worker", "committed", (0, 1), (3, 2, 0, 0, 1, 1), ("accepted", "committed", "partial_only", "available", "interrupted_population_retained")),
    ("committed-publication-failure", (5, 6), "publication", "publication_failed", (0, 1), (3, 2, 1, 0, 0, 0), ("accepted", "committed", "failed", "denied", "not_applicable")),
))
# fmt: on
_PROBE_BY_ID = {probe.identifier: probe for probe in _PROBES}


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, allow_nan=False, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")


def _identity(kind: str, value: object) -> str:
    return hashlib.sha256(
        _canonical_bytes({"kind": kind, "schema_version": 1, "value": value})
    ).hexdigest()


def _text(record: Mapping[str, object], key: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"Lifecycle probe {key} must be a non-empty string")
    return value


def _count(record: Mapping[str, object], key: str) -> int:
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Lifecycle probe {key} must be a non-negative integer")
    return value


@dataclass(frozen=True, slots=True)
class LifecycleOperationCounts:
    attempted: int
    completed: int
    failed: int
    cancelled: int
    interrupted: int
    unstarted: int

    def __post_init__(self) -> None:
        values = tuple(asdict(self).values())
        if any(isinstance(value, bool) or value < 0 for value in values):
            raise ValueError("Lifecycle operation counts must be non-negative integers")
        if sum(values[1:5]) > self.attempted:
            raise ValueError("Lifecycle terminal counts exceed attempted operations")

    def to_record(self) -> dict[str, int]:
        return cast("dict[str, int]", asdict(self))

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> LifecycleOperationCounts:
        names = tuple(item.name for item in fields(cls))
        if set(record) != set(names):
            raise ValueError("Lifecycle counts have unknown or missing fields")
        return cls(*(_count(record, name) for name in names))


@dataclass(frozen=True, slots=True)
class LifecycleProbeResult:
    """One compact typed terminal artifact emitted by a probe occurrence."""

    probe_id: str
    source_commit: str
    lockfile_hash: str
    lane_identity: str
    attestation_identity: str
    environment_identity: str
    case_identity: str
    execution_specification_identity: str
    occurrence_identity: str
    starting_revision: int
    ending_revision: int
    stage: _Stage
    terminal: _Terminal
    counts: LifecycleOperationCounts
    acceptance: _Acceptance
    commit: _Commit
    publication: _Publication
    successor: _Successor
    partial_evidence: _PartialEvidence
    identity: str = field(init=False)

    # fmt: off
    def __post_init__(self) -> None:
        text_values = (
            self.probe_id, self.source_commit, self.lockfile_hash, self.lane_identity,
            self.attestation_identity, self.environment_identity, self.case_identity,
            self.execution_specification_identity, self.occurrence_identity,
        )
        if any(not value for value in text_values):
            raise ValueError("Lifecycle result identities must be non-empty")
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in (self.starting_revision, self.ending_revision)):
            raise ValueError("Lifecycle revisions must be non-negative integers")
        allowed = tuple({getattr(probe, name) for probe in _PROBES} for name in ("stage", "terminal")) + tuple({probe.dispositions[index] for probe in _PROBES} for index in range(5))
        actual = (
            self.stage, self.terminal, self.acceptance, self.commit,
            self.publication, self.successor, self.partial_evidence,
        )
        if any(value not in choices for value, choices in zip(actual, allowed, strict=True)):
            raise ValueError("Lifecycle result contains an unknown terminal value")
        object.__setattr__(self, "identity", _identity("lifecycle-probe-result-v1", self._payload()))

    def _payload(self) -> dict[str, object]:
        result = {item.name: getattr(self, item.name) for item in fields(self) if item.name != "identity"}
        result["counts"] = self.counts.to_record()
        return result

    def to_record(self) -> dict[str, object]:
        return {
            "artifact_type": "migration_core_lifecycle_probe_result",
            "schema_version": 1,
            **self._payload(),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> LifecycleProbeResult:
        names = {item.name for item in fields(cls)}
        if set(record) != names | {"artifact_type", "schema_version"} or record.get("artifact_type") != "migration_core_lifecycle_probe_result" or record.get("schema_version") != 1:
            raise ValueError("Malformed lifecycle probe result")
        counts = record.get("counts")
        if not isinstance(counts, Mapping):
            raise TypeError("Lifecycle probe counts must be a record")
        kwargs: dict[str, Any] = {name: _text(record, name) for name in names - {"identity", "counts", "starting_revision", "ending_revision"}}
        result = cls(
            **kwargs,
            starting_revision=_count(record, "starting_revision"),
            ending_revision=_count(record, "ending_revision"),
            counts=LifecycleOperationCounts.from_record(cast("Mapping[str, object]", counts)),
        )
        if record.get("identity") != result.identity:
            raise ValueError("Lifecycle probe result identity does not match payload")
        return result
    # fmt: on


@dataclass(frozen=True, slots=True)
class LifecycleProbeCapture:
    records: tuple[LifecycleProbeResult, ...]
    capture_version: str = CAPTURE_VERSION
    identity: str = field(init=False)

    # fmt: off
    def __post_init__(self) -> None:
        identifiers = tuple(record.probe_id for record in self.records)
        if self.capture_version != CAPTURE_VERSION or len(set(identifiers)) != len(identifiers):
            raise ValueError("Lifecycle capture version or probe identities are invalid")
        object.__setattr__(self, "identity", _identity("lifecycle-probe-capture-v1", self._payload()))

    def _payload(self) -> dict[str, object]:
        return {
            "capture_version": self.capture_version,
            "records": [record.to_record() for record in self.records],
        }

    def to_bytes(self) -> bytes:
        return _canonical_bytes({"schema_version": 1, **self._payload(), "identity": self.identity})

    @classmethod
    def from_bytes(cls, content: bytes) -> LifecycleProbeCapture:
        raw = json.loads(content)
        if not isinstance(raw, Mapping) or set(raw) != {"capture_version", "identity", "records", "schema_version"} or raw.get("schema_version") != 1:
            raise ValueError("Malformed lifecycle probe capture")
        records = raw.get("records")
        if not isinstance(records, list) or not all(isinstance(item, Mapping) for item in records):
            raise ValueError("Lifecycle probe capture records must be a list")
        capture = cls(tuple(LifecycleProbeResult.from_record(cast("Mapping[str, object]", item)) for item in records))
        if raw.get("identity") != capture.identity:
            raise ValueError("Lifecycle probe capture identity does not match payload")
        return capture
    # fmt: on


def eligible_failure_requirements(
    capture: LifecycleProbeCapture,
    *,
    source_commit: str,
    lockfile_hash: str,
    lane_identity: str,
    attestation_identity: str,
    environment_identity: str,
) -> tuple[str, ...]:
    """Map only fixed terminal semantics—not claims or stored verdicts—to coverage."""

    # fmt: off
    records = {record.probe_id: record for record in capture.records}
    if set(records) != set(_PROBE_BY_ID) or not all(
        (record.source_commit, record.lockfile_hash, record.lane_identity, record.attestation_identity, record.environment_identity)
        == (source_commit, lockfile_hash, lane_identity, attestation_identity, environment_identity)
        for record in records.values()
    ):
        return ()

    def matches(probe: _ProbeDefinition) -> bool:
        result = records[probe.identifier]
        return (result.stage, result.terminal, (result.starting_revision, result.ending_revision), tuple(asdict(result.counts).values()), (result.acceptance, result.commit, result.publication, result.successor, result.partial_evidence)) == (probe.stage, probe.terminal, probe.revisions, probe.counts, probe.dispositions)

    return tuple(
        requirement
        for index, requirement in enumerate(FAILURE_REQUIREMENTS)
        if all(matches(probe) for probe in _PROBES if index in probe.requirement_indexes)
    )
    # fmt: on


__all__ = [
    "FAILURE_REQUIREMENTS",
    "LifecycleOperationCounts",
    "LifecycleProbeCapture",
    "LifecycleProbeResult",
    "eligible_failure_requirements",
]
