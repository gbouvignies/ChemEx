"""Compact lifecycle evidence and fixed semantic predicates for #592."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import cast

from chemex.baselines import (
    CanonicalBaselineValue,
    CaseDefinition,
    ExecutionSpecification,
    Occurrence,
)
from chemex.optimize.method_step import MethodStepOutcome

FAILURE_REQUIREMENTS = tuple(
    f"migration-core.failure.{name}"
    for name in (
        "construction",
        "execution",
        "materialization",
        "commit",
        "partial-stochastic-evidence",
        "workflow-stop",
        "publication",
    )
)
CAPTURE_VERSION = "migration-core-lifecycle-probes-v1"
_PROBE_IDS = (
    "construction-dependent-stop",
    "primary-execution-failure",
    "aggregate-materialization-failure",
    "aggregate-acceptance-stop",
    "accepted-commit-rejected",
    "cancelled-resampling-partial-publication",
    "worker-resampling-failure",
    "requested-resampling-summary-failure",
    "interrupted-resampling-partial-publication",
    "committed-publication-failure",
)
_FACT_KEYS = {
    "state",
    "construction",
    "method_step",
    "primary",
    "checkpoint",
    "commit",
    "resampling",
    "summary_failure_category",
    "successor",
}


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, allow_nan=False, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")


def _identity(kind: str, value: object) -> str:
    return hashlib.sha256(
        _canonical_bytes({"kind": kind, "schema_version": 1, "value": value})
    ).hexdigest()


def _record(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"Lifecycle probe {name} must be a record")
    return cast("Mapping[str, object]", value)


@dataclass(frozen=True, slots=True)
class LifecycleProbeResult:
    """One probe's existing baseline lineage and owner-structured runtime facts."""

    probe_id: str
    case: CaseDefinition
    specification: ExecutionSpecification
    occurrence: Occurrence
    facts: CanonicalBaselineValue
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        expected_inputs = tuple(sorted(member.identity for member in self.case.inputs))
        if self.probe_id not in _PROBE_IDS:
            raise ValueError("Lifecycle result contains an unknown probe")
        if (
            self.specification.case_identity != self.case.identity
            or self.occurrence.execution_specification_identity
            != self.specification.identity
            or self.occurrence.case_identity != self.case.identity
            or self.occurrence.actual_implementation_identity
            != self.specification.implementation.identity
            or self.occurrence.lane_reference != self.specification.lane_reference
            or self.occurrence.input_member_identities != expected_inputs
            or self.occurrence.lifecycle != "REQUESTED"
        ):
            raise ValueError("Lifecycle occurrence lineage is inconsistent")
        facts = self.runtime_facts
        policy = self.specification.policy.to_record_value()
        workflow = self.specification.workflow.to_record_value()
        settings = self.specification.execution_settings.to_record_value()
        state = facts.get("state")
        if (
            set(facts) != _FACT_KEYS
            or not isinstance(policy, Mapping)
            or policy.get("probe_id") != self.probe_id
            or not isinstance(workflow, Mapping)
            or not isinstance(settings, Mapping)
            or not isinstance(settings.get("environment_identity"), str)
            or not isinstance(state, Mapping)
            or workflow.get("starting_revision") != state.get("starting_revision")
        ):
            raise ValueError("Lifecycle execution specification is incompatible")
        outcome = self.method_step_outcome
        if outcome is not None and (
            workflow.get("semantic_identity") != outcome.workflow_identity
            or workflow.get("binding_identity") != outcome.workflow_binding_identity
        ):
            raise ValueError("Lifecycle outcome belongs to another workflow")
        object.__setattr__(
            self, "identity", _identity("lifecycle-probe-result-v1", self._payload())
        )

    @property
    def runtime_facts(self) -> Mapping[str, object]:
        value = self.facts.to_record_value()
        if not isinstance(value, Mapping):
            raise TypeError("Lifecycle runtime facts must be a record")
        return cast("Mapping[str, object]", value)

    @property
    def method_step_outcome(self) -> MethodStepOutcome | None:
        value = self.runtime_facts.get("method_step")
        return (
            None
            if value is None
            else MethodStepOutcome.from_record(dict(_record(value, "method step")))
        )

    @property
    def environment_identity(self) -> str:
        settings = self.specification.execution_settings.to_record_value()
        if not isinstance(settings, Mapping) or not isinstance(
            settings.get("environment_identity"), str
        ):
            raise TypeError("Lifecycle environment identity is unavailable")
        typed = cast("Mapping[str, object]", settings)
        return cast("str", typed["environment_identity"])

    def _payload(self) -> dict[str, object]:
        return {
            "probe_id": self.probe_id,
            "case": self.case.to_record(),
            "specification": self.specification.to_record(),
            "occurrence": self.occurrence.to_record(),
            "facts": self.facts.to_record_value(),
        }

    def to_record(self) -> dict[str, object]:
        return {
            "artifact_type": "migration_core_lifecycle_probe_result",
            "schema_version": 1,
            **self._payload(),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> LifecycleProbeResult:
        if (
            set(record)
            != {
                "artifact_type",
                "schema_version",
                "probe_id",
                "case",
                "specification",
                "occurrence",
                "facts",
                "identity",
            }
            or record.get("artifact_type") != "migration_core_lifecycle_probe_result"
            or record.get("schema_version") != 1
        ):
            raise ValueError("Malformed lifecycle probe result")
        case = CaseDefinition.from_record(_record(record.get("case"), "case"))
        specification = ExecutionSpecification.from_record(
            _record(record.get("specification"), "specification")
        )
        occurrence = Occurrence.from_record(
            _record(record.get("occurrence"), "occurrence"), specification
        )
        probe_id = record.get("probe_id")
        if not isinstance(probe_id, str):
            raise TypeError("Lifecycle probe ID must be a string")
        result = cls(
            probe_id,
            case,
            specification,
            occurrence,
            CanonicalBaselineValue.from_record(record.get("facts"), "runtime facts"),
        )
        if record.get("identity") != result.identity:
            raise ValueError("Lifecycle probe result identity does not match payload")
        return result


@dataclass(frozen=True, slots=True)
class LifecycleProbeCapture:
    records: tuple[LifecycleProbeResult, ...]
    capture_version: str = CAPTURE_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        identifiers = tuple(record.probe_id for record in self.records)
        if self.capture_version != CAPTURE_VERSION or len(set(identifiers)) != len(
            identifiers
        ):
            raise ValueError(
                "Lifecycle capture version or probe identities are invalid"
            )
        object.__setattr__(
            self, "identity", _identity("lifecycle-probe-capture-v1", self._payload())
        )

    def _payload(self) -> dict[str, object]:
        return {
            "capture_version": self.capture_version,
            "records": [record.to_record() for record in self.records],
        }

    def to_bytes(self) -> bytes:
        return _canonical_bytes(
            {"schema_version": 1, **self._payload(), "identity": self.identity}
        )

    @classmethod
    def from_bytes(cls, content: bytes) -> LifecycleProbeCapture:
        raw = json.loads(content)
        if (
            not isinstance(raw, Mapping)
            or set(raw) != {"capture_version", "identity", "records", "schema_version"}
            or raw.get("schema_version") != 1
            or not isinstance(raw.get("records"), list)
        ):
            raise ValueError("Malformed lifecycle probe capture")
        records = cast("list[object]", raw["records"])
        capture = cls(
            tuple(
                LifecycleProbeResult.from_record(_record(item, "result"))
                for item in records
            )
        )
        if raw.get("identity") != capture.identity:
            raise ValueError("Lifecycle probe capture identity does not match payload")
        return capture


@dataclass(frozen=True, slots=True)
class _Expectation:
    probe_id: str
    requirements: tuple[int, ...]
    ending_revision: int
    outcome: tuple[str, str, bool, bool, bool] | None
    construction: tuple[str, tuple[str, ...], tuple[str, ...]] | None
    primary: tuple[str, tuple[str, ...]] | None
    checkpoint: tuple[tuple[str, ...], str | None]
    commit: tuple[str, str | None] | None
    resampling: (
        tuple[str, tuple[int, ...], tuple[int, int, int], tuple[str, ...]] | None
    )
    derivation: str | None
    summary_failure: str | None
    successor: tuple[bool, tuple[str, ...]]
    publication_failed: bool = False


# Each count below belongs only to ResamplingEvidence; there is no cross-stage total.
# fmt: off
_EXPECTATIONS = tuple(_Expectation(*row) for row in (
    ("construction-dependent-stop", (0, 5), 0, None, ("DirectTrfConstructionError", ("first",), ()), None, ((), None), None, None, None, None, (False, ())),
    ("primary-execution-failure", (1, 6), 0, ("failed", "execution_failure", False, False, True), None, ("execution_failure", ("not_started",)), ((), None), None, None, None, None, (True, ())),
    ("aggregate-materialization-failure", (2, 5), 0, ("failed", "decomposition_validation_failure", False, False, False), None, ("decomposition_validation_failure", ("succeeded",)), ((), None), None, None, None, None, (True, ())),
    ("aggregate-acceptance-stop", (3, 5, 6), 0, ("failed", "accepted", True, False, False), None, ("accepted", ()), (("aggregate_accepted",), "RuntimeError"), None, None, None, None, (True, ())),
    ("accepted-commit-rejected", (3, 5, 6), 0, ("accepted_uncommitted", "accepted", True, True, False), None, ("accepted", ("succeeded",)), ((), None), ("failed", "incompatible_state"), None, None, None, (True, ())),
    ("cancelled-resampling-partial-publication", (4, 6), 1, ("committed", "accepted", True, True, True), None, ("accepted", ("succeeded",)), ((), None), ("committed", None), ("cancelled", (0, 1), (0, 0, 0), ("not_started", "not_started")), "cancelled", "insufficient_successful_coverage", (False, ("dependent",))),
    ("worker-resampling-failure", (1, 4), 1, ("committed", "accepted", True, True, True), None, ("accepted", ("succeeded",)), ((), None), ("committed", None), ("completed", (), (2, 1, 1), ("failed", "succeeded")), "failed", "non_finite_summary_arithmetic", (False, ("dependent",))),
    ("requested-resampling-summary-failure", (4,), 1, ("committed", "accepted", True, True, True), None, ("accepted", ("succeeded",)), ((), None), ("committed", None), ("completed", (), (2, 2, 0), ("succeeded", "succeeded")), "failed", "qualification_summary_failure", (False, ("dependent",))),
    ("interrupted-resampling-partial-publication", (4, 5), 1, ("committed", "accepted", True, True, True), None, ("accepted", ("succeeded",)), ((), None), ("committed", None), ("interrupted", (1,), (1, 0, 0), ("interrupted", "not_started")), "interrupted", "insufficient_successful_coverage", (False, ("dependent",))),
    ("committed-publication-failure", (5, 6), 1, ("publication_failed", "accepted", True, True, False), None, ("accepted", ("succeeded",)), ((), None), ("committed", None), None, None, None, (True, ()), True),
))
# fmt: on


def _owned_facts(expected: _Expectation) -> dict[str, object]:
    construction = expected.construction
    primary = expected.primary
    commit = expected.commit
    resampling = expected.resampling
    return {
        "state": {"starting_revision": 0, "ending_revision": expected.ending_revision},
        "construction": None
        if construction is None
        else {
            "failure_category": construction[0],
            "constructor_entries": list(construction[1]),
            "executor_entries": list(construction[2]),
        },
        "primary": None
        if primary is None
        else {
            "terminal": primary[0],
            "component_dispositions": list(primary[1]),
        },
        "checkpoint": {
            "entries": list(expected.checkpoint[0]),
            "failure_category": expected.checkpoint[1],
        },
        "commit": None
        if commit is None
        else {
            "terminal": commit[0],
            "failure_category": commit[1],
        },
        "resampling": None
        if resampling is None
        else {
            "terminal": resampling[0],
            "unstarted_ordinals": list(resampling[1]),
            "completed_count": resampling[2][0],
            "successful_count": resampling[2][1],
            "failed_count": resampling[2][2],
            "replicate_dispositions": list(resampling[3]),
        },
        "summary_failure_category": expected.summary_failure,
        "successor": {
            "denied": expected.successor[0],
            "downstream_entries": list(expected.successor[1]),
        },
    }


def _matches(record: LifecycleProbeResult, expected: _Expectation) -> bool:
    facts = dict(record.runtime_facts)
    facts.pop("method_step")
    if facts != _owned_facts(expected):
        return False
    outcome = record.method_step_outcome
    if expected.outcome is None:
        return outcome is None
    if outcome is None:
        return False
    lifecycle, primary, accepted, committed, published = expected.outcome
    if (
        outcome.lifecycle.value != lifecycle
        or outcome.primary_terminal != primary
        or (outcome.accepted_result_identity is not None) != accepted
        or (outcome.commit_operation_identity is not None) != committed
        or (outcome.publication_identity is not None) != published
        or bool(outcome.publication_failure) != expected.publication_failed
    ):
        return False
    if expected.derivation is None:
        return not outcome.derivations
    return (
        len(outcome.derivations) == 1
        and outcome.derivations[0].stage == "resampling"
        and outcome.derivations[0].disposition.value == expected.derivation
        and outcome.derivations[0].operation_identity is not None
    )


def eligible_failure_requirements(
    capture: LifecycleProbeCapture,
    *,
    source_commit: str,
    lockfile_hash: str,
    lane_identity: str,
    attestation_identity: str | None,
    environment_identity: str,
) -> tuple[str, ...]:
    """Map reconstructed lineage and owner-local facts to fixed requirements."""

    records = {record.probe_id: record for record in capture.records}
    if set(records) != set(_PROBE_IDS):
        return ()

    def provenance(record: LifecycleProbeResult) -> bool:
        source = record.case.source_authority
        settings = record.specification.execution_settings.to_record_value()
        return (
            source.source_commit == source_commit
            and source.lockfile_hash == lockfile_hash
            and record.specification.lane_reference == lane_identity
            and record.occurrence.lane_attestation_identity == attestation_identity
            and isinstance(settings, Mapping)
            and settings.get("environment_identity") == environment_identity
            and settings.get("workers") == settings.get("native_threads") == 1
        )

    if not all(provenance(record) for record in records.values()):
        return ()
    valid = {
        expected.probe_id: _matches(records[expected.probe_id], expected)
        for expected in _EXPECTATIONS
    }
    return tuple(
        requirement
        for index, requirement in enumerate(FAILURE_REQUIREMENTS)
        if all(
            valid[expected.probe_id]
            for expected in _EXPECTATIONS
            if index in expected.requirements
        )
    )


__all__ = [
    "FAILURE_REQUIREMENTS",
    "LifecycleProbeCapture",
    "LifecycleProbeResult",
    "eligible_failure_requirements",
]
