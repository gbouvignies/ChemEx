"""Compact serialization, process, cache, and replay evidence for issue #593."""

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
    ResultBundle,
)

CAPTURE_VERSION = "migration-core-operational-replay-v1"
OPERATIONAL_REQUIREMENTS = frozenset(
    {
        "migration-core.operational.cache",
        "migration-core.operational.deterministic-replay",
        "migration-core.operational.multiprocessing",
        "migration-core.operational.serialization",
        "migration-core.resampling.serial-two-worker-replay",
    }
)
_FACT_KEYS = {
    "cache",
    "deterministic_replay",
    "fail_closed",
    "multiprocessing",
    "serialization",
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
        raise TypeError(f"Operational replay {name} must be a record")
    return cast("Mapping[str, object]", value)


@dataclass(frozen=True, slots=True)
class OperationalReplayCapture:
    """One canonical operational probe with baseline lineage and compact facts."""

    case: CaseDefinition
    specification: ExecutionSpecification
    occurrence: Occurrence
    bundle: ResultBundle
    facts: CanonicalBaselineValue
    capture_version: str = CAPTURE_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        expected_inputs = tuple(sorted(member.identity for member in self.case.inputs))
        if (
            self.capture_version != CAPTURE_VERSION
            or self.case.name != "migration-core-operational-replay"
            or self.specification.case_identity != self.case.identity
            or self.occurrence.execution_specification_identity
            != self.specification.identity
            or self.occurrence.case_identity != self.case.identity
            or self.occurrence.actual_implementation_identity
            != self.specification.implementation.identity
            or self.occurrence.lane_reference != self.specification.lane_reference
            or self.occurrence.input_member_identities != expected_inputs
            or self.occurrence.lifecycle != "SUCCEEDED"
            or self.occurrence.result_bundle_identity != self.bundle.identity
            or self.bundle.occurrence_identity != self.occurrence.identity
            or self.bundle.execution_specification_identity
            != self.specification.identity
            or self.bundle.implementation != self.specification.implementation
            or self.specification.roles
            != ("qualification:migration-core-operational-replay",)
            or self.specification.claims != ("typed-runtime-facts",)
        ):
            raise ValueError("Operational replay lineage is inconsistent")
        facts = self.runtime_facts
        facts_content = _canonical_bytes(facts)
        if (
            set(facts) != _FACT_KEYS
            or len(self.bundle.members) != 1
            or self.bundle.members[0].role != "operational-replay-facts.json"
            or self.bundle.members[0].content_hash
            != hashlib.sha256(facts_content).hexdigest()
            or self.bundle.members[0].size != len(facts_content)
        ):
            raise ValueError("Operational replay facts are incomplete")
        object.__setattr__(
            self,
            "identity",
            _identity("operational-replay-capture-v1", self._payload()),
        )

    @property
    def runtime_facts(self) -> Mapping[str, object]:
        return _record(self.facts.to_record_value(), "facts")

    def _payload(self) -> dict[str, object]:
        return {
            "capture_version": self.capture_version,
            "case": self.case.to_record(),
            "specification": self.specification.to_record(),
            "occurrence": self.occurrence.to_record(),
            "bundle": self.bundle.to_record(),
            "facts": self.facts.to_record_value(),
        }

    def to_bytes(self) -> bytes:
        return _canonical_bytes(
            {"schema_version": 1, **self._payload(), "identity": self.identity}
        )

    @classmethod
    def from_bytes(cls, content: bytes) -> OperationalReplayCapture:
        record = json.loads(content)
        if not isinstance(record, Mapping) or set(record) != {
            "schema_version",
            "capture_version",
            "case",
            "specification",
            "occurrence",
            "bundle",
            "facts",
            "identity",
        }:
            raise ValueError("Malformed operational replay capture")
        typed = cast("Mapping[str, object]", record)
        if typed.get("schema_version") != 1:
            raise ValueError("Unsupported operational replay capture schema")
        case = CaseDefinition.from_record(_record(typed.get("case"), "case"))
        specification = ExecutionSpecification.from_record(
            _record(typed.get("specification"), "specification")
        )
        occurrence = Occurrence.from_record(
            _record(typed.get("occurrence"), "occurrence"), specification
        )
        bundle = ResultBundle.from_record(_record(typed.get("bundle"), "bundle"))
        capture_version = typed.get("capture_version")
        if not isinstance(capture_version, str):
            raise TypeError("Operational replay capture version must be a string")
        capture = cls(
            case,
            specification,
            occurrence,
            bundle,
            CanonicalBaselineValue.from_record(typed.get("facts"), "operational facts"),
            capture_version,
        )
        if typed.get("identity") != capture.identity:
            raise ValueError("Operational replay identity does not match payload")
        return capture


def eligible_operational_requirements(
    capture: OperationalReplayCapture,
    *,
    source_commit: str,
    lockfile_hash: str,
    lane_identity: str,
    attestation_identity: str,
    environment_identity: str,
) -> frozenset[str]:
    """Return exact claims only when lineage and every compact fact qualify."""

    settings = capture.specification.execution_settings.to_record_value()
    facts = capture.runtime_facts
    serialization = _record(facts.get("serialization"), "serialization facts")
    multiprocessing = _record(facts.get("multiprocessing"), "process facts")
    cache = _record(facts.get("cache"), "cache facts")
    replay = _record(facts.get("deterministic_replay"), "replay facts")
    fail_closed = _record(facts.get("fail_closed"), "failure facts")
    if (
        capture.case.source_authority.source_commit != source_commit
        or capture.case.source_authority.lockfile_hash != lockfile_hash
        or capture.specification.lane_reference != lane_identity
        or capture.occurrence.lane_attestation_identity != attestation_identity
        or not isinstance(settings, Mapping)
        or settings.get("environment_identity") != environment_identity
        or settings.get("workers") != 1
        or settings.get("native_threads") != 1
        or serialization.get("plan_round_trip") is not True
        or serialization.get("frame_round_trip") is not True
        or serialization.get("result_round_trip") is not True
        or multiprocessing.get("fresh_worker_ownership") is not True
        or multiprocessing.get("request_order_independent") is not True
        or multiprocessing.get("interruption_suppressed_result") is not True
        or cache.get("hits") != 1
        or cache.get("misses") != 2
        or cache.get("changed_frame_invalidated") is not True
        or replay.get("serial_spawn_identity_match") is not True
        or replay.get("spawn_replay_identity_match") is not True
        or fail_closed.get("malformed_request_rejected") is not True
        or fail_closed.get("stale_request_rejected") is not True
    ):
        return frozenset()
    return OPERATIONAL_REQUIREMENTS
