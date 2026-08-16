# ruff: noqa: I001

import dataclasses
import json
from pathlib import Path

# fmt: off
import pytest

from chemex import migration_core
from chemex.baselines import CaseSourceAuthority
from chemex.migration_core_lifecycle import (
    FAILURE_REQUIREMENTS,
    LifecycleProbeCapture,
    eligible_failure_requirements,
)
from chemex.numerical_lanes import LiveLaneAuthority, RuntimeEnvironment, canonical_lanes
from tests.qualification.capture_migration_core_lifecycle import (
    capture_lifecycle_probes,
    observe_lifecycle_probes,
)
FIXTURE = Path(__file__).parent / "fixtures/migration_core_lifecycle_probes_v1.json"
SOURCE_COMMIT = "0" * 40
LOCKFILE_HASH = "1" * 64


def _live_authority(monkeypatch: pytest.MonkeyPatch) -> LiveLaneAuthority:
    lane = canonical_lanes()[0]
    environment = RuntimeEnvironment(lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: environment),
    )
    return lane.attest_current_process(lane.semantics.image_digest)


@pytest.fixture
def lifecycle_capture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> LifecycleProbeCapture:
    authority = _live_authority(monkeypatch)
    observations = observe_lifecycle_probes(
        tmp_path,
        authority=authority,
        source_commit=SOURCE_COMMIT,
        lockfile_hash=LOCKFILE_HASH,
    )
    return capture_lifecycle_probes(observations)


def test_capture_uses_real_baseline_lineage_and_owner_local_facts(
    lifecycle_capture: LifecycleProbeCapture,
) -> None:
    attestation = lifecycle_capture.records[0].occurrence.lane_attestation_identity
    lane = lifecycle_capture.records[0].specification.lane_reference
    environment = lifecycle_capture.records[0].environment_identity
    assert eligible_failure_requirements(
        lifecycle_capture,
        source_commit=SOURCE_COMMIT,
        lockfile_hash=LOCKFILE_HASH,
        lane_identity=lane,
        attestation_identity=attestation,
        environment_identity=environment,
    ) == FAILURE_REQUIREMENTS

    for record in lifecycle_capture.records:
        assert len(record.case.identity) == 64
        assert len(record.specification.identity) == 64
        assert len(record.occurrence.identity) == 64
        assert record.occurrence.identity != record.occurrence.attempt_token
        assert record.specification.case_identity == record.case.identity
        assert record.occurrence.execution_specification_identity == record.specification.identity
        assert record.occurrence.case_identity == record.case.identity
        assert record.occurrence.actual_implementation_identity == record.specification.implementation.identity
        assert record.occurrence.input_member_identities == tuple(
            sorted(member.identity for member in record.case.inputs)
        )
        assert record.occurrence.lane_reference == record.specification.lane_reference

    raw = json.loads(lifecycle_capture.to_bytes())
    assert all("counts" not in record for record in raw["records"])
    assert all("counts" not in record["facts"] for record in raw["records"])
    assert LifecycleProbeCapture.from_bytes(lifecycle_capture.to_bytes()) == lifecycle_capture


def test_fixed_predicates_reconstruct_exact_typed_owner_facts(
    lifecycle_capture: LifecycleProbeCapture,
) -> None:
    records = {record.probe_id: record for record in lifecycle_capture.records}
    construction = records["construction-dependent-stop"]
    assert construction.runtime_facts["construction"] == {
        "failure_category": "DirectTrfConstructionError",
        "constructor_entries": ["first"],
        "executor_entries": [],
    }
    assert construction.method_step_outcome is None

    primary = records["primary-execution-failure"]
    assert primary.method_step_outcome is not None
    assert primary.method_step_outcome.primary_terminal == "execution_failure"
    assert primary.runtime_facts["primary"] == {
        "terminal": "execution_failure",
        "component_dispositions": ["not_started"],
    }

    commit = records["accepted-commit-rejected"]
    assert commit.runtime_facts["commit"] == {
        "terminal": "failed",
        "failure_category": "incompatible_state",
    }

    cancelled = records["cancelled-resampling-partial-publication"]
    assert cancelled.runtime_facts["resampling"] == {
        "terminal": "cancelled", "unstarted_ordinals": [0, 1],
        "completed_count": 0, "successful_count": 0, "failed_count": 0,
        "replicate_dispositions": ["not_started", "not_started"],
    }

    worker = records["worker-resampling-failure"]
    assert worker.runtime_facts["resampling"] == {
        "terminal": "completed", "unstarted_ordinals": [],
        "completed_count": 2, "successful_count": 1, "failed_count": 1,
        "replicate_dispositions": ["failed", "succeeded"],
    }

    interrupted = records["interrupted-resampling-partial-publication"]
    assert interrupted.runtime_facts["resampling"] == {
        "terminal": "interrupted", "unstarted_ordinals": [1],
        "completed_count": 1, "successful_count": 0, "failed_count": 0,
        "replicate_dispositions": ["interrupted", "not_started"],
    }


def test_stale_flat_counter_capture_and_altered_facts_are_ineligible(
    lifecycle_capture: LifecycleProbeCapture,
) -> None:
    retained = LifecycleProbeCapture.from_bytes(FIXTURE.read_bytes())
    current = migration_core.migration_core_current_release_selection()
    authority = migration_core.migration_core_authority_selection()
    assert retained.identity == current.lifecycle_probe_identity
    assert eligible_failure_requirements(
        retained,
        source_commit=current.lifecycle_source_commit,
        lockfile_hash=current.lifecycle_lockfile_hash,
        lane_identity=authority.lane_identity,
        attestation_identity=authority.attestation_identity,
        environment_identity=authority.environment_identity,
    ) == FAILURE_REQUIREMENTS

    raw = json.loads(lifecycle_capture.to_bytes())
    raw["records"][0].update(claims=list(FAILURE_REQUIREMENTS), passed=True)
    with pytest.raises(ValueError, match="Malformed lifecycle probe result"):
        LifecycleProbeCapture.from_bytes(json.dumps(raw).encode())

    facts = dict(lifecycle_capture.records[0].runtime_facts)
    facts["successor"] = {"denied": True, "downstream_entries": ["dependent"]}
    changed = dataclasses.replace(
        lifecycle_capture.records[0],
        facts=type(lifecycle_capture.records[0].facts).from_value(facts),
    )
    altered = LifecycleProbeCapture((changed, *lifecycle_capture.records[1:]))
    first = lifecycle_capture.records[0]
    eligible = eligible_failure_requirements(
        altered,
        source_commit=SOURCE_COMMIT,
        lockfile_hash=LOCKFILE_HASH,
        lane_identity=first.specification.lane_reference,
        attestation_identity=first.occurrence.lane_attestation_identity,
        environment_identity=first.environment_identity,
    )
    assert FAILURE_REQUIREMENTS[0] not in eligible
    assert FAILURE_REQUIREMENTS[5] not in eligible

    original = lifecycle_capture.records[0]
    foreign_case = dataclasses.replace(
        original.case,
        source_authority=CaseSourceAuthority("foreign-source", LOCKFILE_HASH),
    )
    foreign_specification = dataclasses.replace(
        original.specification, case_identity=foreign_case.identity
    )
    foreign_occurrence = dataclasses.replace(
        original.occurrence,
        execution_specification_identity=foreign_specification.identity,
        case_identity=foreign_case.identity,
    )
    foreign = dataclasses.replace(
        original,
        case=foreign_case,
        specification=foreign_specification,
        occurrence=foreign_occurrence,
    )
    replaced = LifecycleProbeCapture((foreign, *lifecycle_capture.records[1:]))
    assert eligible_failure_requirements(
        replaced,
        source_commit=SOURCE_COMMIT,
        lockfile_hash=LOCKFILE_HASH,
        lane_identity=original.specification.lane_reference,
        attestation_identity=original.occurrence.lane_attestation_identity,
        environment_identity=original.environment_identity,
    ) == ()


def test_current_selection_retains_exact_lifecycle_evidence() -> None:
    current = migration_core.migration_core_current_release_selection()
    assert current.selection_version == "migration-core-current-release-v3"
    assert current.lifecycle_probe_identity is not None
    assert (
        current.lifecycle_source_commit
        == "7a180c4c609dd8683838d41cf5c411ca9dfd3a9d"
    )
    result = migration_core.compile_current_phased_migration_core_status(
        Path(__file__).parents[1]
    )
    assert len(result.eligible_requirement_ids) == 30
    assert len(result.uncovered_requirement_ids) == 21
    assert result.compiler_status == "FAILED_CLOSED"
# fmt: on
