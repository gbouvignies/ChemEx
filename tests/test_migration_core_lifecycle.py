# ruff: noqa: I001

import dataclasses
import json
from pathlib import Path

# fmt: off
import pytest

import chemex.optimize.method_step as step
import chemex.optimize.native_resampling as resampling
from chemex.migration_core_lifecycle import (
    FAILURE_REQUIREMENTS,
    LifecycleOperationCounts,
    LifecycleProbeCapture,
    LifecycleProbeResult,
    eligible_failure_requirements,
)
from chemex.optimize.direct_trf import FitCommitTerminal
from tests.qualification.capture_migration_core_lifecycle import (
    ProbeObservation,
    observe_lifecycle_probes,
)
FIXTURE = Path(__file__).parent / "fixtures/migration_core_lifecycle_probes_v1.json"
def _project(observed: ProbeObservation) -> dict[str, object]:  # noqa: C901
    if observed.construction_failure:
        return {"stage": "construction", "terminal": "failed", "counts": LifecycleOperationCounts(1, 0, 1, 0, 0, 3), "acceptance": "not_started", "commit": "not_started", "publication": "suppressed", "successor": "denied", "partial_evidence": "not_applicable"}
    outcome = observed.outcome
    assert outcome is not None
    derivation = outcome.derivations[0] if outcome.derivations else None
    evidence = next((item for item in derivation.artifacts if isinstance(item, resampling.ResamplingEvidence)), None) if derivation else None
    if outcome.publication_failure:
        stage = "publication"
    elif derivation and (derivation.disposition is step.DerivationDisposition.INTERRUPTED or evidence and evidence.failed_count):
        stage = "worker"
    elif derivation:
        stage = "requested_evidence"
    elif outcome.commit_operation:
        stage = "commit"
    elif outcome.primary_terminal == "accepted":
        stage = "acceptance"
    elif outcome.primary_terminal == "decomposition_validation_failure":
        stage = "materialization"
    else:
        stage = "primary_execution"
    acceptance = "accepted" if outcome.accepted_result else "failed" if stage == "acceptance" else "not_started"
    commit = "not_started" if outcome.commit_operation is None else "committed" if outcome.commit_operation.terminal is FitCommitTerminal.COMMITTED else "rejected"
    publication = "failed" if outcome.publication_failure else "suppressed" if outcome.publication is None else "diagnostics_only" if outcome.lifecycle is step.MethodStepLifecycle.FAILED else "partial_only"
    partial = "not_applicable"
    if derivation and evidence:
        states = tuple(item.disposition.value for item in evidence.outcomes)
        succeeded = states.count("succeeded")
        if derivation.disposition is step.DerivationDisposition.CANCELLED:
            partial, counts = "cancelled_population_retained", LifecycleOperationCounts(3, 2, 0, 1, 0, states.count("not_started"))
        elif derivation.disposition is step.DerivationDisposition.INTERRUPTED:
            partial, counts = "interrupted_population_retained", LifecycleOperationCounts(2 + len(states) - states.count("not_started"), 2 + succeeded, 0, 0, states.count("interrupted"), states.count("not_started"))
        elif evidence.failed_count:
            partial, counts = "worker_failure_retained", LifecycleOperationCounts(2 + len(states), 2 + succeeded, evidence.failed_count, 0, 0, 0)
        else:
            assert any(isinstance(item, resampling.SummaryFailure) for item in derivation.artifacts)
            partial, counts = "summary_failure_retained", LifecycleOperationCounts(3 + len(states), 2 + succeeded, 1, 0, 0, 0)
    else:
        counts = {
            "primary_execution": LifecycleOperationCounts(1, 0, 1, 0, 0, 3),
            "materialization": LifecycleOperationCounts(2, 1, 1, 0, 0, 2),
            "acceptance": LifecycleOperationCounts(2, 1, 1, 0, 0, 2),
            "commit": LifecycleOperationCounts(3, 2, 1, 0, 0, 1),
            "publication": LifecycleOperationCounts(3, 2, 1, 0, 0, 0),
        }[stage]
    return {"stage": stage, "terminal": outcome.lifecycle.value, "counts": counts, "acceptance": acceptance, "commit": commit, "publication": publication, "successor": "available" if observed.successor_snapshot else "denied", "partial_evidence": partial}
def _capture(observations: tuple[ProbeObservation, ...]) -> LifecycleProbeCapture:
    common = {"source_commit": "0" * 40, "lockfile_hash": "1" * 64, "lane_identity": "2" * 64, "attestation_identity": "3" * 64, "environment_identity": "4" * 64, "case_identity": "5" * 64, "execution_specification_identity": "6" * 64, "occurrence_identity": "7" * 64}
    return LifecycleProbeCapture(tuple(LifecycleProbeResult(item.probe_id, starting_revision=item.starting_snapshot.revision, ending_revision=item.ending_snapshot.revision, **common, **_project(item)) for item in observations))
def test_capture_reruns_real_behavior_and_fixed_validator(tmp_path: Path) -> None:
    if not FIXTURE.exists():
        pytest.skip("canonical lifecycle capture is selected only in evidence commit")
    frozen = LifecycleProbeCapture.from_bytes(FIXTURE.read_bytes())
    observed = observe_lifecycle_probes(tmp_path)
    names = ("stage", "terminal", "counts", "acceptance", "commit", "publication", "successor", "partial_evidence")
    assert [tuple(_project(item)[name] for name in names) for item in observed] == [tuple(getattr(item, name) for name in names) for item in frozen.records]
    first = frozen.records[0]
    assert eligible_failure_requirements(frozen, source_commit=first.source_commit, lockfile_hash=first.lockfile_hash, lane_identity=first.lane_identity, attestation_identity=first.attestation_identity, environment_identity=first.environment_identity) == FAILURE_REQUIREMENTS
def test_claims_and_wrong_facts_do_not_grant_eligibility(tmp_path: Path) -> None:
    capture = _capture(observe_lifecycle_probes(tmp_path))
    raw = json.loads(capture.to_bytes())
    raw["records"][0].update(claims=list(FAILURE_REQUIREMENTS), passed=True)
    with pytest.raises(ValueError, match="Malformed lifecycle probe result"):
        LifecycleProbeCapture.from_bytes(json.dumps(raw).encode())
    altered = LifecycleProbeCapture((dataclasses.replace(capture.records[0], successor="available"), *capture.records[1:]))
    eligible = eligible_failure_requirements(altered, source_commit="0" * 40, lockfile_hash="1" * 64, lane_identity="2" * 64, attestation_identity="3" * 64, environment_identity="4" * 64)
    assert FAILURE_REQUIREMENTS[0] not in eligible
    assert FAILURE_REQUIREMENTS[5] not in eligible
# fmt: on
