"""Truth-backed numerical baseline probes for migration qualification (#591)."""

from __future__ import annotations

import copy
import dataclasses
import math
import pickle
from collections.abc import Callable
from types import SimpleNamespace
from typing import cast

import numpy as np
import pytest

from chemex.numerical_lanes import (
    LiveLaneAuthority,
    RuntimeEnvironment,
    canonical_lanes,
)
from chemex.optimize.numerical_probes import (
    CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    AuthoritativeObjectiveAccounting,
    DeProbeEvidence,
    FiniteDifferenceProbeEvidence,
    GridOrderingProbeEvidence,
    GridProbeEvidence,
    NumericalProbeArtifact,
    NumericalProbeBaseline,
    NumericalProbeDefinition,
    NumericalProbeExecutionError,
    NumericalProbeQualificationError,
    SolverProbeEvidence,
    SpectralRiskProbeEvidence,
    numerical_probe_definitions,
    run_numerical_probe,
    run_numerical_probe_baseline,
)
from chemex.typing import Array


def _live_authority(
    monkeypatch: pytest.MonkeyPatch,
    lane_index: int = 0,
) -> LiveLaneAuthority:
    lane = canonical_lanes()[lane_index]
    environment = RuntimeEnvironment(lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: environment),
    )
    return lane.attest_current_process(lane.semantics.image_digest)


def test_probe_catalogue_freezes_required_truth_and_execution_contracts() -> None:
    definitions = numerical_probe_definitions()

    assert tuple(definition.probe_id for definition in definitions) == (
        "trf-routine-quadratic-v1",
        "trf-difficult-rosenbrock-v1",
        "grid-27-seed-coverage-v1",
        "grid-candidate-ordering-v1",
        "de-bounded-search-v1",
        "finite-difference-reliability-v1",
        "trf-active-bound-v1",
        "rank-deficient-linearization-v1",
        "ill-conditioned-linearization-v1",
    )
    assert {definition.category for definition in definitions} == {
        "TRF_ROUTINE",
        "TRF_DIFFICULT",
        "GRID",
        "GRID_ORDERING",
        "DE_SEARCH",
        "FINITE_DIFFERENCE",
        "BOUNDS",
        "RANK",
        "CONDITIONING",
    }
    assert all(definition.derivation.source_issue == 591 for definition in definitions)
    assert all(definition.derivation.source_revision for definition in definitions)
    assert all(definition.derivation.reduction for definition in definitions)
    assert all(definition.truth.reference for definition in definitions)
    assert all(definition.eligible_claims for definition in definitions)
    assert all(definition.failure_expectations for definition in definitions)
    assert all(definition.budget.to_record_value() for definition in definitions)
    assert all(definition.policy.to_record_value() for definition in definitions)
    assert all(
        NumericalProbeDefinition.from_record(definition.to_record()) == definition
        for definition in definitions
    )
    assert len({definition.identity for definition in definitions}) == len(definitions)

    altered = dataclasses.replace(
        definitions[0],
        eligible_claims=("altered-claim",),
    )
    with pytest.raises(ValueError, match="predeclared"):
        run_numerical_probe(altered)


@pytest.mark.parametrize(
    ("probe_id", "expected"),
    (
        ("trf-routine-quadratic-v1", (1.5, -0.5)),
        ("trf-difficult-rosenbrock-v1", (1.0, 1.0)),
    ),
)
def test_trf_probes_replay_with_authoritative_request_accounting(
    probe_id: str,
    expected: tuple[float, ...],
) -> None:
    definition = next(
        item for item in numerical_probe_definitions() if item.probe_id == probe_id
    )

    first = run_numerical_probe(definition)
    replay = run_numerical_probe(definition)

    assert first.terminal == "CONVERGED"
    assert first.risks == ()
    assert first.satisfied_claims == definition.eligible_claims
    assert isinstance(first.evidence, SolverProbeEvidence)
    assert first.evidence.accepted == pytest.approx(expected, rel=0.0, abs=1.0e-8)
    accounting = first.evidence.authoritative_accounting
    assert accounting.requests_received == accounting.requests_completed
    assert accounting.cache_hits + accounting.cache_misses == (
        accounting.requests_completed
    )
    assert accounting.workflow_requests == 1
    assert accounting.materialization_requests == 1
    assert accounting.requests_received == 2
    assert first.evidence.backend_diagnostics.diagnostic_evaluations > (
        accounting.requests_received
    )
    # SciPy 1.15 has no least_squares callback API, so its absence is explicit.
    assert first.evidence.backend_diagnostics.callback_calls == 0
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )
    assert first.identity == replay.identity
    assert NumericalProbeArtifact.from_record(first.to_record(), definition) == first


def test_finite_difference_probe_retains_truth_steps_and_request_fingerprint() -> None:
    definition = next(
        item
        for item in numerical_probe_definitions()
        if item.probe_id == "finite-difference-reliability-v1"
    )

    first = run_numerical_probe(definition)
    replay = run_numerical_probe(definition)

    assert first.terminal == "RELIABLE"
    assert isinstance(first.evidence, FiniteDifferenceProbeEvidence)
    assert first.evidence.estimate == pytest.approx(
        first.evidence.truth,
        rel=0.0,
        abs=first.evidence.absolute_tolerance,
    )
    expected_steps = tuple(
        (0.5 + nominal_step) - 0.5
        for nominal_step in (-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4)
    )
    assert first.evidence.actual_steps == expected_steps
    assert first.evidence.actual_steps[3] != 1.0e-4
    assert first.evidence.estimate == math.fsum(
        weight * value
        for weight, value in zip(
            first.evidence.weights,
            first.evidence.sampled_values,
            strict=True,
        )
    )
    assert first.evidence.authoritative_accounting.requests_received == 5
    assert first.evidence.backend_diagnostics.nfev == 0
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )


def test_bound_rank_and_conditioning_probes_fail_closed_against_exact_truth() -> None:
    definitions = numerical_probe_definitions()
    selected = {
        definition.probe_id: definition
        for definition in definitions
        if definition.probe_id
        in {
            "trf-active-bound-v1",
            "rank-deficient-linearization-v1",
            "ill-conditioned-linearization-v1",
        }
    }

    boundary = run_numerical_probe(selected["trf-active-bound-v1"])
    rank = run_numerical_probe(selected["rank-deficient-linearization-v1"])
    conditioning = run_numerical_probe(selected["ill-conditioned-linearization-v1"])

    assert isinstance(boundary.evidence, SolverProbeEvidence)
    assert boundary.evidence.accepted == pytest.approx((1.0,), rel=0.0, abs=1.0e-8)
    assert boundary.evidence.active_mask == (1,)
    assert boundary.risks == ("ACTIVE_BOUND",)

    assert isinstance(rank.evidence, SpectralRiskProbeEvidence)
    assert rank.terminal == "RISK_DETECTED"
    assert rank.evidence.rank == 1
    assert rank.evidence.dimension == 2
    assert rank.risks == ("RANK_DEFICIENT",)

    assert isinstance(conditioning.evidence, SpectralRiskProbeEvidence)
    assert conditioning.terminal == "RISK_DETECTED"
    assert conditioning.evidence.rank == 2
    assert conditioning.evidence.condition == pytest.approx(
        1.0e8,
        rel=1.0e-12,
        abs=0.0,
    )
    assert conditioning.evidence.condition_limit == pytest.approx(
        1.0e6,
        rel=0.0,
        abs=0.0,
    )
    assert conditioning.risks == ("ILL_CONDITIONED",)

    for definition, artifact in zip(
        (
            selected["trf-active-bound-v1"],
            selected["rank-deficient-linearization-v1"],
            selected["ill-conditioned-linearization-v1"],
        ),
        (boundary, rank, conditioning),
        strict=True,
    ):
        assert (
            NumericalProbeArtifact.from_record(artifact.to_record(), definition)
            == artifact
        )


def test_grid_probe_retains_all_physical_seeds_and_observed_candidate_order() -> None:
    definition = numerical_probe_definitions()[2]

    first = run_numerical_probe(definition)
    replay = run_numerical_probe(definition)

    assert first.terminal == "SELECTED"
    assert isinstance(first.evidence, GridProbeEvidence)
    assert tuple(seed.ordinal for seed in first.evidence.seeds) == tuple(range(27))
    assert len({seed.seed for seed in first.evidence.seeds}) == 27
    assert first.evidence.ordered_seed_ordinals == tuple(
        seed.ordinal
        for seed in sorted(
            first.evidence.seeds,
            key=lambda seed: (
                seed.solver.chi_square,
                seed.solver.accepted,
                seed.ordinal,
            ),
        )
    )
    assert (
        first.evidence.selected_seed_ordinal
        == (first.evidence.ordered_seed_ordinals[0])
    )
    selected = first.evidence.seeds[first.evidence.selected_seed_ordinal]
    assert selected.solver.accepted == pytest.approx(
        (1.0, -1.0, 0.5), rel=0.0, abs=1.0e-8
    )
    assert first.evidence.objective_accounting.requests_received == sum(
        seed.solver.authoritative_accounting.requests_received
        for seed in first.evidence.seeds
    )
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )
    assert first.identity == replay.identity


def test_grid_ordering_claim_uses_predeclared_truth_not_observed_candidates() -> None:
    definition = next(
        item
        for item in numerical_probe_definitions()
        if item.probe_id == "grid-candidate-ordering-v1"
    )

    artifact = run_numerical_probe(definition)

    assert artifact.terminal == "ORDERED"
    assert isinstance(artifact.evidence, GridOrderingProbeEvidence)
    assert artifact.evidence.expected_order == (1, 3, 2, 0)
    assert artifact.evidence.observed_order == artifact.evidence.expected_order
    for wrong_order in ((1, 0, 2, 3), (3, 1, 2, 0)):
        wrong_evidence = dataclasses.replace(
            artifact.evidence,
            observed_order=wrong_order,
            selected_ordinal=wrong_order[0],
        )
        wrong = dataclasses.replace(artifact, evidence=wrong_evidence)
        with pytest.raises(ValueError, match="ordering"):
            wrong.require_qualification(definition)

    reordered = copy.deepcopy(artifact.to_record())
    raw_evidence = cast("dict[str, object]", reordered["evidence"])
    candidates = cast("list[object]", raw_evidence["candidates"])
    candidates[0], candidates[1] = candidates[1], candidates[0]
    with pytest.raises(ValueError, match="ordinal order"):
        NumericalProbeArtifact.from_record(reordered, definition)


def test_de_probe_retains_seeded_population_order_and_separate_trf_polish() -> None:
    definition = next(
        item
        for item in numerical_probe_definitions()
        if item.probe_id == "de-bounded-search-v1"
    )

    first = run_numerical_probe(definition)
    replay = run_numerical_probe(definition)

    assert first.terminal == "POLISHED"
    assert isinstance(first.evidence, DeProbeEvidence)
    assert first.evidence.root_seed == 591
    assert len(first.evidence.final_population) == 5
    assert first.evidence.ordered_population_indices == tuple(
        candidate.population_index
        for candidate in sorted(
            first.evidence.final_population,
            key=lambda candidate: (candidate.objective, candidate.population_index),
        )
    )
    assert (
        first.evidence.selected_population_index
        == (first.evidence.ordered_population_indices[0])
    )
    assert first.evidence.authoritative_search_accounting.requests_received > 0
    assert (
        first.evidence.authoritative_search_accounting.materialization_requests
        == len(first.evidence.final_population)
    )
    assert first.evidence.search_diagnostics.nfev > 0
    assert first.evidence.search_diagnostics.callback_calls == 6
    assert first.evidence.search_diagnostics.iterations == 6
    for candidate in first.evidence.final_population:
        displacement = candidate.vector[0] - 0.25
        expected_objective = displacement**2 + 0.25 * math.sin(5.0 * displacement) ** 2
        assert candidate.objective == pytest.approx(
            expected_objective,
            rel=0.0,
            abs=1.0e-16,
        )
    assert first.evidence.polish.accepted == pytest.approx((0.25,), rel=0.0, abs=1.0e-8)
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )
    assert first.identity == replay.identity


def test_backend_success_cannot_mint_truth_claims(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = numerical_probe_definitions()[0]

    def false_success(*_args: object, **_kwargs: object) -> object:
        return SimpleNamespace(
            x=np.asarray((0.0, 0.0)),
            active_mask=np.asarray((0, 0)),
            nfev=1,
            njev=0,
            success=True,
        )

    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        false_success,
    )

    artifact = run_numerical_probe(definition)

    assert artifact.terminal == "TRUTH_MISMATCH"
    assert artifact.satisfied_claims == ()
    with pytest.raises(ValueError, match="qualify"):
        NumericalProbeBaseline((definition,), (artifact,))


def test_backend_call_pattern_does_not_change_authoritative_request_accounting(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = numerical_probe_definitions()[0]

    def backend(call_count: int) -> Callable[..., object]:
        def solve(
            fun: Callable[[Array], Array],
            start: object,
            **_kwargs: object,
        ) -> object:
            for _index in range(call_count):
                fun(np.asarray(start))
            return SimpleNamespace(
                x=np.asarray((1.5, -0.5)),
                active_mask=np.asarray((0, 0)),
                nfev=call_count,
                njev=1,
                success=True,
            )

        return solve

    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        backend(1),
    )
    first = run_numerical_probe(definition)
    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        backend(7),
    )
    second = run_numerical_probe(definition)

    assert isinstance(first.evidence, SolverProbeEvidence)
    assert isinstance(second.evidence, SolverProbeEvidence)
    assert first.evidence.authoritative_accounting == (
        second.evidence.authoritative_accounting
    )
    assert first.evidence.trajectory_fingerprint == (
        second.evidence.trajectory_fingerprint
    )
    assert first.evidence.backend_diagnostics.nfev == 1
    assert second.evidence.backend_diagnostics.nfev == 7
    assert first.evidence.backend_diagnostics.diagnostic_evaluations == 1
    assert second.evidence.backend_diagnostics.diagnostic_evaluations == 7


def test_invalid_authoritative_request_has_typed_reconciled_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = numerical_probe_definitions()[0]

    def invalid_result(*_args: object, **_kwargs: object) -> object:
        return SimpleNamespace(
            x=np.asarray((float("nan"), 0.0)),
            active_mask=np.asarray((0, 0)),
            nfev=1,
            njev=0,
            success=True,
        )

    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        invalid_result,
    )

    with pytest.raises(NumericalProbeExecutionError) as captured:
        run_numerical_probe(definition)

    error = captured.value
    assert error.terminal == "REQUEST_FAILED"
    assert error.failure_kind == "INVALID_INPUT"
    assert error.accounting.requests_received == 2
    assert error.accounting.requests_failed == 1
    assert error.accounting.requests_completed == 1
    assert error.accounting.requests_refused == 0


def test_evaluator_failure_after_issued_request_is_typed_and_reconciled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = next(
        item
        for item in numerical_probe_definitions()
        if item.probe_id == "finite-difference-reliability-v1"
    )

    def evaluator_failure(_value: float) -> float:
        raise ArithmeticError

    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.math.exp",
        evaluator_failure,
    )

    with pytest.raises(NumericalProbeExecutionError) as captured:
        run_numerical_probe(definition)

    error = captured.value
    assert error.terminal == "REQUEST_FAILED"
    assert error.failure_kind == "EVALUATION_FAILURE"
    assert error.accounting.requests_received == 1
    assert error.accounting.requests_failed == 1
    assert error.accounting.requests_completed == 0


def test_nonfinite_authoritative_result_is_typed_and_never_qualifies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = next(
        item
        for item in numerical_probe_definitions()
        if item.probe_id == "finite-difference-reliability-v1"
    )
    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.math.exp",
        lambda _value: float("inf"),
    )

    with pytest.raises(NumericalProbeExecutionError) as captured:
        run_numerical_probe(definition)

    assert captured.value.terminal == "REQUEST_FAILED"
    assert captured.value.failure_kind == "NON_FINITE_RESULT"
    assert captured.value.accounting.requests_received == 1
    assert captured.value.accounting.requests_failed == 1


def test_authoritative_request_outcomes_form_a_closed_accounting_vocabulary() -> None:
    refused = AuthoritativeObjectiveAccounting(1, 0, 0, 0, 1, ())
    failed = AuthoritativeObjectiveAccounting(
        0,
        1,
        0,
        0,
        0,
        ("EVALUATION_FAILURE",),
    )

    assert refused.requests_received == refused.requests_refused == 1
    assert failed.requests_received == failed.requests_failed == 1
    assert AuthoritativeObjectiveAccounting.from_record(refused.to_record()) == refused
    with pytest.raises(ValueError, match="reconcile"):
        AuthoritativeObjectiveAccounting(1, 0, 0, 0, 0, ())


def test_probe_baseline_rejects_evidence_from_another_category() -> None:
    baseline = run_numerical_probe_baseline()
    wrong = dataclasses.replace(
        baseline.artifacts[0],
        evidence=baseline.artifacts[3].evidence,
    )

    with pytest.raises(TypeError, match="category"):
        NumericalProbeBaseline(
            baseline.definitions,
            (wrong, *baseline.artifacts[1:]),
        )


def test_probe_baseline_manifest_replays_round_trips_and_rejects_tampering() -> None:
    first = run_numerical_probe_baseline()
    replay = run_numerical_probe_baseline()

    assert len(first.definitions) == 9
    assert tuple(artifact.probe_id for artifact in first.artifacts) == tuple(
        definition.probe_id for definition in first.definitions
    )
    assert first.manifest_identity == replay.manifest_identity
    assert first.identity == replay.identity
    assert first.historical_qualification == "CAPTURE_ONLY"
    assert first.historically_satisfied_claims == ()
    assert NumericalProbeBaseline.from_record(first.to_record()) == first
    with pytest.raises(TypeError, match="live authority"):
        run_numerical_probe_baseline(expected_manifest_identity=first.manifest_identity)

    tampered = copy.deepcopy(first.to_record())
    artifacts = tampered["artifacts"]
    assert isinstance(artifacts, list)
    artifact = cast("dict[str, object]", artifacts[0])
    evidence = cast("dict[str, object]", artifact["evidence"])
    accounting = cast("dict[str, object]", evidence["objective_accounting"])
    accounting["fresh_evaluations"] = cast("int", accounting["fresh_evaluations"]) + 1

    with pytest.raises(ValueError, match="accounting|payload|reconcile"):
        NumericalProbeBaseline.from_record(tampered)


def test_serialized_backend_grid_trajectory_and_manifest_tampering_is_rejected() -> (
    None
):
    baseline = run_numerical_probe_baseline()
    records: list[dict[str, object]] = []

    backend = copy.deepcopy(baseline.to_record())
    backend_artifacts = cast("list[dict[str, object]]", backend["artifacts"])
    backend_evidence = cast("dict[str, object]", backend_artifacts[0]["evidence"])
    diagnostics = cast("dict[str, object]", backend_evidence["backend_diagnostics"])
    diagnostics["nfev"] = cast("int", diagnostics["nfev"]) + 1
    records.append(backend)

    trajectory = copy.deepcopy(baseline.to_record())
    trajectory_artifacts = cast("list[dict[str, object]]", trajectory["artifacts"])
    trajectory_evidence = cast("dict[str, object]", trajectory_artifacts[0]["evidence"])
    trajectory_evidence["trajectory_fingerprint"] = "0" * 64
    records.append(trajectory)

    grid = copy.deepcopy(baseline.to_record())
    grid_artifacts = cast("list[dict[str, object]]", grid["artifacts"])
    grid_index = next(
        index
        for index, definition in enumerate(baseline.definitions)
        if definition.probe_id == "grid-candidate-ordering-v1"
    )
    grid_evidence = cast("dict[str, object]", grid_artifacts[grid_index]["evidence"])
    grid_evidence["expected_order"] = [1, 2, 3, 0]
    records.append(grid)

    manifest = copy.deepcopy(baseline.to_record())
    manifest["manifest_identity"] = "0" * 64
    records.append(manifest)

    for record in records:
        with pytest.raises(ValueError):
            NumericalProbeBaseline.from_record(record)


def test_matching_external_manifest_is_recorded_but_remains_evidence_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capture = run_numerical_probe_baseline()
    assert capture.manifest_identity == CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
    authority = _live_authority(monkeypatch)

    authority_without_reference = run_numerical_probe_baseline(authority)
    assert authority_without_reference.historical_qualification == "CAPTURE_ONLY"
    assert authority_without_reference.historically_satisfied_claims == ()

    matched = run_numerical_probe_baseline(
        authority,
        expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    )

    assert matched.historical_qualification == "REFERENCE_MATCHED"
    assert matched.historically_satisfied_claims == ("trajectory-manifest-compatible",)
    tampered_role = copy.deepcopy(matched.to_record())
    tampered_role["observed_lane_role"] = "PYTHON_COMPATIBILITY"
    with pytest.raises(ValueError, match="identity"):
        NumericalProbeBaseline.from_record(tampered_role)
    with pytest.raises(NumericalProbeQualificationError) as mismatch:
        capture.require_live_qualification(
            authority,
            expected_manifest_identity="0" * 64,
            required_lane_role="CANONICAL_NUMERICAL",
        )
    assert mismatch.value.terminal == "EXPECTED_REFERENCE_MISMATCH"


def test_live_run_cannot_mint_a_new_manifest_from_backend_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capture = run_numerical_probe_baseline()
    first_artifact = capture.artifacts[0]
    assert isinstance(first_artifact.evidence, SolverProbeEvidence)
    diagnostics = dataclasses.replace(
        first_artifact.evidence.backend_diagnostics,
        nfev=first_artifact.evidence.backend_diagnostics.nfev + 1,
    )
    drifted_evidence = dataclasses.replace(
        first_artifact.evidence,
        backend_diagnostics=diagnostics,
    )
    drifted_artifact = dataclasses.replace(
        first_artifact,
        evidence=drifted_evidence,
    )
    drifted = NumericalProbeBaseline(
        capture.definitions,
        (drifted_artifact, *capture.artifacts[1:]),
    )
    authority = _live_authority(monkeypatch)

    assert drifted.manifest_identity != CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
    with pytest.raises(NumericalProbeQualificationError) as rejected:
        drifted.require_live_qualification(
            authority,
            expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
            required_lane_role="CANONICAL_NUMERICAL",
        )

    assert rejected.value.terminal == "MANIFEST_MISMATCH"


def test_fabricated_serialized_lane_facts_cannot_gain_live_authority(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capture = run_numerical_probe_baseline()
    fabricated = NumericalProbeBaseline(
        capture.definitions,
        capture.artifacts,
        "a" * 64,
        "CANONICAL_NUMERICAL",
        "b" * 64,
        "c" * 64,
        CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
        "REFERENCE_MATCHED",
    )
    restored = NumericalProbeBaseline.from_record(fabricated.to_record())
    authority = _live_authority(monkeypatch)

    with pytest.raises(NumericalProbeQualificationError) as rejected:
        restored.require_live_qualification(
            authority,
            expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
            required_lane_role="CANONICAL_NUMERICAL",
        )

    assert rejected.value.terminal == "LANE_MISMATCH"


def test_foreign_live_authority_does_not_match_historical_lane_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    canonical = _live_authority(monkeypatch)
    matched = run_numerical_probe_baseline(
        canonical,
        expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    )
    foreign_lane = dataclasses.replace(
        canonical_lanes()[0],
        name="foreign-canonical-numerical-lane",
    )
    foreign_environment = RuntimeEnvironment(foreign_lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(
            lambda cls, image_digest, provenance_path=None: foreign_environment
        ),
    )
    foreign = foreign_lane.attest_current_process(foreign_lane.semantics.image_digest)

    with pytest.raises(NumericalProbeQualificationError) as rejected:
        matched.require_live_qualification(
            foreign,
            expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
            required_lane_role="CANONICAL_NUMERICAL",
        )

    assert rejected.value.terminal == "LANE_MISMATCH"


def test_serialized_probe_baseline_never_recreates_live_lane_authority(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    authority = _live_authority(monkeypatch)
    baseline = run_numerical_probe_baseline(
        authority,
        expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    )
    restored = pickle.loads(pickle.dumps(baseline))  # noqa: S301 - trusted self-record

    assert restored.historical_qualification == "REFERENCE_MATCHED"
    assert restored.observed_lane_identity == authority.lane_identity
    assert not isinstance(restored, LiveLaneAuthority)
    restored.require_live_qualification(
        authority,
        expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
        required_lane_role="CANONICAL_NUMERICAL",
    )

    with pytest.raises(TypeError, match="live"):
        restored.require_live_qualification(  # type: ignore[arg-type]
            restored,
            expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
            required_lane_role="CANONICAL_NUMERICAL",
        )
    with pytest.raises(TypeError):
        copy.copy(authority)
    with pytest.raises(TypeError):
        pickle.dumps(authority)


def test_foreign_or_compatibility_authority_cannot_qualify_canonical_capture(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capture = run_numerical_probe_baseline()
    compatibility = _live_authority(monkeypatch, 1)

    with pytest.raises(NumericalProbeQualificationError) as rejected:
        capture.require_live_qualification(
            compatibility,
            expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
            required_lane_role="CANONICAL_NUMERICAL",
        )

    assert rejected.value.terminal == "LANE_ROLE_MISMATCH"
