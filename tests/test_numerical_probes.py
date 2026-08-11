"""Truth-backed numerical baseline probes for migration qualification (#591)."""

from __future__ import annotations

import copy
import dataclasses
import math
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
    DeProbeEvidence,
    FiniteDifferenceProbeEvidence,
    GridProbeEvidence,
    NumericalProbeArtifact,
    NumericalProbeBaseline,
    NumericalProbeDefinition,
    NumericalProbeExecutionError,
    SolverProbeEvidence,
    SpectralRiskProbeEvidence,
    numerical_probe_definitions,
    run_numerical_probe,
    run_numerical_probe_baseline,
)
from chemex.typing import Array


def _live_authority(monkeypatch: pytest.MonkeyPatch) -> LiveLaneAuthority:
    lane = canonical_lanes()[0]
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
        "grid-27-seed-ordering-v1",
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
    accounting = first.evidence.objective_accounting
    assert accounting.requests_received == accounting.requests_completed
    assert accounting.cache_hits + accounting.cache_misses == (
        accounting.requests_completed
    )
    assert accounting.materialization_requests == 1
    assert accounting.cache_hits >= 1
    assert accounting.workflow_requests > first.evidence.backend_diagnostics.nfev
    # SciPy 1.15 has no least_squares callback API, so its absence is explicit.
    assert first.evidence.backend_diagnostics.callback_calls == 0
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )
    assert first.identity == replay.identity
    assert NumericalProbeArtifact.from_record(first.to_record(), definition) == first


def test_finite_difference_probe_retains_truth_steps_and_request_fingerprint() -> None:
    definition = numerical_probe_definitions()[4]

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
    assert first.evidence.objective_accounting.requests_received == 5
    assert first.evidence.backend_diagnostics.nfev == 0
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )


def test_bound_rank_and_conditioning_probes_fail_closed_against_exact_truth() -> None:
    definitions = numerical_probe_definitions()

    boundary = run_numerical_probe(definitions[5])
    rank = run_numerical_probe(definitions[6])
    conditioning = run_numerical_probe(definitions[7])

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
        definitions[5:],
        (boundary, rank, conditioning),
        strict=True,
    ):
        assert (
            NumericalProbeArtifact.from_record(artifact.to_record(), definition)
            == artifact
        )


def test_grid_probe_retains_all_physical_seeds_and_canonical_candidate_order() -> None:
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
            key=lambda seed: (seed.solver.chi_square, seed.ordinal),
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
        seed.solver.objective_accounting.requests_received
        for seed in first.evidence.seeds
    )
    assert first.evidence.trajectory_fingerprint == (
        replay.evidence.trajectory_fingerprint
    )
    assert first.identity == replay.identity


def test_de_probe_retains_seeded_population_order_and_separate_trf_polish() -> None:
    definition = numerical_probe_definitions()[3]

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
    assert first.evidence.search_accounting.requests_received > 0
    assert first.evidence.search_accounting.materialization_requests == len(
        first.evidence.final_population
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


def test_budget_exhaustion_retains_refused_request_and_trajectory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = numerical_probe_definitions()[0]

    def exhauster(refused: tuple[float, float]) -> Callable[..., object]:
        def exhaust_budget(
            fun: Callable[[Array], Array],
            _start: object,
            **_kwargs: object,
        ) -> object:
            for _index in range(64):
                fun(np.asarray((0.0, 0.0)))
            fun(np.asarray(refused))
            raise AssertionError("Budget should terminate before backend completion")

        return exhaust_budget

    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        exhauster((1.0, 1.0)),
    )

    with pytest.raises(NumericalProbeExecutionError) as first:
        run_numerical_probe(definition)
    monkeypatch.setattr(
        "chemex.optimize.numerical_probes.least_squares",
        exhauster((2.0, 2.0)),
    )
    with pytest.raises(NumericalProbeExecutionError) as second:
        run_numerical_probe(definition)

    assert first.value.terminal == "BUDGET_EXHAUSTED"
    assert first.value.accounting.requests_completed == 64
    assert first.value.accounting.requests_refused == 1
    assert len(first.value.trajectory_fingerprint) == 64
    assert first.value.trajectory_fingerprint != second.value.trajectory_fingerprint


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

    assert len(first.definitions) == 8
    assert tuple(artifact.probe_id for artifact in first.artifacts) == tuple(
        definition.probe_id for definition in first.definitions
    )
    assert first.manifest_identity == replay.manifest_identity
    assert first.identity == replay.identity
    assert first.authority_role == "UNQUALIFIED_DIAGNOSTIC"
    assert first.lane_reference.startswith("unqualified-")
    assert NumericalProbeBaseline.from_record(first.to_record()) == first
    assert (
        run_numerical_probe_baseline(expected_manifest_identity=first.manifest_identity)
        == first
    )
    with pytest.raises(ValueError, match="differs"):
        run_numerical_probe_baseline(expected_manifest_identity="0" * 64)

    tampered = copy.deepcopy(first.to_record())
    artifacts = tampered["artifacts"]
    assert isinstance(artifacts, list)
    artifact = cast("dict[str, object]", artifacts[0])
    evidence = cast("dict[str, object]", artifact["evidence"])
    accounting = cast("dict[str, object]", evidence["objective_accounting"])
    accounting["cache_hits"] = cast("int", accounting["cache_hits"]) + 1

    with pytest.raises(ValueError, match="accounting|payload"):
        NumericalProbeBaseline.from_record(tampered)


def test_probe_baseline_retains_live_lane_authority(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    authority = _live_authority(monkeypatch)

    baseline = run_numerical_probe_baseline(authority)

    assert baseline.authority_role == "QUALIFIED_LANE"
    assert baseline.lane_reference == authority.lane_identity
    assert baseline.lane_attestation_identity == authority.identity
    assert baseline.environment_identity == authority.environment_identity
    assert NumericalProbeBaseline.from_record(baseline.to_record()) == baseline
