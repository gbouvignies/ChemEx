"""Behavioral qualification tests for native resampling evidence (#599).

The public seams are deterministic replicate planning and data generation,
evidence-only projected optimization, and scope-atomic summary derivation.  The
production lmfit runner remains outside this qualification harness.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable

import numpy as np
import pytest

from chemex.optimize.direct_trf import AcceptedFitResult
from chemex.optimize.native_resampling import (
    ClaimState,
    OptimizationStrategy,
    ProjectedOptimizationFailure,
    ProjectedOptimizationResult,
    ProjectedOptimizationSuccess,
    ReplicateDisposition,
    ReplicateRequest,
    ResamplingConstructionError,
    ResamplingDraw,
    ResamplingLifecycle,
    ResamplingPlan,
    ResamplingPopulation,
    ResamplingScheme,
    SummaryTerminal,
    execute_resampling_evidence,
    generate_resampling_draw,
    summarize_resampling_evidence,
)
from chemex.runtime import ExecutionSettings


def _accepted_anchor() -> AcceptedFitResult:
    """Use the exact authoritative constructor without fabricating fit authority."""
    from chemex.evaluation.native import EvaluationResult, ResolvedEvaluationValues

    result = EvaluationResult(
        plan_identity="accepted-plan",
        parameterization_identity="accepted-parameterization",
        evaluator_compatibility_identity="accepted-evaluator",
        resolved_values=ResolvedEvaluationValues(
            "accepted-parameterization",
            "accepted-program",
            (("A", 1.0), ("B", 2.0)),
        ),
        unscaled_calculations=np.array([1.0]),
        normalized_calculations=np.array([1.0]),
        residuals=np.array([0.0]),
        profiles=(),
    )
    return AcceptedFitResult.for_qualification(
        occurrence_identity="accepted-occurrence",
        problem_identity="accepted-problem",
        invocation_identity="accepted-invocation",
        execution_identity="accepted-execution",
        materialization_identity="accepted-materialization",
        parameterization_identity="accepted-parameterization",
        evaluator_parameterization_identity="accepted-parameterization",
        source_occurrence_identity="source-occurrence",
        source_revision=4,
        controlled_ids=("A", "B"),
        vector=(1.0, 2.0),
        chi_square=0.0,
        evaluation_result=result,
        commit_scope=("A", "B"),
        commit_items=(("A", 1.0), ("B", 2.0)),
        origin_context_identity="accepted-origin",
    )


def _population() -> ResamplingPopulation:
    return ResamplingPopulation(
        source_dataset_identity="qualification-dataset",
        observed=(10.0, 20.0, 30.0, 40.0, 50.0),
        calculated=(11.0, 19.0, 31.0, 39.0, 49.0),
        standard_errors=(1.0, 2.0, 1.5, 0.5, 2.5),
        mask=(True, True, True, True, False),
        references=(True, True, False, False, False),
        nucleus_groups=("N1", "N1", "N2", "N3", "N3"),
    )


def _plan(
    accepted: AcceptedFitResult,
    scheme: ResamplingScheme,
    *,
    count: int = 4,
) -> ResamplingPlan:
    return ResamplingPlan.for_accepted(
        accepted,
        scheme=scheme,
        replicate_count=count,
        root_seed=0x1234_5678_90AB_CDEF,
        source_dataset_identity="qualification-dataset",
        output_scope=("A", "B"),
        optimization_projection_identity="direct-trf-projection-v1",
        minimum_successful_count=2,
    )


def _successful_projection(
    offset: float = 0.0,
) -> Callable[[ReplicateRequest, ResamplingDraw], ProjectedOptimizationSuccess]:
    def execute(
        request: ReplicateRequest,
        draw: ResamplingDraw,
    ) -> ProjectedOptimizationSuccess:
        ordinal = request.ordinal
        return ProjectedOptimizationSuccess(
            transformed_data_identity=draw.identity,
            evaluation_plan_identity=f"plan-{ordinal}",
            problem_identity=f"problem-{ordinal}",
            invocation_identity=f"invocation-{ordinal}",
            execution_identity=f"execution-{ordinal}",
            strategy=OptimizationStrategy.DIRECT_TRF,
            component_outcome_identities=(f"component-{ordinal}",),
            resolved_items=(("A", offset + ordinal + 1.0), ("B", 10.0 - ordinal)),
            chi_square=float(ordinal + 1),
        )

    return execute


def test_replicate_plan_uses_canonical_ordinals_and_identity_derived_seeds() -> None:
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.MONTE_CARLO)

    assert tuple(request.ordinal for request in plan.replicates) == (0, 1, 2, 3)
    assert len({request.identity for request in plan.replicates}) == 4
    assert len({request.seed for request in plan.replicates}) == 4
    assert plan == _plan(accepted, ResamplingScheme.MONTE_CARLO)

    changed = ResamplingPlan.for_accepted(
        accepted,
        scheme=ResamplingScheme.MONTE_CARLO,
        replicate_count=4,
        root_seed=0x1234_5678_90AB_CDEE,
        source_dataset_identity="qualification-dataset",
        output_scope=("A", "B"),
        optimization_projection_identity="direct-trf-projection-v1",
        minimum_successful_count=2,
    )
    assert tuple(request.seed for request in changed.replicates) != tuple(
        request.seed for request in plan.replicates
    )


def test_execution_rejects_a_different_source_dataset_lineage() -> None:
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.MONTE_CARLO)
    mismatched = dataclasses.replace(
        _population(),
        source_dataset_identity="other-dataset",
    )

    with pytest.raises(ResamplingConstructionError, match="source dataset"):
        execute_resampling_evidence(
            accepted,
            plan,
            mismatched,
            _successful_projection(),
        )


def test_seeded_generation_preserves_mc_bs_and_bsn_scientific_schemes() -> None:
    population = _population()

    mc = generate_resampling_draw(
        population,
        ResamplingScheme.MONTE_CARLO,
        seed=101,
    )
    assert mc == generate_resampling_draw(
        population,
        ResamplingScheme.MONTE_CARLO,
        seed=101,
    )
    assert mc.source_indices == tuple(range(5))
    assert mc.standard_errors == population.standard_errors
    assert mc.observations != population.calculated

    bootstrap = generate_resampling_draw(
        population,
        ResamplingScheme.BOOTSTRAP,
        seed=202,
    )
    masked_targets = tuple(
        index for index, value in enumerate(population.mask) if value
    )
    assert tuple(
        population.references[source]
        for source in bootstrap.source_indices[: len(masked_targets)]
    ) == tuple(population.references[target] for target in masked_targets)
    assert bootstrap.source_indices[-1] == 4
    assert bootstrap.observations[-1] == population.observed[-1]

    nucleus = generate_resampling_draw(
        population,
        ResamplingScheme.NUCLEUS_BOOTSTRAP,
        seed=303,
    )
    assert len(nucleus.sampled_nucleus_groups) == 3
    assert set(nucleus.sampled_nucleus_groups) <= {"N1", "N2", "N3"}
    expected_sources = tuple(
        index
        for group in nucleus.sampled_nucleus_groups
        for index, source_group in enumerate(population.nucleus_groups)
        if source_group == group
    )
    assert nucleus.source_indices == expected_sources


def test_resampling_execution_is_evidence_only_scope_complete_and_worker_stable() -> (
    None
):
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.MONTE_CARLO)
    population = _population()

    serial = execute_resampling_evidence(
        accepted,
        plan,
        population,
        _successful_projection(),
        execution=ExecutionSettings(workers=1),
    )
    parallel = execute_resampling_evidence(
        accepted,
        plan,
        population,
        _successful_projection(),
        execution=ExecutionSettings(workers=3),
    )

    assert serial.evidence is not None
    assert parallel.evidence is not None
    assert serial.evidence.identity == parallel.evidence.identity
    assert tuple(outcome.ordinal for outcome in serial.evidence.outcomes) == (
        0,
        1,
        2,
        3,
    )
    assert all(
        outcome.disposition is ReplicateDisposition.SUCCEEDED
        for outcome in serial.evidence.outcomes
    )
    assert all(outcome.success is not None for outcome in serial.evidence.outcomes)
    assert all(
        outcome.success.resolved_items == tuple(outcome.success.resolved_items)
        and tuple(param_id for param_id, _value in outcome.success.resolved_items)
        == plan.output_scope
        for outcome in serial.evidence.outcomes
        if outcome.success is not None
    )
    assert not hasattr(serial.evidence, "accepted_result")
    assert not hasattr(serial.evidence, "commit_authority")


def test_failed_replicates_remain_typed_and_summary_uses_one_common_scope() -> None:
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.BOOTSTRAP, count=3)
    population = _population()

    def project(
        request: ReplicateRequest,
        draw: ResamplingDraw,
    ) -> ProjectedOptimizationResult:
        if request.ordinal == 1:
            return ProjectedOptimizationFailure(
                transformed_data_identity=draw.identity,
                disposition=ReplicateDisposition.FAILED,
                category="solver_unsuccessful",
                message="declared local budget exhausted",
                evaluation_plan_identity="plan-1",
                problem_identity="problem-1",
                execution_identity="execution-1",
            )
        return _successful_projection()(request, draw)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        population,
        project,
        execution=ExecutionSettings(workers=2),
    )
    assert operation.evidence is not None
    evidence = operation.evidence
    assert tuple(outcome.disposition for outcome in evidence.outcomes) == (
        ReplicateDisposition.SUCCEEDED,
        ReplicateDisposition.FAILED,
        ReplicateDisposition.SUCCEEDED,
    )
    assert evidence.completed_count == 3
    assert evidence.successful_count == 2
    assert evidence.failed_count == 1
    assert evidence.claim("MINIMUM_SUCCESSFUL_COVERAGE") is ClaimState.SATISFIED
    assert evidence.claim("COMPLETE_SCOPE_SUCCESS_ROWS") is ClaimState.SATISFIED

    summary_outcome = summarize_resampling_evidence(evidence)
    assert summary_outcome.terminal is SummaryTerminal.COMPLETED
    assert summary_outcome.summary is not None
    summary = summary_outcome.summary
    assert summary.included_ordinals == (0, 2)
    assert summary.claim("COMMON_SCOPE_INCLUSION") is ClaimState.SATISFIED
    assert summary.claim("MINIMUM_SUCCESSFUL_COVERAGE") is ClaimState.SATISFIED
    assert tuple(item.ordinal for item in summary.exclusions) == (1,)
    assert all(row.ordinal in summary.included_ordinals for row in summary.samples)
    assert all(
        tuple(param_id for param_id, _value in row.items) == plan.output_scope
        for row in summary.samples
    )
    assert all(
        entry.value is None or np.isfinite(entry.value)
        for entry in summary.correlations
    )


def test_partial_evidence_and_summary_publish_only_under_declared_contract() -> None:
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.NUCLEUS_BOOTSTRAP, count=4)
    population = _population()

    calls = 0

    def cancel_after_one() -> bool:
        return calls >= 1

    def project(
        request: ReplicateRequest,
        draw: ResamplingDraw,
    ) -> ProjectedOptimizationSuccess:
        nonlocal calls
        calls += 1
        return _successful_projection()(request, draw)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        population,
        project,
        execution=ExecutionSettings(workers=1),
        cancellation_probe=cancel_after_one,
    )
    assert operation.evidence is not None
    assert operation.evidence.completed_count == 1
    assert operation.evidence.successful_count == 1
    assert operation.evidence.lifecycle is ResamplingLifecycle.CANCELLED
    assert (
        operation.evidence.claim("MINIMUM_SUCCESSFUL_COVERAGE") is ClaimState.VIOLATED
    )
    assert operation.unstarted_ordinals == (1, 2, 3)

    summary_outcome = summarize_resampling_evidence(operation.evidence)
    assert summary_outcome.terminal is SummaryTerminal.INSUFFICIENT_COVERAGE
    assert summary_outcome.summary is None
    assert summary_outcome.failure is not None

    initially_cancelled = execute_resampling_evidence(
        accepted,
        plan,
        population,
        project,
        execution=ExecutionSettings(workers=1),
        cancellation_probe=lambda: True,
    )
    assert initially_cancelled.evidence is None
    assert initially_cancelled.unstarted_ordinals == (0, 1, 2, 3)


def test_success_payload_rejects_nan_and_non_atomic_scope() -> None:
    with pytest.raises(ResamplingConstructionError, match="complete output scope"):
        ProjectedOptimizationSuccess(
            transformed_data_identity="draw",
            evaluation_plan_identity="plan",
            problem_identity="problem",
            invocation_identity="invocation",
            execution_identity="execution",
            strategy=OptimizationStrategy.DIRECT_TRF,
            component_outcome_identities=("component",),
            resolved_items=(("A", 1.0),),
            chi_square=1.0,
            expected_output_scope=("A", "B"),
        )
    with pytest.raises(ResamplingConstructionError, match="finite"):
        ProjectedOptimizationSuccess(
            transformed_data_identity="draw",
            evaluation_plan_identity="plan",
            problem_identity="problem",
            invocation_identity="invocation",
            execution_identity="execution",
            strategy=OptimizationStrategy.DIRECT_TRF,
            component_outcome_identities=("component",),
            resolved_items=(("A", np.nan), ("B", 2.0)),
            chi_square=1.0,
            expected_output_scope=("A", "B"),
        )


def test_non_finite_summary_arithmetic_returns_typed_unavailability() -> None:
    accepted = _accepted_anchor()
    plan = _plan(accepted, ResamplingScheme.MONTE_CARLO, count=2)

    def project(
        request: ReplicateRequest,
        draw: ResamplingDraw,
    ) -> ProjectedOptimizationSuccess:
        value = -1.0e308 if request.ordinal == 0 else 1.0e308
        return ProjectedOptimizationSuccess(
            transformed_data_identity=draw.identity,
            evaluation_plan_identity=f"plan-{request.ordinal}",
            problem_identity=f"problem-{request.ordinal}",
            invocation_identity=f"invocation-{request.ordinal}",
            execution_identity=f"execution-{request.ordinal}",
            strategy=OptimizationStrategy.DIRECT_TRF,
            component_outcome_identities=(f"component-{request.ordinal}",),
            resolved_items=(("A", value), ("B", value)),
            chi_square=1.0,
        )

    operation = execute_resampling_evidence(
        accepted,
        plan,
        _population(),
        project,
    )
    assert operation.evidence is not None
    outcome = summarize_resampling_evidence(operation.evidence)
    assert outcome.terminal is SummaryTerminal.SOURCE_INVALID
    assert outcome.summary is None
    assert outcome.failure is not None
    assert outcome.failure.category == "non_finite_summary_arithmetic"
