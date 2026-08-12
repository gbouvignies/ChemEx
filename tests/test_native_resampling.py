"""Behavioral qualification tests for native resampling evidence (#599)."""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
from collections.abc import Callable
from threading import Event, Lock, get_ident
from typing import Any

import numpy as np
import pytest

from chemex.evaluation.native import (
    EvaluationPlan,
    EvaluationResult,
    ProfilePlan,
    ResolvedEvaluationValues,
)
from chemex.optimize.direct_trf import AcceptedFitResult, OptimizationProblem
from chemex.optimize.native_resampling import (
    ClaimState,
    CorrelationAvailability,
    OptimizationStrategy,
    ProjectedOptimizationFailure,
    ProjectedOptimizationResult,
    ProjectedOptimizationSuccess,
    ReplicateDisposition,
    ReplicateExecutionPlan,
    ReplicateRequest,
    ReplicateStage,
    ResamplingConstructionError,
    ResamplingDatasetManifest,
    ResamplingLifecycle,
    ResamplingPlan,
    ResamplingScheme,
    ResamplingSummary,
    ResamplingSummaryOutcome,
    ResamplingSummaryPolicy,
    SummaryTerminal,
    execute_resampling_evidence,
    generate_resampling_draw,
    summarize_resampling_evidence,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ConstraintProgram,
    ParameterRole,
    ScientificFunctionBinder,
)
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.runtime import ExecutionSettings


def _evaluation_identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"schema": 1, "kind": kind, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _profile_plan(
    name: str,
    *,
    profile_ordinal: int,
    offset: int,
    observed: tuple[float, ...],
    errors: tuple[float, ...],
    mask: tuple[bool, ...],
) -> ProfilePlan:
    source_identity = f"source-{name}"
    return ProfilePlan(
        _evaluation_identity(
            "profile-plan",
            (source_identity, 0, profile_ordinal, offset),
        ),
        source_identity,
        0,
        profile_ordinal,
        offset,
        (("a", "A"), ("b", "B")),
        True,
        observed,
        errors,
        mask,
        f"kernel-{name}",
        f"kernel-config-{name}",
        f"spectrometer-{name}",
        f"positions-{name}",
        (len(observed),),
    )


def _native_context() -> tuple[
    AcceptedFitResult,
    ResamplingDatasetManifest,
    OptimizationProblem,
    ActiveParameterization,
]:
    binder = ScientificFunctionBinder("qualification", {})
    program = ConstraintProgram(
        "parameter-model",
        "model",
        "definitions",
        "configuration",
        binder.identity,
        ("A", "B"),
        ("A", "B"),
        (),
        (),
        (),
    )
    parameterization = ActiveParameterization(
        program,
        binder,
        "source-occurrence",
        4,
        (("A", ParameterRole.FIT), ("B", ParameterRole.FIT)),
    )
    evaluation_plan = EvaluationPlan(
        parameterization.evaluator_identity,
        program.fingerprint,
        (
            _profile_plan(
                "P1",
                profile_ordinal=0,
                offset=0,
                observed=(10.0, 20.0, 30.0, 40.0),
                errors=(1.0, 2.0, 1.5, 0.5),
                mask=(True, True, True, True),
            ),
            _profile_plan(
                "P2",
                profile_ordinal=1,
                offset=4,
                observed=(50.0,),
                errors=(2.5,),
                mask=(False,),
            ),
        ),
        ("A", "B"),
    )
    dataset = ResamplingDatasetManifest(
        evaluation_plan,
        (11.0, 19.0, 31.0, 39.0, 49.0),
        (True, True, False, False, False),
        ("N1", "N1", "N2", "N3", "N3"),
        (
            "position=100",
            "position=200",
            "position=300",
            "position=400",
            "position=500",
        ),
    )
    snapshot = AnalysisValuesSnapshot(
        "source-occurrence",
        "model",
        "definitions",
        "configuration",
        4,
        (("A", 1.0), ("B", 2.0)),
    )
    problem = OptimizationProblem(
        evaluation_plan.identity,
        parameterization.identity,
        parameterization.evaluator_identity,
        program.fingerprint,
        "configuration",
        snapshot,
        (("A", 1.0), ("B", 2.0)),
        ("A", "B"),
        (),
        (1.0, 2.0),
        (-100.0, -100.0),
        (100.0, 100.0),
        ("A", "B"),
    )
    result = EvaluationResult(
        evaluation_plan.identity,
        parameterization.evaluator_identity,
        parameterization.evaluator_identity,
        ResolvedEvaluationValues(
            parameterization.evaluator_identity,
            program.fingerprint,
            (("A", 1.0), ("B", 2.0)),
        ),
        np.array(dataset.calculated),
        np.array(dataset.calculated),
        np.array((1.0, -0.5, 2.0 / 3.0, -2.0)),
        (),
    )
    accepted = AcceptedFitResult.for_qualification(
        occurrence_identity="accepted-occurrence",
        problem_identity=problem.identity,
        invocation_identity="accepted-invocation",
        execution_identity="accepted-execution",
        materialization_identity="accepted-materialization",
        parameterization_identity=parameterization.identity,
        evaluator_parameterization_identity=parameterization.evaluator_identity,
        source_occurrence_identity=snapshot.occurrence_identity,
        source_revision=snapshot.revision,
        controlled_ids=problem.controlled_ids,
        vector=problem.start,
        chi_square=0.0,
        evaluation_result=result,
        commit_scope=problem.commit_scope,
        commit_items=problem.independent_items,
        origin_context_identity="accepted-origin",
    )
    return accepted, dataset, problem, parameterization


def _plan(
    scheme: ResamplingScheme,
    *,
    count: int = 4,
    minimum: int = 2,
) -> tuple[AcceptedFitResult, ResamplingPlan]:
    accepted, dataset, problem, parameterization = _native_context()
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=dataset,
        source_problem=problem,
        parameterization=parameterization,
        scheme=scheme,
        replicate_count=count,
        replicate_structural_identities=tuple(
            f"qualification-replicate-{name}"
            for name in ("alpha", "beta", "gamma", "delta")[:count]
        ),
        replicate_component_identities=tuple(
            (f"component-{ordinal}",) for ordinal in range(count)
        ),
        root_seed=0x1234_5678_90AB_CDEF,
        output_scope=("A", "B"),
        output_units=("arbitrary", "arbitrary"),
        minimum_successful_count=minimum,
        strategy=OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "80"),),
    )
    return accepted, plan


def _success(
    prepared: ReplicateExecutionPlan,
    *,
    constant_b: bool = False,
) -> ProjectedOptimizationSuccess:
    ordinal = prepared.request.ordinal
    return ProjectedOptimizationSuccess.for_prepared(
        prepared,
        candidate_vector=(float(ordinal + 1), 5.0 if constant_b else 10.0 - ordinal),
        chi_square=float(ordinal + 1),
    )


def _successful_factory(
    *, constant_b: bool = False
) -> Callable[
    [ReplicateExecutionPlan],
    Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess],
]:
    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        return lambda exact: _success(exact, constant_b=constant_b)

    return factory


def _request(
    dataset: ResamplingDatasetManifest,
    scheme: ResamplingScheme,
    seed: int,
) -> ReplicateRequest:
    return ReplicateRequest(
        "generation-plan",
        dataset.identity,
        scheme,
        "chemex-mc-bs-bsn-v1",
        0,
        "generation-replicate",
        "generation-projection",
        ("generation-component",),
        seed,
    )


def test_replicate_plan_uses_canonical_ordinals_and_identity_derived_seeds() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    assert tuple(request.ordinal for request in plan.replicates) == (0, 1, 2, 3)
    assert len({request.identity for request in plan.replicates}) == 4
    assert len({request.seed for request in plan.replicates}) == 4
    assert plan.source_dataset_identity == plan.dataset.identity

    reordered = ResamplingPlan.for_accepted(
        accepted,
        dataset=plan.dataset,
        source_problem=plan.source_problem,
        parameterization=plan.parameterization,
        scheme=plan.scheme,
        replicate_count=plan.replicate_count,
        replicate_structural_identities=tuple(
            reversed(plan.replicate_structural_identities)
        ),
        replicate_component_identities=tuple(
            reversed(plan.replicate_component_identities)
        ),
        root_seed=plan.root_seed,
        output_scope=plan.output_scope,
        output_units=plan.output_units,
        minimum_successful_count=plan.minimum_successful_count,
        strategy=plan.strategy,
        strategy_settings=plan.strategy_settings,
    )
    assert reordered.identity == plan.identity
    assert tuple(item.seed for item in reordered.replicates) == tuple(
        item.seed for item in plan.replicates
    )


def test_altered_source_calculation_cannot_rebind_the_accepted_native_result() -> None:
    accepted, dataset, problem, parameterization = _native_context()
    altered = dataclasses.replace(
        dataset,
        calculated=(12.0, *dataset.calculated[1:]),
    )

    with pytest.raises(ResamplingConstructionError, match="accepted native lineage"):
        ResamplingPlan.for_accepted(
            accepted,
            dataset=altered,
            source_problem=problem,
            parameterization=parameterization,
            scheme=ResamplingScheme.MONTE_CARLO,
            replicate_count=1,
            replicate_structural_identities=("altered-source",),
            replicate_component_identities=(("component",),),
            root_seed=7,
            output_scope=("A", "B"),
            output_units=("arbitrary", "arbitrary"),
            minimum_successful_count=1,
        )


def test_source_dataset_record_rejects_an_altered_observation_with_stale_identity() -> (
    None
):
    _accepted, dataset, _problem, _parameterization = _native_context()
    first_profile = dataset.evaluation_plan.profiles[0]
    changed_plan = dataclasses.replace(
        dataset.evaluation_plan,
        profiles=(
            dataclasses.replace(
                first_profile,
                experimental_values=(
                    first_profile.experimental_values[0] + 0.5,
                    *first_profile.experimental_values[1:],
                ),
            ),
            *dataset.evaluation_plan.profiles[1:],
        ),
    )
    altered = dataclasses.replace(dataset, evaluation_plan=changed_plan)

    assert altered.identity != dataset.identity
    record = altered.to_record()
    record["identity"] = dataset.identity
    with pytest.raises(ValueError, match="identity"):
        ResamplingDatasetManifest.from_record(record)


def test_seeded_generation_preserves_mc_bs_and_bsn_scientific_schemes() -> None:
    _accepted, dataset, _problem, _parameterization = _native_context()

    mc = generate_resampling_draw(
        dataset,
        _request(dataset, ResamplingScheme.MONTE_CARLO, 101),
    )
    np.testing.assert_allclose(
        mc.observations,
        (
            10.209847500036986,
            14.930749036336255,
            31.904952620387146,
            39.372147264939954,
            48.225783000334324,
        ),
        rtol=0.0,
        atol=1.0e-14,
    )
    assert mc.source_indices == tuple(range(5))

    bootstrap = generate_resampling_draw(
        dataset,
        _request(dataset, ResamplingScheme.BOOTSTRAP, 202),
    )
    assert bootstrap.source_indices == (0, 0, 2, 2, 4)
    assert all(
        dataset.profile_blocks[target] == dataset.profile_blocks[source]
        for target, source in enumerate(bootstrap.source_indices)
    )
    assert all(
        dataset.references[target] == dataset.references[source]
        for target, source in enumerate(bootstrap.source_indices)
    )

    nucleus = generate_resampling_draw(
        dataset,
        _request(dataset, ResamplingScheme.NUCLEUS_BOOTSTRAP, 303),
    )
    assert nucleus.sampled_nucleus_groups == ("N2", "N1", "N3")
    assert nucleus.source_indices == (2, 0, 1, 3, 4)


def test_altered_generated_observation_cannot_retain_or_rebind_draw_identity() -> None:
    _accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    request = plan.replicates[0]
    draw = generate_resampling_draw(plan.dataset, request)
    altered = dataclasses.replace(
        draw,
        observations=(draw.observations[0] + 1.0, *draw.observations[1:]),
    )

    assert altered.identity != draw.identity
    with pytest.raises(ResamplingConstructionError, match="canonical seeded"):
        ReplicateExecutionPlan.prepare(plan, request, altered)


@pytest.mark.parametrize(
    ("field_name", "foreign_value"),
    (
        ("accepted_result_identity", "foreign-accepted-result"),
        ("accepted_occurrence_identity", "foreign-accepted-occurrence"),
        ("source_snapshot_occurrence_identity", "foreign-source-snapshot"),
        ("source_snapshot_revision", 99),
        ("evaluation_plan_identity", "foreign-evaluation-plan"),
        ("problem_identity", "foreign-problem"),
        ("parameterization_identity", "foreign-parameterization"),
        ("invocation_identity", "foreign-invocation"),
        ("execution_identity", "foreign-execution"),
        ("workflow_identity", "foreign-workflow"),
        ("strategy", OptimizationStrategy.DE_DIRECT_TRF),
        ("component_identities", ("foreign-component",)),
        ("component_outcome_identities", ("foreign-component-outcome",)),
        ("candidate_identity", "foreign-candidate"),
        ("materialization_identity", "foreign-materialization"),
        ("evaluation_identity", "foreign-evaluation"),
        ("prepared_replicate_identity", "foreign-prepared-replicate"),
        ("requested_output_units", ("foreign-unit", "foreign-unit")),
    ),
)
def test_foreign_native_lineage_never_enters_successful_evidence(
    field_name: str,
    foreign_value: object,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        def execute(exact: ReplicateExecutionPlan) -> ProjectedOptimizationSuccess:
            return dataclasses.replace(
                _success(exact),
                **{field_name: foreign_value},
            )

        return execute

    operation = execute_resampling_evidence(accepted, plan, factory)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(outcome.success is None for outcome in operation.evidence.outcomes)
    assert all(
        outcome.failure is not None
        and outcome.failure.stage is ReplicateStage.VALIDATING
        for outcome in operation.evidence.outcomes
    )


def test_incomplete_requested_scope_never_enters_successful_evidence() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        return lambda exact: dataclasses.replace(
            _success(exact),
            requested_output_scope=("A",),
            requested_output_units=("arbitrary",),
        )

    operation = execute_resampling_evidence(accepted, plan, factory)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "requested_scope_lineage_mismatch"
        for outcome in operation.evidence.outcomes
    )


def test_foreign_transformed_data_with_correct_values_and_dimensions_fails_closed() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.BOOTSTRAP, count=2)

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        return lambda exact: dataclasses.replace(
            _success(exact),
            transformed_data_identity="foreign-draw-with-identical-values",
        )

    operation = execute_resampling_evidence(accepted, plan, factory)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "transformed_data_lineage_mismatch"
        for outcome in operation.evidence.outcomes
    )


def test_candidate_is_freshly_resolved_and_no_fit_or_commit_authority_is_exposed() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    before = accepted
    prepared_problems: list[OptimizationProblem] = []

    def factory(
        prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        prepared_problems.append(prepared.problem)
        return lambda exact: _success(exact)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        factory,
    )

    assert accepted is before
    assert operation.evidence is not None
    assert tuple(
        outcome.success.resolved_items
        for outcome in operation.evidence.outcomes
        if outcome.success is not None
    ) == ((("A", 1.0), ("B", 10.0)), (("A", 2.0), ("B", 9.0)))
    assert all(
        not outcome.success or not hasattr(outcome.success, "accepted_result")
        for outcome in operation.evidence.outcomes
    )
    assert all(
        not outcome.success or not hasattr(outcome.success, "commit_authority")
        for outcome in operation.evidence.outcomes
    )
    assert all(not problem.acceptance_authority for problem in prepared_problems)


def test_out_of_bounds_candidate_fails_fresh_problem_resolution() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        return lambda exact: dataclasses.replace(
            _success(exact),
            candidate_vector=(101.0, 5.0),
        )

    operation = execute_resampling_evidence(accepted, plan, factory)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "candidate_resolution_failure"
        and outcome.stage is ReplicateStage.VALIDATING
        for outcome in operation.evidence.outcomes
    )


class _SingleOwnerExecutor:
    def __init__(
        self,
        prepared: ReplicateExecutionPlan,
        owners: list[tuple[object, int, int]],
        owner_lock: Lock,
        releases: tuple[Event, ...] | None,
    ) -> None:
        self.prepared = prepared
        self.owners = owners
        self.owner_lock = owner_lock
        self.releases = releases
        self.used = False

    def __call__(self, exact: ReplicateExecutionPlan) -> ProjectedOptimizationSuccess:
        if self.used or exact is not self.prepared:
            raise RuntimeError("executor/workspace reused across replicates")
        self.used = True
        ordinal = exact.request.ordinal
        if self.releases is not None:
            if ordinal + 1 < len(self.releases):
                assert self.releases[ordinal + 1].wait(timeout=1.0)
            self.releases[ordinal].set()
        with self.owner_lock:
            self.owners.append((self, ordinal, get_ident()))
        return _success(exact)


def _isolated_factory(
    owners: list[tuple[object, int, int]],
    owner_lock: Lock,
    releases: tuple[Event, ...] | None = None,
) -> Callable[[ReplicateExecutionPlan], _SingleOwnerExecutor]:
    return lambda prepared: _SingleOwnerExecutor(
        prepared,
        owners,
        owner_lock,
        releases,
    )


def test_serial_and_reordered_multi_worker_execution_use_fresh_single_owners() -> None:
    accepted, plan = _plan(ResamplingScheme.NUCLEUS_BOOTSTRAP)
    serial_owners: list[tuple[object, int, int]] = []
    parallel_owners: list[tuple[object, int, int]] = []
    serial = execute_resampling_evidence(
        accepted,
        plan,
        _isolated_factory(serial_owners, Lock()),
        execution=ExecutionSettings(workers=1),
    )
    releases = tuple(Event() for _ in plan.replicates)
    parallel = execute_resampling_evidence(
        accepted,
        plan,
        _isolated_factory(parallel_owners, Lock(), releases),
        execution=ExecutionSettings(workers=4),
    )

    assert serial.evidence is not None
    assert parallel.evidence is not None
    assert serial.evidence.identity == parallel.evidence.identity
    assert len({item[0] for item in serial_owners}) == plan.replicate_count
    assert len({item[0] for item in parallel_owners}) == plan.replicate_count
    assert tuple(item[1] for item in parallel_owners) == (3, 2, 1, 0)
    assert tuple(outcome.ordinal for outcome in parallel.evidence.outcomes) == (
        0,
        1,
        2,
        3,
    )
    assert tuple(item.seed for item in serial.evidence.plan.replicates) == tuple(
        item.seed for item in parallel.evidence.plan.replicates
    )


def test_one_executor_failure_does_not_corrupt_other_replicate_evidence() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    def factory(
        prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        if prepared.request.ordinal == 1:
            return lambda _exact: (_ for _ in ()).throw(RuntimeError("owner failed"))
        return lambda exact: _success(exact)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        factory,
        execution=ExecutionSettings(workers=3),
    )

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 3
    assert operation.evidence.outcomes[1].failure is not None
    assert (
        operation.evidence.outcomes[1].failure.category == "projected_execution_failure"
    )


@pytest.mark.parametrize("workers", (1, 3))
def test_shared_executor_instance_is_rejected_by_single_owner_contract(
    workers: int,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    def shared(exact: ReplicateExecutionPlan) -> ProjectedOptimizationSuccess:
        return _success(exact)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        lambda _prepared: shared,
        execution=ExecutionSettings(workers=workers),
    )

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 1
    assert (
        sum(
            outcome.failure is not None
            and outcome.failure.category == "executor_creation_failure"
            for outcome in operation.evidence.outcomes
        )
        == plan.replicate_count - 1
    )


def test_partial_evidence_preserves_unstarted_ordinals_and_valid_summary() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)
    calls = 0

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        def execute(exact: ReplicateExecutionPlan) -> ProjectedOptimizationSuccess:
            nonlocal calls
            calls += 1
            return _success(exact)

        return execute

    operation = execute_resampling_evidence(
        accepted,
        plan,
        factory,
        cancellation_probe=lambda: calls >= 2,
    )

    assert operation.evidence is not None
    assert operation.evidence.completed_count == 2
    assert operation.evidence.lifecycle is ResamplingLifecycle.CANCELLED
    assert operation.unstarted_ordinals == (2, 3)
    assert tuple(outcome.ordinal for outcome in operation.evidence.outcomes) == (
        0,
        1,
        2,
        3,
    )
    assert all(
        operation.evidence.outcomes[index].disposition
        is ReplicateDisposition.NOT_STARTED
        for index in operation.unstarted_ordinals
    )
    summary = summarize_resampling_evidence(operation.evidence)
    assert summary.terminal is SummaryTerminal.COMPLETED
    assert summary.summary is not None
    assert summary.summary.unstarted_ordinals == (2, 3)


def test_summary_numerical_references_and_zero_variance_are_canonical() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=3)
    operation = execute_resampling_evidence(
        accepted,
        plan,
        _successful_factory(constant_b=True),
    )
    assert operation.evidence is not None
    outcome = summarize_resampling_evidence(operation.evidence)
    assert outcome.summary is not None
    summary = outcome.summary
    distributions = {item.parameter_id: item for item in summary.distributions}
    # Small exact values use stable binary64 reductions; this absolute tolerance is
    # roughly one ulp here and detects any change of estimator or percentile policy.
    assert distributions["A"].mean == pytest.approx(2.0, rel=0.0, abs=1.0e-15)
    assert distributions["A"].standard_deviation == pytest.approx(
        1.0, rel=0.0, abs=1.0e-15
    )
    assert distributions["A"].median == pytest.approx(2.0, rel=0.0, abs=1.0e-15)
    assert distributions["A"].percentile_95_lower == pytest.approx(
        1.05, rel=0.0, abs=1.0e-15
    )
    assert distributions["A"].percentile_95_upper == pytest.approx(
        2.95, rel=0.0, abs=1.0e-15
    )
    covariance = {
        (item.parameter_a, item.parameter_b): item.value for item in summary.covariance
    }
    assert covariance[("A", "A")] == pytest.approx(1.0, rel=0.0, abs=1.0e-15)
    assert covariance[("A", "B")] == pytest.approx(0.0, rel=0.0, abs=1.0e-15)
    b_entries = tuple(
        entry
        for entry in summary.correlations
        if entry.parameter_a == "B" or entry.parameter_b == "B"
    )
    assert all(
        entry.availability is CorrelationAvailability.ZERO_VARIANCE
        and entry.value is None
        for entry in b_entries
    )


def test_summary_record_rejects_every_independent_tampering_vector() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=3)
    operation = execute_resampling_evidence(accepted, plan, _successful_factory())
    assert operation.evidence is not None
    policy = ResamplingSummaryPolicy()
    outcome = summarize_resampling_evidence(operation.evidence, policy)
    assert outcome.summary is not None
    summary = outcome.summary
    evidence = operation.evidence

    def mutate(path: tuple[object, ...], value: object) -> None:
        record = copy.deepcopy(summary.to_record())
        target: Any = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ResamplingConstructionError, match="canonical"):
            ResamplingSummary.from_record(record, evidence, policy)

    mutate(("included_ordinals", 0), 2)
    mutate(("samples", 0, "items", 0, 1), (9.0).hex())
    mutate(("output_scope", 0), "B")
    mutate(("output_units", 0), "foreign-unit")
    mutate(("distributions", 0, "percentile_95_lower"), (0.0).hex())
    mutate(("covariance", 0, "value"), (99.0).hex())
    mutate(("correlations", 0, "value"), (0.5).hex())
    mutate(("correlations", 0, "availability"), "zero_variance")
    mutate(("evidence_identity",), "foreign-evidence")
    mutate(("policy_identity",), "foreign-policy")
    mutate(("claims", 0, "state"), "violated")

    changed_policy = ResamplingSummaryPolicy(percentile_lower=5.0)
    with pytest.raises(ResamplingConstructionError, match="canonical"):
        ResamplingSummary.from_record(
            summary.to_record(),
            evidence,
            changed_policy,
        )
    with pytest.raises(ResamplingConstructionError, match="evidence identity"):
        ResamplingSummaryOutcome(
            SummaryTerminal.COMPLETED,
            summary=summary,
            claimed_evidence_identity="foreign-evidence",
        )


def test_summary_constructor_cannot_accept_caller_supplied_arrays_or_claims() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan, _successful_factory())
    assert operation.evidence is not None
    summary = summarize_resampling_evidence(operation.evidence).summary
    assert summary is not None

    with pytest.raises(ResamplingConstructionError, match="source evidence"):
        dataclasses.replace(summary, _construction_key=object())


@pytest.mark.parametrize(
    ("failure", "expected_stage", "draw_present"),
    (
        (KeyboardInterrupt(), ReplicateStage.GENERATING, False),
        (RuntimeError("generation failed"), ReplicateStage.GENERATING, False),
    ),
)
def test_pre_draw_generation_terminal_is_typed_and_factory_is_not_created(
    monkeypatch: pytest.MonkeyPatch,
    failure: BaseException,
    expected_stage: ReplicateStage,
    draw_present: bool,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)
    factory_calls = 0

    def fail_generation(
        _dataset: ResamplingDatasetManifest,
        _request: ReplicateRequest,
    ) -> None:
        raise failure

    def factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationSuccess]:
        nonlocal factory_calls
        factory_calls += 1
        return _success

    monkeypatch.setattr(
        "chemex.optimize.native_resampling.generate_resampling_draw",
        fail_generation,
    )
    operation = execute_resampling_evidence(accepted, plan, factory)

    assert factory_calls == 0
    assert operation.evidence is not None
    first = operation.evidence.outcomes[0]
    assert first.stage is expected_stage
    assert (first.draw_identity is not None) is draw_present
    assert first.disposition is (
        ReplicateDisposition.INTERRUPTED
        if isinstance(failure, KeyboardInterrupt)
        else ReplicateDisposition.FAILED
    )
    assert tuple(outcome.ordinal for outcome in operation.evidence.outcomes) == (
        0,
        1,
        2,
        3,
    )
    expected_remaining = (
        ReplicateDisposition.NOT_STARTED
        if isinstance(failure, KeyboardInterrupt)
        else ReplicateDisposition.FAILED
    )
    assert all(
        outcome.disposition is expected_remaining
        for outcome in operation.evidence.outcomes[1:]
    )


def test_interruption_after_draw_before_executor_creation_is_typed() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    def interrupt_factory(
        _prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationResult]:
        raise KeyboardInterrupt

    operation = execute_resampling_evidence(accepted, plan, interrupt_factory)

    assert operation.evidence is not None
    first = operation.evidence.outcomes[0]
    assert first.disposition is ReplicateDisposition.INTERRUPTED
    assert first.stage is ReplicateStage.EXECUTOR_CREATING
    assert first.draw_identity is not None
    assert all(
        outcome.disposition is ReplicateDisposition.NOT_STARTED
        for outcome in operation.evidence.outcomes[1:]
    )


@pytest.mark.parametrize("workers", (1, 3))
def test_optimization_interruption_preserves_complete_ordinal_accounting(
    workers: int,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    def factory(
        prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationResult]:
        if prepared.request.ordinal == 0:
            return lambda _exact: (_ for _ in ()).throw(KeyboardInterrupt())
        return lambda exact: _success(exact)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        factory,
        execution=ExecutionSettings(workers=workers),
    )

    assert operation.evidence is not None
    assert tuple(outcome.ordinal for outcome in operation.evidence.outcomes) == (
        0,
        1,
        2,
        3,
    )
    interrupted = operation.evidence.outcomes[0]
    assert interrupted.disposition is ReplicateDisposition.INTERRUPTED
    assert interrupted.stage is ReplicateStage.EXECUTING
    assert interrupted.draw_identity is not None
    assert len(operation.evidence.outcomes) == plan.replicate_count


def test_projected_failure_remains_typed_and_summary_uses_common_scope() -> None:
    accepted, plan = _plan(ResamplingScheme.BOOTSTRAP, count=3)

    def factory(
        prepared: ReplicateExecutionPlan,
    ) -> Callable[[ReplicateExecutionPlan], ProjectedOptimizationResult]:
        if prepared.request.ordinal != 1:
            return lambda exact: _success(exact)

        def fail(exact: ReplicateExecutionPlan) -> ProjectedOptimizationFailure:
            return ProjectedOptimizationFailure(
                exact.draw.identity,
                ReplicateDisposition.FAILED,
                "solver_unsuccessful",
                "declared local budget exhausted",
                optimization_projection_identity=(
                    exact.request.optimization_projection_identity
                ),
                stage=ReplicateStage.EXECUTING,
                prepared_replicate_identity=exact.identity,
                evaluation_plan_identity=exact.evaluation_plan.identity,
                problem_identity=exact.problem.identity,
                invocation_identity=exact.invocation_identity,
            )

        return fail

    operation = execute_resampling_evidence(
        accepted,
        plan,
        factory,
        execution=ExecutionSettings(workers=2),
    )
    assert operation.evidence is not None
    evidence = operation.evidence
    assert evidence.successful_count == 2
    assert evidence.failed_count == 1
    assert evidence.claim("MINIMUM_SUCCESSFUL_COVERAGE") is ClaimState.SATISFIED

    outcome = summarize_resampling_evidence(evidence)
    assert outcome.summary is not None
    assert outcome.summary.included_ordinals == (0, 2)
    assert tuple(item.ordinal for item in outcome.summary.exclusions) == (1,)
