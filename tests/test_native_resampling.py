"""Behavioral qualification tests for native resampling evidence (#599)."""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
from collections.abc import Callable
from copy import deepcopy
from threading import Event, Lock
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pydantic import BaseModel

import chemex.optimize.native_resampling as resampling_module
from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
    ResolvedEvaluationValues,
)
from chemex.nmr.spectrometer import Spectrometer
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
    ValidatedReplicateSuccess,
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
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.printers.data import Printer
from chemex.runtime import ExecutionSettings


class _KernelSettings(BaseModel):
    kind: str = "linear-test-kernel"


class _LinearSpectrometer:
    def __init__(self, name: str) -> None:
        self.spin_system = SpinSystem.from_name(name)
        self.values = {"a": 0.0, "b": 0.0}

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _LinearSpectrometer:
        return deepcopy(self)

    def native_kernel_descriptor(self) -> dict[str, object]:
        return {"kind": "linear-test-spectrometer"}


class _LinearPulseSequence:
    settings = _KernelSettings()

    def calculate(self, spectrometer: _LinearSpectrometer, data: Data) -> np.ndarray:
        return spectrometer.values["a"] + spectrometer.values["b"] * np.asarray(
            data.metadata, dtype=np.float64
        )

    def is_reference(self, metadata: np.ndarray) -> np.ndarray:
        return metadata < 10.0


def _profile(
    name: str,
    observed: tuple[float, ...],
    errors: tuple[float, ...],
    metadata: tuple[float, ...],
    mask: tuple[bool, ...],
) -> Profile:
    data = Data(
        exp=np.asarray(observed, dtype=np.float64),
        err=np.asarray(errors, dtype=np.float64),
        metadata=np.asarray(metadata, dtype=np.float64),
    )
    data.mask = np.asarray(mask, dtype=np.bool_)
    return Profile(
        data,
        cast("Spectrometer", _LinearSpectrometer(name)),
        cast("PulseSequence", _LinearPulseSequence()),
        {"a": "A", "b": "B"},
        cast("Printer", None),
        is_scaled=False,
    )


def _native_context() -> tuple[
    AcceptedFitResult,
    ResamplingDatasetManifest,
    OptimizationProblem,
    ActiveParameterization,
    EvaluationEngine,
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
    snapshot = AnalysisValuesSnapshot(
        "source-occurrence",
        "model",
        "definitions",
        "configuration",
        4,
        (("A", 1.0), ("B", 2.0)),
    )
    profiles = (
        _profile(
            "1N",
            (10.0, 20.0, 30.0, 40.0),
            (1.0, 2.0, 1.5, 0.5),
            (5.0, 9.0, 15.0, 19.0),
            (True, True, True, True),
        ),
        _profile("2N", (50.0,), (2.5,), (24.0,), (False,)),
    )
    experiments = cast("Any", (SimpleNamespace(profiles=profiles),))
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    evaluation_plan = engine.plan
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(snapshot),
    )
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
    dataset = ResamplingDatasetManifest(
        evaluation_plan,
        tuple(float(value) for value in result.normalized_calculations),
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
    return accepted, dataset, problem, parameterization, engine


def _plan(
    scheme: ResamplingScheme,
    *,
    count: int = 4,
    minimum: int = 2,
) -> tuple[AcceptedFitResult, ResamplingPlan]:
    accepted, dataset, problem, parameterization, engine = _native_context()
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=dataset,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
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
    projected: ProjectedOptimizationResult,
) -> ProjectedOptimizationSuccess:
    assert isinstance(projected, ProjectedOptimizationSuccess)
    assert projected.prepared_replicate_identity == prepared.identity
    return projected


def _successful_factory() -> Callable[
    [ReplicateExecutionPlan, ProjectedOptimizationResult],
    ProjectedOptimizationSuccess,
]:
    return _success


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
        source_engine=plan.source_engine,
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
    accepted, dataset, problem, parameterization, engine = _native_context()
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
            source_engine=engine,
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
    _accepted, dataset, _problem, _parameterization, _engine = _native_context()
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
    _accepted, dataset, _problem, _parameterization, _engine = _native_context()

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
        ("prepared_replicate_identity", "foreign-prepared-replicate"),
        ("controlled_ids", ("B", "A")),
    ),
)
def test_foreign_native_lineage_never_enters_successful_evidence(
    field_name: str,
    foreign_value: object,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def hook(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(projected, **{field_name: foreign_value})

    operation = execute_resampling_evidence(accepted, plan, candidate_test_hook=hook)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(outcome.success is None for outcome in operation.evidence.outcomes)
    assert all(
        outcome.failure is not None
        and outcome.failure.stage is ReplicateStage.VALIDATING
        for outcome in operation.evidence.outcomes
    )


def test_incomplete_or_reordered_controlled_scope_never_enters_successful_evidence() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def hook(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(projected, controlled_ids=("B", "A"))

    operation = execute_resampling_evidence(accepted, plan, candidate_test_hook=hook)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "controlled_coordinates_lineage_mismatch"
        for outcome in operation.evidence.outcomes
    )


def test_foreign_transformed_data_with_correct_values_and_dimensions_fails_closed() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.BOOTSTRAP, count=2)

    def hook(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(
            projected, transformed_data_identity="foreign-draw-with-identical-values"
        )

    operation = execute_resampling_evidence(accepted, plan, candidate_test_hook=hook)

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

    def hook(
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        prepared_problems.append(prepared.problem)
        return projected

    operation = execute_resampling_evidence(
        accepted,
        plan,
        candidate_test_hook=hook,
    )

    assert accepted is before
    assert operation.evidence is not None
    assert operation.evidence.successful_count == 2
    assert all(
        outcome.success is not None and outcome.success.fresh_evaluation_count == 1
        for outcome in operation.evidence.outcomes
    )
    assert all(
        not outcome.success or not hasattr(outcome.success, "accepted_result")
        for outcome in operation.evidence.outcomes
    )
    assert all(
        not outcome.success or not hasattr(outcome.success, "commit_authority")
        for outcome in operation.evidence.outcomes
    )
    assert all(not problem.acceptance_authority for problem in prepared_problems)


def test_backend_chi_square_is_diagnostic_and_fresh_chi_square_is_authoritative() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=1, minimum=1)
    backend_value = 9.87654321e123

    def corrupt_backend_diagnostic(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(projected, backend_chi_square=backend_value)

    operation = execute_resampling_evidence(
        accepted,
        plan,
        candidate_test_hook=corrupt_backend_diagnostic,
    )

    assert operation.evidence is not None
    success = operation.evidence.outcomes[0].success
    assert success is not None
    assert success.backend_chi_square == backend_value
    assert success.backend_chi_square_agrees is False
    assert success.chi_square != backend_value
    assert np.isfinite(success.chi_square)
    assert success.fresh_evaluation_count == 1


def test_foreign_candidate_lineage_fails_before_fresh_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    original = EvaluationEngine.resampled
    bound_engines: list[EvaluationEngine] = []

    def record_binding(
        engine: EvaluationEngine,
        evaluation_plan: EvaluationPlan,
        bindings: object,
    ) -> EvaluationEngine:
        result = original(engine, evaluation_plan, cast("Any", bindings))
        bound_engines.append(result)
        return result

    def foreign_problem(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(projected, problem_identity="foreign-problem")

    monkeypatch.setattr(EvaluationEngine, "resampled", record_binding)
    operation = execute_resampling_evidence(
        accepted,
        plan,
        candidate_test_hook=foreign_problem,
    )

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert len(bound_engines) == plan.replicate_count


def test_fresh_validation_evaluator_failure_is_typed_and_excluded(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    original = EvaluationEngine.resampled
    binding_count = 0

    def fail_validation_binding(
        engine: EvaluationEngine,
        evaluation_plan: EvaluationPlan,
        bindings: object,
    ) -> EvaluationEngine:
        nonlocal binding_count
        binding_count += 1
        result = original(engine, evaluation_plan, cast("Any", bindings))
        if binding_count % 2 == 0:

            def fail_new_evaluator() -> object:
                raise RuntimeError("fresh validation workspace failed")

            result.new_evaluator = fail_new_evaluator  # type: ignore[method-assign]
        return result

    monkeypatch.setattr(EvaluationEngine, "resampled", fail_validation_binding)
    operation = execute_resampling_evidence(accepted, plan)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "fresh_evaluation_failure"
        for outcome in operation.evidence.outcomes
    )
    assert summarize_resampling_evidence(operation.evidence).summary is None


def test_each_eligible_candidate_receives_exactly_one_fresh_validation_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=3)
    original_resampled = EvaluationEngine.resampled
    original_new_evaluator = EvaluationEngine.new_evaluator
    binding_count = 0
    validation_engines: list[EvaluationEngine] = []
    validation_engine_ids: set[int] = set()
    evaluation_counts: dict[int, int] = {}

    def identify_validation_engine(
        engine: EvaluationEngine,
        evaluation_plan: EvaluationPlan,
        bindings: object,
    ) -> EvaluationEngine:
        nonlocal binding_count
        binding_count += 1
        result = original_resampled(engine, evaluation_plan, cast("Any", bindings))
        if binding_count % 2 == 0:
            # Keep engines alive so CPython cannot reuse their IDs during the test.
            validation_engines.append(result)
            validation_engine_ids.add(id(result))
        return result

    def count_validation_evaluation(engine: EvaluationEngine) -> object:
        evaluator = original_new_evaluator(engine)
        if id(engine) in validation_engine_ids:
            original_evaluate = evaluator.evaluate

            def counted(frame: EvaluationFrame) -> object:
                evaluation_counts[id(engine)] = evaluation_counts.get(id(engine), 0) + 1
                return original_evaluate(frame)

            evaluator.evaluate = counted  # type: ignore[method-assign]
        return evaluator

    monkeypatch.setattr(EvaluationEngine, "resampled", identify_validation_engine)
    monkeypatch.setattr(EvaluationEngine, "new_evaluator", count_validation_evaluation)
    operation = execute_resampling_evidence(accepted, plan)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == plan.replicate_count
    assert tuple(evaluation_counts.values()) == (1, 1, 1)


def test_out_of_bounds_candidate_fails_fresh_problem_resolution() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)

    def hook(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        assert isinstance(projected, ProjectedOptimizationSuccess)
        return dataclasses.replace(projected, candidate_vector=(101.0, 5.0))

    operation = execute_resampling_evidence(accepted, plan, candidate_test_hook=hook)

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 0
    assert all(
        outcome.failure is not None
        and outcome.failure.category == "candidate_vector_lineage_mismatch"
        and outcome.stage is ReplicateStage.VALIDATING
        for outcome in operation.evidence.outcomes
    )


def test_serial_and_reordered_multi_worker_execution_use_fresh_single_owners(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.NUCLEUS_BOOTSTRAP)
    engines: list[EvaluationEngine] = []
    evaluators: list[object] = []
    engine_lock = Lock()
    original_resampled = EvaluationEngine.resampled
    original_new_evaluator = EvaluationEngine.new_evaluator

    def record_resampled(
        engine: EvaluationEngine,
        evaluation_plan: EvaluationPlan,
        bindings: object,
    ) -> EvaluationEngine:
        result = original_resampled(engine, evaluation_plan, cast("Any", bindings))
        with engine_lock:
            engines.append(result)
        return result

    def record_evaluator(engine: EvaluationEngine) -> object:
        evaluator = original_new_evaluator(engine)
        with engine_lock:
            evaluators.append(evaluator)
        return evaluator

    monkeypatch.setattr(EvaluationEngine, "resampled", record_resampled)
    monkeypatch.setattr(EvaluationEngine, "new_evaluator", record_evaluator)
    serial = execute_resampling_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=1),
    )
    serial_engine_count = len(engines)
    serial_evaluator_count = len(evaluators)
    releases = tuple(Event() for _ in plan.replicates)
    completion_order: list[int] = []

    def reverse_completion(
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        ordinal = prepared.request.ordinal
        if ordinal + 1 < len(releases):
            assert releases[ordinal + 1].wait(timeout=5.0)
        releases[ordinal].set()
        completion_order.append(ordinal)
        return projected

    parallel = execute_resampling_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=4),
        candidate_test_hook=reverse_completion,
    )

    assert serial.evidence is not None
    assert parallel.evidence is not None
    assert serial.evidence.identity == parallel.evidence.identity
    # Each replicate owns optimizer and validation engines; recursive evidence
    # qualification independently rebinds the immutable compatibility descriptor.
    assert serial_engine_count == 3 * plan.replicate_count
    assert len(engines) == 6 * plan.replicate_count
    assert len({id(item) for item in engines}) == len(engines)
    assert serial_evaluator_count >= 2 * plan.replicate_count
    assert len({id(item) for item in evaluators}) == len(evaluators)
    assert tuple(completion_order) == (3, 2, 1, 0)
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

    def hook(
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        if prepared.request.ordinal == 1:
            raise RuntimeError("owner failed")
        return projected

    operation = execute_resampling_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=3),
        candidate_test_hook=hook,
    )

    assert operation.evidence is not None
    assert operation.evidence.successful_count == 3
    assert operation.evidence.outcomes[1].failure is not None
    assert (
        operation.evidence.outcomes[1].failure.category == "projected_execution_failure"
    )


def test_fresh_lambdas_wrapping_shared_backend_cannot_enter_authoritative_path() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)
    shared_backend = object()
    wrappers = tuple(lambda: shared_backend for _ in plan.replicates)

    with pytest.raises(TypeError, match="executor_factory"):
        execute_resampling_evidence(
            accepted,
            plan,
            executor_factory=lambda prepared: wrappers[prepared.request.ordinal],
        )


def test_partial_evidence_preserves_unstarted_ordinals_and_valid_summary() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)
    calls = 0

    def hook(
        _prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        nonlocal calls
        calls += 1
        return projected

    operation = execute_resampling_evidence(
        accepted,
        plan,
        cancellation_probe=lambda: calls >= 2,
        candidate_test_hook=hook,
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


def test_summary_numerical_references_and_zero_variance_are_canonical(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=3)

    def fixed_candidate(
        prepared: ReplicateExecutionPlan,
        _engine: EvaluationEngine,
    ) -> ProjectedOptimizationResult:
        ordinal = prepared.request.ordinal
        return ProjectedOptimizationSuccess.for_prepared(
            prepared,
            candidate_vector=(float(ordinal + 1), 5.0),
            backend_chi_square=float(ordinal + 1),
            execution_identity=f"test-execution-{ordinal}",
        )

    monkeypatch.setattr(resampling_module, "_execute_native_candidate", fixed_candidate)
    operation = execute_resampling_evidence(accepted, plan)
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
        1.0, rel=0.0, abs=1.0e-15
    )
    assert distributions["A"].percentile_95_upper == pytest.approx(
        3.0, rel=0.0, abs=1.0e-15
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
    operation = execute_resampling_evidence(accepted, plan)
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
    with pytest.raises(TypeError, match="evidence_identity"):
        ResamplingSummaryOutcome(
            SummaryTerminal.COMPLETED,
            summary=summary,
            evidence_identity="foreign-evidence",
        )


def test_summary_constructor_cannot_accept_caller_supplied_arrays_or_claims() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    summary = summarize_resampling_evidence(operation.evidence).summary
    assert summary is not None

    with pytest.raises(TypeError, match="init=False"):
        dataclasses.replace(
            summary,
            covariance=(dataclasses.replace(summary.covariance[0], value=999.0),),
        )
    for field_name in (
        "included_ordinals",
        "samples",
        "distributions",
        "correlations",
        "output_scope",
        "output_units",
        "claims",
    ):
        with pytest.raises(TypeError, match="init=False"):
            dataclasses.replace(
                summary,
                **{field_name: getattr(summary, field_name)},
            )

    changed_policy = dataclasses.replace(summary.policy, percentile_lower=5.0)
    recomputed = dataclasses.replace(summary, policy=changed_policy)
    assert recomputed.policy_identity == changed_policy.identity
    assert recomputed.identity != summary.identity


def test_validated_success_cannot_be_retargeted_and_rehashed_into_summary() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    source_values = accepted.evaluation_result.resolved_values.ordered_items()
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    outcome = operation.evidence.outcomes[0]
    assert outcome.success is not None

    for field_name, value in (
        ("resolved_items", (("A", 999.0), ("B", 999.0))),
        ("chi_square", 999.0),
        ("materialization_identity", "forged-materialization"),
        ("evaluation_identity", "forged-evaluation"),
        ("output_scope", ("B", "A")),
        ("output_units", ("wrong", "wrong")),
        ("claims", ()),
        ("identity", "forged-success"),
    ):
        with pytest.raises(TypeError):
            dataclasses.replace(outcome.success, **{field_name: value})

    assert accepted.evaluation_result.resolved_values.ordered_items() == source_values


def test_summary_revalidates_success_after_post_construction_mutation() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    success = operation.evidence.outcomes[0].success
    assert success is not None
    altered_snapshot = ResolvedEvaluationValues(
        success.evaluation_result.resolved_values.parameterization_identity,
        success.evaluation_result.resolved_values.program_identity,
        (("A", 999.0), ("B", 999.0)),
    )
    object.__setattr__(
        success,
        "evaluation_result",
        dataclasses.replace(
            success.evaluation_result,
            resolved_values=altered_snapshot,
        ),
    )
    forged_outcome = dataclasses.replace(
        operation.evidence.outcomes[0], success=success
    )
    with pytest.raises(ResamplingConstructionError):
        dataclasses.replace(
            operation.evidence,
            outcomes=(forged_outcome, *operation.evidence.outcomes[1:]),
        )

    summary = summarize_resampling_evidence(operation.evidence)

    assert summary.terminal is SummaryTerminal.SOURCE_INVALID
    assert summary.summary is None
    assert summary.failure is not None
    assert summary.failure.category == "invalid_source_evidence"


def test_validated_success_record_round_trip_recomputes_every_projection() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    success = operation.evidence.outcomes[0].success
    assert success is not None

    restored = ValidatedReplicateSuccess.from_record(
        success.to_record(),
        plan,
        plan.replicates[0],
    )
    restored_outcome = dataclasses.replace(
        operation.evidence.outcomes[0], success=restored
    )
    restored_evidence = dataclasses.replace(
        operation.evidence,
        outcomes=(restored_outcome, *operation.evidence.outcomes[1:]),
    )

    assert restored.to_record() == success.to_record()
    assert summarize_resampling_evidence(restored_evidence).terminal is (
        SummaryTerminal.COMPLETED
    )


def test_validated_success_record_rejects_tampered_checked_projections() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    success = operation.evidence.outcomes[0].success
    assert success is not None

    def reject(path: tuple[object, ...], value: object) -> None:
        record = copy.deepcopy(success.to_record())
        target: Any = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        record["identity"] = hashlib.sha256(
            json.dumps(record, sort_keys=True).encode()
        ).hexdigest()
        with pytest.raises(ResamplingConstructionError):
            ValidatedReplicateSuccess.from_record(
                record,
                plan,
                plan.replicates[0],
            )

    reject(("resolved_items", 0, 1), (999.0).hex())
    reject(("chi_square",), (999.0).hex())
    reject(("materialization_identity",), "forged-materialization")
    reject(("evaluation_result", "identity"), "forged-evaluation")
    reject(("output_scope", 0), "B")
    reject(("output_units", 0), "wrong-unit")
    reject(("claims", 0, "state"), "violated")
    reject(("projected", "problem_identity"), "foreign-problem")


def test_success_record_rejects_rehashed_altered_complete_resolved_snapshot() -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    success = operation.evidence.outcomes[0].success
    assert success is not None
    altered_snapshot = ResolvedEvaluationValues(
        success.evaluation_result.resolved_values.parameterization_identity,
        success.evaluation_result.resolved_values.program_identity,
        (("A", 999.0), ("B", 999.0)),
    )
    altered_evaluation = dataclasses.replace(
        success.evaluation_result,
        resolved_values=altered_snapshot,
    )
    record = copy.deepcopy(success.to_record())
    record["evaluation_result"] = altered_evaluation.to_record()
    record["resolved_items"] = [["A", (999.0).hex()], ["B", (999.0).hex()]]
    record["identity"] = hashlib.sha256(
        json.dumps(record, sort_keys=True).encode()
    ).hexdigest()

    with pytest.raises(ResamplingConstructionError):
        ValidatedReplicateSuccess.from_record(
            record,
            plan,
            plan.replicates[0],
        )


def test_equivalent_looking_manual_success_without_exact_derivation_fails_closed() -> (
    None
):
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO, count=2)
    operation = execute_resampling_evidence(accepted, plan)
    assert operation.evidence is not None
    first = operation.evidence.outcomes[0]
    second = operation.evidence.outcomes[1]
    assert first.success is not None
    assert second.success is not None
    manual = object.__new__(ValidatedReplicateSuccess)
    object.__setattr__(manual, "prepared", first.success.prepared)
    object.__setattr__(manual, "projected", second.success.projected)
    object.__setattr__(manual, "evaluation_result", first.success.evaluation_result)
    forged_outcome = dataclasses.replace(first, success=manual)

    with pytest.raises(ResamplingConstructionError):
        dataclasses.replace(
            operation.evidence,
            outcomes=(forged_outcome, *operation.evidence.outcomes[1:]),
        )


@pytest.mark.parametrize(
    ("failure", "expected_stage", "draw_present"),
    (
        (KeyboardInterrupt(), ReplicateStage.GENERATING, False),
        (RuntimeError("generation failed"), ReplicateStage.GENERATING, False),
    ),
)
def test_pre_draw_generation_terminal_is_typed_and_optimizer_is_not_created(
    monkeypatch: pytest.MonkeyPatch,
    failure: BaseException,
    expected_stage: ReplicateStage,
    draw_present: bool,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)

    def fail_generation(
        _dataset: ResamplingDatasetManifest,
        _request: ReplicateRequest,
    ) -> None:
        raise failure

    monkeypatch.setattr(
        "chemex.optimize.native_resampling.generate_resampling_draw",
        fail_generation,
    )
    operation = execute_resampling_evidence(accepted, plan)

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


def test_interruption_after_draw_before_executor_creation_is_typed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan(ResamplingScheme.MONTE_CARLO)
    original = EvaluationEngine.resampled
    calls = 0

    def interrupt_binding(
        engine: EvaluationEngine,
        evaluation_plan: EvaluationPlan,
        bindings: object,
    ) -> EvaluationEngine:
        nonlocal calls
        calls += 1
        if calls == 1:
            raise KeyboardInterrupt
        return original(engine, evaluation_plan, cast("Any", bindings))

    monkeypatch.setattr(EvaluationEngine, "resampled", interrupt_binding)
    operation = execute_resampling_evidence(accepted, plan)

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

    def hook(
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        if prepared.request.ordinal == 0:
            raise KeyboardInterrupt
        return projected

    operation = execute_resampling_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=workers),
        candidate_test_hook=hook,
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

    def hook(
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationResult,
    ) -> ProjectedOptimizationResult:
        if prepared.request.ordinal != 1:
            return projected
        return ProjectedOptimizationFailure(
            prepared.draw.identity,
            ReplicateDisposition.FAILED,
            "solver_unsuccessful",
            "declared local budget exhausted",
            optimization_projection_identity=(
                prepared.request.optimization_projection_identity
            ),
            stage=ReplicateStage.EXECUTING,
            prepared_replicate_identity=prepared.identity,
            evaluation_plan_identity=prepared.evaluation_plan.identity,
            problem_identity=prepared.problem.identity,
            invocation_identity=prepared.invocation_identity,
        )

    operation = execute_resampling_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=2),
        candidate_test_hook=hook,
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
