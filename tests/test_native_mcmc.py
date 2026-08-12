"""Behavioral qualification tests for native fixed-topology MCMC evidence (#600)."""

from __future__ import annotations

import dataclasses
from copy import deepcopy
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pydantic import BaseModel

import chemex.optimize.native_mcmc as native_mcmc
from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.nmr.spectrometer import Spectrometer
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    OptimizationProblem,
)
from chemex.optimize.native_mcmc import (
    CalibratedMcmcPolicy,
    CalibrationCandidateMcmcPolicy,
    EnsembleState,
    ExpertMcmcPolicy,
    InitializationKind,
    McmcCalibrationReference,
    McmcConstructionError,
    McmcDiagnosticReason,
    McmcDiagnosticStatus,
    McmcEvidence,
    McmcEvidenceLifecycle,
    McmcExecutionStage,
    McmcInitializationOutcome,
    McmcOperationTerminal,
    McmcPlan,
    McmcPolicyKind,
    McmcTrajectoryClaim,
    PosteriorSampleDisposition,
    ProposalKind,
    ResolvedMcmcPolicy,
    build_bounded_latin_hypercube,
    derive_mcmc_diagnostics,
    derive_mcmc_operation_diagnostics,
    derive_posterior_sample_evidence,
    derive_posterior_summary,
    derive_retained_sample_view,
    execute_mcmc_evidence,
)
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ConstraintProgram,
    ParameterRole,
    ScientificFunctionBinder,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.printers.data import Printer
from chemex.typing import Array


class _KernelSettings(BaseModel):
    kind: str = "linear-test-kernel"


class _LinearSpectrometer:
    def __init__(self) -> None:
        self.spin_system = SpinSystem.from_name("1N")
        self.values = {"a": 0.0, "b": 0.0}

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _LinearSpectrometer:
        return deepcopy(self)

    def native_kernel_descriptor(self) -> dict[str, object]:
        return {"kind": "linear-test-spectrometer"}


class _LinearPulseSequence:
    settings = _KernelSettings()

    def calculate(self, spectrometer: _LinearSpectrometer, data: Data) -> Array:
        return spectrometer.values["a"] + spectrometer.values["b"] * np.asarray(
            data.metadata,
            dtype=np.float64,
        )

    def is_reference(self, metadata: Array) -> Array:
        return metadata < 0.0


def _native_context(
    *, fit_b: bool = True
) -> tuple[
    AcceptedFitResult,
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
        (
            ("A", ParameterRole.FIT),
            ("B", ParameterRole.FIT if fit_b else ParameterRole.FIX),
        ),
    )
    snapshot = AnalysisValuesSnapshot(
        "source-occurrence",
        "model",
        "definitions",
        "configuration",
        4,
        (("A", 0.5), ("B", 1.5)),
    )
    data = Data(
        exp=np.asarray((1.0, 3.0, 5.0), dtype=np.float64),
        err=np.ones(3, dtype=np.float64),
        metadata=np.asarray((0.0, 1.0, 2.0), dtype=np.float64),
    )
    profile = Profile(
        data,
        cast("Spectrometer", _LinearSpectrometer()),
        cast("PulseSequence", _LinearPulseSequence()),
        {"a": "A", "b": "B"},
        cast("Printer", None),
        is_scaled=False,
    )
    experiments = cast("Any", (SimpleNamespace(profiles=(profile,)),))
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    controlled_ids = ("A", "B") if fit_b else ("A",)
    held_items = () if fit_b else (("B", 1.5),)
    start = (0.5, 1.5) if fit_b else (0.5,)
    lower = (0.0, 1.0) if fit_b else (0.0,)
    upper = (2.0, 3.0) if fit_b else (2.0,)
    accepted_vector = (1.0, 2.0) if fit_b else (1.0,)
    problem = OptimizationProblem(
        engine.plan.identity,
        parameterization.identity,
        parameterization.evaluator_identity,
        program.fingerprint,
        "configuration",
        snapshot,
        (("A", 0.5), ("B", 1.5)),
        controlled_ids,
        held_items,
        start,
        lower,
        upper,
        ("A", "B"),
    )
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame(accepted_vector, parameterization),
    )
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
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
        vector=accepted_vector,
        chi_square=float(np.sum(result.residuals**2)),
        evaluation_result=result,
        commit_scope=problem.commit_scope,
        commit_items=result.resolved_values.ordered_items(),
        origin_context_identity="accepted-origin",
    )
    return accepted, problem, parameterization, engine


def _resolved_policy() -> ResolvedMcmcPolicy:
    return ExpertMcmcPolicy(
        burn_steps=2,
        retained_steps=3,
        walkers=8,
        expert_provenance="qualification-test-expert",
    ).resolve(dimension=2, root_seed=1234)


def _same_semantics_foreign_occurrence(
    accepted: AcceptedFitResult,
) -> AcceptedFitResult:
    return AcceptedFitResult.for_qualification(
        occurrence_identity="foreign-accepted-occurrence",
        problem_identity=accepted.problem_identity,
        invocation_identity=accepted.invocation_identity,
        execution_identity=accepted.execution_identity,
        materialization_identity=accepted.materialization_identity,
        parameterization_identity=accepted.parameterization_identity,
        evaluator_parameterization_identity=(
            accepted.evaluator_parameterization_identity
        ),
        source_occurrence_identity=accepted.source_occurrence_identity,
        source_revision=accepted.source_revision,
        controlled_ids=accepted.controlled_ids,
        vector=accepted.vector,
        chi_square=accepted.chi_square,
        evaluation_result=accepted.evaluation_result,
        commit_scope=accepted.commit_scope,
        commit_items=accepted.commit_items,
        origin_context_identity=accepted.origin_context_identity,
    )


def test_arbitrary_settings_cannot_mint_calibrated_authority() -> None:
    with pytest.raises(McmcConstructionError, match="repository-frozen calibration"):
        ResolvedMcmcPolicy(
            kind=McmcPolicyKind.CALIBRATED,
            policy_version="caller-labelled-calibrated",
            dimension=2,
            walkers=12,
            burn_steps=20,
            retained_steps=40,
            root_seed=1234,
            provenance_identity="caller-labelled-settings",
            qualification_dimension_range=(1, 3),
        )


def test_foreign_calibration_reference_cannot_authorize_policy() -> None:
    reference = McmcCalibrationReference(
        calibration_identity="future-calibration-release",
        baseline_release_identity="future-baseline-release",
        numerical_lane_requirement="canonical-python-3.13-x86-64",
        policy_version="bounded-analytic-v1",
        minimum_dimension=1,
        maximum_dimension=3,
        walkers=12,
        burn_steps=20,
        retained_steps=40,
    )
    policy = CalibratedMcmcPolicy(reference, authority=object())

    with pytest.raises(McmcConstructionError, match="calibration reference"):
        policy.resolve(dimension=2, root_seed=0x1234_5678_90AB_CDEF)


def test_provisional_calibration_candidate_is_explicitly_unqualified() -> None:
    policy = CalibrationCandidateMcmcPolicy(
        candidate_identity="candidate-bounded-analytic-v1",
        minimum_dimension=1,
        maximum_dimension=3,
        walkers=12,
        burn_steps=20,
        retained_steps=40,
    )

    resolved = policy.resolve(dimension=2, root_seed=1234)

    assert resolved.kind is McmcPolicyKind.CALIBRATION_CANDIDATE
    assert not resolved.has_calibrated_adequacy
    assert resolved.qualification_dimension_range == (1, 3)
    assert resolved.dimension == 2
    assert resolved.walkers == 12
    assert resolved.burn_steps == 20
    assert resolved.retained_steps == 40
    assert resolved.total_steps == 60
    assert resolved.thin == 1
    assert resolved.initialization is InitializationKind.BOUNDED_LATIN_HYPERCUBE
    assert resolved.proposal is ProposalKind.STRETCH
    assert resolved.proposal_scale == 2.0
    assert resolved.objective_request_budget == 12 * 61
    assert resolved.root_seed == 1234
    assert resolved.initialization_seed != resolved.sampler_seed
    assert resolved.identity


def test_expert_policy_only_overrides_predeclared_topology() -> None:
    policy = ExpertMcmcPolicy(
        burn_steps=7,
        retained_steps=13,
        walkers=8,
        expert_provenance="user-declared-expert-settings",
    )

    resolved = policy.resolve(dimension=2, root_seed=99)

    assert resolved.kind is McmcPolicyKind.EXPERT
    assert resolved.policy_version == "expert-v1"
    assert resolved.provenance_identity == "user-declared-expert-settings"
    assert not resolved.has_calibrated_adequacy
    assert resolved.qualification_dimension_range is None
    assert resolved.burn_steps == 7
    assert resolved.retained_steps == 13
    assert resolved.walkers == 8
    assert resolved.initialization is InitializationKind.BOUNDED_LATIN_HYPERCUBE
    assert resolved.proposal is ProposalKind.STRETCH
    assert resolved.thin == 1
    assert resolved.objective_request_budget == 8 * 21


def test_calibration_candidate_fails_closed_outside_its_declared_stratum() -> None:
    policy = CalibrationCandidateMcmcPolicy(
        candidate_identity="bounded-analytic-candidate-v1",
        minimum_dimension=1,
        maximum_dimension=2,
        walkers=8,
        burn_steps=10,
        retained_steps=20,
    )

    with pytest.raises(McmcConstructionError, match="does not cover"):
        policy.resolve(dimension=3, root_seed=1)


def test_initial_ensemble_is_a_seeded_bounded_latin_hypercube() -> None:
    lower = (-2.0, 10.0)
    upper = (2.0, 22.0)

    first = build_bounded_latin_hypercube(
        lower,
        upper,
        walkers=8,
        seed=0xCAFE,
    )
    second = build_bounded_latin_hypercube(
        lower,
        upper,
        walkers=8,
        seed=0xCAFE,
    )
    positions = np.asarray(first)

    assert first == second
    assert positions.shape == (8, 2)
    assert np.all(positions > np.asarray(lower))
    assert np.all(positions < np.asarray(upper))
    strata = np.floor(
        (positions - np.asarray(lower)) / (np.asarray(upper) - np.asarray(lower)) * 8
    ).astype(int)
    for dimension in range(2):
        assert sorted(strata[:, dimension]) == list(range(8))


@pytest.mark.parametrize(
    ("lower", "upper"),
    (
        ((-np.inf,), (1.0,)),
        ((0.0,), (np.inf,)),
        ((1.0,), (1.0,)),
        ((2.0,), (1.0,)),
    ),
)
def test_initial_ensemble_rejects_nonfinite_or_empty_bounds(
    lower: tuple[float, ...],
    upper: tuple[float, ...],
) -> None:
    with pytest.raises(McmcConstructionError, match="finite strictly ordered bounds"):
        build_bounded_latin_hypercube(lower, upper, walkers=4, seed=1)


@pytest.mark.parametrize("seed", (True, -1, 1 << 64, 1.5, np.nan))
def test_initial_ensemble_rejects_noncanonical_seed(seed: object) -> None:
    with pytest.raises(McmcConstructionError, match="unsigned 64-bit seed"):
        build_bounded_latin_hypercube((0.0,), (1.0,), walkers=4, seed=cast("int", seed))


def test_initial_ensemble_rejects_transposed_or_unrepresentably_narrow_bounds() -> None:
    with pytest.raises(McmcConstructionError, match="finite strictly ordered bounds"):
        build_bounded_latin_hypercube(
            cast("tuple[float, ...]", ((0.0, 1.0),)),
            cast("tuple[float, ...]", ((2.0, 3.0),)),
            walkers=4,
            seed=1,
        )
    lower = 1.0
    upper = float(np.nextafter(lower, np.inf))
    with pytest.raises(McmcConstructionError, match="closed bound"):
        build_bounded_latin_hypercube((lower,), (upper,), walkers=4, seed=1)


def test_plan_binds_initial_ensemble_to_exact_accepted_native_lineage() -> None:
    accepted, problem, parameterization, engine = _native_context()

    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    assert plan.accepted_result_identity == accepted.identity
    assert plan.anchor.accepted_occurrence_identity == accepted.occurrence_identity
    assert (
        plan.anchor.source_occurrence_identity
        == problem.source_snapshot.occurrence_identity
    )
    assert plan.anchor.source_revision == problem.source_snapshot.revision
    assert plan.anchor.evaluation_plan_identity == engine.plan.identity
    assert plan.anchor.parameterization_identity == parameterization.identity
    assert plan.anchor.held_items == problem.held_items
    assert (
        plan.anchor.affine_feasibility_identity == problem.affine_feasibility_identity
    )
    assert accepted.vector != problem.start
    assert plan.problem_identity == problem.identity
    assert plan.coordinate_ids == ("A", "B")
    assert plan.lower_bounds == (0.0, 1.0)
    assert plan.upper_bounds == (2.0, 3.0)
    assert len(plan.initial_ensemble) == plan.policy.walkers
    assert all(
        lower < value < upper
        for row in plan.initial_ensemble
        for value, lower, upper in zip(
            row,
            plan.lower_bounds,
            plan.upper_bounds,
            strict=True,
        )
    )
    assert plan.identity


def test_plan_rejects_foreign_problem_parameterization_and_coordinate_order() -> None:
    accepted, problem, parameterization, engine = _native_context()
    policy = _resolved_policy()
    units = (
        ("A", ParameterUnit.DIMENSIONLESS),
        ("B", ParameterUnit.DIMENSIONLESS),
    )
    foreign_problem = dataclasses.replace(problem, start=(0.75, 1.5))
    foreign_parameterization = dataclasses.replace(parameterization, source_revision=5)

    with pytest.raises(McmcConstructionError, match="exact accepted native lineage"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=foreign_problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=units,
        )
    with pytest.raises(McmcConstructionError, match="exact accepted native lineage"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=problem,
            parameterization=foreign_parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=units,
        )
    with pytest.raises(McmcConstructionError, match="coordinate units"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=tuple(reversed(units)),
        )


def test_complete_chain_is_primary_and_flat_samples_are_derived_views() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    source_vector = accepted.vector
    source_items = tuple(problem.source_snapshot.items())

    operation = execute_mcmc_evidence(accepted, plan)
    replay = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    evidence = operation.evidence
    assert replay.evidence is not None
    assert replay.evidence.identity == evidence.identity
    assert evidence.lifecycle is McmcEvidenceLifecycle.COMPLETED
    assert evidence.plan_identity == plan.identity
    assert evidence.accepted_result_identity == accepted.identity
    assert evidence.completed_transition_count == plan.policy.total_steps
    assert tuple(state.ordinal for state in evidence.states) == tuple(
        range(plan.policy.total_steps + 1)
    )
    assert evidence.states[0].positions == plan.initial_ensemble
    assert all(len(state.positions) == plan.policy.walkers for state in evidence.states)
    assert all(
        len(position) == plan.policy.dimension
        for state in evidence.states
        for position in state.positions
    )
    assert all(
        len(state.log_densities) == plan.policy.walkers for state in evidence.states
    )
    assert evidence.objective_request_count == plan.policy.objective_request_budget
    assert accepted.vector == source_vector
    assert tuple(problem.source_snapshot.items()) == source_items

    initial = np.asarray(evidence.states[0].positions)
    metadata = np.asarray((0.0, 1.0, 2.0))
    observed = np.asarray((1.0, 3.0, 5.0))
    expected_initial_log_density = -0.5 * np.sum(
        (initial[:, 0, np.newaxis] + initial[:, 1, np.newaxis] * metadata - observed)
        ** 2,
        axis=1,
    )
    # These small, well-conditioned binary64 calculations use only affine
    # arithmetic and a three-term sum. 1e-14 is roughly 45 machine eps at unit
    # scale: enough for operation-order differences while still detecting a
    # changed seeded trajectory or log-density convention.
    np.testing.assert_allclose(
        evidence.states[0].log_densities,
        expected_initial_log_density,
        rtol=0.0,
        atol=1.0e-14,
    )
    # A local seeded replay is ordinary capture evidence only. It is not a
    # canonical trajectory claim without the #588 live-lane authority and an
    # independently frozen reference.
    assert replay.evidence.states == evidence.states
    assert evidence.trajectory_claim is McmcTrajectoryClaim.ORDINARY_CAPTURE
    assert not evidence.canonical_lane_qualified

    retained = derive_retained_sample_view(evidence)
    assert retained.source_evidence_identity == evidence.identity
    assert retained.coordinate_ids == plan.coordinate_ids
    assert retained.selected_state_ordinals == (3, 4, 5)
    assert retained.sample_indices[0] == (3, 0)
    assert retained.sample_indices[-1] == (5, 7)
    assert len(retained.samples) == plan.policy.walkers * plan.policy.retained_steps
    assert len(retained.log_densities) == len(retained.samples)
    assert retained.is_complete

    diagnostics = derive_mcmc_diagnostics(evidence)
    assert diagnostics.source_evidence_identity == evidence.identity
    assert diagnostics.state_ordinals == (1, 2, 3, 4, 5)
    assert diagnostics.walker_ordinals == tuple(range(plan.policy.walkers))
    assert diagnostics.acceptance_fractions is not None
    assert diagnostics.mean_acceptance_fraction is not None
    assert len(diagnostics.acceptance_fractions) == plan.policy.walkers
    assert 0.0 <= diagnostics.mean_acceptance_fraction <= 1.0


def test_interruption_preserves_only_a_contiguous_complete_state_prefix() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    def interrupt_after_second_transition(state: EnsembleState) -> None:
        if state.ordinal == 2:
            raise KeyboardInterrupt

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        state_observer=interrupt_after_second_transition,
    )

    assert operation.terminal is McmcOperationTerminal.INTERRUPTED
    assert operation.failure_category == "interrupted"
    assert operation.evidence is not None
    evidence = operation.evidence
    assert evidence.lifecycle is McmcEvidenceLifecycle.PARTIAL
    assert tuple(state.ordinal for state in evidence.states) == (0, 1, 2)
    assert evidence.completed_transition_count == 2
    assert evidence.objective_request_count == plan.policy.walkers * 3
    retained = derive_retained_sample_view(evidence)
    assert retained.selected_state_ordinals == ()
    assert not retained.is_complete


def test_primary_chain_evidence_round_trips_and_rejects_tampering() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    evidence = operation.evidence

    record = evidence.to_record()
    assert operation.raw_capture is not None
    restored = McmcEvidence.from_record(
        record,
        plan=plan,
        raw_capture=operation.raw_capture,
    )

    assert restored.identity == evidence.identity
    assert restored.states == evidence.states
    assert restored.objective_request_count == evidence.objective_request_count

    tampered = deepcopy(record)
    states = cast("list[dict[str, object]]", tampered["states"])
    positions = cast("list[list[str]]", states[1]["positions"])
    positions[0][0] = float(float.fromhex(positions[0][0]) + 0.1).hex()
    with pytest.raises(McmcConstructionError, match="identity"):
        McmcEvidence.from_record(
            tampered,
            plan=plan,
            raw_capture=operation.raw_capture,
        )


def test_execution_rejects_same_semantics_from_foreign_accepted_occurrence() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    foreign = _same_semantics_foreign_occurrence(accepted)
    assert foreign.identity == accepted.identity
    assert foreign.occurrence_identity != accepted.occurrence_identity

    with pytest.raises(McmcConstructionError, match="accepted occurrence"):
        execute_mcmc_evidence(foreign, plan)


def test_fresh_validation_rejects_tampered_backend_log_density() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    first = operation.raw_capture.states[0]
    tampered_first = dataclasses.replace(
        first,
        log_densities=(first.log_densities[0] + 1.0, *first.log_densities[1:]),
    )
    with pytest.raises(McmcConstructionError, match="transition"):
        dataclasses.replace(
            operation.raw_capture,
            states=(tampered_first, *operation.raw_capture.states[1:]),
        )


def test_fresh_validator_rejects_backend_originated_log_density_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan_context()
    original = native_mcmc._ensemble_state

    def wrong_backend_density(
        ordinal: int,
        positions: Array,
        log_densities: Array,
        accepted_mask: Array | None,
    ) -> EnsembleState:
        state = original(ordinal, positions, log_densities, accepted_mask)
        if ordinal != 0:
            return state
        return EnsembleState(
            state.ordinal,
            state.positions,
            (state.log_densities[0] + 1.0, *state.log_densities[1:]),
            state.accepted,
        )

    monkeypatch.setattr(native_mcmc, "_ensemble_state", wrong_backend_density)

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.FAILED
    assert operation.evidence is None
    assert operation.validation is not None
    assert operation.validation.failures[0].category == "backend_log_density_mismatch"


def _plan_context() -> tuple[AcceptedFitResult, McmcPlan]:
    accepted, problem, parameterization, engine = _native_context()
    return accepted, McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )


def test_cancellation_before_initialization_freezes_zero_state_operation() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()
    cancellation.cancel()

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
    )

    assert operation.terminal is McmcOperationTerminal.CANCELLED
    assert operation.complete_state_count == 0
    assert operation.stage is McmcExecutionStage.BEFORE_INITIALIZATION
    assert operation.evidence is None
    assert operation.raw_capture is not None
    assert operation.raw_capture.policy_identity == plan.policy.identity
    assert operation.raw_capture.root_seed == plan.policy.root_seed
    assert (
        operation.raw_capture.initialization_outcome
        is McmcInitializationOutcome.CANCELLED
    )
    diagnostics = derive_mcmc_operation_diagnostics(operation)
    assert diagnostics.status is McmcDiagnosticStatus.UNAVAILABLE
    assert diagnostics.reason is McmcDiagnosticReason.NO_TRANSITIONS


@pytest.mark.parametrize(
    ("target_stage", "terminal", "category"),
    (
        (
            McmcExecutionStage.INITIALIZING,
            McmcOperationTerminal.INTERRUPTED,
            "interrupted",
        ),
        (
            McmcExecutionStage.BEFORE_TRANSITION,
            McmcOperationTerminal.CANCELLED,
            "cancelled",
        ),
        (
            McmcExecutionStage.DURING_TRANSITION,
            McmcOperationTerminal.INTERRUPTED,
            "interrupted",
        ),
    ),
)
def test_zero_or_initial_state_terminal_checkpoints_are_truthful(
    target_stage: McmcExecutionStage,
    terminal: McmcOperationTerminal,
    category: str,
) -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def stop(stage: McmcExecutionStage, _count: int) -> None:
        if stage is not target_stage:
            return
        if terminal is McmcOperationTerminal.CANCELLED:
            cancellation.cancel()
        else:
            raise KeyboardInterrupt

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=stop,
    )

    assert operation.terminal is terminal
    assert operation.failure_category == category
    assert operation.stage is target_stage
    expected_count = 0 if target_stage is McmcExecutionStage.INITIALIZING else 1
    assert operation.complete_state_count == expected_count


def test_initialization_failure_and_later_failure_preserve_only_complete_prefix() -> (
    None
):
    accepted, plan = _plan_context()

    def fail_initialization(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.INITIALIZING:
            raise RuntimeError("initialization failure")

    initial_failure = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=fail_initialization,
    )
    assert initial_failure.terminal is McmcOperationTerminal.FAILED
    assert initial_failure.complete_state_count == 0
    assert initial_failure.evidence is None
    assert initial_failure.raw_capture is not None
    assert (
        initial_failure.raw_capture.initialization_outcome
        is McmcInitializationOutcome.FAILED
    )

    def fail_after_two_transitions(stage: McmcExecutionStage, count: int) -> None:
        if stage is McmcExecutionStage.AFTER_COMPLETE_STATE and count == 3:
            raise RuntimeError("later failure")

    later_failure = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=fail_after_two_transitions,
    )
    assert later_failure.terminal is McmcOperationTerminal.FAILED
    assert later_failure.complete_state_count == 3
    assert later_failure.evidence is not None
    assert tuple(state.ordinal for state in later_failure.evidence.states) == (0, 1, 2)


def test_cancellation_after_complete_state_and_before_final_assembly() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def cancel_after_one_transition(stage: McmcExecutionStage, count: int) -> None:
        if stage is McmcExecutionStage.AFTER_COMPLETE_STATE and count == 2:
            cancellation.cancel()

    partial = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=cancel_after_one_transition,
    )
    assert partial.terminal is McmcOperationTerminal.CANCELLED
    assert partial.complete_state_count == 2
    assert partial.evidence is not None
    assert tuple(state.ordinal for state in partial.evidence.states) == (0, 1)

    final_cancellation = CancellationToken()

    def cancel_final(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.BEFORE_FINAL_ASSEMBLY:
            final_cancellation.cancel()

    final = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=final_cancellation,
        checkpoint_observer=cancel_final,
    )
    assert final.terminal is McmcOperationTerminal.CANCELLED
    assert final.complete_state_count == plan.policy.total_steps + 1
    assert final.evidence is not None
    assert final.evidence.lifecycle is McmcEvidenceLifecycle.PARTIAL


def test_zero_transition_acceptance_is_typed_unavailable() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def cancel_before_transition(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.BEFORE_TRANSITION:
            cancellation.cancel()

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=cancel_before_transition,
    )
    assert operation.evidence is not None

    diagnostics = derive_mcmc_diagnostics(operation.evidence)

    assert diagnostics.status is McmcDiagnosticStatus.UNAVAILABLE
    assert diagnostics.reason is McmcDiagnosticReason.NO_TRANSITIONS
    assert diagnostics.acceptance_fractions is None
    assert diagnostics.mean_acceptance_fraction is None


def test_complete_scope_posterior_outcomes_and_summary_are_canonical() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)

    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    summary = derive_posterior_summary(posterior)

    assert (
        tuple((item.state_ordinal, item.walker_ordinal) for item in posterior.outcomes)
        == selection.sample_indices
    )
    assert all(
        item.disposition is PosteriorSampleDisposition.SUCCESS
        and item.resolved_items is not None
        and tuple(key for key, _value in item.resolved_items) == ("A", "B")
        for item in posterior.outcomes
    )
    assert summary.included_labels == selection.sample_indices
    assert not summary.excluded_labels
    assert len(summary.parameter_summaries) == 2
    assert all(
        item.credible_interval_name == "equal_tailed_95_percent"
        and item.posterior_standard_deviation >= 0.0
        and item.autocorrelation_time.status is McmcDiagnosticStatus.UNAVAILABLE
        for item in summary.parameter_summaries
    )
    assert len(summary.covariance) == 4
    assert len(summary.correlations) == 4

    with pytest.raises(McmcConstructionError, match="witness"):
        dataclasses.replace(summary, included_labels=())


def test_posterior_samples_combine_exact_held_values_before_scope_projection() -> None:
    accepted, problem, parameterization, engine = _native_context(fit_b=False)
    policy = ExpertMcmcPolicy(
        burn_steps=1,
        retained_steps=2,
        walkers=4,
        expert_provenance="held-value-qualification",
    ).resolve(dimension=1, root_seed=9876)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(("A", ParameterUnit.DIMENSIONLESS),),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)

    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    assert all(
        outcome.resolved_items is not None and dict(outcome.resolved_items)["B"] == 1.5
        for outcome in posterior.outcomes
    )


def test_primary_tampering_cannot_reach_any_derived_posterior_artifact() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    evidence = operation.evidence
    with pytest.raises(McmcConstructionError, match="contiguous"):
        dataclasses.replace(
            evidence,
            states=(evidence.states[0], *reversed(evidence.states[1:])),
        )


def test_all_operation_and_posterior_paths_preserve_source_values_bitwise() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    before = problem.source_snapshot.to_json()
    completed = execute_mcmc_evidence(accepted, plan)
    cancellation = CancellationToken()
    cancellation.cancel()
    cancelled = execute_mcmc_evidence(accepted, plan, cancellation=cancellation)

    def interrupt(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.INITIALIZING:
            raise KeyboardInterrupt

    interrupted = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=interrupt,
    )
    assert completed.evidence is not None
    selection = derive_retained_sample_view(completed.evidence)
    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    derive_posterior_summary(posterior)

    assert cancelled.terminal is McmcOperationTerminal.CANCELLED
    assert interrupted.terminal is McmcOperationTerminal.INTERRUPTED
    assert problem.source_snapshot.to_json() == before
    assert accepted.vector == (1.0, 2.0)
