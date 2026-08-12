"""Behavioral qualification tests for native fixed-topology MCMC evidence (#600)."""

from __future__ import annotations

from copy import deepcopy
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pydantic import BaseModel

from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.nmr.spectrometer import Spectrometer
from chemex.optimize.direct_trf import AcceptedFitResult, OptimizationProblem
from chemex.optimize.native_mcmc import (
    CalibratedMcmcPolicy,
    EnsembleState,
    ExpertMcmcPolicy,
    InitializationKind,
    McmcConstructionError,
    McmcEvidence,
    McmcEvidenceLifecycle,
    McmcOperationTerminal,
    McmcPlan,
    McmcPolicyKind,
    ProposalKind,
    ResolvedMcmcPolicy,
    build_bounded_latin_hypercube,
    derive_mcmc_diagnostics,
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


def _native_context() -> tuple[
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
        (("A", ParameterRole.FIT), ("B", ParameterRole.FIT)),
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
    problem = OptimizationProblem(
        engine.plan.identity,
        parameterization.identity,
        parameterization.evaluator_identity,
        program.fingerprint,
        "configuration",
        snapshot,
        (("A", 0.5), ("B", 1.5)),
        ("A", "B"),
        (),
        (0.5, 1.5),
        (0.0, 1.0),
        (2.0, 3.0),
        ("A", "B"),
    )
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame((1.0, 2.0), parameterization),
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
        vector=(1.0, 2.0),
        chi_square=0.0,
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
    ).resolve(dimension=2, root_seed=1234)


def test_calibrated_policy_resolves_the_complete_prospective_run() -> None:
    policy = CalibratedMcmcPolicy(
        policy_version="bounded-analytic-v1",
        minimum_dimension=1,
        maximum_dimension=3,
        walkers=12,
        burn_steps=20,
        retained_steps=40,
    )

    resolved = policy.resolve(dimension=2, root_seed=0x1234_5678_90AB_CDEF)

    assert resolved.kind is McmcPolicyKind.CALIBRATED
    assert resolved.policy_version == "bounded-analytic-v1"
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
    assert resolved.root_seed == 0x1234_5678_90AB_CDEF
    assert resolved.initialization_seed != resolved.sampler_seed
    assert resolved.identity


def test_expert_policy_only_overrides_predeclared_topology() -> None:
    policy = ExpertMcmcPolicy(
        burn_steps=7,
        retained_steps=13,
        walkers=8,
    )

    resolved = policy.resolve(dimension=2, root_seed=99)

    assert resolved.kind is McmcPolicyKind.EXPERT
    assert resolved.policy_version == "expert-v1"
    assert resolved.qualification_dimension_range is None
    assert resolved.burn_steps == 7
    assert resolved.retained_steps == 13
    assert resolved.walkers == 8
    assert resolved.initialization is InitializationKind.BOUNDED_LATIN_HYPERCUBE
    assert resolved.proposal is ProposalKind.STRETCH
    assert resolved.thin == 1
    assert resolved.objective_request_budget == 8 * 21


def test_calibrated_policy_fails_closed_outside_its_qualified_stratum() -> None:
    policy = CalibratedMcmcPolicy(
        policy_version="bounded-analytic-v1",
        minimum_dimension=1,
        maximum_dimension=2,
        walkers=8,
        burn_steps=10,
        retained_steps=20,
    )

    with pytest.raises(McmcConstructionError, match="No calibrated MCMC policy"):
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
    np.testing.assert_allclose(
        evidence.states[-1].positions,
        (
            (0.22193955880807614, 1.2886666343140305),
            (0.9121794119730182, 2.2579167141311447),
            (1.8830633726465182, 1.5610464161494841),
            (1.383428732177615, 1.1821818901130026),
            (1.1957090816652476, 1.6603544159276828),
            (1.2211689734419329, 2.26568078108725),
            (1.1292648992933922, 1.9097405704768293),
            (0.029980653967193294, 2.8508821421336843),
        ),
        rtol=0.0,
        atol=1.0e-14,
    )

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
    restored = McmcEvidence.from_record(record, plan=plan)

    assert restored.identity == evidence.identity
    assert restored.states == evidence.states
    assert restored.objective_request_count == evidence.objective_request_count

    tampered = deepcopy(record)
    states = cast("list[dict[str, object]]", tampered["states"])
    positions = cast("list[list[str]]", states[1]["positions"])
    positions[0][0] = float(float.fromhex(positions[0][0]) + 0.1).hex()
    with pytest.raises(McmcConstructionError, match="identity"):
        McmcEvidence.from_record(tampered, plan=plan)
