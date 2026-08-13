"""Behavioral qualification for the closed native method-step composer (#641).

The public seam is an immutable ``MethodStepWorkflow`` executed against the
exact session-owned AnalysisValues occurrence.  Tests observe only its typed
outcome and the explicit successor-state gate.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace
from typing import TypedDict
from unittest.mock import patch

import numpy as np
import pytest

from chemex.baselines import (
    LegacyObservationImplementation,
    Occurrence,
    ResultBundle,
    ResultMember,
)
from chemex.configuration.methods import (
    McmcSettings,
    Method,
    Selection,
    Statistics,
    read_methods,
)
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.native_provenance import (
    BaselineReference,
    ProvenanceEnvironment,
    WorkflowProvenance,
)
from chemex.optimize.de_direct_trf import (
    DeCoordinateSemantics,
    DeDirectTrfInvocation,
    DeDirectTrfOutcome,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    DirectTrfConstructionError,
    MaterializationTerminal,
    OptimizationProblem,
    RootMaterializationFailure,
    TerminalFailure,
)
from chemex.optimize.grid_direct_trf import GridDirectTrfInvocation
from chemex.optimize.grouped_direct_trf import (
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfOutcome,
)
from chemex.optimize.grouped_grid_direct_trf import GroupedGridDirectTrfOutcome
from chemex.optimize.method_step import (
    DerivationDisposition,
    EvaluationPurpose,
    McmcDerivationRequest,
    MethodStepCheckpoint,
    MethodStepLifecycle,
    MethodStepOutcome,
    MethodStepPublicationRequest,
    MethodStepRunProvenanceRequest,
    MethodStepStrategy,
    MethodStepWorkflow,
    ResamplingDerivationRequest,
    UncertaintyDerivationRequest,
    execute_method_step,
    require_successor_state,
)
from chemex.optimize.native_mcmc import ExpertMcmcPolicy, McmcEvidence
from chemex.optimize.native_resampling import (
    OperationTerminal as ResamplingOperationTerminal,
)
from chemex.optimize.native_resampling import (
    OptimizationStrategy,
    ResamplingEvidence,
    ResamplingOperation,
    ResamplingPlan,
    ResamplingScheme,
    ResamplingSummaryOutcome,
    ResamplingSummaryPolicy,
    SummaryFailure,
    SummaryTerminal,
    execute_resampling_evidence,
)
from chemex.optimize.uncertainty import (
    ParameterUnit,
    ResidualVarianceScaling,
    UncertaintyPolicy,
    compile_constraint_linearization_capabilities,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    SealedParameterModel,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.run_info import capture_native_inputs
from chemex.runtime import AnalysisSession, ExecutionSettings

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


class _OptimizationInputs(TypedDict):
    starting_snapshot: AnalysisValuesSnapshot
    parameter_model: SealedParameterModel
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    method: Method
    problem: OptimizationProblem
    decomposition: FitDecomposition


def _evaluation_workflow(
    *,
    purpose: EvaluationPurpose = EvaluationPurpose.NO_OPTIMIZATION_REQUIRED,
) -> tuple[AnalysisSession, MethodStepWorkflow]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = (
        Method(fit=["PB"], fix=["R1A_A", "KEX_AB"])
        if purpose is EvaluationPurpose.EVALUATE_ONLY
        else Method(fix=["R1A_A", "PB", "KEX_AB"])
    )
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    if purpose is EvaluationPurpose.NO_OBJECTIVE_DATA:
        for experiment in experiments:
            for profile in experiment.profiles:
                profile.is_scaled = False
                profile.data.mask = np.zeros(profile.data.size, dtype=np.bool_)
                profile.data.mark_dirty()
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    workflow = MethodStepWorkflow.for_evaluation(
        starting_snapshot=session.analysis_values.snapshot(),
        parameter_model=parameter_model,
        parameterization=parameterization,
        engine=engine,
        method=method,
        purpose=purpose,
    )
    return session, workflow


def _optimization_inputs(
    *,
    all_profiles: bool = False,
    method: Method | None = None,
    stabilize_data: bool = False,
) -> tuple[AnalysisSession, _OptimizationInputs]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(
            include=(None if all_profiles else [SpinSystem.from_name("G2N-HN")]),
            exclude=None,
        ),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = read_methods([METHOD])["DEFAULT"] if method is None else method
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    if stabilize_data:
        frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
        )
        initial = engine.new_evaluator().evaluate(frame)
        assert isinstance(initial, EvaluationResult)
        offset = 0
        for experiment in experiments:
            for profile in experiment.profiles:
                stop = offset + profile.data.size
                profile.data.exp = np.asarray(
                    initial.normalized_calculations[offset:stop],
                    dtype=np.float64,
                ).copy()
                profile.data.mark_dirty()
                offset = stop
        engine = EvaluationEngine.from_experiments(experiments, parameterization)
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    starting = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        starting,
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    return session, {
        "starting_snapshot": starting,
        "parameter_model": parameter_model,
        "parameterization": parameterization,
        "engine": engine,
        "method": method,
        "problem": problem,
        "decomposition": decomposition,
    }


def _direct_workflow(
    *,
    all_profiles: bool = False,
    stabilize_data: bool = False,
    derivations: tuple[
        UncertaintyDerivationRequest
        | ResamplingDerivationRequest
        | McmcDerivationRequest,
        ...,
    ] = (),
) -> tuple[AnalysisSession, MethodStepWorkflow]:
    session, inputs = _optimization_inputs(
        all_profiles=all_profiles,
        stabilize_data=stabilize_data,
    )
    decomposition = inputs["decomposition"]
    assert isinstance(decomposition, FitDecomposition)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    workflow = MethodStepWorkflow.for_optimization(
        **inputs,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=invocation,
        derivations=derivations,
    )
    return session, workflow


def _grid_workflow(
    *,
    all_profiles: bool = False,
) -> tuple[AnalysisSession, MethodStepWorkflow]:
    session, inputs = _optimization_inputs(all_profiles=all_profiles)
    problem = inputs["problem"]
    assert isinstance(problem, OptimizationProblem)
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((problem.controlled_ids[0], (problem.start[0],)),),
        objective_request_budget=80,
    )
    return session, MethodStepWorkflow.for_optimization(
        **inputs,
        strategy=MethodStepStrategy.GRID_DIRECT_TRF,
        invocation=invocation,
    )


def _de_workflow() -> tuple[AnalysisSession, MethodStepWorkflow]:
    session, inputs = _optimization_inputs()
    problem = inputs["problem"]
    assert isinstance(problem, OptimizationProblem)
    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=(
            (
                problem.controlled_ids[0],
                0.5 * problem.start[0],
                1.5 * problem.start[0],
                DeCoordinateSemantics.LINEAR,
            ),
        ),
        root_seed=641,
        de_objective_request_budget=30,
        polish_objective_request_budget=80,
        population_multiplier=4,
        maximum_generations=5,
    )
    return session, MethodStepWorkflow.for_optimization(
        **inputs,
        strategy=MethodStepStrategy.DE_DIRECT_TRF,
        invocation=invocation,
    )


def _uncertainty_request(workflow: MethodStepWorkflow) -> UncertaintyDerivationRequest:
    controlled_id = workflow.problem.controlled_ids[0] if workflow.problem else ""
    policy = UncertaintyPolicy(
        calibration_identity="issue-641-qualification-v1",
        numerical_compatibility_requirement=(
            "binary64-scipy-economical-svd-fixed-pairwise-v1"
        ),
        coordinate_scales=((controlled_id, 1.0),),
        coordinate_units=((controlled_id, ParameterUnit.RATE_PER_SECOND),),
        residual_variance_scaling=(
            ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE
        ),
        relative_step_tolerance=1.0e-4,
        roundoff_multiplier=64.0,
        smaller_step_extent=8,
        larger_step_extent=8,
        svd_driver="gesdd",
        rank_absolute_tolerance=0.0,
        rank_relative_tolerance=1.0e-12,
        weak_relative_tolerance=1.0e-6,
        singular_value_cluster_relative_tolerance=1.0e-10,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1",
    )
    capabilities = compile_constraint_linearization_capabilities(
        workflow.parameterization,
        (),
        (),
    )
    return UncertaintyDerivationRequest(
        policy,
        compiled_capabilities=capabilities,
        resolved_environment_identity="issue-641-local-environment",
    )


def _resampling_request(workflow: MethodStepWorkflow) -> ResamplingDerivationRequest:
    size = workflow.engine.plan.observation_count
    scope = workflow.problem.commit_scope if workflow.problem else ()
    return ResamplingDerivationRequest(
        references=(False,) * size,
        nucleus_groups=("G2N",) * size,
        observation_descriptors=tuple(f"ordinal={index}" for index in range(size)),
        scheme=ResamplingScheme.MONTE_CARLO,
        replicate_count=2,
        replicate_structural_identities=("replicate-alpha", "replicate-beta"),
        replicate_component_identities=(("component-alpha",), ("component-beta",)),
        root_seed=641,
        output_scope=scope,
        output_units=("native",) * len(scope),
        minimum_successful_count=1,
        strategy=OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "80"),),
        summary_policy=ResamplingSummaryPolicy(),
    )


def _publication_request(path: Path) -> MethodStepPublicationRequest:
    requested = Occurrence(
        "a" * 64,
        "b" * 64,
        "c" * 64,
        "unqualified-local-lane-v1",
        None,
        ("d" * 64,),
        "issue-641-publication-baseline",
    )
    bundle = ResultBundle.create(
        requested.identity,
        requested.execution_specification_identity,
        LegacyObservationImplementation(),
        (ResultMember("result", "e" * 64, 1),),
    )
    occurrence = requested.succeeded(bundle)
    return MethodStepPublicationRequest(
        path,
        ProvenanceEnvironment.from_current_process(),
        (
            BaselineReference.from_occurrence(occurrence),
            BaselineReference.from_result_bundle(bundle),
        ),
    )


def test_evaluation_only_retains_objective_without_fit_acceptance_or_commit() -> None:
    session, workflow = _evaluation_workflow()
    starting = session.analysis_values.snapshot()

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE, (
        outcome.publication_failure
    )
    assert isinstance(outcome.evaluation_result, EvaluationResult)
    assert outcome.evaluation_result.residuals.size > 0
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None
    assert session.analysis_values.snapshot() == starting
    assert require_successor_state(outcome, session.analysis_values) is (
        workflow.starting_snapshot
    )


def test_pre_cancelled_evaluation_never_creates_an_evaluator() -> None:
    session, workflow = _evaluation_workflow()
    starting = session.analysis_values.snapshot()
    cancellation = CancellationToken()
    cancellation.cancel()

    with patch.object(
        workflow.engine,
        "new_evaluator",
        side_effect=AssertionError("evaluate-only cancellation gate was crossed"),
    ) as evaluator_factory:
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
            cancellation=cancellation,
        )

    assert outcome.lifecycle is MethodStepLifecycle.CANCELLED
    assert outcome.primary_terminal == "cancelled"
    assert outcome.evaluation_result is None
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None
    assert session.analysis_values.snapshot() == starting
    evaluator_factory.assert_not_called()
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_no_objective_data_is_typed_and_does_not_evaluate_or_continue() -> None:
    session, workflow = _evaluation_workflow(
        purpose=EvaluationPurpose.NO_OBJECTIVE_DATA
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.NO_OBJECTIVE_DATA
    assert outcome.evaluation_result is None
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_no_objective_data_publishes_typed_diagnostics(tmp_path: Path) -> None:
    session, base = _evaluation_workflow(purpose=EvaluationPurpose.NO_OBJECTIVE_DATA)
    output = tmp_path / "no-objective"
    workflow = dataclasses.replace(base, publication=_publication_request(output))

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.NO_OBJECTIVE_DATA
    assert outcome.publication is not None, outcome.publication_failure
    assert (output / "Diagnostics" / "outcome.json").is_file()
    assert not (output / "Parameters").exists()


def test_explicit_evaluate_only_with_fit_roles_can_publish(tmp_path: Path) -> None:
    session, base = _evaluation_workflow(purpose=EvaluationPurpose.EVALUATE_ONLY)
    output = tmp_path / "evaluate-only"
    workflow = dataclasses.replace(base, publication=_publication_request(output))

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE, (
        outcome.publication_failure
    )
    assert outcome.publication is not None, outcome.publication_failure
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None
    assert require_successor_state(outcome, session.analysis_values).revision == 0


def test_serialized_outcome_cannot_recreate_successor_authority() -> None:
    session, workflow = _evaluation_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    historical = type(outcome).from_record(outcome.to_record())

    assert historical.to_record() == outcome.to_record()
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(historical, session.analysis_values)


def test_workflow_rejects_a_starting_snapshot_from_another_occurrence() -> None:
    _session, workflow = _evaluation_workflow()
    foreign_session, _foreign_workflow = _evaluation_workflow()

    with pytest.raises(ValueError, match="starting state"):
        MethodStepWorkflow.for_evaluation(
            starting_snapshot=foreign_session.analysis_values.snapshot(),
            parameter_model=workflow.parameter_model,
            parameterization=workflow.parameterization,
            engine=workflow.engine,
            method=read_methods([METHOD])["DEFAULT"],
            purpose=EvaluationPurpose.NO_OPTIMIZATION_REQUIRED,
        )


def test_one_component_direct_trf_uses_decomposition_aggregate_and_one_commit() -> None:
    session, workflow = _direct_workflow()
    assert workflow.decomposition is not None
    assert len(workflow.decomposition.components) == 1

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert isinstance(outcome.primary_execution, GroupedDirectTrfOutcome)
    assert len(outcome.primary_execution.components) == 1
    assert outcome.accepted_result is outcome.primary_execution.accepted_result
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.receipt is not None
    assert outcome.commit_operation.receipt.new_revision == 1
    successor = require_successor_state(outcome, session.analysis_values)
    assert successor.revision == 1


def test_mutated_different_budget_invocation_cannot_execute() -> None:
    session, workflow = _direct_workflow()
    decomposition = workflow.decomposition
    assert decomposition is not None
    replacement = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(81,) * len(decomposition.components),
    )
    original_identity = workflow.semantic_identity
    object.__setattr__(workflow, "invocation", replacement)

    with pytest.raises(ValueError, match="integrity"):
        execute_method_step(workflow, analysis_values=session.analysis_values)

    assert workflow.semantic_identity == original_identity
    assert session.analysis_values.snapshot().revision == 0


def test_recomputed_workflow_hash_cannot_repair_stale_invocation_child() -> None:
    session, workflow = _direct_workflow()
    invocation = workflow.invocation
    assert isinstance(invocation, GroupedDirectTrfInvocation)
    stale_child = invocation.component_invocations[0]
    object.__setattr__(
        stale_child,
        "objective_request_budget",
        stale_child.objective_request_budget + 1,
    )
    repaired_child = dataclasses.replace(stale_child)
    repaired_invocation = dataclasses.replace(
        invocation,
        component_invocations=(
            repaired_child,
            *invocation.component_invocations[1:],
        ),
    )
    repaired_workflow = dataclasses.replace(workflow, invocation=repaired_invocation)
    object.__setattr__(
        workflow,
        "semantic_identity",
        repaired_workflow.semantic_identity,
    )
    object.__setattr__(
        workflow,
        "binding_identity",
        repaired_workflow.binding_identity,
    )

    with pytest.raises(ValueError, match="child integrity"):
        execute_method_step(workflow, analysis_values=session.analysis_values)
    assert session.analysis_values.snapshot().revision == 0


def test_dataclass_replacement_rederives_invocation_identity() -> None:
    session, workflow = _direct_workflow()
    decomposition = workflow.decomposition
    assert decomposition is not None
    replacement = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(81,) * len(decomposition.components),
    )

    replaced = dataclasses.replace(workflow, invocation=replacement)

    assert replaced.semantic_identity != workflow.semantic_identity
    assert replaced.binding_identity != workflow.binding_identity
    outcome = execute_method_step(replaced, analysis_values=session.analysis_values)
    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED


def test_stale_caller_supplied_semantic_identity_cannot_execute() -> None:
    session, original = _direct_workflow()
    decomposition = original.decomposition
    assert decomposition is not None
    replacement = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(81,) * len(decomposition.components),
    )
    reconstructed = dataclasses.replace(original, invocation=replacement)
    object.__setattr__(
        reconstructed,
        "semantic_identity",
        original.semantic_identity,
    )

    with pytest.raises(ValueError, match="workflow integrity"):
        execute_method_step(reconstructed, analysis_values=session.analysis_values)
    assert session.analysis_values.snapshot().revision == 0


def test_mutated_committed_lifecycle_cannot_grant_successor() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    object.__setattr__(outcome, "lifecycle", MethodStepLifecycle.FAILED)

    with pytest.raises(ValueError, match="integrity"):
        require_successor_state(outcome, session.analysis_values)

    with pytest.raises(ValueError, match="integrity"):
        outcome.to_record()


def _recompute_outcome_record_identity(record: dict[str, object]) -> str:
    payload = dict(record)
    payload.pop("identity", None)
    return hashlib.sha256(
        json.dumps(
            ("native-method-step-outcome-v1", payload),
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
    ).hexdigest()


def _repair_live_outcome_record_identity(outcome: MethodStepOutcome) -> None:
    payload = outcome._record_payload()
    object.__setattr__(
        outcome,
        "record_identity",
        hashlib.sha256(
            json.dumps(
                ("native-method-step-outcome-v1", payload),
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            ).encode("ascii")
        ).hexdigest(),
    )


def test_serialized_lifecycle_tampering_survives_no_rehashed_envelope() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    record = outcome.to_record()
    record["lifecycle"] = MethodStepLifecycle.FAILED.value
    record["identity"] = _recompute_outcome_record_identity(record)

    with pytest.raises(ValueError, match="lifecycle integrity"):
        type(outcome).from_record(record)


def test_mutated_commit_operation_cannot_grant_successor() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    operation = outcome.commit_operation
    assert operation is not None
    tampered = dataclasses.replace(
        operation,
        occurrence_identity="tampered-commit-occurrence",
        _occurrence_witness=None,
    )
    object.__setattr__(outcome, "commit_operation", tampered)
    object.__setattr__(outcome, "commit_operation_identity", tampered.identity)
    _repair_live_outcome_record_identity(outcome)

    with pytest.raises(ValueError, match="integrity"):
        require_successor_state(outcome, session.analysis_values)


def test_mutated_successor_revision_cannot_grant_successor() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    operation = outcome.commit_operation
    assert operation is not None
    successor = operation.committed_snapshot
    assert successor is not None
    object.__setattr__(successor, "revision", successor.revision + 1)

    with pytest.raises(ValueError, match="integrity"):
        require_successor_state(outcome, session.analysis_values)


def test_mutated_source_workflow_cannot_grant_successor() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    _foreign_session, foreign = _direct_workflow()
    object.__setattr__(outcome, "source_workflow", foreign)

    with pytest.raises(ValueError, match="integrity"):
        require_successor_state(outcome, session.analysis_values)


def test_mutated_outcome_starting_state_cannot_grant_successor() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    object.__setattr__(outcome, "starting_state_identity", "stale-starting-state")
    _repair_live_outcome_record_identity(outcome)

    with pytest.raises(ValueError, match="integrity"):
        require_successor_state(outcome, session.analysis_values)


def test_operational_worker_and_thread_settings_are_not_semantic(
    tmp_path: Path,
) -> None:
    _session, base = _evaluation_workflow(purpose=EvaluationPurpose.EVALUATE_ONLY)
    mcmc_one = McmcSettings(steps=8, burn=2, walkers=6, seed=641, workers=1)
    mcmc_four = mcmc_one.model_copy(update={"workers": 4})
    one = dataclasses.replace(
        base,
        method=base.method.model_copy(update={"statistics": Statistics(mcmc=mcmc_one)}),
    )
    four = dataclasses.replace(
        base,
        method=base.method.model_copy(
            update={"statistics": Statistics(mcmc=mcmc_four)}
        ),
    )

    assert one.semantic_identity == four.semantic_identity
    assert one.binding_identity != four.binding_identity
    request = _publication_request(tmp_path / "unused")
    common = {
        "parameterization": one.parameterization,
        "plan": one.engine.plan,
        "semantic_workflow_identity": one.semantic_identity,
        "grouping_topology": (("aggregate", one.parameterization.independent_ids),),
        "policies": (),
        "budgets": (),
        "seeds": (),
        "execution": ExecutionSettings(),
        "environment": request.environment,
        "baseline_references": request.baseline_references,
    }
    provenance_one = WorkflowProvenance.create_method_step(method=one.method, **common)
    provenance_four = WorkflowProvenance.create_method_step(
        method=four.method, **common
    )
    assert provenance_one.execution_identity != provenance_four.execution_identity
    assert provenance_one.identity != provenance_four.identity

    _direct_session, direct = _direct_workflow()
    invocation = direct.invocation
    assert isinstance(invocation, GroupedDirectTrfInvocation)
    with_one_thread = GroupedDirectTrfInvocation(
        invocation.root_problem_identity,
        invocation.decomposition_identity,
        tuple(
            dataclasses.replace(
                item,
                execution_settings=ExecutionSettings(native_threads=1),
            )
            for item in invocation.component_invocations
        ),
    )
    with_four_threads = GroupedDirectTrfInvocation(
        invocation.root_problem_identity,
        invocation.decomposition_identity,
        tuple(
            dataclasses.replace(
                item,
                execution_settings=ExecutionSettings(native_threads=4),
            )
            for item in invocation.component_invocations
        ),
    )
    thread_one = dataclasses.replace(direct, invocation=with_one_thread)
    thread_four = dataclasses.replace(direct, invocation=with_four_threads)
    assert thread_one.semantic_identity == thread_four.semantic_identity
    assert thread_one.binding_identity != thread_four.binding_identity
    decomposition = direct.decomposition
    assert decomposition is not None
    thread_grouping = tuple(
        (component.identity, component.controlled_ids)
        for component in decomposition.components
    )
    thread_common = {
        "parameterization": direct.parameterization,
        "plan": direct.engine.plan,
        "method": direct.method,
        "semantic_workflow_identity": thread_one.semantic_identity,
        "grouping_topology": thread_grouping,
        "policies": (),
        "budgets": (),
        "seeds": (),
        "environment": request.environment,
        "baseline_references": request.baseline_references,
    }
    thread_provenance_one = WorkflowProvenance.create_method_step(
        execution=ExecutionSettings(native_threads=1),
        **thread_common,
    )
    thread_provenance_four = WorkflowProvenance.create_method_step(
        execution=ExecutionSettings(native_threads=4),
        **thread_common,
    )
    assert (
        thread_provenance_one.execution_identity
        != thread_provenance_four.execution_identity
    )


def test_scientific_mcmc_setting_remains_semantic() -> None:
    _session, base = _evaluation_workflow(purpose=EvaluationPurpose.EVALUATE_ONLY)
    first = dataclasses.replace(
        base,
        method=base.method.model_copy(
            update={
                "statistics": Statistics(
                    mcmc=McmcSettings(
                        steps=8,
                        burn=2,
                        walkers=6,
                        seed=641,
                        workers=1,
                    )
                )
            }
        ),
    )
    changed = dataclasses.replace(
        first,
        method=first.method.model_copy(
            update={
                "statistics": Statistics(
                    mcmc=McmcSettings(
                        steps=10,
                        burn=2,
                        walkers=6,
                        seed=641,
                        workers=1,
                    )
                )
            }
        ),
    )

    assert first.semantic_identity != changed.semantic_identity


def test_optimization_workflow_rejects_foreign_or_zero_component_decomposition() -> (
    None
):
    _session, workflow = _direct_workflow()
    _foreign_session, foreign = _direct_workflow()
    assert workflow.problem is not None
    assert workflow.decomposition is not None
    assert workflow.invocation is not None
    assert foreign.decomposition is not None

    with pytest.raises(DirectTrfConstructionError, match="decomposition"):
        dataclasses.replace(workflow.decomposition, components=())

    with pytest.raises(ValueError, match="decomposition"):
        MethodStepWorkflow.for_optimization(
            starting_snapshot=workflow.starting_snapshot,
            parameter_model=workflow.parameter_model,
            parameterization=workflow.parameterization,
            engine=workflow.engine,
            method=workflow.method,
            problem=workflow.problem,
            decomposition=foreign.decomposition,
            strategy=MethodStepStrategy.DIRECT_TRF,
            invocation=workflow.invocation,
        )


def test_multi_component_direct_trf_keeps_components_non_authoritative() -> None:
    session, workflow = _direct_workflow(all_profiles=True)
    assert workflow.decomposition is not None
    assert len(workflow.decomposition.components) > 1

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert isinstance(outcome.primary_execution, GroupedDirectTrfOutcome)
    assert len(outcome.primary_execution.components) == len(
        workflow.decomposition.components
    )
    assert all(
        component.candidate is not None
        for component in outcome.primary_execution.components
    )
    assert outcome.primary_execution.accepted_result is outcome.accepted_result


@pytest.mark.parametrize("all_profiles", [False, True])
def test_grid_then_trf_always_uses_grouped_decomposition(all_profiles: bool) -> None:
    session, workflow = _grid_workflow(all_profiles=all_profiles)

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert isinstance(outcome.primary_execution, GroupedGridDirectTrfOutcome)
    assert outcome.accepted_result is outcome.primary_execution.accepted_result
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.receipt is not None


def test_de_then_trf_routes_existing_seeded_search_and_commits_once() -> None:
    session, workflow = _de_workflow()

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert isinstance(outcome.primary_execution, DeDirectTrfOutcome)
    assert outcome.accepted_result is outcome.primary_execution.accepted_result
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.receipt is not None


def test_cancellation_before_primary_starts_blocks_commit_and_successor() -> None:
    session, workflow = _direct_workflow()
    cancellation = CancellationToken()
    cancellation.cancel()

    outcome = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        cancellation=cancellation,
    )

    assert outcome.lifecycle is MethodStepLifecycle.CANCELLED
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_stale_commit_retains_accepted_diagnostics_and_blocks_successor() -> None:
    session, workflow = _direct_workflow()

    def make_stale(checkpoint: MethodStepCheckpoint) -> None:
        if checkpoint is not MethodStepCheckpoint.AGGREGATE_ACCEPTED:
            return
        current = session.analysis_values.snapshot()
        session.analysis_values.commit(
            dict(current.items()),
            expected=current,
            scope=tuple(current),
        )

    outcome = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        checkpoint_observer=make_stale,
    )

    assert outcome.lifecycle is MethodStepLifecycle.ACCEPTED_UNCOMMITTED
    assert outcome.accepted_result is not None
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.failure is not None
    assert outcome.commit_operation.failure.category.value == "stale_revision"
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_stale_first_step_cannot_compile_a_required_second_step() -> None:
    session, workflow = _direct_workflow()

    def make_stale(checkpoint: MethodStepCheckpoint) -> None:
        if checkpoint is MethodStepCheckpoint.AGGREGATE_ACCEPTED:
            current = session.analysis_values.snapshot()
            session.analysis_values.commit(
                dict(current.items()),
                expected=current,
                scope=tuple(current),
            )

    first = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        checkpoint_observer=make_stale,
    )
    second_compilation_count = 0

    def compile_second() -> None:
        nonlocal second_compilation_count
        require_successor_state(first, session.analysis_values)
        second_compilation_count += 1

    with pytest.raises(ValueError, match="successor authority"):
        compile_second()
    assert second_compilation_count == 0


def test_optimizer_non_convergence_is_typed_and_never_committed() -> None:
    session, workflow = _direct_workflow()

    def unsuccessful(*_args: object, **_kwargs: object) -> object:
        problem = workflow.problem
        assert problem is not None
        return SimpleNamespace(
            x=np.asarray(problem.start),
            fun=np.ones(workflow.engine.plan.retained_observation_count),
            status=0,
            success=False,
            message="qualified non-convergence",
            nfev=1,
            njev=0,
            cost=1.0,
            optimality=1.0,
            active_mask=np.zeros(len(problem.start), dtype=np.int64),
        )

    with patch("chemex.optimize.direct_trf.least_squares", unsuccessful):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert outcome.lifecycle is MethodStepLifecycle.FAILED
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None


def test_requested_uncertainty_runs_only_after_commit_and_cannot_mutate_central_values() -> (
    None
):
    session, base = _direct_workflow(stabilize_data=True)
    request = _uncertainty_request(base)
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        derivations=(request,),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert len(outcome.derivations) == 1
    assert outcome.derivations[0].stage == "uncertainty"
    assert outcome.derivations[0].disposition in {
        DerivationDisposition.COMPLETED,
        DerivationDisposition.FAILED,
    }
    assert outcome.commit_operation is not None
    assert (
        outcome.commit_operation.committed_snapshot
        == session.analysis_values.snapshot()
    )


def test_cancellation_after_commit_stops_requested_derivations_without_rollback() -> (
    None
):
    session, base = _direct_workflow()
    request = _uncertainty_request(base)
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        derivations=(request,),
    )
    cancellation = CancellationToken()

    def cancel_after_commit(checkpoint: MethodStepCheckpoint) -> None:
        if checkpoint is MethodStepCheckpoint.COMMIT_COMPLETED:
            cancellation.cancel()

    outcome = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        cancellation=cancellation,
        checkpoint_observer=cancel_after_commit,
    )

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.receipt is not None
    assert outcome.derivations[0].disposition is (
        DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP
    )
    assert require_successor_state(outcome, session.analysis_values).revision == 1


def test_cancellation_after_acceptance_blocks_commit_and_all_derivations() -> None:
    session, base = _direct_workflow()
    request = _uncertainty_request(base)
    workflow = dataclasses.replace(base, derivations=(request,))
    cancellation = CancellationToken()

    def cancel_at_acceptance(checkpoint: MethodStepCheckpoint) -> None:
        if checkpoint is MethodStepCheckpoint.AGGREGATE_ACCEPTED:
            cancellation.cancel()

    outcome = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        cancellation=cancellation,
        checkpoint_observer=cancel_at_acceptance,
    )

    assert outcome.lifecycle is MethodStepLifecycle.CANCELLED
    assert outcome.accepted_result is not None
    assert outcome.commit_operation is None
    assert session.analysis_values.snapshot().revision == 0
    assert outcome.derivations[0].disposition is (
        DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP
    )
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_independent_resampling_runs_after_uncertainty_branch_failure() -> None:
    session, base = _direct_workflow(stabilize_data=True)
    valid_uncertainty = _uncertainty_request(base)
    incompatible_policy = dataclasses.replace(
        valid_uncertainty.policy,
        coordinate_scales=(("__FOREIGN", 1.0),),
        coordinate_units=(("__FOREIGN", ParameterUnit.RATE_PER_SECOND),),
    )
    uncertainty = dataclasses.replace(
        valid_uncertainty,
        policy=incompatible_policy,
    )
    resampling = _resampling_request(base)
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        derivations=(uncertainty, resampling),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert [item.stage for item in outcome.derivations] == [
        "uncertainty",
        "resampling",
    ]
    assert outcome.derivations[0].disposition is DerivationDisposition.FAILED
    assert outcome.derivations[1].disposition is DerivationDisposition.COMPLETED
    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED


def test_interrupted_derivation_stops_later_requested_stages() -> None:
    session, base = _direct_workflow(stabilize_data=True)
    resampling = _resampling_request(base)
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    policy = ExpertMcmcPolicy(
        burn_steps=1,
        retained_steps=2,
        walkers=max(6, 2 * len(base.problem.controlled_ids)),
        expert_provenance="issue-641-interruption",
    ).resolve(dimension=len(base.problem.controlled_ids), root_seed=641)
    mcmc = McmcDerivationRequest(
        policy,
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS)
            for param_id in base.problem.controlled_ids
        ),
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS)
            for param_id in base.problem.commit_scope
        ),
    )
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        derivations=(resampling, mcmc),
    )
    interrupted = SimpleNamespace(
        terminal=ResamplingOperationTerminal.INTERRUPTED,
        identity="qualified-resampling-interruption",
        evidence=None,
    )

    with (
        patch(
            "chemex.optimize.method_step.execute_resampling_evidence",
            return_value=interrupted,
        ),
        patch("chemex.optimize.method_step.execute_mcmc_evidence") as execute_mcmc,
    ):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert [item.disposition for item in outcome.derivations] == [
        DerivationDisposition.INTERRUPTED,
        DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP,
    ]
    execute_mcmc.assert_not_called()


def test_committed_publication_retains_partial_derivation_evidence(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow(stabilize_data=True)
    request = _resampling_request(base)
    output = tmp_path / "partial"
    workflow = dataclasses.replace(
        base,
        derivations=(request,),
        publication=_publication_request(output),
    )

    def cancel_resampling(
        accepted: AcceptedFitResult,
        plan: ResamplingPlan,
        **_kwargs: object,
    ) -> ResamplingOperation:
        return execute_resampling_evidence(
            accepted,
            plan,
            cancellation_probe=lambda: True,
        )

    with patch(
        "chemex.optimize.method_step.execute_resampling_evidence",
        side_effect=cancel_resampling,
    ):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED, (
        outcome.publication_failure
    )
    assert outcome.derivations[0].disposition is DerivationDisposition.CANCELLED
    assert outcome.publication is not None
    assert (
        output / "PartialEvidence" / "Resampling" / "MC" / "evidence.json"
    ).is_file()
    assert (output / "Parameters" / "restart.toml").is_file()


def test_resampling_summary_failure_is_preserved_and_published(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow(stabilize_data=True)
    request = _resampling_request(base)
    output = tmp_path / "summary-failure"
    workflow = dataclasses.replace(
        base,
        derivations=(request,),
        publication=_publication_request(output),
    )
    captured_failure: SummaryFailure | None = None

    def fail_summary(
        evidence: ResamplingEvidence,
        _policy: ResamplingSummaryPolicy,
    ) -> ResamplingSummaryOutcome:
        nonlocal captured_failure
        captured_failure = SummaryFailure(
            evidence.identity,
            "qualification_summary_failure",
            "Qualification forced the derived summary to fail",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.INSUFFICIENT_COVERAGE,
            failure=captured_failure,
        )

    with patch(
        "chemex.optimize.method_step.summarize_resampling_evidence",
        side_effect=fail_summary,
    ) as summarize:
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert captured_failure is not None
    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED, (
        outcome.publication_failure
    )
    assert outcome.derivations[0].disposition is DerivationDisposition.FAILED
    assert captured_failure.identity in outcome.derivations[0].artifact_identities
    assert summarize.call_count == 1
    assert (
        output / "PartialEvidence" / "Resampling" / "MC" / "summary-failure.json"
    ).is_file()


@pytest.mark.parametrize(
    ("posterior_error", "expected_disposition"),
    [
        (
            ValueError("qualification posterior failure"),
            DerivationDisposition.FAILED,
        ),
        (KeyboardInterrupt(), DerivationDisposition.INTERRUPTED),
    ],
)
def test_mcmc_posterior_stop_retains_completed_primary_evidence(
    tmp_path: Path,
    posterior_error: BaseException,
    expected_disposition: DerivationDisposition,
) -> None:
    session, inputs = _optimization_inputs(
        method=Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]),
        stabilize_data=True,
    )
    problem = inputs["problem"]
    decomposition = inputs["decomposition"]
    assert isinstance(problem, OptimizationProblem)
    assert isinstance(decomposition, FitDecomposition)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    policy = ExpertMcmcPolicy(
        burn_steps=1,
        retained_steps=2,
        walkers=max(6, 2 * len(problem.controlled_ids)),
        expert_provenance="issue-641-posterior-failure",
    ).resolve(dimension=len(problem.controlled_ids), root_seed=641)
    request = McmcDerivationRequest(
        policy,
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS)
            for param_id in problem.controlled_ids
        ),
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS) for param_id in problem.commit_scope
        ),
    )
    output = tmp_path / "mcmc-posterior-failure"
    workflow = MethodStepWorkflow.for_optimization(
        **inputs,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=invocation,
        derivations=(request,),
        publication=_publication_request(output),
    )

    with patch(
        "chemex.optimize.method_step.derive_retained_sample_view",
        side_effect=posterior_error,
    ):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    derivation = outcome.derivations[0]
    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED, (
        outcome.publication_failure
    )
    assert derivation.disposition is expected_disposition
    assert isinstance(derivation.artifacts[0], McmcEvidence)
    assert type(posterior_error).__name__ in derivation.message
    assert (output / "PartialEvidence" / "MCMC" / "evidence.json").is_file()


def test_requested_mcmc_retains_exact_committed_fit_lineage() -> None:
    session, inputs = _optimization_inputs(
        method=Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]),
        stabilize_data=True,
    )
    problem = inputs["problem"]
    decomposition = inputs["decomposition"]
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    policy = ExpertMcmcPolicy(
        burn_steps=1,
        retained_steps=2,
        walkers=6,
        expert_provenance="issue-641-explicit-mcmc",
    ).resolve(dimension=len(problem.controlled_ids), root_seed=641)
    mcmc = McmcDerivationRequest(
        policy,
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS)
            for param_id in problem.controlled_ids
        ),
        tuple(
            (param_id, ParameterUnit.DIMENSIONLESS) for param_id in problem.commit_scope
        ),
    )
    workflow = MethodStepWorkflow.for_optimization(
        **inputs,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=invocation,
        derivations=(mcmc,),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert outcome.derivations[0].stage == "mcmc"
    assert outcome.derivations[0].disposition is DerivationDisposition.COMPLETED
    assert len(outcome.derivations[0].artifact_identities) == 3


def test_successful_method_step_publication_uses_atomic_native_step_root(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow()
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    output = tmp_path / "step-0001"
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        publication=_publication_request(output),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert outcome.publication is not None, outcome.publication_failure
    outcome.publication.require_exact_live_publication()
    assert outcome.publication.path == output
    assert (output / "fit-manifest.toml").is_file()
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert (output / "Parameters" / "restart.toml").is_file()


def test_publication_failure_never_rolls_back_committed_fit(tmp_path: Path) -> None:
    session, base = _direct_workflow()
    assert base.problem is not None
    assert base.decomposition is not None
    assert base.invocation is not None
    output = tmp_path / "occupied"
    output.mkdir()
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=base.starting_snapshot,
        parameter_model=base.parameter_model,
        parameterization=base.parameterization,
        engine=base.engine,
        method=base.method,
        problem=base.problem,
        decomposition=base.decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=base.invocation,
        publication=_publication_request(output),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.PUBLICATION_FAILED
    assert outcome.commit_operation is not None
    assert outcome.commit_operation.receipt is not None
    assert session.analysis_values.snapshot().revision == 1
    assert "FileExistsError" in outcome.publication_failure
    with pytest.raises(ValueError, match="successor authority"):
        require_successor_state(outcome, session.analysis_values)


def test_evaluation_only_publication_uses_existing_native_reporting(
    tmp_path: Path,
) -> None:
    session, base = _evaluation_workflow()
    output = tmp_path / "evaluation"
    workflow = dataclasses.replace(base, publication=_publication_request(output))

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE, (
        outcome.publication_failure
    )
    assert outcome.publication is not None, outcome.publication_failure
    assert (output / "fit-manifest.toml").is_file()
    assert not (output / "Parameters" / "restart.toml").exists()
    assert require_successor_state(outcome, session.analysis_values).revision == 0


def test_failed_primary_publishes_diagnostics_only(tmp_path: Path) -> None:
    session, base = _direct_workflow()
    output = tmp_path / "failed"
    workflow = dataclasses.replace(base, publication=_publication_request(output))

    with patch.object(workflow.engine, "project_profiles", side_effect=RuntimeError):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert outcome.lifecycle is MethodStepLifecycle.FAILED
    assert outcome.publication is not None, outcome.publication_failure
    assert (output / "Diagnostics" / "outcome.json").is_file()
    assert not (output / "Parameters").exists()
    assert not (output / "Statistics").exists()


def test_stale_commit_publication_suppresses_fitted_and_restart_output(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow()
    output = tmp_path / "stale"
    workflow = dataclasses.replace(base, publication=_publication_request(output))

    def make_stale(checkpoint: MethodStepCheckpoint) -> None:
        if checkpoint is MethodStepCheckpoint.AGGREGATE_ACCEPTED:
            current = session.analysis_values.snapshot()
            session.analysis_values.commit(
                dict(current.items()),
                expected=current,
                scope=tuple(current),
            )

    outcome = execute_method_step(
        workflow,
        analysis_values=session.analysis_values,
        checkpoint_observer=make_stale,
    )

    assert outcome.lifecycle is MethodStepLifecycle.ACCEPTED_UNCOMMITTED
    assert outcome.publication is not None, outcome.publication_failure
    assert (output / "Diagnostics" / "outcome.json").is_file()
    assert not (output / "Parameters").exists()


def test_interruption_during_primary_starts_no_later_stage() -> None:
    session, workflow = _direct_workflow()

    with patch.object(
        workflow.engine, "project_profiles", side_effect=KeyboardInterrupt
    ):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert outcome.lifecycle is MethodStepLifecycle.INTERRUPTED
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None


def test_aggregate_materialization_failure_never_accepts_or_commits() -> None:
    session, workflow = _direct_workflow()
    calls = 0

    def fail_aggregate(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        return RootMaterializationFailure(
            MaterializationTerminal.FAILURE,
            1,
            workflow.engine.compatibility_identity,
            0,
            1,
            TerminalFailure("qualified_aggregate_failure"),
        )

    with patch(
        "chemex.optimize.grouped_direct_trf.direct_trf_owner.materialize_root_candidate",
        fail_aggregate,
    ):
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert calls == 1
    assert outcome.lifecycle is MethodStepLifecycle.FAILED
    assert outcome.accepted_result is None
    assert outcome.commit_operation is None


def test_two_steps_compile_second_from_exact_committed_successor() -> None:
    session, first = _direct_workflow(stabilize_data=True)
    first_outcome = execute_method_step(
        first,
        analysis_values=session.analysis_values,
    )
    successor = require_successor_state(first_outcome, session.analysis_values)
    parameterization = session.compile_parameterization(
        first.method,
        set(first.parameterization.scope_ids),
    )
    engine = first.engine.rebind_parameterization(parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None
    assert parameter_model is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        successor,
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    second = MethodStepWorkflow.for_optimization(
        starting_snapshot=successor,
        parameter_model=parameter_model,
        parameterization=parameterization,
        engine=engine,
        method=first.method,
        problem=problem,
        decomposition=decomposition,
        strategy=MethodStepStrategy.DIRECT_TRF,
        invocation=GroupedDirectTrfInvocation.for_decomposition(
            decomposition,
            objective_request_budgets=(80,) * len(decomposition.components),
        ),
    )

    second_outcome = execute_method_step(
        second,
        analysis_values=session.analysis_values,
    )

    assert second.starting_snapshot is successor
    assert second_outcome.lifecycle is MethodStepLifecycle.COMMITTED
    assert (
        require_successor_state(second_outcome, session.analysis_values).revision == 2
    )


def test_method_step_publication_can_write_existing_run_information(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow()
    output = tmp_path / "Output"
    captured = capture_native_inputs(
        Namespace(
            experiments=[EXPERIMENT],
            parameters=[PARAMETERS],
            method=[METHOD],
            output=output,
        ),
        working_directory=ROOT,
    )
    run_request = MethodStepRunProvenanceRequest(
        output,
        "issue-641-native-method-step",
        captured,
        base.starting_snapshot,
        working_directory=ROOT,
    )
    publication = dataclasses.replace(
        _publication_request(output / "STEP-0001"),
        run_provenance=run_request,
    )
    workflow = dataclasses.replace(base, publication=publication)

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED, (
        outcome.publication_failure
    )
    assert outcome.run_provenance is not None
    assert outcome.publication is not None
    outcome.publication.require_exact_live_publication()
    assert outcome.run_provenance_identity == (
        outcome.run_provenance.execution_occurrence_identity
    )
    assert (output / "run_info" / "run.toml").is_file()


def test_run_information_failure_retains_successful_step_publication(
    tmp_path: Path,
) -> None:
    session, base = _direct_workflow()
    output = tmp_path / "Output"
    (output / "run_info").mkdir(parents=True)
    captured = capture_native_inputs(
        Namespace(
            experiments=[EXPERIMENT],
            parameters=[PARAMETERS],
            method=[METHOD],
            output=output,
        ),
        working_directory=ROOT,
    )
    run_request = MethodStepRunProvenanceRequest(
        output,
        "issue-641-run-info-failure",
        captured,
        base.starting_snapshot,
        working_directory=ROOT,
    )
    workflow = dataclasses.replace(
        base,
        publication=dataclasses.replace(
            _publication_request(output / "STEP-0001"),
            run_provenance=run_request,
        ),
    )

    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)

    assert outcome.lifecycle is MethodStepLifecycle.PUBLICATION_FAILED
    assert outcome.publication is not None
    assert outcome.publication_identity == outcome.publication.identity
    outcome.publication.require_exact_live_publication()
    assert (output / "STEP-0001" / "fit-manifest.toml").is_file()
    assert "FileExistsError" in outcome.publication_failure


def test_semantic_identity_excludes_publication_occurrence_facts(
    tmp_path: Path,
) -> None:
    _session, base = _direct_workflow()
    first = dataclasses.replace(
        base,
        publication=_publication_request(tmp_path / "one"),
    )
    second = dataclasses.replace(
        base,
        publication=_publication_request(tmp_path / "two"),
    )

    assert first.semantic_identity == second.semantic_identity
    assert first.binding_identity != second.binding_identity


def test_semantic_identity_includes_strategy_and_derivation_choices() -> None:
    _session, direct = _direct_workflow()
    problem = direct.problem
    assert problem is not None
    grid = dataclasses.replace(
        direct,
        strategy=MethodStepStrategy.GRID_DIRECT_TRF,
        invocation=GridDirectTrfInvocation.for_problem(
            problem,
            axes=((problem.controlled_ids[0], (problem.start[0],)),),
            objective_request_budget=80,
        ),
    )
    with_uncertainty = dataclasses.replace(
        direct,
        derivations=(_uncertainty_request(direct),),
    )

    assert direct.semantic_identity != grid.semantic_identity
    assert direct.semantic_identity != with_uncertainty.semantic_identity


def test_historical_outcome_rejects_lifecycle_tampering() -> None:
    session, workflow = _direct_workflow()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    record = outcome.to_record()
    record["lifecycle"] = MethodStepLifecycle.FAILED.value

    with pytest.raises(ValueError, match="integrity"):
        type(outcome).from_record(record)


def test_native_method_step_never_calls_legacy_run_methods() -> None:
    session, workflow = _direct_workflow()

    with patch(
        "chemex.optimize.fitting.run_methods",
        side_effect=AssertionError("legacy execution was called"),
    ) as legacy:
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )

    assert outcome.lifecycle is MethodStepLifecycle.COMMITTED
    legacy.assert_not_called()
