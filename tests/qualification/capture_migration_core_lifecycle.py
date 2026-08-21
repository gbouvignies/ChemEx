"""Execute #592 probes and return observation-only production facts."""

# ruff: noqa: I001

from __future__ import annotations

# fmt: off
import dataclasses
import hashlib
import json
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

import numpy as np

import chemex.optimize.direct_trf as direct
import chemex.optimize.method_step as step
import chemex.optimize.native_resampling as resampling
from chemex.baselines import (
    CanonicalBaselineValue, CaseDefinition, CaseSourceAuthority,
    ExecutionSpecification, InputMember, LegacyObservationImplementation,
    Occurrence, ResultBundle, ResultMember,
)
from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.migration_core_lifecycle import LifecycleProbeCapture, LifecycleProbeResult
from chemex.native_provenance import BaselineReference, ProvenanceEnvironment
from chemex.numerical_lanes import LaneAttestation, LiveLaneAuthority
from chemex.optimize.grouped_direct_trf import FitDecomposition, GroupedDirectTrfInvocation
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.runtime import AnalysisSession
ROOT = Path(__file__).parents[2]
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
FROZEN_LEGACY_IMPLEMENTATION = LegacyObservationImplementation(
    package_version="2026.6.1",
    lmfit_version="1.3.4",
    source_manifest_hash="46977ae3475257d7c2561d7f53ddf5b36cc6b221a7aa9f43e73403140278f846",
)
@dataclass(frozen=True, slots=True)
class ProbeObservation:
    """Raw objects and entry facts emitted by one production-seam execution."""

    probe_id: str
    starting_snapshot: AnalysisValuesSnapshot
    ending_snapshot: AnalysisValuesSnapshot
    case: CaseDefinition
    specification: ExecutionSpecification
    occurrence: Occurrence
    outcome: step.MethodStepOutcome | None = None
    construction_failure: direct.DirectTrfConstructionError | None = None
    constructor_entries: tuple[str, ...] = ()
    executor_entries: tuple[str, ...] = ()
    successor_snapshot: AnalysisValuesSnapshot | None = None
    successor_denial: ValueError | None = None
    checkpoint_entries: tuple[str, ...] = ()
    checkpoint_failures: tuple[str, ...] = ()
    downstream_entries: tuple[str, ...] = ()
    outcomes: tuple[step.MethodStepOutcome, ...] = ()
def _base_workflow(
    *, no_objective: bool = False, stabilize: bool = False
) -> tuple[AnalysisSession, step.MethodStepWorkflow]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT], Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None), session=session
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = Method(fix=["R1A_A", "PB", "KEX_AB"]) if no_objective else read_methods([METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    if no_objective:
        for experiment in experiments:
            for profile in experiment.profiles:
                profile.is_scaled = False
                profile.data.mask[:] = False
                profile.data.mark_dirty()
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    if stabilize:
        frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
        )
        evaluated = engine.new_evaluator().evaluate(frame)
        assert isinstance(evaluated, EvaluationResult)
        offset = 0
        for experiment in experiments:
            for profile in experiment.profiles:
                stop = offset + profile.data.size
                profile.data.exp = np.asarray(evaluated.normalized_calculations[offset:stop]).copy()
                profile.data.mark_dirty()
                offset = stop
        engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None and parameter_model is not None
    starting = session.analysis_values.snapshot()
    if no_objective:
        return session, step.MethodStepWorkflow.for_evaluation(
            starting_snapshot=starting, parameter_model=parameter_model,
            parameterization=parameterization, engine=engine, method=method,
            purpose=step.EvaluationPurpose.NO_OBJECTIVE_DATA,
        )
    problem = direct.OptimizationProblem.from_native(engine.plan, parameterization, configuration, starting)
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    return session, step.MethodStepWorkflow.for_optimization(
        starting_snapshot=starting, parameter_model=parameter_model,
        parameterization=parameterization, engine=engine, method=method,
        problem=problem, decomposition=decomposition,
        strategy=step.MethodStepStrategy.DIRECT_TRF,
        invocation=GroupedDirectTrfInvocation.for_decomposition(
            decomposition, objective_request_budgets=(80,) * len(decomposition.components)
        ),
    )
def _resampling_request(workflow: step.MethodStepWorkflow) -> step.ResamplingDerivationRequest:
    size = workflow.engine.plan.observation_count
    assert workflow.problem is not None
    scope = workflow.problem.commit_scope
    return step.ResamplingDerivationRequest(
        references=(False,) * size, nucleus_groups=("G2N",) * size,
        observation_descriptors=tuple(f"ordinal={index}" for index in range(size)),
        scheme=resampling.ResamplingScheme.MONTE_CARLO, replicate_count=2,
        replicate_structural_identities=("replicate-alpha", "replicate-beta"),
        replicate_component_identities=(("component-alpha",), ("component-beta",)),
        root_seed=592, output_scope=scope, output_units=("native",) * len(scope),
        minimum_successful_count=1, strategy=resampling.OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "80"),),
        summary_policy=resampling.ResamplingSummaryPolicy(),
    )
def _member(role: str, path: Path) -> InputMember:
    content = path.read_bytes()
    return InputMember(role, hashlib.sha256(content).hexdigest(), len(content))


def _publication(
    path: Path, source_commit: str, lockfile_hash: str
) -> step.MethodStepPublicationRequest:
    implementation = FROZEN_LEGACY_IMPLEMENTATION
    case = CaseDefinition.create(
        "migration-core-lifecycle-publication-reference",
        CaseSourceAuthority(source_commit, lockfile_hash),
        {"purpose": "publication-provenance-context"},
        (_member("experiment", EXPERIMENT),),
    )
    specification = ExecutionSpecification.create(
        case, implementation,
        workflow={"purpose": "publication-provenance-context"},
        lane_reference="unqualified-local-lane-v1", policy={}, budget={}, seed=None,
        execution_settings={}, artifact_inventory={"roles": ["result"]},
        roles=("qualification:publication-reference",), claims=("context-only",),
    )
    requested = Occurrence.requested(
        specification, case, f"publication-reference:{uuid4().hex}"
    )
    content = json.dumps(requested.to_record(), sort_keys=True).encode("ascii")
    bundle = ResultBundle.create(
        requested.identity, specification.identity, implementation,
        (ResultMember("result", hashlib.sha256(content).hexdigest(), len(content)),),
    )
    return step.MethodStepPublicationRequest(
        path, ProvenanceEnvironment.from_current_process(),
        (BaselineReference.from_occurrence(requested.succeeded(bundle)), BaselineReference.from_result_bundle(bundle)),
    )


def _lineage(
    probe_id: str, workflow: step.MethodStepWorkflow, authority: LiveLaneAuthority,
    source_commit: str, lockfile_hash: str, *, required_steps: tuple[str, ...] = (),
) -> tuple[CaseDefinition, ExecutionSpecification, Occurrence]:
    attestation = LaneAttestation.from_record(authority.to_record())
    implementation = FROZEN_LEGACY_IMPLEMENTATION
    case = CaseDefinition.create(
        f"migration-core-lifecycle:{probe_id}",
        CaseSourceAuthority(source_commit, lockfile_hash),
        {"probe_id": probe_id, "required_steps": list(required_steps)},
        tuple(
            _member(role, path)
            for role, path in (
                ("capture-runner", Path(__file__)), ("experiment", EXPERIMENT),
                ("method", METHOD), ("parameters", PARAMETERS),
            )
        ),
    )
    specification = ExecutionSpecification.create(
        case, implementation,
        workflow={
            "semantic_identity": workflow.semantic_identity,
            "binding_identity": workflow.binding_identity,
            "evaluation_plan_identity": workflow.engine.plan.identity,
            "parameterization_identity": workflow.parameterization.identity,
            "starting_occurrence_identity": workflow.starting_snapshot.occurrence_identity,
            "starting_revision": workflow.starting_snapshot.revision,
            "required_steps": list(required_steps),
        },
        lane_reference=attestation.lane_identity,
        policy={"probe_id": probe_id},
        budget={
            "components": 0 if workflow.decomposition is None else len(workflow.decomposition.components),
            "replicates": next(
                (item.replicate_count for item in workflow.derivations if isinstance(item, step.ResamplingDerivationRequest)), 0
            ),
        },
        seed=next(
            (item.root_seed for item in workflow.derivations if isinstance(item, step.ResamplingDerivationRequest)), None
        ),
        execution_settings={
            "environment_identity": attestation.environment_identity,
            "workers": attestation.workers,
            "native_threads": attestation.native_threads,
        },
        artifact_inventory={"owner_records": ["method_step", "primary", "commit", "resampling", "successor"]},
        roles=("qualification:migration-core-lifecycle",),
        claims=("typed-runtime-facts",),
    )
    occurrence = Occurrence.requested(
        specification, case, f"{probe_id}:{uuid4().hex}", authority
    )
    return case, specification, occurrence
def _observed(
    probe_id: str, session: AnalysisSession, starting: AnalysisValuesSnapshot,
    outcome: step.MethodStepOutcome, lineage: tuple[CaseDefinition, ExecutionSpecification, Occurrence],
    checkpoint_entries: tuple[str, ...] = (), checkpoint_failures: tuple[str, ...] = (),
) -> ProbeObservation:
    downstream: list[str] = []
    try:
        successor = step.require_successor_state(outcome, session.analysis_values)
    except ValueError as denial:
        return ProbeObservation(
            probe_id, starting, session.analysis_values.snapshot(), *lineage, outcome,
            successor_denial=denial, checkpoint_entries=checkpoint_entries,
            checkpoint_failures=checkpoint_failures,
        )
    downstream.append("dependent")
    return ProbeObservation(
        probe_id, starting, session.analysis_values.snapshot(), *lineage, outcome,
        successor_snapshot=successor, checkpoint_entries=checkpoint_entries,
        checkpoint_failures=checkpoint_failures,
        downstream_entries=tuple(downstream),
    )
def observe_two_required_steps(
    *, no_objective: bool, authority: LiveLaneAuthority,
    source_commit: str, lockfile_hash: str,
) -> ProbeObservation:
    """Run one unconditional two-step production composition path."""

    session, template = _base_workflow(
        no_objective=no_objective, stabilize=not no_objective
    )
    starting = session.analysis_values.snapshot()
    lineage = _lineage(
        "construction-dependent-stop", template, authority, source_commit,
        lockfile_hash, required_steps=("first", "second"),
    )
    constructors: list[str] = []
    executors: list[str] = []
    outcomes: list[step.MethodStepOutcome] = []

    def compose_and_execute(label: str, snapshot: AnalysisValuesSnapshot) -> step.MethodStepOutcome:
        constructors.append(label)
        configuration = session.parameter_factory.sealed_configuration
        parameter_model = session.parameter_factory.sealed_parameter_model
        assert configuration is not None and parameter_model is not None
        parameterization = session.compile_parameterization(template.method, set(template.parameterization.scope_ids))
        engine = template.engine.rebind_parameterization(parameterization)
        problem = direct.OptimizationProblem.from_native(engine.plan, parameterization, configuration, snapshot)
        decomposition = FitDecomposition.from_root(problem, parameterization, engine)
        workflow = step.MethodStepWorkflow.for_optimization(
            starting_snapshot=snapshot, parameter_model=parameter_model,
            parameterization=parameterization, engine=engine, method=template.method,
            problem=problem, decomposition=decomposition,
            strategy=step.MethodStepStrategy.DIRECT_TRF,
            invocation=GroupedDirectTrfInvocation.for_decomposition(
                decomposition, objective_request_budgets=(80,) * len(decomposition.components)
            ),
        )
        executors.append(label)
        return step.execute_method_step(workflow, analysis_values=session.analysis_values)

    try:
        first = compose_and_execute("first", starting)
        outcomes.append(first)
        outcomes.append(compose_and_execute("second", step.require_successor_state(first, session.analysis_values)))
    except direct.DirectTrfConstructionError as failure:
        return ProbeObservation(
            "construction-dependent-stop", starting, session.analysis_values.snapshot(), *lineage,
            construction_failure=failure, constructor_entries=tuple(constructors),
            executor_entries=tuple(executors), outcomes=tuple(outcomes),
        )
    return ProbeObservation(
        "successful-two-step", starting, session.analysis_values.snapshot(), *lineage,
        constructor_entries=tuple(constructors), executor_entries=tuple(executors),
        outcomes=tuple(outcomes),
    )
def _method_probe(  # noqa: C901
    probe_id: str, root: Path, authority: LiveLaneAuthority,
    source_commit: str, lockfile_hash: str,
) -> ProbeObservation:
    session, workflow = _base_workflow(stabilize="resampling" in probe_id)
    starting = session.analysis_values.snapshot()
    context = nullcontext()
    kwargs: dict[str, object] = {}
    checkpoint_entries: list[str] = []
    checkpoint_failures: list[str] = []
    if probe_id == "primary-execution-failure":
        workflow = dataclasses.replace(workflow, publication=_publication(root, source_commit, lockfile_hash))
        context = patch.object(workflow.engine, "project_profiles", side_effect=RuntimeError("qualified primary failure"))
    elif probe_id == "aggregate-materialization-failure":
        failure = direct.RootMaterializationFailure(
            direct.MaterializationTerminal.FAILURE, 1, workflow.engine.compatibility_identity,
            0, 1, direct.TerminalFailure("qualified materialization failure"),
        )
        context = patch("chemex.optimize.grouped_direct_trf.direct_trf_owner.materialize_root_candidate", return_value=failure)
    elif probe_id == "aggregate-acceptance-stop":
        def fail_acceptance(checkpoint: step.MethodStepCheckpoint) -> None:
            if checkpoint is step.MethodStepCheckpoint.AGGREGATE_ACCEPTED:
                checkpoint_entries.append(checkpoint.value)
                failure = RuntimeError("qualified acceptance failure")
                checkpoint_failures.append(type(failure).__name__)
                raise failure
        kwargs["checkpoint_observer"] = fail_acceptance
    elif probe_id == "accepted-commit-rejected":
        configuration = session.parameter_factory.sealed_configuration
        assert configuration is not None
        foreign = AnalysisValues()
        foreign.initialize(starting.model_identity, configuration)
        def reject(accepted: direct.AcceptedFitResult, authority: object, *, problem: direct.OptimizationProblem, **other: object):
            return direct.execute_fit_commit(
                accepted, authority, problem=problem,
                parameterization=other["parameterization"], analysis_values=foreign,
            )
        context = patch("chemex.optimize.method_step.execute_fit_commit", side_effect=reject)
    elif "resampling" in probe_id:
        workflow = dataclasses.replace(
            workflow, derivations=(_resampling_request(workflow),),
            publication=_publication(root, source_commit, lockfile_hash)
        )
        if probe_id == "requested-resampling-summary-failure":
            def fail_summary(evidence: resampling.ResamplingEvidence, _policy: resampling.ResamplingSummaryPolicy):
                failure = resampling.SummaryFailure(
                    evidence.identity, "qualification_summary_failure", "Qualification forced summary failure"
                )
                return resampling.ResamplingSummaryOutcome(resampling.SummaryTerminal.INSUFFICIENT_COVERAGE, failure=failure)
            context = patch("chemex.optimize.method_step.summarize_resampling_evidence", side_effect=fail_summary)
        else:
            mode = probe_id
            def run(accepted: direct.AcceptedFitResult, plan: resampling.ResamplingPlan, **_kwargs: object):
                if mode == "cancelled-resampling-partial-publication":
                    return resampling.execute_resampling_evidence(accepted, plan, cancellation_probe=lambda: True)
                def hook(prepared: object, projected: object):
                    if prepared.request.ordinal == 0:
                        if mode == "worker-resampling-failure":
                            raise RuntimeError("qualified worker failure")
                        raise KeyboardInterrupt
                    return projected
                return resampling.execute_resampling_evidence(accepted, plan, candidate_test_hook=hook)
            context = patch("chemex.optimize.method_step.execute_resampling_evidence", side_effect=run)
    elif probe_id == "committed-publication-failure":
        root.mkdir()
        workflow = dataclasses.replace(workflow, publication=_publication(root, source_commit, lockfile_hash))
    lineage = _lineage(
        probe_id, workflow, authority, source_commit, lockfile_hash
    )
    with context:
        outcome = step.execute_method_step(workflow, analysis_values=session.analysis_values, **kwargs)
    return _observed(
        probe_id, session, starting, outcome, lineage,
        tuple(checkpoint_entries), tuple(checkpoint_failures),
    )


def observe_lifecycle_probes(
    root: Path, *, authority: LiveLaneAuthority,
    source_commit: str, lockfile_hash: str,
) -> tuple[ProbeObservation, ...]:
    """Run fixed operational scenarios without consulting semantic expectations."""

    identifiers = (
        "primary-execution-failure", "aggregate-materialization-failure",
        "aggregate-acceptance-stop", "accepted-commit-rejected",
        "cancelled-resampling-partial-publication", "worker-resampling-failure",
        "requested-resampling-summary-failure",
        "interrupted-resampling-partial-publication", "committed-publication-failure",
    )
    return (observe_two_required_steps(
        no_objective=True, authority=authority, source_commit=source_commit,
        lockfile_hash=lockfile_hash,
    ),) + tuple(
        _method_probe(
            identifier, root / identifier, authority, source_commit, lockfile_hash
        ) for identifier in identifiers
    )


def _facts(observed: ProbeObservation) -> CanonicalBaselineValue:
    outcome = observed.outcome
    primary = None
    commit = None
    resampling_facts = None
    summary_failure = None
    if outcome is not None:
        grouped = outcome.primary_execution
        primary = {
            "terminal": outcome.primary_terminal,
            "component_dispositions": [] if grouped is None else [
                item.disposition.value for item in grouped.components
            ],
        }
        if outcome.commit_operation is not None:
            operation = outcome.commit_operation
            commit = {
                "terminal": operation.terminal.value,
                "failure_category": None if operation.failure is None else operation.failure.category.value,
            }
        if outcome.derivations:
            derivation = outcome.derivations[0]
            operation = derivation.operation
            if isinstance(operation, resampling.ResamplingOperation):
                evidence = operation.evidence
                assert evidence is not None
                resampling_facts = {
                    "terminal": operation.terminal.value,
                    "unstarted_ordinals": list(operation.unstarted_ordinals),
                    "completed_count": evidence.completed_count,
                    "successful_count": evidence.successful_count,
                    "failed_count": evidence.failed_count,
                    "replicate_dispositions": [item.disposition.value for item in evidence.outcomes],
                }
            failure = next(
                (item for item in derivation.artifacts if isinstance(item, resampling.SummaryFailure)), None
            )
            summary_failure = None if failure is None else failure.category
    return CanonicalBaselineValue.from_value({
        "state": {
            "starting_revision": observed.starting_snapshot.revision,
            "ending_revision": observed.ending_snapshot.revision,
        },
        "construction": None if observed.construction_failure is None else {
            "failure_category": type(observed.construction_failure).__name__,
            "constructor_entries": list(observed.constructor_entries),
            "executor_entries": list(observed.executor_entries),
        },
        "method_step": None if outcome is None else outcome.to_record(),
        "primary": primary,
        "checkpoint": {
            "entries": list(observed.checkpoint_entries),
            "failure_category": observed.checkpoint_failures[0] if observed.checkpoint_failures else None,
        },
        "commit": commit,
        "resampling": resampling_facts,
        "summary_failure_category": summary_failure,
        "successor": {
            "denied": observed.successor_denial is not None,
            "downstream_entries": list(observed.downstream_entries),
        },
    })


def capture_lifecycle_probes(
    observations: tuple[ProbeObservation, ...],
) -> LifecycleProbeCapture:
    return LifecycleProbeCapture(tuple(
        LifecycleProbeResult(
            item.probe_id, item.case, item.specification, item.occurrence, _facts(item)
        )
        for item in observations
    ))
# fmt: on
