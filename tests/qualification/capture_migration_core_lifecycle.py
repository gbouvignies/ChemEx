"""Execute #592 probes and return observation-only production facts."""

# ruff: noqa: I001

from __future__ import annotations

# fmt: off
import dataclasses
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch

import numpy as np

import chemex.optimize.direct_trf as direct
import chemex.optimize.method_step as step
import chemex.optimize.native_resampling as resampling
from chemex.baselines import LegacyObservationImplementation, Occurrence, ResultBundle, ResultMember
from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.native_provenance import BaselineReference, ProvenanceEnvironment
from chemex.optimize.grouped_direct_trf import FitDecomposition, GroupedDirectTrfInvocation
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.runtime import AnalysisSession
ROOT = Path(__file__).parents[2]
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
@dataclass(frozen=True, slots=True)
class ProbeObservation:
    """Raw objects and entry facts emitted by one production-seam execution."""

    probe_id: str
    starting_snapshot: AnalysisValuesSnapshot
    ending_snapshot: AnalysisValuesSnapshot
    outcome: step.MethodStepOutcome | None = None
    construction_failure: direct.DirectTrfConstructionError | None = None
    constructor_entries: tuple[str, ...] = ()
    executor_entries: tuple[str, ...] = ()
    successor_snapshot: AnalysisValuesSnapshot | None = None
    successor_denial: ValueError | None = None
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
def _publication(path: Path) -> step.MethodStepPublicationRequest:
    requested = Occurrence(
        "a" * 64, "b" * 64, "c" * 64, "unqualified-local-lane-v1", None,
        ("d" * 64,), "issue-592-lifecycle-probe",
    )
    bundle = ResultBundle.create(
        requested.identity, requested.execution_specification_identity,
        LegacyObservationImplementation(), (ResultMember("result", "e" * 64, 1),),
    )
    return step.MethodStepPublicationRequest(
        path, ProvenanceEnvironment.from_current_process(),
        (BaselineReference.from_occurrence(requested.succeeded(bundle)), BaselineReference.from_result_bundle(bundle)),
    )
def _observed(
    probe_id: str, session: AnalysisSession, starting: AnalysisValuesSnapshot,
    outcome: step.MethodStepOutcome,
) -> ProbeObservation:
    try:
        successor = step.require_successor_state(outcome, session.analysis_values)
    except ValueError as denial:
        return ProbeObservation(
            probe_id, starting, session.analysis_values.snapshot(), outcome,
            successor_denial=denial,
        )
    return ProbeObservation(
        probe_id, starting, session.analysis_values.snapshot(), outcome,
        successor_snapshot=successor,
    )
def observe_two_required_steps(*, no_objective: bool) -> ProbeObservation:
    """Run one unconditional two-step production composition path."""

    session, template = _base_workflow(
        no_objective=no_objective, stabilize=not no_objective
    )
    starting = session.analysis_values.snapshot()
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
            "construction-dependent-stop", starting, session.analysis_values.snapshot(),
            construction_failure=failure, constructor_entries=tuple(constructors),
            executor_entries=tuple(executors), outcomes=tuple(outcomes),
        )
    return ProbeObservation(
        "successful-two-step", starting, session.analysis_values.snapshot(),
        constructor_entries=tuple(constructors), executor_entries=tuple(executors),
        outcomes=tuple(outcomes),
    )
def _method_probe(probe_id: str, root: Path) -> ProbeObservation:  # noqa: C901
    session, workflow = _base_workflow(stabilize="resampling" in probe_id)
    starting = session.analysis_values.snapshot()
    context = nullcontext()
    kwargs: dict[str, object] = {}
    if probe_id == "primary-execution-failure":
        workflow = dataclasses.replace(workflow, publication=_publication(root))
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
                raise RuntimeError("qualified acceptance failure")
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
            workflow, derivations=(_resampling_request(workflow),), publication=_publication(root)
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
        workflow = dataclasses.replace(workflow, publication=_publication(root))
    with context:
        outcome = step.execute_method_step(workflow, analysis_values=session.analysis_values, **kwargs)
    return _observed(probe_id, session, starting, outcome)
def observe_lifecycle_probes(root: Path) -> tuple[ProbeObservation, ...]:
    """Run fixed operational scenarios without consulting semantic expectations."""

    identifiers = (
        "primary-execution-failure", "aggregate-materialization-failure",
        "aggregate-acceptance-stop", "accepted-commit-rejected",
        "cancelled-resampling-partial-publication", "worker-resampling-failure",
        "requested-resampling-summary-failure",
        "interrupted-resampling-partial-publication", "committed-publication-failure",
    )
    return (observe_two_required_steps(no_objective=True),) + tuple(
        _method_probe(identifier, root / identifier) for identifier in identifiers
    )
# fmt: on
