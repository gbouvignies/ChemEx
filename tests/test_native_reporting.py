"""Behavioral qualification tests for native step-root reporting (#608).

The seam is the published method-step directory.  Tests deliberately exercise
real native evaluation, acceptance, and commit artifacts while treating the
renderer implementation as private.
"""

from __future__ import annotations

import json
import tomllib
from argparse import Namespace
from dataclasses import replace
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Literal

import numpy as np
import pytest

import chemex.run_info as run_info_module
from chemex.baselines import (
    LegacyObservationImplementation,
    Occurrence,
    ResultBundle,
    ResultMember,
)
from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.native_provenance import (
    ArtifactReference,
    BaselineReference,
    BudgetRecord,
    PolicyRecord,
    ProvenanceEnvironment,
    PublishedStepReference,
    SeedRecord,
    WorkflowProvenance,
)
from chemex.optimize import native_reporting
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CommitReceipt,
    DirectTrfExecution,
    DirectTrfInvocation,
    FitCommitOperation,
    FitCommitTerminal,
    OptimizationProblem,
    canonical_chi_square,
    committed_values_identity,
    execute_direct_trf,
    execute_fit_commit,
)
from chemex.optimize.native_mcmc import (
    ExpertMcmcPolicy,
    McmcExecutionStage,
    McmcPlan,
    derive_posterior_sample_evidence,
    derive_posterior_summary,
    derive_retained_sample_view,
    execute_mcmc_evidence,
)
from chemex.optimize.native_reporting import (
    CommittedFitPublication,
    ComponentDiagnostic,
    EvaluationPublication,
    McmcPublication,
    ResamplingPublication,
    SuppressedPublication,
    publish_native_results,
)
from chemex.optimize.native_resampling import (
    OperationTerminal,
    ResamplingDatasetManifest,
    ResamplingEvidence,
    ResamplingLifecycle,
    ResamplingPlan,
    ResamplingScheme,
    ResamplingSummaryOutcome,
    SummaryFailure,
    SummaryTerminal,
    execute_resampling_evidence,
    summarize_resampling_evidence,
)
from chemex.optimize.uncertainty import (
    ParameterUnit,
    ResidualVarianceScaling,
    UncertaintyEvidence,
    UncertaintyPolicy,
    compile_constraint_linearization_capabilities,
    derive_uncertainty_evidence,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.spin_system import SpinSystem
from chemex.run_info import (
    CapturedNativeInputs,
    NativeRunInformation,
    RunInformationKind,
    capture_native_inputs,
    classify_run_information,
    read_native_run_information,
    write_native_run_info,
)
from chemex.runtime import AnalysisSession, ExecutionSettings

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


def _captured_inputs(output: Path) -> CapturedNativeInputs:
    return capture_native_inputs(
        Namespace(
            experiments=[EXPERIMENT],
            parameters=[PARAMETERS],
            method=[METHOD],
            output=output,
        ),
        working_directory=ROOT,
    )


def _workflow_provenance(
    *,
    parameterization: ActiveParameterization,
    plan: EvaluationPlan,
    method: Method,
    invocation: DirectTrfInvocation | None = None,
    execution: DirectTrfExecution | None = None,
    evaluation_result: EvaluationResult | None = None,
    uncertainty: UncertaintyEvidence | None = None,
    resampling: tuple[ResamplingPublication, ...] = (),
    mcmc: McmcPublication | None = None,
) -> WorkflowProvenance:
    assert (execution is not None) != (evaluation_result is not None)
    native_execution = execution if execution is not None else evaluation_result
    assert native_execution is not None
    requested = Occurrence(
        "a" * 64,
        "b" * 64,
        "c" * 64,
        "unqualified-local-lane-v1",
        None,
        ("d" * 64,),
        "issue-609-baseline-attempt",
    )
    bundle = ResultBundle.create(
        requested.identity,
        requested.execution_specification_identity,
        LegacyObservationImplementation(),
        (ResultMember("result", "e" * 64, 1),),
    )
    occurrence = requested.succeeded(bundle)
    policies: list[PolicyRecord] = []
    budgets: list[BudgetRecord] = []
    seeds: list[SeedRecord] = []
    if invocation is not None or execution is not None:
        assert invocation is not None
        assert execution is not None
        policies.append(PolicyRecord("primary_direct_trf", invocation.identity))
        budgets.append(
            BudgetRecord(
                "primary_direct_trf_objective_requests",
                invocation.objective_request_budget,
                execution.counters.objective_requests_accepted,
            )
        )
    if uncertainty is not None:
        policies.append(PolicyRecord("uncertainty", uncertainty.source_policy.identity))
    for item in resampling:
        resampling_plan = item.evidence.plan
        name = f"resampling_{resampling_plan.scheme.value}"
        policies.append(PolicyRecord(name, resampling_plan.identity))
        budgets.append(
            BudgetRecord(
                f"{name}_replicates",
                resampling_plan.replicate_count,
                item.evidence.completed_count,
            )
        )
        objective_budget = int(
            dict(resampling_plan.strategy_settings).get(
                "objective_request_budget", "100"
            )
        )
        budgets.append(
            BudgetRecord(
                f"{name}_objective_requests",
                objective_budget * resampling_plan.replicate_count,
                None,
            )
        )
        seeds.append(
            SeedRecord(
                f"{name}_root",
                resampling_plan.root_seed,
                resampling_plan.seed_policy_version,
            )
        )
    if mcmc is not None:
        policy = mcmc.evidence.plan.policy
        policies.append(PolicyRecord("mcmc", policy.identity))
        budgets.append(
            BudgetRecord(
                "mcmc_objective_requests",
                policy.objective_request_budget,
                mcmc.evidence.objective_request_count,
            )
        )
        seeds.append(
            SeedRecord("mcmc_root", policy.root_seed, policy.seed_policy_version)
        )
    return WorkflowProvenance.create(
        parameterization=parameterization,
        plan=plan,
        method=method,
        invocation=invocation,
        native_execution=native_execution,
        policies=tuple(policies),
        budgets=tuple(budgets),
        seeds=tuple(seeds),
        environment=ProvenanceEnvironment(
            chemex_version="2026.8",
            python_version="3.14.5",
            python_implementation="CPython",
            platform="macOS-arm64",
            numpy_version="2.5.1",
            scipy_version="1.17.0",
            emcee_version="3.1.6",
            numerical_libraries=(("blas", "accelerate", "unknown"),),
        ),
        baseline_references=(
            BaselineReference.from_occurrence(occurrence),
            BaselineReference.from_result_bundle(bundle),
        ),
    )


def _committed_fit(
    method: Method | None = None,
    *,
    with_uncertainty: bool = False,
    with_resampling: bool = False,
    with_mcmc: bool = False,
) -> CommittedFitPublication:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    selected_method = read_methods([METHOD])["DEFAULT"] if method is None else method
    parameterization = session.compile_parameterization(
        selected_method,
        experiments.param_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    if with_resampling or with_mcmc:
        source = session.analysis_values.snapshot()
        initial_frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            parameterization.frame_from_snapshot(source),
        )
        initial = engine.new_evaluator().evaluate(initial_frame)
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
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None
    assert parameter_model is not None
    starting_snapshot = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        starting_snapshot,
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=80,
        execution_settings=ExecutionSettings(workers=1, native_threads=1),
    )
    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
    )
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    accepted = outcome.accepted_result
    commit_operation = execute_fit_commit(
        accepted,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )
    assert commit_operation.terminal is FitCommitTerminal.COMMITTED
    assert commit_operation.receipt is not None
    assert commit_operation.committed_snapshot is not None
    receipt = commit_operation.receipt
    committed = commit_operation.committed_snapshot
    uncertainty = None
    if with_uncertainty:
        controlled_id = problem.controlled_ids[0]
        derived_id = "__R1A_B_G2N_H_800_0MHZ"
        policy = UncertaintyPolicy(
            calibration_identity="issue-608-reporting-qualification-v1",
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
            affine_feasibility_policy=("canonical-root-affine-halfspace-zeta-gt-3-v1"),
        )
        compiled = compile_constraint_linearization_capabilities(
            parameterization,
            (controlled_id, derived_id),
            (),
        )
        uncertainty = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=policy,
            constrained_scope=(controlled_id, derived_id),
            constrained_units=(
                (controlled_id, ParameterUnit.RATE_PER_SECOND),
                (derived_id, ParameterUnit.RATE_PER_SECOND),
            ),
            constrained_scales=((controlled_id, 1.0), (derived_id, 1.0)),
            compiled_constraint_linearization=compiled,
            resolved_environment_identity="issue-608-local-environment",
        )
    resampling: tuple[ResamplingPublication, ...] = ()
    if with_resampling:
        observation_count = engine.plan.observation_count
        dataset = ResamplingDatasetManifest(
            engine.plan,
            tuple(
                float(value)
                for value in accepted.evaluation_result.normalized_calculations
            ),
            tuple(False for _index in range(observation_count)),
            tuple("G2N" for _index in range(observation_count)),
            tuple(f"ordinal={index}" for index in range(observation_count)),
        )
        plan = ResamplingPlan.for_accepted(
            accepted,
            dataset=dataset,
            source_problem=problem,
            parameterization=parameterization,
            source_engine=engine,
            scheme=ResamplingScheme.MONTE_CARLO,
            replicate_count=2,
            replicate_structural_identities=("replicate-alpha", "replicate-beta"),
            replicate_component_identities=(("component",), ("component",)),
            root_seed=608,
            output_scope=parameterization.scope_ids,
            output_units=tuple("native" for _id in parameterization.scope_ids),
            minimum_successful_count=2,
            strategy_settings=(("objective_request_budget", "80"),),
        )
        operation = execute_resampling_evidence(accepted, plan)
        assert operation.evidence is not None
        summary = summarize_resampling_evidence(operation.evidence)
        resampling = (ResamplingPublication(operation.evidence, summary),)
    mcmc = None
    if with_mcmc:
        policy = ExpertMcmcPolicy(
            burn_steps=1,
            retained_steps=3,
            walkers=4,
            expert_provenance="issue-608-reporting-qualification",
        ).resolve(dimension=len(problem.controlled_ids), root_seed=608)
        plan = McmcPlan.for_accepted(
            accepted,
            source_problem=problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=tuple(
                (param_id, ParameterUnit.DIMENSIONLESS)
                for param_id in problem.controlled_ids
            ),
        )
        operation = execute_mcmc_evidence(accepted, plan)
        assert operation.evidence is not None
        retained = derive_retained_sample_view(operation.evidence)
        posterior = derive_posterior_sample_evidence(retained, plan.coordinate_units)
        summary = derive_posterior_summary(posterior)
        mcmc = McmcPublication(operation.evidence, posterior, summary)
    return CommittedFitPublication(
        engine.plan,
        selected_method,
        parameterization,
        parameter_model,
        starting_snapshot,
        problem,
        invocation,
        outcome.execution,
        accepted,
        receipt,
        committed,
        commit_operation,
        session.analysis_values,
        provenance=_workflow_provenance(
            parameterization=parameterization,
            plan=engine.plan,
            method=selected_method,
            invocation=invocation,
            execution=outcome.execution,
            uncertainty=uncertainty,
            resampling=resampling,
            mcmc=mcmc,
        ),
        uncertainty=uncertainty,
        resampling=resampling,
        mcmc=mcmc,
    )


def _evaluation_only(method: Method | None = None) -> EvaluationPublication:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = Method(fix=["R1A_A", "PB", "KEX_AB"]) if method is None else method
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
    return EvaluationPublication(
        engine.plan,
        method,
        parameterization,
        result,
        provenance=_workflow_provenance(
            parameterization=parameterization,
            plan=engine.plan,
            method=method,
            evaluation_result=result,
        ),
    )


def _sequential_committed_steps(
    output: Path,
) -> tuple[
    tuple[CommittedFitPublication, CommittedFitPublication],
    tuple[PublishedStepReference, PublishedStepReference],
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    publications: list[CommittedFitPublication] = []
    steps: list[PublishedStepReference] = []
    for method in (
        Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]),
        Method(fit=["KEX_AB"], fix=["R1A_A", "PB"]),
    ):
        parameterization = session.compile_parameterization(
            method,
            experiments.param_ids,
        )
        engine = EvaluationEngine.from_experiments(experiments, parameterization)
        starting = session.analysis_values.snapshot()
        problem = OptimizationProblem.from_native(
            engine.plan,
            parameterization,
            configuration,
            starting,
        )
        invocation = DirectTrfInvocation.for_problem(
            problem,
            objective_request_budget=80,
        )
        outcome = execute_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )
        assert outcome.accepted_result is not None
        assert outcome.commit_authority is not None
        accepted = outcome.accepted_result
        operation = execute_fit_commit(
            accepted,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
        assert operation.receipt is not None
        assert operation.committed_snapshot is not None
        publications.append(
            CommittedFitPublication(
                engine.plan,
                method,
                parameterization,
                parameter_model,
                starting,
                problem,
                invocation,
                outcome.execution,
                accepted,
                operation.receipt,
                operation.committed_snapshot,
                operation,
                session.analysis_values,
                provenance=_workflow_provenance(
                    parameterization=parameterization,
                    plan=engine.plan,
                    method=method,
                    invocation=invocation,
                    execution=outcome.execution,
                ),
            )
        )
        steps.append(
            publish_native_results(
                output / f"STEP-{len(publications)}",
                publications[-1],
            )
        )
    return (publications[0], publications[1]), (steps[0], steps[1])


def _failed_commit_operation() -> tuple[
    FitCommitOperation,
    EvaluationPlan,
    ActiveParameterization,
    OptimizationProblem,
    DirectTrfInvocation,
    DirectTrfExecution,
    WorkflowProvenance,
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = read_methods([METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None
    assert parameter_model is not None
    source = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, source
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=80,
    )
    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
    )
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    session.analysis_values.commit(
        {param_id: source[param_id] for param_id in problem.commit_scope},
        expected=source,
        scope=problem.commit_scope,
    )
    operation = execute_fit_commit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )
    assert operation.terminal is FitCommitTerminal.FAILED
    return (
        operation,
        engine.plan,
        parameterization,
        problem,
        invocation,
        outcome.execution,
        _workflow_provenance(
            parameterization=parameterization,
            plan=engine.plan,
            method=method,
            invocation=invocation,
            execution=outcome.execution,
        ),
    )


def _publication_from_direct_commit_and_fabricated_evidence() -> (
    CommittedFitPublication
):
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    method = read_methods([METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None
    assert parameter_model is not None
    source = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, source
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=80,
    )
    outcome = execute_direct_trf(problem, invocation, parameterization, engine)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame(problem.start, parameterization),
    )
    evaluation = engine.new_evaluator().evaluate(frame)
    assert isinstance(evaluation, EvaluationResult)
    accepted = AcceptedFitResult.for_qualification(
        occurrence_identity="issue-608-direct-commit-attack",
        problem_identity=problem.identity,
        invocation_identity="issue-608-fabricated-invocation",
        execution_identity="issue-608-fabricated-execution",
        materialization_identity="issue-608-fabricated-materialization",
        parameterization_identity=parameterization.identity,
        evaluator_parameterization_identity=parameterization.evaluator_identity,
        source_occurrence_identity=source.occurrence_identity,
        source_revision=source.revision,
        controlled_ids=problem.controlled_ids,
        vector=problem.start,
        chi_square=canonical_chi_square(evaluation.residuals),
        evaluation_result=evaluation,
        commit_scope=problem.commit_scope,
        commit_items=evaluation.resolved_values.ordered_items(),
        origin_context_identity="issue-608-fabricated-origin",
    )
    committed = session.analysis_values.commit(
        dict(accepted.commit_items),
        expected=source,
        scope=problem.commit_scope,
    )
    receipt = CommitReceipt(
        accepted.occurrence_identity,
        accepted.identity,
        problem.identity,
        source.revision,
        committed.revision,
        problem.commit_scope,
        committed_values_identity(committed, problem.commit_scope),
        committed.model_identity,
        committed.configuration_identity,
    )
    fabricated_operation = FitCommitOperation(
        "fabricated-commit-operation",
        accepted.identity,
        accepted.occurrence_identity,
        problem.identity,
        FitCommitTerminal.COMMITTED,
        receipt=receipt,
        committed_snapshot=committed,
    )
    return CommittedFitPublication(
        engine.plan,
        method,
        parameterization,
        parameter_model,
        source,
        problem,
        invocation,
        outcome.execution,
        accepted,
        receipt,
        committed,
        fabricated_operation,
        session.analysis_values,
        provenance=_workflow_provenance(
            parameterization=parameterization,
            plan=engine.plan,
            method=method,
            invocation=invocation,
            execution=outcome.execution,
        ),
    )


def _tree_bytes(path: Path) -> dict[str, bytes]:
    return {
        str(item.relative_to(path)): item.read_bytes()
        for item in sorted(path.rglob("*"))
        if item.is_file()
    }


def _staging_residue(parent: Path, destination_name: str) -> tuple[Path, ...]:
    return tuple(parent.glob(f".{destination_name}.tmp-*"))


def test_committed_native_fit_publishes_only_aggregate_step_root(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"

    publish_native_results(output, publication)

    assert {path.name for path in output.iterdir()} == {
        "Data",
        "Parameters",
        "Plots",
        "Statistics",
        "fit-manifest.toml",
        "statistics.toml",
    }
    assert (output / "Data" / "profile_0000.tsv").is_file()
    assert (output / "Plots" / "profiles.pdf").stat().st_size > 0
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()

    fitted = (output / "Parameters" / "fitted.toml").read_text()
    assert publication.accepted.controlled_ids[0] in fitted
    assert "stderr" not in fitted.lower()
    assert "error" not in fitted.lower()
    assert "±" not in fitted

    restart_path = output / "Parameters" / "restart.toml"
    restart_defaults = read_defaults([restart_path])
    assert len(restart_defaults) == len(publication.parameterization.independent_ids)
    assert all(
        setting.min is not None and setting.max is not None
        for _name, setting in restart_defaults
    )

    manifest = tomllib.loads((output / "fit-manifest.toml").read_text())
    assert manifest["schema_version"] == 1
    assert manifest["lifecycle"] == "committed"
    assert manifest["authority"] == "committed_fit"
    assert manifest["accepted_result_identity"] == publication.accepted.identity
    assert manifest["commit_receipt_identity"] == publication.commit_receipt.identity
    assert manifest["starting_state_kind"] == "workflow_starting_provenance"
    assert manifest["starting_revision"] == 0
    assert manifest["workflow"]["identity"] == publication.provenance.workflow_identity
    assert manifest["method"]["identity"] == publication.provenance.method_identity
    assert (
        manifest["method"]["normalized_sha256"]
        == sha256(publication.provenance.normalized_method_text.encode()).hexdigest()
    )
    assert manifest["selection"] == {
        "model": "2st",
        "include": [profile.identity for profile in publication.plan.profiles],
        "exclude": [],
    }
    assert manifest["execution"] == {"workers": 1, "native_threads": 1}
    assert manifest["budgets"][0] == {
        "name": "primary_direct_trf_objective_requests",
        "limit": publication.primary_invocation.objective_request_budget,
        "used": publication.primary_execution.counters.objective_requests_accepted,
    }
    assert "seeds" not in manifest
    assert manifest["environment"]["numpy_version"] == "2.5.1"
    assert manifest["baseline_references"] == [
        {
            "kind": "occurrence",
            "identity": publication.provenance.baseline_references[0].identity,
            "occurrence_identity": publication.provenance.baseline_references[
                0
            ].occurrence_identity,
            "result_bundle_identity": publication.provenance.baseline_references[
                0
            ].result_bundle_identity,
        },
        {
            "kind": "bundle",
            "identity": publication.provenance.baseline_references[1].identity,
            "occurrence_identity": publication.provenance.baseline_references[
                1
            ].occurrence_identity,
            "result_bundle_identity": publication.provenance.baseline_references[
                1
            ].result_bundle_identity,
        },
    ]
    assert manifest["artifacts"]["Parameters/restart.toml"]["role"] == (
        "committed_restart_state"
    )
    assert manifest["artifacts"]["Parameters/fitted.toml"]["role"] == (
        "report_only_fitted_values"
    )
    for relative_path, artifact in manifest["artifacts"].items():
        artifact_path = output / relative_path
        assert artifact["sha256"] == sha256(artifact_path.read_bytes()).hexdigest()

    statistics = tomllib.loads((output / "statistics.toml").read_text())
    residual_count = len(publication.accepted.evaluation_result.residuals)
    controlled_count = len(publication.accepted.controlled_ids)
    normalization_count = sum(
        profile.is_scaled for profile in publication.plan.profiles
    )
    degrees_of_freedom = residual_count - controlled_count - normalization_count
    assert statistics["retained residual count"] == residual_count
    assert statistics["controlled parameter count"] == controlled_count
    assert statistics["profiled normalization count"] == normalization_count
    assert statistics["nominal residual degrees of freedom"] == degrees_of_freedom
    assert statistics["chi-square"] == pytest.approx(publication.accepted.chi_square)
    assert statistics["reduced-chi-square"] == pytest.approx(
        publication.accepted.chi_square / degrees_of_freedom
    )


def test_many_components_are_diagnostic_views_without_authoritative_copies(
    tmp_path: Path,
) -> None:
    base = _committed_fit(
        Method(fit=["PB", "R1A_A"], fix=["KEX_AB"]),
    )
    publication = CommittedFitPublication(
        base.plan,
        base.method,
        base.parameterization,
        base.parameter_model,
        base.starting_snapshot,
        base.primary_problem,
        base.primary_invocation,
        base.primary_execution,
        base.accepted,
        base.commit_receipt,
        base.committed_snapshot,
        base.commit_operation,
        base.analysis_values,
        provenance=base.provenance,
        components=(
            ComponentDiagnostic(
                "component-alpha",
                "succeeded",
                (base.accepted.controlled_ids[0],),
                local_chi_square=1.25,
            ),
            ComponentDiagnostic(
                "component-beta",
                "succeeded",
                base.accepted.controlled_ids[1:],
                local_chi_square=2.5,
            ),
        ),
    )
    output = tmp_path / "STEP"

    publish_native_results(output, publication)

    components = json.loads((output / "Components" / "index.json").read_text())
    manifest = tomllib.loads((output / "fit-manifest.toml").read_text())
    assert manifest["artifacts"]["Components/index.json"]["role"] == (
        "diagnostic_provenance"
    )
    assert [item["ordinal"] for item in components["components"]] == [0, 1]
    assert [item["identity"] for item in components["components"]] == [
        "component-alpha",
        "component-beta",
    ]
    assert components["authority"] == "diagnostic_only"
    manifest = tomllib.loads((output / "fit-manifest.toml").read_text())
    assert [item["identity"] for item in manifest["components"]] == [
        "component-alpha",
        "component-beta",
    ]
    assert all(
        item["location"] == "Components/index.json"
        and item["authority"] == "diagnostic_only"
        for item in manifest["components"]
    )
    assert not list((output / "Components").rglob("Parameters"))
    assert not list((output / "Components").rglob("Data"))
    assert not list((output / "Components").rglob("Plots"))
    assert not list((output / "Components").rglob("Statistics"))
    assert (output / "Parameters" / "fitted.toml").is_file()


def test_committed_publication_rejects_a_receipt_without_its_snapshot_witness(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    reconstructed_receipt = replace(publication.commit_receipt)

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, commit_receipt=reconstructed_receipt),
        )
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_a_manually_reconstructed_receipt(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    receipt = publication.commit_receipt
    reconstructed = CommitReceipt(
        receipt.accepted_occurrence_identity,
        receipt.accepted_result_identity,
        receipt.problem_identity,
        receipt.old_revision,
        receipt.new_revision,
        receipt.scope,
        receipt.committed_value_identity,
        receipt.model_identity,
        receipt.configuration_identity,
    )

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, commit_receipt=reconstructed),
        )
    assert not (tmp_path / "STEP").exists()


def test_direct_analysis_values_commit_cannot_fabricate_publication_authority(
    tmp_path: Path,
) -> None:
    publication = _publication_from_direct_commit_and_fabricated_evidence()

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(tmp_path / "STEP", publication)
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_another_occurrence_receipt(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    other = _committed_fit()

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, commit_receipt=other.commit_receipt),
        )
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_a_stale_committed_revision(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    current = publication.analysis_values.snapshot()
    publication.analysis_values.commit(
        {param_id: current[param_id] for param_id in publication.commit_receipt.scope},
        expected=current,
        scope=publication.commit_receipt.scope,
    )

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(tmp_path / "STEP", publication)
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_a_copied_successful_commit_operation(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    copied_operation = replace(publication.commit_operation)

    with pytest.raises(ValueError, match="exact successful commit operation"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, commit_operation=copied_operation),
        )
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_forged_execution_budget_provenance(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    source = publication.provenance.budgets[0]
    forged_provenance = replace(
        publication.provenance,
        budgets=(BudgetRecord(source.name, source.limit, source.used - 1),),
    )

    with pytest.raises(ValueError, match="execution record"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, provenance=forged_provenance),
        )
    assert not (tmp_path / "STEP").exists()


def test_committed_publication_rejects_a_foreign_normalized_method(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    other = _committed_fit(Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]))
    forged_provenance = replace(
        publication.provenance,
        normalized_method_text=other.provenance.normalized_method_text,
    )

    with pytest.raises(ValueError, match="method and selection provenance"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, provenance=forged_provenance),
        )
    assert not (tmp_path / "STEP").exists()


def test_normalized_method_rejects_noncanonical_nested_payload() -> None:
    publication = _committed_fit()
    with pytest.raises(ValueError, match="Normalized method semantics"):
        replace(
            publication.provenance,
            normalized_method_text=(
                'FORMAT_VERSION = 2\n[METHOD]\nCANONICAL_JSON = "{}"\n'
            ),
        )


def test_workflow_provenance_rejects_substituted_baseline_occurrence() -> None:
    publication = _committed_fit()
    occurrence, bundle = publication.provenance.baseline_references
    with pytest.raises(ValueError, match="linked occurrence/bundle"):
        replace(
            publication.provenance,
            baseline_references=(replace(occurrence, identity="f" * 64), bundle),
        )


@pytest.mark.parametrize(
    "change",
    ("workflow", "grouping", "selection", "policy", "budget", "seed", "execution"),
)
def test_execution_identity_covers_every_non_method_execution_choice(
    change: str,
) -> None:
    publication = _committed_fit()
    base = publication.provenance
    changes: dict[str, object]
    if change == "workflow":
        changes = {"workflow_type": "grid"}
    elif change == "grouping":
        changes = {"grouping_topology": (("separate", base.grouping_topology[0][1]),)}
    elif change == "selection":
        changes = {"selection": replace(base.selection, include=("foreign",))}
    elif change == "policy":
        changes = {"policies": (PolicyRecord("primary_direct_trf", "foreign"),)}
    elif change == "budget":
        budget = base.budgets[0]
        changes = {
            "budgets": (BudgetRecord(budget.name, budget.limit + 1, budget.used),)
        }
    elif change == "seed":
        changes = {"seeds": (SeedRecord("root", 609, "seed-policy"),)}
    else:
        changes = {"execution": ExecutionSettings(workers=3, native_threads=1)}

    changed = replace(base, **changes)

    assert changed.method_identity == base.method_identity
    assert changed.execution_identity != base.execution_identity


@pytest.mark.parametrize(
    ("changes", "description"),
    (
        ({"workflow_type": "grid"}, "GRID workflow"),
        ({"workflow_type": "de"}, "DE workflow"),
        (
            {"grouping_topology": (("foreign-component", ("FOREIGN",)),)},
            "grouping topology",
        ),
        (
            {"execution": ExecutionSettings(workers=999, native_threads=1)},
            "worker count",
        ),
        (
            {"execution": ExecutionSettings(workers=1, native_threads=999)},
            "native thread count",
        ),
    ),
)
def test_committed_publication_rejects_recomputed_false_execution_provenance(
    tmp_path: Path,
    changes: dict[str, object],
    description: str,
) -> None:
    publication = _committed_fit()
    forged_provenance = replace(publication.provenance, **changes)

    with pytest.raises(ValueError, match="execution provenance"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, provenance=forged_provenance),
        )
    assert not (tmp_path / "STEP").exists(), description


def test_committed_publication_rejects_unsuccessful_component_diagnostics(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    component = ComponentDiagnostic(
        "failed-component",
        "failed",
        publication.accepted.controlled_ids,
    )

    with pytest.raises(ValueError, match="exactly partition"):
        publish_native_results(
            tmp_path / "STEP",
            replace(publication, components=(component,)),
        )


def test_zero_component_evaluation_publishes_no_fitted_estimate(tmp_path: Path) -> None:
    publication = _evaluation_only()
    output = tmp_path / "STEP"

    publish_native_results(output, publication)

    assert (output / "Data" / "profile_0000.tsv").is_file()
    assert (output / "Plots" / "profiles.pdf").is_file()
    assert (output / "Statistics").is_dir()
    assert not (output / "Parameters" / "fitted.toml").exists()
    assert not (output / "Components").exists()
    manifest = tomllib.loads((output / "fit-manifest.toml").read_text())
    assert manifest["lifecycle"] == "successful_no_state_change"
    assert manifest["authority"] == "evaluation_only"
    assert "accepted_result_identity" not in manifest
    assert "commit_receipt_identity" not in manifest
    assert "Parameters/restart.toml" not in manifest["artifacts"]
    statistics = tomllib.loads((output / "statistics.toml").read_text())
    residual_count = len(publication.evaluation_result.residuals)
    normalization_count = sum(
        profile.is_scaled for profile in publication.plan.profiles
    )
    assert statistics["controlled parameter count"] == 0
    assert statistics["nominal residual degrees of freedom"] == (
        residual_count - normalization_count
    )


def test_covariance_and_constrained_evidence_are_separate_from_central_reports(
    tmp_path: Path,
) -> None:
    publication = _committed_fit(with_uncertainty=True)
    output = tmp_path / "STEP"

    publish_native_results(output, publication)

    covariance = json.loads(
        (output / "Statistics" / "Covariance" / "evidence.json").read_text()
    )
    constrained = json.loads(
        (output / "Statistics" / "Constrained" / "evidence.json").read_text()
    )
    assert covariance["artifact_type"] == "native_covariance_evidence"
    assert constrained["artifact_type"] == "native_constrained_evidence"
    assert covariance["accepted_result_identity"] == publication.accepted.identity
    assert constrained["accepted_result_identity"] == publication.accepted.identity

    central_reports = "\n".join(
        path.read_text() for path in (output / "Parameters").glob("*.toml")
    ).lower()
    assert "stderr" not in central_reports
    assert "standard error" not in central_reports
    assert "±" not in central_reports


def test_resampling_and_posterior_evidence_keep_family_specific_semantics(
    tmp_path: Path,
) -> None:
    resampling_publication = _committed_fit(with_resampling=True)
    mcmc_publication = _committed_fit(
        Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]),
        with_mcmc=True,
    )
    resampling_output = tmp_path / "RESAMPLING"
    mcmc_output = tmp_path / "MCMC"

    publish_native_results(resampling_output, resampling_publication)
    publish_native_results(mcmc_output, mcmc_publication)

    resampling = json.loads(
        (
            resampling_output / "Statistics" / "Resampling" / "MC" / "summary.json"
        ).read_text()
    )
    posterior = json.loads(
        (mcmc_output / "Statistics" / "MCMC" / "posterior-summary.json").read_text()
    )
    assert resampling["artifact_type"] == "native_resampling_summary"
    assert resampling["scheme"] == "mc"
    resampling_manifest = tomllib.loads(
        (resampling_output / "fit-manifest.toml").read_text()
    )
    assert resampling_manifest["evidence"][0]["family"] == "resampling:mc"
    assert resampling_manifest["evidence"][0]["validity"] == (
        "required_claims_satisfied"
    )
    assert all("standard_deviation" in item for item in resampling["distributions"])
    assert posterior["artifact_type"] == "native_posterior_summary"
    assert all(
        "posterior_standard_deviation" in item
        for item in posterior["parameter_summaries"]
    )
    rendered = json.dumps((resampling, posterior)).lower()
    assert "stderr" not in rendered
    assert (
        resampling_output / "Statistics" / "Resampling" / "MC" / "evidence.json"
    ).is_file()
    assert (mcmc_output / "Statistics" / "MCMC" / "evidence.json").is_file()


@pytest.mark.parametrize(
    ("terminal", "truncate", "expected_lifecycle"),
    (
        (OperationTerminal.COMPLETED, True, ResamplingLifecycle.PARTIAL),
        (OperationTerminal.CANCELLED, False, ResamplingLifecycle.CANCELLED),
        (OperationTerminal.INTERRUPTED, False, ResamplingLifecycle.INTERRUPTED),
    ),
)
def test_incomplete_resampling_is_rejected_from_ordinary_statistics(
    tmp_path: Path,
    terminal: OperationTerminal,
    truncate: bool,
    expected_lifecycle: ResamplingLifecycle,
) -> None:
    publication = _committed_fit(with_resampling=True)
    source = publication.resampling[0].evidence
    outcomes = source.outcomes[:1] if truncate else source.outcomes
    incomplete = ResamplingEvidence(
        source.plan,
        source.population_identity,
        outcomes,
        terminal,
    )
    assert incomplete.lifecycle is expected_lifecycle
    summary = summarize_resampling_evidence(incomplete)

    with pytest.raises(ValueError, match="completed resampling evidence"):
        publish_native_results(
            tmp_path / "STEP",
            replace(
                publication,
                resampling=(ResamplingPublication(incomplete, summary),),
            ),
        )
    assert not (tmp_path / "STEP").exists()


def test_failed_resampling_summary_is_rejected_from_ordinary_statistics(
    tmp_path: Path,
) -> None:
    publication = _committed_fit(with_resampling=True)
    evidence = publication.resampling[0].evidence
    failure = SummaryFailure(
        evidence.identity,
        "failed_summary",
        "summary qualification failed",
    )
    failed = ResamplingSummaryOutcome(
        SummaryTerminal.SOURCE_INVALID,
        failure=failure,
    )

    with pytest.raises(ValueError, match="completed resampling evidence"):
        publish_native_results(
            tmp_path / "STEP",
            replace(
                publication,
                resampling=(ResamplingPublication(evidence, failed),),
            ),
        )
    assert not (tmp_path / "STEP").exists()


@pytest.mark.parametrize(
    ("terminal", "truncate", "suppressed_lifecycle"),
    (
        (OperationTerminal.COMPLETED, True, "failed"),
        (OperationTerminal.CANCELLED, False, "cancelled"),
        (OperationTerminal.INTERRUPTED, False, "interrupted"),
    ),
)
def test_incomplete_resampling_is_retained_only_as_partial_evidence(
    tmp_path: Path,
    terminal: OperationTerminal,
    truncate: bool,
    suppressed_lifecycle: Literal["failed", "cancelled", "interrupted"],
) -> None:
    publication = _committed_fit(with_resampling=True)
    source = publication.resampling[0].evidence
    outcomes = source.outcomes[:1] if truncate else source.outcomes
    incomplete = ResamplingEvidence(
        source.plan,
        source.population_identity,
        outcomes,
        terminal,
    )
    partial_resampling = ResamplingPublication(
        incomplete,
        summarize_resampling_evidence(incomplete),
    )
    suppressed = SuppressedPublication(
        suppressed_lifecycle,
        incomplete,
        publication.plan,
        publication.method,
        publication.parameterization,
        provenance=_workflow_provenance(
            parameterization=publication.parameterization,
            plan=publication.plan,
            method=publication.method,
            invocation=publication.primary_invocation,
            execution=publication.primary_execution,
            resampling=(partial_resampling,),
        ),
        primary_invocation=publication.primary_invocation,
        primary_execution=publication.primary_execution,
        primary_problem=publication.primary_problem,
        accepted_result_identity=publication.accepted.identity,
        accepted_occurrence_identity=publication.accepted.occurrence_identity,
        partial_resampling=(incomplete,),
    )

    publish_native_results(tmp_path / "STEP", suppressed)

    assert (tmp_path / "STEP" / "PartialEvidence" / "Resampling").is_dir()
    assert not (tmp_path / "STEP" / "Statistics").exists()


@pytest.mark.parametrize(
    "alteration",
    (
        "included_ordinals",
        "covariance",
        "correlations",
        "percentile",
        "claims",
        "source_identity",
    ),
)
def test_altered_resampling_summary_is_rejected_before_staging(
    tmp_path: Path,
    alteration: str,
) -> None:
    publication = _committed_fit(with_resampling=True)
    summary = publication.resampling[0].summary.summary
    assert summary is not None
    if alteration == "included_ordinals":
        object.__setattr__(
            summary, "included_ordinals", (*summary.included_ordinals, 99)
        )
    elif alteration == "covariance":
        object.__setattr__(summary, "covariance", summary.covariance[:-1])
    elif alteration == "correlations":
        object.__setattr__(summary, "correlations", summary.correlations[:-1])
    elif alteration == "percentile":
        distribution = summary.distributions[0]
        object.__setattr__(
            distribution,
            "percentile_95_upper",
            distribution.percentile_95_upper + 1.0,
        )
    elif alteration == "claims":
        object.__setattr__(summary, "claims", summary.claims[:-1])
    else:
        object.__setattr__(summary, "evidence_identity", "foreign-evidence")

    with pytest.raises(ValueError, match="differs from canonical source evidence"):
        publish_native_results(tmp_path / "STEP", publication)
    assert not (tmp_path / "STEP").exists()


def test_round_tripped_resampling_summary_remains_publishable(tmp_path: Path) -> None:
    publication = _committed_fit(with_resampling=True)
    item = publication.resampling[0]
    summary = item.summary.summary
    assert summary is not None
    restored = type(summary).from_record(
        summary.to_record(),
        item.evidence,
        summary.policy,
    )
    outcome = replace(item.summary, summary=restored)

    publish_native_results(
        tmp_path / "STEP",
        replace(
            publication,
            resampling=(ResamplingPublication(item.evidence, outcome),),
        ),
    )

    assert (tmp_path / "STEP" / "Statistics" / "Resampling").is_dir()


def test_existing_destination_collision_preserves_authoritative_tree(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"
    (output / "existing").mkdir(parents=True)
    (output / "existing" / "sentinel.bin").write_bytes(b"authoritative-old-tree")
    before = _tree_bytes(output)

    with pytest.raises(FileExistsError, match="destination exists"):
        publish_native_results(output, publication)

    assert _tree_bytes(output) == before
    assert not _staging_residue(tmp_path, output.name)


@pytest.mark.parametrize("destination_kind", ["directory", "file", "symlink"])
def test_destination_appearing_after_staging_preserves_competing_object(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    destination_kind: Literal["directory", "file", "symlink"],
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"
    symlink_target = tmp_path / "competing-target"
    original_write_components = native_reporting.write_components

    def finish_staging_then_create_destination(
        path: Path,
        components: tuple[ComponentDiagnostic, ...],
    ) -> None:
        original_write_components(path, components)
        if destination_kind == "directory":
            sentinel = output / "competing" / "sentinel.bin"
            sentinel.parent.mkdir(parents=True)
            sentinel.write_bytes(b"competing-authoritative-tree")
        elif destination_kind == "file":
            output.write_bytes(b"competing-authoritative-file")
        else:
            symlink_target.mkdir()
            (symlink_target / "sentinel.bin").write_bytes(
                b"competing-authoritative-symlink-target"
            )
            output.symlink_to(symlink_target, target_is_directory=True)

    monkeypatch.setattr(
        native_reporting,
        "write_components",
        finish_staging_then_create_destination,
    )

    with pytest.raises(FileExistsError, match="destination exists"):
        publish_native_results(output, publication)

    if destination_kind == "directory":
        assert {str(item.relative_to(output)) for item in output.rglob("*")} == {
            "competing",
            "competing/sentinel.bin",
        }
        assert _tree_bytes(output) == {
            "competing/sentinel.bin": b"competing-authoritative-tree"
        }
    elif destination_kind == "file":
        assert output.read_bytes() == b"competing-authoritative-file"
    else:
        assert output.is_symlink()
        assert output.readlink() == symlink_target
        assert {
            str(item.relative_to(symlink_target)) for item in symlink_target.rglob("*")
        } == {"sentinel.bin"}
        assert _tree_bytes(symlink_target) == {
            "sentinel.bin": b"competing-authoritative-symlink-target"
        }
    assert not _staging_residue(tmp_path, output.name)


def test_mid_render_failure_removes_staging_without_publishing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"

    def fail_statistics(*_args: object, **_kwargs: object) -> None:
        raise OSError("injected statistics rendering failure")

    monkeypatch.setattr(native_reporting, "write_statistics", fail_statistics)

    with pytest.raises(OSError, match="injected statistics rendering failure"):
        publish_native_results(output, publication)

    assert not output.exists()
    assert not _staging_residue(tmp_path, output.name)


def test_final_rename_failure_removes_complete_staging_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"

    def fail_rename(_source: Path, _destination: Path) -> None:
        raise OSError("injected atomic rename failure")

    monkeypatch.setattr(native_reporting, "publish_directory_noreplace", fail_rename)

    with pytest.raises(OSError, match="injected atomic rename failure"):
        publish_native_results(output, publication)

    assert not output.exists()
    assert not _staging_residue(tmp_path, output.name)


def test_final_publication_primitive_receives_the_unchanged_validated_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"
    original = native_reporting.publish_directory_noreplace
    candidate_at_atomic_boundary: dict[str, bytes] = {}

    def observe_atomic_boundary(staging: Path, destination: Path) -> None:
        candidate_at_atomic_boundary.update(_tree_bytes(staging))
        original(staging, destination)

    monkeypatch.setattr(
        native_reporting,
        "publish_directory_noreplace",
        observe_atomic_boundary,
    )

    publish_native_results(output, publication)

    assert candidate_at_atomic_boundary
    assert _tree_bytes(output) == candidate_at_atomic_boundary
    assert not _staging_residue(tmp_path, output.name)


@pytest.mark.parametrize(
    "relative_path",
    (
        "statistics.toml",
        "Parameters/restart.toml",
        "Parameters/restart-provenance.json",
        "Parameters/fitted.toml",
    ),
)
def test_staged_artifact_mutation_after_manifest_prevents_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative_path: str,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"
    original = native_reporting._write_committed_manifest

    def write_then_mutate(
        staging: Path,
        source: CommittedFitPublication,
        occurrence_identity: str,
        restart: object,
    ) -> object:
        result = original(staging, source, occurrence_identity, restart)
        (staging / relative_path).write_bytes(b"mutated after manifest\n")
        return result

    monkeypatch.setattr(
        native_reporting,
        "_write_committed_manifest",
        write_then_mutate,
    )

    with pytest.raises(ValueError, match="changed after hashing"):
        publish_native_results(output, publication)
    assert not output.exists()
    assert not _staging_residue(tmp_path, output.name)


def test_uncommitted_workflow_suppresses_output_and_downstream_evidence(
    tmp_path: Path,
) -> None:
    committed = _committed_fit(
        Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]),
        with_mcmc=True,
    )
    assert committed.mcmc is not None
    cancellation = CancellationToken()

    def cancel_after_initialization(
        stage: McmcExecutionStage,
        _complete_state_count: int,
    ) -> None:
        if stage is McmcExecutionStage.BEFORE_TRANSITION:
            cancellation.cancel()

    operation = execute_mcmc_evidence(
        committed.accepted,
        committed.mcmc.evidence.plan,
        cancellation=cancellation,
        checkpoint_observer=cancel_after_initialization,
    )
    assert operation.evidence is not None
    foreign_occurrence = SuppressedPublication(
        "cancelled",
        operation.evidence,
        committed.plan,
        committed.method,
        committed.parameterization,
        provenance=_workflow_provenance(
            parameterization=committed.parameterization,
            plan=committed.plan,
            method=committed.method,
            invocation=committed.primary_invocation,
            execution=committed.primary_execution,
            mcmc=replace(committed.mcmc, evidence=operation.evidence),
        ),
        primary_invocation=committed.primary_invocation,
        primary_execution=committed.primary_execution,
        primary_problem=committed.primary_problem,
        accepted_result_identity=committed.accepted.identity,
        accepted_occurrence_identity="foreign-accepted-occurrence",
        partial_mcmc=operation.evidence,
    )
    with pytest.raises(ValueError, match="genuine partial MCMC"):
        publish_native_results(tmp_path / "FOREIGN", foreign_occurrence)

    (
        failed_commit,
        failed_plan,
        failed_parameterization,
        failed_problem,
        failed_invocation,
        failed_execution,
        failed_provenance,
    ) = _failed_commit_operation()
    forged_commit = replace(failed_commit, occurrence_identity="forged-occurrence")
    forged = SuppressedPublication(
        "accepted_uncommitted",
        forged_commit,
        failed_plan,
        read_methods([METHOD])["DEFAULT"],
        failed_parameterization,
        provenance=failed_provenance,
        primary_invocation=failed_invocation,
        primary_execution=failed_execution,
        primary_problem=failed_problem,
        accepted_result_identity=forged_commit.accepted_result_identity,
        accepted_occurrence_identity=forged_commit.accepted_occurrence_identity,
    )
    with pytest.raises(ValueError, match="operation provenance disagree"):
        publish_native_results(tmp_path / "FORGED", forged)

    invalid = SuppressedPublication(
        "accepted_uncommitted",
        failed_commit,
        failed_plan,
        read_methods([METHOD])["DEFAULT"],
        failed_parameterization,
        provenance=failed_provenance,
        primary_invocation=failed_invocation,
        primary_execution=failed_execution,
        primary_problem=failed_problem,
        accepted_result_identity=failed_commit.accepted_result_identity,
        accepted_occurrence_identity=failed_commit.accepted_occurrence_identity,
        partial_mcmc=operation.evidence,
    )
    output = tmp_path / "STEP"

    with pytest.raises(ValueError, match="cannot start downstream evidence"):
        publish_native_results(output, invalid)
    assert not output.exists()

    publication = SuppressedPublication(
        "accepted_uncommitted",
        failed_commit,
        failed_plan,
        read_methods([METHOD])["DEFAULT"],
        failed_parameterization,
        provenance=failed_provenance,
        primary_invocation=failed_invocation,
        primary_execution=failed_execution,
        primary_problem=failed_problem,
        accepted_result_identity=failed_commit.accepted_result_identity,
        accepted_occurrence_identity=failed_commit.accepted_occurrence_identity,
    )

    foreign = _committed_fit(Method(fit=["PB"], fix=["R1A_A", "KEX_AB"]))
    with pytest.raises(ValueError, match="method and selection provenance"):
        publish_native_results(
            tmp_path / "FOREIGN-SCOPE",
            replace(publication, parameterization=foreign.parameterization),
        )

    published_step = publish_native_results(output, publication)

    with pytest.raises(ValueError, match="exact live publication"):
        replace(
            published_step,
            lifecycle="committed",
            authority="committed_fit",
        )
    with pytest.raises(ValueError, match="exact live publication"):
        replace(
            published_step,
            artifacts=(
                *published_step.artifacts,
                ArtifactReference(
                    "Parameters/restart.toml",
                    "committed_restart_state",
                    "0" * 64,
                ),
            ),
        )

    assert {path.name for path in output.iterdir()} == {
        "Diagnostics",
        "fit-manifest.toml",
    }
    diagnostics = json.loads((output / "Diagnostics" / "outcome.json").read_text())
    assert diagnostics["operation"]["identity"] == failed_commit.identity
    assert diagnostics["operation"]["terminal"] == "failed"
    manifest = tomllib.loads((output / "fit-manifest.toml").read_text())
    assert manifest["lifecycle"] == "accepted_uncommitted"
    assert manifest["authority"] == "diagnostic_only"
    assert manifest["operation_identity"] == failed_commit.identity
    assert (
        manifest["artifacts"]["Diagnostics/outcome.json"]["sha256"]
        == sha256((output / "Diagnostics" / "outcome.json").read_bytes()).hexdigest()
    )
    assert not (output / "PartialEvidence").exists()
    for normal_name in (
        "Parameters",
        "Data",
        "Plots",
        "Statistics",
        "statistics.toml",
        "Components",
    ):
        assert not (output / normal_name).exists()

    assert published_step.independent_ids == failed_parameterization.independent_ids
    write_native_run_info(
        Namespace(
            experiments=[EXPERIMENT],
            parameters=[PARAMETERS],
            method=[METHOD],
            output=tmp_path,
        ),
        NativeRunInformation(
            invocation_identity="issue-609-suppressed-only-run",
            inputs=_captured_inputs(tmp_path),
            parameter_model=committed.parameter_model,
            starting_snapshot=committed.starting_snapshot,
            steps=(published_step,),
        ),
        working_directory=ROOT,
    )
    run = tomllib.loads((tmp_path / "run_info" / "run.toml").read_text())
    assert run["steps"][0]["lifecycle"] == "accepted_uncommitted"


def test_native_run_information_archives_replay_complete_provenance(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    published_step = publish_native_results(output / "STEP", publication)
    args = Namespace(
        experiments=[EXPERIMENT],
        parameters=[PARAMETERS],
        method=[METHOD],
        output=output,
    )

    write_native_run_info(
        args,
        NativeRunInformation(
            invocation_identity="issue-609-invocation",
            inputs=_captured_inputs(output),
            parameter_model=publication.parameter_model,
            starting_snapshot=publication.starting_snapshot,
            steps=(published_step,),
        ),
        argv=("chemex", "fit", "--model", "2st"),
        working_directory=ROOT,
        timestamp=datetime(2026, 8, 12, 18, 0, tzinfo=UTC),
    )

    run_info = output / "run_info"
    run = tomllib.loads((run_info / "run.toml").read_text())
    assert run["schema_version"] == 2
    assert run["schema_kind"] == "native_product_run_information"
    assert run["invocation_identity"] == "issue-609-invocation"
    assert run["starting_state"]["kind"] == "starting_independent_state"
    assert run["starting_state"]["revision"] == 0
    assert (
        run["normalized_methods"][0]["sha256"]
        == sha256(publication.provenance.normalized_method_text.encode()).hexdigest()
    )
    assert run["workflow_records"]["path"] == "workflows.json"
    assert "lmfit" not in (run_info / "run.toml").read_text().lower()
    assert "numdifftools" not in (run_info / "run.toml").read_text().lower()

    copied_experiment = run["inputs"]["experiments"][0]
    archived_experiment = run_info / copied_experiment["copied_path"]
    assert archived_experiment.read_bytes() == EXPERIMENT.read_bytes()
    assert copied_experiment["sha256"] == sha256(EXPERIMENT.read_bytes()).hexdigest()

    starting = read_defaults([run_info / "parameters_used.toml"])
    restart = read_defaults([output / "STEP" / "Parameters" / "restart.toml"])
    assert len(starting) == len(publication.parameterization.independent_ids)
    assert len(restart) == len(publication.parameterization.independent_ids)
    assert (run_info / "parameters_used.toml").read_bytes() != (
        output / "STEP" / "Parameters" / "restart.toml"
    ).read_bytes()

    workflows = json.loads((run_info / "workflows.json").read_text())
    assert workflows["workflows"][0]["identity"] == (
        publication.provenance.workflow_identity
    )
    assert workflows["workflows"][0]["outcome"]["lifecycle"] == "committed"
    assert workflows["workflows"][0]["manifest"] == {
        "path": "STEP/fit-manifest.toml",
        "identity": published_step.manifest_identity,
        "sha256": published_step.manifest_sha256,
    }
    assert workflows["workflows"][0]["artifacts"]

    historical = tmp_path / "historical-v1.toml"
    historical.write_text("schema_version = 1\n", encoding="utf-8")
    assert classify_run_information(historical) is RunInformationKind.HISTORICAL_V1
    assert classify_run_information(run_info / "run.toml") is (
        RunInformationKind.NATIVE_V2
    )

    original_tree = _tree_bytes(run_info)
    with pytest.raises(FileExistsError, match="destination exists"):
        write_native_run_info(
            args,
            NativeRunInformation(
                invocation_identity="issue-609-second-invocation",
                inputs=_captured_inputs(output),
                parameter_model=publication.parameter_model,
                starting_snapshot=publication.starting_snapshot,
                steps=(published_step,),
            ),
            working_directory=ROOT,
        )
    assert _tree_bytes(run_info) == original_tree
    assert not list(output.glob(".run_info-*"))


def test_native_run_information_archives_each_method_section(
    tmp_path: Path,
) -> None:
    run_state = _committed_fit()
    first = _evaluation_only(Method(fix=["R1A_A", "PB", "KEX_AB"]))
    second = _evaluation_only(Method(fitmethod="trf", fix=["R1A_A", "PB", "KEX_AB"]))
    output = tmp_path / "Output"
    steps = (
        publish_native_results(output / "STEP-1", first),
        publish_native_results(output / "STEP-2", second),
    )
    write_native_run_info(
        Namespace(output=output),
        NativeRunInformation(
            invocation_identity="issue-609-multi-method",
            inputs=_captured_inputs(output),
            parameter_model=run_state.parameter_model,
            starting_snapshot=run_state.starting_snapshot,
            steps=steps,
        ),
        working_directory=ROOT,
    )

    run = tomllib.loads((output / "run_info" / "run.toml").read_text())
    assert len(run["normalized_methods"]) == 2
    assert {item["method_identity"] for item in run["normalized_methods"]} == {
        first.provenance.method_identity,
        second.provenance.method_identity,
    }
    assert read_native_run_information(
        output / "run_info"
    ).published_step_identities == (
        steps[0].identity,
        steps[1].identity,
    )


def test_native_run_information_walks_sequential_committed_revisions(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    publications, steps = _sequential_committed_steps(output)
    first = publications[0]
    assert steps[0].lifecycle == steps[1].lifecycle == "committed"
    write_native_run_info(
        Namespace(output=output),
        NativeRunInformation(
            invocation_identity="issue-609-sequential-commits",
            inputs=_captured_inputs(output),
            parameter_model=first.parameter_model,
            starting_snapshot=first.starting_snapshot,
            steps=steps,
        ),
        working_directory=ROOT,
    )

    assert read_native_run_information(
        output / "run_info"
    ).published_step_identities == (
        steps[0].identity,
        steps[1].identity,
    )


def test_native_run_information_rejects_changed_step_artifacts_without_replacement(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    published_step = publish_native_results(output / "STEP", publication)
    run_info = output / "run_info"
    run_info.mkdir()
    existing = run_info / "run.toml"
    existing.write_text("schema_version = 1\n", encoding="utf-8")
    (output / "STEP" / "statistics.toml").write_text(
        'tampered = "after publication"\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="artifact changed"):
        write_native_run_info(
            Namespace(
                experiments=[EXPERIMENT],
                parameters=[PARAMETERS],
                method=[METHOD],
                output=output,
            ),
            NativeRunInformation(
                invocation_identity="issue-609-tampered-invocation",
                inputs=_captured_inputs(output),
                parameter_model=publication.parameter_model,
                starting_snapshot=publication.starting_snapshot,
                steps=(published_step,),
            ),
            working_directory=ROOT,
        )

    assert existing.read_text(encoding="utf-8") == "schema_version = 1\n"
    assert not list(output.glob(".run_info-*"))


def test_native_run_information_rejects_unlisted_step_artifact(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    (output / "STEP" / "foreign.bin").write_bytes(b"not catalogued")

    with pytest.raises(ValueError, match="catalogue is incomplete"):
        write_native_run_info(
            Namespace(output=output),
            NativeRunInformation(
                invocation_identity="issue-609-extra-artifact",
                inputs=_captured_inputs(output),
                parameter_model=publication.parameter_model,
                starting_snapshot=publication.starting_snapshot,
                steps=(step,),
            ),
            working_directory=ROOT,
        )
    assert not (output / "run_info").exists()


@pytest.mark.parametrize("field", ("lifecycle", "authority", "artifacts"))
def test_native_run_information_rejects_forged_step_references(
    tmp_path: Path,
    field: str,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    changes = {
        "lifecycle": {"lifecycle": "failed"},
        "authority": {"authority": "diagnostic_only"},
        "artifacts": {"artifacts": step.artifacts[:-1]},
    }
    with pytest.raises(ValueError, match="exact live publication"):
        replace(step, **changes[field])
    assert not (output / "run_info").exists()


def test_native_run_information_rejects_reconstructed_live_step_reference(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    assert not hasattr(type(step), "from_live_publication")
    with pytest.raises(ValueError, match="exact live publication"):
        type(step)(
            object(),
            step.publication_occurrence_identity,
            step.path,
            step.lifecycle,
            step.authority,
            step.manifest_identity,
            step.manifest_sha256,
            step.provenance,
            step.independent_ids,
            step.artifacts,
            _witness=None,
        )


def test_native_run_information_rejects_manifest_selection_with_new_file_hash(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    manifest_path = output / "STEP" / "fit-manifest.toml"
    manifest_path.write_text(
        manifest_path.read_text().replace(
            'model = "2st"',
            'model = "FOREIGN"',
            1,
        )
    )
    with pytest.raises(ValueError, match="exact live publication"):
        replace(
            step,
            manifest_sha256=sha256(manifest_path.read_bytes()).hexdigest(),
        )


def test_native_run_information_rejects_relabel_and_foreign_workflow_reference(
    tmp_path: Path,
) -> None:
    first = _committed_fit()
    second = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", first)
    for changes in (
        {"lifecycle": "failed", "authority": "diagnostic_only"},
        {"provenance": second.provenance},
        {"artifacts": step.artifacts[:-2]},
    ):
        with pytest.raises(ValueError, match="exact live publication"):
            replace(step, **changes)


def test_native_run_information_rejects_foreign_equivalent_publication(
    tmp_path: Path,
) -> None:
    run_occurrence = _committed_fit()
    foreign_occurrence = _committed_fit()
    output = tmp_path / "Output"
    foreign_step = publish_native_results(output / "STEP", foreign_occurrence)

    with pytest.raises(ValueError, match="another run occurrence"):
        write_native_run_info(
            Namespace(output=output),
            NativeRunInformation(
                invocation_identity="issue-609-foreign-publication",
                inputs=_captured_inputs(output),
                parameter_model=run_occurrence.parameter_model,
                starting_snapshot=run_occurrence.starting_snapshot,
                steps=(foreign_step,),
            ),
            working_directory=ROOT,
        )
    assert not (output / "run_info").exists()


def test_committed_restart_provenance_binds_exact_commit_occurrence(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "STEP"
    publish_native_results(output, publication)

    record = json.loads((output / "Parameters" / "restart-provenance.json").read_text())
    assert record["starting_occurrence_identity"] == (
        publication.starting_snapshot.occurrence_identity
    )
    assert record["starting_revision"] == publication.starting_snapshot.revision
    assert record["accepted_result_identity"] == publication.accepted.identity
    assert record["accepted_occurrence_identity"] == (
        publication.accepted.occurrence_identity
    )
    assert record["commit_operation_identity"] == publication.commit_operation.identity
    assert record["commit_occurrence_identity"] == (
        publication.commit_operation.occurrence_identity
    )
    assert record["committed_occurrence_identity"] == (
        publication.committed_snapshot.occurrence_identity
    )
    assert record["committed_revision"] == publication.committed_snapshot.revision
    assert [item["id"] for item in record["independent_items"]] == list(
        publication.parameterization.independent_ids
    )


def test_durable_native_run_rejects_child_tampering_with_new_outer_hash(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    write_native_run_info(
        Namespace(output=output),
        NativeRunInformation(
            invocation_identity="issue-609-durable-tamper",
            inputs=_captured_inputs(output),
            parameter_model=publication.parameter_model,
            starting_snapshot=publication.starting_snapshot,
            steps=(step,),
        ),
        working_directory=ROOT,
    )
    run_info = output / "run_info"
    historical = read_native_run_information(run_info)
    assert historical.published_step_identities == (step.identity,)

    workflow_path = run_info / "workflows.json"
    workflows = json.loads(workflow_path.read_text())
    workflows["workflows"][0]["published_step_identity"] = "0" * 64
    workflow_path.write_text(json.dumps(workflows, sort_keys=True) + "\n")
    run_path = run_info / "run.toml"
    run_text = run_path.read_text()
    parsed = tomllib.loads(run_text)
    old_workflow_hash = parsed["workflow_records"]["sha256"]
    new_workflow_hash = sha256(workflow_path.read_bytes()).hexdigest()
    run_text = run_text.replace(old_workflow_hash, new_workflow_hash)
    parsed = tomllib.loads(run_text)
    old_record_identity = parsed["record_identity"]
    new_record_identity = run_info_module._native_run_record_identity(parsed)
    run_path.write_text(run_text.replace(old_record_identity, new_record_identity))

    with pytest.raises(ValueError, match="step lineage"):
        read_native_run_information(run_info)


@pytest.mark.parametrize(
    "relative_path",
    (
        "parameters_used.toml",
        "normalized/method-0001.toml",
        "workflows.json",
        "run.toml",
        "inputs/experiments/800mhz.toml",
    ),
)
def test_staged_run_info_mutation_prevents_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative_path: str,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    original = run_info_module.read_native_run_information

    def mutate_then_validate(path: Path) -> object:
        (path / relative_path).write_bytes(b"mutated before final publication\n")
        return original(path)

    monkeypatch.setattr(
        run_info_module,
        "read_native_run_information",
        mutate_then_validate,
    )
    with pytest.raises((ValueError, KeyError)):
        write_native_run_info(
            Namespace(output=output),
            NativeRunInformation(
                invocation_identity="issue-609-staged-run-mutation",
                inputs=_captured_inputs(output),
                parameter_model=publication.parameter_model,
                starting_snapshot=publication.starting_snapshot,
                steps=(step,),
            ),
            working_directory=ROOT,
        )
    assert not (output / "run_info").exists()
    assert not list(output.glob(".run_info-*"))


def test_native_run_information_uses_occurrence_start_input_bytes(
    tmp_path: Path,
) -> None:
    publication = _committed_fit()
    output = tmp_path / "Output"
    step = publish_native_results(output / "STEP", publication)
    first = tmp_path / "one" / "input.toml"
    second = tmp_path / "two" / "input.toml"
    target_a = tmp_path / "target-a.toml"
    target_b = tmp_path / "target-b.toml"
    link = tmp_path / "linked.toml"
    traversal_like = tmp_path / "..input.toml"
    for path, content in (
        (first, b'A = "first"\n'),
        (second, b'A = "second"\n'),
        (target_a, b'A = "symlink-a"\n'),
        (target_b, b'A = "symlink-b"\n'),
        (traversal_like, b'A = "safe-name"\n'),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    link.symlink_to(target_a)
    capture_args = Namespace(
        experiments=[first, second, link, traversal_like],
        parameters=None,
        method=None,
        output=output,
    )
    captured = capture_native_inputs(capture_args, working_directory=tmp_path)

    first.write_bytes(b'A = "mutated"\n')
    second.unlink()
    link.unlink()
    link.symlink_to(target_b)
    write_native_run_info(
        Namespace(output=output),
        NativeRunInformation(
            invocation_identity="issue-609-input-capture",
            inputs=captured,
            parameter_model=publication.parameter_model,
            starting_snapshot=publication.starting_snapshot,
            steps=(step,),
        ),
        working_directory=tmp_path,
    )

    run = tomllib.loads((output / "run_info" / "run.toml").read_text())
    archived = {
        item["provided_path"]: (output / "run_info" / item["copied_path"]).read_bytes()
        for item in run["inputs"]["experiments"]
    }
    assert archived[str(first)] == b'A = "first"\n'
    assert archived[str(second)] == b'A = "second"\n'
    assert archived[str(link)] == b'A = "symlink-a"\n'
    assert archived[str(traversal_like)] == b'A = "safe-name"\n'
    copied_names = [
        Path(item["copied_path"]).name for item in run["inputs"]["experiments"]
    ]
    assert len(copied_names) == len(set(copied_names))
    assert all(name not in {"", ".", ".."} for name in copied_names)
