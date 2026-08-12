"""Behavioral qualification tests for native step-root reporting (#608).

The seam is the published method-step directory.  Tests deliberately exercise
real native evaluation, acceptance, and commit artifacts while treating the
renderer implementation as private.
"""

from __future__ import annotations

import json
import tomllib
from dataclasses import replace
from pathlib import Path
from typing import Literal

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.optimize import native_reporting
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CommitReceipt,
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
    UncertaintyPolicy,
    compile_constraint_linearization_capabilities,
    derive_uncertainty_evidence,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


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
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    outcome = execute_direct_trf(
        problem,
        DirectTrfInvocation.for_problem(problem, objective_request_budget=80),
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
        parameterization,
        accepted,
        receipt,
        committed,
        commit_operation,
        session.analysis_values,
        uncertainty=uncertainty,
        resampling=resampling,
        mcmc=mcmc,
    )


def _evaluation_only() -> EvaluationPublication:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    parameterization = session.compile_parameterization(
        Method(fix=["R1A_A", "PB", "KEX_AB"]), experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
    return EvaluationPublication(engine.plan, parameterization, result)


def _failed_commit_operation() -> FitCommitOperation:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values()
    parameterization = session.compile_parameterization(
        read_methods([METHOD])["DEFAULT"], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    source = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, source
    )
    outcome = execute_direct_trf(
        problem,
        DirectTrfInvocation.for_problem(problem, objective_request_budget=80),
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
    return operation


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
    parameterization = session.compile_parameterization(
        read_methods([METHOD])["DEFAULT"], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    source = session.analysis_values.snapshot()
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, source
    )
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
        parameterization,
        accepted,
        receipt,
        committed,
        fabricated_operation,
        session.analysis_values,
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
        base.parameterization,
        base.accepted,
        base.commit_receipt,
        base.committed_snapshot,
        base.commit_operation,
        base.analysis_values,
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
    assert [item["ordinal"] for item in components["components"]] == [0, 1]
    assert [item["identity"] for item in components["components"]] == [
        "component-alpha",
        "component-beta",
    ]
    assert components["authority"] == "diagnostic_only"
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
    suppressed = SuppressedPublication(
        suppressed_lifecycle,
        incomplete,
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

    def fail_replace(_source: object, _destination: object) -> None:
        raise OSError("injected atomic rename failure")

    monkeypatch.setattr(native_reporting.os, "replace", fail_replace)

    with pytest.raises(OSError, match="injected atomic rename failure"):
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
        accepted_result_identity=committed.accepted.identity,
        accepted_occurrence_identity="foreign-accepted-occurrence",
        partial_mcmc=operation.evidence,
    )
    with pytest.raises(ValueError, match="genuine partial MCMC"):
        publish_native_results(tmp_path / "FOREIGN", foreign_occurrence)

    failed_commit = _failed_commit_operation()
    forged_commit = replace(failed_commit, occurrence_identity="forged-occurrence")
    forged = SuppressedPublication(
        "accepted_uncommitted",
        forged_commit,
        accepted_result_identity=forged_commit.accepted_result_identity,
        accepted_occurrence_identity=forged_commit.accepted_occurrence_identity,
    )
    with pytest.raises(ValueError, match="operation provenance disagree"):
        publish_native_results(tmp_path / "FORGED", forged)

    invalid = SuppressedPublication(
        "accepted_uncommitted",
        failed_commit,
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
        accepted_result_identity=failed_commit.accepted_result_identity,
        accepted_occurrence_identity=failed_commit.accepted_occurrence_identity,
    )

    publish_native_results(output, publication)

    assert {path.name for path in output.iterdir()} == {"Diagnostics"}
    diagnostics = json.loads((output / "Diagnostics" / "outcome.json").read_text())
    assert diagnostics["operation"]["identity"] == failed_commit.identity
    assert diagnostics["operation"]["terminal"] == "failed"
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
