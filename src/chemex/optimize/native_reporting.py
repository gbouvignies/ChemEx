"""Atomic step-root publication for native qualification results (#608).

This is an isolated qualification seam.  The production CLI remains on its
legacy reporting path until the native migration gates promote these artifacts.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal
from uuid import uuid4

from chemex.atomic import publish_directory_noreplace
from chemex.evaluation.native import EvaluationPlan, EvaluationResult
from chemex.native_provenance import (
    ArtifactReference,
    ArtifactRole,
    BudgetRecord,
    NativeProvenanceError,
    PolicyRecord,
    PublishedStepReference,
    SeedRecord,
    WorkflowProvenance,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CommitReceipt,
    DirectTrfExecution,
    DirectTrfInvocation,
    DirectTrfTerminal,
    FitCommitOperation,
    FitCommitTerminal,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
    canonical_chi_square,
    committed_values_identity,
    fit_commit_operation_is_authoritative,
)
from chemex.optimize.native_mcmc import (
    McmcEvidence,
    PosteriorSampleEvidence,
    PosteriorSummary,
)
from chemex.optimize.native_resampling import (
    ClaimState,
    ResamplingEvidence,
    ResamplingLifecycle,
    ResamplingSummaryOutcome,
    SummaryTerminal,
)
from chemex.optimize.uncertainty import UncertaintyEvidence
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.plotters.native_reporting import write_native_plots
from chemex.printers.native_reporting import (
    write_components,
    write_data,
    write_mcmc,
    write_parameter_reports,
    write_partial_mcmc,
    write_partial_resampling,
    write_resampling,
    write_restart_parameters,
    write_statistics,
    write_suppressed_outcome,
    write_uncertainty,
)


class NativePublicationError(ValueError):
    """Raised before publication when authoritative source artifacts disagree."""


@dataclass(frozen=True, slots=True)
class ComponentDiagnostic:
    """One explicitly non-authoritative fit-component execution view."""

    identity: str
    disposition: Literal[
        "succeeded",
        "failed",
        "execution_failure",
        "cancelled",
        "interrupted",
        "not_started",
    ]
    controlled_ids: tuple[str, ...]
    local_chi_square: float | None = None

    def __post_init__(self) -> None:
        if (
            not self.identity
            or not self.controlled_ids
            or len(set(self.controlled_ids)) != len(self.controlled_ids)
        ):
            raise NativePublicationError(
                "Component diagnostics require an identity and controlled scope"
            )
        if self.local_chi_square is not None and (
            not math.isfinite(self.local_chi_square) or self.local_chi_square < 0.0
        ):
            raise NativePublicationError(
                "Component diagnostic chi-square must be finite and non-negative"
            )


@dataclass(frozen=True, slots=True)
class ResamplingPublication:
    """One typed resampling family and its scope-atomic summary outcome."""

    evidence: ResamplingEvidence
    summary: ResamplingSummaryOutcome


@dataclass(frozen=True, slots=True)
class McmcPublication:
    """Primary chain topology plus explicitly derived posterior evidence."""

    evidence: McmcEvidence
    posterior_samples: PosteriorSampleEvidence
    summary: PosteriorSummary


@dataclass(frozen=True, slots=True)
class CommittedFitPublication:
    """Exact committed anchor from which ordinary fitted output may be rendered."""

    plan: EvaluationPlan
    parameterization: ActiveParameterization
    parameter_model: SealedParameterModel
    starting_snapshot: AnalysisValuesSnapshot
    primary_problem: OptimizationProblem
    primary_invocation: DirectTrfInvocation
    primary_execution: DirectTrfExecution
    accepted: AcceptedFitResult
    commit_receipt: CommitReceipt
    committed_snapshot: AnalysisValuesSnapshot
    commit_operation: FitCommitOperation
    analysis_values: AnalysisValues = field(repr=False, compare=False)
    provenance: WorkflowProvenance = field(kw_only=True)
    components: tuple[ComponentDiagnostic, ...] = ()
    uncertainty: UncertaintyEvidence | None = None
    resampling: tuple[ResamplingPublication, ...] = ()
    mcmc: McmcPublication | None = None


@dataclass(frozen=True, slots=True)
class EvaluationPublication:
    """Authoritative evaluate-only result with no estimate or commit semantics."""

    plan: EvaluationPlan
    parameterization: ActiveParameterization
    evaluation_result: EvaluationResult
    provenance: WorkflowProvenance = field(kw_only=True)


@dataclass(frozen=True, slots=True)
class SuppressedPublication:
    """Failure/uncommitted provenance with only genuine partial evidence."""

    lifecycle: Literal[
        "failed",
        "accepted_uncommitted",
        "cancelled",
        "interrupted",
    ]
    operation: (
        FitCommitOperation | DirectTrfExecution | ResamplingEvidence | McmcEvidence
    )
    plan: EvaluationPlan
    parameterization: ActiveParameterization
    provenance: WorkflowProvenance = field(kw_only=True)
    primary_invocation: DirectTrfInvocation | None = field(
        default=None,
        kw_only=True,
    )
    primary_execution: DirectTrfExecution | None = field(
        default=None,
        kw_only=True,
    )
    primary_problem: OptimizationProblem | None = field(default=None, kw_only=True)
    accepted_result_identity: str | None = None
    accepted_occurrence_identity: str | None = None
    components: tuple[ComponentDiagnostic, ...] = ()
    partial_resampling: tuple[ResamplingEvidence, ...] = ()
    partial_mcmc: McmcEvidence | None = None


type NativePublication = (
    CommittedFitPublication | EvaluationPublication | SuppressedPublication
)


def _validate_exact_provenance_records(
    provenance: WorkflowProvenance,
    policies: tuple[PolicyRecord, ...],
    budgets: tuple[BudgetRecord, ...],
    seeds: tuple[SeedRecord, ...],
) -> None:
    """Reject missing, forged, or unrelated execution records."""
    if (
        provenance.policies != policies
        or provenance.budgets != budgets
        or provenance.seeds != seeds
    ):
        raise NativePublicationError(
            "Workflow provenance differs from exact execution records"
        )


def _primary_execution_records(
    invocation: DirectTrfInvocation,
    execution: DirectTrfExecution,
) -> tuple[PolicyRecord, BudgetRecord]:
    if (
        execution.problem_identity != invocation.problem_identity
        or execution.invocation_identity != invocation.identity
    ):
        raise NativePublicationError(
            "Primary invocation and execution provenance disagree"
        )
    return (
        PolicyRecord("primary_direct_trf", invocation.identity),
        BudgetRecord(
            "primary_direct_trf_objective_requests",
            invocation.objective_request_budget,
            execution.counters.objective_requests_accepted,
        ),
    )


def _stochastic_execution_records(
    resampling: tuple[ResamplingEvidence, ...],
    mcmc: McmcEvidence | None,
) -> tuple[tuple[PolicyRecord, ...], tuple[BudgetRecord, ...], tuple[SeedRecord, ...]]:
    policies: list[PolicyRecord] = []
    budgets: list[BudgetRecord] = []
    seeds: list[SeedRecord] = []
    for evidence in resampling:
        plan = evidence.plan
        name = f"resampling_{plan.scheme.value}"
        policies.append(PolicyRecord(name, plan.identity))
        budgets.append(
            BudgetRecord(
                f"{name}_replicates",
                plan.replicate_count,
                evidence.completed_count,
            )
        )
        objective_budget = int(
            dict(plan.strategy_settings).get("objective_request_budget", "100")
        )
        budgets.append(
            BudgetRecord(
                f"{name}_objective_requests",
                objective_budget * plan.replicate_count,
                None,
            )
        )
        seeds.append(
            SeedRecord(f"{name}_root", plan.root_seed, plan.seed_policy_version)
        )
    if mcmc is None:
        return tuple(policies), tuple(budgets), tuple(seeds)
    policy = mcmc.plan.policy
    policies.append(PolicyRecord("mcmc", policy.identity))
    budgets.append(
        BudgetRecord(
            "mcmc_objective_requests",
            policy.objective_request_budget,
            mcmc.objective_request_count,
        )
    )
    seeds.append(
        SeedRecord("mcmc_root", policy.root_seed, policy.seed_policy_version),
    )
    return tuple(policies), tuple(budgets), tuple(seeds)


def _validate_provenance_context(
    provenance: WorkflowProvenance,
    plan: EvaluationPlan,
    parameterization: ActiveParameterization,
) -> None:
    try:
        provenance.validate_execution_context(parameterization, plan)
    except NativeProvenanceError as error:
        raise NativePublicationError(
            "Workflow method and selection provenance differ from execution"
        ) from error


def _validate_optional_primary_execution(
    publication: SuppressedPublication,
) -> tuple[PolicyRecord | None, BudgetRecord | None]:
    invocation = publication.primary_invocation
    execution = publication.primary_execution
    if invocation is None and execution is None:
        return None, None
    if invocation is None or execution is None:
        raise NativePublicationError(
            "Suppressed primary provenance requires invocation and execution"
        )
    return _primary_execution_records(invocation, execution)


def _validate_suppressed_execution_lineage(
    publication: SuppressedPublication,
) -> tuple[PolicyRecord | None, BudgetRecord | None]:
    _validate_provenance_context(
        publication.provenance,
        publication.plan,
        publication.parameterization,
    )
    primary_policy, primary_budget = _validate_optional_primary_execution(publication)
    primary_problem = publication.primary_problem
    if primary_policy is not None:
        if (
            primary_problem is None
            or primary_problem.parameterization_identity
            != publication.parameterization.identity
            or primary_problem.evaluator_parameterization_identity
            != publication.parameterization.evaluator_identity
            or primary_problem.evaluation_plan_identity != publication.plan.identity
            or publication.primary_invocation is None
            or publication.primary_invocation.problem_identity
            != primary_problem.identity
        ):
            raise NativePublicationError(
                "Suppressed restart scope differs from primary execution lineage"
            )
    elif primary_problem is not None:
        raise NativePublicationError(
            "Suppressed primary problem requires its invocation and execution"
        )
    if isinstance(publication.operation, (FitCommitOperation, DirectTrfExecution)) and (
        primary_problem is None
        or primary_problem.identity != publication.operation.problem_identity
    ):
        raise NativePublicationError(
            "Suppressed operation differs from primary problem lineage"
        )
    return primary_policy, primary_budget


def _validate_evaluation_artifacts(
    plan: EvaluationPlan,
    parameterization: ActiveParameterization,
    result: EvaluationResult,
) -> None:
    if (
        result.plan_identity != plan.identity
        or result.parameterization_identity != parameterization.evaluator_identity
        or plan.parameterization_identity != parameterization.evaluator_identity
        or plan.constraint_program_identity != parameterization.program.fingerprint
    ):
        raise NativePublicationError(
            "Evaluation result, parameterization, and plan are incompatible"
        )
    try:
        EvaluationPlan.from_record(plan.to_record())
        EvaluationResult.from_record(result.to_record(), plan)
    except (ArithmeticError, TypeError, ValueError) as error:
        raise NativePublicationError(
            "Native evaluation artifacts fail structural validation"
        ) from error


def _validate_evaluation_only(publication: EvaluationPublication) -> None:
    _validate_provenance_context(
        publication.provenance,
        publication.plan,
        publication.parameterization,
    )
    _validate_exact_provenance_records(publication.provenance, (), (), ())
    if any(
        publication.parameterization.role(param_id) is ParameterRole.FIT
        for param_id in publication.parameterization.independent_ids
    ):
        raise NativePublicationError(
            "Evaluate-only publication cannot contain controlled parameters"
        )
    _validate_evaluation_artifacts(
        publication.plan,
        publication.parameterization,
        publication.evaluation_result,
    )


def _validate_suppressed(publication: SuppressedPublication) -> None:
    has_partial_evidence = bool(
        publication.partial_resampling or publication.partial_mcmc is not None
    )
    if (
        publication.lifecycle == "accepted_uncommitted"
        and (
            not publication.accepted_result_identity
            or not publication.accepted_occurrence_identity
        )
    ) or (
        has_partial_evidence
        and (
            not publication.accepted_result_identity
            or not publication.accepted_occurrence_identity
        )
    ):
        raise NativePublicationError(
            "Suppressed publication requires typed lifecycle provenance"
        )
    operation = publication.operation
    primary_policy, primary_budget = _validate_suppressed_execution_lineage(publication)
    commit_failure = (
        isinstance(operation, FitCommitOperation)
        and fit_commit_operation_is_authoritative(operation)
        and operation.terminal is FitCommitTerminal.FAILED
        and operation.accepted_result_identity == publication.accepted_result_identity
        and operation.accepted_occurrence_identity
        == publication.accepted_occurrence_identity
    )
    direct_terminal = (
        operation.terminal if isinstance(operation, DirectTrfExecution) else None
    )
    stochastic_terminal = (
        operation.terminal.value
        if isinstance(operation, McmcEvidence)
        else operation.lifecycle.value
        if isinstance(operation, ResamplingEvidence)
        else None
    )
    provenance_matches = (
        (publication.lifecycle == "accepted_uncommitted" and commit_failure)
        or (
            publication.lifecycle == "failed"
            and (
                direct_terminal
                not in (None, DirectTrfTerminal.CONVERGED, DirectTrfTerminal.CANCELLED)
                or stochastic_terminal == "partial"
            )
        )
        or (
            publication.lifecycle == "cancelled"
            and (
                direct_terminal is DirectTrfTerminal.CANCELLED
                or stochastic_terminal == "cancelled"
            )
        )
        or (
            publication.lifecycle == "interrupted"
            and (
                direct_terminal is DirectTrfTerminal.INTERRUPTED
                or stochastic_terminal == "interrupted"
            )
        )
    )
    if not provenance_matches:
        raise NativePublicationError(
            "Suppressed lifecycle and atomic operation provenance disagree"
        )
    if publication.lifecycle == "accepted_uncommitted" and has_partial_evidence:
        raise NativePublicationError(
            "Accepted-but-uncommitted workflows cannot start downstream evidence"
        )
    component_ids = tuple(item.identity for item in publication.components)
    if len(set(component_ids)) != len(component_ids):
        raise NativePublicationError(
            "Suppressed component diagnostics require unique identities"
        )
    schemes: set[str] = set()
    for evidence in publication.partial_resampling:
        evidence.validate_integrity()
        scheme = evidence.plan.scheme.value
        if (
            evidence.lifecycle.value == "completed"
            or scheme in schemes
            or evidence.plan.accepted_result_identity
            != publication.accepted_result_identity
            or evidence.plan.accepted_occurrence_identity
            != publication.accepted_occurrence_identity
            or evidence.plan.parameterization is not publication.parameterization
            or evidence.plan.source_engine.plan.identity != publication.plan.identity
        ):
            raise NativePublicationError(
                "Suppressed workflows may retain only genuine partial resampling"
            )
        schemes.add(scheme)
    if publication.partial_mcmc is not None:
        publication.partial_mcmc.validate_integrity()
        if (
            publication.partial_mcmc.lifecycle.value == "completed"
            or publication.partial_mcmc.accepted_result_identity
            != publication.accepted_result_identity
            or publication.partial_mcmc.plan.accepted_occurrence_identity
            != publication.accepted_occurrence_identity
            or publication.partial_mcmc.plan.parameterization
            is not publication.parameterization
            or publication.partial_mcmc.plan.source_engine.plan.identity
            != publication.plan.identity
        ):
            raise NativePublicationError(
                "Suppressed workflows may retain only genuine partial MCMC evidence"
            )
    stochastic_policies, stochastic_budgets, stochastic_seeds = (
        _stochastic_execution_records(
            publication.partial_resampling,
            publication.partial_mcmc,
        )
    )
    _validate_exact_provenance_records(
        publication.provenance,
        (() if primary_policy is None else (primary_policy,)) + stochastic_policies,
        (() if primary_budget is None else (primary_budget,)) + stochastic_budgets,
        stochastic_seeds,
    )


def _validate_committed_parameter_state(
    publication: CommittedFitPublication,
) -> None:
    accepted = publication.accepted
    committed = publication.committed_snapshot
    parameter_model = publication.parameter_model
    if (
        parameter_model.identity
        != publication.parameterization.program.parameter_model_identity
        or parameter_model.model_identity != committed.model_identity
        or parameter_model.definitions.identity != committed.definitions_identity
        or parameter_model.configuration.identity != committed.configuration_identity
    ):
        raise NativePublicationError(
            "Committed publication parameter model differs from committed state"
        )
    starting = publication.starting_snapshot
    if (
        starting.occurrence_identity != accepted.source_occurrence_identity
        or starting.revision != accepted.source_revision
        or starting.model_identity != committed.model_identity
        or starting.definitions_identity != committed.definitions_identity
        or starting.configuration_identity != committed.configuration_identity
    ):
        raise NativePublicationError(
            "Committed publication starting state differs from accepted provenance"
        )


def _validate_committed_fit(publication: CommittedFitPublication) -> None:
    accepted = publication.accepted
    receipt = publication.commit_receipt
    committed = publication.committed_snapshot
    operation = publication.commit_operation
    current = publication.analysis_values.snapshot()
    _validate_committed_parameter_state(publication)
    _validate_provenance_context(
        publication.provenance,
        publication.plan,
        publication.parameterization,
    )
    _primary_execution_records(
        publication.primary_invocation,
        publication.primary_execution,
    )
    if not accepted_occurrence_is_authoritative(accepted):
        raise NativePublicationError(
            "Committed publication requires an authoritative accepted occurrence"
        )
    if (
        not fit_commit_operation_is_authoritative(operation)
        or operation.terminal is not FitCommitTerminal.COMMITTED
        or operation.receipt is not receipt
        or operation.committed_snapshot is not committed
        or operation.accepted_result_identity != accepted.identity
        or operation.accepted_occurrence_identity != accepted.occurrence_identity
        or operation.problem_identity != accepted.problem_identity
        or publication.primary_problem.identity != accepted.problem_identity
        or publication.primary_problem.parameterization_identity
        != publication.parameterization.identity
        or publication.primary_problem.evaluator_parameterization_identity
        != publication.parameterization.evaluator_identity
        or publication.primary_problem.evaluation_plan_identity
        != publication.plan.identity
        or publication.primary_invocation.problem_identity != accepted.problem_identity
        or publication.primary_execution.problem_identity != accepted.problem_identity
        or publication.primary_execution.invocation_identity
        != accepted.invocation_identity
        or publication.primary_execution.identity != accepted.execution_identity
        or current != committed
    ):
        raise NativePublicationError(
            "Committed publication requires the exact successful commit operation"
        )
    if (
        accepted.problem_identity != receipt.problem_identity
        or accepted.identity != receipt.accepted_result_identity
        or accepted.occurrence_identity != receipt.accepted_occurrence_identity
        or accepted.source_revision != receipt.old_revision
        or receipt.new_revision != receipt.old_revision + 1
        or accepted.commit_scope != receipt.scope
        or accepted.evaluator_parameterization_identity
        != publication.parameterization.evaluator_identity
        or committed.occurrence_identity != accepted.source_occurrence_identity
        or committed.revision != receipt.new_revision
        or committed.model_identity != receipt.model_identity
        or committed.configuration_identity != receipt.configuration_identity
        or committed.model_identity
        != publication.parameterization.program.model_identity
        or committed.definitions_identity
        != publication.parameterization.program.definitions_identity
        or committed.configuration_identity
        != publication.parameterization.program.configuration_identity
        or receipt.committed_value_identity
        != committed_values_identity(committed, receipt.scope)
    ):
        raise NativePublicationError(
            "Commit receipt does not belong to the accepted fit result"
        )
    if accepted.parameterization_identity != publication.parameterization.identity:
        raise NativePublicationError(
            "Accepted fit, parameterization, and evaluation plan are incompatible"
        )
    _validate_evaluation_artifacts(
        publication.plan,
        publication.parameterization,
        accepted.evaluation_result,
    )
    if accepted.chi_square != canonical_chi_square(
        accepted.evaluation_result.residuals
    ):
        raise NativePublicationError(
            "Accepted chi-square differs from the retained residual vector"
        )
    expected_commit_items = tuple(
        (param_id, accepted.evaluation_result.resolved_values[param_id])
        for param_id in accepted.commit_scope
    )
    committed_items = tuple(
        (param_id, committed[param_id]) for param_id in accepted.commit_scope
    )
    if (
        accepted.commit_items != expected_commit_items
        or committed_items != accepted.commit_items
        or accepted.commit_scope
        != tuple(
            param_id
            for param_id, _value in accepted.evaluation_result.resolved_values.ordered_items()
        )
    ):
        raise NativePublicationError(
            "Committed snapshot is not the accepted aggregate resolved state"
        )
    expected_controlled_ids = tuple(
        param_id
        for param_id in publication.parameterization.independent_ids
        if publication.parameterization.role(param_id) is ParameterRole.FIT
    )
    if accepted.controlled_ids != expected_controlled_ids:
        raise NativePublicationError(
            "Accepted control scope differs from the active FIT coordinates"
        )
    component_ids = tuple(item.identity for item in publication.components)
    component_controls = tuple(
        param_id for item in publication.components for param_id in item.controlled_ids
    )
    if publication.components and (
        len(set(component_ids)) != len(component_ids)
        or len(set(component_controls)) != len(component_controls)
        or set(component_controls) != set(accepted.controlled_ids)
        or any(item.disposition != "succeeded" for item in publication.components)
    ):
        raise NativePublicationError(
            "Component diagnostics must exactly partition the accepted control scope"
        )
    _validate_typed_evidence(publication)


def _validate_typed_evidence(publication: CommittedFitPublication) -> None:
    """Fail closed when any requested evidence has a foreign anchor."""
    accepted = publication.accepted
    uncertainty = publication.uncertainty
    if uncertainty is not None and (
        uncertainty.accepted_anchor is not accepted
        or uncertainty.accepted_result_identity != accepted.identity
        or uncertainty.accepted_occurrence_identity != accepted.occurrence_identity
        or uncertainty.source_parameterization is not publication.parameterization
        or uncertainty.source_engine.plan.identity != publication.plan.identity
    ):
        raise NativePublicationError(
            "Uncertainty evidence belongs to another accepted fit occurrence"
        )
    schemes: set[str] = set()
    for item in publication.resampling:
        item.evidence.validate_integrity()
        scheme = item.evidence.plan.scheme.value
        summary = item.summary.summary
        if (
            scheme in schemes
            or item.evidence.lifecycle is not ResamplingLifecycle.COMPLETED
            or item.summary.terminal is not SummaryTerminal.COMPLETED
            or summary is None
            or summary.evidence is not item.evidence
            or item.evidence.plan.accepted_result_identity != accepted.identity
            or item.evidence.plan.accepted_occurrence_identity
            != accepted.occurrence_identity
            or item.summary.evidence_identity != item.evidence.identity
        ):
            raise NativePublicationError(
                "Ordinary statistics require exact completed resampling evidence"
            )
        required_claims = (
            "STRUCTURAL_INTEGRITY",
            "INTENDED_POPULATION_TERMINAL",
            "MINIMUM_SUCCESSFUL_COVERAGE",
            "COMPLETE_SCOPE_SUCCESS_ROWS",
        )
        if any(
            item.evidence.claim(name) is not ClaimState.SATISFIED
            for name in required_claims
        ):
            raise NativePublicationError(
                "Ordinary statistics require satisfied resampling validity claims"
            )
        summary.validate_integrity()
        schemes.add(scheme)
    mcmc = publication.mcmc
    if mcmc is not None:
        mcmc.evidence.validate_integrity()
        mcmc.posterior_samples.validate_integrity()
        mcmc.summary.validate_integrity()
        if (
            mcmc.evidence.plan.accepted is not accepted
            or mcmc.evidence.accepted_result_identity != accepted.identity
            or mcmc.evidence.plan.accepted_occurrence_identity
            != accepted.occurrence_identity
            or mcmc.posterior_samples.selection.source_evidence is not mcmc.evidence
            or mcmc.summary.source is not mcmc.posterior_samples
        ):
            raise NativePublicationError(
                "MCMC evidence belongs to another accepted fit occurrence"
            )
    stochastic_policies, stochastic_budgets, stochastic_seeds = (
        _stochastic_execution_records(
            tuple(item.evidence for item in publication.resampling),
            None if mcmc is None else mcmc.evidence,
        )
    )
    primary_policy, primary_budget = _primary_execution_records(
        publication.primary_invocation,
        publication.primary_execution,
    )
    uncertainty_policies = (
        ()
        if uncertainty is None
        else (PolicyRecord("uncertainty", uncertainty.source_policy.identity),)
    )
    _validate_exact_provenance_records(
        publication.provenance,
        (primary_policy, *uncertainty_policies, *stochastic_policies),
        (primary_budget, *stochastic_budgets),
        stochastic_seeds,
    )


def _suppressed_operation_record(
    operation: FitCommitOperation
    | DirectTrfExecution
    | ResamplingEvidence
    | McmcEvidence,
) -> dict[str, object]:
    if isinstance(operation, FitCommitOperation):
        return operation.to_record()
    if isinstance(operation, McmcEvidence):
        return {"artifact_type": "native_mcmc_evidence", **operation.to_record()}
    if isinstance(operation, ResamplingEvidence):
        return {
            "artifact_type": "native_resampling_operation",
            "identity": operation.identity,
            "terminal": operation.lifecycle.value,
            "lifecycle": operation.lifecycle.value,
            "accepted_result_identity": operation.plan.accepted_result_identity,
            "accepted_occurrence_identity": (
                operation.plan.accepted_occurrence_identity
            ),
        }
    return {
        "artifact_type": "native_direct_trf_execution",
        "identity": operation.identity,
        "occurrence_identity": operation.occurrence_identity,
        "problem_identity": operation.problem_identity,
        "invocation_identity": operation.invocation_identity,
        "terminal": operation.terminal.value,
        "failure": (
            None
            if operation.failure is None
            else {
                "identity": operation.failure.identity,
                "category": operation.failure.category,
                "message": operation.failure.message,
            }
        ),
    }


def _render_suppressed(path: Path, publication: SuppressedPublication) -> None:
    path.mkdir()
    diagnostics = path / "Diagnostics"
    diagnostics.mkdir()
    write_suppressed_outcome(
        diagnostics / "outcome.json",
        lifecycle=publication.lifecycle,
        operation_record=_suppressed_operation_record(publication.operation),
        accepted_result_identity=publication.accepted_result_identity,
        accepted_occurrence_identity=publication.accepted_occurrence_identity,
        components=publication.components,
    )
    if not publication.partial_resampling and publication.partial_mcmc is None:
        return
    partial = path / "PartialEvidence"
    partial.mkdir()
    if publication.partial_resampling:
        resampling = partial / "Resampling"
        resampling.mkdir()
        for evidence in publication.partial_resampling:
            family = resampling / evidence.plan.scheme.value.upper()
            family.mkdir()
            write_partial_resampling(family / "evidence.json", evidence)
    if publication.partial_mcmc is not None:
        mcmc = partial / "MCMC"
        mcmc.mkdir()
        write_partial_mcmc(mcmc / "evidence.json", publication.partial_mcmc)


def _render_committed_fit(
    path: Path,
    publication: CommittedFitPublication,
) -> None:
    path.mkdir()
    write_parameter_reports(
        path / "Parameters",
        publication.parameterization,
        publication.accepted.evaluation_result,
    )
    write_restart_parameters(
        path / "Parameters" / "restart.toml",
        publication.parameter_model,
        publication.parameterization,
        publication.committed_snapshot,
    )
    write_data(path / "Data", publication.plan, publication.accepted.evaluation_result)
    write_native_plots(
        path / "Plots", publication.plan, publication.accepted.evaluation_result
    )
    (path / "Statistics").mkdir()
    write_statistics(
        path,
        publication.plan,
        publication.accepted.evaluation_result,
        len(publication.accepted.controlled_ids),
    )
    write_uncertainty(path / "Statistics", publication.uncertainty)
    write_resampling(
        path / "Statistics",
        tuple((item.evidence, item.summary) for item in publication.resampling),
    )
    mcmc = publication.mcmc
    write_mcmc(
        path / "Statistics",
        None if mcmc is None else mcmc.evidence,
        None if mcmc is None else mcmc.posterior_samples,
        None if mcmc is None else mcmc.summary,
    )
    write_components(path / "Components", publication.components)


def _artifact_role(relative_path: str) -> ArtifactRole:
    if relative_path == "Parameters/restart.toml":
        return "committed_restart_state"
    if relative_path == "Parameters/fitted.toml":
        return "report_only_fitted_values"
    if relative_path.startswith("Diagnostics/"):
        return "diagnostic_provenance"
    if relative_path.startswith("PartialEvidence/"):
        return "partial_evidence"
    return "product_output"


def _write_manifest(
    path: Path,
    fields: tuple[tuple[str, str | int], ...],
    provenance: WorkflowProvenance,
    *,
    components: tuple[ComponentDiagnostic, ...] = (),
    evidence: tuple[tuple[str, str, str, str], ...] = (),
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    native_threads = (
        "auto"
        if provenance.execution.native_threads is None
        else provenance.execution.native_threads
    )
    artifacts = tuple(
        ArtifactReference(
            str(item.relative_to(path)),
            _artifact_role(str(item.relative_to(path))),
            hashlib.sha256(item.read_bytes()).hexdigest(),
        )
        for item in sorted(path.rglob("*"))
        if item.is_file()
    )
    manifest_identity = hashlib.sha256(
        json.dumps(
            {
                "schema_version": 1,
                "fields": fields,
                "provenance_identity": provenance.identity,
                "artifacts": [
                    (item.path, item.role, item.sha256) for item in artifacts
                ],
                "components": [
                    (item.identity, item.disposition, item.controlled_ids)
                    for item in components
                ],
                "evidence": evidence,
            },
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode()
    ).hexdigest()
    lines = [
        "# Atomic native method-step manifest.",
        "schema_version = 1",
        f"manifest_identity = {json.dumps(manifest_identity)}",
        *(f"{name} = {json.dumps(value)}" for name, value in fields),
        "",
        "[workflow]",
        f"identity = {json.dumps(provenance.workflow_identity)}",
        f"provenance_identity = {json.dumps(provenance.identity)}",
        "",
        "[method]",
        f"identity = {json.dumps(provenance.method_identity)}",
        f"parameterization_identity = {json.dumps(provenance.parameterization_identity)}",
        f"evaluation_plan_identity = {json.dumps(provenance.evaluation_plan_identity)}",
        f"normalized = {json.dumps(provenance.normalized_method_text)}",
        f"normalized_sha256 = {json.dumps(provenance.normalized_method_sha256)}",
        "",
        "[selection]",
        f"model = {json.dumps(provenance.selection.model_name)}",
        f"include = {json.dumps(list(provenance.selection.include))}",
        f"exclude = {json.dumps(list(provenance.selection.exclude))}",
        "",
        "[execution]",
        f"workers = {provenance.execution.workers}",
        f"native_threads = {json.dumps(native_threads)}"
        if isinstance(native_threads, str)
        else f"native_threads = {native_threads}",
        "",
        "[environment]",
        *(
            f"{name} = {json.dumps(value)}"
            for name, value in provenance.environment.to_record().items()
            if name != "numerical_libraries"
        ),
        "",
    ]
    for kind, name, library_version in provenance.environment.numerical_libraries:
        lines.extend(
            (
                "[[environment.numerical_libraries]]",
                f"kind = {json.dumps(kind)}",
                f"name = {json.dumps(name)}",
                f"version = {json.dumps(library_version)}",
                "",
            )
        )
    for policy in provenance.policies:
        lines.extend(
            (
                "[[policies]]",
                f"name = {json.dumps(policy.name)}",
                f"identity = {json.dumps(policy.identity)}",
                "",
            )
        )
    for budget in provenance.budgets:
        lines.extend(
            (
                "[[budgets]]",
                f"name = {json.dumps(budget.name)}",
                f"limit = {budget.limit}",
            )
        )
        if budget.used is not None:
            lines.append(f"used = {budget.used}")
        lines.append("")
    for seed in provenance.seeds:
        lines.extend(
            (
                "[[seeds]]",
                f"name = {json.dumps(seed.name)}",
                f"value = {seed.value}",
                f"policy_identity = {json.dumps(seed.policy_identity)}",
                "",
            )
        )
    for reference in provenance.baseline_references:
        lines.extend(
            (
                "[[baseline_references]]",
                f"kind = {json.dumps(reference.kind)}",
                f"identity = {json.dumps(reference.identity)}",
                "",
            )
        )
    for component in components:
        lines.extend(
            (
                "[[components]]",
                f"identity = {json.dumps(component.identity)}",
                f"disposition = {json.dumps(component.disposition)}",
                f"controlled_ids = {json.dumps(list(component.controlled_ids))}",
                'location = "Components/index.json"',
                'authority = "diagnostic_only"',
                "",
            )
        )
    for family, identity, validity, location in evidence:
        lines.extend(
            (
                "[[evidence]]",
                f"family = {json.dumps(family)}",
                f"identity = {json.dumps(identity)}",
                f"validity = {json.dumps(validity)}",
                f"location = {json.dumps(location)}",
                "",
            )
        )
    for artifact in artifacts:
        lines.extend(
            (
                f"[artifacts.{json.dumps(artifact.path)}]",
                f"role = {json.dumps(artifact.role)}",
                f"sha256 = {json.dumps(artifact.sha256)}",
                "",
            )
        )
    manifest_path = path / "fit-manifest.toml"
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    return (
        manifest_identity,
        hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        artifacts,
    )


def _write_committed_manifest(
    path: Path,
    publication: CommittedFitPublication,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    evidence: list[tuple[str, str, str, str]] = []
    if publication.uncertainty is not None:
        evidence.append(
            (
                "covariance",
                publication.uncertainty.identity,
                "typed_claims_in_artifact",
                "Statistics/Covariance/evidence.json",
            )
        )
        if publication.uncertainty.requested_output_scope:
            evidence.append(
                (
                    "constrained",
                    publication.uncertainty.identity,
                    "typed_claims_in_artifact",
                    "Statistics/Constrained/evidence.json",
                )
            )
    for item in publication.resampling:
        scheme = item.evidence.plan.scheme.value
        evidence.append(
            (
                f"resampling:{scheme}",
                item.evidence.identity,
                "required_claims_satisfied",
                f"Statistics/Resampling/{scheme.upper()}/evidence.json",
            )
        )
    if publication.mcmc is not None:
        evidence.append(
            (
                "mcmc",
                publication.mcmc.evidence.identity,
                "integrity_validated",
                "Statistics/MCMC/evidence.json",
            )
        )
    return _write_manifest(
        path,
        (
            ("lifecycle", "committed"),
            ("authority", "committed_fit"),
            ("starting_state_kind", "workflow_starting_provenance"),
            (
                "starting_occurrence_identity",
                publication.starting_snapshot.occurrence_identity,
            ),
            ("starting_revision", publication.starting_snapshot.revision),
            ("accepted_result_identity", publication.accepted.identity),
            (
                "accepted_occurrence_identity",
                publication.accepted.occurrence_identity,
            ),
            ("commit_receipt_identity", publication.commit_receipt.identity),
        ),
        publication.provenance,
        components=publication.components,
        evidence=tuple(evidence),
    )


def _write_evaluation_manifest(
    path: Path,
    publication: EvaluationPublication,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    return _write_manifest(
        path,
        (
            ("lifecycle", "successful_no_state_change"),
            ("authority", "evaluation_only"),
        ),
        publication.provenance,
    )


def _write_suppressed_manifest(
    path: Path,
    publication: SuppressedPublication,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    fields: list[tuple[str, str | int]] = [
        ("lifecycle", publication.lifecycle),
        ("authority", "diagnostic_only"),
        ("operation_identity", publication.operation.identity),
    ]
    if publication.accepted_result_identity is not None:
        fields.append(
            ("accepted_result_identity", publication.accepted_result_identity)
        )
    if publication.accepted_occurrence_identity is not None:
        fields.append(
            (
                "accepted_occurrence_identity",
                publication.accepted_occurrence_identity,
            )
        )
    evidence = tuple(
        (
            f"partial_resampling:{item.plan.scheme.value}",
            item.identity,
            item.lifecycle.value,
            f"PartialEvidence/Resampling/{item.plan.scheme.value.upper()}/evidence.json",
        )
        for item in publication.partial_resampling
    )
    if publication.partial_mcmc is not None:
        evidence += (
            (
                "partial_mcmc",
                publication.partial_mcmc.identity,
                publication.partial_mcmc.lifecycle.value,
                "PartialEvidence/MCMC/evidence.json",
            ),
        )
    return _write_manifest(
        path,
        tuple(fields),
        publication.provenance,
        components=publication.components,
        evidence=evidence,
    )


def _render_evaluation(path: Path, publication: EvaluationPublication) -> None:
    path.mkdir()
    write_parameter_reports(
        path / "Parameters",
        publication.parameterization,
        publication.evaluation_result,
    )
    write_data(path / "Data", publication.plan, publication.evaluation_result)
    write_native_plots(path / "Plots", publication.plan, publication.evaluation_result)
    (path / "Statistics").mkdir()
    write_statistics(path, publication.plan, publication.evaluation_result, 0)


def publish_native_results(
    path: Path,
    publication: NativePublication,
) -> PublishedStepReference:
    """Validate and atomically publish one native method-step result tree."""
    if isinstance(publication, CommittedFitPublication):
        _validate_committed_fit(publication)
    elif isinstance(publication, EvaluationPublication):
        _validate_evaluation_only(publication)
    else:
        _validate_suppressed(publication)
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(f"Native publication destination exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.parent / f".{destination.name}.tmp-{uuid4().hex}"
    try:
        if isinstance(publication, CommittedFitPublication):
            _render_committed_fit(staging, publication)
            manifest = _write_committed_manifest(staging, publication)
            lifecycle = "committed"
            authority = "committed_fit"
            independent_ids = publication.parameterization.independent_ids
        elif isinstance(publication, EvaluationPublication):
            _render_evaluation(staging, publication)
            manifest = _write_evaluation_manifest(staging, publication)
            lifecycle = "successful_no_state_change"
            authority = "evaluation_only"
            independent_ids = publication.parameterization.independent_ids
        else:
            _render_suppressed(staging, publication)
            manifest = _write_suppressed_manifest(staging, publication)
            lifecycle = publication.lifecycle
            authority = "diagnostic_only"
            independent_ids = publication.parameterization.independent_ids
        publish_directory_noreplace(staging, destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    manifest_identity, manifest_sha256, artifacts = manifest
    return PublishedStepReference(
        destination,
        lifecycle,
        authority,
        manifest_identity,
        manifest_sha256,
        publication.provenance,
        independent_ids,
        artifacts,
    )
