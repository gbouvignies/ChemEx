"""Atomic step-root publication for native qualification results (#608).

This is an isolated qualification seam.  The production CLI remains on its
legacy reporting path until the native migration gates promote these artifacts.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal
from uuid import uuid4

from chemex.atomic import publish_directory_noreplace
from chemex.configuration.methods import Method
from chemex.evaluation.native import EvaluationFailure, EvaluationPlan, EvaluationResult
from chemex.native_provenance import (
    ArtifactReference,
    ArtifactRole,
    BudgetRecord,
    CommittedRestartRecord,
    NativeProvenanceError,
    PolicyRecord,
    PublishedStepReference,
    SeedRecord,
    WorkflowProvenance,
    _published_step_from_successful_native_publication,
    native_step_manifest_identity,
    validate_native_step_manifest_bytes,
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
    SummaryFailure,
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
    write_evaluated_parameters,
    write_mcmc,
    write_parameter_reports,
    write_partial_mcmc,
    write_partial_resampling,
    write_resampling,
    write_resampling_summary_failure,
    write_restart_parameters,
    write_statistics,
    write_suppressed_outcome,
    write_uncertainty,
)
from chemex.runtime import ExecutionSettings


class NativePublicationError(ValueError):
    """Raised before publication when authoritative source artifacts disagree."""


type MethodStepPrimaryTerminal = Literal[
    "accepted",
    "component_failure",
    "execution_failure",
    "decomposition_validation_failure",
    "no_eligible_candidate",
    "materialization_failure",
    "search_unsuccessful",
    "polish_unsuccessful",
    "cancelled",
    "interrupted",
]


@dataclass(frozen=True, slots=True)
class MethodStepPrimaryRecord:
    """Exact aggregate-primary view used by the existing reporting seam."""

    semantic_workflow_identity: str
    problem_identity: str
    invocation_identity: str
    execution_identity: str
    aggregate_execution_identity: str
    terminal: MethodStepPrimaryTerminal
    grouping_topology: tuple[tuple[str, tuple[str, ...]], ...]
    execution_settings: ExecutionSettings
    policy: PolicyRecord
    budget: BudgetRecord
    seeds: tuple[SeedRecord, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            not self.semantic_workflow_identity
            or not self.problem_identity
            or not self.invocation_identity
            or not self.execution_identity
            or not self.aggregate_execution_identity
            or self.terminal
            not in {
                "accepted",
                "component_failure",
                "execution_failure",
                "decomposition_validation_failure",
                "no_eligible_candidate",
                "materialization_failure",
                "search_unsuccessful",
                "polish_unsuccessful",
                "cancelled",
                "interrupted",
            }
            or not self.grouping_topology
        ):
            raise NativePublicationError(
                "Method-step primary reporting lineage is incomplete"
            )
        object.__setattr__(
            self,
            "identity",
            hashlib.sha256(
                json.dumps(
                    (
                        "native-method-step-primary-reporting-v1",
                        self.semantic_workflow_identity,
                        self.problem_identity,
                        self.invocation_identity,
                        self.execution_identity,
                        self.aggregate_execution_identity,
                        self.terminal,
                        self.grouping_topology,
                        (
                            self.execution_settings.workers,
                            self.execution_settings.native_threads,
                        ),
                        (self.policy.name, self.policy.identity),
                        (self.budget.name, self.budget.limit, self.budget.used),
                        tuple(
                            (item.name, item.value, item.policy_identity)
                            for item in self.seeds
                        ),
                    ),
                    ensure_ascii=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        )


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
class ResamplingSummaryFailurePublication:
    """Completed primary resampling evidence whose derived summary failed."""

    evidence: ResamplingEvidence
    failure: SummaryFailure


@dataclass(frozen=True, slots=True)
class McmcPublication:
    """Primary chain topology plus explicitly derived posterior evidence."""

    evidence: McmcEvidence
    posterior_samples: PosteriorSampleEvidence
    summary: PosteriorSummary


@dataclass(frozen=True, slots=True, weakref_slot=True)
class CommittedFitPublication:
    """Exact committed anchor from which ordinary fitted output may be rendered."""

    plan: EvaluationPlan
    method: Method
    parameterization: ActiveParameterization
    parameter_model: SealedParameterModel
    starting_snapshot: AnalysisValuesSnapshot
    primary_problem: OptimizationProblem
    primary_invocation: DirectTrfInvocation | MethodStepPrimaryRecord
    primary_execution: DirectTrfExecution | MethodStepPrimaryRecord
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
    partial_resampling: tuple[ResamplingEvidence, ...] = ()
    resampling_summary_failures: tuple[ResamplingSummaryFailurePublication, ...] = ()
    partial_mcmc: McmcEvidence | None = None


@dataclass(frozen=True, slots=True, weakref_slot=True)
class EvaluationPublication:
    """Authoritative evaluate-only result with no estimate or commit semantics."""

    plan: EvaluationPlan
    method: Method
    parameterization: ActiveParameterization
    evaluation_result: EvaluationResult
    provenance: WorkflowProvenance = field(kw_only=True)
    method_step_semantic_identity: str | None = field(default=None, kw_only=True)
    allow_controlled: bool = field(default=False, kw_only=True)


@dataclass(frozen=True, slots=True, weakref_slot=True)
class NoObjectivePublication:
    """Typed no-objective occurrence with diagnostic-only method-step output."""

    plan: EvaluationPlan
    method: Method
    parameterization: ActiveParameterization
    semantic_workflow_identity: str
    provenance: WorkflowProvenance = field(kw_only=True)


@dataclass(frozen=True, slots=True, weakref_slot=True)
class SuppressedPublication:
    """Failure/uncommitted provenance with only genuine partial evidence."""

    lifecycle: Literal[
        "failed",
        "accepted_uncommitted",
        "cancelled",
        "interrupted",
    ]
    operation: (
        FitCommitOperation
        | DirectTrfExecution
        | MethodStepPrimaryRecord
        | EvaluationFailure
        | ResamplingEvidence
        | McmcEvidence
    )
    plan: EvaluationPlan
    method: Method
    parameterization: ActiveParameterization
    provenance: WorkflowProvenance = field(kw_only=True)
    primary_invocation: DirectTrfInvocation | MethodStepPrimaryRecord | None = field(
        default=None,
        kw_only=True,
    )
    primary_execution: DirectTrfExecution | MethodStepPrimaryRecord | None = field(
        default=None,
        kw_only=True,
    )
    primary_problem: OptimizationProblem | None = field(default=None, kw_only=True)
    method_step_semantic_identity: str | None = field(default=None, kw_only=True)
    accepted_result_identity: str | None = None
    accepted_occurrence_identity: str | None = None
    components: tuple[ComponentDiagnostic, ...] = ()
    partial_resampling: tuple[ResamplingEvidence, ...] = ()
    partial_mcmc: McmcEvidence | None = None


type NativePublication = (
    CommittedFitPublication
    | EvaluationPublication
    | NoObjectivePublication
    | SuppressedPublication
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
    invocation: DirectTrfInvocation | MethodStepPrimaryRecord,
    execution: DirectTrfExecution | MethodStepPrimaryRecord,
) -> tuple[PolicyRecord, BudgetRecord]:
    if isinstance(invocation, MethodStepPrimaryRecord):
        if execution is not invocation:
            raise NativePublicationError(
                "Method-step reporting requires one exact aggregate primary record"
            )
        return invocation.policy, invocation.budget
    if not isinstance(execution, DirectTrfExecution):
        raise NativePublicationError("Direct TRF reporting execution is incompatible")
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


def _primary_execution_seeds(
    invocation: DirectTrfInvocation | MethodStepPrimaryRecord,
) -> tuple[SeedRecord, ...]:
    """Return aggregate-primary seeds, if the selected strategy owns any."""
    return invocation.seeds if isinstance(invocation, MethodStepPrimaryRecord) else ()


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


def publication_provenance_records(
    primary: MethodStepPrimaryRecord,
    *,
    uncertainty: UncertaintyEvidence | None = None,
    resampling: tuple[ResamplingEvidence, ...] = (),
    mcmc: McmcEvidence | None = None,
) -> tuple[tuple[PolicyRecord, ...], tuple[BudgetRecord, ...], tuple[SeedRecord, ...]]:
    """Return the exact records expected by native committed publication."""
    stochastic_policies, stochastic_budgets, stochastic_seeds = (
        _stochastic_execution_records(resampling, mcmc)
    )
    uncertainty_policies = (
        ()
        if uncertainty is None
        else (PolicyRecord("uncertainty", uncertainty.source_policy.identity),)
    )
    return (
        (primary.policy, *uncertainty_policies, *stochastic_policies),
        (primary.budget, *stochastic_budgets),
        (*primary.seeds, *stochastic_seeds),
    )


def _validate_provenance_context(
    provenance: WorkflowProvenance,
    plan: EvaluationPlan,
    parameterization: ActiveParameterization,
    method: Method,
    invocation: DirectTrfInvocation | MethodStepPrimaryRecord | None,
    native_execution: DirectTrfExecution | MethodStepPrimaryRecord | EvaluationResult,
) -> None:
    try:
        if isinstance(invocation, MethodStepPrimaryRecord):
            if native_execution is not invocation:
                raise NativeProvenanceError(
                    "Method-step publication uses another aggregate primary"
                )
            provenance.validate_method_step_context(
                parameterization=parameterization,
                plan=plan,
                method=method,
                semantic_workflow_identity=invocation.semantic_workflow_identity,
                grouping_topology=invocation.grouping_topology,
                execution=invocation.execution_settings,
            )
        else:
            if isinstance(native_execution, MethodStepPrimaryRecord):
                raise NativeProvenanceError(
                    "Direct publication cannot use method-step primary evidence"
                )
            provenance.validate_execution_context(
                parameterization,
                plan,
                method,
                invocation,
                native_execution,
            )
    except NativeProvenanceError as error:
        raise NativePublicationError(
            "Workflow method and selection provenance, or execution provenance, "
            "differ from execution"
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
    if isinstance(publication.operation, EvaluationFailure):
        failure = publication.operation
        if (
            publication.primary_invocation is not None
            or publication.primary_execution is not None
            or publication.primary_problem is not None
            or publication.method_step_semantic_identity is None
            or failure.plan_identity != publication.plan.identity
            or failure.parameterization_identity
            != publication.parameterization.evaluator_identity
        ):
            raise NativePublicationError(
                "Suppressed evaluation failure has incompatible lineage"
            )
        EvaluationFailure.from_record(failure.to_record(), publication.plan)
        try:
            publication.provenance.validate_method_step_context(
                parameterization=publication.parameterization,
                plan=publication.plan,
                method=publication.method,
                semantic_workflow_identity=publication.method_step_semantic_identity,
                grouping_topology=(
                    ("aggregate", publication.parameterization.independent_ids),
                ),
                execution=ExecutionSettings(),
            )
        except NativeProvenanceError as error:
            raise NativePublicationError(
                "Suppressed evaluation provenance differs from its method step"
            ) from error
        return None, None
    if publication.primary_execution is None:
        raise NativePublicationError(
            "Suppressed provenance requires the exact primary execution"
        )
    _validate_provenance_context(
        publication.provenance,
        publication.plan,
        publication.parameterization,
        publication.method,
        publication.primary_invocation,
        publication.primary_execution,
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
    if isinstance(
        publication.operation,
        (FitCommitOperation, DirectTrfExecution, MethodStepPrimaryRecord),
    ) and (
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
    if publication.method_step_semantic_identity is None:
        _validate_provenance_context(
            publication.provenance,
            publication.plan,
            publication.parameterization,
            publication.method,
            None,
            publication.evaluation_result,
        )
    else:
        try:
            publication.provenance.validate_method_step_context(
                parameterization=publication.parameterization,
                plan=publication.plan,
                method=publication.method,
                semantic_workflow_identity=publication.method_step_semantic_identity,
                grouping_topology=(
                    ("aggregate", publication.parameterization.independent_ids),
                ),
                execution=ExecutionSettings(),
            )
        except NativeProvenanceError as error:
            raise NativePublicationError(
                "Evaluation publication differs from its method-step context"
            ) from error
    _validate_exact_provenance_records(publication.provenance, (), (), ())
    if not publication.allow_controlled and any(
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


def _validate_no_objective(publication: NoObjectivePublication) -> None:
    try:
        publication.provenance.validate_method_step_context(
            parameterization=publication.parameterization,
            plan=publication.plan,
            method=publication.method,
            semantic_workflow_identity=publication.semantic_workflow_identity,
            grouping_topology=(
                ("aggregate", publication.parameterization.independent_ids),
            ),
            execution=ExecutionSettings(),
        )
    except NativeProvenanceError as error:
        raise NativePublicationError(
            "No-objective publication differs from its method-step context"
        ) from error
    _validate_exact_provenance_records(publication.provenance, (), (), ())
    if publication.plan.retained_observation_count:
        raise NativePublicationError("No-objective publication retained objective data")


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
    method_step_terminal = (
        operation.terminal if isinstance(operation, MethodStepPrimaryRecord) else None
    )
    evaluation_failed = isinstance(operation, EvaluationFailure)
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
                or method_step_terminal
                not in (None, "accepted", "cancelled", "interrupted")
                or evaluation_failed
                or stochastic_terminal == "partial"
            )
        )
        or (
            publication.lifecycle == "cancelled"
            and (
                direct_terminal is DirectTrfTerminal.CANCELLED
                or method_step_terminal == "cancelled"
                or stochastic_terminal == "cancelled"
            )
        )
        or (
            publication.lifecycle == "interrupted"
            and (
                direct_terminal is DirectTrfTerminal.INTERRUPTED
                or method_step_terminal == "interrupted"
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
        (
            ()
            if publication.primary_invocation is None
            else _primary_execution_seeds(publication.primary_invocation)
        )
        + stochastic_seeds,
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
        publication.method,
        publication.primary_invocation,
        publication.primary_execution,
    )
    _primary_execution_records(
        publication.primary_invocation,
        publication.primary_execution,
    )
    if not accepted_occurrence_is_authoritative(accepted):
        raise NativePublicationError(
            "Committed publication requires an authoritative accepted occurrence"
        )
    primary_invocation = publication.primary_invocation
    primary_execution = publication.primary_execution
    if isinstance(primary_invocation, MethodStepPrimaryRecord):
        primary_matches = (
            primary_execution is primary_invocation
            and primary_invocation.terminal == "accepted"
            and primary_invocation.problem_identity == accepted.problem_identity
            and primary_invocation.invocation_identity == accepted.invocation_identity
            and primary_invocation.execution_identity == accepted.execution_identity
        )
    else:
        primary_matches = (
            isinstance(primary_execution, DirectTrfExecution)
            and primary_invocation.problem_identity == accepted.problem_identity
            and primary_execution.problem_identity == accepted.problem_identity
            and primary_execution.invocation_identity == accepted.invocation_identity
            and primary_execution.identity == accepted.execution_identity
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
        or not primary_matches
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


def _validate_typed_evidence(  # noqa: C901 - closed evidence families
    publication: CommittedFitPublication,
) -> None:
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
    ordinary_schemes = {
        item.evidence.plan.scheme.value for item in publication.resampling
    }
    for evidence in publication.partial_resampling:
        evidence.validate_integrity()
        if (
            evidence.plan.scheme.value in ordinary_schemes
            or evidence.plan.accepted_result_identity != accepted.identity
            or evidence.plan.accepted_occurrence_identity
            != accepted.occurrence_identity
            or evidence.plan.parameterization is not publication.parameterization
            or evidence.plan.source_engine.plan.identity != publication.plan.identity
        ):
            raise NativePublicationError(
                "Partial resampling evidence belongs to another committed fit"
            )
    partial_schemes = {
        evidence.plan.scheme.value for evidence in publication.partial_resampling
    }
    failed_summary_schemes: set[str] = set()
    for item in publication.resampling_summary_failures:
        item.evidence.validate_integrity()
        scheme = item.evidence.plan.scheme.value
        if (
            item.evidence.lifecycle is not ResamplingLifecycle.COMPLETED
            or item.failure.source_evidence_identity != item.evidence.identity
            or scheme in ordinary_schemes
            or scheme in partial_schemes
            or scheme in failed_summary_schemes
            or item.evidence.plan.accepted_result_identity != accepted.identity
            or item.evidence.plan.accepted_occurrence_identity
            != accepted.occurrence_identity
            or item.evidence.plan.parameterization is not publication.parameterization
            or item.evidence.plan.source_engine.plan.identity
            != publication.plan.identity
        ):
            raise NativePublicationError(
                "Resampling summary failure belongs to another committed fit"
            )
        failed_summary_schemes.add(scheme)
    if publication.partial_mcmc is not None:
        publication.partial_mcmc.validate_integrity()
        if (
            mcmc is not None
            or publication.partial_mcmc.accepted_result_identity != accepted.identity
            or publication.partial_mcmc.plan.accepted_occurrence_identity
            != accepted.occurrence_identity
            or publication.partial_mcmc.plan.parameterization
            is not publication.parameterization
            or publication.partial_mcmc.plan.source_engine.plan.identity
            != publication.plan.identity
        ):
            raise NativePublicationError(
                "Partial MCMC evidence belongs to another committed fit"
            )
    stochastic_policies, stochastic_budgets, stochastic_seeds = (
        _stochastic_execution_records(
            (
                *(item.evidence for item in publication.resampling),
                *publication.partial_resampling,
                *(item.evidence for item in publication.resampling_summary_failures),
            ),
            (publication.partial_mcmc if mcmc is None else mcmc.evidence),
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
        (*_primary_execution_seeds(publication.primary_invocation), *stochastic_seeds),
    )


def _suppressed_operation_record(
    operation: FitCommitOperation
    | DirectTrfExecution
    | MethodStepPrimaryRecord
    | EvaluationFailure
    | ResamplingEvidence
    | McmcEvidence,
) -> dict[str, object]:
    if isinstance(operation, FitCommitOperation):
        return operation.to_record()
    if isinstance(operation, EvaluationFailure):
        return {
            "artifact_type": "native_evaluation_failure",
            **operation.to_record(),
        }
    if isinstance(operation, MethodStepPrimaryRecord):
        return {
            "artifact_type": "native_method_step_primary",
            "identity": operation.identity,
            "problem_identity": operation.problem_identity,
            "invocation_identity": operation.invocation_identity,
            "execution_identity": operation.execution_identity,
            "aggregate_execution_identity": operation.aggregate_execution_identity,
            "terminal": operation.terminal,
        }
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
    write_components(path / "Components", publication.components)
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
) -> CommittedRestartRecord:
    path.mkdir()
    write_parameter_reports(
        path / "Parameters",
        publication.parameterization,
        publication.accepted.evaluation_result,
        publication.uncertainty,
    )
    write_restart_parameters(
        path / "Parameters" / "restart.toml",
        publication.parameter_model,
        publication.parameterization,
        publication.committed_snapshot,
    )
    restart_path = path / "Parameters" / "restart.toml"
    restart = CommittedRestartRecord(
        publication.starting_snapshot.occurrence_identity,
        publication.starting_snapshot.revision,
        publication.accepted.identity,
        publication.accepted.occurrence_identity,
        publication.commit_operation.identity,
        publication.commit_operation.occurrence_identity,
        publication.committed_snapshot.occurrence_identity,
        publication.committed_snapshot.revision,
        publication.provenance.workflow_identity,
        publication.primary_problem.identity,
        publication.parameterization.identity,
        publication.parameter_model.identity,
        publication.committed_snapshot.configuration_identity,
        tuple(
            (
                param_id,
                float(publication.committed_snapshot[param_id]).hex(),
                float(
                    publication.parameter_model.configuration[param_id].lower_bound
                ).hex(),
                float(
                    publication.parameter_model.configuration[param_id].upper_bound
                ).hex(),
            )
            for param_id in publication.parameterization.independent_ids
        ),
        hashlib.sha256(restart_path.read_bytes()).hexdigest(),
    )
    (path / "Parameters" / "restart-provenance.json").write_text(
        json.dumps(
            restart.to_record(),
            allow_nan=False,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
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
    if (
        publication.partial_resampling
        or publication.resampling_summary_failures
        or publication.partial_mcmc is not None
    ):
        partial = path / "PartialEvidence"
        partial.mkdir()
        if publication.partial_resampling or publication.resampling_summary_failures:
            resampling_path = partial / "Resampling"
            resampling_path.mkdir()
            for evidence in publication.partial_resampling:
                family = resampling_path / evidence.plan.scheme.value.upper()
                family.mkdir()
                write_partial_resampling(family / "evidence.json", evidence)
            for item in publication.resampling_summary_failures:
                family = resampling_path / item.evidence.plan.scheme.value.upper()
                family.mkdir()
                write_partial_resampling(family / "evidence.json", item.evidence)
                write_resampling_summary_failure(
                    family / "summary-failure.json",
                    item.failure,
                )
        if publication.partial_mcmc is not None:
            mcmc_path = partial / "MCMC"
            mcmc_path.mkdir()
            write_partial_mcmc(mcmc_path / "evidence.json", publication.partial_mcmc)
    write_components(path / "Components", publication.components)
    return restart


def _artifact_role(relative_path: str) -> ArtifactRole:
    if relative_path == "Parameters/restart.toml":
        return "committed_restart_state"
    if relative_path == "Parameters/restart-provenance.json":
        return "committed_restart_state"
    if relative_path == "Parameters/fitted.toml":
        return "report_only_fitted_values"
    if relative_path.startswith(("Diagnostics/", "Components/")):
        return "diagnostic_provenance"
    if relative_path.startswith("PartialEvidence/"):
        return "partial_evidence"
    return "product_output"


def _write_manifest(  # noqa: C901 - deterministic TOML sections
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
    lines = [
        "# Atomic native method-step manifest.",
        "schema_version = 1",
        'manifest_identity = ""',
        *(f"{name} = {json.dumps(value)}" for name, value in fields),
        "",
        "[workflow]",
        f"identity = {json.dumps(provenance.workflow_identity)}",
        f"execution_identity = {json.dumps(provenance.execution_identity)}",
        f"provenance_identity = {json.dumps(provenance.identity)}",
        f"type = {json.dumps(provenance.workflow_type)}",
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
    for group_identity, controlled_ids in provenance.grouping_topology:
        lines.extend(
            (
                "[[workflow.groups]]",
                f"identity = {json.dumps(group_identity)}",
                f"controlled_ids = {json.dumps(list(controlled_ids))}",
                "",
            )
        )
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
                f"occurrence_identity = {json.dumps(reference.occurrence_identity)}",
                f"result_bundle_identity = {json.dumps(reference.result_bundle_identity)}",
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
    provisional = "\n".join(lines)
    manifest_identity = native_step_manifest_identity(tomllib.loads(provisional))
    lines[2] = f"manifest_identity = {json.dumps(manifest_identity)}"
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    validate_native_step_manifest_bytes(manifest_path.read_bytes())
    return (
        manifest_identity,
        hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        artifacts,
    )


def _write_committed_manifest(
    path: Path,
    publication: CommittedFitPublication,
    publication_occurrence_identity: str,
    restart: CommittedRestartRecord,
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
    for item in publication.partial_resampling:
        scheme = item.plan.scheme.value
        evidence.append(
            (
                f"partial_resampling:{scheme}",
                item.identity,
                item.lifecycle.value,
                f"PartialEvidence/Resampling/{scheme.upper()}/evidence.json",
            )
        )
    for item in publication.resampling_summary_failures:
        scheme = item.evidence.plan.scheme.value
        evidence.extend(
            (
                (
                    f"partial_resampling:{scheme}",
                    item.evidence.identity,
                    "summary_failed",
                    f"PartialEvidence/Resampling/{scheme.upper()}/evidence.json",
                ),
                (
                    f"resampling_summary_failure:{scheme}",
                    item.failure.identity,
                    item.failure.category,
                    f"PartialEvidence/Resampling/{scheme.upper()}/summary-failure.json",
                ),
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
    if publication.partial_mcmc is not None:
        evidence.append(
            (
                "partial_mcmc",
                publication.partial_mcmc.identity,
                publication.partial_mcmc.lifecycle.value,
                "PartialEvidence/MCMC/evidence.json",
            )
        )
    return _write_manifest(
        path,
        (
            ("lifecycle", "committed"),
            ("authority", "committed_fit"),
            ("publication_occurrence_identity", publication_occurrence_identity),
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
            ("commit_operation_identity", publication.commit_operation.identity),
            (
                "commit_occurrence_identity",
                publication.commit_operation.occurrence_identity,
            ),
            (
                "committed_occurrence_identity",
                publication.committed_snapshot.occurrence_identity,
            ),
            ("committed_revision", publication.committed_snapshot.revision),
            ("problem_identity", publication.primary_problem.identity),
            ("parameterization_identity", publication.parameterization.identity),
            ("restart_record_identity", restart.identity),
        ),
        publication.provenance,
        components=publication.components,
        evidence=tuple(evidence),
    )


def _write_evaluation_manifest(
    path: Path,
    publication: EvaluationPublication,
    publication_occurrence_identity: str,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    return _write_manifest(
        path,
        (
            ("lifecycle", "successful_no_state_change"),
            ("authority", "evaluation_only"),
            ("publication_occurrence_identity", publication_occurrence_identity),
        ),
        publication.provenance,
    )


def _write_suppressed_manifest(
    path: Path,
    publication: SuppressedPublication,
    publication_occurrence_identity: str,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    fields: list[tuple[str, str | int]] = [
        ("lifecycle", publication.lifecycle),
        ("authority", "diagnostic_only"),
        ("publication_occurrence_identity", publication_occurrence_identity),
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
    if publication.allow_controlled:
        write_evaluated_parameters(
            path / "Parameters",
            publication.parameterization,
            publication.evaluation_result,
        )
    else:
        write_parameter_reports(
            path / "Parameters",
            publication.parameterization,
            publication.evaluation_result,
        )
    write_data(path / "Data", publication.plan, publication.evaluation_result)
    write_native_plots(path / "Plots", publication.plan, publication.evaluation_result)
    (path / "Statistics").mkdir()
    write_statistics(path, publication.plan, publication.evaluation_result, 0)


def _render_no_objective(path: Path) -> None:
    path.mkdir()
    diagnostics = path / "Diagnostics"
    diagnostics.mkdir()
    write_suppressed_outcome(
        diagnostics / "outcome.json",
        lifecycle="no_objective_data",
        operation_record={
            "artifact_type": "native_no_objective_data",
            "terminal": "no_objective_data",
        },
        accepted_result_identity=None,
        accepted_occurrence_identity=None,
        components=(),
    )


def _write_no_objective_manifest(
    path: Path,
    publication: NoObjectivePublication,
    publication_occurrence_identity: str,
) -> tuple[str, str, tuple[ArtifactReference, ...]]:
    return _write_manifest(
        path,
        (
            ("lifecycle", "no_objective_data"),
            ("authority", "evaluation_only"),
            ("publication_occurrence_identity", publication_occurrence_identity),
        ),
        publication.provenance,
    )


def _validate_and_publish_staged_step(
    staging: Path,
    destination: Path,
    manifest: tuple[str, str, tuple[ArtifactReference, ...]],
    expected_restart: CommittedRestartRecord | None,
) -> None:
    """Seal the complete staged #608 tree immediately before atomic publication."""
    manifest_identity, manifest_sha256, artifacts = manifest
    manifest_path = staging / "fit-manifest.toml"
    manifest_bytes = manifest_path.read_bytes()
    record = validate_native_step_manifest_bytes(manifest_bytes)
    if (
        record.get("manifest_identity") != manifest_identity
        or hashlib.sha256(manifest_bytes).hexdigest() != manifest_sha256
    ):
        raise NativePublicationError(
            "Staged native manifest changed before publication"
        )
    expected = {item.path: item for item in artifacts}
    actual_paths = {
        str(item.relative_to(staging))
        for item in staging.rglob("*")
        if item.is_file() and item != manifest_path
    }
    if actual_paths != set(expected):
        raise NativePublicationError("Staged native artifact catalogue changed")
    for relative_path, artifact in expected.items():
        if (
            hashlib.sha256((staging / relative_path).read_bytes()).hexdigest()
            != artifact.sha256
        ):
            raise NativePublicationError("Staged native artifact changed after hashing")
    if record.get("lifecycle") == "committed":
        restart_path = staging / "Parameters" / "restart.toml"
        restart_record = CommittedRestartRecord.from_record(
            json.loads(
                (staging / "Parameters" / "restart-provenance.json").read_text(
                    encoding="utf-8"
                )
            )
        )
        if (
            expected_restart is None
            or restart_record != expected_restart
            or restart_record.identity != record.get("restart_record_identity")
            or restart_record.restart_sha256
            != hashlib.sha256(restart_path.read_bytes()).hexdigest()
        ):
            raise NativePublicationError("Committed restart provenance changed")
    # Final ownership boundary: every ChemEx renderer and manifest writer has
    # completed above. Keep the validation adjacent to this atomic primitive;
    # staging is private and no application callback may run in between.
    publish_directory_noreplace(staging, destination)


def publish_native_results(
    path: Path,
    publication: NativePublication,
) -> PublishedStepReference:
    """Validate and atomically publish one native method-step result tree."""
    if isinstance(publication, CommittedFitPublication):
        _validate_committed_fit(publication)
    elif isinstance(publication, EvaluationPublication):
        _validate_evaluation_only(publication)
    elif isinstance(publication, NoObjectivePublication):
        _validate_no_objective(publication)
    else:
        _validate_suppressed(publication)
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(f"Native publication destination exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.parent / f".{destination.name}.tmp-{uuid4().hex}"
    publication_occurrence_identity = uuid4().hex
    restart: CommittedRestartRecord | None = None
    try:
        if isinstance(publication, CommittedFitPublication):
            restart = _render_committed_fit(staging, publication)
            manifest = _write_committed_manifest(
                staging,
                publication,
                publication_occurrence_identity,
                restart,
            )
            lifecycle = "committed"
            authority = "committed_fit"
            independent_ids = publication.parameterization.independent_ids
        elif isinstance(publication, EvaluationPublication):
            _render_evaluation(staging, publication)
            manifest = _write_evaluation_manifest(
                staging,
                publication,
                publication_occurrence_identity,
            )
            lifecycle = "successful_no_state_change"
            authority = "evaluation_only"
            independent_ids = publication.parameterization.independent_ids
        elif isinstance(publication, NoObjectivePublication):
            _render_no_objective(staging)
            manifest = _write_no_objective_manifest(
                staging,
                publication,
                publication_occurrence_identity,
            )
            lifecycle = "no_objective_data"
            authority = "evaluation_only"
            independent_ids = publication.parameterization.independent_ids
        else:
            _render_suppressed(staging, publication)
            manifest = _write_suppressed_manifest(
                staging,
                publication,
                publication_occurrence_identity,
            )
            lifecycle = publication.lifecycle
            authority = "diagnostic_only"
            independent_ids = publication.parameterization.independent_ids
        _validate_and_publish_staged_step(staging, destination, manifest, restart)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    manifest_identity, manifest_sha256, artifacts = manifest
    return _published_step_from_successful_native_publication(
        publication,
        publication_occurrence_identity,
        destination,
        lifecycle,
        authority,
        manifest_identity,
        manifest_sha256,
        publication.provenance,
        independent_ids,
        artifacts,
    )
