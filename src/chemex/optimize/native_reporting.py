"""Atomic step-root publication for native qualification results (#608).

This is an isolated qualification seam.  The production CLI remains on its
legacy reporting path until the native migration gates promote these artifacts.
"""

from __future__ import annotations

import ctypes
import errno
import math
import os
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal
from uuid import uuid4

from chemex.evaluation.native import EvaluationPlan, EvaluationResult
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CommitReceipt,
    DirectTrfExecution,
    DirectTrfTerminal,
    FitCommitOperation,
    FitCommitTerminal,
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
from chemex.parameters.parameterization import ActiveParameterization, ParameterRole
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
    accepted: AcceptedFitResult
    commit_receipt: CommitReceipt
    committed_snapshot: AnalysisValuesSnapshot
    commit_operation: FitCommitOperation
    analysis_values: AnalysisValues = field(repr=False, compare=False)
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
    accepted_result_identity: str | None = None
    accepted_occurrence_identity: str | None = None
    components: tuple[ComponentDiagnostic, ...] = ()
    partial_resampling: tuple[ResamplingEvidence, ...] = ()
    partial_mcmc: McmcEvidence | None = None


type NativePublication = (
    CommittedFitPublication | EvaluationPublication | SuppressedPublication
)

_AT_FDCWD = -100
_RENAME_NOREPLACE = 1
_RENAME_EXCL = 4


def _atomic_publish_noreplace(staging: Path, destination: Path) -> None:
    """Atomically rename a staged directory without replacing any destination."""
    if sys.platform == "win32":
        try:
            os.rename(staging, destination)
        except FileExistsError as error:
            raise FileExistsError(
                errno.EEXIST,
                "Native publication destination exists",
                destination,
            ) from error
        return

    source_bytes = os.fsencode(staging)
    destination_bytes = os.fsencode(destination)
    if sys.platform == "linux":
        function_name = "renameat2"
        argument_types = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        arguments = (
            _AT_FDCWD,
            source_bytes,
            _AT_FDCWD,
            destination_bytes,
            _RENAME_NOREPLACE,
        )
    elif sys.platform == "darwin":
        function_name = "renamex_np"
        argument_types = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
        arguments = (source_bytes, destination_bytes, _RENAME_EXCL)
    else:
        raise OSError(
            errno.ENOTSUP,
            "Atomic no-replace directory rename is unavailable on this platform",
            destination,
        )

    libc = ctypes.CDLL(None, use_errno=True)
    try:
        rename_noreplace = getattr(libc, function_name)
    except AttributeError as error:
        raise OSError(
            errno.ENOSYS,
            f"{function_name} is unavailable for atomic native publication",
            destination,
        ) from error
    rename_noreplace.argtypes = argument_types
    rename_noreplace.restype = ctypes.c_int
    result = rename_noreplace(*arguments)
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        raise FileExistsError(
            errno.EEXIST,
            "Native publication destination exists",
            destination,
        )
    raise OSError(error_number, os.strerror(error_number), destination)


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
        ):
            raise NativePublicationError(
                "Suppressed workflows may retain only genuine partial MCMC evidence"
            )


def _validate_committed_fit(publication: CommittedFitPublication) -> None:
    accepted = publication.accepted
    receipt = publication.commit_receipt
    committed = publication.committed_snapshot
    operation = publication.commit_operation
    current = publication.analysis_values.snapshot()
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
) -> None:
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
        elif isinstance(publication, EvaluationPublication):
            _render_evaluation(staging, publication)
        else:
            _render_suppressed(staging, publication)
        _atomic_publish_noreplace(staging, destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
