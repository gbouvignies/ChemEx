"""Closed composition of one native method step.

The production deterministic dispatcher composes already-compiled native
artifacts here while preserving the live AnalysisValues occurrence boundary.
"""

from __future__ import annotations

import hashlib
import json
import weakref
from argparse import Namespace
from collections.abc import Callable
from dataclasses import dataclass, field, fields, is_dataclass, replace
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal, Protocol, Self, cast
from uuid import uuid4

from chemex.configuration.methods import Method
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.native_provenance import (
    BaselineReference,
    BudgetRecord,
    PolicyRecord,
    ProvenanceEnvironment,
    PublishedStepReference,
    SeedRecord,
    WorkflowProvenance,
    published_step_reference_identity,
)
from chemex.optimize.de_direct_trf import (
    AcceptedDeDirectTrfResult,
    DeDirectTrfInterrupted,
    DeDirectTrfInvocation,
    DeDirectTrfOutcome,
    execute_de_direct_trf,
    validate_de_commit_lineage,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    DirectTrfInvocation,
    FinalResidualJacobianEvidence,
    FitCommitOperation,
    FitCommitTerminal,
    LiveFitCommitAuthority,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
    execute_fit_commit,
    fit_commit_authority_is_authoritative,
    fit_commit_operation_is_authoritative,
)
from chemex.optimize.grid_direct_trf import (
    AcceptedGridDirectTrfResult,
    GridDirectTrfInterrupted,
    GridDirectTrfInvocation,
    GridDirectTrfOutcome,
    validate_grid_commit_lineage,
)
from chemex.optimize.grouped_direct_trf import (
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfOutcome,
    execute_grouped_direct_trf,
)
from chemex.optimize.grouped_grid_direct_trf import (
    GroupedGridDirectTrfOutcome,
    GroupedGridSeedOutcome,
    execute_grouped_grid_direct_trf,
)
from chemex.optimize.native_mcmc import (
    McmcEvidence,
    McmcOperationTerminal,
    McmcPlan,
    PosteriorSampleEvidence,
    PosteriorSummary,
    ResolvedMcmcPolicy,
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
    MethodStepPrimaryRecord,
    MethodStepPrimaryTerminal,
    NoObjectivePublication,
    ResamplingPublication,
    ResamplingSummaryFailurePublication,
    SuppressedPublication,
    publication_provenance_records,
    publish_native_results,
)
from chemex.optimize.native_resampling import (
    OperationTerminal as ResamplingOperationTerminal,
)
from chemex.optimize.native_resampling import (
    OptimizationStrategy,
    ResamplingDatasetManifest,
    ResamplingEvidence,
    ResamplingPlan,
    ResamplingScheme,
    ResamplingSummaryPolicy,
    SummaryFailure,
    execute_resampling_evidence,
    summarize_resampling_evidence,
)
from chemex.optimize.progress import ContextualProgressObserver
from chemex.optimize.uncertainty import (
    CompiledConstraintLinearizationCapabilities,
    ParameterUnit,
    UncertaintyEvidence,
    UncertaintyPolicy,
    derive_uncertainty_evidence,
)
from chemex.optimize.uncertainty import (
    OperationTerminal as UncertaintyOperationTerminal,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.run_info import (
    CapturedNativeInputs,
    NativeRunInformation,
    write_native_run_info,
)
from chemex.runtime import ExecutionSettings


class EvaluationPurpose(StrEnum):
    """Closed meanings of an evaluation-only method step."""

    EVALUATE_ONLY = "evaluate_only"
    NO_OPTIMIZATION_REQUIRED = "no_optimization_required"
    NO_OBJECTIVE_DATA = "no_objective_data"


class MethodStepLifecycle(StrEnum):
    """Closed terminal lifecycle of one native method-step occurrence."""

    SUCCESSFUL_NO_STATE_CHANGE = "successful_no_state_change"
    NO_OBJECTIVE_DATA = "no_objective_data"
    COMMITTED = "committed"
    ACCEPTED_UNCOMMITTED = "accepted_uncommitted"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    PUBLICATION_FAILED = "publication_failed"


class MethodStepStrategy(StrEnum):
    """Closed primary optimization strategies supported by the composer."""

    DIRECT_TRF = "direct_trf"
    GRID_DIRECT_TRF = "grid_direct_trf"
    DE_DIRECT_TRF = "de_direct_trf"


class MethodStepCheckpoint(StrEnum):
    """Occurrence-only observation points that cannot change workflow identity."""

    AGGREGATE_ACCEPTED = "aggregate_accepted"
    COMMIT_COMPLETED = "commit_completed"


class DerivationDisposition(StrEnum):
    """Closed terminal truth for a requested post-commit evidence stage."""

    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    BLOCKED_BY_PREREQUISITE = "blocked_by_prerequisite"
    NOT_STARTED_BY_WORKFLOW_STOP = "not_started_by_workflow_stop"


@dataclass(frozen=True, slots=True)
class UncertaintyDerivationRequest:
    """Exact covariance and optional constrained-propagation request."""

    policy: UncertaintyPolicy
    constrained_scope: tuple[str, ...] = ()
    constrained_units: tuple[tuple[str, ParameterUnit], ...] = ()
    constrained_scales: tuple[tuple[str, float], ...] = ()
    compiled_capabilities: CompiledConstraintLinearizationCapabilities | None = None
    resolved_environment_identity: str = ""
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.resolved_environment_identity:
            raise ValueError("Uncertainty request requires an environment identity")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-uncertainty-request-v1",
                (
                    self.policy.identity,
                    self.constrained_scope,
                    tuple((key, value.value) for key, value in self.constrained_units),
                    tuple(
                        (key, float(value).hex())
                        for key, value in self.constrained_scales
                    ),
                    None
                    if self.compiled_capabilities is None
                    else self.compiled_capabilities.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ResamplingDerivationRequest:
    """Exact deterministic resampling population requested after commit."""

    references: tuple[bool, ...]
    nucleus_groups: tuple[str, ...]
    observation_descriptors: tuple[str, ...]
    scheme: ResamplingScheme
    replicate_count: int
    replicate_structural_identities: tuple[str, ...]
    replicate_component_identities: tuple[tuple[str, ...], ...]
    root_seed: int
    output_scope: tuple[str, ...]
    output_units: tuple[str, ...]
    minimum_successful_count: int
    strategy: OptimizationStrategy
    strategy_settings: tuple[tuple[str, str], ...]
    summary_policy: ResamplingSummaryPolicy
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-resampling-request-v1",
                (
                    self.references,
                    self.nucleus_groups,
                    self.observation_descriptors,
                    self.scheme.value,
                    self.replicate_count,
                    self.replicate_structural_identities,
                    self.replicate_component_identities,
                    self.root_seed,
                    self.output_scope,
                    self.output_units,
                    self.minimum_successful_count,
                    self.strategy.value,
                    self.strategy_settings,
                    self.summary_policy.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class McmcDerivationRequest:
    """Exact MCMC execution and typed posterior-summary request."""

    policy: ResolvedMcmcPolicy
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    output_units: tuple[tuple[str, ParameterUnit], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-mcmc-request-v1",
                (
                    self.policy.identity,
                    tuple((key, value.value) for key, value in self.coordinate_units),
                    tuple((key, value.value) for key, value in self.output_units),
                ),
            ),
        )


type MethodStepDerivationRequest = (
    UncertaintyDerivationRequest | ResamplingDerivationRequest | McmcDerivationRequest
)


@dataclass(frozen=True, slots=True)
class DerivationOutcome:
    """One requested stage disposition referencing its existing source artifacts."""

    request_identity: str
    stage: str
    disposition: DerivationDisposition
    operation_identity: str | None = None
    artifact_identities: tuple[str, ...] = ()
    message: str = ""
    operation: object | None = field(default=None, repr=False, compare=False)
    artifacts: tuple[object, ...] = field(default=(), repr=False, compare=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.operation is not None and getattr(self.operation, "identity", None) != (
            self.operation_identity
        ):
            raise ValueError("Derivation operation identity differs from its artifact")
        if self.artifacts and tuple(
            getattr(item, "identity", None) for item in self.artifacts
        ) != (self.artifact_identities):
            raise ValueError("Derivation artifact identities differ from their sources")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-derivation-outcome-v1",
                (
                    self.request_identity,
                    self.stage,
                    self.disposition.value,
                    self.operation_identity,
                    self.artifact_identities,
                    self.message,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class MethodStepPublicationRequest:
    """Exact native step-root publication destination and provenance context."""

    path: Path
    environment: ProvenanceEnvironment
    baseline_references: tuple[BaselineReference, ...]
    run_provenance: MethodStepRunProvenanceRequest | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        destination = Path(self.path).resolve()
        object.__setattr__(self, "path", destination)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-publication-request-v1",
                (
                    str(destination),
                    self.environment.identity,
                    tuple(
                        (
                            item.kind,
                            item.identity,
                            item.occurrence_identity,
                            item.result_bundle_identity,
                        )
                        for item in self.baseline_references
                    ),
                    None
                    if self.run_provenance is None
                    else self.run_provenance.identity,
                    "native-step-root-no-clobber-v1",
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class MethodStepRunProvenanceRequest:
    """Occurrence-owned #609 run-information request captured before execution."""

    output_directory: Path
    invocation_identity: str
    inputs: CapturedNativeInputs
    starting_snapshot: AnalysisValuesSnapshot
    prior_steps: tuple[PublishedStepReference, ...] = ()
    argv: tuple[str, ...] = ()
    working_directory: Path = field(default_factory=Path.cwd)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        output = Path(self.output_directory).resolve()
        working = Path(self.working_directory).resolve()
        if not self.invocation_identity.strip():
            raise ValueError("Run-provenance invocation identity cannot be empty")
        object.__setattr__(self, "output_directory", output)
        object.__setattr__(self, "working_directory", working)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-method-step-run-provenance-request-v1",
                (
                    str(output),
                    self.invocation_identity,
                    self.inputs.identity,
                    _snapshot_identity(self.starting_snapshot),
                    tuple(item.identity for item in self.prior_steps),
                    self.argv,
                    str(working),
                ),
            ),
        )


type MethodStepInvocation = (
    GroupedDirectTrfInvocation | GridDirectTrfInvocation | DeDirectTrfInvocation
)
type MethodStepPrimaryExecution = (
    GroupedDirectTrfOutcome
    | GroupedGridDirectTrfOutcome
    | GridDirectTrfOutcome
    | DeDirectTrfOutcome
)


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        (kind, record),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


class _IdentifiedArtifact(Protocol):
    @property
    def identity(self) -> str: ...


def _snapshot_identity(snapshot: AnalysisValuesSnapshot) -> str:
    return hashlib.sha256(snapshot.to_json().encode("ascii")).hexdigest()


def _validate_snapshot_integrity(snapshot: AnalysisValuesSnapshot) -> None:
    try:
        restored = AnalysisValuesSnapshot.from_json(snapshot.to_json())
    except Exception as error:
        raise ValueError("Method-step snapshot integrity validation failed") from error
    if restored != snapshot or tuple(snapshot.items()) != tuple(restored.items()):
        raise ValueError("Method-step snapshot integrity validation failed")


def _validate_rederived_identity(
    value: object,
    *,
    identity_attribute: str = "identity",
) -> None:
    """Reconstruct one frozen declarative child and compare its current identity."""
    try:
        expected = replace(cast("Any", value))
    except Exception as error:
        raise ValueError(
            f"Method-step child integrity validation failed for {type(value).__name__}"
        ) from error
    if getattr(expected, identity_attribute) != getattr(value, identity_attribute):
        raise ValueError(
            f"Method-step child integrity validation failed for {type(value).__name__}"
        )


def _validate_recursive_dataclass_children(  # noqa: C901 - closed evidence tree
    value: object,
    *,
    seen: set[int] | None = None,
) -> None:
    """Reconstruct the closed native evidence tree while skipping live authorities."""
    visited = set() if seen is None else seen
    if id(value) in visited:
        return
    visited.add(id(value))
    if isinstance(value, tuple):
        for item in value:
            _validate_recursive_dataclass_children(item, seen=visited)
        return
    if isinstance(value, AnalysisValuesSnapshot):
        _validate_snapshot_integrity(value)
        return
    if isinstance(value, FinalResidualJacobianEvidence):
        value.validate_integrity()
        return
    if isinstance(value, (AcceptedFitResult, FitCommitOperation)) or not is_dataclass(
        value
    ):
        return
    for declared in fields(value):
        if "witness" in declared.name or "authority" in declared.name:
            continue
        _validate_recursive_dataclass_children(
            getattr(value, declared.name),
            seen=visited,
        )
    if type(value).__module__.startswith("chemex."):
        reconstructed = replace(cast("Any", value))
        for identity_attribute in ("identity", "fingerprint"):
            if hasattr(value, identity_attribute) and getattr(
                reconstructed, identity_attribute
            ) != getattr(value, identity_attribute):
                raise ValueError(
                    "Method-step recursive child integrity validation failed for "
                    f"{type(value).__name__}"
                )


def _semantic_method_record(method: Method) -> dict[str, object]:
    """Return the closed scientific/lifecycle Method fields used by #641."""
    # Profile selection has already been materialized in the evaluation plan,
    # whose identity captures the selected profiles and observations.  Exclude
    # the parsed SpinSystem objects here just as normalized method provenance
    # does, avoiding a second and non-canonical serialization of that state.
    record = cast(
        "dict[str, object]",
        method.model_dump(mode="json", exclude={"include", "exclude"}),
    )
    statistics = record.get("statistics")
    if isinstance(statistics, dict):
        mcmc = statistics.get("mcmc")
        if isinstance(mcmc, dict):
            # Worker allocation changes execution, not the requested posterior.
            mcmc.pop("workers", None)
    return record


def _operational_method_record(method: Method) -> dict[str, object]:
    """Return Method-owned execution controls excluded from semantic identity."""
    mcmc = None if method.statistics is None else method.statistics.mcmc
    return {
        "mcmc_workers": None if mcmc is None else mcmc.workers,
    }


def _method_reconstruction_record(method: Method) -> dict[str, object]:
    """Return validation input without serializing SpinSystem internals."""
    record = cast(
        "dict[str, object]",
        method.model_dump(exclude={"include", "exclude"}, exclude_none=True),
    )
    for field_name in ("include", "exclude"):
        selection = getattr(method, field_name)
        if isinstance(selection, list):
            record[field_name] = [str(item) for item in selection]
        elif selection is not None:
            record[field_name] = selection
    return record


def _direct_policy_semantics(invocation: DirectTrfInvocation) -> tuple[object, ...]:
    return (
        invocation.objective_request_budget,
        tuple(float(value).hex() for value in invocation.x_scale),
        tuple(
            None if value is None else float(value).hex()
            for value in (
                invocation.ftol,
                invocation.xtol,
                invocation.gtol,
            )
        ),
    )


def _optimization_semantics(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    strategy: MethodStepStrategy,
    invocation: MethodStepInvocation,
) -> dict[str, object]:
    """Value/lifecycle semantics with occurrence and runtime facts removed."""
    problem_record = (
        problem.evaluation_plan_identity,
        problem.evaluator_parameterization_identity,
        problem.constraint_program_identity,
        problem.configuration_identity,
        tuple((key, float(value).hex()) for key, value in problem.independent_items),
        problem.controlled_ids,
        tuple((key, float(value).hex()) for key, value in problem.held_items),
        tuple(float(value).hex() for value in problem.start),
        tuple(float(value).hex() for value in problem.lower_bounds),
        tuple(float(value).hex() for value in problem.upper_bounds),
        problem.commit_scope,
        problem.affine_feasibility_identity,
        problem.scalarization_version,
    )
    decomposition_record = tuple(
        (component.controlled_ids, component.root_profile_indices)
        for component in decomposition.components
    )
    if isinstance(invocation, GroupedDirectTrfInvocation):
        invocation_record: object = tuple(
            _direct_policy_semantics(item) for item in invocation.component_invocations
        )
    elif isinstance(invocation, GridDirectTrfInvocation):
        invocation_record = (
            tuple(axis.identity for axis in invocation.axes),
            invocation.objective_request_budget,
            tuple(float(value).hex() for value in invocation.x_scale),
            tuple(
                None if value is None else float(value).hex()
                for value in (invocation.ftol, invocation.xtol, invocation.gtol)
            ),
        )
    else:
        invocation_record = (
            tuple(item.identity for item in invocation.search_coordinates),
            invocation.root_seed,
            invocation.population.identity,
            invocation.de_objective_request_budget,
            invocation.polish_objective_request_budget,
            invocation.mutation,
            float(invocation.recombination).hex(),
            float(invocation.tol).hex(),
            float(invocation.atol).hex(),
            tuple(float(value).hex() for value in invocation.polish_x_scale),
            tuple(
                None if value is None else float(value).hex()
                for value in (
                    invocation.polish_ftol,
                    invocation.polish_xtol,
                    invocation.polish_gtol,
                )
            ),
        )
    return {
        "kind": "optimization",
        "strategy": strategy.value,
        "problem": problem_record,
        "decomposition": decomposition_record,
        "invocation": invocation_record,
    }


@dataclass(frozen=True, slots=True)
class MethodStepWorkflow:
    """One immutable exact native method-step composition root."""

    starting_snapshot: AnalysisValuesSnapshot = field(repr=False, compare=False)
    parameter_model: SealedParameterModel = field(repr=False, compare=False)
    parameterization: ActiveParameterization = field(repr=False, compare=False)
    engine: EvaluationEngine = field(repr=False, compare=False)
    method: Method = field(repr=False, compare=False)
    evaluation_purpose: EvaluationPurpose | None = None
    problem: OptimizationProblem | None = field(default=None, repr=False, compare=False)
    decomposition: FitDecomposition | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    strategy: MethodStepStrategy | None = None
    invocation: MethodStepInvocation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    derivations: tuple[MethodStepDerivationRequest, ...] = ()
    publication: MethodStepPublicationRequest | None = None
    semantic_identity: str = field(init=False)
    binding_identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - closed construction boundary
        snapshot = self.starting_snapshot
        program = self.parameterization.program
        if (
            self.parameter_model.identity != program.parameter_model_identity
            or snapshot.occurrence_identity != self.parameterization.occurrence_identity
            or snapshot.revision != self.parameterization.source_revision
            or snapshot.model_identity != program.model_identity
            or snapshot.definitions_identity != program.definitions_identity
            or snapshot.configuration_identity != program.configuration_identity
            or self.engine.plan.parameterization_identity
            != self.parameterization.evaluator_identity
            or self.engine.plan.constraint_program_identity != program.fingerprint
        ):
            raise ValueError("Method-step workflow has an incompatible starting state")
        retained = self.engine.plan.retained_observation_count
        if self.evaluation_purpose is not None:
            if any(
                value is not None
                for value in (
                    self.problem,
                    self.decomposition,
                    self.strategy,
                    self.invocation,
                )
            ):
                raise ValueError("Evaluation-only workflow cannot own optimization")
            fitted = tuple(
                param_id
                for param_id in self.parameterization.independent_ids
                if self.parameterization.role(param_id) is ParameterRole.FIT
            )
            if self.evaluation_purpose is EvaluationPurpose.NO_OBJECTIVE_DATA:
                if retained:
                    raise ValueError("NO_OBJECTIVE_DATA requires an empty objective")
            elif retained < 1:
                raise ValueError(
                    "Evaluation-only workflow requires retained objective data"
                )
            elif (
                self.evaluation_purpose is EvaluationPurpose.NO_OPTIMIZATION_REQUIRED
                and fitted
            ):
                raise ValueError(
                    "NO_OPTIMIZATION_REQUIRED cannot retain controlled coordinates"
                )
            primary_record: dict[str, object] = {
                "kind": "evaluation_only",
                "purpose": self.evaluation_purpose.value,
            }
        else:
            if (
                self.problem is None
                or self.decomposition is None
                or self.strategy is None
                or self.invocation is None
            ):
                raise ValueError("Optimization workflow is missing its primary inputs")
            if retained < 1:
                raise ValueError("Optimization workflow requires objective data")
            if self.problem.source_snapshot is not snapshot:
                raise ValueError("Optimization problem has another starting state")
            try:
                exact_decomposition = FitDecomposition.from_root(
                    self.problem,
                    self.parameterization,
                    self.engine,
                )
            except Exception as error:
                raise ValueError(
                    f"Method-step optimization decomposition is invalid: {error}"
                ) from error
            if (
                not self.decomposition.components
                or exact_decomposition.identity != self.decomposition.identity
                or self.decomposition.root_problem_identity != self.problem.identity
            ):
                raise ValueError(
                    "Method-step optimization requires its exact non-empty decomposition"
                )
            if self.strategy is MethodStepStrategy.DIRECT_TRF:
                direct_invocation = cast(
                    "GroupedDirectTrfInvocation",
                    self.invocation,
                )
                supported = isinstance(self.invocation, GroupedDirectTrfInvocation)
                supported = supported and (
                    direct_invocation.decomposition_identity
                    == self.decomposition.identity
                    and direct_invocation.root_problem_identity == self.problem.identity
                )
                if (
                    supported
                    and len(
                        {
                            item.execution_settings
                            for item in direct_invocation.component_invocations
                        }
                    )
                    != 1
                ):
                    raise ValueError(
                        "Grouped Direct TRF requires one resolved execution setting"
                    )
            elif self.strategy is MethodStepStrategy.GRID_DIRECT_TRF:
                supported = isinstance(self.invocation, GridDirectTrfInvocation)
                supported = supported and (
                    self.invocation.root_problem_identity == self.problem.identity
                )
            else:
                supported = isinstance(self.invocation, DeDirectTrfInvocation)
                supported = supported and (
                    self.invocation.root_problem_identity == self.problem.identity
                    and len(self.decomposition.components) == 1
                )
            if not supported:
                raise ValueError(
                    "Primary strategy and invocation are unsupported or incompatible"
                )
            primary_record = _optimization_semantics(
                self.problem,
                self.decomposition,
                self.strategy,
                self.invocation,
            )
        ordered_derivations = tuple(
            sorted(
                self.derivations,
                key=lambda item: {
                    "uncertainty": 0,
                    "resampling": 1,
                    "mcmc": 2,
                }[_derivation_stage(item)],
            )
        )
        object.__setattr__(self, "derivations", ordered_derivations)
        derivation_kinds = tuple(
            _derivation_stage(item) for item in ordered_derivations
        )
        if len(set(derivation_kinds)) != len(derivation_kinds):
            raise ValueError("Method-step derivation kinds must be unique")
        if self.evaluation_purpose is not None and self.derivations:
            raise ValueError("Evaluation-only workflow cannot request fit derivations")
        for item in self.derivations:
            if isinstance(item, ResamplingDerivationRequest) and any(
                len(values) != self.engine.plan.observation_count
                for values in (
                    item.references,
                    item.nucleus_groups,
                    item.observation_descriptors,
                )
            ):
                raise ValueError(
                    "Resampling request metadata must cover the exact EvaluationPlan"
                )
        semantic_identity, binding_identity = _rederive_workflow_identities(
            self,
            primary_record=primary_record,
            derivation_kinds=derivation_kinds,
        )
        object.__setattr__(self, "semantic_identity", semantic_identity)
        object.__setattr__(self, "binding_identity", binding_identity)

    def validate_integrity(self) -> None:
        """Rederive this workflow and every identity-bearing child."""
        _validate_method_step_workflow(self)

    @classmethod
    def for_evaluation(
        cls,
        *,
        starting_snapshot: AnalysisValuesSnapshot,
        parameter_model: SealedParameterModel,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
        method: Method,
        purpose: EvaluationPurpose,
        publication: MethodStepPublicationRequest | None = None,
    ) -> Self:
        """Bind an explicit evaluation-only meaning to exact native inputs."""
        return cls(
            starting_snapshot,
            parameter_model,
            parameterization,
            engine,
            method,
            evaluation_purpose=purpose,
            publication=publication,
        )

    @classmethod
    def for_optimization(
        cls,
        *,
        starting_snapshot: AnalysisValuesSnapshot,
        parameter_model: SealedParameterModel,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
        method: Method,
        problem: OptimizationProblem,
        decomposition: FitDecomposition,
        strategy: MethodStepStrategy,
        invocation: MethodStepInvocation,
        derivations: tuple[MethodStepDerivationRequest, ...] = (),
        publication: MethodStepPublicationRequest | None = None,
    ) -> Self:
        """Bind one explicit supported strategy to its exact decomposition."""
        return cls(
            starting_snapshot,
            parameter_model,
            parameterization,
            engine,
            method,
            problem=problem,
            decomposition=decomposition,
            strategy=strategy,
            invocation=invocation,
            derivations=derivations,
            publication=publication,
        )


def _rederive_workflow_identities(
    workflow: MethodStepWorkflow,
    *,
    primary_record: dict[str, object] | None = None,
    derivation_kinds: tuple[str, ...] | None = None,
) -> tuple[str, str]:
    """Compute the two workflow identities solely from current declarative fields."""
    if primary_record is None:
        if workflow.evaluation_purpose is not None:
            primary_record = {
                "kind": "evaluation_only",
                "purpose": workflow.evaluation_purpose.value,
            }
        else:
            primary_record = _optimization_semantics(
                cast("OptimizationProblem", workflow.problem),
                cast("FitDecomposition", workflow.decomposition),
                cast("MethodStepStrategy", workflow.strategy),
                cast("MethodStepInvocation", workflow.invocation),
            )
    if derivation_kinds is None:
        derivation_kinds = tuple(
            _derivation_stage(item) for item in workflow.derivations
        )
    semantic_identity = _identity(
        "native-method-step-workflow-semantics-v1",
        {
            "primary": primary_record,
            "parameter_model_identity": workflow.parameter_model.identity,
            "parameterization_program": workflow.parameterization.program.fingerprint,
            "evaluation_plan_identity": workflow.engine.plan.identity,
            "method": _semantic_method_record(workflow.method),
            "publication": (
                None
                if workflow.publication is None
                else "native-step-root-no-clobber-v1"
            ),
            "derivations": tuple(
                (kind, item.identity)
                for kind, item in zip(
                    derivation_kinds,
                    workflow.derivations,
                    strict=True,
                )
            ),
            "failure_policy": "closed-method-step-failure-v1",
        },
    )
    snapshot = workflow.starting_snapshot
    binding_identity = _identity(
        "native-method-step-workflow-binding-v1",
        {
            "semantic_identity": semantic_identity,
            "starting_state_identity": _snapshot_identity(snapshot),
            "occurrence_identity": snapshot.occurrence_identity,
            "revision": snapshot.revision,
            "parameterization_identity": workflow.parameterization.identity,
            "operational_method": _operational_method_record(workflow.method),
            "invocation_identity": (
                None if workflow.invocation is None else workflow.invocation.identity
            ),
            "uncertainty_environments": tuple(
                item.resolved_environment_identity
                for item in workflow.derivations
                if isinstance(item, UncertaintyDerivationRequest)
            ),
            "publication_request_identity": (
                None if workflow.publication is None else workflow.publication.identity
            ),
        },
    )
    return semantic_identity, binding_identity


def _validate_derivation_request_integrity(
    request: MethodStepDerivationRequest,
) -> None:
    if isinstance(request, UncertaintyDerivationRequest):
        _validate_rederived_identity(request.policy)
        if request.compiled_capabilities is not None:
            _validate_rederived_identity(request.compiled_capabilities)
    elif isinstance(request, ResamplingDerivationRequest):
        _validate_rederived_identity(request.summary_policy)
    else:
        _validate_rederived_identity(request.policy)
    _validate_rederived_identity(request)


def _validate_invocation_integrity(invocation: MethodStepInvocation) -> None:
    if isinstance(invocation, GroupedDirectTrfInvocation):
        for component in invocation.component_invocations:
            _validate_rederived_identity(component)
    elif isinstance(invocation, GridDirectTrfInvocation):
        for axis in invocation.axes:
            _validate_rederived_identity(axis)
        for seed in invocation.seeds:
            if seed.problem is not None:
                _validate_snapshot_integrity(seed.problem.source_snapshot)
                _validate_rederived_identity(seed.problem)
            if seed.invocation is not None:
                _validate_rederived_identity(seed.invocation)
            _validate_rederived_identity(seed)
    else:
        _validate_snapshot_integrity(invocation.root_problem.source_snapshot)
        _validate_snapshot_integrity(invocation.search_problem.source_snapshot)
        _validate_rederived_identity(invocation.root_problem)
        _validate_rederived_identity(invocation.search_problem)
        for coordinate in invocation.search_coordinates:
            _validate_rederived_identity(coordinate)
        _validate_rederived_identity(invocation.population)
    _validate_rederived_identity(invocation)


def _validate_publication_request_integrity(
    publication: MethodStepPublicationRequest,
) -> None:
    _validate_rederived_identity(publication.environment)
    for reference in publication.baseline_references:
        _validate_rederived_identity(reference)
    if publication.run_provenance is not None:
        run = publication.run_provenance
        _validate_snapshot_integrity(run.starting_snapshot)
        _validate_rederived_identity(run.inputs)
        for reference in run.prior_steps:
            reference.require_exact_live_publication()
        _validate_rederived_identity(run)
    _validate_rederived_identity(publication)


def _validate_method_step_workflow(  # noqa: C901 - closed recursive workflow tree
    workflow: MethodStepWorkflow,
) -> None:
    """Canonically rederive one workflow without trusting construction-time hashes."""
    _validate_snapshot_integrity(workflow.starting_snapshot)
    _validate_rederived_identity(workflow.parameter_model)
    _validate_rederived_identity(
        workflow.parameterization.program,
        identity_attribute="fingerprint",
    )
    _validate_rederived_identity(workflow.parameterization)
    try:
        method_record = _method_reconstruction_record(workflow.method)
        canonical_method = Method.model_validate(method_record)
    except Exception as error:
        raise ValueError("Method-step Method integrity validation failed") from error
    if _method_reconstruction_record(canonical_method) != method_record:
        raise ValueError("Method-step Method integrity validation failed")
    try:
        restored_plan = type(workflow.engine.plan).from_record(
            workflow.engine.plan.to_record()
        )
    except Exception as error:
        raise ValueError(
            "Method-step evaluation-plan integrity validation failed"
        ) from error
    if restored_plan.identity != workflow.engine.plan.identity:
        raise ValueError("Method-step evaluation-plan integrity validation failed")
    if workflow.engine._parameterization is not workflow.parameterization:
        raise ValueError("Method-step engine integrity validation failed")
    if workflow.problem is not None:
        _validate_snapshot_integrity(workflow.problem.source_snapshot)
        _validate_rederived_identity(workflow.problem)
    if workflow.decomposition is not None:
        for component in workflow.decomposition.components:
            _validate_snapshot_integrity(component.problem.source_snapshot)
            _validate_rederived_identity(component.problem)
            _validate_rederived_identity(component)
        _validate_rederived_identity(workflow.decomposition.partition_proof)
        _validate_rederived_identity(workflow.decomposition)
    if workflow.invocation is not None:
        _validate_invocation_integrity(workflow.invocation)
    for request in workflow.derivations:
        _validate_derivation_request_integrity(request)
    if workflow.publication is not None:
        _validate_publication_request_integrity(workflow.publication)
    retained = workflow.engine.plan.retained_observation_count
    if workflow.evaluation_purpose is not None:
        if any(
            item is not None
            for item in (
                workflow.problem,
                workflow.decomposition,
                workflow.strategy,
                workflow.invocation,
            )
        ) or (
            workflow.evaluation_purpose is not EvaluationPurpose.NO_OBJECTIVE_DATA
            and retained < 1
        ):
            raise ValueError("Method-step workflow integrity validation failed")
    else:
        problem = workflow.problem
        decomposition = workflow.decomposition
        strategy = workflow.strategy
        invocation = workflow.invocation
        if (
            problem is None
            or decomposition is None
            or strategy is None
            or invocation is None
            or retained < 1
            or problem.source_snapshot is not workflow.starting_snapshot
            or decomposition.root_problem_identity != problem.identity
            or decomposition.root_plan_identity != workflow.engine.plan.identity
        ):
            raise ValueError("Method-step workflow integrity validation failed")
        if strategy is MethodStepStrategy.DIRECT_TRF:
            compatible = (
                isinstance(invocation, GroupedDirectTrfInvocation)
                and invocation.root_problem_identity == problem.identity
                and invocation.decomposition_identity == decomposition.identity
            )
        elif strategy is MethodStepStrategy.GRID_DIRECT_TRF:
            compatible = (
                isinstance(invocation, GridDirectTrfInvocation)
                and invocation.root_problem_identity == problem.identity
            )
        else:
            compatible = (
                isinstance(invocation, DeDirectTrfInvocation)
                and invocation.root_problem_identity == problem.identity
                and len(decomposition.components) == 1
            )
        if not compatible:
            raise ValueError("Method-step workflow integrity validation failed")
    ordered = tuple(
        sorted(
            workflow.derivations,
            key=lambda item: {"uncertainty": 0, "resampling": 1, "mcmc": 2}[
                _derivation_stage(item)
            ],
        )
    )
    derivation_kinds = tuple(_derivation_stage(item) for item in workflow.derivations)
    if (
        ordered != workflow.derivations
        or len(set(derivation_kinds)) != len(derivation_kinds)
        or (workflow.evaluation_purpose is not None and workflow.derivations)
    ):
        raise ValueError("Method-step workflow integrity validation failed")
    semantic_identity, binding_identity = _rederive_workflow_identities(workflow)
    if (
        semantic_identity != workflow.semantic_identity
        or binding_identity != workflow.binding_identity
    ):
        raise ValueError("Method-step workflow integrity validation failed")


@dataclass(frozen=True, slots=True)
class _ExecutionStartWorkflowContext:
    semantic_identity: str
    binding_identity: str
    starting_snapshot: AnalysisValuesSnapshot
    starting_snapshot_object: AnalysisValuesSnapshot
    parameter_model: SealedParameterModel
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    method: Method
    problem: OptimizationProblem | None
    decomposition: FitDecomposition | None
    strategy: MethodStepStrategy | None
    invocation: MethodStepInvocation | None
    derivations: tuple[MethodStepDerivationRequest, ...]
    publication: MethodStepPublicationRequest | None
    evaluation_purpose: EvaluationPurpose | None


def _pin_execution_start_workflow(
    workflow: MethodStepWorkflow,
) -> _ExecutionStartWorkflowContext:
    """Capture an execution-owned context after full recursive validation."""
    workflow.validate_integrity()
    semantic_identity, binding_identity = _rederive_workflow_identities(workflow)
    return _ExecutionStartWorkflowContext(
        semantic_identity,
        binding_identity,
        AnalysisValuesSnapshot.from_json(workflow.starting_snapshot.to_json()),
        workflow.starting_snapshot,
        workflow.parameter_model,
        workflow.parameterization,
        workflow.engine,
        workflow.method,
        workflow.problem,
        workflow.decomposition,
        workflow.strategy,
        workflow.invocation,
        workflow.derivations,
        workflow.publication,
        workflow.evaluation_purpose,
    )


def _validate_execution_start_workflow(
    workflow: MethodStepWorkflow,
    expected: _ExecutionStartWorkflowContext,
) -> None:
    """Reject callback retargeting, including a fully rehashed replacement."""
    workflow.validate_integrity()
    semantic_identity, binding_identity = _rederive_workflow_identities(workflow)
    if (
        semantic_identity != expected.semantic_identity
        or binding_identity != expected.binding_identity
        or workflow.starting_snapshot is not expected.starting_snapshot_object
        or workflow.starting_snapshot != expected.starting_snapshot
        or _snapshot_identity(workflow.starting_snapshot)
        != _snapshot_identity(expected.starting_snapshot)
        or workflow.parameter_model is not expected.parameter_model
        or workflow.parameterization is not expected.parameterization
        or workflow.engine is not expected.engine
        or workflow.method is not expected.method
        or workflow.problem is not expected.problem
        or workflow.decomposition is not expected.decomposition
        or workflow.strategy is not expected.strategy
        or workflow.invocation is not expected.invocation
        or workflow.derivations is not expected.derivations
        or workflow.publication is not expected.publication
        or workflow.evaluation_purpose is not expected.evaluation_purpose
    ):
        raise ValueError("Method-step execution-start workflow integrity failed")


@dataclass(frozen=True, slots=True, weakref_slot=True, eq=False)
class MethodStepOutcome:
    """Immutable occurrence evidence; live transition authority stays external."""

    workflow_identity: str
    workflow_binding_identity: str
    occurrence_identity: str
    starting_state_identity: str
    lifecycle: MethodStepLifecycle
    primary_strategy: MethodStepStrategy | None = None
    primary_terminal: str | None = None
    primary_execution_identity: str | None = None
    evaluation_result_identity: str | None = None
    evaluation_failure_identity: str | None = None
    accepted_result_identity: str | None = None
    commit_operation_identity: str | None = None
    derivations: tuple[DerivationOutcome, ...] = ()
    publication_identity: str | None = None
    publication_failure: str = ""
    run_provenance_identity: str | None = None
    successor_state_identity: str | None = None
    evaluation_result: EvaluationResult | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    evaluation_failure: EvaluationFailure | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    primary_execution: MethodStepPrimaryExecution | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    accepted_result: AcceptedFitResult | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    commit_operation: FitCommitOperation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    publication: PublishedStepReference | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    run_provenance: NativeRunInformation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    source_workflow: MethodStepWorkflow | None = field(
        default=None,
        repr=False,
        compare=False,
        kw_only=True,
    )
    record_identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.evaluation_result is not None and (
            self.evaluation_result_identity != self.evaluation_result.identity
        ):
            raise ValueError("Method-step evaluation identity differs from its result")
        if self.evaluation_failure is not None and (
            self.evaluation_failure_identity != self.evaluation_failure.identity
        ):
            raise ValueError("Method-step evaluation identity differs from its failure")
        if self.accepted_result is not None and (
            self.accepted_result_identity != self.accepted_result.identity
        ):
            raise ValueError("Method-step accepted identity differs from its result")
        if self.commit_operation is not None and (
            self.commit_operation_identity != self.commit_operation.identity
        ):
            raise ValueError("Method-step commit identity differs from its operation")
        if self.publication is not None and (
            self.publication_identity != self.publication.identity
        ):
            raise ValueError("Method-step publication identity differs from its source")
        if self.run_provenance is not None and (
            self.run_provenance_identity
            != self.run_provenance.execution_occurrence_identity
        ):
            raise ValueError("Method-step run provenance differs from its source")
        object.__setattr__(
            self,
            "record_identity",
            _identity("native-method-step-outcome-v1", self._record_payload()),
        )
        _validate_method_step_outcome(self)

    def validate_integrity(self) -> None:
        """Rederive this outcome and every available live child."""
        _validate_method_step_outcome(self)

    def _record_payload(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "workflow_identity": self.workflow_identity,
            "workflow_binding_identity": self.workflow_binding_identity,
            "occurrence_identity": self.occurrence_identity,
            "starting_state_identity": self.starting_state_identity,
            "lifecycle": self.lifecycle.value,
            "primary_strategy": (
                None if self.primary_strategy is None else self.primary_strategy.value
            ),
            "primary_terminal": self.primary_terminal,
            "primary_execution_identity": self.primary_execution_identity,
            "evaluation_result_identity": self.evaluation_result_identity,
            "evaluation_failure_identity": self.evaluation_failure_identity,
            "accepted_result_identity": self.accepted_result_identity,
            "commit_operation_identity": self.commit_operation_identity,
            "derivations": [
                {
                    "identity": item.identity,
                    "request_identity": item.request_identity,
                    "stage": item.stage,
                    "disposition": item.disposition.value,
                    "operation_identity": item.operation_identity,
                    "artifact_identities": list(item.artifact_identities),
                    "message": item.message,
                }
                for item in self.derivations
            ],
            "publication_identity": self.publication_identity,
            "publication_failure": self.publication_failure,
            "run_provenance_identity": self.run_provenance_identity,
            "successor_state_identity": self.successor_state_identity,
        }

    def to_record(self) -> dict[str, object]:
        """Serialize historical evidence without any process-local authority."""
        self.validate_integrity()
        return {**self._record_payload(), "identity": self.record_identity}

    @classmethod
    def from_record(cls, record: object) -> MethodStepOutcome:
        """Restore historical evidence without recreating successor authority."""
        if not isinstance(record, dict):
            raise TypeError("Method-step outcome record must be a mapping")
        typed_record = cast("dict[str, object]", record)
        required = {
            "schema_version",
            "workflow_identity",
            "workflow_binding_identity",
            "occurrence_identity",
            "starting_state_identity",
            "lifecycle",
            "primary_strategy",
            "primary_terminal",
            "primary_execution_identity",
            "evaluation_result_identity",
            "evaluation_failure_identity",
            "accepted_result_identity",
            "commit_operation_identity",
            "derivations",
            "publication_identity",
            "publication_failure",
            "run_provenance_identity",
            "successor_state_identity",
            "identity",
        }
        if set(typed_record) != required or typed_record.get("schema_version") != 1:
            raise ValueError("Malformed method-step outcome record")
        if not isinstance(typed_record.get("derivations"), list):
            raise TypeError("Historical method-step derivations must be a list")
        derivation_records = cast("list[object]", typed_record["derivations"])
        derivations: list[DerivationOutcome] = []
        for raw in derivation_records:
            if not isinstance(raw, dict):
                raise TypeError("Historical derivation outcome must be a mapping")
            item = cast("dict[str, object]", raw)
            if set(item) != {
                "identity",
                "request_identity",
                "stage",
                "disposition",
                "operation_identity",
                "artifact_identities",
                "message",
            } or not isinstance(item["artifact_identities"], list):
                raise ValueError("Malformed historical derivation outcome")
            derivation = DerivationOutcome(
                str(item["request_identity"]),
                str(item["stage"]),
                DerivationDisposition(str(item["disposition"])),
                cast("str | None", item["operation_identity"]),
                tuple(str(value) for value in item["artifact_identities"]),
                str(item["message"]),
            )
            if item["identity"] != derivation.identity:
                raise ValueError("Historical derivation identity differs from payload")
            derivations.append(derivation)
        outcome = cls(
            str(typed_record["workflow_identity"]),
            str(typed_record["workflow_binding_identity"]),
            str(typed_record["occurrence_identity"]),
            str(typed_record["starting_state_identity"]),
            MethodStepLifecycle(str(typed_record["lifecycle"])),
            (
                None
                if typed_record["primary_strategy"] is None
                else MethodStepStrategy(str(typed_record["primary_strategy"]))
            ),
            cast("str | None", typed_record["primary_terminal"]),
            cast("str | None", typed_record["primary_execution_identity"]),
            cast("str | None", typed_record["evaluation_result_identity"]),
            cast("str | None", typed_record["evaluation_failure_identity"]),
            cast("str | None", typed_record["accepted_result_identity"]),
            cast("str | None", typed_record["commit_operation_identity"]),
            tuple(derivations),
            cast("str | None", typed_record["publication_identity"]),
            str(typed_record["publication_failure"]),
            cast("str | None", typed_record["run_provenance_identity"]),
            cast("str | None", typed_record["successor_state_identity"]),
        )
        if typed_record["identity"] != outcome.record_identity:
            raise ValueError("Method-step outcome identity differs from its payload")
        outcome.validate_integrity()
        return outcome


def _validate_derivation_outcome_integrity(outcome: DerivationOutcome) -> None:
    expected = DerivationOutcome(
        outcome.request_identity,
        outcome.stage,
        outcome.disposition,
        outcome.operation_identity,
        outcome.artifact_identities,
        outcome.message,
        operation=outcome.operation,
        artifacts=outcome.artifacts,
    )
    if expected.identity != outcome.identity:
        raise ValueError("Method-step derivation outcome integrity validation failed")


def _validate_published_step_integrity(
    reference: PublishedStepReference,
    *,
    workflow: MethodStepWorkflow | None = None,
) -> None:
    _validate_recursive_dataclass_children(reference.provenance)
    if workflow is not None:
        reference.provenance.validate_method_step_context(
            parameterization=workflow.parameterization,
            plan=workflow.engine.plan,
            method=workflow.method,
            semantic_workflow_identity=workflow.semantic_identity,
            grouping_topology=reference.provenance.grouping_topology,
            execution=reference.provenance.execution,
        )
    expected_identity = published_step_reference_identity(
        publication_occurrence_identity=reference.publication_occurrence_identity,
        lifecycle=reference.lifecycle,
        authority=reference.authority,
        manifest_identity=reference.manifest_identity,
        manifest_sha256=reference.manifest_sha256,
        provenance=reference.provenance,
        independent_ids=reference.independent_ids,
        artifacts=reference.artifacts,
    )
    if expected_identity != reference.identity:
        raise ValueError("Method-step publication integrity failed")
    reference.require_exact_live_publication()


def _validate_outcome_lifecycle(outcome: MethodStepOutcome) -> None:
    lifecycle = outcome.lifecycle
    has_evaluation = outcome.evaluation_result_identity is not None
    has_evaluation_failure = outcome.evaluation_failure_identity is not None
    has_primary = outcome.primary_execution_identity is not None
    has_accepted = outcome.accepted_result_identity is not None
    has_commit = outcome.commit_operation_identity is not None
    has_successor = outcome.successor_state_identity is not None
    if has_evaluation and has_evaluation_failure:
        raise ValueError(
            "Method-step outcome integrity has conflicting evaluation truth"
        )
    if lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE:
        valid = (
            has_evaluation
            and has_successor
            and not any((has_evaluation_failure, has_primary, has_accepted, has_commit))
        )
    elif lifecycle is MethodStepLifecycle.NO_OBJECTIVE_DATA:
        valid = not any(
            (
                has_evaluation,
                has_evaluation_failure,
                has_primary,
                has_accepted,
                has_commit,
                has_successor,
            )
        )
    elif lifecycle is MethodStepLifecycle.COMMITTED:
        valid = (
            has_primary
            and has_accepted
            and has_commit
            and has_successor
            and outcome.primary_terminal == "accepted"
            and not has_evaluation
            and not has_evaluation_failure
        )
    elif lifecycle is MethodStepLifecycle.ACCEPTED_UNCOMMITTED:
        valid = (
            has_primary
            and has_accepted
            and has_commit
            and not has_successor
            and not has_evaluation
            and not has_evaluation_failure
        )
    elif lifecycle is MethodStepLifecycle.PUBLICATION_FAILED:
        committed_failure = (
            has_primary
            and has_accepted
            and has_commit
            and has_successor
            and outcome.primary_terminal == "accepted"
        )
        evaluation_failure = has_evaluation and not any(
            (has_primary, has_accepted, has_commit, has_successor)
        )
        valid = bool(outcome.publication_failure) and (
            committed_failure or evaluation_failure
        )
    else:
        valid = not has_successor and not has_commit
        if lifecycle is MethodStepLifecycle.FAILED:
            valid = valid and (has_evaluation_failure or has_primary)
        elif lifecycle in {
            MethodStepLifecycle.CANCELLED,
            MethodStepLifecycle.INTERRUPTED,
        }:
            valid = valid and outcome.primary_terminal == lifecycle.value
    if not valid:
        raise ValueError("Method-step outcome lifecycle integrity validation failed")


def _validate_method_step_outcome(  # noqa: C901 - closed recursive outcome tree
    outcome: MethodStepOutcome,
) -> None:
    """Canonically rederive one outcome without trusting its stored outer hash."""
    for derivation in outcome.derivations:
        _validate_derivation_outcome_integrity(derivation)
    if outcome.record_identity != _identity(
        "native-method-step-outcome-v1",
        outcome._record_payload(),
    ):
        raise ValueError("Method-step outcome integrity validation failed")
    _validate_outcome_lifecycle(outcome)
    workflow = outcome.source_workflow
    if workflow is None:
        return
    workflow.validate_integrity()
    if (
        outcome.workflow_identity != workflow.semantic_identity
        or outcome.workflow_binding_identity != workflow.binding_identity
        or outcome.starting_state_identity
        != _snapshot_identity(workflow.starting_snapshot)
    ):
        raise ValueError("Method-step outcome source-workflow integrity failed")
    if outcome.evaluation_result is not None:
        if (
            outcome.evaluation_result.identity != outcome.evaluation_result_identity
            or outcome.evaluation_result.plan_identity != workflow.engine.plan.identity
            or outcome.evaluation_result.parameterization_identity
            != workflow.parameterization.evaluator_identity
        ):
            raise ValueError("Method-step evaluation-result integrity failed")
        type(outcome.evaluation_result).from_record(
            outcome.evaluation_result.to_record(),
            workflow.engine.plan,
        )
    if outcome.evaluation_failure is not None and (
        outcome.evaluation_failure.identity != outcome.evaluation_failure_identity
    ):
        raise ValueError("Method-step evaluation-failure integrity failed")
    if outcome.primary_execution is not None:
        _validate_recursive_dataclass_children(outcome.primary_execution)
        expected_primary_identity = _primary_execution_identity(
            outcome.primary_execution
        )
        if expected_primary_identity != outcome.primary_execution_identity:
            raise ValueError("Method-step primary-outcome integrity failed")
        replace(outcome.primary_execution)
    if outcome.accepted_result is not None:
        accepted = outcome.accepted_result
        if accepted.final_residual_jacobian is not None:
            accepted.final_residual_jacobian.validate_integrity()
        reconstructed = replace(accepted, occurrence_witness=None)
        if (
            reconstructed.identity != accepted.identity
            or accepted.identity != outcome.accepted_result_identity
            or not accepted_occurrence_is_authoritative(accepted)
            or accepted.source_occurrence_identity
            != workflow.starting_snapshot.occurrence_identity
            or accepted.source_revision != workflow.starting_snapshot.revision
        ):
            raise ValueError("Method-step accepted-result integrity failed")
        if (
            outcome.primary_execution is None
            or outcome.primary_execution.accepted_result is not accepted
        ):
            raise ValueError("Method-step aggregate acceptance integrity failed")
    if outcome.commit_operation is not None:
        operation = outcome.commit_operation
        if operation.receipt is not None:
            _validate_rederived_identity(operation.receipt)
        if operation.failure is not None:
            _validate_rederived_identity(operation.failure)
        if operation.committed_snapshot is not None:
            _validate_snapshot_integrity(operation.committed_snapshot)
        reconstructed = replace(operation, _occurrence_witness=None)
        if (
            reconstructed.identity != operation.identity
            or operation.identity != outcome.commit_operation_identity
            or not fit_commit_operation_is_authoritative(operation)
            or outcome.accepted_result is None
            or operation.accepted_result_identity != outcome.accepted_result.identity
            or operation.accepted_occurrence_identity
            != outcome.accepted_result.occurrence_identity
        ):
            raise ValueError("Method-step commit-operation integrity failed")
    if outcome.publication is not None:
        if outcome.publication.identity != outcome.publication_identity:
            raise ValueError("Method-step publication integrity failed")
        _validate_published_step_integrity(outcome.publication, workflow=workflow)
    if outcome.run_provenance is not None:
        run = outcome.run_provenance
        _validate_rederived_identity(run.inputs)
        _validate_rederived_identity(run.parameter_model)
        _validate_snapshot_integrity(run.starting_snapshot)
        for step in run.steps:
            _validate_published_step_integrity(step)
        reconstructed = replace(run)
        if (
            reconstructed.execution_occurrence_identity
            != run.execution_occurrence_identity
            or run.execution_occurrence_identity != outcome.run_provenance_identity
            or run.parameter_model.identity != workflow.parameter_model.identity
            or (
                outcome.publication is not None
                and run.steps[-1] is not outcome.publication
            )
        ):
            raise ValueError("Method-step run-provenance integrity failed")
    if outcome.successor_state_identity is not None:
        expected_successor = (
            outcome.commit_operation.committed_snapshot
            if outcome.commit_operation is not None
            else workflow.starting_snapshot
        )
        if (
            expected_successor is None
            or _snapshot_identity(expected_successor)
            != outcome.successor_state_identity
        ):
            raise ValueError("Method-step successor-state integrity failed")


@dataclass(frozen=True, slots=True)
class _SuccessorAuthorityBinding:
    outcome: weakref.ReferenceType[MethodStepOutcome]
    workflow: MethodStepWorkflow
    snapshot: AnalysisValuesSnapshot
    outcome_record_identity: str
    workflow_binding_identity: str
    snapshot_identity: str
    occurrence_identity: str
    revision: int


_SUCCESSOR_AUTHORITIES: dict[
    int,
    _SuccessorAuthorityBinding,
] = {}
_PUBLICATION_SOURCES: dict[
    int,
    tuple[weakref.ReferenceType[PublishedStepReference], object],
] = {}


def _grant_successor_authority(
    outcome: MethodStepOutcome,
    snapshot: AnalysisValuesSnapshot,
) -> None:
    outcome.validate_integrity()
    workflow = outcome.source_workflow
    if workflow is None or outcome.lifecycle not in {
        MethodStepLifecycle.COMMITTED,
        MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE,
    }:
        raise ValueError("Method-step outcome cannot mint successor authority")
    _validate_snapshot_integrity(snapshot)
    key = id(outcome)

    def remove(_reference: weakref.ReferenceType[MethodStepOutcome]) -> None:
        _SUCCESSOR_AUTHORITIES.pop(key, None)

    _SUCCESSOR_AUTHORITIES[key] = _SuccessorAuthorityBinding(
        weakref.ref(outcome, remove),
        workflow,
        snapshot,
        outcome.record_identity,
        workflow.binding_identity,
        _snapshot_identity(snapshot),
        snapshot.occurrence_identity,
        snapshot.revision,
    )


def _retain_publication_source(
    reference: PublishedStepReference,
    source: object,
) -> None:
    key = id(reference)

    def remove(_reference: weakref.ReferenceType[PublishedStepReference]) -> None:
        _PUBLICATION_SOURCES.pop(key, None)

    _PUBLICATION_SOURCES[key] = (weakref.ref(reference, remove), source)


def require_successor_state(
    outcome: MethodStepOutcome,
    analysis_values: AnalysisValues,
) -> AnalysisValuesSnapshot:
    """Return the exact live successor or reject historical/failed evidence."""
    try:
        outcome.validate_integrity()
    except (TypeError, ValueError) as error:
        raise ValueError(
            "Method-step outcome integrity blocks successor authority"
        ) from error
    binding = _SUCCESSOR_AUTHORITIES.get(id(outcome))
    current = analysis_values.snapshot()
    if (
        binding is None
        or binding.outcome() is not outcome
        or outcome.source_workflow is not binding.workflow
        or binding.outcome_record_identity != outcome.record_identity
        or binding.workflow_binding_identity != binding.workflow.binding_identity
        or outcome.lifecycle
        not in {
            MethodStepLifecycle.COMMITTED,
            MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE,
        }
        or _snapshot_identity(binding.snapshot) != binding.snapshot_identity
        or binding.snapshot.occurrence_identity != binding.occurrence_identity
        or binding.snapshot.revision != binding.revision
        or current != binding.snapshot
        or current.occurrence_identity != binding.occurrence_identity
        or current.revision != binding.revision
    ):
        raise ValueError("Method-step outcome has no live successor authority")
    return binding.snapshot


def execute_method_step(  # noqa: C901 - closed evaluation/optimization lifecycle
    workflow: MethodStepWorkflow,
    *,
    analysis_values: AnalysisValues,
    cancellation: CancellationToken | None = None,
    checkpoint_observer: Callable[[MethodStepCheckpoint], None] | None = None,
    commit_completed_observer: (
        Callable[[AcceptedFitResult, FitCommitOperation], None] | None
    ) = None,
    progress_observer: ContextualProgressObserver | None = None,
) -> MethodStepOutcome:
    """Execute one exact native workflow occurrence."""
    execution_start = _pin_execution_start_workflow(workflow)
    current = analysis_values.snapshot()
    if current != workflow.starting_snapshot:
        raise ValueError("Method-step workflow no longer has its exact starting state")
    starting_state_identity = _snapshot_identity(workflow.starting_snapshot)
    if workflow.evaluation_purpose is None:
        return _execute_optimization(
            workflow,
            execution_start=execution_start,
            analysis_values=analysis_values,
            cancellation=cancellation,
            checkpoint_observer=checkpoint_observer,
            commit_completed_observer=commit_completed_observer,
            progress_observer=progress_observer,
        )
    if cancellation is not None and cancellation.is_cancelled:
        return MethodStepOutcome(
            workflow.semantic_identity,
            workflow.binding_identity,
            uuid4().hex,
            starting_state_identity,
            MethodStepLifecycle.CANCELLED,
            primary_terminal="cancelled",
            source_workflow=workflow,
        )
    if workflow.evaluation_purpose is EvaluationPurpose.NO_OBJECTIVE_DATA:
        publication: PublishedStepReference | None = None
        run_provenance: NativeRunInformation | None = None
        publication_failure = ""
        if workflow.publication is not None:
            try:
                publication = _publish_no_objective_method_step(workflow)
                run_provenance = _publish_run_provenance(workflow, publication)
            except Exception as error:  # noqa: BLE001 - typed outcome stays truthful
                publication_failure = f"{type(error).__name__}: {error}"
        return MethodStepOutcome(
            workflow.semantic_identity,
            workflow.binding_identity,
            uuid4().hex,
            starting_state_identity,
            MethodStepLifecycle.NO_OBJECTIVE_DATA,
            publication_identity=(
                None if publication is None else publication.identity
            ),
            publication_failure=publication_failure,
            run_provenance_identity=(
                None
                if run_provenance is None
                else run_provenance.execution_occurrence_identity
            ),
            publication=publication,
            run_provenance=run_provenance,
            source_workflow=workflow,
        )
    frame = EvaluationFrame.from_lifecycle_frame(
        workflow.parameterization,
        workflow.parameterization.frame_from_snapshot(workflow.starting_snapshot),
    )
    result = workflow.engine.new_evaluator().evaluate(frame)
    if isinstance(result, EvaluationFailure):
        publication: PublishedStepReference | None = None
        run_provenance: NativeRunInformation | None = None
        publication_failure = ""
        if workflow.publication is not None:
            try:
                publication = _publish_evaluation_failure_method_step(workflow, result)
                run_provenance = _publish_run_provenance(workflow, publication)
            except Exception as error:  # noqa: BLE001 - failure evidence remains truthful
                publication_failure = f"{type(error).__name__}: {error}"
        return MethodStepOutcome(
            workflow.semantic_identity,
            workflow.binding_identity,
            uuid4().hex,
            starting_state_identity,
            MethodStepLifecycle.FAILED,
            evaluation_failure_identity=result.identity,
            publication_identity=(
                None if publication is None else publication.identity
            ),
            publication_failure=publication_failure,
            run_provenance_identity=(
                None
                if run_provenance is None
                else run_provenance.execution_occurrence_identity
            ),
            evaluation_failure=result,
            publication=publication,
            run_provenance=run_provenance,
            source_workflow=workflow,
        )
    publication: PublishedStepReference | None = None
    run_provenance: NativeRunInformation | None = None
    publication_failure = ""
    lifecycle = MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE
    if workflow.publication is not None:
        try:
            publication = _publish_evaluation_method_step(workflow, result)
            run_provenance = _publish_run_provenance(workflow, publication)
        except Exception as error:  # noqa: BLE001 - evaluated evidence remains truthful
            lifecycle = MethodStepLifecycle.PUBLICATION_FAILED
            publication_failure = f"{type(error).__name__}: {error}"
    outcome = MethodStepOutcome(
        workflow.semantic_identity,
        workflow.binding_identity,
        uuid4().hex,
        starting_state_identity,
        lifecycle,
        evaluation_result_identity=result.identity,
        publication_identity=(None if publication is None else publication.identity),
        publication_failure=publication_failure,
        run_provenance_identity=(
            None
            if run_provenance is None
            else run_provenance.execution_occurrence_identity
        ),
        successor_state_identity=(
            starting_state_identity
            if lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE
            else None
        ),
        evaluation_result=result,
        publication=publication,
        run_provenance=run_provenance,
        source_workflow=workflow,
    )
    if lifecycle is MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE:
        _grant_successor_authority(outcome, workflow.starting_snapshot)
    return outcome


def _primary_execution_identity(execution: MethodStepPrimaryExecution) -> str:
    if isinstance(execution, GroupedGridDirectTrfOutcome):
        return _identity(
            "native-method-step-grouped-grid-execution",
            (
                execution.terminal.value,
                tuple(item.identity for item in execution.attempts),
                None
                if execution.accepted_result is None
                else execution.accepted_result.identity,
                None if execution.failure is None else execution.failure.identity,
            ),
        )
    if isinstance(execution, GroupedDirectTrfOutcome):
        return _identity(
            "native-method-step-grouped-direct-execution",
            (
                execution.terminal.value,
                tuple(item.identity for item in execution.components),
                None
                if execution.accepted_result is None
                else execution.accepted_result.identity,
                None if execution.failure is None else execution.failure.identity,
            ),
        )
    if isinstance(execution, GridDirectTrfOutcome):
        return _identity(
            "native-method-step-grid-execution",
            (
                execution.terminal.value,
                tuple(item.identity for item in execution.attempts),
                None
                if execution.accepted_result is None
                else execution.accepted_result.identity,
                None if execution.failure is None else execution.failure.identity,
            ),
        )
    return _identity(
        "native-method-step-de-execution",
        (
            execution.terminal.value,
            execution.search.identity,
            execution.accounting.identity,
            None
            if execution.accepted_result is None
            else execution.accepted_result.identity,
            None if execution.failure is None else execution.failure.identity,
        ),
    )


def _validate_fit_commit_boundary(
    workflow: MethodStepWorkflow,
    expected: _ExecutionStartWorkflowContext,
    execution: MethodStepPrimaryExecution,
    primary_identity: str,
    accepted: AcceptedFitResult,
    authority: LiveFitCommitAuthority,
    analysis_values: AnalysisValues,
) -> None:
    """Prove exact workflow, acceptance, and authority lineage before commit."""
    _validate_execution_start_workflow(workflow, expected)
    problem = expected.problem
    invocation = expected.invocation
    strategy = expected.strategy
    current = analysis_values.snapshot()
    if (
        problem is None
        or invocation is None
        or strategy is None
        or current != expected.starting_snapshot
        or current.occurrence_identity != expected.starting_snapshot.occurrence_identity
        or current.revision != expected.starting_snapshot.revision
        or problem is not workflow.problem
        or problem.source_snapshot is not workflow.starting_snapshot
        or _primary_execution_identity(execution) != primary_identity
        or execution.accepted_result is not accepted
        or execution.commit_authority is not authority
        or not accepted_occurrence_is_authoritative(accepted)
        or accepted.problem_identity != problem.identity
        or accepted.parameterization_identity != expected.parameterization.identity
        or accepted.evaluator_parameterization_identity
        != expected.parameterization.evaluator_identity
        or accepted.source_occurrence_identity
        != expected.starting_snapshot.occurrence_identity
        or accepted.source_revision != expected.starting_snapshot.revision
        or accepted.controlled_ids != problem.controlled_ids
        or accepted.commit_scope != problem.commit_scope
        or accepted.evaluation_result.plan_identity != expected.engine.plan.identity
    ):
        raise ValueError("Method-step fit-commit boundary integrity failed")
    reconstructed = replace(accepted, occurrence_witness=None)
    if reconstructed.identity != accepted.identity:
        raise ValueError("Method-step accepted-result integrity failed before commit")
    if strategy is MethodStepStrategy.DIRECT_TRF:
        if accepted.invocation_identity != invocation.identity:
            raise ValueError(
                "Method-step Direct invocation lineage changed before commit"
            )
    elif strategy is MethodStepStrategy.GRID_DIRECT_TRF:
        typed = cast("AcceptedGridDirectTrfResult", accepted)
        if typed.workflow_invocation is not invocation:
            raise ValueError(
                "Method-step GRID invocation lineage changed before commit"
            )
        validate_grid_commit_lineage(typed, problem)
    else:
        typed_de = cast("AcceptedDeDirectTrfResult", accepted)
        if typed_de.workflow_invocation is not invocation:
            raise ValueError("Method-step DE invocation lineage changed before commit")
        validate_de_commit_lineage(typed_de, problem)
    if not fit_commit_authority_is_authoritative(accepted, authority, problem):
        raise ValueError("Method-step fit-commit authority pairing is no longer live")


def _precommit_failure_outcome(
    expected: _ExecutionStartWorkflowContext,
    primary_identity: str,
    accepted: AcceptedFitResult,
    *,
    interrupted: bool = False,
) -> MethodStepOutcome:
    """Retain non-authoritative acceptance diagnostics after a fenced stop."""
    disposition = (
        DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP
        if interrupted
        else DerivationDisposition.BLOCKED_BY_PREREQUISITE
    )
    derivations = tuple(
        DerivationOutcome(
            item.identity,
            _derivation_stage(item),
            disposition,
            message="Fit commit boundary validation did not permit transition",
        )
        for item in expected.derivations
    )
    return MethodStepOutcome(
        expected.semantic_identity,
        expected.binding_identity,
        uuid4().hex,
        _snapshot_identity(expected.starting_snapshot),
        (
            MethodStepLifecycle.INTERRUPTED
            if interrupted
            else MethodStepLifecycle.FAILED
        ),
        expected.strategy,
        "interrupted" if interrupted else "accepted",
        primary_identity,
        accepted_result_identity=accepted.identity,
        derivations=derivations,
    )


def _execute_optimization(  # noqa: C901 - closed primary transition
    workflow: MethodStepWorkflow,
    *,
    execution_start: _ExecutionStartWorkflowContext,
    analysis_values: AnalysisValues,
    cancellation: CancellationToken | None,
    checkpoint_observer: Callable[[MethodStepCheckpoint], None] | None,
    commit_completed_observer: (
        Callable[[AcceptedFitResult, FitCommitOperation], None] | None
    ),
    progress_observer: ContextualProgressObserver | None,
) -> MethodStepOutcome:
    _validate_execution_start_workflow(workflow, execution_start)
    problem = cast("OptimizationProblem", workflow.problem)
    decomposition = cast("FitDecomposition", workflow.decomposition)
    strategy = cast("MethodStepStrategy", workflow.strategy)
    invocation = cast("MethodStepInvocation", workflow.invocation)
    token = CancellationToken() if cancellation is None else cancellation
    try:
        if strategy is MethodStepStrategy.DIRECT_TRF:
            execution: MethodStepPrimaryExecution = execute_grouped_direct_trf(
                problem,
                decomposition,
                cast("GroupedDirectTrfInvocation", invocation),
                workflow.parameterization,
                workflow.engine,
                cancellation=token,
                progress_observer=progress_observer,
            )
        elif strategy is MethodStepStrategy.GRID_DIRECT_TRF:
            execution = execute_grouped_grid_direct_trf(
                problem,
                decomposition,
                cast("GridDirectTrfInvocation", invocation),
                workflow.parameterization,
                workflow.engine,
                cancellation=token,
                progress_observer=progress_observer,
            )
        else:
            execution = execute_de_direct_trf(
                problem,
                cast("DeDirectTrfInvocation", invocation),
                workflow.parameterization,
                workflow.engine,
                cancellation=token,
            )
    except GridDirectTrfInterrupted as error:
        execution = error.outcome
    except DeDirectTrfInterrupted as error:
        execution = error.outcome
    primary_identity = _primary_execution_identity(execution)
    terminal = execution.terminal
    accepted = execution.accepted_result
    authority = execution.commit_authority
    if accepted is None or authority is None:
        if terminal.value == "cancelled":
            lifecycle = MethodStepLifecycle.CANCELLED
        elif terminal.value == "interrupted":
            lifecycle = MethodStepLifecycle.INTERRUPTED
        else:
            lifecycle = MethodStepLifecycle.FAILED
        derivations = (
            _stopped_derivations(workflow.derivations)
            if lifecycle
            in {MethodStepLifecycle.CANCELLED, MethodStepLifecycle.INTERRUPTED}
            else tuple(
                DerivationOutcome(
                    item.identity,
                    _derivation_stage(item),
                    DerivationDisposition.BLOCKED_BY_PREREQUISITE,
                    message="Successful aggregate acceptance prerequisite was unavailable",
                )
                for item in workflow.derivations
            )
        )
        publication: PublishedStepReference | None = None
        run_provenance: NativeRunInformation | None = None
        publication_failure = ""
        if workflow.publication is not None:
            try:
                publication = _publish_suppressed_method_step(
                    workflow,
                    execution,
                    primary_identity,
                    lifecycle,
                )
                run_provenance = _publish_run_provenance(workflow, publication)
            except Exception as error:  # noqa: BLE001 - primary evidence remains truthful
                publication_failure = f"{type(error).__name__}: {error}"
        return MethodStepOutcome(
            workflow.semantic_identity,
            workflow.binding_identity,
            uuid4().hex,
            _snapshot_identity(workflow.starting_snapshot),
            lifecycle,
            strategy,
            terminal.value,
            primary_identity,
            derivations=derivations,
            publication_identity=(
                None if publication is None else publication.identity
            ),
            publication_failure=publication_failure,
            run_provenance_identity=(
                None
                if run_provenance is None
                else run_provenance.execution_occurrence_identity
            ),
            primary_execution=execution,
            publication=publication,
            run_provenance=run_provenance,
            source_workflow=workflow,
        )
    if checkpoint_observer is not None:
        try:
            checkpoint_observer(MethodStepCheckpoint.AGGREGATE_ACCEPTED)
        except KeyboardInterrupt:
            return _precommit_failure_outcome(
                execution_start,
                primary_identity,
                accepted,
                interrupted=True,
            )
        except Exception:  # noqa: BLE001 - observer failure is a typed stop
            return _precommit_failure_outcome(
                execution_start,
                primary_identity,
                accepted,
            )
    try:
        _validate_fit_commit_boundary(
            workflow,
            execution_start,
            execution,
            primary_identity,
            accepted,
            authority,
            analysis_values,
        )
    except Exception:  # noqa: BLE001 - integrity failure is a typed stop
        return _precommit_failure_outcome(
            execution_start,
            primary_identity,
            accepted,
        )
    if token.is_cancelled:
        derivations = _stopped_derivations(workflow.derivations)
        publication: PublishedStepReference | None = None
        run_provenance: NativeRunInformation | None = None
        publication_failure = ""
        if workflow.publication is not None:
            try:
                publication = _publish_suppressed_method_step(
                    workflow,
                    execution,
                    primary_identity,
                    MethodStepLifecycle.CANCELLED,
                    accepted=accepted,
                    terminal="cancelled",
                )
                run_provenance = _publish_run_provenance(workflow, publication)
            except Exception as error:  # noqa: BLE001 - accepted evidence stays truthful
                publication_failure = f"{type(error).__name__}: {error}"
        return MethodStepOutcome(
            workflow.semantic_identity,
            workflow.binding_identity,
            uuid4().hex,
            _snapshot_identity(workflow.starting_snapshot),
            MethodStepLifecycle.CANCELLED,
            strategy,
            "cancelled",
            primary_identity,
            accepted_result_identity=accepted.identity,
            derivations=derivations,
            publication_identity=(
                None if publication is None else publication.identity
            ),
            publication_failure=publication_failure,
            run_provenance_identity=(
                None
                if run_provenance is None
                else run_provenance.execution_occurrence_identity
            ),
            primary_execution=execution,
            accepted_result=accepted,
            publication=publication,
            run_provenance=run_provenance,
            source_workflow=workflow,
        )
    operation = execute_fit_commit(
        accepted,
        authority,
        problem=problem,
        parameterization=workflow.parameterization,
        analysis_values=analysis_values,
    )
    if checkpoint_observer is not None:
        checkpoint_observer(MethodStepCheckpoint.COMMIT_COMPLETED)
        _validate_execution_start_workflow(workflow, execution_start)
    if commit_completed_observer is not None:
        try:
            commit_completed_observer(accepted, operation)
        except KeyboardInterrupt:
            token.cancel()
        except Exception:  # noqa: BLE001, S110 - reporting is non-scientific
            pass
    committed = operation.terminal is FitCommitTerminal.COMMITTED
    successor = operation.committed_snapshot if committed else None
    derivations = (
        _execute_derivations(workflow, accepted, token)
        if committed
        else tuple(
            DerivationOutcome(
                item.identity,
                _derivation_stage(item),
                DerivationDisposition.BLOCKED_BY_PREREQUISITE,
                message="Successful fit commit prerequisite was unavailable",
            )
            for item in workflow.derivations
        )
    )
    publication: PublishedStepReference | None = None
    run_provenance: NativeRunInformation | None = None
    publication_failure = ""
    lifecycle = (
        MethodStepLifecycle.COMMITTED
        if committed
        else MethodStepLifecycle.ACCEPTED_UNCOMMITTED
    )
    if workflow.publication is not None:
        try:
            if committed:
                publication = _publish_committed_method_step(
                    workflow,
                    execution,
                    primary_identity,
                    accepted,
                    operation,
                    derivations,
                    analysis_values,
                )
            else:
                publication = _publish_suppressed_method_step(
                    workflow,
                    execution,
                    primary_identity,
                    MethodStepLifecycle.ACCEPTED_UNCOMMITTED,
                    operation=operation,
                    accepted=accepted,
                )
            run_provenance = _publish_run_provenance(workflow, publication)
        except Exception as error:  # noqa: BLE001 - commit remains authoritative
            if committed:
                lifecycle = MethodStepLifecycle.PUBLICATION_FAILED
            publication_failure = f"{type(error).__name__}: {error}"
    outcome = MethodStepOutcome(
        workflow.semantic_identity,
        workflow.binding_identity,
        uuid4().hex,
        _snapshot_identity(workflow.starting_snapshot),
        lifecycle,
        strategy,
        terminal.value,
        primary_identity,
        accepted_result_identity=accepted.identity,
        commit_operation_identity=operation.identity,
        derivations=derivations,
        publication_identity=(None if publication is None else publication.identity),
        publication_failure=publication_failure,
        run_provenance_identity=(
            None
            if run_provenance is None
            else run_provenance.execution_occurrence_identity
        ),
        successor_state_identity=(
            None if successor is None else _snapshot_identity(successor)
        ),
        primary_execution=execution,
        accepted_result=accepted,
        commit_operation=operation,
        publication=publication,
        run_provenance=run_provenance,
        source_workflow=workflow,
    )
    if successor is not None and lifecycle is MethodStepLifecycle.COMMITTED:
        _grant_successor_authority(outcome, successor)
    return outcome


def _method_step_primary_record(
    workflow: MethodStepWorkflow,
    execution: MethodStepPrimaryExecution,
    aggregate_execution_identity: str,
    accepted: AcceptedFitResult | None,
    *,
    terminal: MethodStepPrimaryTerminal | None = None,
) -> MethodStepPrimaryRecord:
    decomposition = cast("FitDecomposition", workflow.decomposition)
    invocation = cast("MethodStepInvocation", workflow.invocation)
    grouping = tuple(
        (component.identity, component.controlled_ids)
        for component in decomposition.components
    )
    if isinstance(invocation, GroupedDirectTrfInvocation):
        limit = sum(
            item.objective_request_budget for item in invocation.component_invocations
        )
        settings = invocation.component_invocations[0].execution_settings
        used = None
        seeds: tuple[SeedRecord, ...] = ()
    elif isinstance(invocation, GridDirectTrfInvocation):
        limit = invocation.objective_request_budget * len(invocation.seeds)
        settings = ExecutionSettings()
        used = None
        seeds = ()
    else:
        de_execution = cast("DeDirectTrfOutcome", execution)
        limit = (
            invocation.de_objective_request_budget
            + invocation.polish_objective_request_budget
        )
        settings = ExecutionSettings()
        used = de_execution.accounting.de_counters.objective_requests_accepted + (
            0
            if de_execution.accounting.polish_counters is None
            else de_execution.accounting.polish_counters.objective_requests_accepted
        )
        seeds = (
            SeedRecord("primary_de_root", invocation.root_seed, invocation.identity),
        )
    strategy = cast("MethodStepStrategy", workflow.strategy)
    return MethodStepPrimaryRecord(
        workflow.semantic_identity,
        cast("OptimizationProblem", workflow.problem).identity,
        invocation.identity if accepted is None else accepted.invocation_identity,
        (
            aggregate_execution_identity
            if accepted is None
            else accepted.execution_identity
        ),
        aggregate_execution_identity,
        execution.terminal.value if terminal is None else terminal,
        grouping,
        settings,
        PolicyRecord(f"primary_{strategy.value}", invocation.identity),
        BudgetRecord(f"primary_{strategy.value}_objective_requests", limit, used),
        seeds,
    )


def _components_for_publication(
    workflow: MethodStepWorkflow,
    execution: MethodStepPrimaryExecution,
) -> tuple[ComponentDiagnostic, ...]:
    decomposition = cast("FitDecomposition", workflow.decomposition)
    if isinstance(execution, GroupedDirectTrfOutcome):
        return tuple(
            ComponentDiagnostic(
                item.identity,
                item.disposition.value,
                item.controlled_ids,
                None if item.candidate is None else item.candidate.chi_square,
            )
            for item in execution.components
        )
    if isinstance(execution, GroupedGridDirectTrfOutcome):
        selected = None
        if execution.selection is not None:
            selected = execution.selection.selected_record
        components = (
            tuple(
                component
                for attempt in execution.attempts
                for component in attempt.components
            )
            if selected is None
            else cast("GroupedGridSeedOutcome", selected).components
        )
        return tuple(
            ComponentDiagnostic(
                item.identity,
                item.disposition.value,
                item.controlled_ids,
                None if item.candidate is None else item.candidate.chi_square,
            )
            for item in components
        )
    terminal = execution.terminal.value
    disposition = (
        "succeeded"
        if terminal == "accepted"
        else terminal
        if terminal in {"cancelled", "interrupted"}
        else "failed"
    )
    return tuple(
        ComponentDiagnostic(
            component.identity,
            disposition,
            component.controlled_ids,
        )
        for component in decomposition.components
    )


def _publication_evidence(
    workflow: MethodStepWorkflow,
    derivations: tuple[DerivationOutcome, ...],
) -> tuple[
    UncertaintyEvidence | None,
    tuple[ResamplingPublication, ...],
    McmcPublication | None,
    tuple[ResamplingEvidence, ...],
    tuple[ResamplingSummaryFailurePublication, ...],
    McmcEvidence | None,
]:
    uncertainty: UncertaintyEvidence | None = None
    resampling: list[ResamplingPublication] = []
    mcmc: McmcPublication | None = None
    partial_resampling: list[ResamplingEvidence] = []
    resampling_summary_failures: list[ResamplingSummaryFailurePublication] = []
    partial_mcmc: McmcEvidence | None = None
    requests = {item.identity: item for item in workflow.derivations}
    for outcome in derivations:
        request = requests[outcome.request_identity]
        if isinstance(request, UncertaintyDerivationRequest):
            if outcome.artifacts:
                uncertainty = cast("UncertaintyEvidence", outcome.artifacts[0])
        elif isinstance(request, ResamplingDerivationRequest):
            if not outcome.artifacts:
                continue
            evidence = cast("ResamplingEvidence", outcome.artifacts[0])
            if outcome.disposition is DerivationDisposition.COMPLETED:
                summary = summarize_resampling_evidence(
                    evidence,
                    request.summary_policy,
                )
                resampling.append(
                    ResamplingPublication(
                        evidence,
                        summary,
                    )
                )
            elif (
                evidence.lifecycle.value == "completed"
                and len(outcome.artifacts) > 1
                and isinstance(outcome.artifacts[1], SummaryFailure)
            ):
                resampling_summary_failures.append(
                    ResamplingSummaryFailurePublication(
                        evidence,
                        outcome.artifacts[1],
                    )
                )
            else:
                partial_resampling.append(evidence)
        else:
            if not outcome.artifacts:
                continue
            evidence = cast("McmcEvidence", outcome.artifacts[0])
            if outcome.disposition is DerivationDisposition.COMPLETED:
                posterior = cast("PosteriorSampleEvidence", outcome.artifacts[1])
                summary = cast("PosteriorSummary", outcome.artifacts[2])
                mcmc = McmcPublication(evidence, posterior, summary)
            else:
                partial_mcmc = evidence
    return (
        uncertainty,
        tuple(resampling),
        mcmc,
        tuple(partial_resampling),
        tuple(resampling_summary_failures),
        partial_mcmc,
    )


def _method_step_provenance(
    workflow: MethodStepWorkflow,
    primary: MethodStepPrimaryRecord | None,
    *,
    policies: tuple[PolicyRecord, ...] = (),
    budgets: tuple[BudgetRecord, ...] = (),
    seeds: tuple[SeedRecord, ...] = (),
) -> WorkflowProvenance:
    request = cast("MethodStepPublicationRequest", workflow.publication)
    grouping = (
        (("aggregate", workflow.parameterization.independent_ids),)
        if primary is None
        else primary.grouping_topology
    )
    execution = ExecutionSettings() if primary is None else primary.execution_settings
    return WorkflowProvenance.create_method_step(
        parameterization=workflow.parameterization,
        plan=workflow.engine.plan,
        method=workflow.method,
        semantic_workflow_identity=workflow.semantic_identity,
        grouping_topology=grouping,
        policies=policies,
        budgets=budgets,
        seeds=seeds,
        execution=execution,
        environment=request.environment,
        baseline_references=request.baseline_references,
    )


def _publish_evaluation_method_step(
    workflow: MethodStepWorkflow,
    result: EvaluationResult,
) -> PublishedStepReference:
    workflow.validate_integrity()
    request = cast("MethodStepPublicationRequest", workflow.publication)
    provenance = _method_step_provenance(workflow, None)
    publication = EvaluationPublication(
        workflow.engine.plan,
        workflow.method,
        workflow.parameterization,
        result,
        provenance=provenance,
        method_step_semantic_identity=workflow.semantic_identity,
        allow_controlled=(
            workflow.evaluation_purpose is EvaluationPurpose.EVALUATE_ONLY
        ),
    )
    reference = publish_native_results(
        request.path,
        publication,
    )
    _retain_publication_source(reference, publication)
    return reference


def _publish_evaluation_failure_method_step(
    workflow: MethodStepWorkflow,
    failure: EvaluationFailure,
) -> PublishedStepReference:
    workflow.validate_integrity()
    request = cast("MethodStepPublicationRequest", workflow.publication)
    publication = SuppressedPublication(
        "failed",
        failure,
        workflow.engine.plan,
        workflow.method,
        workflow.parameterization,
        provenance=_method_step_provenance(workflow, None),
        method_step_semantic_identity=workflow.semantic_identity,
    )
    reference = publish_native_results(request.path, publication)
    _retain_publication_source(reference, publication)
    return reference


def _publish_no_objective_method_step(
    workflow: MethodStepWorkflow,
) -> PublishedStepReference:
    workflow.validate_integrity()
    request = cast("MethodStepPublicationRequest", workflow.publication)
    publication = NoObjectivePublication(
        workflow.engine.plan,
        workflow.method,
        workflow.parameterization,
        workflow.semantic_identity,
        provenance=_method_step_provenance(workflow, None),
    )
    reference = publish_native_results(request.path, publication)
    _retain_publication_source(reference, publication)
    return reference


def _publish_suppressed_method_step(
    workflow: MethodStepWorkflow,
    execution: MethodStepPrimaryExecution,
    aggregate_execution_identity: str,
    lifecycle: MethodStepLifecycle,
    *,
    operation: FitCommitOperation | None = None,
    accepted: AcceptedFitResult | None = None,
    terminal: MethodStepPrimaryTerminal | None = None,
) -> PublishedStepReference:
    workflow.validate_integrity()
    request = cast("MethodStepPublicationRequest", workflow.publication)
    primary = _method_step_primary_record(
        workflow,
        execution,
        aggregate_execution_identity,
        accepted,
        terminal=terminal,
    )
    provenance = _method_step_provenance(
        workflow,
        primary,
        policies=(primary.policy,),
        budgets=(primary.budget,),
        seeds=primary.seeds,
    )
    suppressed_lifecycle = cast(
        "Literal['failed', 'accepted_uncommitted', 'cancelled', 'interrupted']",
        lifecycle.value,
    )
    publication = SuppressedPublication(
        suppressed_lifecycle,
        primary if operation is None else operation,
        workflow.engine.plan,
        workflow.method,
        workflow.parameterization,
        provenance=provenance,
        primary_invocation=primary,
        primary_execution=primary,
        primary_problem=cast("OptimizationProblem", workflow.problem),
        accepted_result_identity=(None if accepted is None else accepted.identity),
        accepted_occurrence_identity=(
            None if accepted is None else accepted.occurrence_identity
        ),
        components=_components_for_publication(workflow, execution),
    )
    reference = publish_native_results(
        request.path,
        publication,
    )
    _retain_publication_source(reference, publication)
    return reference


def _publish_run_provenance(
    workflow: MethodStepWorkflow,
    publication: PublishedStepReference,
) -> NativeRunInformation | None:
    workflow.validate_integrity()
    publication_request = cast("MethodStepPublicationRequest", workflow.publication)
    request = publication_request.run_provenance
    if request is None:
        return None
    run = NativeRunInformation(
        request.invocation_identity,
        request.inputs,
        workflow.parameter_model,
        request.starting_snapshot,
        (*request.prior_steps, publication),
    )
    write_native_run_info(
        Namespace(output=request.output_directory),
        run,
        argv=request.argv,
        working_directory=request.working_directory,
    )
    return run


def _publish_committed_method_step(
    workflow: MethodStepWorkflow,
    execution: MethodStepPrimaryExecution,
    aggregate_execution_identity: str,
    accepted: AcceptedFitResult,
    operation: FitCommitOperation,
    derivations: tuple[DerivationOutcome, ...],
    analysis_values: AnalysisValues,
) -> PublishedStepReference:
    workflow.validate_integrity()
    request = cast("MethodStepPublicationRequest", workflow.publication)
    if operation.receipt is None or operation.committed_snapshot is None:
        raise ValueError("Committed publication requires exact commit artifacts")
    primary = _method_step_primary_record(
        workflow,
        execution,
        aggregate_execution_identity,
        accepted,
    )
    (
        uncertainty,
        resampling,
        mcmc,
        partial_resampling,
        resampling_summary_failures,
        partial_mcmc,
    ) = _publication_evidence(workflow, derivations)
    policies, budgets, seeds = publication_provenance_records(
        primary,
        uncertainty=uncertainty,
        resampling=(
            *(item.evidence for item in resampling),
            *partial_resampling,
            *(item.evidence for item in resampling_summary_failures),
        ),
        mcmc=(partial_mcmc if mcmc is None else mcmc.evidence),
    )
    provenance = WorkflowProvenance.create_method_step(
        parameterization=workflow.parameterization,
        plan=workflow.engine.plan,
        method=workflow.method,
        semantic_workflow_identity=workflow.semantic_identity,
        grouping_topology=primary.grouping_topology,
        policies=policies,
        budgets=budgets,
        seeds=seeds,
        execution=primary.execution_settings,
        environment=request.environment,
        baseline_references=request.baseline_references,
    )
    publication = CommittedFitPublication(
        workflow.engine.plan,
        workflow.method,
        workflow.parameterization,
        workflow.parameter_model,
        workflow.starting_snapshot,
        cast("OptimizationProblem", workflow.problem),
        primary,
        primary,
        accepted,
        operation.receipt,
        operation.committed_snapshot,
        operation,
        analysis_values,
        provenance=provenance,
        components=_components_for_publication(workflow, execution),
        uncertainty=uncertainty,
        resampling=resampling,
        mcmc=mcmc,
        partial_resampling=partial_resampling,
        resampling_summary_failures=resampling_summary_failures,
        partial_mcmc=partial_mcmc,
    )
    reference = publish_native_results(
        request.path,
        publication,
    )
    _retain_publication_source(reference, publication)
    return reference


def _derivation_stage(request: MethodStepDerivationRequest) -> str:
    if isinstance(request, UncertaintyDerivationRequest):
        return "uncertainty"
    if isinstance(request, ResamplingDerivationRequest):
        return "resampling"
    return "mcmc"


def _stopped_derivations(
    requests: tuple[MethodStepDerivationRequest, ...],
) -> tuple[DerivationOutcome, ...]:
    return tuple(
        DerivationOutcome(
            item.identity,
            _derivation_stage(item),
            DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP,
            message="Workflow stop was observed before the requested stage started",
        )
        for item in requests
    )


def _execute_derivations(
    workflow: MethodStepWorkflow,
    accepted: AcceptedFitResult,
    cancellation: CancellationToken,
) -> tuple[DerivationOutcome, ...]:
    workflow.validate_integrity()
    outcomes: list[DerivationOutcome] = []
    for index, request in enumerate(workflow.derivations):
        workflow.validate_integrity()
        if cancellation.is_cancelled:
            outcomes.extend(_stopped_derivations(workflow.derivations[index:]))
            break
        try:
            if isinstance(request, UncertaintyDerivationRequest):
                outcome = _execute_uncertainty(
                    workflow, accepted, request, cancellation
                )
            elif isinstance(request, ResamplingDerivationRequest):
                outcome = _execute_resampling(workflow, accepted, request, cancellation)
            else:
                outcome = _execute_mcmc(workflow, accepted, request, cancellation)
        except KeyboardInterrupt:
            outcomes.append(
                DerivationOutcome(
                    request.identity,
                    _derivation_stage(request),
                    DerivationDisposition.INTERRUPTED,
                    message="Derivation was interrupted",
                )
            )
            outcomes.extend(_stopped_derivations(workflow.derivations[index + 1 :]))
            break
        except Exception as error:  # noqa: BLE001 - branch-local typed evidence
            outcomes.append(
                DerivationOutcome(
                    request.identity,
                    _derivation_stage(request),
                    DerivationDisposition.FAILED,
                    message=f"{type(error).__name__}: {error}",
                )
            )
            continue
        outcomes.append(outcome)
        if outcome.disposition in {
            DerivationDisposition.CANCELLED,
            DerivationDisposition.INTERRUPTED,
        }:
            outcomes.extend(_stopped_derivations(workflow.derivations[index + 1 :]))
            break
    return tuple(outcomes)


def _execute_uncertainty(
    workflow: MethodStepWorkflow,
    accepted: AcceptedFitResult,
    request: UncertaintyDerivationRequest,
    cancellation: CancellationToken,
) -> DerivationOutcome:
    workflow.validate_integrity()
    problem = cast("OptimizationProblem", workflow.problem)

    def cancellation_probe() -> UncertaintyOperationTerminal | None:
        return (
            UncertaintyOperationTerminal.CANCELLED
            if cancellation.is_cancelled
            else None
        )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=workflow.parameterization,
        engine=workflow.engine,
        policy=request.policy,
        constrained_scope=request.constrained_scope,
        constrained_units=request.constrained_units,
        constrained_scales=request.constrained_scales,
        compiled_constraint_linearization=request.compiled_capabilities,
        cancellation_probe=cancellation_probe,
        resolved_environment_identity=request.resolved_environment_identity,
    )
    terminals = {item.terminal for item in evidence.operations}
    disposition = (
        DerivationDisposition.CANCELLED
        if UncertaintyOperationTerminal.CANCELLED in terminals
        else DerivationDisposition.INTERRUPTED
        if UncertaintyOperationTerminal.INTERRUPTED in terminals
        else DerivationDisposition.FAILED
        if evidence.failures
        else DerivationDisposition.COMPLETED
    )
    return DerivationOutcome(
        request.identity,
        "uncertainty",
        disposition,
        artifact_identities=(evidence.identity,),
        artifacts=(evidence,),
    )


def _execute_resampling(
    workflow: MethodStepWorkflow,
    accepted: AcceptedFitResult,
    request: ResamplingDerivationRequest,
    cancellation: CancellationToken,
) -> DerivationOutcome:
    workflow.validate_integrity()
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=ResamplingDatasetManifest(
            workflow.engine.plan,
            tuple(
                float(value)
                for value in accepted.evaluation_result.normalized_calculations
            ),
            request.references,
            request.nucleus_groups,
            request.observation_descriptors,
        ),
        source_problem=cast("OptimizationProblem", workflow.problem),
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        scheme=request.scheme,
        replicate_count=request.replicate_count,
        replicate_structural_identities=request.replicate_structural_identities,
        replicate_component_identities=request.replicate_component_identities,
        root_seed=request.root_seed,
        output_scope=request.output_scope,
        output_units=request.output_units,
        minimum_successful_count=request.minimum_successful_count,
        strategy=request.strategy,
        strategy_settings=request.strategy_settings,
    )
    operation = execute_resampling_evidence(
        accepted,
        plan,
        cancellation_probe=lambda: cancellation.is_cancelled,
    )
    artifacts: list[_IdentifiedArtifact] = []
    if operation.evidence is not None:
        artifacts.append(operation.evidence)
        summary = summarize_resampling_evidence(
            operation.evidence,
            request.summary_policy,
        )
        summary_payload = (
            summary.summary if summary.summary is not None else summary.failure
        )
        if summary_payload is not None:
            artifacts.append(summary_payload)
        summary_failed = summary.terminal.value != "completed"
    else:
        summary_failed = False
    disposition = (
        DerivationDisposition.CANCELLED
        if operation.terminal is ResamplingOperationTerminal.CANCELLED
        else DerivationDisposition.INTERRUPTED
        if operation.terminal is ResamplingOperationTerminal.INTERRUPTED
        else DerivationDisposition.FAILED
        if summary_failed
        else DerivationDisposition.COMPLETED
    )
    return DerivationOutcome(
        request.identity,
        "resampling",
        disposition,
        operation.identity,
        tuple(item.identity for item in artifacts),
        operation=operation,
        artifacts=tuple(artifacts),
    )


def _execute_mcmc(
    workflow: MethodStepWorkflow,
    accepted: AcceptedFitResult,
    request: McmcDerivationRequest,
    cancellation: CancellationToken,
) -> DerivationOutcome:
    workflow.validate_integrity()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=cast("OptimizationProblem", workflow.problem),
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        policy=request.policy,
        coordinate_units=request.coordinate_units,
    )
    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
    )
    artifacts: list[_IdentifiedArtifact] = []
    if operation.evidence is not None:
        artifacts.append(operation.evidence)
    if operation.terminal is McmcOperationTerminal.COMPLETED:
        evidence = cast("McmcEvidence", operation.evidence)
        try:
            selection = derive_retained_sample_view(
                evidence,
                cancellation=cancellation,
            )
            posterior: PosteriorSampleEvidence = derive_posterior_sample_evidence(
                selection,
                request.output_units,
                cancellation=cancellation,
            )
            summary: PosteriorSummary = derive_posterior_summary(
                posterior,
                cancellation=cancellation,
            )
        except KeyboardInterrupt:
            return DerivationOutcome(
                request.identity,
                "mcmc",
                DerivationDisposition.INTERRUPTED,
                operation.identity,
                (evidence.identity,),
                "KeyboardInterrupt: MCMC posterior derivation was interrupted",
                operation=operation,
                artifacts=(evidence,),
            )
        except Exception as error:  # noqa: BLE001 - retain completed primary evidence
            return DerivationOutcome(
                request.identity,
                "mcmc",
                (
                    DerivationDisposition.CANCELLED
                    if cancellation.is_cancelled
                    else DerivationDisposition.FAILED
                ),
                operation.identity,
                (evidence.identity,),
                f"{type(error).__name__}: {error}",
                operation=operation,
                artifacts=(evidence,),
            )
        artifacts.extend((posterior, summary))
        disposition = DerivationDisposition.COMPLETED
    elif operation.terminal is McmcOperationTerminal.CANCELLED:
        disposition = DerivationDisposition.CANCELLED
    elif operation.terminal is McmcOperationTerminal.INTERRUPTED:
        disposition = DerivationDisposition.INTERRUPTED
    else:
        disposition = DerivationDisposition.FAILED
    return DerivationOutcome(
        request.identity,
        "mcmc",
        disposition,
        operation.identity,
        tuple(item.identity for item in artifacts),
        operation=operation,
        artifacts=tuple(artifacts),
    )


__all__ = [
    "DerivationDisposition",
    "DerivationOutcome",
    "EvaluationPurpose",
    "McmcDerivationRequest",
    "MethodStepCheckpoint",
    "MethodStepLifecycle",
    "MethodStepOutcome",
    "MethodStepPublicationRequest",
    "MethodStepRunProvenanceRequest",
    "MethodStepStrategy",
    "MethodStepWorkflow",
    "ResamplingDerivationRequest",
    "UncertaintyDerivationRequest",
    "execute_method_step",
    "require_successor_state",
]
