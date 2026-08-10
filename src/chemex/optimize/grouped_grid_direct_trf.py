"""Exact fit-component GRID to native Direct TRF qualification (#596).

This isolated qualification path composes the grouped Direct TRF contract with
the Cartesian GRID contract.  It is intentionally not wired into production
dispatch.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import cast
from uuid import uuid4

from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CandidateMaterialization,
    CandidateSummary,
    CommitReceipt,
    DirectTrfConstructionError,
    DirectTrfInvocation,
    GridSeedProblemDerivation,
    LiveFitCommitAuthority,
    MaterializationTerminal,
    MaterializedDirectTrfCandidate,
    OptimizationProblem,
    TerminalFailure,
    _accept_materialized_fit_for_derived_workflow,
    _grant_derived_fit_commit_authority,
    canonical_chi_square,
    commit_accepted_fit,
)
from chemex.optimize.grid_direct_trf import (
    GridCandidate,
    GridCoordinate,
    GridDirectTrfInvocation,
    GridSeed,
    GridSeedDisposition,
    _validate_grid_context,
)
from chemex.optimize.grouped_direct_trf import (
    ComponentDisposition,
    FitComponentOutcome,
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfOutcome,
    GroupedDirectTrfTerminal,
    _reconstruct_fresh_aggregate,
    _root_projection_matches,
    _validate_component_projections,
    execute_direct_trf_components,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.values import AnalysisValues

_SCHEMA_VERSION = 1
_WORKFLOW_VERSION = "native-grouped-cartesian-grid-direct-trf-v1"
_AGGREGATE_ORDER_VERSION = "aggregate-chi-square-vector-seed-ordinal-v1"


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class GroupedGridSeedOutcome:
    """One seed's component evidence and fresh non-authoritative aggregate."""

    seed_identity: str
    seed_ordinal: int
    axis_items: tuple[tuple[str, GridCoordinate], ...]
    start: tuple[GridCoordinate, ...]
    disposition: GridSeedDisposition
    components: tuple[FitComponentOutcome, ...]
    seed_decomposition: FitDecomposition | None = field(
        default=None,
        repr=False,
        compare=False,
        kw_only=True,
    )
    component_invocation: GroupedDirectTrfInvocation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    execution_identity: str | None = None
    objective: float | None = None
    candidate: GridCandidate | None = field(default=None, repr=False, compare=False)
    failure: TerminalFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        eligible = self.disposition is GridSeedDisposition.ELIGIBLE
        if eligible != (self.candidate is not None and self.objective is not None):
            raise ValueError("Only an eligible grouped GRID seed has an aggregate")
        if eligible and (
            self.seed_decomposition is None
            or self.component_invocation is None
            or self.execution_identity is None
            or self.candidate is None
            or self.candidate.seed_identity != self.seed_identity
            or self.candidate.seed_ordinal != self.seed_ordinal
            or self.candidate.objective != self.objective
            or self.candidate.candidate.execution_identity != self.execution_identity
            or any(
                component.disposition is not ComponentDisposition.SUCCEEDED
                for component in self.components
            )
        ):
            raise ValueError(
                "Eligible grouped GRID seed lacks its exact aggregate evidence"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grouped-grid-seed-outcome",
                (
                    self.seed_identity,
                    self.seed_ordinal,
                    self.disposition.value,
                    tuple(component.identity for component in self.components),
                    None
                    if self.seed_decomposition is None
                    else self.seed_decomposition.identity,
                    None
                    if self.component_invocation is None
                    else self.component_invocation.identity,
                    self.execution_identity,
                    None if self.candidate is None else self.candidate.identity,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class GroupedGridSelection:
    """Canonical ordering over fresh whole-problem aggregate candidates."""

    selected_seed_identity: str
    selected_seed_ordinal: int
    selected_candidate_identity: str
    eligible_candidate_identities: tuple[str, ...]
    candidate_records: tuple[GroupedGridSeedOutcome, ...] = field(
        repr=False,
        compare=False,
    )
    ordering_policy: str = _AGGREGATE_ORDER_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        selected = tuple(
            record
            for record in self.candidate_records
            if record.seed_identity == self.selected_seed_identity
            and record.seed_ordinal == self.selected_seed_ordinal
            and record.candidate is not None
            and record.candidate.identity == self.selected_candidate_identity
        )
        candidate_identities = tuple(
            cast("GridCandidate", record.candidate).identity
            for record in self.candidate_records
        )
        if (
            self.ordering_policy != _AGGREGATE_ORDER_VERSION
            or not self.candidate_records
            or any(
                record.disposition is not GridSeedDisposition.ELIGIBLE
                or record.candidate is None
                for record in self.candidate_records
            )
            or candidate_identities != self.eligible_candidate_identities
            or self.eligible_candidate_identities[0] != self.selected_candidate_identity
            or len(selected) != 1
            or selected[0].identity != self.candidate_records[0].identity
        ):
            raise ValueError(
                "Grouped GRID selection requires one canonical aggregate winner"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grouped-grid-selection",
                (
                    self.ordering_policy,
                    self.selected_seed_identity,
                    self.selected_seed_ordinal,
                    self.selected_candidate_identity,
                    self.eligible_candidate_identities,
                ),
            ),
        )

    @classmethod
    def from_candidate_records(
        cls,
        records: Sequence[GroupedGridSeedOutcome],
    ) -> GroupedGridSelection:
        ordered = tuple(
            sorted(
                records,
                key=lambda record: cast(
                    "GridCandidate", record.candidate
                ).ordering_key(),
            )
        )
        if not ordered or any(record.candidate is None for record in ordered):
            raise ValueError("Grouped GRID selection requires aggregate candidates")
        selected = ordered[0]
        candidate = cast("GridCandidate", selected.candidate)
        return cls(
            selected.seed_identity,
            selected.seed_ordinal,
            candidate.identity,
            tuple(
                cast("GridCandidate", record.candidate).identity for record in ordered
            ),
            ordered,
        )

    @property
    def selected_record(self) -> GroupedGridSeedOutcome:
        return self.candidate_records[0]


@dataclass(frozen=True, slots=True)
class GroupedGridSelectionProvenance:
    """Complete authority lineage from one aggregate GRID winner."""

    workflow_invocation_identity: str
    root_problem_identity: str
    decomposition_identity: str
    selection_identity: str
    selected_seed_identity: str
    selected_seed_ordinal: int
    selected_outcome_identity: str
    aggregate_candidate_identity: str
    aggregate_materialization_identity: str
    accepted_materialization_identity: str
    accepted_evaluation_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grouped-grid-selection-provenance",
                (
                    _WORKFLOW_VERSION,
                    self.workflow_invocation_identity,
                    self.root_problem_identity,
                    self.decomposition_identity,
                    self.selection_identity,
                    self.selected_seed_identity,
                    self.selected_seed_ordinal,
                    self.selected_outcome_identity,
                    self.aggregate_candidate_identity,
                    self.aggregate_materialization_identity,
                    self.accepted_materialization_identity,
                    self.accepted_evaluation_identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class AcceptedGroupedGridDirectTrfResult(AcceptedFitResult):
    """Accepted root evidence carrying exact grouped-GRID selection lineage."""

    grouped_grid_provenance: GroupedGridSelectionProvenance
    workflow_invocation: GridDirectTrfInvocation = field(repr=False, compare=False)
    decomposition: FitDecomposition = field(repr=False, compare=False)
    selection: GroupedGridSelection = field(repr=False, compare=False)
    selected_seed: GridSeed = field(repr=False, compare=False)
    selected_outcome: GroupedGridSeedOutcome = field(repr=False, compare=False)
    fresh_candidate: MaterializedDirectTrfCandidate = field(
        repr=False,
        compare=False,
    )

    @classmethod
    def from_accepted(
        cls,
        accepted: AcceptedFitResult,
        provenance: GroupedGridSelectionProvenance,
        workflow_invocation: GridDirectTrfInvocation,
        decomposition: FitDecomposition,
        selection: GroupedGridSelection,
        selected_seed: GridSeed,
        selected_outcome: GroupedGridSeedOutcome,
        fresh_candidate: MaterializedDirectTrfCandidate,
    ) -> AcceptedGroupedGridDirectTrfResult:
        return cls(
            accepted.occurrence_identity,
            accepted.problem_identity,
            accepted.invocation_identity,
            accepted.execution_identity,
            accepted.materialization_identity,
            accepted.parameterization_identity,
            accepted.evaluator_parameterization_identity,
            accepted.source_occurrence_identity,
            accepted.source_revision,
            accepted.controlled_ids,
            accepted.vector,
            accepted.chi_square,
            accepted.evaluation_result,
            accepted.commit_scope,
            accepted.commit_items,
            accepted.origin_context_identity,
            provenance,
            workflow_invocation,
            decomposition,
            selection,
            selected_seed,
            selected_outcome,
            fresh_candidate,
        )


class GroupedGridDirectTrfTerminal(StrEnum):
    """Closed outcomes of one grouped GRID to Direct TRF occurrence."""

    ACCEPTED = "accepted"
    NO_ELIGIBLE_CANDIDATE = "no_eligible_candidate"
    EXECUTION_FAILURE = "execution_failure"
    MATERIALIZATION_FAILURE = "materialization_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


@dataclass(frozen=True, slots=True)
class GroupedGridDirectTrfOutcome:
    """Complete occurrence; only the selected root aggregate has authority."""

    terminal: GroupedGridDirectTrfTerminal
    attempts: tuple[GroupedGridSeedOutcome, ...]
    selection: GroupedGridSelection | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    accepted_result: AcceptedGroupedGridDirectTrfResult | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    commit_authority: LiveFitCommitAuthority | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    failure: TerminalFailure | None = None

    def __post_init__(self) -> None:
        accepted = self.terminal is GroupedGridDirectTrfTerminal.ACCEPTED
        if accepted != (
            self.accepted_result is not None and self.commit_authority is not None
        ):
            raise ValueError("Only accepted grouped GRID execution has authority")
        if not accepted and (
            self.accepted_result is not None or self.commit_authority is not None
        ):
            raise ValueError("Unaccepted grouped GRID execution cannot have authority")
        if accepted and self.selection is None:
            raise ValueError("Accepted grouped GRID execution requires a selection")


def _validate_context(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> None:
    _validate_grid_context(problem, invocation, parameterization, engine)
    expected = FitDecomposition.from_root(problem, parameterization, engine)
    if (
        not problem.acceptance_authority
        or decomposition.root_problem_identity != problem.identity
        or decomposition.identity != expected.identity
        or invocation.root_problem_identity != problem.identity
    ):
        raise DirectTrfConstructionError(
            "Grouped GRID context is not rooted in one complete problem"
        )


def _component_invocation(
    root_problem: OptimizationProblem,
    seed_decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
) -> GroupedDirectTrfInvocation:
    root_scales = dict(
        zip(root_problem.controlled_ids, invocation.x_scale, strict=True)
    )
    component_invocations = tuple(
        DirectTrfInvocation.for_problem(
            component.problem,
            objective_request_budget=invocation.objective_request_budget,
            x_scale=tuple(
                root_scales[param_id] for param_id in component.controlled_ids
            ),
            ftol=invocation.ftol,
            xtol=invocation.xtol,
            gtol=invocation.gtol,
        )
        for component in seed_decomposition.components
    )
    return GroupedDirectTrfInvocation(
        seed_decomposition.root_problem_identity,
        seed_decomposition.identity,
        component_invocations,
    )


def _fresh_root_candidate(
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    vector: tuple[float, ...],
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> MaterializedDirectTrfCandidate:
    lifecycle = problem.lifecycle_frame(vector, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    evaluator = engine.new_evaluator()
    evaluated = evaluator.evaluate(frame)
    if isinstance(evaluated, EvaluationFailure):
        raise DirectTrfConstructionError(
            f"Grouped GRID root materialization failed: {evaluated.message}"
        )
    return _materialized_root_candidate(
        problem,
        invocation_identity,
        execution_identity,
        vector,
        evaluator,
        evaluated,
    )


def _materialized_root_candidate(
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    vector: tuple[float, ...],
    evaluator: BoundEvaluator,
    evaluated: EvaluationResult,
) -> MaterializedDirectTrfCandidate:
    """Record one already evaluated root vector as fresh candidate evidence."""
    chi_square = canonical_chi_square(evaluated.residuals)
    statistics = evaluator.cache_statistics
    summary = CandidateSummary(vector, chi_square, 0)
    materialization = CandidateMaterialization(
        uuid4().hex,
        problem.identity,
        invocation_identity,
        execution_identity,
        summary,
        MaterializationTerminal.SUCCESS,
        1,
        evaluator.compatibility_identity,
        evaluated.identity,
        statistics.hits,
        statistics.misses,
    )
    return MaterializedDirectTrfCandidate(
        problem.identity,
        invocation_identity,
        execution_identity,
        materialization,
        vector,
        chi_square,
        evaluated,
    )


def _failed_component_seed_outcome(
    seed: GridSeed,
    seed_decomposition: FitDecomposition,
    components: tuple[FitComponentOutcome, ...],
    component_invocation: GroupedDirectTrfInvocation,
) -> GroupedGridSeedOutcome:
    dispositions = tuple(component.disposition for component in components)
    if ComponentDisposition.CANCELLED in dispositions:
        disposition = GridSeedDisposition.CANCELLED
    elif ComponentDisposition.INTERRUPTED in dispositions:
        disposition = GridSeedDisposition.INTERRUPTED
    elif ComponentDisposition.EXECUTION_FAILURE in dispositions:
        disposition = GridSeedDisposition.IMPLEMENTATION_FAILURE
    else:
        disposition = GridSeedDisposition.NON_CONVERGED
    failure = next(
        (component.failure for component in components if component.failure),
        TerminalFailure("component_failure", "A grouped GRID component failed"),
    )
    return GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        disposition,
        components,
        component_invocation,
        failure=failure,
        seed_decomposition=seed_decomposition,
    )


def _failed_aggregate_seed_outcome(
    seed: GridSeed,
    seed_decomposition: FitDecomposition,
    components: tuple[FitComponentOutcome, ...],
    component_invocation: GroupedDirectTrfInvocation,
    outcome: GroupedDirectTrfOutcome,
) -> GroupedGridSeedOutcome:
    disposition = {
        GroupedDirectTrfTerminal.CANCELLED: GridSeedDisposition.CANCELLED,
        GroupedDirectTrfTerminal.INTERRUPTED: GridSeedDisposition.INTERRUPTED,
        GroupedDirectTrfTerminal.EXECUTION_FAILURE: (
            GridSeedDisposition.IMPLEMENTATION_FAILURE
        ),
    }.get(outcome.terminal, GridSeedDisposition.MATERIALIZATION_FAILURE)
    return GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        disposition,
        components,
        component_invocation,
        failure=outcome.failure,
        seed_decomposition=seed_decomposition,
    )


def _execute_seed(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
    seed: GridSeed,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    token: CancellationToken,
) -> GroupedGridSeedOutcome:
    if seed.rejection is not None:
        disposition = (
            GridSeedDisposition.UNUSABLE_VALUE
            if seed.rejection.category == "grid_seed_unusable_value"
            else GridSeedDisposition.OUT_OF_BOUNDS
        )
        return GroupedGridSeedOutcome(
            seed.identity,
            seed.ordinal,
            seed.axis_items,
            seed.start,
            disposition,
            (),
            failure=seed.rejection,
        )
    if seed.problem is None:
        raise DirectTrfConstructionError("Feasible grouped GRID seed lacks a problem")
    seed_decomposition = FitDecomposition.from_root(
        seed.problem,
        parameterization,
        engine,
        cancellation=token,
    )
    if (
        seed_decomposition.identity != decomposition.identity
        or tuple(component.identity for component in seed_decomposition.components)
        != tuple(component.identity for component in decomposition.components)
        or seed.problem.source_snapshot is not problem.source_snapshot
    ):
        raise DirectTrfConstructionError(
            "Grouped GRID seed differs from the exact root decomposition"
        )
    component_invocation = _component_invocation(
        problem,
        seed_decomposition,
        invocation,
    )
    components = execute_direct_trf_components(
        seed_decomposition,
        component_invocation,
        parameterization,
        engine,
        cancellation=token,
    )
    if any(
        component.disposition is not ComponentDisposition.SUCCEEDED
        for component in components
    ):
        return _failed_component_seed_outcome(
            seed,
            seed_decomposition,
            components,
            component_invocation,
        )
    projections = _validate_component_projections(
        seed_decomposition,
        component_invocation,
        engine,
        components,
        token,
    )
    if isinstance(projections, GroupedDirectTrfOutcome):
        return _failed_aggregate_seed_outcome(
            seed,
            seed_decomposition,
            components,
            component_invocation,
            projections,
        )
    fresh = _reconstruct_fresh_aggregate(
        problem,
        seed_decomposition,
        parameterization,
        engine,
        components,
        token,
    )
    if isinstance(fresh, GroupedDirectTrfOutcome):
        return _failed_aggregate_seed_outcome(
            seed,
            seed_decomposition,
            components,
            component_invocation,
            fresh,
        )
    vector, evaluator, aggregate = fresh
    if any(
        not _root_projection_matches(aggregate, engine, component, projection)
        for component, projection in zip(
            seed_decomposition.components,
            projections,
            strict=True,
        )
    ):
        return GroupedGridSeedOutcome(
            seed.identity,
            seed.ordinal,
            seed.axis_items,
            seed.start,
            GridSeedDisposition.MATERIALIZATION_FAILURE,
            components,
            component_invocation,
            failure=TerminalFailure(
                "decomposition_projection_mismatch",
                "Fresh grouped GRID aggregate differs from component evidence",
            ),
            seed_decomposition=seed_decomposition,
        )
    execution_identity = _identity(
        "native-grouped-grid-seed-execution",
        (
            invocation.identity,
            seed.identity,
            component_invocation.identity,
            tuple(component.identity for component in components),
            aggregate.identity,
        ),
    )
    materialized = _materialized_root_candidate(
        problem,
        invocation.identity,
        execution_identity,
        vector,
        evaluator,
        aggregate,
    )
    candidate = GridCandidate(seed.identity, seed.ordinal, materialized)
    return GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.ELIGIBLE,
        components,
        component_invocation,
        execution_identity,
        candidate.objective,
        candidate,
        seed_decomposition=seed_decomposition,
    )


def _not_started_seed(
    seed: GridSeed,
    decomposition: FitDecomposition,
) -> GroupedGridSeedOutcome:
    components = tuple(
        FitComponentOutcome(
            component.identity,
            component.controlled_ids,
            ComponentDisposition.NOT_STARTED,
        )
        for component in decomposition.components
    )
    return GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.NOT_STARTED,
        components,
    )


def _seed_exception(
    seed: GridSeed,
    decomposition: FitDecomposition,
    error: BaseException,
    *,
    interrupted: bool = False,
) -> GroupedGridSeedOutcome:
    components = tuple(
        FitComponentOutcome(
            component.identity,
            component.controlled_ids,
            ComponentDisposition.NOT_STARTED,
        )
        for component in decomposition.components
    )
    return GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        (
            GridSeedDisposition.INTERRUPTED
            if interrupted
            else GridSeedDisposition.IMPLEMENTATION_FAILURE
        ),
        components,
        failure=TerminalFailure(
            "interrupted" if interrupted else "grouped_grid_seed_exception",
            f"{type(error).__name__}: {error}",
        ),
    )


def _select(
    attempts: Sequence[GroupedGridSeedOutcome],
) -> tuple[GroupedGridSelection, GroupedGridSeedOutcome]:
    eligible = tuple(
        attempt
        for attempt in attempts
        if attempt.disposition is GridSeedDisposition.ELIGIBLE
        and attempt.candidate is not None
    )
    selection = GroupedGridSelection.from_candidate_records(eligible)
    return selection, selection.selected_record


def _component_lineage_matches(
    accepted: AcceptedGroupedGridDirectTrfResult,
    selected: GroupedGridSeedOutcome,
    seed: GridSeed,
    problem: OptimizationProblem,
) -> bool:
    seed_problem = seed.problem
    seed_decomposition = selected.seed_decomposition
    component_invocation = selected.component_invocation
    if (
        seed_problem is None
        or seed_problem.source_snapshot != problem.source_snapshot
        or not isinstance(seed_problem.derivation, GridSeedProblemDerivation)
        or seed_problem.derivation.root_problem_identity != problem.identity
        or seed_problem.derivation.seed_identity != seed.coordinate_identity
        or seed_problem.derivation.seed_ordinal != seed.ordinal
        or seed_decomposition is None
        or seed_decomposition.root_problem_identity != seed_problem.identity
        or seed_decomposition.identity != accepted.decomposition.identity
        or component_invocation is None
        or component_invocation.root_problem_identity != seed_problem.identity
        or component_invocation.decomposition_identity
        != accepted.decomposition.identity
    ):
        return False
    root_components = accepted.decomposition.components
    seed_components = seed_decomposition.components
    direct_invocations = component_invocation.component_invocations
    if not (
        len(root_components)
        == len(seed_components)
        == len(selected.components)
        == len(direct_invocations)
    ):
        return False
    root_scales = dict(
        zip(
            problem.controlled_ids,
            accepted.workflow_invocation.x_scale,
            strict=True,
        )
    )
    for root_component, seed_component, outcome, direct in zip(
        root_components,
        seed_components,
        selected.components,
        direct_invocations,
        strict=True,
    ):
        component_candidate = outcome.candidate
        if component_candidate is None:
            return False
        materialization = component_candidate.materialization
        if (
            seed_component.identity != root_component.identity
            or seed_component.controlled_ids != root_component.controlled_ids
            or outcome.component_identity != root_component.identity
            or outcome.controlled_ids != root_component.controlled_ids
            or outcome.disposition is not ComponentDisposition.SUCCEEDED
            or outcome.execution_identity != component_candidate.execution_identity
            or direct.problem_identity != seed_component.problem.identity
            or component_candidate.problem_identity != direct.problem_identity
            or component_candidate.invocation_identity != direct.identity
            or direct.objective_request_budget
            != accepted.workflow_invocation.objective_request_budget
            or direct.x_scale
            != tuple(
                root_scales[param_id] for param_id in root_component.controlled_ids
            )
            or direct.ftol != accepted.workflow_invocation.ftol
            or direct.xtol != accepted.workflow_invocation.xtol
            or direct.gtol != accepted.workflow_invocation.gtol
            or materialization.problem_identity != component_candidate.problem_identity
            or materialization.invocation_identity
            != component_candidate.invocation_identity
            or materialization.execution_identity
            != component_candidate.execution_identity
            or materialization.evaluation_identity
            != component_candidate.evaluation_result.identity
            or component_candidate.evaluation_result.plan_identity
            != seed_component.problem.evaluation_plan_identity
        ):
            return False
    return True


def _root_candidate_matches(
    candidate: MaterializedDirectTrfCandidate,
    problem: OptimizationProblem,
    invocation_identity: str,
) -> bool:
    materialization = candidate.materialization
    return (
        candidate.problem_identity == problem.identity
        and candidate.invocation_identity == invocation_identity
        and materialization.problem_identity == problem.identity
        and materialization.invocation_identity == invocation_identity
        and materialization.execution_identity == candidate.execution_identity
        and materialization.evaluation_identity == candidate.evaluation_result.identity
        and candidate.evaluation_result.plan_identity
        == problem.evaluation_plan_identity
        and candidate.evaluation_result.parameterization_identity
        == problem.evaluator_parameterization_identity
        and tuple(candidate.evaluation_result.resolved_values) == problem.commit_scope
    )


def _aggregate_lineage_matches(
    accepted: AcceptedGroupedGridDirectTrfResult,
    selected: GroupedGridSeedOutcome,
    candidate: GridCandidate,
    fresh: MaterializedDirectTrfCandidate,
    problem: OptimizationProblem,
) -> bool:
    aggregate = candidate.candidate
    invocation_identity = accepted.workflow_invocation.identity
    return (
        candidate.seed_identity == selected.seed_identity
        and candidate.seed_ordinal == selected.seed_ordinal
        and selected.execution_identity == aggregate.execution_identity
        and _root_candidate_matches(aggregate, problem, invocation_identity)
        and _root_candidate_matches(fresh, problem, invocation_identity)
        and fresh.vector == candidate.vector
        and fresh.chi_square == candidate.objective
    )


def _validate_acceptance_lineage(
    accepted: AcceptedGroupedGridDirectTrfResult,
    problem: OptimizationProblem,
) -> None:
    provenance = accepted.grouped_grid_provenance
    selected = accepted.selected_outcome
    candidate = selected.candidate
    fresh = accepted.fresh_candidate
    try:
        seed = accepted.workflow_invocation.seeds[provenance.selected_seed_ordinal]
    except (AttributeError, IndexError, TypeError):
        seed = None
    if (
        accepted.workflow_invocation.root_problem_identity != problem.identity
        or len(accepted.workflow_invocation.x_scale) != len(problem.controlled_ids)
        or accepted.decomposition.root_problem_identity != problem.identity
        or accepted.decomposition.identity != provenance.decomposition_identity
        or provenance.root_problem_identity != problem.identity
        or provenance.workflow_invocation_identity
        != accepted.workflow_invocation.identity
        or provenance.selection_identity != accepted.selection.identity
        or seed is None
        or seed.identity != accepted.selected_seed.identity
        or provenance.selected_seed_identity != seed.identity
        or provenance.selected_seed_ordinal != seed.ordinal
        or selected.seed_identity != seed.identity
        or selected.seed_ordinal != seed.ordinal
        or accepted.selection.selected_seed_identity != seed.identity
        or accepted.selection.selected_seed_ordinal != seed.ordinal
        or not _component_lineage_matches(accepted, selected, seed, problem)
        or selected.identity != provenance.selected_outcome_identity
        or accepted.selection.selected_record.identity != selected.identity
        or candidate is None
        or candidate.identity != provenance.aggregate_candidate_identity
        or candidate.candidate.materialization.identity
        != provenance.aggregate_materialization_identity
        or not _aggregate_lineage_matches(
            accepted,
            selected,
            candidate,
            fresh,
            problem,
        )
        or fresh.materialization.identity
        != provenance.accepted_materialization_identity
        or fresh.evaluation_result.identity != provenance.accepted_evaluation_identity
        or accepted.problem_identity != problem.identity
        or accepted.vector != fresh.vector
        or accepted.chi_square != fresh.chi_square
        or accepted.materialization_identity != fresh.materialization.identity
        or accepted.evaluation_result.identity != fresh.evaluation_result.identity
        or accepted.origin_context_identity != provenance.identity
    ):
        raise DirectTrfConstructionError(
            "Accepted grouped GRID result lacks canonical aggregate authority"
        )


def _execute_all_seeds(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    token: CancellationToken,
) -> tuple[GroupedGridSeedOutcome, ...] | GroupedGridDirectTrfOutcome:
    attempted: list[GroupedGridSeedOutcome] = []
    for index, seed in enumerate(invocation.seeds):
        if token.is_cancelled:
            attempted.extend(
                _not_started_seed(item, decomposition)
                for item in invocation.seeds[index:]
            )
            return GroupedGridDirectTrfOutcome(
                GroupedGridDirectTrfTerminal.CANCELLED,
                tuple(attempted),
            )
        try:
            attempt = _execute_seed(
                problem,
                decomposition,
                invocation,
                seed,
                parameterization,
                engine,
                token,
            )
        except KeyboardInterrupt as error:
            attempt = _seed_exception(
                seed,
                decomposition,
                error,
                interrupted=True,
            )
        except Exception as error:  # noqa: BLE001 - orchestration fails closed
            attempt = _seed_exception(seed, decomposition, error)
        attempted.append(attempt)
        if attempt.disposition in {
            GridSeedDisposition.CANCELLED,
            GridSeedDisposition.INTERRUPTED,
            GridSeedDisposition.IMPLEMENTATION_FAILURE,
            GridSeedDisposition.BACKEND_FAILURE,
        }:
            attempted.extend(
                _not_started_seed(item, decomposition)
                for item in invocation.seeds[index + 1 :]
            )
            terminal = (
                GroupedGridDirectTrfTerminal.CANCELLED
                if attempt.disposition is GridSeedDisposition.CANCELLED
                else GroupedGridDirectTrfTerminal.INTERRUPTED
                if attempt.disposition is GridSeedDisposition.INTERRUPTED
                else GroupedGridDirectTrfTerminal.EXECUTION_FAILURE
            )
            return GroupedGridDirectTrfOutcome(
                terminal,
                tuple(attempted),
                failure=attempt.failure,
            )
    return tuple(attempted)


def commit_grouped_grid_accepted_fit(
    accepted: AcceptedGroupedGridDirectTrfResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Validate grouped-GRID lineage, then use the generic commit seam."""
    _validate_acceptance_lineage(accepted, problem)
    return commit_accepted_fit(
        accepted,
        authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=analysis_values,
    )


def execute_grouped_grid_direct_trf(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> GroupedGridDirectTrfOutcome:
    """Run exact components per seed and accept only the best root aggregate."""
    _validate_context(problem, decomposition, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    seed_execution = _execute_all_seeds(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        token,
    )
    if isinstance(seed_execution, GroupedGridDirectTrfOutcome):
        return seed_execution
    attempts = seed_execution
    if not any(
        attempt.disposition is GridSeedDisposition.ELIGIBLE for attempt in attempts
    ):
        return GroupedGridDirectTrfOutcome(
            GroupedGridDirectTrfTerminal.NO_ELIGIBLE_CANDIDATE,
            attempts,
        )
    selection, selected = _select(attempts)
    aggregate_candidate = cast("GridCandidate", selected.candidate)
    execution_identity = _identity(
        "native-grouped-grid-selection-execution",
        (
            invocation.identity,
            decomposition.identity,
            selection.identity,
            selected.identity,
            aggregate_candidate.identity,
        ),
    )
    if token.is_cancelled:
        return GroupedGridDirectTrfOutcome(
            GroupedGridDirectTrfTerminal.CANCELLED,
            attempts,
            selection,
        )
    try:
        fresh = _fresh_root_candidate(
            problem,
            invocation.identity,
            execution_identity,
            aggregate_candidate.vector,
            parameterization,
            engine,
        )
    except KeyboardInterrupt:
        return GroupedGridDirectTrfOutcome(
            GroupedGridDirectTrfTerminal.INTERRUPTED,
            attempts,
            selection,
        )
    except Exception as error:  # noqa: BLE001 - root promotion fails closed
        return GroupedGridDirectTrfOutcome(
            GroupedGridDirectTrfTerminal.MATERIALIZATION_FAILURE,
            attempts,
            selection,
            failure=TerminalFailure(
                "grouped_grid_selection_materialization_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    if token.is_cancelled:
        return GroupedGridDirectTrfOutcome(
            GroupedGridDirectTrfTerminal.CANCELLED,
            attempts,
            selection,
        )
    provenance = GroupedGridSelectionProvenance(
        invocation.identity,
        problem.identity,
        decomposition.identity,
        selection.identity,
        selected.seed_identity,
        selected.seed_ordinal,
        selected.identity,
        aggregate_candidate.identity,
        aggregate_candidate.candidate.materialization.identity,
        fresh.materialization.identity,
        fresh.evaluation_result.identity,
    )
    base_accepted = _accept_materialized_fit_for_derived_workflow(
        problem=problem,
        invocation_identity=invocation.identity,
        execution_identity=execution_identity,
        materialization=fresh.materialization,
        vector=fresh.vector,
        chi_square=fresh.chi_square,
        evaluation_result=fresh.evaluation_result,
        authority_context_identity=provenance.identity,
    )
    selected_seed = invocation.seeds[selected.seed_ordinal]
    accepted = AcceptedGroupedGridDirectTrfResult.from_accepted(
        base_accepted,
        provenance,
        invocation,
        decomposition,
        selection,
        selected_seed,
        selected,
        fresh,
    )
    _validate_acceptance_lineage(accepted, problem)
    authority = _grant_derived_fit_commit_authority(accepted, problem)
    return GroupedGridDirectTrfOutcome(
        GroupedGridDirectTrfTerminal.ACCEPTED,
        attempts,
        selection,
        accepted,
        authority,
    )
