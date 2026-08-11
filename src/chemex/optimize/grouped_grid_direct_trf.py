"""Exact fit-component GRID to native Direct TRF qualification (#596).

This isolated qualification path composes the grouped Direct TRF contract with
the Cartesian GRID contract.  It is intentionally not wired into production
dispatch.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from chemex.evaluation.native import EvaluationEngine
from chemex.optimize import direct_trf as direct_trf_owner
from chemex.optimize.direct_trf import (
    CancellationToken,
    CandidateMaterialization,
    CommitReceipt,
    ComponentProblemDerivation,
    DirectTrfConstructionError,
    DirectTrfInvocation,
    GridSeedProblemDerivation,
    LiveFitCommitAuthority,
    OptimizationProblem,
    TerminalFailure,
)
from chemex.optimize.grid_direct_trf import (
    AcceptedGridDirectTrfResult,
    GridCandidate,
    GridCoordinate,
    GridDirectTrfInvocation,
    GridDirectTrfOutcome,
    GridDirectTrfTerminal,
    GridSeed,
    GridSeedDisposition,
    GridSelection,
    _finalize_grid_candidates,
    _materialized_grid_candidate_matches,
    _validate_grid_context,
    commit_grid_accepted_fit,
)
from chemex.optimize.grouped_direct_trf import (
    ComponentDisposition,
    FitComponentOutcome,
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfOutcome,
    GroupedDirectTrfTerminal,
    _aggregate_vector,
    _grouped_materialization_failure,
    _root_projection_matches,
    _validate_component_projections,
    execute_direct_trf_components,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.values import AnalysisValues

_SCHEMA_VERSION = 1


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
    root_decomposition_identity: str = field(kw_only=True)
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
            not self.root_decomposition_identity
            or self.seed_decomposition is None
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
                    self.root_decomposition_identity,
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

    def validate_for_grid_selection(
        self,
        problem: OptimizationProblem,
        invocation: GridDirectTrfInvocation,
        seed: GridSeed,
    ) -> None:
        """Prove exact grouped component ownership before GRID selection."""
        _validate_grouped_grid_seed(self, problem, invocation, seed)


@dataclass(frozen=True, slots=True)
class GroupedGridDirectTrfOutcome:
    """Thin grouped-evidence view over canonical #595 GRID lifecycle state."""

    attempts: tuple[GroupedGridSeedOutcome, ...]
    grid_outcome: GridDirectTrfOutcome = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        if tuple(record.identity for record in self.grid_outcome.attempts) != tuple(
            record.identity for record in self.attempts
        ):
            raise ValueError("Grouped evidence differs from canonical GRID outcome")

    @property
    def terminal(self) -> GridDirectTrfTerminal:
        return self.grid_outcome.terminal

    @property
    def selection(self) -> GridSelection | None:
        return self.grid_outcome.selection

    @property
    def accepted_result(self) -> AcceptedGridDirectTrfResult | None:
        return self.grid_outcome.accepted_result

    @property
    def commit_authority(self) -> LiveFitCommitAuthority | None:
        return self.grid_outcome.commit_authority

    @property
    def failure(self) -> TerminalFailure | None:
        return self.grid_outcome.failure


GroupedGridDirectTrfTerminal = GridDirectTrfTerminal


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


def _failed_component_seed_outcome(
    seed: GridSeed,
    root_decomposition_identity: str,
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
        root_decomposition_identity=root_decomposition_identity,
        seed_decomposition=seed_decomposition,
    )


def _failed_aggregate_seed_outcome(
    seed: GridSeed,
    root_decomposition_identity: str,
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
        root_decomposition_identity=root_decomposition_identity,
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
            root_decomposition_identity=decomposition.identity,
        )
    if seed.problem is None:
        raise DirectTrfConstructionError("Feasible grouped GRID seed lacks a problem")
    problem.validate_derived_problem(seed.problem)
    seed_decomposition = FitDecomposition.from_root(
        seed.problem,
        parameterization,
        engine,
        cancellation=token,
    )
    if seed_decomposition.identity != decomposition.identity or tuple(
        component.identity for component in seed_decomposition.components
    ) != tuple(component.identity for component in decomposition.components):
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
            decomposition.identity,
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
            decomposition.identity,
            seed_decomposition,
            components,
            component_invocation,
            projections,
        )
    vector_or_failure = _aggregate_vector(
        problem,
        seed_decomposition,
        components,
    )
    if isinstance(vector_or_failure, TerminalFailure):
        return _failed_aggregate_seed_outcome(
            seed,
            decomposition.identity,
            seed_decomposition,
            components,
            component_invocation,
            GroupedDirectTrfOutcome(
                GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
                components,
                failure=vector_or_failure,
            ),
        )
    vector = vector_or_failure
    materialized = direct_trf_owner.materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=vector,
        invocation_identity=invocation.identity,
        execution_identity=lambda aggregate: _identity(
            "native-grouped-grid-seed-execution",
            (
                invocation.identity,
                seed.identity,
                component_invocation.identity,
                tuple(component.identity for component in components),
                aggregate.identity,
            ),
        ),
        cancellation=token,
    )
    if isinstance(materialized, direct_trf_owner.RootMaterializationFailure):
        grouped_failure = _grouped_materialization_failure(materialized, components)
        return _failed_aggregate_seed_outcome(
            seed,
            decomposition.identity,
            seed_decomposition,
            components,
            component_invocation,
            grouped_failure,
        )
    if isinstance(materialized, CandidateMaterialization):
        return _failed_aggregate_seed_outcome(
            seed,
            decomposition.identity,
            seed_decomposition,
            components,
            component_invocation,
            GroupedDirectTrfOutcome(
                GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
                components,
                failure=materialized.failure,
            ),
        )
    aggregate = materialized.evaluation_result
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
            root_decomposition_identity=decomposition.identity,
            seed_decomposition=seed_decomposition,
        )
    candidate = GridCandidate(seed.identity, seed.ordinal, materialized)
    record = GroupedGridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.ELIGIBLE,
        components,
        component_invocation,
        materialized.execution_identity,
        candidate.objective,
        candidate,
        root_decomposition_identity=decomposition.identity,
        seed_decomposition=seed_decomposition,
    )
    record.validate_for_grid_selection(problem, invocation, seed)
    return record


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
        root_decomposition_identity=decomposition.identity,
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
        root_decomposition_identity=decomposition.identity,
    )


def _validate_grouped_grid_seed(
    selected: GroupedGridSeedOutcome,
    problem: OptimizationProblem,
    invocation: GridDirectTrfInvocation,
    seed: GridSeed,
) -> None:
    """Reject foreign or cross-seed component evidence before selection."""
    seed_problem = seed.problem
    seed_decomposition = selected.seed_decomposition
    component_invocation = selected.component_invocation
    aggregate = selected.candidate
    materialized_aggregate = None if aggregate is None else aggregate.candidate
    aggregate_evaluation = (
        None
        if materialized_aggregate is None
        else materialized_aggregate.evaluation_result
    )
    seed_derivation = None if seed_problem is None else seed_problem.derivation
    if seed_problem is None:
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        )
    try:
        problem.validate_derived_problem(seed_problem)
    except DirectTrfConstructionError as error:
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        ) from error
    if (
        invocation.root_problem_identity != problem.identity
        or selected.seed_identity != seed.identity
        or selected.seed_ordinal != seed.ordinal
        or selected.axis_items != seed.axis_items
        or selected.start != seed.start
        or selected.disposition is not GridSeedDisposition.ELIGIBLE
        or not isinstance(seed_derivation, GridSeedProblemDerivation)
        or seed_derivation.root_problem_identity != problem.identity
        or seed_derivation.seed_identity != seed.coordinate_identity
        or seed_derivation.seed_ordinal != seed.ordinal
        or seed_derivation.axis_items != seed.axis_items
        or seed_derivation.start != seed.start
        or seed_problem.controlled_ids != problem.controlled_ids
        or seed_decomposition is None
        or seed_decomposition.root_problem_identity != seed_problem.identity
        or seed_decomposition.root_plan_identity
        != seed_problem.evaluation_plan_identity
        or seed_decomposition.identity != selected.root_decomposition_identity
        or component_invocation is None
        or component_invocation.root_problem_identity != seed_problem.identity
        or component_invocation.decomposition_identity != seed_decomposition.identity
        or aggregate is None
        or materialized_aggregate is None
        or aggregate.seed_identity != seed.identity
        or aggregate.seed_ordinal != seed.ordinal
        or selected.execution_identity != materialized_aggregate.execution_identity
        or selected.objective != aggregate.objective
        or aggregate_evaluation is None
        or not _materialized_grid_candidate_matches(
            materialized_aggregate,
            problem,
            invocation.identity,
            materialized_aggregate.execution_identity,
        )
    ):
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        )
    seed_components = seed_decomposition.components
    direct_invocations = component_invocation.component_invocations
    if not (
        len(seed_components) == len(selected.components) == len(direct_invocations)
        and seed_components
    ):
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        )
    root_scales = dict(
        zip(
            problem.controlled_ids,
            invocation.x_scale,
            strict=True,
        )
    )
    expected_execution_identity = _identity(
        "native-grouped-grid-seed-execution",
        (
            invocation.identity,
            seed.identity,
            component_invocation.identity,
            tuple(component.identity for component in selected.components),
            aggregate_evaluation.identity,
        ),
    )
    if selected.execution_identity != expected_execution_identity:
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        )
    component_assignments: dict[str, float] = {}
    for seed_component, outcome, direct in zip(
        seed_components,
        selected.components,
        direct_invocations,
        strict=True,
    ):
        component_candidate = outcome.candidate
        component_problem = seed_component.problem
        component_derivation = component_problem.derivation
        seed_coordinates = dict(
            zip(seed_problem.controlled_ids, seed_problem.start, strict=True)
        )
        try:
            seed_problem.validate_derived_problem(component_problem)
        except DirectTrfConstructionError as error:
            raise DirectTrfConstructionError(
                "Grouped GRID seed lacks exact component ownership"
            ) from error
        if (
            component_candidate is None
            or outcome.component_identity != seed_component.identity
            or outcome.controlled_ids != seed_component.controlled_ids
            or outcome.disposition is not ComponentDisposition.SUCCEEDED
            or outcome.execution_identity != component_candidate.execution_identity
            or not isinstance(component_derivation, ComponentProblemDerivation)
            or component_derivation.root_problem_identity != seed_problem.identity
            or component_derivation.component_identity != seed_component.identity
            or component_derivation.controlled_ids != seed_component.controlled_ids
            or component_problem.start
            != tuple(
                seed_coordinates[param_id] for param_id in seed_component.controlled_ids
            )
            or direct.problem_identity != component_problem.identity
            or component_candidate.problem_identity != direct.problem_identity
            or component_candidate.invocation_identity != direct.identity
            or direct.objective_request_budget != invocation.objective_request_budget
            or direct.x_scale
            != tuple(
                root_scales[param_id] for param_id in seed_component.controlled_ids
            )
            or direct.ftol != invocation.ftol
            or direct.xtol != invocation.xtol
            or direct.gtol != invocation.gtol
            or not _materialized_grid_candidate_matches(
                component_candidate,
                component_problem,
                direct.identity,
                component_candidate.execution_identity,
            )
        ):
            raise DirectTrfConstructionError(
                "Grouped GRID seed lacks exact component ownership"
            )
        component_assignments.update(
            zip(
                seed_component.controlled_ids,
                component_candidate.vector,
                strict=True,
            )
        )
    if materialized_aggregate.vector != tuple(
        component_assignments[param_id] for param_id in problem.controlled_ids
    ):
        raise DirectTrfConstructionError(
            "Grouped GRID seed lacks exact component ownership"
        )


def _lineage_failure(
    record: GroupedGridSeedOutcome,
    error: BaseException,
) -> GroupedGridSeedOutcome:
    """Strip aggregate eligibility from one invalid grouped seed record."""
    return GroupedGridSeedOutcome(
        record.seed_identity,
        record.seed_ordinal,
        record.axis_items,
        record.start,
        GridSeedDisposition.IMPLEMENTATION_FAILURE,
        record.components,
        record.component_invocation,
        failure=TerminalFailure(
            "grouped_grid_seed_lineage_failure",
            f"{type(error).__name__}: {error}",
        ),
        root_decomposition_identity=record.root_decomposition_identity,
        seed_decomposition=record.seed_decomposition,
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
        try:
            cancelled = token.is_cancelled
        except KeyboardInterrupt:
            attempted.extend(
                _not_started_seed(item, decomposition)
                for item in invocation.seeds[index:]
            )
            records = tuple(attempted)
            return GroupedGridDirectTrfOutcome(
                records,
                GridDirectTrfOutcome(
                    GridDirectTrfTerminal.INTERRUPTED,
                    records,
                ),
            )
        if cancelled:
            attempted.extend(
                _not_started_seed(item, decomposition)
                for item in invocation.seeds[index:]
            )
            records = tuple(attempted)
            return GroupedGridDirectTrfOutcome(
                records,
                GridDirectTrfOutcome(
                    GridDirectTrfTerminal.CANCELLED,
                    records,
                ),
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
                GridDirectTrfTerminal.CANCELLED
                if attempt.disposition is GridSeedDisposition.CANCELLED
                else GridDirectTrfTerminal.INTERRUPTED
                if attempt.disposition is GridSeedDisposition.INTERRUPTED
                else GridDirectTrfTerminal.EXECUTION_FAILURE
            )
            records = tuple(attempted)
            return GroupedGridDirectTrfOutcome(
                records,
                GridDirectTrfOutcome(
                    terminal,
                    records,
                    failure=attempt.failure,
                ),
            )
    return tuple(attempted)


def _validate_grouped_attempts(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GridDirectTrfInvocation,
    attempts: tuple[GroupedGridSeedOutcome, ...],
) -> GroupedGridDirectTrfOutcome | None:
    """Fail closed before canonical GRID sees malformed grouped evidence."""
    records = list(attempts)
    failure: BaseException | None = None
    failed_index: int | None = None
    if len(records) != len(invocation.seeds):
        failure = DirectTrfConstructionError(
            "Grouped GRID occurrence lacks its complete canonical seed table"
        )
    else:
        for index, (record, seed) in enumerate(
            zip(records, invocation.seeds, strict=True)
        ):
            if (
                record.seed_identity != seed.identity
                or record.seed_ordinal != seed.ordinal
                or record.root_decomposition_identity != decomposition.identity
            ):
                failure = DirectTrfConstructionError(
                    "Grouped GRID occurrence contains foreign seed evidence"
                )
                failed_index = index
                break
            if record.disposition is GridSeedDisposition.ELIGIBLE:
                try:
                    record.validate_for_grid_selection(problem, invocation, seed)
                except Exception as error:  # noqa: BLE001 - evidence fails closed
                    failure = error
                    failed_index = index
                    break
    if failure is None:
        return None
    if failed_index is not None:
        records[failed_index] = _lineage_failure(records[failed_index], failure)
    frozen = tuple(records)
    return GroupedGridDirectTrfOutcome(
        frozen,
        GridDirectTrfOutcome(
            GridDirectTrfTerminal.EXECUTION_FAILURE,
            frozen,
            failure=TerminalFailure(
                "grouped_grid_occurrence_lineage_failure",
                f"{type(failure).__name__}: {failure}",
            ),
        ),
    )


def commit_grouped_grid_accepted_fit(
    accepted: AcceptedGridDirectTrfResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Validate grouped-GRID lineage, then use the generic commit seam."""
    selected = accepted.selected_outcome
    if not isinstance(selected, GroupedGridSeedOutcome):
        raise DirectTrfConstructionError(
            "Accepted grouped GRID result lacks grouped seed evidence"
        )
    selected.validate_for_grid_selection(
        problem,
        accepted.workflow_invocation,
        accepted.selected_seed,
    )
    return commit_grid_accepted_fit(
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
    invalid = _validate_grouped_attempts(
        problem,
        decomposition,
        invocation,
        attempts,
    )
    if invalid is not None:
        return invalid
    grid_outcome = _finalize_grid_candidates(
        problem,
        invocation,
        attempts,
        parameterization,
        engine,
        token,
    )
    return GroupedGridDirectTrfOutcome(attempts, grid_outcome)
