"""Exact profiled chi-square GRID execution and surface reconstruction."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import combinations, product

import numpy as np

from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CandidateMaterialization,
    DirectTrfCandidateTerminal,
    DirectTrfInterrupted,
    DirectTrfInvocation,
    DirectTrfScalePolicy,
    LiveFitCommitAuthority,
    MaterializedDirectTrfCandidate,
    OptimizationProblem,
    RootMaterializationFailure,
    TerminalFailure,
    _accept_materialized_fit_for_derived_workflow,
    _grant_derived_fit_commit_authority,
    canonical_chi_square,
    execute_direct_trf_candidate,
    materialize_root_candidate,
)
from chemex.optimize.grouped_direct_trf import _profile_dependencies
from chemex.optimize.progress import ContextualProgressObserver, FitProgressContext
from chemex.parameters.parameterization import ActiveParameterization
from chemex.typing import Array


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": 1},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


class ProfiledGridConstructionError(ValueError):
    """The requested GRID cannot be represented as an exact profiled product."""


class ProfiledGridPointStatus(StrEnum):
    """Truthful terminal status for one factor-local profiled point."""

    SUCCESS = "success"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    NOT_STARTED = "not_started"


@dataclass(frozen=True, slots=True)
class ProfiledGridFactor:
    """One exact objective factor after treating GRID coordinates as held."""

    ordinal: int
    profile_indices: tuple[int, ...]
    grid_ids: tuple[str, ...]
    nuisance_ids: tuple[str, ...]
    profile_keys: tuple[tuple[int, int], ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.ordinal < 0:
            raise ProfiledGridConstructionError(
                "GRID factor ordinal cannot be negative"
            )
        if not self.profile_indices:
            raise ProfiledGridConstructionError(
                "GRID factor requires retained profiles"
            )
        if len(set(self.grid_ids)) != len(self.grid_ids) or len(
            set(self.nuisance_ids)
        ) != len(self.nuisance_ids):
            raise ProfiledGridConstructionError(
                "GRID factor coordinates must be unique"
            )
        if set(self.grid_ids) & set(self.nuisance_ids):
            raise ProfiledGridConstructionError(
                "GRID and nuisance factor coordinates must be disjoint"
            )
        if self.profile_keys and len(self.profile_keys) != len(self.profile_indices):
            raise ProfiledGridConstructionError(
                "GRID factor profile keys differ from its profile indices"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "profiled-grid-factor",
                (
                    self.ordinal,
                    self.profile_indices,
                    self.grid_ids,
                    self.nuisance_ids,
                    self.profile_keys,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ProfiledGridPoint:
    """One raw factor-grid point, including unsuccessful evidence."""

    ordinal: int
    axis_items: tuple[tuple[str, float], ...]
    status: ProfiledGridPointStatus
    chi_square: float | None = None
    nuisance_items: tuple[tuple[str, float], ...] = ()
    objective_evaluations: int = 0
    failure: str | None = None

    def __post_init__(self) -> None:
        successful = self.status is ProfiledGridPointStatus.SUCCESS
        if successful != (self.chi_square is not None):
            raise ProfiledGridConstructionError(
                "Only a successful GRID point may expose chi-square"
            )
        chi_square = self.chi_square
        if successful and (
            chi_square is None
            or not math.isfinite(chi_square)
            or self.objective_evaluations < 1
        ):
            raise ProfiledGridConstructionError(
                "Successful GRID point requires finite objective evidence"
            )
        if not successful and self.nuisance_items:
            raise ProfiledGridConstructionError(
                "Unsuccessful GRID point cannot expose nuisance values"
            )


@dataclass(frozen=True, slots=True)
class ProfiledGridFactorResult:
    """Retained raw point evidence for one exact objective factor."""

    factor: ProfiledGridFactor
    points: tuple[ProfiledGridPoint, ...]

    def __post_init__(self) -> None:
        if not self.points:
            raise ProfiledGridConstructionError("GRID factor has no evaluated points")
        expected_axes = self.factor.grid_ids
        if any(
            tuple(param_id for param_id, _value in point.axis_items) != expected_axes
            or tuple(param_id for param_id, _value in point.nuisance_items)
            != self.factor.nuisance_ids
            for point in self.points
            if point.status is ProfiledGridPointStatus.SUCCESS
        ):
            raise ProfiledGridConstructionError(
                "GRID factor point coordinates differ from their factor"
            )


@dataclass(frozen=True, slots=True)
class ProfiledGridSurface:
    """One exact numerical marginal surface."""

    axis_ids: tuple[str, ...]
    axis_values: tuple[tuple[float, ...], ...]
    chi_square: Array = field(repr=False, compare=False)


@dataclass(frozen=True, slots=True)
class ProfiledGridSelection:
    """One coherent exact joint GRID and conditional nuisance solution."""

    grid_items: tuple[tuple[str, float], ...]
    nuisance_items: tuple[tuple[str, float], ...]
    chi_square: float
    factor_point_ordinals: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ProfiledGridAggregate:
    """Exact aggregate selection plus reusable one- and two-dimensional surfaces."""

    selection: ProfiledGridSelection
    profiles_1d: Mapping[str, ProfiledGridSurface]
    profiles_2d: Mapping[tuple[str, str], ProfiledGridSurface]


class ProfiledGridTerminal(StrEnum):
    ACCEPTED = "accepted"
    NO_USABLE_FACTOR = "no_usable_factor"
    MATERIALIZATION_FAILURE = "materialization_failure"
    EXECUTION_FAILURE = "execution_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


@dataclass(frozen=True, slots=True)
class ProfiledGridOutcome:
    terminal: ProfiledGridTerminal
    factors: tuple[ProfiledGridFactorResult, ...]
    aggregate: ProfiledGridAggregate | None = None
    accepted_result: AcceptedFitResult | None = field(
        default=None, repr=False, compare=False
    )
    commit_authority: LiveFitCommitAuthority | None = field(
        default=None, repr=False, compare=False
    )
    failure: TerminalFailure | None = None

    def __post_init__(self) -> None:
        accepted = self.terminal is ProfiledGridTerminal.ACCEPTED
        if accepted != (
            self.aggregate is not None
            and self.accepted_result is not None
            and self.commit_authority is not None
        ):
            raise ProfiledGridConstructionError(
                "Only accepted profiled GRID execution has aggregate authority"
            )


def discover_profiled_grid_factors(  # noqa: C901 - exact factor proof is kept local
    profile_dependencies: Sequence[frozenset[str]],
    retained_profiles: Sequence[bool],
    profile_keys: Sequence[tuple[int, int]] | None = None,
    *,
    grid_ids: Sequence[str],
    controlled_ids: Sequence[str],
) -> tuple[ProfiledGridFactor, ...]:
    """Find exact factors after removing globally shared held GRID coupling."""
    dependencies = tuple(profile_dependencies)
    retained = tuple(retained_profiles)
    keys = (
        tuple((index, index) for index in range(len(dependencies)))
        if profile_keys is None
        else tuple(profile_keys)
    )
    ordered_grid = tuple(grid_ids)
    ordered_controls = tuple(controlled_ids)
    if len(dependencies) != len(retained) or len(dependencies) != len(keys):
        raise ProfiledGridConstructionError(
            "Profile dependency and retained-objective records differ"
        )
    if not set(ordered_grid).issubset(ordered_controls):
        raise ProfiledGridConstructionError(
            "GRID axes must be final controlled coordinates"
        )
    objective_indices = tuple(index for index, keep in enumerate(retained) if keep)
    variable_indices = tuple(
        index for index in objective_indices if dependencies[index]
    )
    if not variable_indices:
        raise ProfiledGridConstructionError(
            "GRID requires a retained objective dependency"
        )
    observed = set().union(*(dependencies[index] for index in variable_indices))
    missing = set(ordered_controls) - observed
    if missing:
        raise ProfiledGridConstructionError(
            "Controlled coordinates lack a retained objective dependency: "
            + ", ".join(sorted(missing))
        )
    global_grid = set(ordered_grid).intersection(
        *(dependencies[index] for index in variable_indices)
    )

    parent = {index: index for index in variable_indices}

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left_position, left in enumerate(variable_indices):
        for right in variable_indices[left_position + 1 :]:
            if (dependencies[left] & dependencies[right]) - global_grid:
                left_root = find(left)
                right_root = find(right)
                if left_root != right_root:
                    parent[right_root] = left_root

    groups: dict[int, list[int]] = {}
    for index in variable_indices:
        groups.setdefault(find(index), []).append(index)
    ordered_groups = sorted(groups.values(), key=lambda indices: tuple(indices))
    constant_indices = tuple(
        index for index in objective_indices if not dependencies[index]
    )
    if constant_indices:
        ordered_groups.append(list(constant_indices))

    factors: list[ProfiledGridFactor] = []
    grid_set = set(ordered_grid)
    for ordinal, indices in enumerate(ordered_groups):
        factor_dependencies = set().union(*(dependencies[index] for index in indices))
        factors.append(
            ProfiledGridFactor(
                ordinal,
                tuple(indices),
                tuple(
                    param_id
                    for param_id in ordered_grid
                    if param_id in factor_dependencies
                ),
                tuple(
                    param_id
                    for param_id in ordered_controls
                    if param_id not in grid_set and param_id in factor_dependencies
                ),
                tuple(keys[index] for index in indices),
            )
        )
    return tuple(factors)


def _factor_point_failure(
    ordinal: int,
    axis_items: tuple[tuple[str, float], ...],
    *,
    status: ProfiledGridPointStatus = ProfiledGridPointStatus.FAILED,
    evaluations: int = 0,
    failure: object,
) -> ProfiledGridPoint:
    return ProfiledGridPoint(
        ordinal,
        axis_items,
        status,
        objective_evaluations=evaluations,
        failure=str(failure),
    )


def _evaluate_only_factor_point(
    problem: OptimizationProblem,
    factor: ProfiledGridFactor,
    child_engine: EvaluationEngine,
    parameterization: ActiveParameterization,
    point_ordinal: int,
    axis_items: tuple[tuple[str, float], ...],
) -> ProfiledGridPoint:
    coordinate_values = dict(zip(problem.controlled_ids, problem.start, strict=True))
    coordinate_values.update(axis_items)
    vector = tuple(coordinate_values[param_id] for param_id in problem.controlled_ids)
    lifecycle = problem.lifecycle_frame(vector, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    evaluated = child_engine.new_evaluator().evaluate(frame)
    if isinstance(evaluated, EvaluationFailure):
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            failure=f"{evaluated.category}: {evaluated.message}",
            evaluations=1,
        )
    try:
        chi_square = canonical_chi_square(evaluated.residuals)
    except (TypeError, ValueError) as error:
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            failure=error,
            evaluations=1,
        )
    return ProfiledGridPoint(
        point_ordinal,
        axis_items,
        ProfiledGridPointStatus.SUCCESS,
        chi_square,
        (),
        1,
    )


def _fit_factor_point(
    problem: OptimizationProblem,
    factor: ProfiledGridFactor,
    child_engine: EvaluationEngine,
    parameterization: ActiveParameterization,
    point_ordinal: int,
    axis_items: tuple[tuple[str, float], ...],
    *,
    objective_request_budget: int,
    cancellation: CancellationToken,
    progress_observer: ContextualProgressObserver | None,
    factor_count: int,
    point_count: int,
) -> ProfiledGridPoint:
    child = problem.derive_profiled_grid_point(
        factor_identity=factor.identity,
        point_ordinal=point_ordinal,
        projected_plan_identity=child_engine.plan.identity,
        grid_items=axis_items,
        controlled_ids=factor.nuisance_ids,
    )
    invocation = DirectTrfInvocation.for_problem(
        child,
        objective_request_budget=objective_request_budget,
    )
    context = FitProgressContext(
        factor.ordinal + 1,
        factor_count,
        factor.nuisance_ids,
        point_ordinal + 1,
        point_count,
    )
    local_observer = (
        None
        if progress_observer is None
        else lambda event: progress_observer(context, event)
    )
    try:
        outcome = execute_direct_trf_candidate(
            child,
            invocation,
            parameterization,
            child_engine,
            cancellation=cancellation,
            progress_observer=local_observer,
        )
    except DirectTrfInterrupted as error:
        evaluations = error.execution.counters.objective_evaluations_completed
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            status=ProfiledGridPointStatus.INTERRUPTED,
            evaluations=evaluations,
            failure="Direct TRF interrupted",
        )
    evaluations = outcome.execution.counters.objective_evaluations_completed
    if outcome.terminal is DirectTrfCandidateTerminal.CANCELLED:
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            status=ProfiledGridPointStatus.CANCELLED,
            evaluations=evaluations,
            failure=outcome.execution.failure or "cancelled",
        )
    if outcome.terminal is not DirectTrfCandidateTerminal.SUCCESS:
        failure = (
            outcome.materialization.failure
            if outcome.materialization is not None
            else outcome.execution.failure
        )
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            evaluations=evaluations,
            failure=failure or outcome.terminal.value,
        )
    candidate = outcome.candidate
    if candidate is None:
        return _factor_point_failure(
            point_ordinal,
            axis_items,
            evaluations=evaluations,
            failure="Successful nuisance TRF lacks a candidate",
        )
    return ProfiledGridPoint(
        point_ordinal,
        axis_items,
        ProfiledGridPointStatus.SUCCESS,
        candidate.chi_square,
        tuple(zip(factor.nuisance_ids, candidate.vector, strict=True)),
        evaluations,
    )


def _validate_selected_factor_evidence(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    materialized: MaterializedDirectTrfCandidate,
    factor_results: Sequence[ProfiledGridFactorResult],
    aggregate: ProfiledGridAggregate,
) -> TerminalFailure | None:
    """Check selected factor evidence against one fresh complete root frame."""
    lifecycle = problem.lifecycle_frame(materialized.vector, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    factor_total = 0.0
    for result, selected_ordinal in zip(
        factor_results,
        aggregate.selection.factor_point_ordinals,
        strict=True,
    ):
        selected = next(
            (point for point in result.points if point.ordinal == selected_ordinal),
            None,
        )
        if selected is None or selected.chi_square is None:
            return TerminalFailure(
                "grid_factor_validation_failure",
                f"Selected point is unavailable for factor {result.factor.ordinal + 1}",
            )
        evaluated = (
            engine.project_profiles(result.factor.profile_indices)
            .new_evaluator()
            .evaluate(frame)
        )
        if isinstance(evaluated, EvaluationFailure):
            return TerminalFailure(
                "grid_factor_validation_failure",
                f"Fresh factor {result.factor.ordinal + 1} evaluation failed: "
                f"{evaluated.category}: {evaluated.message}",
                evaluated,
            )
        try:
            fresh_chi_square = canonical_chi_square(evaluated.residuals)
        except (TypeError, ValueError) as error:
            return TerminalFailure(
                "grid_factor_validation_failure",
                f"Fresh factor {result.factor.ordinal + 1} chi-square failed: {error}",
            )
        if not math.isclose(
            fresh_chi_square,
            selected.chi_square,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        ):
            return TerminalFailure(
                "grid_factor_validation_failure",
                f"Fresh factor {result.factor.ordinal + 1} chi-square "
                f"{fresh_chi_square!r} differs from selected evidence "
                f"{selected.chi_square!r}",
            )
        factor_total += fresh_chi_square
    if not math.isclose(
        factor_total,
        aggregate.selection.chi_square,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ) or not math.isclose(
        materialized.chi_square,
        aggregate.selection.chi_square,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        return TerminalFailure(
            "grid_aggregate_validation_failure",
            "Fresh complete-root chi-square, selected factor sum, and aggregate "
            "profile minimum differ",
        )
    return None


def execute_profiled_grid(  # noqa: C901 - closed factor/point lifecycle dispatcher
    problem: OptimizationProblem,
    axes: Mapping[str, Sequence[float]],
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    objective_request_budget: int,
    cancellation: CancellationToken | None = None,
    progress_observer: ContextualProgressObserver | None = None,
) -> ProfiledGridOutcome:
    """Evaluate exact factor-local profiled grids and accept one root aggregate."""
    problem.validate_parameterization(parameterization)
    if engine.plan.identity != problem.evaluation_plan_identity:
        raise ProfiledGridConstructionError(
            "Profiled GRID evaluator belongs to another root problem"
        )
    ordered_axes = {
        param_id: tuple(float(value) for value in values)
        for param_id, values in axes.items()
    }
    grid_ids = tuple(ordered_axes)
    if not grid_ids or not set(grid_ids).issubset(problem.controlled_ids):
        raise ProfiledGridConstructionError(
            "Profiled GRID axes must be active final FIT coordinates"
        )
    dependencies = _profile_dependencies(problem, parameterization, engine)
    factors = discover_profiled_grid_factors(
        dependencies,
        tuple(
            bool(profile.retained_observation_indices)
            for profile in engine.plan.profiles
        ),
        tuple(
            (profile.experiment_ordinal, profile.profile_ordinal)
            for profile in engine.plan.profiles
        ),
        grid_ids=grid_ids,
        controlled_ids=problem.controlled_ids,
    )
    if problem.affine_half_spaces or problem.affine_equalities:
        retained_indices = tuple(
            index
            for index, profile in enumerate(engine.plan.profiles)
            if profile.retained_observation_indices
        )
        factors = (
            ProfiledGridFactor(
                0,
                retained_indices,
                grid_ids,
                tuple(
                    param_id
                    for param_id in problem.controlled_ids
                    if param_id not in grid_ids
                ),
                tuple(
                    (
                        engine.plan.profiles[index].experiment_ordinal,
                        engine.plan.profiles[index].profile_ordinal,
                    )
                    for index in retained_indices
                ),
            ),
        )
    token = CancellationToken() if cancellation is None else cancellation
    factor_results: list[ProfiledGridFactorResult] = []
    for factor in factors:
        child_engine = engine.project_profiles(factor.profile_indices)
        coordinate_values = tuple(
            ordered_axes[param_id] for param_id in factor.grid_ids
        )
        point_count = math.prod(len(values) for values in coordinate_values) or 1
        assignments = product(*coordinate_values)
        points: list[ProfiledGridPoint] = []
        for point_ordinal, values in enumerate(assignments):
            if token.is_cancelled:
                points.append(
                    ProfiledGridPoint(
                        point_ordinal,
                        tuple(zip(factor.grid_ids, values, strict=True)),
                        ProfiledGridPointStatus.NOT_STARTED,
                        failure="cancelled before execution",
                    )
                )
                factor_results.append(ProfiledGridFactorResult(factor, tuple(points)))
                return ProfiledGridOutcome(
                    ProfiledGridTerminal.CANCELLED,
                    tuple(factor_results),
                    failure=TerminalFailure("cancelled", "GRID execution cancelled"),
                )
            axis_items = tuple(zip(factor.grid_ids, values, strict=True))
            try:
                point = (
                    _fit_factor_point(
                        problem,
                        factor,
                        child_engine,
                        parameterization,
                        point_ordinal,
                        axis_items,
                        objective_request_budget=objective_request_budget,
                        cancellation=token,
                        progress_observer=progress_observer,
                        factor_count=len(factors),
                        point_count=point_count,
                    )
                    if factor.nuisance_ids
                    else _evaluate_only_factor_point(
                        problem,
                        factor,
                        child_engine,
                        parameterization,
                        point_ordinal,
                        axis_items,
                    )
                )
            except KeyboardInterrupt:
                point = _factor_point_failure(
                    point_ordinal,
                    axis_items,
                    status=ProfiledGridPointStatus.INTERRUPTED,
                    failure="KeyboardInterrupt",
                )
            except Exception as error:  # noqa: BLE001 - one point fails independently
                point = _factor_point_failure(
                    point_ordinal,
                    axis_items,
                    failure=f"{type(error).__name__}: {error}",
                )
            points.append(point)
            if point.status is ProfiledGridPointStatus.INTERRUPTED:
                factor_results.append(ProfiledGridFactorResult(factor, tuple(points)))
                return ProfiledGridOutcome(
                    ProfiledGridTerminal.INTERRUPTED,
                    tuple(factor_results),
                    failure=TerminalFailure(
                        "interrupted", "GRID execution interrupted"
                    ),
                )
        result = ProfiledGridFactorResult(factor, tuple(points))
        factor_results.append(result)
        if not any(
            point.status is ProfiledGridPointStatus.SUCCESS for point in result.points
        ):
            return ProfiledGridOutcome(
                ProfiledGridTerminal.NO_USABLE_FACTOR,
                tuple(factor_results),
                failure=TerminalFailure(
                    "no_usable_factor",
                    f"GRID factor {factor.ordinal + 1} has no usable point",
                ),
            )

    try:
        aggregate = aggregate_profiled_grids(ordered_axes, factor_results)
    except ProfiledGridConstructionError as error:
        return ProfiledGridOutcome(
            ProfiledGridTerminal.EXECUTION_FAILURE,
            tuple(factor_results),
            failure=TerminalFailure("grid_aggregation_failure", str(error)),
        )
    selected_values = dict(aggregate.selection.grid_items)
    selected_values.update(aggregate.selection.nuisance_items)
    vector = tuple(
        selected_values.get(param_id, start)
        for param_id, start in zip(problem.controlled_ids, problem.start, strict=True)
    )
    invocation_identity = _identity(
        "profiled-grid-invocation",
        (
            problem.identity,
            tuple((param_id, ordered_axes[param_id]) for param_id in grid_ids),
            tuple(factor.identity for factor in factors),
            objective_request_budget,
            (DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM.value),
        ),
    )
    execution_identity = _identity(
        "profiled-grid-execution",
        (
            invocation_identity,
            tuple(
                (
                    result.factor.identity,
                    tuple(
                        (
                            point.ordinal,
                            point.axis_items,
                            point.status.value,
                            point.chi_square,
                            point.nuisance_items,
                        )
                        for point in result.points
                    ),
                )
                for result in factor_results
            ),
            aggregate.selection.grid_items,
            aggregate.selection.nuisance_items,
        ),
    )
    materialized = materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=vector,
        invocation_identity=invocation_identity,
        execution_identity=execution_identity,
        cancellation=token,
    )
    if not isinstance(materialized, MaterializedDirectTrfCandidate):
        failure = materialized.failure
        terminal = (
            ProfiledGridTerminal.CANCELLED
            if isinstance(
                materialized, (CandidateMaterialization, RootMaterializationFailure)
            )
            and materialized.terminal.value == "cancelled"
            else ProfiledGridTerminal.INTERRUPTED
            if isinstance(
                materialized, (CandidateMaterialization, RootMaterializationFailure)
            )
            and materialized.terminal.value == "interrupted"
            else ProfiledGridTerminal.MATERIALIZATION_FAILURE
        )
        return ProfiledGridOutcome(
            terminal,
            tuple(factor_results),
            aggregate,
            failure=failure,
        )
    validation_failure = _validate_selected_factor_evidence(
        problem,
        parameterization,
        engine,
        materialized,
        factor_results,
        aggregate,
    )
    if validation_failure is not None:
        return ProfiledGridOutcome(
            ProfiledGridTerminal.EXECUTION_FAILURE,
            tuple(factor_results),
            aggregate,
            failure=validation_failure,
        )
    authority_context = _identity(
        "profiled-grid-aggregate-acceptance",
        (
            invocation_identity,
            execution_identity,
            aggregate.selection.grid_items,
            aggregate.selection.nuisance_items,
            materialized.materialization.identity,
            materialized.evaluation_result.identity,
        ),
    )
    accepted = _accept_materialized_fit_for_derived_workflow(
        problem=problem,
        invocation_identity=invocation_identity,
        execution_identity=execution_identity,
        materialization=materialized.materialization,
        vector=materialized.vector,
        chi_square=materialized.chi_square,
        evaluation_result=materialized.evaluation_result,
        authority_context_identity=authority_context,
    )
    authority = _grant_derived_fit_commit_authority(accepted, problem)
    return ProfiledGridOutcome(
        ProfiledGridTerminal.ACCEPTED,
        tuple(factor_results),
        aggregate,
        accepted,
        authority,
    )


def _best_factor_point(
    result: ProfiledGridFactorResult,
    fixed: Mapping[str, float],
) -> ProfiledGridPoint | None:
    best: ProfiledGridPoint | None = None
    for point in result.points:
        if point.status is not ProfiledGridPointStatus.SUCCESS:
            continue
        objective = point.chi_square
        if objective is None:
            continue
        coordinates = dict(point.axis_items)
        if any(
            axis_id in fixed and coordinates[axis_id] != fixed[axis_id]
            for axis_id in result.factor.grid_ids
        ):
            continue
        best_objective = None if best is None else best.chi_square
        if best_objective is None or objective < best_objective:
            best = point
    return best


def _shared_local_axes(
    axis_ids: tuple[str, ...],
    factors: Sequence[ProfiledGridFactorResult],
) -> tuple[str, ...]:
    objective_factors = tuple(result for result in factors if result.factor.grid_ids)
    if not objective_factors:
        raise ProfiledGridConstructionError("GRID factors expose no GRID coordinate")
    incidence = {
        axis_id: tuple(
            index
            for index, result in enumerate(objective_factors)
            if axis_id in result.factor.grid_ids
        )
        for axis_id in axis_ids
    }
    missing = tuple(axis_id for axis_id, owners in incidence.items() if not owners)
    if missing:
        raise ProfiledGridConstructionError(
            "GRID coordinates lack objective factors: " + ", ".join(missing)
        )
    shared = tuple(
        axis_id
        for axis_id, owners in incidence.items()
        if len(owners) == len(objective_factors)
    )
    invalid = tuple(
        axis_id
        for axis_id, owners in incidence.items()
        if axis_id not in shared and len(owners) != 1
    )
    if invalid:
        raise ProfiledGridConstructionError(
            "Exact GRID factorization requires shared-global or factor-local axes; "
            "overlapping axes: " + ", ".join(invalid)
        )
    return shared


def _profile_surface_pairs(
    axis_ids: tuple[str, ...],
    factors: Sequence[ProfiledGridFactor],
) -> tuple[tuple[str, str], ...]:
    """Select pairs that coexist in at least one exact GRID factor."""
    axis_order = {axis_id: index for index, axis_id in enumerate(axis_ids)}
    pairs: set[tuple[str, str]] = set()
    for factor in factors:
        for left, right in combinations(factor.grid_ids, 2):
            pair = (
                (left, right) if axis_order[left] < axis_order[right] else (right, left)
            )
            pairs.add(pair)
    return tuple(
        sorted(
            pairs,
            key=lambda pair: (axis_order[pair[0]], axis_order[pair[1]]),
        )
    )


def aggregate_profiled_grids(  # noqa: C901 - exact shared/local reduction
    axes: Mapping[str, Sequence[float]],
    factors: Sequence[ProfiledGridFactorResult],
) -> ProfiledGridAggregate:
    """Recombine exact factor grids without constructing local-axis root products."""
    ordered_axes = {
        param_id: tuple(float(value) for value in values)
        for param_id, values in axes.items()
    }
    axis_ids = tuple(ordered_axes)
    factor_results = tuple(factors)
    shared_ids = _shared_local_axes(axis_ids, factor_results)
    value_indices = {
        axis_id: {value: index for index, value in enumerate(values)}
        for axis_id, values in ordered_axes.items()
    }
    factor_tables: list[Array] = []
    for result in factor_results:
        shape = tuple(len(ordered_axes[axis_id]) for axis_id in result.factor.grid_ids)
        table = np.full(shape, np.inf, dtype=np.float64)
        seen: set[tuple[int, ...]] = set()
        for point in result.points:
            try:
                indices = tuple(
                    value_indices[axis_id][value] for axis_id, value in point.axis_items
                )
            except KeyError as error:
                raise ProfiledGridConstructionError(
                    "GRID factor point lies outside the declared aggregate axes"
                ) from error
            if indices in seen:
                raise ProfiledGridConstructionError(
                    "GRID factor contains a duplicate concrete point"
                )
            seen.add(indices)
            if point.status is ProfiledGridPointStatus.SUCCESS:
                chi_square = point.chi_square
                if chi_square is None:  # guarded by ProfiledGridPoint invariants
                    raise ProfiledGridConstructionError(
                        "Successful GRID factor point lacks chi-square"
                    )
                table[indices] = chi_square
        factor_tables.append(table)

    def minimize_with_fixed(
        fixed: Mapping[str, float],
    ) -> tuple[float, tuple[ProfiledGridPoint, ...], dict[str, float]] | None:
        remaining_shared = tuple(
            axis_id for axis_id in shared_ids if axis_id not in fixed
        )
        best: tuple[float, tuple[ProfiledGridPoint, ...], dict[str, float]] | None = (
            None
        )
        for shared_values in product(
            *(ordered_axes[axis_id] for axis_id in remaining_shared)
        ):
            assignment = dict(fixed)
            assignment.update(zip(remaining_shared, shared_values, strict=True))
            selected_points: list[ProfiledGridPoint] = []
            total = 0.0
            for result in factor_results:
                point = _best_factor_point(result, assignment)
                if point is None:
                    break
                selected_points.append(point)
                objective = point.chi_square
                if objective is None:
                    break
                total += objective
            else:
                candidate = (total, tuple(selected_points), assignment)
                if best is None or total < best[0]:
                    best = candidate
        return best

    selected = minimize_with_fixed({})
    if selected is None:
        raise ProfiledGridConstructionError(
            "No coherent usable point exists across all GRID factors"
        )
    selected_chi_square, selected_points, selected_assignment = selected
    selected_grid = dict(selected_assignment)
    nuisance: dict[str, float] = {}
    for point in selected_points:
        selected_grid.update(
            {
                param_id: value
                for param_id, value in point.axis_items
                if param_id not in shared_ids
            }
        )
        for param_id, value in point.nuisance_items:
            if param_id in nuisance:
                raise ProfiledGridConstructionError(
                    f"Nuisance coordinate {param_id} appears in multiple GRID factors"
                )
            nuisance[param_id] = value
    if set(selected_grid) != set(axis_ids):
        raise ProfiledGridConstructionError(
            "Coherent GRID selection does not cover every concrete axis"
        )

    def surface(selection: tuple[str, ...]) -> ProfiledGridSurface:
        values = tuple(ordered_axes[axis_id] for axis_id in selection)
        chisqr = np.full(tuple(len(item) for item in values), np.inf, dtype=np.float64)
        selected_local = tuple(
            axis_id for axis_id in selection if axis_id not in shared_ids
        )
        local_shape = tuple(len(ordered_axes[axis_id]) for axis_id in selected_local)
        shared_domains = tuple(
            tuple(enumerate(ordered_axes[axis_id])) for axis_id in shared_ids
        )
        for indexed_shared in product(*shared_domains):
            shared_indices = {
                axis_id: indexed[0]
                for axis_id, indexed in zip(shared_ids, indexed_shared, strict=True)
            }
            total = np.zeros(local_shape, dtype=np.float64)
            for result, table in zip(factor_results, factor_tables, strict=True):
                remaining_ids: list[str] = []
                indexer: list[int | slice] = []
                for axis_id in result.factor.grid_ids:
                    if axis_id in shared_indices:
                        indexer.append(shared_indices[axis_id])
                    else:
                        indexer.append(slice(None))
                        remaining_ids.append(axis_id)
                contribution = np.asarray(table[tuple(indexer)], dtype=np.float64)
                reduce_axes = tuple(
                    index
                    for index, axis_id in enumerate(remaining_ids)
                    if axis_id not in selected_local
                )
                if reduce_axes:
                    contribution = np.min(contribution, axis=reduce_axes)
                retained_ids = tuple(
                    axis_id for axis_id in remaining_ids if axis_id in selected_local
                )
                reshape = tuple(
                    len(ordered_axes[axis_id]) if axis_id in retained_ids else 1
                    for axis_id in selected_local
                )
                total += contribution.reshape(reshape)
            target = tuple(
                shared_indices[axis_id] if axis_id in shared_indices else slice(None)
                for axis_id in selection
            )
            chisqr[target] = np.minimum(chisqr[target], total)
        return ProfiledGridSurface(selection, values, chisqr)

    profiles_1d = {axis_id: surface((axis_id,)) for axis_id in axis_ids}
    profiles_2d = {
        selection: surface(selection)
        for selection in _profile_surface_pairs(
            axis_ids,
            tuple(result.factor for result in factor_results),
        )
    }
    return ProfiledGridAggregate(
        ProfiledGridSelection(
            tuple((param_id, selected_grid[param_id]) for param_id in axis_ids),
            tuple(nuisance.items()),
            selected_chi_square,
            tuple(point.ordinal for point in selected_points),
        ),
        profiles_1d,
        profiles_2d,
    )
