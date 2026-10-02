"""Bounded, non-authoritative factor-local preparation of one TRF start."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from chemex.evaluation.native import EvaluationEngine
from chemex.optimize.direct_trf import (
    CancellationToken,
    DirectTrfConstructionError,
    MaterializedDirectTrfCandidate,
    OptimizationProblem,
    materialize_root_candidate,
)
from chemex.optimize.profiled_grid import (
    ProfiledGridConstructionError,
    ProfiledGridFactor,
    ProfiledGridPoint,
    ProfiledGridPointStatus,
    _fit_factor_point,
    build_profiled_factors,
)
from chemex.parameters.parameterization import ActiveParameterization


def preparation_request_budget(nvars: int, *, alternate: bool = False) -> int:
    """Limit actual objective requests separately from authoritative refinement."""
    return min(600 * (nvars + 1), 500 if alternate else 3000)


def mirrored_endpoints(
    endpoint: tuple[tuple[str, float], ...], mirror_ids: tuple[str, ...]
) -> tuple[tuple[tuple[str, float], ...], ...]:
    """Enumerate signs of nonzero fitted magnitudes inside one proven factor."""
    selected = tuple(key for key, _ in endpoint if key in mirror_ids)
    if len(selected) > 2:
        raise DirectTrfConstructionError(
            "TRY_DW_SIGNS supports at most two selected DW parameters within the same coupled fit; the sign search was skipped"
        )
    active = tuple(key for key, value in endpoint if key in selected and value != 0.0)
    values = dict(endpoint)
    starts = []
    for signs in product((-1.0, 1.0), repeat=len(active)):
        updates = {
            key: sign * abs(values[key])
            for key, sign in zip(active, signs, strict=True)
        }
        start = tuple((key, updates.get(key, value)) for key, value in endpoint)
        if start != endpoint:
            starts.append(start)
    return tuple(starts)


@dataclass(frozen=True, slots=True)
class ProfilePreparationResult:
    """A validated start and unsuccessful-attempt evidence; never an accepted fit."""

    vector: tuple[float, ...]
    attempts: tuple[ProfiledGridPoint, ...]
    fallback_factors: tuple[int, ...]
    root_fallback: bool = False

    @property
    def incomplete(self) -> bool:
        return self.root_fallback or any(
            point.status is not ProfiledGridPointStatus.SUCCESS
            for point in self.attempts
        )


def prepare_profile_start(
    problem: OptimizationProblem,
    hold_ids: tuple[str, ...],
    mirror_ids: tuple[str, ...],
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> ProfilePreparationResult:
    """Prepare one state, keeping original factor values when nuisance TRF fails.

    No local endpoint carries acceptance or commit authority. The caller must
    perform the ordinary complete/grouped TRF, including normal acceptance.
    """
    problem.validate_parameterization(parameterization)
    if (
        not problem.acceptance_authority
        or engine.plan.identity != problem.evaluation_plan_identity
        or not {*hold_ids, *mirror_ids}.issubset(problem.controlled_ids)
        or set(hold_ids) & set(mirror_ids)
    ):
        raise DirectTrfConstructionError(
            "Preparation requires its complete root objective and disjoint FIT subsets"
        )
    token = cancellation or CancellationToken()
    values = dict(zip(problem.controlled_ids, problem.start, strict=True))
    original = dict(values)
    attempts: list[ProfiledGridPoint] = []
    fallback: list[int] = []
    try:
        factors = build_profiled_factors(problem, hold_ids, parameterization, engine)
    except ProfiledGridConstructionError:
        return ProfilePreparationResult(problem.start, (), (), True)
    for factor in factors:
        if not factor.nuisance_ids:
            continue
        child_engine = engine.project_profiles(factor.profile_indices)
        held = tuple((key, original[key]) for key in factor.grid_ids)

        best = _fit_preparation_factor(
            problem,
            factor,
            child_engine,
            parameterization,
            held,
            tuple((key, original[key]) for key in factor.nuisance_ids),
            token,
            attempts,
            len(factors),
        )
        if best.status is not ProfiledGridPointStatus.SUCCESS:
            fallback.append(factor.ordinal)
            continue
        # Always branch from the normal successful endpoint, never from an
        # earlier alternate or from numerical safety bounds.
        for start in _factor_branch_starts(best, mirror_ids, attempts):
            point = _fit_preparation_factor(
                problem,
                factor,
                child_engine,
                parameterization,
                held,
                start,
                token,
                attempts,
                len(factors),
                alternate=True,
            )
            if (
                point.status is ProfiledGridPointStatus.SUCCESS
                and point.chi_square is not None
                and best.chi_square is not None
                and point.chi_square < best.chi_square
            ):
                best = point
        values.update(best.nuisance_items)
    vector = tuple(values[key] for key in problem.controlled_ids)
    materialized = materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=vector,
        invocation_identity="profile-preparation",
        execution_identity=problem.identity,
        cancellation=token,
    )
    if token.is_cancelled:
        raise KeyboardInterrupt("Profile preparation interrupted")
    if isinstance(materialized, MaterializedDirectTrfCandidate):
        return ProfilePreparationResult(vector, tuple(attempts), tuple(fallback))
    # Optional work must never suppress the existing authoritative fit path.
    return ProfilePreparationResult(
        problem.start, tuple(attempts), tuple(fallback), True
    )


def _fit_preparation_factor(
    problem: OptimizationProblem,
    factor: ProfiledGridFactor,
    engine: EvaluationEngine,
    parameterization: ActiveParameterization,
    held: tuple[tuple[str, float], ...],
    start: tuple[tuple[str, float], ...],
    token: CancellationToken,
    attempts: list[ProfiledGridPoint],
    factor_count: int,
    *,
    alternate: bool = False,
) -> ProfiledGridPoint:
    ordinal = len(attempts)
    try:
        point = _fit_factor_point(
            problem,
            factor,
            engine,
            parameterization,
            ordinal,
            held,
            objective_request_budget=preparation_request_budget(
                len(factor.nuisance_ids), alternate=alternate
            ),
            cancellation=token,
            progress_observer=None,
            factor_count=factor_count,
            point_count=1,
            nuisance_start_items=start,
        )
    except DirectTrfConstructionError as error:
        point = ProfiledGridPoint(
            ordinal, held, ProfiledGridPointStatus.FAILED, failure=str(error)
        )
    if point.status in {
        ProfiledGridPointStatus.CANCELLED,
        ProfiledGridPointStatus.INTERRUPTED,
    }:
        raise KeyboardInterrupt("Profile preparation interrupted")
    attempts.append(point)
    return point


def _factor_branch_starts(
    normal: ProfiledGridPoint,
    mirror_ids: tuple[str, ...],
    attempts: list[ProfiledGridPoint],
) -> tuple[tuple[tuple[str, float], ...], ...]:
    """Reject unqualified combinatorics without discarding a valid normal fit."""
    try:
        return mirrored_endpoints(normal.nuisance_items, mirror_ids)
    except DirectTrfConstructionError as error:
        attempts.append(
            ProfiledGridPoint(
                len(attempts),
                normal.axis_items,
                ProfiledGridPointStatus.FAILED,
                failure=str(error),
            )
        )
        return ()
