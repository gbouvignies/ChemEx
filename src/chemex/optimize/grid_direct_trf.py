"""Single-component Cartesian GRID to native Direct TRF qualification (#595).

This module is intentionally not wired into production dispatch.  It plans
physical-coordinate Cartesian seeds rooted in one immutable native problem;
execution and selection are added at the same public qualification seam.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import product
from numbers import Real
from typing import Protocol, cast

from chemex.evaluation.native import EvaluationEngine
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CandidateMaterialization,
    CommitReceipt,
    DirectTrfCandidateOutcome,
    DirectTrfCandidateTerminal,
    DirectTrfConstructionError,
    DirectTrfInterrupted,
    DirectTrfInvocation,
    DirectTrfTerminal,
    GridSeedProblemDerivation,
    GridSelectionProvenance,
    LiveFitCommitAuthority,
    MaterializationTerminal,
    MaterializedDirectTrfCandidate,
    OptimizationProblem,
    TerminalFailure,
    _accept_materialized_fit_for_derived_workflow,
    _grant_derived_fit_commit_authority,
    _materialize_derived_direct_trf_candidate_for_root,
    canonical_chi_square,
    commit_accepted_fit,
    execute_direct_trf_candidate,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.values import AnalysisValues

_SCHEMA_VERSION = 1
_GRID_WORKFLOW_VERSION = "native-cartesian-grid-direct-trf-v1"
_GRID_CANDIDATE_ORDER_VERSION = "chi-square-vector-seed-ordinal-v1"


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _float_token(value: float) -> str:
    return float(value).hex()


@dataclass(frozen=True, slots=True)
class UnusableGridValue:
    """Deterministic evidence for one seed-local unusable axis value."""

    category: str
    type_identity: str
    value_token: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-unusable-value",
                (self.category, self.type_identity, self.value_token),
            ),
        )


type GridCoordinate = float | UnusableGridValue


def _grid_coordinate(value: object) -> GridCoordinate:
    if isinstance(value, UnusableGridValue):
        return value
    type_identity = f"{type(value).__module__}.{type(value).__qualname__}"
    if not isinstance(value, bool) and isinstance(value, Real):
        try:
            result = float(value)
        except (OverflowError, TypeError, ValueError):
            return UnusableGridValue("non_real", type_identity, type_identity)
        if math.isfinite(result):
            return 0.0 if result == 0.0 else result
        category = (
            "nan"
            if math.isnan(result)
            else "positive_infinity"
            if result > 0.0
            else "negative_infinity"
        )
        return UnusableGridValue(category, type_identity, result.hex())
    if isinstance(value, str):
        canonical_value = value
    elif isinstance(value, complex):
        canonical_value = f"{float(value.real).hex()}:{float(value.imag).hex()}"
    elif value is None:
        canonical_value = "<none>"
    else:
        try:
            canonical_value = json.dumps(
                value,
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
        except (TypeError, ValueError):
            canonical_value = type_identity
    return UnusableGridValue("non_real", type_identity, canonical_value)


def _coordinate_token(value: GridCoordinate) -> tuple[str, str]:
    if isinstance(value, UnusableGridValue):
        return "unusable", value.identity
    return "finite", _float_token(value)


@dataclass(frozen=True, slots=True)
class GridAxis:
    """One declared physical-coordinate axis with canonical binary64 values."""

    param_id: str
    values: tuple[GridCoordinate, ...]
    declaration_ordinal: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.param_id:
            raise DirectTrfConstructionError("GRID axis parameter ID cannot be empty")
        if isinstance(self.declaration_ordinal, bool) or self.declaration_ordinal < 0:
            raise DirectTrfConstructionError(
                "GRID declaration ordinal must be a non-negative integer"
            )
        values = tuple(_grid_coordinate(value) for value in self.values)
        if not values:
            raise DirectTrfConstructionError("GRID axes cannot be empty")
        object.__setattr__(self, "values", values)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-axis",
                (
                    self.param_id,
                    tuple(_coordinate_token(value) for value in values),
                    self.declaration_ordinal,
                ),
            ),
        )


def _seed_identity(
    root_problem_identity: str,
    ordinal: int,
    axis_items: tuple[tuple[str, GridCoordinate], ...],
    start: tuple[GridCoordinate, ...],
) -> str:
    return _identity(
        "native-grid-seed",
        (
            root_problem_identity,
            ordinal,
            tuple(
                (param_id, _coordinate_token(value)) for param_id, value in axis_items
            ),
            tuple(_coordinate_token(value) for value in start),
        ),
    )


@dataclass(frozen=True, slots=True)
class GridSeed:
    """One canonical seed, including typed rejection before local TRF."""

    root_problem_identity: str
    ordinal: int
    axis_items: tuple[tuple[str, GridCoordinate], ...]
    start: tuple[GridCoordinate, ...]
    problem: OptimizationProblem | None = field(repr=False, compare=False)
    invocation: DirectTrfInvocation | None = field(repr=False, compare=False)
    rejection: TerminalFailure | None = None
    coordinate_identity: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        planned = self.problem is not None and self.invocation is not None
        if planned == (self.rejection is not None):
            raise ValueError(
                "A GRID seed must have either one local attempt or one rejection"
            )
        coordinate_identity = _seed_identity(
            self.root_problem_identity,
            self.ordinal,
            self.axis_items,
            self.start,
        )
        if planned:
            problem = self.problem
            invocation = self.invocation
            derivation = problem.derivation
            if (
                not isinstance(derivation, GridSeedProblemDerivation)
                or derivation.root_problem_identity != self.root_problem_identity
                or derivation.seed_identity != coordinate_identity
                or derivation.seed_ordinal != self.ordinal
                or derivation.axis_items != self.axis_items
                or derivation.start != self.start
                or invocation.problem_identity != problem.identity
            ):
                raise DirectTrfConstructionError(
                    "GRID seed child problem or invocation differs from its seed"
                )
        object.__setattr__(self, "coordinate_identity", coordinate_identity)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-seed-plan",
                (
                    coordinate_identity,
                    None if self.problem is None else self.problem.identity,
                    None if self.invocation is None else self.invocation.identity,
                    None if self.rejection is None else self.rejection.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class GridDirectTrfInvocation:
    """Immutable Cartesian seed plan and uniform local Direct TRF policy."""

    root_problem_identity: str
    axes: tuple[GridAxis, ...]
    seeds: tuple[GridSeed, ...]
    objective_request_budget: int
    x_scale: tuple[float, ...]
    ftol: float | None
    xtol: float | None
    gtol: float | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.axes or not self.seeds:
            raise DirectTrfConstructionError("GRID requires at least one seed")
        if any(
            seed.root_problem_identity != self.root_problem_identity
            for seed in self.seeds
        ):
            raise DirectTrfConstructionError(
                "GRID seeds belong to another root problem"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-direct-trf-invocation",
                (
                    _GRID_WORKFLOW_VERSION,
                    self.root_problem_identity,
                    tuple(axis.identity for axis in self.axes),
                    tuple(seed.identity for seed in self.seeds),
                    self.objective_request_budget,
                    tuple(_float_token(value) for value in self.x_scale),
                    tuple(
                        None if value is None else _float_token(value)
                        for value in (self.ftol, self.xtol, self.gtol)
                    ),
                ),
            ),
        )

    @classmethod
    def for_problem(
        cls,
        problem: OptimizationProblem,
        *,
        axes: Sequence[tuple[str, Sequence[object]]],
        objective_request_budget: int,
        x_scale: float | Sequence[float] | None = None,
        ftol: float | None = 1.0e-8,
        xtol: float | None = 1.0e-8,
        gtol: float | None = 1.0e-8,
    ) -> GridDirectTrfInvocation:
        """Plan stable-ID axes and lexicographic seeds in physical coordinates."""
        if not problem.acceptance_authority:
            raise DirectTrfConstructionError("GRID requires a complete root problem")
        declared = tuple(
            GridAxis(
                param_id,
                tuple(_grid_coordinate(value) for value in values),
                ordinal,
            )
            for ordinal, (param_id, values) in enumerate(axes)
        )
        if not declared:
            raise DirectTrfConstructionError("GRID requires at least one axis")
        axis_ids = tuple(axis.param_id for axis in declared)
        if len(set(axis_ids)) != len(axis_ids):
            raise DirectTrfConstructionError("Duplicate GRID axis parameter ID")
        unknown = set(axis_ids) - set(problem.controlled_ids)
        if unknown:
            raise DirectTrfConstructionError(
                "GRID axes must target controlled external coordinates: "
                + ", ".join(sorted(unknown))
            )
        canonical_axes = tuple(
            sorted(declared, key=lambda axis: axis.param_id.encode("utf-8"))
        )
        template = DirectTrfInvocation.for_problem(
            problem,
            objective_request_budget=objective_request_budget,
            x_scale=x_scale,
            ftol=ftol,
            xtol=xtol,
            gtol=gtol,
        )
        coordinate_indices = {
            param_id: index for index, param_id in enumerate(problem.controlled_ids)
        }
        seeds: list[GridSeed] = []
        for ordinal, values in enumerate(
            product(*(axis.values for axis in canonical_axes))
        ):
            axis_items = tuple(
                (axis.param_id, value)
                for axis, value in zip(canonical_axes, values, strict=True)
            )
            start_list: list[GridCoordinate] = [float(value) for value in problem.start]
            for param_id, value in axis_items:
                start_list[coordinate_indices[param_id]] = value
            start = tuple(start_list)
            seed_identity = _seed_identity(
                problem.identity,
                ordinal,
                axis_items,
                start,
            )
            unusable = tuple(
                param_id
                for param_id, value in axis_items
                if isinstance(value, UnusableGridValue)
            )
            if unusable:
                seeds.append(
                    GridSeed(
                        problem.identity,
                        ordinal,
                        axis_items,
                        start,
                        None,
                        None,
                        TerminalFailure(
                            "grid_seed_unusable_value",
                            "GRID seed has an unusable physical coordinate: "
                            + ", ".join(unusable),
                        ),
                    )
                )
                continue
            finite_axis_items = cast(
                "tuple[tuple[str, float], ...]",
                axis_items,
            )
            finite_start = cast("tuple[float, ...]", start)
            outside = tuple(
                param_id
                for param_id, value in finite_axis_items
                if not (
                    problem.lower_bounds[coordinate_indices[param_id]]
                    <= value
                    <= problem.upper_bounds[coordinate_indices[param_id]]
                )
            )
            if outside:
                seeds.append(
                    GridSeed(
                        problem.identity,
                        ordinal,
                        axis_items,
                        start,
                        None,
                        None,
                        TerminalFailure(
                            "grid_seed_out_of_bounds",
                            "GRID seed is outside physical bounds: "
                            + ", ".join(outside),
                        ),
                    )
                )
                continue
            derivation = GridSeedProblemDerivation(
                problem.identity,
                seed_identity,
                ordinal,
                finite_axis_items,
                problem.controlled_ids,
                finite_start,
            )
            child = OptimizationProblem(
                problem.evaluation_plan_identity,
                problem.parameterization_identity,
                problem.evaluator_parameterization_identity,
                problem.constraint_program_identity,
                problem.configuration_identity,
                problem.source_snapshot,
                problem.independent_items,
                problem.controlled_ids,
                problem.held_items,
                finite_start,
                problem.lower_bounds,
                problem.upper_bounds,
                problem.commit_scope,
                derivation,
                problem.scalarization_version,
            )
            child_invocation = DirectTrfInvocation.for_problem(
                child,
                objective_request_budget=template.objective_request_budget,
                x_scale=template.x_scale,
                ftol=template.ftol,
                xtol=template.xtol,
                gtol=template.gtol,
            )
            seeds.append(
                GridSeed(
                    problem.identity,
                    ordinal,
                    axis_items,
                    start,
                    child,
                    child_invocation,
                )
            )
        return cls(
            problem.identity,
            canonical_axes,
            tuple(seeds),
            template.objective_request_budget,
            template.x_scale,
            template.ftol,
            template.xtol,
            template.gtol,
        )


class GridSeedDisposition(StrEnum):
    """Closed disposition for every canonical seed occurrence."""

    ELIGIBLE = "eligible"
    UNUSABLE_VALUE = "unusable_value"
    OUT_OF_BOUNDS = "out_of_bounds"
    PREFLIGHT_INVALID = "preflight_invalid"
    INVALID_TRIAL = "invalid_trial"
    BUDGET_EXHAUSTED = "budget_exhausted"
    NON_CONVERGED = "non_converged"
    BACKEND_FAILURE = "backend_failure"
    IMPLEMENTATION_FAILURE = "implementation_failure"
    MATERIALIZATION_FAILURE = "materialization_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    NOT_STARTED = "not_started"


@dataclass(frozen=True, slots=True)
class GridCandidate:
    """One workflow-local eligible candidate under the canonical total order."""

    seed_identity: str
    seed_ordinal: int
    candidate: MaterializedDirectTrfCandidate = field(repr=False, compare=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-candidate",
                (
                    _GRID_CANDIDATE_ORDER_VERSION,
                    self.seed_identity,
                    self.seed_ordinal,
                    self.candidate.identity,
                ),
            ),
        )

    @property
    def vector(self) -> tuple[float, ...]:
        return self.candidate.vector

    @property
    def objective(self) -> float:
        return self.candidate.chi_square

    def ordering_key(self) -> tuple[float, tuple[float, ...], int]:
        """Return χ², final external vector, then seed ordinal."""
        return self.objective, self.vector, self.seed_ordinal


class GridCandidateRecord(Protocol):
    """Narrow evidence contract consumed by canonical GRID finalization."""

    @property
    def seed_identity(self) -> str: ...

    @property
    def seed_ordinal(self) -> int: ...

    @property
    def axis_items(self) -> tuple[tuple[str, GridCoordinate], ...]: ...

    @property
    def start(self) -> tuple[GridCoordinate, ...]: ...

    @property
    def disposition(self) -> GridSeedDisposition: ...

    @property
    def execution_identity(self) -> str | None: ...

    @property
    def objective(self) -> float | None: ...

    @property
    def candidate(self) -> GridCandidate | None: ...

    @property
    def failure(self) -> TerminalFailure | None: ...

    @property
    def identity(self) -> str: ...

    def validate_for_grid_selection(
        self,
        problem: OptimizationProblem,
        invocation: GridDirectTrfInvocation,
        seed: GridSeed,
    ) -> None:
        """Raise when this candidate is not exact evidence for ``seed``."""
        ...


@dataclass(frozen=True, slots=True)
class GridSeedOutcome:
    """Compact exact table/provenance row for one canonical seed."""

    seed_identity: str
    seed_ordinal: int
    axis_items: tuple[tuple[str, GridCoordinate], ...]
    start: tuple[GridCoordinate, ...]
    disposition: GridSeedDisposition
    execution_identity: str | None = None
    objective: float | None = None
    candidate: GridCandidate | None = field(default=None, repr=False, compare=False)
    direct_outcome: DirectTrfCandidateOutcome | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    failure: TerminalFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        eligible = self.disposition is GridSeedDisposition.ELIGIBLE
        candidate = self.candidate
        if eligible != (self.candidate is not None and self.objective is not None):
            raise ValueError("Only an eligible GRID seed may expose a candidate")
        if candidate is not None and (
            candidate.seed_identity != self.seed_identity
            or candidate.seed_ordinal != self.seed_ordinal
            or candidate.objective != self.objective
        ):
            raise ValueError("GRID seed candidate provenance is inconsistent")
        if eligible and (
            candidate is None
            or self.direct_outcome is None
            or self.direct_outcome.terminal is not DirectTrfCandidateTerminal.SUCCESS
            or self.direct_outcome.candidate is not candidate.candidate
            or self.execution_identity != self.direct_outcome.execution.identity
        ):
            raise ValueError(
                "Eligible GRID seed lacks its exact Direct TRF candidate result"
            )
        if self.disposition is GridSeedDisposition.NOT_STARTED and (
            self.execution_identity is not None or self.direct_outcome is not None
        ):
            raise ValueError("An unstarted GRID seed cannot expose attempt evidence")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-seed-outcome",
                (
                    self.seed_identity,
                    self.seed_ordinal,
                    tuple(
                        (param_id, _coordinate_token(value))
                        for param_id, value in self.axis_items
                    ),
                    tuple(_coordinate_token(value) for value in self.start),
                    self.disposition.value,
                    self.execution_identity,
                    None if self.objective is None else _float_token(self.objective),
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
        """Prove exact single-component seed lineage before selection."""
        _validate_direct_grid_record(self, problem, invocation, seed)


def _validate_direct_grid_record(
    record: GridSeedOutcome,
    problem: OptimizationProblem,
    invocation: GridDirectTrfInvocation,
    seed: GridSeed,
) -> None:
    """Validate #595 child/candidate evidence at GRID's ownership boundary."""
    candidate = record.candidate
    direct = record.direct_outcome
    child_problem = seed.problem
    child_invocation = seed.invocation
    materialized = None if candidate is None else candidate.candidate
    materialization = None if materialized is None else materialized.materialization
    derivation = None if child_problem is None else child_problem.derivation
    unchanged_root_semantics = (
        child_problem is not None
        and child_problem.evaluation_plan_identity == problem.evaluation_plan_identity
        and child_problem.parameterization_identity == problem.parameterization_identity
        and child_problem.evaluator_parameterization_identity
        == problem.evaluator_parameterization_identity
        and child_problem.constraint_program_identity
        == problem.constraint_program_identity
        and child_problem.configuration_identity == problem.configuration_identity
        and child_problem.source_snapshot == problem.source_snapshot
        and child_problem.independent_items == problem.independent_items
        and child_problem.controlled_ids == problem.controlled_ids
        and child_problem.held_items == problem.held_items
        and child_problem.lower_bounds == problem.lower_bounds
        and child_problem.upper_bounds == problem.upper_bounds
        and child_problem.commit_scope == problem.commit_scope
        and child_problem.scalarization_version == problem.scalarization_version
    )
    if (
        invocation.root_problem_identity != problem.identity
        or seed.root_problem_identity != problem.identity
        or seed.identity != record.seed_identity
        or seed.ordinal != record.seed_ordinal
        or record.disposition is not GridSeedDisposition.ELIGIBLE
        or candidate is None
        or candidate.seed_identity != seed.identity
        or candidate.seed_ordinal != seed.ordinal
        or direct is None
        or materialized is None
        or direct.terminal is not DirectTrfCandidateTerminal.SUCCESS
        or direct.execution.terminal is not DirectTrfTerminal.CONVERGED
        or direct.candidate is not materialized
        or direct.materialization is not materialization
        or record.execution_identity != direct.execution.identity
        or child_problem is None
        or child_invocation is None
        or not isinstance(derivation, GridSeedProblemDerivation)
        or derivation.root_problem_identity != problem.identity
        or derivation.seed_identity != seed.coordinate_identity
        or derivation.seed_ordinal != seed.ordinal
        or not unchanged_root_semantics
        or direct.execution.problem_identity != child_problem.identity
        or direct.execution.invocation_identity != child_invocation.identity
        or not _materialized_grid_candidate_matches(
            materialized,
            child_problem,
            child_invocation.identity,
            direct.execution.identity,
        )
    ):
        raise DirectTrfConstructionError(
            "GRID candidate lacks exact canonical seed lineage"
        )


def _materialized_grid_candidate_matches(
    candidate: MaterializedDirectTrfCandidate,
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
) -> bool:
    """Check self-consistent materialized evidence for a known GRID problem."""
    evaluation = candidate.evaluation_result
    materialization = candidate.materialization
    try:
        objective = canonical_chi_square(evaluation.residuals)
        vector = tuple(
            evaluation.resolved_values[param_id] for param_id in problem.controlled_ids
        )
    except Exception:  # noqa: BLE001 - malformed evidence fails closed
        return False
    return (
        candidate.problem_identity == problem.identity
        and candidate.invocation_identity == invocation_identity
        and candidate.execution_identity == execution_identity
        and candidate.vector == vector
        and candidate.chi_square == objective
        and materialization.problem_identity == problem.identity
        and materialization.invocation_identity == invocation_identity
        and materialization.execution_identity == execution_identity
        and materialization.terminal is MaterializationTerminal.SUCCESS
        and materialization.candidate.vector == candidate.vector
        and materialization.candidate.chi_square == candidate.chi_square
        and materialization.evaluation_identity == evaluation.identity
        and materialization.evaluator_compatibility_identity
        == evaluation.evaluator_compatibility_identity
        and evaluation.plan_identity == problem.evaluation_plan_identity
        and evaluation.parameterization_identity
        == problem.evaluator_parameterization_identity
        and tuple(evaluation.resolved_values) == problem.commit_scope
    )


@dataclass(frozen=True, slots=True)
class GridSelection:
    """Deterministic winner and complete eligible-order provenance."""

    selected_seed_identity: str
    selected_seed_ordinal: int
    selected_candidate_identity: str
    eligible_candidate_identities: tuple[str, ...]
    candidate_records: tuple[GridCandidateRecord, ...] = field(
        repr=False,
        compare=False,
    )
    ordering_policy: str = _GRID_CANDIDATE_ORDER_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        records = self.candidate_records
        record_candidate_identities = tuple(
            record.candidate.identity
            for record in records
            if record.candidate is not None
        )
        selected_records = tuple(
            record
            for record in records
            if record.seed_identity == self.selected_seed_identity
            and record.seed_ordinal == self.selected_seed_ordinal
            and record.candidate is not None
            and record.candidate.identity == self.selected_candidate_identity
        )
        if (
            not self.eligible_candidate_identities
            or self.eligible_candidate_identities[0] != self.selected_candidate_identity
            or not records
            or any(
                record.disposition is not GridSeedDisposition.ELIGIBLE
                or record.candidate is None
                for record in records
            )
            or record_candidate_identities != self.eligible_candidate_identities
            or len(selected_records) != 1
            or selected_records[0].identity != records[0].identity
        ):
            raise ValueError(
                "GRID selection must lead with one exact canonical candidate record"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-selection",
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
        *,
        selected_seed_identity: str,
        selected_seed_ordinal: int,
        selected_candidate_identity: str,
        candidate_records: Sequence[GridCandidateRecord],
    ) -> GridSelection:
        """Validate the requested winner against canonical eligible records."""
        records = tuple(candidate_records)
        if not records or any(
            record.disposition is not GridSeedDisposition.ELIGIBLE
            or record.candidate is None
            for record in records
        ):
            raise ValueError(
                "GRID selection requires exact canonical candidate records"
            )
        ordered = tuple(
            sorted(
                records,
                key=lambda record: cast(
                    "GridCandidate", record.candidate
                ).ordering_key(),
            )
        )
        return cls(
            selected_seed_identity,
            selected_seed_ordinal,
            selected_candidate_identity,
            tuple(
                cast("GridCandidate", record.candidate).identity for record in ordered
            ),
            ordered,
        )

    @property
    def selected_record(self) -> GridCandidateRecord:
        """Return the one exact candidate record that won selection."""
        return self.candidate_records[0]


@dataclass(frozen=True, slots=True)
class AcceptedGridDirectTrfResult(AcceptedFitResult):
    """Convenient immutable evidence for one accepted GRID selection."""

    grid_selection_provenance: GridSelectionProvenance
    workflow_invocation: GridDirectTrfInvocation = field(repr=False, compare=False)
    selection: GridSelection = field(repr=False, compare=False)
    selected_seed: GridSeed = field(repr=False, compare=False)
    selected_outcome: GridCandidateRecord = field(repr=False, compare=False)
    fresh_candidate: MaterializedDirectTrfCandidate = field(
        repr=False,
        compare=False,
    )

    @classmethod
    def from_accepted(
        cls,
        accepted: AcceptedFitResult,
        provenance: GridSelectionProvenance,
        workflow_invocation: GridDirectTrfInvocation,
        selection: GridSelection,
        selected_seed: GridSeed,
        selected_outcome: GridCandidateRecord,
        fresh_candidate: MaterializedDirectTrfCandidate,
    ) -> AcceptedGridDirectTrfResult:
        """Attach GRID provenance to already accepted root evidence."""
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
            selection,
            selected_seed,
            selected_outcome,
            fresh_candidate,
        )


def validate_grid_acceptance_lineage(
    provenance: GridSelectionProvenance,
    invocation: GridDirectTrfInvocation,
    selection: GridSelection,
    selected_seed: GridSeed,
    selected_outcome: GridCandidateRecord,
    fresh_candidate: MaterializedDirectTrfCandidate,
    problem: OptimizationProblem,
) -> None:
    """Validate canonical selection and fresh root promotion for any GRID record."""
    candidate = selected_outcome.candidate
    try:
        canonical_seed = invocation.seeds[provenance.selected_seed_ordinal]
    except (AttributeError, IndexError, TypeError):
        canonical_seed = None
    materialized_candidate = None if candidate is None else candidate.candidate
    candidate_materialization = (
        None
        if materialized_candidate is None
        else materialized_candidate.materialization
    )
    fresh_materialization = fresh_candidate.materialization
    try:
        if canonical_seed is not None:
            selected_outcome.validate_for_grid_selection(
                problem,
                invocation,
                canonical_seed,
            )
    except Exception as error:
        raise DirectTrfConstructionError(
            "Accepted GRID result lacks its canonical selection authority"
        ) from error
    if (
        not isinstance(provenance, GridSelectionProvenance)
        or invocation.root_problem_identity != problem.identity
        or provenance.root_problem_identity != problem.identity
        or provenance.workflow_invocation_identity != invocation.identity
        or provenance.selection_identity != selection.identity
        or canonical_seed is None
        or canonical_seed.identity != selected_seed.identity
        or canonical_seed.ordinal != provenance.selected_seed_ordinal
        or provenance.selected_seed_identity != selected_seed.identity
        or selection.selected_seed_identity != selected_seed.identity
        or selection.selected_seed_ordinal != selected_seed.ordinal
        or selection.selected_record.identity != selected_outcome.identity
        or selected_outcome.seed_identity != selected_seed.identity
        or selected_outcome.seed_ordinal != selected_seed.ordinal
        or candidate is None
        or materialized_candidate is None
        or selection.selected_candidate_identity != candidate.identity
        or provenance.grid_candidate_identity != candidate.identity
        or provenance.materialized_candidate_identity != materialized_candidate.identity
        or provenance.candidate_problem_identity
        != materialized_candidate.problem_identity
        or provenance.candidate_invocation_identity
        != materialized_candidate.invocation_identity
        or provenance.candidate_execution_identity
        != materialized_candidate.execution_identity
        or candidate_materialization is None
        or provenance.candidate_materialization_identity
        != candidate_materialization.identity
        or candidate_materialization.problem_identity
        != materialized_candidate.problem_identity
        or candidate_materialization.invocation_identity
        != materialized_candidate.invocation_identity
        or candidate_materialization.execution_identity
        != materialized_candidate.execution_identity
        or candidate_materialization.evaluation_identity
        != materialized_candidate.evaluation_result.identity
        or provenance.accepted_materialization_identity
        != fresh_materialization.identity
        or fresh_candidate.problem_identity != problem.identity
        or fresh_candidate.invocation_identity
        != materialized_candidate.invocation_identity
        or fresh_candidate.execution_identity
        != materialized_candidate.execution_identity
        or fresh_materialization.problem_identity != problem.identity
        or fresh_materialization.invocation_identity
        != materialized_candidate.invocation_identity
        or fresh_materialization.execution_identity
        != materialized_candidate.execution_identity
        or fresh_materialization.terminal is not MaterializationTerminal.SUCCESS
        or fresh_materialization.candidate.vector != fresh_candidate.vector
        or fresh_materialization.candidate.chi_square != fresh_candidate.chi_square
        or fresh_materialization.evaluation_identity
        != fresh_candidate.evaluation_result.identity
        or fresh_materialization.evaluator_compatibility_identity
        != fresh_candidate.evaluation_result.evaluator_compatibility_identity
        or provenance.accepted_evaluation_identity
        != fresh_candidate.evaluation_result.identity
        or fresh_candidate.evaluation_result.plan_identity
        != problem.evaluation_plan_identity
        or fresh_candidate.evaluation_result.parameterization_identity
        != problem.evaluator_parameterization_identity
        or tuple(fresh_candidate.evaluation_result.resolved_values)
        != problem.commit_scope
        or fresh_candidate.vector != candidate.vector
        or fresh_candidate.chi_square != candidate.objective
    ):
        raise DirectTrfConstructionError(
            "Accepted GRID result lacks its canonical selection authority"
        )


def validate_grid_commit_lineage(
    accepted: AcceptedGridDirectTrfResult,
    problem: OptimizationProblem,
) -> None:
    """Revalidate GRID evidence independently of generic commit authority."""
    validate_grid_acceptance_lineage(
        accepted.grid_selection_provenance,
        accepted.workflow_invocation,
        accepted.selection,
        accepted.selected_seed,
        accepted.selected_outcome,
        accepted.fresh_candidate,
        problem,
    )
    fresh = accepted.fresh_candidate
    if (
        accepted.problem_identity != problem.identity
        or accepted.invocation_identity != fresh.invocation_identity
        or accepted.execution_identity != fresh.execution_identity
        or accepted.materialization_identity != fresh.materialization.identity
        or accepted.vector != fresh.vector
        or accepted.chi_square != fresh.chi_square
        or accepted.evaluation_result.identity != fresh.evaluation_result.identity
        or accepted.origin_context_identity
        != accepted.grid_selection_provenance.identity
    ):
        raise DirectTrfConstructionError(
            "Accepted GRID result lacks its canonical selection authority"
        )


def commit_grid_accepted_fit(
    accepted: AcceptedGridDirectTrfResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Validate GRID evidence, then cross the generic Direct TRF commit seam."""
    validate_grid_commit_lineage(accepted, problem)
    return commit_accepted_fit(
        accepted,
        authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=analysis_values,
    )


class GridDirectTrfTerminal(StrEnum):
    """Closed outcomes of one GRID to Direct TRF workflow occurrence."""

    ACCEPTED = "accepted"
    NO_ELIGIBLE_CANDIDATE = "no_eligible_candidate"
    MATERIALIZATION_FAILURE = "materialization_failure"
    EXECUTION_FAILURE = "execution_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


@dataclass(frozen=True, slots=True)
class GridTableRow:
    """Exact per-seed evidence retained for a future GRID table renderer."""

    seed_identity: str
    seed_ordinal: int
    axis_items: tuple[tuple[str, GridCoordinate], ...]
    start: tuple[GridCoordinate, ...]
    objective: float | None
    disposition: GridSeedDisposition
    candidate_identity: str | None
    selected: bool
    selection_identity: str | None


@dataclass(frozen=True, slots=True)
class GridDiagnosticPoint:
    """One exact seed projection for supported one- or two-axis diagnostics."""

    seed_identity: str
    seed_ordinal: int
    coordinates: tuple[tuple[str, GridCoordinate], ...]
    objective: float | None
    disposition: GridSeedDisposition
    candidate_identity: str | None
    selected: bool
    selection_identity: str | None


@dataclass(frozen=True, slots=True)
class GridDirectTrfOutcome:
    """Complete GRID occurrence; only ACCEPTED carries root fit authority."""

    terminal: GridDirectTrfTerminal
    attempts: tuple[GridCandidateRecord, ...]
    selection: GridSelection | None = None
    materialization: CandidateMaterialization | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    accepted_result: AcceptedGridDirectTrfResult | None = field(
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
        accepted = self.terminal is GridDirectTrfTerminal.ACCEPTED
        if accepted != (
            self.accepted_result is not None and self.commit_authority is not None
        ):
            raise ValueError(
                "Only accepted GRID execution has evidence and live authority"
            )
        if not accepted and (
            self.accepted_result is not None or self.commit_authority is not None
        ):
            raise ValueError("Unaccepted GRID execution cannot expose fit authority")
        if (
            self.terminal
            in {
                GridDirectTrfTerminal.ACCEPTED,
                GridDirectTrfTerminal.MATERIALIZATION_FAILURE,
            }
            and self.selection is None
        ):
            raise ValueError("Completed GRID acceptance requires a selection")
        if self.selection is not None:
            canonical_candidates = tuple(
                sorted(
                    (
                        attempt
                        for attempt in self.attempts
                        if attempt.disposition is GridSeedDisposition.ELIGIBLE
                        and attempt.candidate is not None
                    ),
                    key=lambda attempt: cast(
                        "GridCandidate", attempt.candidate
                    ).ordering_key(),
                )
            )
            if tuple(
                record.identity for record in self.selection.candidate_records
            ) != tuple(record.identity for record in canonical_candidates):
                raise ValueError("GRID selection provenance is inconsistent")

    @property
    def table_rows(self) -> tuple[GridTableRow, ...]:
        """Return canonical rows without re-evaluation or objective reduction."""
        selected_record_identity = (
            None if self.selection is None else self.selection.selected_record.identity
        )
        selection_identity = None if self.selection is None else self.selection.identity
        return tuple(
            GridTableRow(
                attempt.seed_identity,
                attempt.seed_ordinal,
                attempt.axis_items,
                attempt.start,
                attempt.objective,
                attempt.disposition,
                None if attempt.candidate is None else attempt.candidate.identity,
                attempt.identity == selected_record_identity,
                selection_identity,
            )
            for attempt in self.attempts
        )

    def diagnostic_points(
        self,
        axis_ids: Sequence[str],
    ) -> tuple[GridDiagnosticPoint, ...]:
        """Project retained rows onto one or two explicitly ordered GRID axes."""
        requested = tuple(axis_ids)
        available = {
            param_id
            for attempt in self.attempts
            for param_id, _value in attempt.axis_items
        }
        if (
            len(requested) not in {1, 2}
            or len(set(requested)) != len(requested)
            or not set(requested).issubset(available)
        ):
            raise DirectTrfConstructionError(
                "GRID diagnostics require one or two unique planned axis IDs"
            )
        return tuple(
            GridDiagnosticPoint(
                row.seed_identity,
                row.seed_ordinal,
                tuple(
                    (param_id, dict(row.axis_items)[param_id]) for param_id in requested
                ),
                row.objective,
                row.disposition,
                row.candidate_identity,
                row.selected,
                row.selection_identity,
            )
            for row in self.table_rows
        )


class GridDirectTrfInterrupted(KeyboardInterrupt):
    """Propagated interruption carrying the frozen GRID workflow outcome."""

    def __init__(self, outcome: GridDirectTrfOutcome) -> None:
        self.outcome = outcome
        super().__init__("Native GRID to Direct TRF was interrupted")


def _validate_grid_context(
    problem: OptimizationProblem,
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> None:
    problem.validate_parameterization(parameterization)
    declared_axes = tuple(
        (axis.param_id, axis.values)
        for axis in sorted(invocation.axes, key=lambda axis: axis.declaration_ordinal)
    )
    expected = GridDirectTrfInvocation.for_problem(
        problem,
        axes=declared_axes,
        objective_request_budget=invocation.objective_request_budget,
        x_scale=invocation.x_scale,
        ftol=invocation.ftol,
        xtol=invocation.xtol,
        gtol=invocation.gtol,
    )
    if (
        invocation.root_problem_identity != problem.identity
        or invocation.identity != expected.identity
        or engine.plan.identity != problem.evaluation_plan_identity
    ):
        raise DirectTrfConstructionError(
            "GRID invocation and evaluator must share one root problem"
        )


def _not_started(seed: GridSeed) -> GridSeedOutcome:
    return GridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.NOT_STARTED,
    )


def _attempt_failure(
    seed: GridSeed,
    disposition: GridSeedDisposition,
    outcome: DirectTrfCandidateOutcome,
) -> GridSeedOutcome:
    materialization_failure = (
        outcome.materialization.failure if outcome.materialization is not None else None
    )
    return GridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        disposition,
        outcome.execution.identity,
        direct_outcome=outcome,
        failure=materialization_failure or outcome.execution.failure,
    )


def _interrupted_attempt(
    seed: GridSeed,
    interrupted: DirectTrfInterrupted,
) -> GridSeedOutcome:
    return GridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.INTERRUPTED,
        interrupted.execution.identity,
        failure=(
            interrupted.materialization.failure
            if interrupted.materialization is not None
            else interrupted.execution.failure
        ),
    )


def _interrupted_outcome(
    attempts: Sequence[GridSeedOutcome],
    remaining: Sequence[GridSeed],
    selection: GridSelection | None = None,
    materialization: CandidateMaterialization | None = None,
) -> GridDirectTrfOutcome:
    return GridDirectTrfOutcome(
        GridDirectTrfTerminal.INTERRUPTED,
        (*attempts, *(_not_started(seed) for seed in remaining)),
        selection,
        materialization,
    )


def _disposition_for_unsuccessful(
    outcome: DirectTrfCandidateOutcome,
) -> GridSeedDisposition:
    if outcome.terminal is DirectTrfCandidateTerminal.MATERIALIZATION_FAILURE:
        return GridSeedDisposition.MATERIALIZATION_FAILURE
    return {
        DirectTrfTerminal.PREFLIGHT_INVALID: GridSeedDisposition.PREFLIGHT_INVALID,
        DirectTrfTerminal.INVALID_TRIAL: GridSeedDisposition.INVALID_TRIAL,
        DirectTrfTerminal.BUDGET_EXHAUSTED: GridSeedDisposition.BUDGET_EXHAUSTED,
        DirectTrfTerminal.NON_CONVERGED: GridSeedDisposition.NON_CONVERGED,
        DirectTrfTerminal.BACKEND_FAILURE: GridSeedDisposition.BACKEND_FAILURE,
        DirectTrfTerminal.IMPLEMENTATION_FAILURE: (
            GridSeedDisposition.IMPLEMENTATION_FAILURE
        ),
    }.get(
        outcome.execution.terminal,
        GridSeedDisposition.IMPLEMENTATION_FAILURE,
    )


def _selection(
    attempts: Sequence[GridCandidateRecord],
) -> tuple[GridSelection, GridCandidateRecord]:
    eligible = sorted(
        (
            attempt
            for attempt in attempts
            if attempt.disposition is GridSeedDisposition.ELIGIBLE
            and attempt.candidate is not None
        ),
        key=lambda attempt: cast("GridCandidate", attempt.candidate).ordering_key(),
    )
    selected = eligible[0]
    candidate = cast("GridCandidate", selected.candidate)
    selection = GridSelection.from_candidate_records(
        selected_seed_identity=selected.seed_identity,
        selected_seed_ordinal=selected.seed_ordinal,
        selected_candidate_identity=candidate.identity,
        candidate_records=eligible,
    )
    return selection, selected


def _grid_finalization_gate(
    cancellation: CancellationToken,
) -> GridDirectTrfTerminal | None:
    """Gate canonical selection against late cancellation or interruption."""
    try:
        return GridDirectTrfTerminal.CANCELLED if cancellation.is_cancelled else None
    except KeyboardInterrupt:
        return GridDirectTrfTerminal.INTERRUPTED


def _finalize_grid_candidates(
    problem: OptimizationProblem,
    invocation: GridDirectTrfInvocation,
    attempts: Sequence[GridCandidateRecord],
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    cancellation: CancellationToken,
) -> GridDirectTrfOutcome:
    """Own canonical validation, selection, promotion, acceptance, and authority."""
    records = tuple(attempts)
    gate = _grid_finalization_gate(cancellation)
    if gate is not None:
        return GridDirectTrfOutcome(gate, records)
    eligible = tuple(
        record
        for record in records
        if record.disposition is GridSeedDisposition.ELIGIBLE
    )
    if not eligible:
        return GridDirectTrfOutcome(
            GridDirectTrfTerminal.NO_ELIGIBLE_CANDIDATE,
            records,
        )
    try:
        for record in eligible:
            seed = invocation.seeds[record.seed_ordinal]
            record.validate_for_grid_selection(problem, invocation, seed)
    except Exception as error:  # noqa: BLE001 - foreign evidence fails closed
        return GridDirectTrfOutcome(
            GridDirectTrfTerminal.EXECUTION_FAILURE,
            records,
            failure=TerminalFailure(
                "grid_candidate_lineage_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    gate = _grid_finalization_gate(cancellation)
    if gate is not None:
        return GridDirectTrfOutcome(gate, records)
    selection, selected = _selection(eligible)
    selected_seed = invocation.seeds[selected.seed_ordinal]
    selected_candidate = cast("GridCandidate", selected.candidate)
    promoted = _materialize_derived_direct_trf_candidate_for_root(
        problem,
        selected_candidate.candidate,
        parameterization,
        engine,
        cancellation=cancellation,
    )
    if isinstance(promoted, CandidateMaterialization):
        terminal = {
            MaterializationTerminal.CANCELLED: GridDirectTrfTerminal.CANCELLED,
            MaterializationTerminal.INTERRUPTED: GridDirectTrfTerminal.INTERRUPTED,
        }.get(
            promoted.terminal,
            GridDirectTrfTerminal.MATERIALIZATION_FAILURE,
        )
        return GridDirectTrfOutcome(
            terminal,
            records,
            selection,
            promoted,
            failure=promoted.failure,
        )
    provenance = GridSelectionProvenance(
        invocation.identity,
        problem.identity,
        selection.identity,
        selected.seed_identity,
        selected.seed_ordinal,
        selected_candidate.identity,
        selected_candidate.candidate.identity,
        selected_candidate.candidate.problem_identity,
        selected_candidate.candidate.invocation_identity,
        selected_candidate.candidate.execution_identity,
        selected_candidate.candidate.materialization.identity,
        promoted.materialization.identity,
        promoted.evaluation_result.identity,
    )
    validate_grid_acceptance_lineage(
        provenance,
        invocation,
        selection,
        selected_seed,
        selected,
        promoted,
        problem,
    )
    root_accepted = _accept_materialized_fit_for_derived_workflow(
        problem=problem,
        invocation_identity=promoted.invocation_identity,
        execution_identity=promoted.execution_identity,
        materialization=promoted.materialization,
        vector=promoted.vector,
        chi_square=promoted.chi_square,
        evaluation_result=promoted.evaluation_result,
        authority_context_identity=provenance.identity,
    )
    accepted = AcceptedGridDirectTrfResult.from_accepted(
        root_accepted,
        provenance,
        invocation,
        selection,
        selected_seed,
        selected,
        promoted,
    )
    validate_grid_commit_lineage(accepted, problem)
    authority = _grant_derived_fit_commit_authority(accepted, problem)
    return GridDirectTrfOutcome(
        GridDirectTrfTerminal.ACCEPTED,
        records,
        selection,
        promoted.materialization,
        accepted,
        authority,
    )


class _GridSeedInterrupted(KeyboardInterrupt):
    def __init__(self, outcome: GridSeedOutcome) -> None:
        self.outcome = outcome
        super().__init__("Native GRID seed was interrupted")


def _seed_execution_exception(seed: GridSeed, error: Exception) -> GridSeedOutcome:
    failure = TerminalFailure(
        "grid_seed_execution_exception",
        f"{type(error).__name__}: {error}",
    )
    return GridSeedOutcome(
        seed.identity,
        seed.ordinal,
        seed.axis_items,
        seed.start,
        GridSeedDisposition.IMPLEMENTATION_FAILURE,
        failure=failure,
    )


def _execute_grid_seed(
    seed: GridSeed,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    token: CancellationToken,
) -> GridSeedOutcome:
    if seed.rejection is not None:
        disposition = (
            GridSeedDisposition.UNUSABLE_VALUE
            if seed.rejection.category == "grid_seed_unusable_value"
            else GridSeedDisposition.OUT_OF_BOUNDS
        )
        return GridSeedOutcome(
            seed.identity,
            seed.ordinal,
            seed.axis_items,
            seed.start,
            disposition,
            failure=seed.rejection,
        )
    if seed.problem is None or seed.invocation is None:
        return _seed_execution_exception(
            seed,
            RuntimeError("A feasible GRID seed lacks its Direct TRF child"),
        )
    try:
        direct = execute_direct_trf_candidate(
            seed.problem,
            seed.invocation,
            parameterization,
            engine,
            cancellation=token,
        )
    except DirectTrfInterrupted as interrupted:
        raise _GridSeedInterrupted(
            _interrupted_attempt(seed, interrupted)
        ) from interrupted
    except KeyboardInterrupt as error:
        raise _GridSeedInterrupted(
            GridSeedOutcome(
                seed.identity,
                seed.ordinal,
                seed.axis_items,
                seed.start,
                GridSeedDisposition.INTERRUPTED,
                failure=TerminalFailure("interrupted", "KeyboardInterrupt"),
            )
        ) from error
    except Exception as error:  # noqa: BLE001 - orchestration fails closed
        return _seed_execution_exception(seed, error)
    if direct.terminal is DirectTrfCandidateTerminal.SUCCESS:
        if direct.candidate is None:
            return _seed_execution_exception(
                seed,
                RuntimeError("Successful Direct TRF seed lacks its candidate"),
            )
        candidate = GridCandidate(seed.identity, seed.ordinal, direct.candidate)
        return GridSeedOutcome(
            seed.identity,
            seed.ordinal,
            seed.axis_items,
            seed.start,
            GridSeedDisposition.ELIGIBLE,
            direct.execution.identity,
            candidate.objective,
            candidate,
            direct,
        )
    if direct.terminal is DirectTrfCandidateTerminal.CANCELLED:
        return _attempt_failure(seed, GridSeedDisposition.CANCELLED, direct)
    return _attempt_failure(seed, _disposition_for_unsuccessful(direct), direct)


def _execute_grid_seeds(
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    token: CancellationToken,
) -> tuple[GridSeedOutcome, ...] | GridDirectTrfOutcome:
    attempts: list[GridSeedOutcome] = []
    for index, seed in enumerate(invocation.seeds):
        remaining = invocation.seeds[index + 1 :]
        if token.is_cancelled:
            return GridDirectTrfOutcome(
                GridDirectTrfTerminal.CANCELLED,
                (*attempts, *(_not_started(item) for item in invocation.seeds[index:])),
            )
        try:
            seed_outcome = _execute_grid_seed(
                seed,
                parameterization,
                engine,
                token,
            )
        except _GridSeedInterrupted as interrupted:
            attempts.append(interrupted.outcome)
            raise GridDirectTrfInterrupted(
                _interrupted_outcome(attempts, remaining)
            ) from interrupted
        attempts.append(seed_outcome)
        if seed_outcome.disposition is GridSeedDisposition.CANCELLED:
            return GridDirectTrfOutcome(
                GridDirectTrfTerminal.CANCELLED,
                (*attempts, *(_not_started(item) for item in remaining)),
            )
        if seed_outcome.disposition in {
            GridSeedDisposition.BACKEND_FAILURE,
            GridSeedDisposition.IMPLEMENTATION_FAILURE,
        }:
            return GridDirectTrfOutcome(
                GridDirectTrfTerminal.EXECUTION_FAILURE,
                (*attempts, *(_not_started(item) for item in remaining)),
                failure=seed_outcome.failure,
            )
    return tuple(attempts)


def execute_grid_direct_trf(
    problem: OptimizationProblem,
    invocation: GridDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> GridDirectTrfOutcome:
    """Run every expected seed failure safely, then accept only the winner."""
    _validate_grid_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    seed_execution = _execute_grid_seeds(
        invocation,
        parameterization,
        engine,
        token,
    )
    if isinstance(seed_execution, GridDirectTrfOutcome):
        return seed_execution
    return _finalize_grid_candidates(
        problem,
        invocation,
        seed_execution,
        parameterization,
        engine,
        token,
    )
