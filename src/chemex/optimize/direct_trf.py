"""Bounded native Direct-TRF qualification path for issue #586.

This is deliberately not a production dispatcher or a generic optimizer layer.
It connects the native parameterization and evaluation contracts to one closed
SciPy ``least_squares(method="trf")`` invocation, fresh result materialization,
and the existing revision-checked Analysis Values commit boundary.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real
from threading import Event
from typing import cast
from uuid import uuid4

import numpy as np
from scipy.optimize import least_squares

from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    IndependentValueFrame,
    ParameterRole,
)
from chemex.parameters.sealed import SealedConfiguration
from chemex.parameters.values import (
    AnalysisValues,
    AnalysisValuesSnapshot,
)
from chemex.typing import Array

_SCHEMA_VERSION = 1
_PROBLEM_VERSION = "native-direct-trf-problem-v1"
_SCALARIZATION_VERSION = "ordered-pairwise-chi-square-v1"
_INVOCATION_VERSION = "scipy-direct-trf-v1"
_CANDIDATE_ORDER_VERSION = "chi-square-vector-ordinal-v1"
_BACKEND_SETTINGS = (
    ("method", "trf"),
    ("jac", "2-point"),
    ("diff_step", None),
    ("tr_solver", "exact"),
    ("tr_options", ()),
    ("loss", "linear"),
    ("f_scale", "0x1.0000000000000p+0"),
    ("verbose", 0),
)


class DirectTrfConstructionError(ValueError):
    """Raised before an execution attempt when native semantics are incompatible."""


class ObjectiveScalarizationError(ValueError):
    """Raised when canonical chi-square cannot be formed as finite binary64."""


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {
            "kind": kind,
            "record": record,
            "schema_version": _SCHEMA_VERSION,
        },
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _finite_binary64(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real binary64 scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return 0.0 if result == 0.0 else result


def _bound_binary64(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real binary64 scalar")
    result = float(value)
    if math.isnan(result):
        raise ValueError(f"{name} cannot be NaN")
    return 0.0 if result == 0.0 else result


def _float_token(value: float) -> str:
    return float(value).hex()


def _vector_tokens(values: Sequence[float]) -> tuple[str, ...]:
    return tuple(_float_token(value) for value in values)


def _canonical_vector(
    values: object,
    *,
    dimension: int,
    name: str,
) -> tuple[float, ...]:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (dimension,):
        raise ValueError(f"{name} must have shape ({dimension},)")
    return tuple(
        _finite_binary64(value, name=f"{name}[{index}]")
        for index, value in enumerate(array)
    )


def _validate_problem_ordering(
    controlled_ids: tuple[str, ...],
    independent_items: tuple[tuple[str, float], ...],
    held_items: tuple[tuple[str, float], ...],
    commit_scope: tuple[str, ...],
) -> None:
    if not controlled_ids:
        raise DirectTrfConstructionError(
            "Direct TRF requires at least one controlled parameter"
        )
    if len(set(controlled_ids)) != len(controlled_ids):
        raise DirectTrfConstructionError("Controlled parameter IDs are not unique")
    if len(set(commit_scope)) != len(commit_scope):
        raise DirectTrfConstructionError("Commit scope IDs are not unique")
    independent_ids = tuple(param_id for param_id, _value in independent_items)
    held_ids = tuple(param_id for param_id, _value in held_items)
    if len(set(independent_ids)) != len(independent_ids):
        raise DirectTrfConstructionError("Independent parameter IDs are not unique")
    canonical_controlled = tuple(
        param_id for param_id in independent_ids if param_id in controlled_ids
    )
    if canonical_controlled != controlled_ids:
        raise DirectTrfConstructionError("Controlled parameter order is not canonical")
    canonical_held = tuple(
        param_id for param_id in independent_ids if param_id not in controlled_ids
    )
    if canonical_held != held_ids:
        raise DirectTrfConstructionError("Held parameter order is not canonical")


def _normalize_problem_bounds(
    dimension: int,
    start: tuple[float, ...],
    lower_bounds: tuple[float, ...],
    upper_bounds: tuple[float, ...],
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    if len(start) != dimension:
        raise DirectTrfConstructionError("Start vector has the wrong dimension")
    if len(lower_bounds) != dimension or len(upper_bounds) != dimension:
        raise DirectTrfConstructionError("Bound vectors have the wrong dimension")
    normalized_start = tuple(
        _finite_binary64(value, name=f"start[{index}]")
        for index, value in enumerate(start)
    )
    lower = tuple(
        _bound_binary64(value, name=f"lower bound[{index}]")
        for index, value in enumerate(lower_bounds)
    )
    upper = tuple(
        _bound_binary64(value, name=f"upper bound[{index}]")
        for index, value in enumerate(upper_bounds)
    )
    for index, (value, minimum, maximum) in enumerate(
        zip(normalized_start, lower, upper, strict=True)
    ):
        if minimum >= maximum:
            raise DirectTrfConstructionError(
                f"Controlled bound interval {index} is empty or fixed"
            )
        if not minimum <= value <= maximum:
            raise DirectTrfConstructionError(
                f"Start coordinate {index} is outside its effective bounds"
            )
    return normalized_start, lower, upper


def canonical_chi_square(residuals: Sequence[float] | Array) -> float:
    """Apply the versioned ordered pairwise binary64 residual scalarization."""
    terms: list[float] = []
    for ordinal, residual in enumerate(residuals):
        try:
            value = _finite_binary64(residual, name=f"residual[{ordinal}]")
        except (TypeError, ValueError) as error:
            raise ObjectiveScalarizationError(str(error)) from error
        square = value * value
        if not math.isfinite(square):
            raise ObjectiveScalarizationError(
                f"residual square is non-finite at ordinal {ordinal}"
            )
        terms.append(0.0 if square == 0.0 else square)
    while len(terms) > 1:
        reduced: list[float] = []
        for index in range(0, len(terms) - 1, 2):
            partial = terms[index] + terms[index + 1]
            if not math.isfinite(partial):
                raise ObjectiveScalarizationError(
                    f"chi-square partial sum is non-finite at level index {index // 2}"
                )
            reduced.append(0.0 if partial == 0.0 else partial)
        if len(terms) % 2:
            reduced.append(terms[-1])
        terms = reduced
    return terms[0] if terms else 0.0


@dataclass(frozen=True, slots=True)
class CandidateSummary:
    """A successfully evaluated candidate under the canonical total order."""

    vector: tuple[float, ...]
    chi_square: float
    ordinal: int

    def __post_init__(self) -> None:
        vector = tuple(
            _finite_binary64(value, name=f"candidate[{index}]")
            for index, value in enumerate(self.vector)
        )
        chi_square = _finite_binary64(self.chi_square, name="candidate chi-square")
        if not vector:
            raise ValueError("Candidate vector cannot be empty")
        if chi_square < 0.0:
            raise ValueError("Candidate chi-square cannot be negative")
        if isinstance(self.ordinal, bool) or self.ordinal < 0:
            raise ValueError("Candidate ordinal must be a non-negative integer")
        object.__setattr__(self, "vector", vector)
        object.__setattr__(self, "chi_square", chi_square)

    def ordering_key(self) -> tuple[float, tuple[float, ...], int]:
        """Return the declared χ², external-vector, ordinal total order."""
        return self.chi_square, self.vector, self.ordinal

    @property
    def identity(self) -> str:
        return _identity(
            "direct-trf-candidate",
            (
                _CANDIDATE_ORDER_VERSION,
                _vector_tokens(self.vector),
                _float_token(self.chi_square),
                self.ordinal,
            ),
        )


@dataclass(frozen=True, slots=True)
class OptimizationProblem:
    """One immutable canonical external-coordinate fit and commit context."""

    evaluation_plan_identity: str
    parameterization_identity: str
    evaluator_parameterization_identity: str
    constraint_program_identity: str
    configuration_identity: str
    source_snapshot: AnalysisValuesSnapshot = field(repr=False, compare=False)
    independent_items: tuple[tuple[str, float], ...]
    controlled_ids: tuple[str, ...]
    held_items: tuple[tuple[str, float], ...]
    start: tuple[float, ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    commit_scope: tuple[str, ...]
    scalarization_version: str = _SCALARIZATION_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_problem_ordering(
            self.controlled_ids,
            self.independent_items,
            self.held_items,
            self.commit_scope,
        )
        normalized_start, lower, upper = _normalize_problem_bounds(
            len(self.controlled_ids),
            self.start,
            self.lower_bounds,
            self.upper_bounds,
        )
        object.__setattr__(self, "start", normalized_start)
        object.__setattr__(self, "lower_bounds", lower)
        object.__setattr__(self, "upper_bounds", upper)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-direct-trf-problem",
                (
                    _PROBLEM_VERSION,
                    self.evaluation_plan_identity,
                    self.parameterization_identity,
                    self.evaluator_parameterization_identity,
                    self.constraint_program_identity,
                    self.configuration_identity,
                    self.source_snapshot.occurrence_identity,
                    self.source_snapshot.revision,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in self.independent_items
                    ),
                    self.controlled_ids,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in self.held_items
                    ),
                    _vector_tokens(normalized_start),
                    _vector_tokens(lower),
                    _vector_tokens(upper),
                    self.commit_scope,
                    self.scalarization_version,
                ),
            ),
        )

    @classmethod
    def from_native(
        cls,
        plan: EvaluationPlan,
        parameterization: ActiveParameterization,
        configuration: SealedConfiguration,
        snapshot: AnalysisValuesSnapshot,
    ) -> OptimizationProblem:
        """Construct exact vector, held values, bounds, and commit scope."""
        if plan.parameterization_identity != parameterization.evaluator_identity:
            raise DirectTrfConstructionError(
                "Evaluation plan belongs to another parameterization"
            )
        if plan.constraint_program_identity != parameterization.program.fingerprint:
            raise DirectTrfConstructionError(
                "Evaluation plan belongs to another constraint program"
            )
        if plan.resolved_ids != parameterization.scope_ids:
            raise DirectTrfConstructionError(
                "Evaluation plan does not cover the complete active scope"
            )
        if plan.retained_observation_count < 1:
            raise DirectTrfConstructionError("Residual objective cannot be empty")
        if configuration.identity != snapshot.configuration_identity:
            raise DirectTrfConstructionError(
                "Configuration and starting snapshot identities differ"
            )
        frame = parameterization.frame_from_snapshot(snapshot)
        independent_items = frame.ordered_items()
        controlled_ids = tuple(
            param_id
            for param_id in parameterization.independent_ids
            if parameterization.role(param_id) is ParameterRole.FIT
        )
        controlled = set(controlled_ids)
        held_items = tuple(
            (param_id, value)
            for param_id, value in independent_items
            if param_id not in controlled
        )
        independent_values = dict(independent_items)
        start = tuple(independent_values[param_id] for param_id in controlled_ids)
        lower = tuple(
            configuration[param_id].lower_bound for param_id in controlled_ids
        )
        upper = tuple(
            configuration[param_id].upper_bound for param_id in controlled_ids
        )
        for param_id, value in held_items:
            item = configuration[param_id]
            if not item.lower_bound <= value <= item.upper_bound:
                raise DirectTrfConstructionError(
                    f"Held parameter {param_id!r} is outside its effective bounds"
                )
        return cls(
            plan.identity,
            parameterization.identity,
            parameterization.evaluator_identity,
            parameterization.program.fingerprint,
            configuration.identity,
            snapshot,
            independent_items,
            controlled_ids,
            held_items,
            start,
            lower,
            upper,
            parameterization.scope_ids,
        )

    def validate_parameterization(
        self, parameterization: ActiveParameterization
    ) -> None:
        """Reject a parameterization outside this problem's lifecycle context."""
        if (
            parameterization.identity != self.parameterization_identity
            or parameterization.evaluator_identity
            != self.evaluator_parameterization_identity
            or parameterization.program.fingerprint != self.constraint_program_identity
        ):
            raise DirectTrfConstructionError(
                "Optimization problem belongs to another parameterization"
            )

    def lifecycle_frame(
        self,
        vector: Sequence[float] | Array,
        parameterization: ActiveParameterization,
    ) -> IndependentValueFrame:
        """Combine one canonical controlled vector with captured held values."""
        self.validate_parameterization(parameterization)
        candidate = _canonical_vector(
            vector,
            dimension=len(self.controlled_ids),
            name="external candidate",
        )
        updates = dict(zip(self.controlled_ids, candidate, strict=True))
        return IndependentValueFrame(
            parameterization.identity,
            parameterization.program.fingerprint,
            self.source_snapshot.occurrence_identity,
            self.source_snapshot.revision,
            self.independent_items,
        ).with_updates(updates)


@dataclass(frozen=True, slots=True)
class DirectTrfInvocation:
    """Closed immutable settings for one atomic SciPy Direct-TRF attempt."""

    problem_identity: str
    objective_request_budget: int
    x_scale: tuple[float, ...]
    ftol: float | None = 1.0e-8
    xtol: float | None = 1.0e-8
    gtol: float | None = 1.0e-8
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.objective_request_budget, bool)
            or not isinstance(self.objective_request_budget, int)
            or self.objective_request_budget < 1
        ):
            raise DirectTrfConstructionError(
                "Objective-request budget must be a positive integer"
            )
        scales = tuple(
            _finite_binary64(value, name=f"x_scale[{index}]")
            for index, value in enumerate(self.x_scale)
        )
        if not scales or any(value <= 0.0 for value in scales):
            raise DirectTrfConstructionError("x_scale entries must be positive")
        epsilon = float(np.finfo(np.float64).eps)
        tolerances: list[float | None] = []
        for name, tolerance in (
            ("ftol", self.ftol),
            ("xtol", self.xtol),
            ("gtol", self.gtol),
        ):
            if tolerance is None:
                tolerances.append(None)
                continue
            value = _finite_binary64(tolerance, name=name)
            if value <= epsilon:
                raise DirectTrfConstructionError(
                    f"{name} must be greater than binary64 machine epsilon"
                )
            tolerances.append(value)
        if all(value is None for value in tolerances):
            raise DirectTrfConstructionError(
                "At least one Direct-TRF convergence tolerance must be enabled"
            )
        object.__setattr__(self, "x_scale", scales)
        object.__setattr__(self, "ftol", tolerances[0])
        object.__setattr__(self, "xtol", tolerances[1])
        object.__setattr__(self, "gtol", tolerances[2])
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-direct-trf-invocation",
                (
                    _INVOCATION_VERSION,
                    self.problem_identity,
                    self.objective_request_budget,
                    _vector_tokens(scales),
                    tuple(
                        None if value is None else _float_token(value)
                        for value in tolerances
                    ),
                    _BACKEND_SETTINGS,
                ),
            ),
        )

    @classmethod
    def for_problem(
        cls,
        problem: OptimizationProblem,
        *,
        objective_request_budget: int,
        x_scale: float | Sequence[float] | None = None,
        ftol: float | None = 1.0e-8,
        xtol: float | None = 1.0e-8,
        gtol: float | None = 1.0e-8,
    ) -> DirectTrfInvocation:
        """Resolve scalar/default scaling into one value per external coordinate."""
        dimension = len(problem.controlled_ids)
        if x_scale is None:
            scales = (1.0,) * dimension
        elif isinstance(x_scale, Real) and not isinstance(x_scale, bool):
            scales = (float(x_scale),) * dimension
        else:
            scales = tuple(cast("Sequence[float]", x_scale))
        if len(scales) != dimension:
            raise DirectTrfConstructionError("x_scale has the wrong dimension")
        return cls(
            problem.identity,
            objective_request_budget,
            scales,
            ftol,
            xtol,
            gtol,
        )


class DirectTrfTerminal(StrEnum):
    """Closed terminal outcomes for the atomic solver attempt."""

    CONVERGED = "converged"
    NON_CONVERGED = "non_converged"
    PREFLIGHT_INVALID = "preflight_invalid"
    BUDGET_EXHAUSTED = "budget_exhausted"
    INVALID_TRIAL = "invalid_trial"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    BACKEND_FAILURE = "backend_failure"
    IMPLEMENTATION_FAILURE = "implementation_failure"


@dataclass(frozen=True, slots=True)
class AttemptCounters:
    """ChemEx-owned objective counters, separate from backend evidence."""

    solver_requests_received: int
    objective_requests_accepted: int
    objective_evaluations_completed: int


@dataclass(frozen=True, slots=True)
class TerminalFailure:
    """Stable sanitized reason for one unsuccessful terminal outcome."""

    category: str
    message: str = ""
    evaluation_failure: EvaluationFailure | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    @property
    def identity(self) -> str:
        return _identity(
            "native-direct-trf-failure",
            (
                self.category,
                self.message,
                None
                if self.evaluation_failure is None
                else self.evaluation_failure.identity,
            ),
        )


@dataclass(frozen=True, slots=True)
class BackendEvidence:
    """Detached closed evidence normalized from SciPy's mutable result."""

    status: int
    success: bool
    message: str
    nfev: int
    njev: int | None
    cost: float
    optimality: float
    active_mask: tuple[int, ...]
    final_vector: tuple[float, ...]
    final_residuals: tuple[float, ...]

    @property
    def identity(self) -> str:
        return _identity(
            "native-direct-trf-backend-evidence",
            (
                self.status,
                self.success,
                self.message,
                self.nfev,
                self.njev,
                _float_token(self.cost),
                _float_token(self.optimality),
                self.active_mask,
                _vector_tokens(self.final_vector),
                _vector_tokens(self.final_residuals),
            ),
        )


@dataclass(frozen=True, slots=True)
class DirectTrfExecution:
    """Immutable record of one begun atomic attempt."""

    occurrence_identity: str = field(compare=False)
    problem_identity: str
    invocation_identity: str
    terminal: DirectTrfTerminal
    counters: AttemptCounters
    preflight_evaluation_identity: str | None
    best_candidate: CandidateSummary | None
    final_candidate: CandidateSummary | None
    backend: BackendEvidence | None
    failure: TerminalFailure | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.terminal is DirectTrfTerminal.CONVERGED and (
            self.final_candidate is None
            or self.backend is None
            or self.failure is not None
        ):
            raise ValueError("Converged Direct TRF lacks complete final evidence")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-direct-trf-execution",
                (
                    self.problem_identity,
                    self.invocation_identity,
                    self.terminal.value,
                    (
                        self.counters.solver_requests_received,
                        self.counters.objective_requests_accepted,
                        self.counters.objective_evaluations_completed,
                    ),
                    self.preflight_evaluation_identity,
                    None
                    if self.best_candidate is None
                    else self.best_candidate.identity,
                    None
                    if self.final_candidate is None
                    else self.final_candidate.identity,
                    None if self.backend is None else self.backend.identity,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )


class DirectTrfInterrupted(KeyboardInterrupt):
    """Propagated interruption carrying the frozen atomic-attempt record."""

    def __init__(
        self,
        execution: DirectTrfExecution,
        materialization: AcceptanceMaterialization | None = None,
    ) -> None:
        self.execution = execution
        self.materialization = materialization
        super().__init__("Native Direct TRF was interrupted")


class MaterializationTerminal(StrEnum):
    """Closed outcomes of fresh authoritative candidate materialization."""

    SUCCESS = "success"
    FAILURE = "failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


@dataclass(frozen=True, slots=True)
class AcceptanceMaterialization:
    """One fresh-workspace aggregate materialization record."""

    occurrence_identity: str = field(compare=False)
    problem_identity: str
    invocation_identity: str
    execution_identity: str
    candidate: CandidateSummary
    terminal: MaterializationTerminal
    evaluation_count: int
    evaluator_compatibility_identity: str
    evaluation_identity: str | None
    cache_hits: int
    cache_misses: int
    failure: TerminalFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.terminal is MaterializationTerminal.SUCCESS:
            if (
                self.evaluation_count != 1
                or self.evaluation_identity is None
                or self.failure is not None
            ):
                raise ValueError("Successful materialization lacks complete evidence")
        elif self.evaluation_identity is not None:
            raise ValueError("Failed materialization cannot expose a usable evaluation")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-direct-trf-materialization",
                (
                    self.problem_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.candidate.identity,
                    self.terminal.value,
                    self.evaluation_count,
                    self.evaluator_compatibility_identity,
                    self.evaluation_identity,
                    self.cache_hits,
                    self.cache_misses,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class AcceptedFitResult:
    """Fresh-materialized immutable fit result eligible for one atomic commit."""

    occurrence_identity: str = field(compare=False)
    problem_identity: str
    invocation_identity: str
    execution_identity: str
    materialization_identity: str
    parameterization_identity: str
    evaluator_parameterization_identity: str
    source_occurrence_identity: str
    source_revision: int
    controlled_ids: tuple[str, ...]
    vector: tuple[float, ...]
    chi_square: float
    evaluation_result: EvaluationResult = field(repr=False, compare=False)
    commit_scope: tuple[str, ...]
    commit_items: tuple[tuple[str, float], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.evaluation_result.residuals.size < 1:
            raise ValueError("Accepted result cannot have an empty residual objective")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-accepted-fit-result",
                (
                    self.problem_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.materialization_identity,
                    self.parameterization_identity,
                    self.evaluator_parameterization_identity,
                    self.source_occurrence_identity,
                    self.source_revision,
                    self.controlled_ids,
                    _vector_tokens(self.vector),
                    _float_token(self.chi_square),
                    self.evaluation_result.identity,
                    self.commit_scope,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in self.commit_items
                    ),
                ),
            ),
        )


class DirectTrfOutcomeTerminal(StrEnum):
    """Top-level outcome after optional materialization, before commit."""

    ACCEPTED = "accepted"
    SOLVER_UNSUCCESSFUL = "solver_unsuccessful"
    MATERIALIZATION_FAILURE = "materialization_failure"
    CANCELLED = "cancelled"


@dataclass(frozen=True, slots=True)
class DirectTrfOutcome:
    """Typed terminal outcome for the bounded qualification slice."""

    terminal: DirectTrfOutcomeTerminal
    execution: DirectTrfExecution
    materialization: AcceptanceMaterialization | None = None
    accepted_result: AcceptedFitResult | None = None

    def __post_init__(self) -> None:
        accepted = self.terminal is DirectTrfOutcomeTerminal.ACCEPTED
        if accepted != (self.accepted_result is not None):
            raise ValueError("Only ACCEPTED may expose an Accepted Fit Result")
        if accepted and (
            self.materialization is None
            or self.materialization.terminal is not MaterializationTerminal.SUCCESS
        ):
            raise ValueError("Accepted outcome lacks successful materialization")
        if self.terminal is DirectTrfOutcomeTerminal.MATERIALIZATION_FAILURE and (
            self.materialization is None
            or self.materialization.terminal is not MaterializationTerminal.FAILURE
        ):
            raise ValueError("Materialization failure lacks its failure record")
        if (
            self.terminal is DirectTrfOutcomeTerminal.CANCELLED
            and self.materialization is not None
            and self.materialization.terminal is not MaterializationTerminal.CANCELLED
        ):
            raise ValueError("Cancelled outcome has incompatible materialization")


@dataclass(frozen=True, slots=True)
class CommitReceipt:
    """Immutable evidence of one successful revision-checked atomic commit."""

    accepted_occurrence_identity: str
    accepted_result_identity: str
    problem_identity: str
    old_revision: int
    new_revision: int
    scope: tuple[str, ...]
    committed_value_identity: str
    model_identity: str
    configuration_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-fit-commit-receipt",
                (
                    self.accepted_occurrence_identity,
                    self.accepted_result_identity,
                    self.problem_identity,
                    self.old_revision,
                    self.new_revision,
                    self.scope,
                    self.committed_value_identity,
                    self.model_identity,
                    self.configuration_identity,
                ),
            ),
        )


class CancellationToken:
    """Occurrence-local cooperative cancellation signal."""

    def __init__(self) -> None:
        self._event = Event()

    @property
    def is_cancelled(self) -> bool:
        return self._event.is_set()

    def cancel(self) -> None:
        self._event.set()


@dataclass(frozen=True, slots=True)
class _CompletedRequest:
    summary: CandidateSummary
    residuals: tuple[float, ...]
    evaluation_identity: str


class _AttemptStop(RuntimeError):
    def __init__(self, terminal: DirectTrfTerminal, failure: TerminalFailure) -> None:
        self.terminal = terminal
        self.failure = failure
        super().__init__(failure.message or failure.category)


class _LiveAttempt:
    def __init__(
        self,
        problem: OptimizationProblem,
        invocation: DirectTrfInvocation,
        parameterization: ActiveParameterization,
        evaluator: BoundEvaluator,
        cancellation: CancellationToken,
    ) -> None:
        self.problem = problem
        self.invocation = invocation
        self.parameterization = parameterization
        self.evaluator = evaluator
        self.cancellation = cancellation
        self.received = 0
        self.accepted = 0
        self.completed = 0
        self.requests: list[_CompletedRequest] = []
        self.best: CandidateSummary | None = None

    @property
    def counters(self) -> AttemptCounters:
        return AttemptCounters(self.received, self.accepted, self.completed)

    def objective(self, solver_vector: Array) -> Array:
        self.received += 1
        if self.cancellation.is_cancelled:
            raise _AttemptStop(
                DirectTrfTerminal.CANCELLED,
                TerminalFailure("cancelled", "Cancellation observed before request"),
            )
        if self.accepted >= self.invocation.objective_request_budget:
            raise _AttemptStop(
                DirectTrfTerminal.BUDGET_EXHAUSTED,
                TerminalFailure(
                    "objective_budget_exhausted",
                    "Objective-request budget exhausted",
                ),
            )
        self.accepted += 1
        try:
            vector = _canonical_vector(
                solver_vector,
                dimension=len(self.problem.controlled_ids),
                name="solver vector",
            )
            if any(
                not lower <= value <= upper
                for value, lower, upper in zip(
                    vector,
                    self.problem.lower_bounds,
                    self.problem.upper_bounds,
                    strict=True,
                )
            ):
                raise ValueError("SciPy supplied an out-of-bounds external candidate")
            lifecycle = self.problem.lifecycle_frame(vector, self.parameterization)
            frame = EvaluationFrame.from_lifecycle_frame(
                self.parameterization,
                lifecycle,
            )
        except Exception as error:
            raise _AttemptStop(
                DirectTrfTerminal.IMPLEMENTATION_FAILURE,
                TerminalFailure(
                    "candidate_contract_failure",
                    f"{type(error).__name__}: {error}",
                ),
            ) from error
        outcome = self.evaluator.evaluate(frame)
        self.completed += 1
        if isinstance(outcome, EvaluationFailure):
            terminal = (
                DirectTrfTerminal.INVALID_TRIAL
                if outcome.validity == "INVALID_TRIAL"
                else DirectTrfTerminal.IMPLEMENTATION_FAILURE
            )
            raise _AttemptStop(
                terminal,
                TerminalFailure(
                    outcome.category,
                    outcome.message,
                    outcome,
                ),
            )
        try:
            residuals = tuple(
                _finite_binary64(value, name=f"residual[{index}]")
                for index, value in enumerate(outcome.residuals)
            )
            chi_square = canonical_chi_square(residuals)
        except (TypeError, ValueError, ObjectiveScalarizationError) as error:
            raise _AttemptStop(
                DirectTrfTerminal.INVALID_TRIAL,
                TerminalFailure(
                    "objective_scalarization_failure",
                    f"{type(error).__name__}: {error}",
                ),
            ) from error
        summary = CandidateSummary(vector, chi_square, self.received)
        request = _CompletedRequest(summary, residuals, outcome.identity)
        self.requests.append(request)
        if self.best is None or summary.ordering_key() < self.best.ordering_key():
            self.best = summary
        if self.cancellation.is_cancelled:
            raise _AttemptStop(
                DirectTrfTerminal.CANCELLED,
                TerminalFailure("cancelled", "Cancellation observed after evaluation"),
            )
        return np.asarray(residuals, dtype=np.float64)


def _preflight(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    evaluator: BoundEvaluator,
) -> EvaluationResult | TerminalFailure:
    try:
        lifecycle = problem.lifecycle_frame(problem.start, parameterization)
        frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    except Exception as error:  # noqa: BLE001 - lifecycle preflight fails closed
        return TerminalFailure(
            "preflight_frame_failure",
            f"{type(error).__name__}: {error}",
        )
    outcome = evaluator.evaluate(frame)
    if isinstance(outcome, EvaluationFailure):
        return TerminalFailure(outcome.category, outcome.message, outcome)
    try:
        canonical_chi_square(outcome.residuals)
    except (TypeError, ValueError, ObjectiveScalarizationError) as error:
        return TerminalFailure(
            "preflight_scalarization_failure",
            f"{type(error).__name__}: {error}",
        )
    return outcome


def _execution(
    occurrence_identity: str,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    live: _LiveAttempt,
    terminal: DirectTrfTerminal,
    preflight_identity: str | None,
    *,
    final_candidate: CandidateSummary | None = None,
    backend: BackendEvidence | None = None,
    failure: TerminalFailure | None = None,
) -> DirectTrfExecution:
    return DirectTrfExecution(
        occurrence_identity,
        problem.identity,
        invocation.identity,
        terminal,
        live.counters,
        preflight_identity,
        live.best,
        final_candidate,
        backend,
        failure,
    )


def _backend_int(value: object, *, name: str, optional: bool = False) -> int | None:
    if value is None and optional:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"SciPy {name} must be an integer")
    result = int(value)
    if result < 0 and name in {"nfev", "njev"}:
        raise ValueError(f"SciPy {name} cannot be negative")
    return result


def _normalize_backend(
    result: object,
    *,
    dimension: int,
    residual_count: int,
) -> BackendEvidence:
    status_value = _backend_int(getattr(result, "status", None), name="status")
    if status_value is None:
        raise TypeError("SciPy status evidence is missing")
    success = getattr(result, "success", None)
    message = getattr(result, "message", None)
    if not isinstance(success, (bool, np.bool_)) or not isinstance(message, str):
        raise TypeError("SciPy success/message evidence is malformed")
    nfev_value = _backend_int(getattr(result, "nfev", None), name="nfev")
    if nfev_value is None:
        raise TypeError("SciPy nfev evidence is missing")
    njev_value = _backend_int(
        getattr(result, "njev", None),
        name="njev",
        optional=True,
    )
    cost = _finite_binary64(getattr(result, "cost", None), name="SciPy cost")
    optimality = _finite_binary64(
        getattr(result, "optimality", None), name="SciPy optimality"
    )
    if cost < 0.0 or optimality < 0.0:
        raise ValueError("SciPy cost and optimality cannot be negative")
    final_vector = _canonical_vector(
        getattr(result, "x", None),
        dimension=dimension,
        name="SciPy final vector",
    )
    final_residuals = _canonical_vector(
        getattr(result, "fun", None),
        dimension=residual_count,
        name="SciPy final residual",
    )
    active_array = np.asarray(getattr(result, "active_mask", None))
    if active_array.shape != (dimension,) or any(
        not isinstance(value, (int, np.integer))
        or isinstance(value, (bool, np.bool_))
        or int(value) not in {-1, 0, 1}
        for value in active_array
    ):
        raise ValueError("SciPy active_mask evidence is malformed")
    return BackendEvidence(
        status_value,
        bool(success),
        message,
        nfev_value,
        njev_value,
        cost,
        optimality,
        tuple(int(value) for value in active_array),
        final_vector,
        final_residuals,
    )


def _matching_request(
    requests: Sequence[_CompletedRequest],
    backend: BackendEvidence,
) -> _CompletedRequest | None:
    for request in reversed(requests):
        if (
            request.summary.vector == backend.final_vector
            and request.residuals == backend.final_residuals
        ):
            return request
    return None


def _failed_outcome(execution: DirectTrfExecution) -> DirectTrfOutcome:
    terminal = (
        DirectTrfOutcomeTerminal.CANCELLED
        if execution.terminal is DirectTrfTerminal.CANCELLED
        else DirectTrfOutcomeTerminal.SOLVER_UNSUCCESSFUL
    )
    return DirectTrfOutcome(terminal, execution)


def _materialization_failure(
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    candidate: CandidateSummary,
    compatibility_identity: str,
    evaluation_count: int,
    cache_hits: int,
    cache_misses: int,
    failure: TerminalFailure,
) -> DirectTrfOutcome:
    materialization = AcceptanceMaterialization(
        uuid4().hex,
        problem.identity,
        invocation.identity,
        execution.identity,
        candidate,
        MaterializationTerminal.FAILURE,
        evaluation_count,
        compatibility_identity,
        None,
        cache_hits,
        cache_misses,
        failure,
    )
    return DirectTrfOutcome(
        DirectTrfOutcomeTerminal.MATERIALIZATION_FAILURE,
        execution,
        materialization,
    )


def _cancelled_materialization(
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    candidate: CandidateSummary,
    compatibility_identity: str,
    evaluation_count: int,
    cache_hits: int,
    cache_misses: int,
) -> DirectTrfOutcome:
    materialization = AcceptanceMaterialization(
        uuid4().hex,
        problem.identity,
        invocation.identity,
        execution.identity,
        candidate,
        MaterializationTerminal.CANCELLED,
        evaluation_count,
        compatibility_identity,
        None,
        cache_hits,
        cache_misses,
        TerminalFailure("cancelled", "Cancellation observed during materialization"),
    )
    return DirectTrfOutcome(
        DirectTrfOutcomeTerminal.CANCELLED,
        execution,
        materialization,
    )


def _interrupted_materialization(
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    candidate: CandidateSummary,
    compatibility_identity: str,
    evaluation_count: int,
    cache_hits: int,
    cache_misses: int,
) -> AcceptanceMaterialization:
    return AcceptanceMaterialization(
        uuid4().hex,
        problem.identity,
        invocation.identity,
        execution.identity,
        candidate,
        MaterializationTerminal.INTERRUPTED,
        evaluation_count,
        compatibility_identity,
        None,
        cache_hits,
        cache_misses,
        TerminalFailure("interrupted", "KeyboardInterrupt during materialization"),
    )


@dataclass(frozen=True, slots=True)
class _ConvergedBackend:
    backend: BackendEvidence
    request: _CompletedRequest


def _validate_execution_context(
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> None:
    if invocation.problem_identity != problem.identity:
        raise DirectTrfConstructionError("Invocation belongs to another problem")
    problem.validate_parameterization(parameterization)
    if engine.plan.identity != problem.evaluation_plan_identity:
        raise DirectTrfConstructionError("Evaluator belongs to another plan")
    if len(invocation.x_scale) != len(problem.controlled_ids):
        raise DirectTrfConstructionError("Invocation has the wrong vector dimension")


def _invoke_least_squares(
    live: _LiveAttempt,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
) -> object:
    return least_squares(
        live.objective,
        np.asarray(problem.start, dtype=np.float64),
        bounds=(
            np.asarray(problem.lower_bounds, dtype=np.float64),
            np.asarray(problem.upper_bounds, dtype=np.float64),
        ),
        method="trf",
        jac="2-point",
        diff_step=None,
        tr_solver="exact",
        tr_options={},
        loss="linear",
        f_scale=1.0,
        x_scale=np.asarray(invocation.x_scale, dtype=np.float64),
        ftol=invocation.ftol,
        xtol=invocation.xtol,
        gtol=invocation.gtol,
        max_nfev=invocation.objective_request_budget,
        verbose=0,
    )


def _terminal_backend_failure(
    backend: BackendEvidence,
) -> tuple[DirectTrfTerminal, TerminalFailure] | None:
    if backend.status == 0 and not backend.success:
        return (
            DirectTrfTerminal.NON_CONVERGED,
            TerminalFailure("non_converged", backend.message),
        )
    if backend.status < 0 and not backend.success:
        return (
            DirectTrfTerminal.BACKEND_FAILURE,
            TerminalFailure("backend_failure", backend.message),
        )
    if not backend.success or backend.status not in {1, 2, 3, 4}:
        return (
            DirectTrfTerminal.IMPLEMENTATION_FAILURE,
            TerminalFailure(
                "inconsistent_backend_terminal",
                "SciPy success and status evidence disagree",
            ),
        )
    return None


def _process_backend_result(
    backend_result: object,
    occurrence_identity: str,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    live: _LiveAttempt,
    preflight_identity: str,
    residual_count: int,
) -> _ConvergedBackend | DirectTrfExecution:
    try:
        backend = _normalize_backend(
            backend_result,
            dimension=len(problem.controlled_ids),
            residual_count=residual_count,
        )
    except Exception as error:  # noqa: BLE001 - malformed third-party result boundary
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.IMPLEMENTATION_FAILURE,
            preflight_identity,
            failure=TerminalFailure(
                "malformed_backend_result",
                f"{type(error).__name__}: {error}",
            ),
        )
    terminal = _terminal_backend_failure(backend)
    if terminal is not None:
        outcome, failure = terminal
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            outcome,
            preflight_identity,
            backend=backend,
            failure=failure,
        )
    request = _matching_request(live.requests, backend)
    if request is None:
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.IMPLEMENTATION_FAILURE,
            preflight_identity,
            backend=backend,
            failure=TerminalFailure(
                "unmatched_backend_candidate",
                "Converged backend candidate lacks exact completed-request evidence",
            ),
        )
    return _ConvergedBackend(backend, request)


def _run_solver_attempt(
    occurrence_identity: str,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    live: _LiveAttempt,
    preflight_identity: str,
    residual_count: int,
) -> _ConvergedBackend | DirectTrfExecution:
    try:
        backend_result = _invoke_least_squares(live, problem, invocation)
    except _AttemptStop as stop:
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            stop.terminal,
            preflight_identity,
            failure=stop.failure,
        )
    except KeyboardInterrupt as error:
        interrupted = _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.INTERRUPTED,
            preflight_identity,
            failure=TerminalFailure("interrupted", "KeyboardInterrupt"),
        )
        raise DirectTrfInterrupted(interrupted) from error
    except Exception as error:  # noqa: BLE001 - undeclared backend failures fail closed
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.IMPLEMENTATION_FAILURE,
            preflight_identity,
            failure=TerminalFailure(
                "unexpected_backend_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    if live.cancellation.is_cancelled:
        return _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.CANCELLED,
            preflight_identity,
            failure=TerminalFailure(
                "cancelled", "Cancellation observed after backend return"
            ),
        )
    return _process_backend_result(
        backend_result,
        occurrence_identity,
        problem,
        invocation,
        live,
        preflight_identity,
        residual_count,
    )


def _validate_materialized_result(
    materialized: EvaluationResult,
    request: _CompletedRequest,
    problem: OptimizationProblem,
    evaluator: BoundEvaluator,
) -> float:
    residuals = tuple(
        _finite_binary64(value, name=f"materialized residual[{index}]")
        for index, value in enumerate(materialized.residuals)
    )
    chi_square = canonical_chi_square(residuals)
    if residuals != request.residuals or chi_square != request.summary.chi_square:
        raise ValueError(
            "Fresh materialization differs from converged request evidence"
        )
    if (
        materialized.identity != request.evaluation_identity
        or materialized.plan_identity != problem.evaluation_plan_identity
        or materialized.parameterization_identity
        != problem.evaluator_parameterization_identity
        or materialized.evaluator_compatibility_identity
        != evaluator.compatibility_identity
        or tuple(materialized.resolved_values) != problem.commit_scope
    ):
        raise ValueError("Fresh materialization identities or scope differ")
    return chi_square


def _materialize_accepted_result(
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    request: _CompletedRequest,
    cancellation: CancellationToken,
) -> DirectTrfOutcome:
    try:
        fresh_evaluator = engine.new_evaluator()
    except KeyboardInterrupt as error:
        materialization = _interrupted_materialization(
            execution,
            problem,
            invocation,
            request.summary,
            engine.compatibility_identity,
            0,
            0,
            0,
        )
        raise DirectTrfInterrupted(execution, materialization) from error
    except Exception as error:  # noqa: BLE001 - fresh binding must fail closed
        return _materialization_failure(
            execution,
            problem,
            invocation,
            request.summary,
            engine.compatibility_identity,
            0,
            0,
            0,
            TerminalFailure(
                "materialization_binding_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    try:
        lifecycle = problem.lifecycle_frame(request.summary.vector, parameterization)
        frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
        materialized = fresh_evaluator.evaluate(frame)
    except KeyboardInterrupt as error:
        statistics = fresh_evaluator.cache_statistics
        materialization = _interrupted_materialization(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
        )
        raise DirectTrfInterrupted(execution, materialization) from error
    except Exception as error:  # noqa: BLE001 - aggregate materialization fails closed
        statistics = fresh_evaluator.cache_statistics
        return _materialization_failure(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "materialization_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    statistics = fresh_evaluator.cache_statistics
    if cancellation.is_cancelled:
        return _cancelled_materialization(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
        )
    if isinstance(materialized, EvaluationFailure):
        return _materialization_failure(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                materialized.category,
                materialized.message,
                materialized,
            ),
        )
    try:
        chi_square = _validate_materialized_result(
            materialized,
            request,
            problem,
            fresh_evaluator,
        )
    except Exception as error:  # noqa: BLE001 - aggregate validation fails closed
        return _materialization_failure(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "materialization_validation_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    if cancellation.is_cancelled:
        return _cancelled_materialization(
            execution,
            problem,
            invocation,
            request.summary,
            fresh_evaluator.compatibility_identity,
            1,
            statistics.hits,
            statistics.misses,
        )
    materialization = AcceptanceMaterialization(
        uuid4().hex,
        problem.identity,
        invocation.identity,
        execution.identity,
        request.summary,
        MaterializationTerminal.SUCCESS,
        1,
        fresh_evaluator.compatibility_identity,
        materialized.identity,
        statistics.hits,
        statistics.misses,
    )
    commit_items = tuple(
        (param_id, materialized.resolved_values[param_id])
        for param_id in problem.commit_scope
    )
    accepted = AcceptedFitResult(
        uuid4().hex,
        problem.identity,
        invocation.identity,
        execution.identity,
        materialization.identity,
        problem.parameterization_identity,
        problem.evaluator_parameterization_identity,
        problem.source_snapshot.occurrence_identity,
        problem.source_snapshot.revision,
        problem.controlled_ids,
        request.summary.vector,
        chi_square,
        materialized,
        problem.commit_scope,
        commit_items,
    )
    return DirectTrfOutcome(
        DirectTrfOutcomeTerminal.ACCEPTED,
        execution,
        materialization,
        accepted,
    )


def execute_direct_trf(
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> DirectTrfOutcome:
    """Execute one bounded Direct-TRF attempt and fresh acceptance materialization."""
    _validate_execution_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    occurrence_identity = uuid4().hex
    evaluator = engine.new_evaluator()
    live = _LiveAttempt(problem, invocation, parameterization, evaluator, token)
    if token.is_cancelled:
        return _failed_outcome(
            _execution(
                occurrence_identity,
                problem,
                invocation,
                live,
                DirectTrfTerminal.CANCELLED,
                None,
                failure=TerminalFailure(
                    "cancelled", "Cancellation observed before preflight"
                ),
            )
        )
    try:
        preflight = _preflight(problem, parameterization, evaluator)
    except KeyboardInterrupt as error:
        interrupted = _execution(
            occurrence_identity,
            problem,
            invocation,
            live,
            DirectTrfTerminal.INTERRUPTED,
            None,
            failure=TerminalFailure("interrupted", "KeyboardInterrupt"),
        )
        raise DirectTrfInterrupted(interrupted) from error
    if isinstance(preflight, TerminalFailure):
        return _failed_outcome(
            _execution(
                occurrence_identity,
                problem,
                invocation,
                live,
                DirectTrfTerminal.PREFLIGHT_INVALID,
                None,
                failure=preflight,
            )
        )
    preflight_identity = preflight.identity
    if token.is_cancelled:
        return _failed_outcome(
            _execution(
                occurrence_identity,
                problem,
                invocation,
                live,
                DirectTrfTerminal.CANCELLED,
                preflight_identity,
                failure=TerminalFailure(
                    "cancelled", "Cancellation observed after preflight"
                ),
            )
        )
    solver = _run_solver_attempt(
        occurrence_identity,
        problem,
        invocation,
        live,
        preflight_identity,
        engine.plan.retained_observation_count,
    )
    if isinstance(solver, DirectTrfExecution):
        return _failed_outcome(solver)
    backend = solver.backend
    request = solver.request
    execution = _execution(
        occurrence_identity,
        problem,
        invocation,
        live,
        DirectTrfTerminal.CONVERGED,
        preflight_identity,
        final_candidate=request.summary,
        backend=backend,
    )

    if token.is_cancelled:
        return DirectTrfOutcome(DirectTrfOutcomeTerminal.CANCELLED, execution)
    return _materialize_accepted_result(
        execution,
        problem,
        invocation,
        parameterization,
        engine,
        request,
        token,
    )


def commit_accepted_fit(
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Commit exactly one validated accepted result at its captured revision."""
    if (
        accepted.problem_identity != problem.identity
        or accepted.parameterization_identity != parameterization.identity
        or accepted.evaluator_parameterization_identity
        != parameterization.evaluator_identity
        or accepted.source_occurrence_identity
        != problem.source_snapshot.occurrence_identity
        or accepted.source_revision != problem.source_snapshot.revision
        or accepted.controlled_ids != problem.controlled_ids
        or accepted.commit_scope != problem.commit_scope
    ):
        raise DirectTrfConstructionError(
            "Accepted result is incompatible with its commit context"
        )
    expected_items = tuple(
        (param_id, accepted.evaluation_result.resolved_values[param_id])
        for param_id in problem.commit_scope
    )
    if accepted.commit_items != expected_items:
        raise DirectTrfConstructionError(
            "Accepted result commit snapshot is not its materialized aggregate"
        )
    committed = analysis_values.commit(
        dict(accepted.commit_items),
        expected=problem.source_snapshot,
        scope=problem.commit_scope,
    )
    committed_value_identity = _identity(
        "native-committed-values",
        tuple(
            (param_id, _float_token(committed[param_id]))
            for param_id in problem.commit_scope
        ),
    )
    return CommitReceipt(
        accepted.occurrence_identity,
        accepted.identity,
        problem.identity,
        problem.source_snapshot.revision,
        committed.revision,
        problem.commit_scope,
        committed_value_identity,
        committed.model_identity,
        committed.configuration_identity,
    )
