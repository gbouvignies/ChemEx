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
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real
from threading import Event, RLock
from typing import SupportsIndex, cast
from uuid import uuid4
from weakref import ReferenceType, WeakKeyDictionary, ref

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
    AnalysisValuesCommitError,
    AnalysisValuesSnapshot,
    IncompatibleAnalysisValuesError,
    InvalidAnalysisValuesCommitError,
    StaleAnalysisValuesError,
)
from chemex.runtime.execution import ExecutionSettings, native_thread_environment
from chemex.typing import Array

_SCHEMA_VERSION = 1
_PROBLEM_VERSION = "native-direct-trf-problem-v1"
_SCALARIZATION_VERSION = "ordered-pairwise-chi-square-v1"
_INVOCATION_VERSION = "scipy-direct-trf-v2"
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
class ComponentProblemDerivation:
    """Exact immutable derivation of one non-authoritative child problem."""

    root_problem_identity: str
    root_affine_feasibility_identity: str
    component_identity: str
    projection_policy: str
    projected_plan_identity: str
    controlled_ids: tuple[str, ...]
    root_start_bindings: tuple[tuple[str, float], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-component-problem-derivation",
                (
                    self.root_problem_identity,
                    self.root_affine_feasibility_identity,
                    self.component_identity,
                    self.projection_policy,
                    self.projected_plan_identity,
                    self.controlled_ids,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in self.root_start_bindings
                    ),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class GridSeedProblemDerivation:
    """Exact lineage for one full-coordinate GRID seed child problem."""

    root_problem_identity: str
    root_affine_feasibility_identity: str
    seed_identity: str
    seed_ordinal: int
    axis_items: tuple[tuple[str, float], ...]
    controlled_ids: tuple[str, ...]
    start: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if isinstance(self.seed_ordinal, bool) or self.seed_ordinal < 0:
            raise DirectTrfConstructionError(
                "GRID seed ordinal must be a non-negative integer"
            )
        axis_items = tuple(
            (param_id, _finite_binary64(value, name=f"GRID axis {param_id!r}"))
            for param_id, value in self.axis_items
        )
        axis_ids = tuple(param_id for param_id, _value in axis_items)
        if (
            not axis_items
            or len(set(axis_ids)) != len(axis_ids)
            or axis_ids
            != tuple(sorted(axis_ids, key=lambda item: item.encode("utf-8")))
            or not set(axis_ids).issubset(self.controlled_ids)
        ):
            raise DirectTrfConstructionError(
                "GRID seed axes must be unique canonical controlled IDs"
            )
        start = tuple(
            _finite_binary64(value, name=f"GRID start[{index}]")
            for index, value in enumerate(self.start)
        )
        if len(start) != len(self.controlled_ids):
            raise DirectTrfConstructionError("GRID seed start has the wrong dimension")
        object.__setattr__(self, "axis_items", axis_items)
        object.__setattr__(self, "start", start)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-seed-problem-derivation",
                (
                    self.root_problem_identity,
                    self.root_affine_feasibility_identity,
                    self.seed_identity,
                    self.seed_ordinal,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in axis_items
                    ),
                    self.controlled_ids,
                    _vector_tokens(start),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class DeSearchProblemDerivation:
    """Exact lineage for one selected-coordinate DE search problem."""

    root_problem_identity: str
    root_affine_feasibility_identity: str
    search_specification_identity: str
    selected_ids: tuple[str, ...]
    captured_held_items: tuple[tuple[str, float], ...]
    start: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            not self.root_problem_identity
            or not self.search_specification_identity
            or not self.selected_ids
            or len(set(self.selected_ids)) != len(self.selected_ids)
            or self.selected_ids
            != tuple(sorted(self.selected_ids, key=lambda item: item.encode("utf-8")))
        ):
            raise DirectTrfConstructionError(
                "DE search derivation requires unique canonical selected IDs"
            )
        held_items = tuple(
            (
                param_id,
                _finite_binary64(value, name=f"DE held value {param_id!r}"),
            )
            for param_id, value in self.captured_held_items
        )
        held_ids = tuple(param_id for param_id, _value in held_items)
        if len(set(held_ids)) != len(held_ids) or set(held_ids) & set(
            self.selected_ids
        ):
            raise DirectTrfConstructionError(
                "DE search derivation has inconsistent held coordinates"
            )
        start = tuple(
            _finite_binary64(value, name=f"DE search start[{index}]")
            for index, value in enumerate(self.start)
        )
        if len(start) != len(self.selected_ids):
            raise DirectTrfConstructionError("DE search start has the wrong dimension")
        object.__setattr__(self, "captured_held_items", held_items)
        object.__setattr__(self, "start", start)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-search-problem-derivation",
                (
                    self.root_problem_identity,
                    self.root_affine_feasibility_identity,
                    self.search_specification_identity,
                    self.selected_ids,
                    tuple(
                        (param_id, _float_token(value))
                        for param_id, value in held_items
                    ),
                    _vector_tokens(start),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class DePolishProblemDerivation:
    """Exact lineage for the one full-coordinate TRF polish after DE."""

    root_problem_identity: str
    root_affine_feasibility_identity: str
    workflow_invocation_identity: str
    search_problem_identity: str
    search_execution_identity: str
    search_candidate_identity: str
    controlled_ids: tuple[str, ...]
    start: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            not self.root_problem_identity
            or not self.workflow_invocation_identity
            or not self.search_problem_identity
            or not self.search_execution_identity
            or not self.search_candidate_identity
            or not self.controlled_ids
        ):
            raise DirectTrfConstructionError(
                "DE polish derivation requires complete transition lineage"
            )
        start = tuple(
            _finite_binary64(value, name=f"DE polish start[{index}]")
            for index, value in enumerate(self.start)
        )
        if len(start) != len(self.controlled_ids):
            raise DirectTrfConstructionError("DE polish start has the wrong dimension")
        object.__setattr__(self, "start", start)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-polish-problem-derivation",
                (
                    self.root_problem_identity,
                    self.root_affine_feasibility_identity,
                    self.workflow_invocation_identity,
                    self.search_problem_identity,
                    self.search_execution_identity,
                    self.search_candidate_identity,
                    self.controlled_ids,
                    _vector_tokens(start),
                ),
            ),
        )


type ProblemDerivation = (
    ComponentProblemDerivation
    | GridSeedProblemDerivation
    | DeSearchProblemDerivation
    | DePolishProblemDerivation
)


def _normalize_affine_restriction(
    restriction_id: str,
    coefficients: tuple[float, ...],
    scalar: float,
    *,
    scalar_name: str,
) -> tuple[tuple[float, ...], float]:
    normalized_coefficients = tuple(
        _finite_binary64(value, name=f"affine coefficient[{index}]")
        for index, value in enumerate(coefficients)
    )
    normalized_scalar = _finite_binary64(scalar, name=scalar_name)
    if not restriction_id or not normalized_coefficients:
        raise DirectTrfConstructionError(
            "Affine restrictions require an ID and full-frame coefficients"
        )
    return normalized_coefficients, normalized_scalar


@dataclass(frozen=True, slots=True)
class AffineHalfSpace:
    """One canonical full independent-frame restriction ``a^T x <= b``."""

    restriction_id: str
    coefficients: tuple[float, ...]
    upper_bound: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        coefficients, upper_bound = _normalize_affine_restriction(
            self.restriction_id,
            self.coefficients,
            self.upper_bound,
            scalar_name="affine upper bound",
        )
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "upper_bound", upper_bound)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-affine-half-space",
                (
                    self.restriction_id,
                    _vector_tokens(coefficients),
                    _float_token(upper_bound),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class AffineEquality:
    """One canonical full independent-frame restriction ``a^T x == b``."""

    restriction_id: str
    coefficients: tuple[float, ...]
    value: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        coefficients, value = _normalize_affine_restriction(
            self.restriction_id,
            self.coefficients,
            self.value,
            scalar_name="affine equality value",
        )
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "value", value)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-affine-equality",
                (
                    self.restriction_id,
                    _vector_tokens(coefficients),
                    _float_token(value),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class CoordinateLineFeasibility:
    """Exact feasible displacement interval along one controlled coordinate."""

    param_id: str
    minimum_displacement: float
    maximum_displacement: float
    lower_limiters: tuple[str, ...]
    upper_limiters: tuple[str, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        minimum = _bound_binary64(
            self.minimum_displacement,
            name="minimum feasible displacement",
        )
        maximum = _bound_binary64(
            self.maximum_displacement,
            name="maximum feasible displacement",
        )
        if not minimum <= 0.0 <= maximum:
            raise DirectTrfConstructionError(
                "Coordinate line feasibility must contain the accepted point"
            )
        object.__setattr__(self, "minimum_displacement", minimum)
        object.__setattr__(self, "maximum_displacement", maximum)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-coordinate-line-feasibility",
                (
                    self.param_id,
                    _float_token(minimum),
                    _float_token(maximum),
                    self.lower_limiters,
                    self.upper_limiters,
                ),
            ),
        )


def _validate_affine_feasibility(
    independent_items: tuple[tuple[str, float], ...],
    half_spaces: tuple[AffineHalfSpace, ...],
    equalities: tuple[AffineEquality, ...],
) -> None:
    restrictions = (*half_spaces, *equalities)
    restriction_ids = tuple(item.restriction_id for item in restrictions)
    if len(set(restriction_ids)) != len(restriction_ids) or any(
        len(item.coefficients) != len(independent_items) for item in restrictions
    ):
        raise DirectTrfConstructionError(
            "Affine feasibility must use unique IDs and full independent-frame coefficients"
        )


def _affine_feasibility_identity_record(
    half_spaces: tuple[AffineHalfSpace, ...],
    equalities: tuple[AffineEquality, ...],
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    if not half_spaces and not equalities:
        return ()
    return (
        ("affine-half-spaces", tuple(item.identity for item in half_spaces)),
        ("affine-equalities", tuple(item.identity for item in equalities)),
    )


def _affine_feasibility_identity(
    half_spaces: tuple[AffineHalfSpace, ...],
    equalities: tuple[AffineEquality, ...],
) -> str:
    return _identity(
        "native-affine-feasibility",
        _affine_feasibility_identity_record(half_spaces, equalities),
    )


def _pairwise_feasibility_sum(values: Sequence[float]) -> float:
    terms = list(values)
    while len(terms) > 1:
        reduced = [
            terms[index] + terms[index + 1] for index in range(0, len(terms) - 1, 2)
        ]
        if len(terms) % 2:
            reduced.append(terms[-1])
        terms = reduced
    result = terms[0] if terms else 0.0
    return 0.0 if result == 0.0 else result


def _full_frame_candidate(
    independent_items: tuple[tuple[str, float], ...],
    controlled_ids: tuple[str, ...],
    candidate: tuple[float, ...],
) -> tuple[float, ...]:
    updates = dict(zip(controlled_ids, candidate, strict=True))
    return tuple(updates.get(param_id, value) for param_id, value in independent_items)


def _affine_total(
    coefficients: tuple[float, ...],
    full_values: tuple[float, ...],
) -> float:
    products = tuple(
        coefficient * value
        for coefficient, value in zip(coefficients, full_values, strict=True)
    )
    return _pairwise_feasibility_sum(products)


def _validate_candidate_feasibility(
    independent_items: tuple[tuple[str, float], ...],
    controlled_ids: tuple[str, ...],
    candidate: tuple[float, ...],
    lower_bounds: tuple[float, ...],
    upper_bounds: tuple[float, ...],
    half_spaces: tuple[AffineHalfSpace, ...],
    equalities: tuple[AffineEquality, ...],
) -> tuple[float, ...]:
    if any(
        not lower <= value <= upper
        for value, lower, upper in zip(
            candidate,
            lower_bounds,
            upper_bounds,
            strict=True,
        )
    ):
        raise DirectTrfConstructionError(
            "External candidate is outside exact box bounds"
        )
    full_values = _full_frame_candidate(
        independent_items,
        controlled_ids,
        candidate,
    )
    for restriction in half_spaces:
        total = _affine_total(restriction.coefficients, full_values)
        slack = restriction.upper_bound - total
        if not math.isfinite(total) or not math.isfinite(slack) or slack < 0.0:
            raise DirectTrfConstructionError(
                f"External candidate violates affine half-space "
                f"{restriction.restriction_id!r}"
            )
    for restriction in equalities:
        total = _affine_total(restriction.coefficients, full_values)
        if not math.isfinite(total) or total != restriction.value:
            raise DirectTrfConstructionError(
                f"External candidate violates affine equality "
                f"{restriction.restriction_id!r}"
            )
    return full_values


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
    derivation: ProblemDerivation | None = None
    scalarization_version: str = _SCALARIZATION_VERSION
    affine_half_spaces: tuple[AffineHalfSpace, ...] = ()
    affine_equalities: tuple[AffineEquality, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - complete problem invariant
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
        _validate_affine_feasibility(
            self.independent_items,
            self.affine_half_spaces,
            self.affine_equalities,
        )
        affine_feasibility_identity = _affine_feasibility_identity(
            self.affine_half_spaces,
            self.affine_equalities,
        )
        _validate_candidate_feasibility(
            self.independent_items,
            self.controlled_ids,
            normalized_start,
            lower,
            upper,
            self.affine_half_spaces,
            self.affine_equalities,
        )
        if isinstance(self.derivation, ComponentProblemDerivation):
            if (
                self.derivation.root_affine_feasibility_identity
                != affine_feasibility_identity
            ):
                raise DirectTrfConstructionError(
                    "Component problem affine feasibility differs from its "
                    "derivation record"
                )
            if (
                self.derivation.projected_plan_identity != self.evaluation_plan_identity
                or self.derivation.controlled_ids != self.controlled_ids
                or self.derivation.root_start_bindings != self.held_items
            ):
                raise DirectTrfConstructionError(
                    "Component problem differs from its derivation record"
                )
        elif isinstance(self.derivation, GridSeedProblemDerivation):
            if (
                self.derivation.root_affine_feasibility_identity
                != affine_feasibility_identity
            ):
                raise DirectTrfConstructionError(
                    "GRID seed affine feasibility differs from its derivation record"
                )
            if (
                self.derivation.controlled_ids != self.controlled_ids
                or self.derivation.start != normalized_start
            ):
                raise DirectTrfConstructionError(
                    "GRID seed problem differs from its derivation record"
                )
        elif isinstance(self.derivation, DeSearchProblemDerivation):
            if (
                self.derivation.root_affine_feasibility_identity
                != affine_feasibility_identity
            ):
                raise DirectTrfConstructionError(
                    "DE search affine feasibility differs from its derivation record"
                )
            if (
                self.derivation.selected_ids != self.controlled_ids
                or self.derivation.captured_held_items != self.held_items
                or self.derivation.start != normalized_start
            ):
                raise DirectTrfConstructionError(
                    "DE search problem differs from its derivation record"
                )
        elif isinstance(self.derivation, DePolishProblemDerivation):
            if (
                self.derivation.root_affine_feasibility_identity
                != affine_feasibility_identity
            ):
                raise DirectTrfConstructionError(
                    "DE polish affine feasibility differs from its derivation record"
                )
            if (
                self.derivation.controlled_ids != self.controlled_ids
                or self.derivation.start != normalized_start
            ):
                raise DirectTrfConstructionError(
                    "DE polish problem differs from its derivation record"
                )
        object.__setattr__(self, "start", normalized_start)
        object.__setattr__(self, "lower_bounds", lower)
        object.__setattr__(self, "upper_bounds", upper)
        if self.derivation is None:
            derivation_record: tuple[tuple[str, str], ...] = ()
        elif isinstance(self.derivation, ComponentProblemDerivation):
            derivation_record = (("derived-fit-component", self.derivation.identity),)
        elif isinstance(self.derivation, GridSeedProblemDerivation):
            derivation_record = (("derived-grid-seed", self.derivation.identity),)
        elif isinstance(self.derivation, DeSearchProblemDerivation):
            derivation_record = (("derived-de-search", self.derivation.identity),)
        else:
            derivation_record = (("derived-de-polish", self.derivation.identity),)
        # Preserve the established box-only problem identity exactly.
        feasibility_record = _affine_feasibility_identity_record(
            self.affine_half_spaces,
            self.affine_equalities,
        )
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
                    *derivation_record,
                    self.scalarization_version,
                    *feasibility_record,
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

    def derive_child(
        self,
        *,
        controlled_ids: tuple[str, ...],
        start: tuple[float, ...],
        derivation: ProblemDerivation,
    ) -> OptimizationProblem:
        """Construct one child while retaining every root-owned invariant."""
        controlled = set(controlled_ids)
        expected_controlled = tuple(
            param_id for param_id in self.controlled_ids if param_id in controlled
        )
        if expected_controlled != controlled_ids:
            raise DirectTrfConstructionError(
                "Derived problem controls a non-canonical root coordinate subset"
            )
        held_items = tuple(
            item for item in self.independent_items if item[0] not in controlled
        )
        root_indices = {
            param_id: index for index, param_id in enumerate(self.controlled_ids)
        }
        evaluation_plan_identity = (
            derivation.projected_plan_identity
            if isinstance(derivation, ComponentProblemDerivation)
            else self.evaluation_plan_identity
        )
        child = OptimizationProblem(
            evaluation_plan_identity,
            self.parameterization_identity,
            self.evaluator_parameterization_identity,
            self.constraint_program_identity,
            self.configuration_identity,
            self.source_snapshot,
            self.independent_items,
            controlled_ids,
            held_items,
            start,
            tuple(
                self.lower_bounds[root_indices[param_id]] for param_id in controlled_ids
            ),
            tuple(
                self.upper_bounds[root_indices[param_id]] for param_id in controlled_ids
            ),
            self.commit_scope,
            derivation,
            self.scalarization_version,
            self.affine_half_spaces,
            self.affine_equalities,
        )
        self.validate_derived_problem(child)
        return child

    def validate_derived_problem(self, child: OptimizationProblem) -> None:
        """Fail closed when a child changes root-owned scientific context."""
        derivation = child.derivation
        controlled = set(child.controlled_ids)
        expected_controlled = tuple(
            param_id for param_id in self.controlled_ids if param_id in controlled
        )
        root_indices = {
            param_id: index for index, param_id in enumerate(self.controlled_ids)
        }
        expected_held = tuple(
            item for item in self.independent_items if item[0] not in controlled
        )
        expected_plan_identity = (
            derivation.projected_plan_identity
            if isinstance(derivation, ComponentProblemDerivation)
            else self.evaluation_plan_identity
        )
        if (
            derivation is None
            or derivation.root_problem_identity != self.identity
            or derivation.root_affine_feasibility_identity
            != self.affine_feasibility_identity
            or expected_controlled != child.controlled_ids
            or child.evaluation_plan_identity != expected_plan_identity
            or child.parameterization_identity != self.parameterization_identity
            or child.evaluator_parameterization_identity
            != self.evaluator_parameterization_identity
            or child.constraint_program_identity != self.constraint_program_identity
            or child.configuration_identity != self.configuration_identity
            or child.source_snapshot is not self.source_snapshot
            or child.independent_items != self.independent_items
            or child.held_items != expected_held
            or child.lower_bounds
            != tuple(
                self.lower_bounds[root_indices[param_id]]
                for param_id in child.controlled_ids
            )
            or child.upper_bounds
            != tuple(
                self.upper_bounds[root_indices[param_id]]
                for param_id in child.controlled_ids
            )
            or child.commit_scope != self.commit_scope
            or child.scalarization_version != self.scalarization_version
            or child.affine_half_spaces != self.affine_half_spaces
            or child.affine_equalities != self.affine_equalities
        ):
            raise DirectTrfConstructionError(
                "Derived problem differs from its root-owned context"
            )

    @property
    def acceptance_authority(self) -> bool:
        """Whether this is the complete root rather than a derived component."""
        return self.derivation is None

    @property
    def affine_feasibility_identity(self) -> str:
        """Return the canonical box-independent affine feasibility identity."""
        return _affine_feasibility_identity(
            self.affine_half_spaces,
            self.affine_equalities,
        )

    def coordinate_line_feasibility(  # noqa: C901 - exact interval intersection
        self,
        vector: Sequence[float] | Array,
        column: int,
    ) -> CoordinateLineFeasibility:
        """Intersect exact box and affine feasibility along one coordinate line."""
        candidate = _canonical_vector(
            vector,
            dimension=len(self.controlled_ids),
            name="coordinate-line center",
        )
        full_values = _validate_candidate_feasibility(
            self.independent_items,
            self.controlled_ids,
            candidate,
            self.lower_bounds,
            self.upper_bounds,
            self.affine_half_spaces,
            self.affine_equalities,
        )
        if not 0 <= column < len(self.controlled_ids):
            raise IndexError("Controlled coordinate column is out of range")
        param_id = self.controlled_ids[column]
        independent_index = next(
            index
            for index, (independent_id, _value) in enumerate(self.independent_items)
            if independent_id == param_id
        )
        minimum = self.lower_bounds[column] - candidate[column]
        maximum = self.upper_bounds[column] - candidate[column]
        lower_limiters = (f"box-lower:{param_id}",) if math.isfinite(minimum) else ()
        upper_limiters = (f"box-upper:{param_id}",) if math.isfinite(maximum) else ()
        for restriction in self.affine_half_spaces:
            coefficient = restriction.coefficients[independent_index]
            if coefficient == 0.0:
                continue
            slack = restriction.upper_bound - _affine_total(
                restriction.coefficients,
                full_values,
            )
            boundary = _bound_binary64(
                slack / coefficient,
                name=f"affine line boundary {restriction.restriction_id!r}",
            )
            limiter = f"affine:{restriction.restriction_id}"
            if coefficient > 0.0:
                if boundary < maximum:
                    maximum = boundary
                    upper_limiters = (limiter,)
                elif boundary == maximum:
                    upper_limiters += (limiter,)
            elif boundary > minimum:
                minimum = boundary
                lower_limiters = (limiter,)
            elif boundary == minimum:
                lower_limiters += (limiter,)
        for restriction in self.affine_equalities:
            coefficient = restriction.coefficients[independent_index]
            if coefficient == 0.0:
                continue
            minimum = 0.0
            maximum = 0.0
            limiter = f"equality:{restriction.restriction_id}"
            lower_limiters = (limiter,)
            upper_limiters = (limiter,)
        return CoordinateLineFeasibility(
            param_id,
            minimum,
            maximum,
            lower_limiters,
            upper_limiters,
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
        _validate_candidate_feasibility(
            self.independent_items,
            self.controlled_ids,
            candidate,
            self.lower_bounds,
            self.upper_bounds,
            self.affine_half_spaces,
            self.affine_equalities,
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
    execution_settings: ExecutionSettings = field(default_factory=ExecutionSettings)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.execution_settings.workers != 1:
            raise DirectTrfConstructionError(
                "One Direct TRF invocation has exactly one ChemEx worker"
            )
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
                    (
                        self.execution_settings.workers,
                        self.execution_settings.native_threads,
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
        execution_settings: ExecutionSettings | None = None,
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
            ExecutionSettings() if execution_settings is None else execution_settings,
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
        materialization: CandidateMaterialization | None = None,
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
class CandidateMaterialization:
    """One fresh-workspace candidate materialization record."""

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
class RootMaterializationFailure:
    """Fresh root lifecycle failure before materialization evidence can exist."""

    terminal: MaterializationTerminal
    evaluation_count: int
    evaluator_compatibility_identity: str
    cache_hits: int
    cache_misses: int
    failure: TerminalFailure


@dataclass(frozen=True, slots=True)
class GridSelectionProvenance:
    """Exact GRID workflow selection that authorizes one root materialization."""

    workflow_invocation_identity: str
    root_problem_identity: str
    selection_identity: str
    selected_seed_identity: str
    selected_seed_ordinal: int
    grid_candidate_identity: str
    materialized_candidate_identity: str
    candidate_problem_identity: str
    candidate_invocation_identity: str
    candidate_execution_identity: str
    candidate_materialization_identity: str
    accepted_materialization_identity: str
    accepted_evaluation_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.selected_seed_ordinal, bool)
            or self.selected_seed_ordinal < 0
        ):
            raise ValueError("Selected GRID seed ordinal must be non-negative")
        identities = (
            self.workflow_invocation_identity,
            self.root_problem_identity,
            self.selection_identity,
            self.selected_seed_identity,
            self.grid_candidate_identity,
            self.materialized_candidate_identity,
            self.candidate_problem_identity,
            self.candidate_invocation_identity,
            self.candidate_execution_identity,
            self.candidate_materialization_identity,
            self.accepted_materialization_identity,
            self.accepted_evaluation_identity,
        )
        if any(not identity for identity in identities):
            raise ValueError("GRID selection provenance identities cannot be empty")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grid-selection-provenance",
                (
                    self.workflow_invocation_identity,
                    self.root_problem_identity,
                    self.selection_identity,
                    self.selected_seed_identity,
                    self.selected_seed_ordinal,
                    self.grid_candidate_identity,
                    self.materialized_candidate_identity,
                    self.candidate_problem_identity,
                    self.candidate_invocation_identity,
                    self.candidate_execution_identity,
                    self.candidate_materialization_identity,
                    self.accepted_materialization_identity,
                    self.accepted_evaluation_identity,
                ),
            ),
        )


class _OpaqueOccurrenceWitness:
    """Noncopyable, nonserializable base for process-local occurrence witnesses."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> _OpaqueOccurrenceWitness:
        raise TypeError("Occurrence witnesses are minted only by native execution")

    def __copy__(self) -> _OpaqueOccurrenceWitness:
        raise TypeError("Occurrence witnesses cannot be copied")

    def __deepcopy__(self, _memo: object) -> _OpaqueOccurrenceWitness:
        raise TypeError("Occurrence witnesses cannot be copied")

    def __reduce__(self) -> tuple[object, ...]:
        raise TypeError("Occurrence witnesses cannot be serialized")

    def __reduce_ex__(self, _protocol: SupportsIndex, /) -> tuple[object, ...]:
        raise TypeError("Occurrence witnesses cannot be serialized")


class _AcceptedOccurrenceWitness(_OpaqueOccurrenceWitness):
    """Opaque process-local witness for one accepted-result occurrence."""

    __slots__ = ()


@dataclass(frozen=True, slots=True)
class _OccurrenceWitnessBinding:
    occurrence_identity: str
    artifact_identity: str | None = None
    artifact_reference: ReferenceType[object] | None = None


class _OccurrenceWitnessRegistry[WitnessT]:
    """Bind opaque witnesses to one occurrence and one immutable artifact."""

    def __init__(self) -> None:
        self._bindings: WeakKeyDictionary[WitnessT, _OccurrenceWitnessBinding] = (
            WeakKeyDictionary()
        )
        self._lock = RLock()

    def mint(self, witness_type: type[WitnessT], occurrence_identity: str) -> WitnessT:
        witness = object.__new__(witness_type)
        with self._lock:
            self._bindings[witness] = _OccurrenceWitnessBinding(occurrence_identity)
        return witness

    def bind(
        self,
        witness: WitnessT,
        occurrence_identity: str,
        artifact_identity: str,
        *,
        artifact_object: object | None = None,
    ) -> bool:
        artifact_reference = None if artifact_object is None else ref(artifact_object)
        with self._lock:
            binding = self._bindings.get(witness)
            if binding is None or binding.occurrence_identity != occurrence_identity:
                return False
            if binding.artifact_identity is None:
                self._bindings[witness] = _OccurrenceWitnessBinding(
                    occurrence_identity,
                    artifact_identity,
                    artifact_reference,
                )
                return True
            return (
                binding.occurrence_identity == occurrence_identity
                and binding.artifact_identity == artifact_identity
                and (
                    binding.artifact_reference is None
                    if artifact_object is None
                    else binding.artifact_reference is not None
                    and binding.artifact_reference() is artifact_object
                )
            )

    def is_bound(
        self,
        witness: WitnessT,
        occurrence_identity: str,
        artifact_identity: str,
        *,
        artifact_object: object | None = None,
    ) -> bool:
        with self._lock:
            binding = self._bindings.get(witness)
            return (
                binding is not None
                and binding.occurrence_identity == occurrence_identity
                and binding.artifact_identity == artifact_identity
                and (
                    binding.artifact_reference is None
                    if artifact_object is None
                    else binding.artifact_reference is not None
                    and binding.artifact_reference() is artifact_object
                )
            )


_ACCEPTED_OCCURRENCE_WITNESSES = _OccurrenceWitnessRegistry[
    _AcceptedOccurrenceWitness
]()


def _mint_accepted_occurrence_witness(
    occurrence_identity: str,
) -> _AcceptedOccurrenceWitness:
    return _ACCEPTED_OCCURRENCE_WITNESSES.mint(
        _AcceptedOccurrenceWitness,
        occurrence_identity,
    )


def _bind_accepted_occurrence_witness(
    witness: _AcceptedOccurrenceWitness,
    occurrence_identity: str,
    accepted_result_identity: str,
) -> bool:
    return _ACCEPTED_OCCURRENCE_WITNESSES.bind(
        witness,
        occurrence_identity,
        accepted_result_identity,
    )


@dataclass(frozen=True, slots=True)
class AcceptedFitResult:
    """Fresh-materialized immutable evidence with no live commit authority."""

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
    origin_context_identity: str
    occurrence_witness: _AcceptedOccurrenceWitness | None = field(
        compare=False,
        repr=False,
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.evaluation_result.residuals.size < 1:
            raise ValueError("Accepted result cannot have an empty residual objective")
        if not self.origin_context_identity:
            raise ValueError("Accepted result requires an origin context identity")
        identity = _identity(
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
                self.origin_context_identity,
            ),
        )
        object.__setattr__(self, "identity", identity)
        if self.occurrence_witness is not None:
            _bind_accepted_occurrence_witness(
                self.occurrence_witness,
                self.occurrence_identity,
                identity,
            )

    @classmethod
    def for_qualification(
        cls,
        *,
        occurrence_identity: str,
        problem_identity: str,
        invocation_identity: str,
        execution_identity: str,
        materialization_identity: str,
        parameterization_identity: str,
        evaluator_parameterization_identity: str,
        source_occurrence_identity: str,
        source_revision: int,
        controlled_ids: tuple[str, ...],
        vector: tuple[float, ...],
        chi_square: float,
        evaluation_result: EvaluationResult,
        commit_scope: tuple[str, ...],
        commit_items: tuple[tuple[str, float], ...],
        origin_context_identity: str,
    ) -> AcceptedFitResult:
        """Construct an authoritative anchor only for isolated qualification."""
        return cls(
            occurrence_identity,
            problem_identity,
            invocation_identity,
            execution_identity,
            materialization_identity,
            parameterization_identity,
            evaluator_parameterization_identity,
            source_occurrence_identity,
            source_revision,
            controlled_ids,
            vector,
            chi_square,
            evaluation_result,
            commit_scope,
            commit_items,
            origin_context_identity,
            occurrence_witness=_mint_accepted_occurrence_witness(occurrence_identity),
        )


def accepted_occurrence_is_authoritative(accepted: AcceptedFitResult) -> bool:
    """Report whether this object retains its exact construction occurrence."""
    if accepted.occurrence_witness is None:
        return False
    return _ACCEPTED_OCCURRENCE_WITNESSES.is_bound(
        accepted.occurrence_witness,
        accepted.occurrence_identity,
        accepted.identity,
    )


class LiveFitCommitAuthority:
    """Opaque process-local token whose authority lives only in Direct TRF."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> LiveFitCommitAuthority:
        raise TypeError(
            "Live fit commit authority is minted only by Direct TRF acceptance"
        )

    def __copy__(self) -> LiveFitCommitAuthority:
        raise TypeError("Live fit commit authority cannot be copied")

    def __deepcopy__(self, _memo: object) -> LiveFitCommitAuthority:
        raise TypeError("Live fit commit authority cannot be copied")

    def __reduce__(self) -> tuple[object, ...]:
        raise TypeError("Live fit commit authority cannot be serialized")

    def __reduce_ex__(self, _protocol: SupportsIndex, /) -> tuple[object, ...]:
        raise TypeError("Live fit commit authority cannot be serialized")


@dataclass(frozen=True, slots=True)
class _CommitAuthorityBinding:
    """Private immutable binding for one exact accepted evidence object."""

    accepted_result: AcceptedFitResult = field(repr=False, compare=False)
    accepted_result_identity: str
    accepted_occurrence_identity: str
    problem_identity: str
    snapshot_context_identity: str
    source_occurrence_identity: str
    source_revision: int
    origin_context_identity: str


_LIVE_FIT_COMMIT_AUTHORITIES: WeakKeyDictionary[
    LiveFitCommitAuthority,
    _CommitAuthorityBinding,
] = WeakKeyDictionary()
_LIVE_FIT_COMMIT_AUTHORITIES_LOCK = RLock()


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
    materialization: CandidateMaterialization | None = None
    accepted_result: AcceptedFitResult | None = None
    commit_authority: LiveFitCommitAuthority | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        accepted = self.terminal is DirectTrfOutcomeTerminal.ACCEPTED
        if accepted != (
            self.accepted_result is not None and self.commit_authority is not None
        ):
            raise ValueError(
                "Only ACCEPTED may expose fit evidence and live commit authority"
            )
        if not accepted and (
            self.accepted_result is not None or self.commit_authority is not None
        ):
            raise ValueError("Unaccepted outcome cannot expose fit authority")
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


class DirectTrfCandidateTerminal(StrEnum):
    """Terminal outcomes for non-authoritative component candidate execution."""

    SUCCESS = "success"
    SOLVER_UNSUCCESSFUL = "solver_unsuccessful"
    MATERIALIZATION_FAILURE = "materialization_failure"
    CANCELLED = "cancelled"


@dataclass(frozen=True, slots=True)
class MaterializedDirectTrfCandidate:
    """Fresh component candidate evidence with no acceptance or commit authority."""

    problem_identity: str
    invocation_identity: str
    execution_identity: str
    materialization: CandidateMaterialization = field(repr=False, compare=False)
    vector: tuple[float, ...]
    chi_square: float
    evaluation_result: EvaluationResult = field(repr=False, compare=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-materialized-direct-trf-candidate",
                (
                    self.problem_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.materialization.identity,
                    _vector_tokens(self.vector),
                    _float_token(self.chi_square),
                    self.evaluation_result.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class DirectTrfCandidateOutcome:
    """Non-authoritative result of one component-local Direct TRF run."""

    terminal: DirectTrfCandidateTerminal
    execution: DirectTrfExecution
    materialization: CandidateMaterialization | None = None
    candidate: MaterializedDirectTrfCandidate | None = None

    def __post_init__(self) -> None:
        succeeded = self.terminal is DirectTrfCandidateTerminal.SUCCESS
        if succeeded != (self.candidate is not None):
            raise ValueError("Only successful component execution exposes a candidate")
        if succeeded and (
            self.materialization is None
            or self.materialization.terminal is not MaterializationTerminal.SUCCESS
        ):
            raise ValueError("Successful component candidate lacks materialization")


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


class FitCommitFailureCategory(StrEnum):
    """Closed failure categories for one attempted fit commit."""

    STALE_REVISION = "stale_revision"
    INCOMPATIBLE_STATE = "incompatible_state"
    INVALID_CANDIDATE = "invalid_candidate"


@dataclass(frozen=True, slots=True)
class FitCommitFailure:
    """Sanitized typed reason that an atomic fit commit did not occur."""

    category: FitCommitFailureCategory
    message: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.message:
            raise ValueError("Fit commit failure requires a message")
        object.__setattr__(
            self,
            "identity",
            _identity("native-fit-commit-failure", (self.category.value, self.message)),
        )


class FitCommitTerminal(StrEnum):
    """Terminal state of one atomic fit commit operation."""

    COMMITTED = "committed"
    FAILED = "failed"


class _FitCommitOperationWitness(_OpaqueOccurrenceWitness):
    """Opaque process-local witness for one fit-commit operation occurrence."""

    __slots__ = ()


_FIT_COMMIT_OPERATION_WITNESSES = _OccurrenceWitnessRegistry[
    _FitCommitOperationWitness
]()


def _mint_fit_commit_operation_witness(
    occurrence_identity: str,
) -> _FitCommitOperationWitness:
    return _FIT_COMMIT_OPERATION_WITNESSES.mint(
        _FitCommitOperationWitness,
        occurrence_identity,
    )


def _bind_fit_commit_operation_witness(
    witness: _FitCommitOperationWitness,
    occurrence_identity: str,
    operation_identity: str,
    *,
    artifact_object: object,
) -> bool:
    return _FIT_COMMIT_OPERATION_WITNESSES.bind(
        witness,
        occurrence_identity,
        operation_identity,
        artifact_object=artifact_object,
    )


@dataclass(frozen=True, slots=True, weakref_slot=True)
class FitCommitOperation:
    """Frozen occurrence-exact result of one attempted atomic fit commit."""

    occurrence_identity: str = field(compare=False)
    accepted_result_identity: str
    accepted_occurrence_identity: str
    problem_identity: str
    terminal: FitCommitTerminal
    receipt: CommitReceipt | None = None
    committed_snapshot: AnalysisValuesSnapshot | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    failure: FitCommitFailure | None = None
    _occurrence_witness: _FitCommitOperationWitness | None = field(
        default=None,
        repr=False,
        compare=False,
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        committed = self.terminal is FitCommitTerminal.COMMITTED
        if (
            not self.occurrence_identity
            or committed != (self.receipt is not None)
            or committed != (self.committed_snapshot is not None)
            or committed == (self.failure is not None)
        ):
            raise ValueError("Fit commit operation lacks its exact terminal payload")
        if self.receipt is not None and (
            self.receipt.accepted_result_identity != self.accepted_result_identity
            or self.receipt.accepted_occurrence_identity
            != self.accepted_occurrence_identity
            or self.receipt.problem_identity != self.problem_identity
        ):
            raise ValueError("Fit commit receipt belongs to another operation anchor")
        identity = _identity(
            "native-fit-commit-operation",
            (
                self.occurrence_identity,
                self.accepted_result_identity,
                self.accepted_occurrence_identity,
                self.problem_identity,
                self.terminal.value,
                None if self.receipt is None else self.receipt.identity,
                None if self.failure is None else self.failure.identity,
            ),
        )
        object.__setattr__(self, "identity", identity)
        if self._occurrence_witness is not None:
            _bind_fit_commit_operation_witness(
                self._occurrence_witness,
                self.occurrence_identity,
                identity,
                artifact_object=self,
            )

    def to_record(self) -> dict[str, object]:
        """Return the typed diagnostic representation of this operation."""
        return {
            "artifact_type": "native_fit_commit_operation",
            "schema_version": 1,
            "identity": self.identity,
            "occurrence_identity": self.occurrence_identity,
            "accepted_result_identity": self.accepted_result_identity,
            "accepted_occurrence_identity": self.accepted_occurrence_identity,
            "problem_identity": self.problem_identity,
            "terminal": self.terminal.value,
            "receipt_identity": None if self.receipt is None else self.receipt.identity,
            "failure": (
                None
                if self.failure is None
                else {
                    "identity": self.failure.identity,
                    "category": self.failure.category.value,
                    "message": self.failure.message,
                }
            ),
        }


def fit_commit_operation_is_authoritative(operation: FitCommitOperation) -> bool:
    """Return whether execution minted this exact commit-operation occurrence."""
    witness = operation._occurrence_witness
    if witness is None:
        return False
    return _FIT_COMMIT_OPERATION_WITNESSES.is_bound(
        witness,
        operation.occurrence_identity,
        operation.identity,
        artifact_object=operation,
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
    execution = invocation.execution_settings
    with native_thread_environment(
        execution.native_thread_env(parallel=execution.is_parallel)
    ):
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


def _materialize_candidate_lifecycle(  # noqa: C901 - ordered lifecycle gate
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    vector: tuple[float, ...],
    invocation_identity: str,
    execution_identity: str | Callable[[EvaluationResult], str],
    candidate_summary: CandidateSummary | None = None,
    expected_request: _CompletedRequest | None = None,
    cancellation: CancellationToken | None = None,
    check_cancellation_before_evaluation: bool,
    check_cancellation_after_validation: bool,
) -> (
    MaterializedDirectTrfCandidate
    | CandidateMaterialization
    | RootMaterializationFailure
):
    """Own fresh binding, evaluation, validation, and immutable evidence."""
    token = CancellationToken() if cancellation is None else cancellation
    static_execution_identity = (
        execution_identity if isinstance(execution_identity, str) else None
    )

    def failed(
        terminal: MaterializationTerminal,
        evaluation_count: int,
        compatibility_identity: str,
        cache_hits: int,
        cache_misses: int,
        failure: TerminalFailure,
    ) -> CandidateMaterialization | RootMaterializationFailure:
        if candidate_summary is None or static_execution_identity is None:
            return RootMaterializationFailure(
                terminal,
                evaluation_count,
                compatibility_identity,
                cache_hits,
                cache_misses,
                failure,
            )
        return CandidateMaterialization(
            uuid4().hex,
            problem.identity,
            invocation_identity,
            static_execution_identity,
            candidate_summary,
            terminal,
            evaluation_count,
            compatibility_identity,
            None,
            cache_hits,
            cache_misses,
            failure,
        )

    if check_cancellation_before_evaluation and token.is_cancelled:
        return failed(
            MaterializationTerminal.CANCELLED,
            0,
            engine.compatibility_identity,
            0,
            0,
            TerminalFailure("cancelled", "Cancellation before root materialization"),
        )
    problem.validate_parameterization(parameterization)
    if engine.plan.identity != problem.evaluation_plan_identity:
        raise DirectTrfConstructionError(
            "Root candidate materialization evaluator belongs to another plan"
        )
    canonical_vector = _canonical_vector(
        vector,
        dimension=len(problem.controlled_ids),
        name="root materialization candidate",
    )
    if candidate_summary is not None and (
        candidate_summary.vector != canonical_vector
        or static_execution_identity is None
    ):
        raise DirectTrfConstructionError(
            "Root materialization summary requires one static execution identity"
        )
    try:
        evaluator = engine.new_evaluator()
    except KeyboardInterrupt:
        return failed(
            MaterializationTerminal.INTERRUPTED,
            0,
            engine.compatibility_identity,
            0,
            0,
            TerminalFailure(
                "interrupted",
                "KeyboardInterrupt during root evaluator binding",
            ),
        )
    except Exception as error:  # noqa: BLE001 - binding failures fail closed
        return failed(
            MaterializationTerminal.FAILURE,
            0,
            engine.compatibility_identity,
            0,
            0,
            TerminalFailure(
                "materialization_binding_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    try:
        lifecycle = problem.lifecycle_frame(canonical_vector, parameterization)
        frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
        evaluated = evaluator.evaluate(frame)
    except KeyboardInterrupt:
        statistics = evaluator.cache_statistics
        return failed(
            MaterializationTerminal.INTERRUPTED,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure("interrupted", "KeyboardInterrupt during materialization"),
        )
    except Exception as error:  # noqa: BLE001 - root evaluation fails closed
        statistics = evaluator.cache_statistics
        return failed(
            MaterializationTerminal.FAILURE,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "materialization_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    statistics = evaluator.cache_statistics
    if token.is_cancelled:
        return failed(
            MaterializationTerminal.CANCELLED,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "cancelled",
                "Cancellation observed during root materialization",
            ),
        )
    if isinstance(evaluated, EvaluationFailure):
        return failed(
            MaterializationTerminal.FAILURE,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(evaluated.category, evaluated.message, evaluated),
        )
    try:
        chi_square = (
            _validate_materialized_result(
                evaluated,
                expected_request,
                problem,
                evaluator,
            )
            if expected_request is not None
            else canonical_chi_square(evaluated.residuals)
        )
        if (
            evaluated.plan_identity != problem.evaluation_plan_identity
            or evaluated.parameterization_identity
            != problem.evaluator_parameterization_identity
            or evaluated.evaluator_compatibility_identity
            != evaluator.compatibility_identity
            or tuple(evaluated.resolved_values) != problem.commit_scope
            or canonical_vector
            != tuple(
                evaluated.resolved_values[param_id]
                for param_id in problem.controlled_ids
            )
        ):
            raise ValueError("Fresh root materialization identities or scope differ")
        if isinstance(execution_identity, str):
            resolved_execution_identity = execution_identity
        else:
            resolved_execution_identity = execution_identity(evaluated)
        summary = candidate_summary or CandidateSummary(
            canonical_vector,
            chi_square,
            0,
        )
        if summary.chi_square != chi_square:
            raise ValueError("Fresh root materialization objective differs")
    except KeyboardInterrupt:
        return failed(
            MaterializationTerminal.INTERRUPTED,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "interrupted",
                "KeyboardInterrupt during materialization validation",
            ),
        )
    except Exception as error:  # noqa: BLE001 - root validation fails closed
        return failed(
            MaterializationTerminal.FAILURE,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "materialization_validation_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    if check_cancellation_after_validation and token.is_cancelled:
        return failed(
            MaterializationTerminal.CANCELLED,
            1,
            evaluator.compatibility_identity,
            statistics.hits,
            statistics.misses,
            TerminalFailure(
                "cancelled",
                "Cancellation observed during root materialization",
            ),
        )
    materialization = CandidateMaterialization(
        uuid4().hex,
        problem.identity,
        invocation_identity,
        resolved_execution_identity,
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
        resolved_execution_identity,
        materialization,
        canonical_vector,
        chi_square,
        evaluated,
    )


def materialize_root_candidate(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    vector: tuple[float, ...],
    invocation_identity: str,
    execution_identity: str | Callable[[EvaluationResult], str],
    candidate_summary: CandidateSummary | None = None,
    expected_request: _CompletedRequest | None = None,
    cancellation: CancellationToken | None = None,
) -> (
    MaterializedDirectTrfCandidate
    | CandidateMaterialization
    | RootMaterializationFailure
):
    """Freshly materialize one validated root vector under Direct-TRF ownership."""
    if not problem.acceptance_authority:
        raise DirectTrfConstructionError(
            "Root candidate materialization requires an authoritative problem"
        )
    return _materialize_candidate_lifecycle(
        problem,
        parameterization,
        engine,
        vector=vector,
        invocation_identity=invocation_identity,
        execution_identity=execution_identity,
        candidate_summary=candidate_summary,
        expected_request=expected_request,
        cancellation=cancellation,
        check_cancellation_before_evaluation=True,
        check_cancellation_after_validation=False,
    )


def _materialize_candidate(
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    request: _CompletedRequest,
    cancellation: CancellationToken,
) -> MaterializedDirectTrfCandidate | DirectTrfOutcome:
    materialized = _materialize_candidate_lifecycle(
        problem,
        parameterization,
        engine,
        vector=request.summary.vector,
        invocation_identity=invocation.identity,
        execution_identity=execution.identity,
        candidate_summary=request.summary,
        expected_request=request,
        cancellation=cancellation,
        check_cancellation_before_evaluation=False,
        check_cancellation_after_validation=True,
    )
    if isinstance(materialized, RootMaterializationFailure):
        raise DirectTrfConstructionError(
            "Direct TRF materialization must retain failure evidence"
        )
    if isinstance(materialized, CandidateMaterialization):
        if materialized.terminal is MaterializationTerminal.INTERRUPTED:
            raise DirectTrfInterrupted(execution, materialized)
        terminal = (
            DirectTrfOutcomeTerminal.CANCELLED
            if materialized.terminal is MaterializationTerminal.CANCELLED
            else DirectTrfOutcomeTerminal.MATERIALIZATION_FAILURE
        )
        return DirectTrfOutcome(terminal, execution, materialized)
    return materialized


def _accepted_fit_evidence(
    *,
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    materialization: CandidateMaterialization,
    vector: tuple[float, ...],
    chi_square: float,
    evaluation_result: EvaluationResult,
    origin_context_identity: str,
) -> AcceptedFitResult:
    """Construct evidence only from one fresh complete root materialization."""
    if not problem.acceptance_authority:
        raise DirectTrfConstructionError(
            "A derived component problem has no acceptance authority"
        )
    if (
        materialization.problem_identity != problem.identity
        or materialization.invocation_identity != invocation_identity
        or materialization.execution_identity != execution_identity
        or materialization.candidate.vector != vector
        or materialization.candidate.chi_square != chi_square
        or materialization.evaluation_identity != evaluation_result.identity
        or evaluation_result.plan_identity != problem.evaluation_plan_identity
        or evaluation_result.parameterization_identity
        != problem.evaluator_parameterization_identity
        or tuple(evaluation_result.resolved_values) != problem.commit_scope
    ):
        raise DirectTrfConstructionError(
            "Fresh materialization is incompatible with root acceptance"
        )
    commit_items = tuple(
        (param_id, evaluation_result.resolved_values[param_id])
        for param_id in problem.commit_scope
    )
    occurrence_identity = uuid4().hex
    return AcceptedFitResult(
        occurrence_identity,
        problem.identity,
        invocation_identity,
        execution_identity,
        materialization.identity,
        problem.parameterization_identity,
        problem.evaluator_parameterization_identity,
        problem.source_snapshot.occurrence_identity,
        problem.source_snapshot.revision,
        problem.controlled_ids,
        vector,
        chi_square,
        evaluation_result,
        problem.commit_scope,
        commit_items,
        origin_context_identity,
        occurrence_witness=_mint_accepted_occurrence_witness(occurrence_identity),
    )


def _snapshot_commit_context_identity(snapshot: AnalysisValuesSnapshot) -> str:
    return _identity(
        "native-analysis-values-commit-context",
        (
            snapshot.occurrence_identity,
            snapshot.model_identity,
            snapshot.definitions_identity,
            snapshot.configuration_identity,
            snapshot.revision,
            tuple(
                (param_id, _float_token(value)) for param_id, value in snapshot.items()
            ),
        ),
    )


def committed_values_identity(
    snapshot: AnalysisValuesSnapshot,
    scope: tuple[str, ...],
) -> str:
    """Return the receipt identity for committed values in an exact scope."""
    return _identity(
        "native-committed-values",
        tuple((param_id, _float_token(snapshot[param_id])) for param_id in scope),
    )


def _mint_fit_commit_authority(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
) -> LiveFitCommitAuthority:
    if (
        accepted.problem_identity != problem.identity
        or accepted.source_occurrence_identity
        != problem.source_snapshot.occurrence_identity
        or accepted.source_revision != problem.source_snapshot.revision
        or not accepted.origin_context_identity
    ):
        raise DirectTrfConstructionError(
            "Fit commit authority evidence differs from its root context"
        )
    authority = object.__new__(LiveFitCommitAuthority)
    binding = _CommitAuthorityBinding(
        accepted,
        accepted.identity,
        accepted.occurrence_identity,
        problem.identity,
        _snapshot_commit_context_identity(problem.source_snapshot),
        problem.source_snapshot.occurrence_identity,
        problem.source_snapshot.revision,
        accepted.origin_context_identity,
    )
    with _LIVE_FIT_COMMIT_AUTHORITIES_LOCK:
        _LIVE_FIT_COMMIT_AUTHORITIES[authority] = binding
    return authority


def _direct_fit_commit_origin(
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    materialization: CandidateMaterialization,
    evaluation_result: EvaluationResult,
) -> str:
    return _identity(
        "native-direct-trf-fit-commit-origin",
        (
            problem.identity,
            invocation_identity,
            execution_identity,
            materialization.identity,
            evaluation_result.identity,
        ),
    )


def accept_materialized_fit(
    *,
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    materialization: CandidateMaterialization,
    vector: tuple[float, ...],
    chi_square: float,
    evaluation_result: EvaluationResult,
) -> tuple[AcceptedFitResult, LiveFitCommitAuthority]:
    """Accept a fresh Direct TRF fit and mint its process-owned capability."""
    origin_context_identity = _direct_fit_commit_origin(
        problem,
        invocation_identity,
        execution_identity,
        materialization,
        evaluation_result,
    )
    accepted = _accepted_fit_evidence(
        problem=problem,
        invocation_identity=invocation_identity,
        execution_identity=execution_identity,
        materialization=materialization,
        vector=vector,
        chi_square=chi_square,
        evaluation_result=evaluation_result,
        origin_context_identity=origin_context_identity,
    )
    authority = _mint_fit_commit_authority(accepted, problem)
    return accepted, authority


def _accept_materialized_fit_for_derived_workflow(
    *,
    problem: OptimizationProblem,
    invocation_identity: str,
    execution_identity: str,
    materialization: CandidateMaterialization,
    vector: tuple[float, ...],
    chi_square: float,
    evaluation_result: EvaluationResult,
    authority_context_identity: str,
) -> AcceptedFitResult:
    """Construct derived evidence after workflow validation, before authority."""
    return _accepted_fit_evidence(
        problem=problem,
        invocation_identity=invocation_identity,
        execution_identity=execution_identity,
        materialization=materialization,
        vector=vector,
        chi_square=chi_square,
        evaluation_result=evaluation_result,
        origin_context_identity=authority_context_identity,
    )


def _grant_derived_fit_commit_authority(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
) -> LiveFitCommitAuthority:
    """Mint authority for one exact, already validated derived evidence object."""
    return _mint_fit_commit_authority(accepted, problem)


def _accepted_outcome(
    candidate: MaterializedDirectTrfCandidate,
    execution: DirectTrfExecution,
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
) -> DirectTrfOutcome:
    accepted, authority = accept_materialized_fit(
        problem=problem,
        invocation_identity=invocation.identity,
        execution_identity=execution.identity,
        materialization=candidate.materialization,
        vector=candidate.vector,
        chi_square=candidate.chi_square,
        evaluation_result=candidate.evaluation_result,
    )
    return DirectTrfOutcome(
        DirectTrfOutcomeTerminal.ACCEPTED,
        execution,
        candidate.materialization,
        accepted,
        authority,
    )


def _execute_direct_trf_attempt(
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    token: CancellationToken,
) -> tuple[DirectTrfExecution, _CompletedRequest] | DirectTrfOutcome:
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
    request = solver.request
    execution = _execution(
        occurrence_identity,
        problem,
        invocation,
        live,
        DirectTrfTerminal.CONVERGED,
        preflight_identity,
        final_candidate=request.summary,
        backend=solver.backend,
    )
    if token.is_cancelled:
        return DirectTrfOutcome(DirectTrfOutcomeTerminal.CANCELLED, execution)
    return execution, request


def execute_direct_trf(
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> DirectTrfOutcome:
    """Execute one bounded Direct-TRF attempt and fresh acceptance materialization."""
    if not problem.acceptance_authority:
        raise DirectTrfConstructionError(
            "A derived component problem has no acceptance authority"
        )
    _validate_execution_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    attempt = _execute_direct_trf_attempt(
        problem,
        invocation,
        parameterization,
        engine,
        token,
    )
    if isinstance(attempt, DirectTrfOutcome):
        return attempt
    execution, request = attempt
    materialized = _materialize_candidate(
        execution,
        problem,
        invocation,
        parameterization,
        engine,
        request,
        token,
    )
    if isinstance(materialized, DirectTrfOutcome):
        return materialized
    return _accepted_outcome(
        materialized,
        execution,
        problem,
        invocation,
    )


def execute_direct_trf_candidate(
    problem: OptimizationProblem,
    invocation: DirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> DirectTrfCandidateOutcome:
    """Run Direct TRF for one component without acceptance or commit authority."""
    _validate_execution_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    attempt = _execute_direct_trf_attempt(
        problem,
        invocation,
        parameterization,
        engine,
        token,
    )
    if isinstance(attempt, DirectTrfOutcome):
        terminal = {
            DirectTrfOutcomeTerminal.SOLVER_UNSUCCESSFUL: (
                DirectTrfCandidateTerminal.SOLVER_UNSUCCESSFUL
            ),
            DirectTrfOutcomeTerminal.MATERIALIZATION_FAILURE: (
                DirectTrfCandidateTerminal.MATERIALIZATION_FAILURE
            ),
            DirectTrfOutcomeTerminal.CANCELLED: DirectTrfCandidateTerminal.CANCELLED,
        }[attempt.terminal]
        return DirectTrfCandidateOutcome(
            terminal,
            attempt.execution,
            attempt.materialization,
        )
    execution, request = attempt
    materialized = _materialize_candidate(
        execution,
        problem,
        invocation,
        parameterization,
        engine,
        request,
        token,
    )
    if isinstance(materialized, DirectTrfOutcome):
        terminal = (
            DirectTrfCandidateTerminal.CANCELLED
            if materialized.terminal is DirectTrfOutcomeTerminal.CANCELLED
            else DirectTrfCandidateTerminal.MATERIALIZATION_FAILURE
        )
        return DirectTrfCandidateOutcome(
            terminal,
            execution,
            materialized.materialization,
        )
    return DirectTrfCandidateOutcome(
        DirectTrfCandidateTerminal.SUCCESS,
        execution,
        materialized.materialization,
        materialized,
    )


def _validate_derived_candidate_for_root(
    problem: OptimizationProblem,
    candidate: MaterializedDirectTrfCandidate,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> _CompletedRequest:
    if not problem.acceptance_authority:
        raise DirectTrfConstructionError(
            "Derived root materialization requires an authoritative root problem"
        )
    problem.validate_parameterization(parameterization)
    if engine.plan.identity != problem.evaluation_plan_identity:
        raise DirectTrfConstructionError(
            "Derived root materialization evaluator belongs to another plan"
        )
    materialization = candidate.materialization
    evaluation = candidate.evaluation_result
    try:
        candidate_chi_square = canonical_chi_square(evaluation.residuals)
    except (TypeError, ValueError, ObjectiveScalarizationError) as error:
        raise DirectTrfConstructionError(
            "Derived candidate has invalid objective evidence"
        ) from error
    if (
        materialization.terminal is not MaterializationTerminal.SUCCESS
        or materialization.problem_identity != candidate.problem_identity
        or materialization.invocation_identity != candidate.invocation_identity
        or materialization.execution_identity != candidate.execution_identity
        or materialization.candidate.vector != candidate.vector
        or materialization.candidate.chi_square != candidate.chi_square
        or materialization.evaluation_identity != evaluation.identity
        or evaluation.plan_identity != problem.evaluation_plan_identity
        or evaluation.parameterization_identity
        != problem.evaluator_parameterization_identity
        or evaluation.evaluator_compatibility_identity != engine.compatibility_identity
        or tuple(evaluation.resolved_values) != problem.commit_scope
        or candidate.vector
        != tuple(
            evaluation.resolved_values[param_id] for param_id in problem.controlled_ids
        )
        or candidate.chi_square != candidate_chi_square
    ):
        raise DirectTrfConstructionError(
            "Derived candidate provenance is incompatible with its root"
        )
    return _CompletedRequest(
        materialization.candidate,
        tuple(float(value) for value in evaluation.residuals),
        evaluation.identity,
    )


def _materialize_derived_direct_trf_candidate_for_root(
    problem: OptimizationProblem,
    candidate: MaterializedDirectTrfCandidate,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> MaterializedDirectTrfCandidate | CandidateMaterialization:
    """Freshly evaluate workflow-validated candidate evidence at its root."""
    request = _validate_derived_candidate_for_root(
        problem,
        candidate,
        parameterization,
        engine,
    )
    materialized = materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=candidate.vector,
        invocation_identity=candidate.invocation_identity,
        execution_identity=candidate.execution_identity,
        candidate_summary=CandidateSummary(
            candidate.vector,
            candidate.chi_square,
            0,
        ),
        expected_request=request,
        cancellation=cancellation,
    )
    if isinstance(materialized, RootMaterializationFailure):
        raise DirectTrfConstructionError(
            "Referenced root materialization must retain failure evidence"
        )
    return materialized


def _consume_fit_commit_authority(
    accepted: AcceptedFitResult,
    authority: LiveFitCommitAuthority,
    problem: OptimizationProblem,
) -> None:
    if not isinstance(authority, LiveFitCommitAuthority):
        raise DirectTrfConstructionError(
            "Accepted result lacks its exact live Direct TRF commit authority"
        )
    with _LIVE_FIT_COMMIT_AUTHORITIES_LOCK:
        binding = _LIVE_FIT_COMMIT_AUTHORITIES.get(authority)
        if (
            binding is None
            or binding.accepted_result is not accepted
            or binding.accepted_result_identity != accepted.identity
            or binding.accepted_occurrence_identity != accepted.occurrence_identity
            or binding.problem_identity != problem.identity
            or binding.snapshot_context_identity
            != _snapshot_commit_context_identity(problem.source_snapshot)
            or binding.source_occurrence_identity
            != problem.source_snapshot.occurrence_identity
            or binding.source_revision != problem.source_snapshot.revision
            or binding.origin_context_identity != accepted.origin_context_identity
        ):
            raise DirectTrfConstructionError(
                "Accepted result lacks its exact live Direct TRF commit authority"
            )
        del _LIVE_FIT_COMMIT_AUTHORITIES[authority]


def commit_accepted_fit(
    accepted: AcceptedFitResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Commit evidence once through its exact process-owned live capability."""
    if not problem.acceptance_authority:
        raise DirectTrfConstructionError(
            "A derived component problem has no commit authority"
        )
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
    _consume_fit_commit_authority(accepted, authority, problem)
    committed = analysis_values.commit(
        dict(accepted.commit_items),
        expected=problem.source_snapshot,
        scope=problem.commit_scope,
    )
    committed_value_identity = committed_values_identity(
        committed,
        problem.commit_scope,
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


def execute_fit_commit(
    accepted: AcceptedFitResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> FitCommitOperation:
    """Attempt one fit commit and freeze either its receipt or typed failure."""
    occurrence_identity = uuid4().hex
    witness = _mint_fit_commit_operation_witness(occurrence_identity)
    try:
        receipt = commit_accepted_fit(
            accepted,
            authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=analysis_values,
        )
    except AnalysisValuesCommitError as error:
        if isinstance(error, StaleAnalysisValuesError):
            category = FitCommitFailureCategory.STALE_REVISION
        elif isinstance(error, IncompatibleAnalysisValuesError):
            category = FitCommitFailureCategory.INCOMPATIBLE_STATE
        elif isinstance(error, InvalidAnalysisValuesCommitError):
            category = FitCommitFailureCategory.INVALID_CANDIDATE
        else:  # pragma: no cover - closed subclasses above are exhaustive today
            raise
        return FitCommitOperation(
            occurrence_identity,
            accepted.identity,
            accepted.occurrence_identity,
            problem.identity,
            FitCommitTerminal.FAILED,
            failure=FitCommitFailure(category, str(error)),
            _occurrence_witness=witness,
        )
    return FitCommitOperation(
        occurrence_identity,
        accepted.identity,
        accepted.occurrence_identity,
        problem.identity,
        FitCommitTerminal.COMMITTED,
        receipt=receipt,
        committed_snapshot=analysis_values.snapshot(),
        _occurrence_witness=witness,
    )
