"""Exact private solver coordinates for model-owned relaxation feasibility."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np

from chemex.parameters.parameterization import (
    ActiveParameterization,
    BinaryExpression,
    CompiledConstraint,
    FunctionExpression,
    IndependentValueFrame,
    LiteralExpression,
    ReferenceExpression,
    UnaryExpression,
)
from chemex.parameters.relaxation import RelaxationPsdBlock
from chemex.typing import Array

# Eigenvalue and Schur checks operate on tiny (2x2/3x3) symmetric binary64
# matrices. These tolerances admit eigvalsh round-off at a semidefinite boundary
# while remaining six orders tighter than the scientific regression tolerances.
_PSD_EIGENVALUE_RTOL = 1.0e-11
_SCHUR_RANGE_RTOL = 1.0e-11
_SCHUR_RANGE_ATOL = 1.0e-12
_DIFFERENTIAL_STEP_RELATIVE = math.sqrt(np.finfo(np.float64).eps)


class ScientificFeasibilityError(ValueError):
    """A public ChemEx state violates a model-owned scientific domain."""


class FeasibleCoordinateConstructionError(ValueError):
    """The active FIT/FIX partition lacks an exact supported chart."""


def _identity(records: object) -> str:
    encoded = json.dumps(records, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _matrix(block: RelaxationPsdBlock, values: Mapping[str, float]) -> Array:
    matrix = np.diag([values[param_id] for param_id in block.diagonal_ids])
    for row, column, param_id in block.off_diagonal_ids:
        matrix[row, column] = matrix[column, row] = values[param_id]
    return matrix


def _validate_blocks(
    blocks: Sequence[RelaxationPsdBlock],
    values: Mapping[str, float],
) -> None:
    for block in blocks:
        matrix = _matrix(block, values)
        scale = max(1.0, float(np.max(np.abs(matrix))))
        smallest = float(np.linalg.eigvalsh(matrix)[0])
        if smallest < -_PSD_EIGENVALUE_RTOL * scale:
            raise ScientificFeasibilityError(
                "Relaxation block is not positive semidefinite: "
                f"domain={block.domain_id!r}, state={block.state!r}, "
                f"diagonal_ids={block.diagonal_ids!r}, "
                f"off_diagonal_ids={block.off_diagonal_ids!r}, "
                f"minimum_eigenvalue={smallest!r}"
            )


def _conditional_interval(
    block: RelaxationPsdBlock,
    cross_id: str,
    values: Mapping[str, float],
) -> tuple[float, float]:
    entries = [
        (row, column)
        for row, column, param_id in block.off_diagonal_ids
        if param_id == cross_id
    ]
    if len(entries) != 1:
        raise FeasibleCoordinateConstructionError(
            f"Relaxation cross-rate {cross_id!r} has an ambiguous block entry"
        )
    row, column = entries[0]
    diagonals = [values[param_id] for param_id in block.diagonal_ids]
    if any(value < 0.0 for value in diagonals):
        return (1.0, -1.0)
    pair_limit = math.sqrt(diagonals[row] * diagonals[column])
    lower, upper = -pair_limit, pair_limit
    if len(diagonals) == 2:
        return lower, upper
    other = 3 - row - column
    pivot = diagonals[other]
    fixed = {
        frozenset((left, right)): values[param_id]
        for left, right, param_id in block.off_diagonal_ids
        if param_id != cross_id
    }
    left_other = fixed.get(frozenset((row, other)), 0.0)
    right_other = fixed.get(frozenset((column, other)), 0.0)
    if pivot == 0.0:
        if left_other != 0.0 or right_other != 0.0:
            return (1.0, -1.0)
        return lower, upper
    left_margin = diagonals[row] * pivot - left_other * left_other
    right_margin = diagonals[column] * pivot - right_other * right_other
    if left_margin < 0.0 or right_margin < 0.0:
        return (1.0, -1.0)
    radius = math.sqrt(max(0.0, left_margin * right_margin)) / pivot
    center = left_other * right_other / pivot
    return max(lower, center - radius), min(upper, center + radius)


def _conditional_diagonal_floor(
    block: RelaxationPsdBlock,
    diagonal_id: str,
    values: Mapping[str, float],
) -> float:
    """Return the exact Schur-complement floor for one free diagonal."""
    index = block.diagonal_ids.index(diagonal_id)
    matrix = _matrix(block, values)
    keep = [position for position in range(len(matrix)) if position != index]
    principal = matrix[np.ix_(keep, keep)]
    scale = max(1.0, float(np.max(np.abs(principal))))
    if float(np.linalg.eigvalsh(principal)[0]) < -_PSD_EIGENVALUE_RTOL * scale:
        return math.inf
    coupling = matrix[index, keep]
    inverse = np.linalg.pinv(principal, hermitian=True)
    if not np.allclose(
        principal @ inverse @ coupling,
        coupling,
        rtol=_SCHUR_RANGE_RTOL,
        atol=_SCHUR_RANGE_ATOL * scale,
    ):
        return math.inf
    return max(0.0, float(coupling @ inverse @ coupling))


@dataclass(frozen=True, slots=True)
class FeasiblePoint:
    """One decoded public state that is safe to send to scientific evaluation."""

    frame: IndependentValueFrame
    vector: tuple[float, ...]


class RateFloorKind(StrEnum):
    """Closed exact two-mode PSD floor forms supported by the private chart."""

    OTHER_DIAGONAL = "other-diagonal"
    DERIVED_AFFINE = "derived-affine"
    REPEATED_DIAGONAL = "repeated-diagonal"


@dataclass(frozen=True, slots=True)
class RateFloor:
    """One exact lower envelope imposed on a controlled public rate."""

    rate_id: str
    kind: RateFloorKind
    other_id: str
    cross_id: str


def _rate_floor(specification: RateFloor, values: Mapping[str, float]) -> float:
    cross = values[specification.cross_id]
    if specification.kind is RateFloorKind.REPEATED_DIAGONAL:
        return abs(cross)
    other = values[specification.other_id]
    if specification.kind is RateFloorKind.OTHER_DIAGONAL:
        if other <= 0.0:
            return 0.0 if cross == 0.0 else math.inf
        return cross * cross / other
    intercept = other - values[specification.rate_id]
    discriminant = math.sqrt(intercept * intercept + 4.0 * cross * cross)
    return 0.5 * (-intercept + discriminant)


@dataclass(frozen=True, slots=True)
class FeasibleCoordinates:
    """Role-aware exact chart over the represented relaxation PSD domains."""

    parameterization: ActiveParameterization = field(repr=False, compare=False)
    base_frame: IndependentValueFrame = field(repr=False, compare=False)
    controlled_ids: tuple[str, ...]
    public_lower_bounds: tuple[float, ...]
    public_upper_bounds: tuple[float, ...]
    rate_excess_ids: tuple[tuple[str, str], ...]
    rate_floors: tuple[RateFloor, ...]
    static_rate_floors: tuple[tuple[str, float], ...]
    cross_rate_ids: tuple[str, ...]
    blocks: tuple[RelaxationPsdBlock, ...]
    solver_start: tuple[float, ...]
    solver_lower_bounds: tuple[float, ...]
    solver_upper_bounds: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                (
                    self.parameterization.evaluator_identity,
                    self.controlled_ids,
                    self.public_lower_bounds,
                    self.public_upper_bounds,
                    self.rate_excess_ids,
                    tuple(
                        (item.rate_id, item.kind, item.other_id, item.cross_id)
                        for item in self.rate_floors
                    ),
                    self.static_rate_floors,
                    self.cross_rate_ids,
                    tuple(block.domain_id for block in self.blocks),
                    self.solver_start,
                    self.solver_lower_bounds,
                    self.solver_upper_bounds,
                )
            ),
        )

    @property
    def has_coordinate_transform(self) -> bool:
        return bool(self.rate_excess_ids or self.rate_floors or self.cross_rate_ids)

    @property
    def is_noop(self) -> bool:
        return (
            not self.rate_excess_ids
            and not self.rate_floors
            and not self.cross_rate_ids
            and self.solver_lower_bounds == self.public_lower_bounds
            and self.solver_upper_bounds == self.public_upper_bounds
        )

    @property
    def supports_box_only_algorithms(self) -> bool:
        """Return whether public and private coordinates differ only by box bounds."""
        return not self.has_coordinate_transform

    @property
    def solver_domain(
        self,
    ) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
        """Return the private start, lower bounds, and upper bounds."""
        return self.solver_start, self.solver_lower_bounds, self.solver_upper_bounds

    @property
    def solver_bounds(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return private lower and upper bounds for box-only consumers."""
        return self.solver_lower_bounds, self.solver_upper_bounds

    @property
    def domain_parameter_groups(self) -> tuple[frozenset[str], ...]:
        """Return parameter-ID hyperedges that must remain in one fit component."""
        return tuple(
            frozenset(
                (
                    *block.diagonal_ids,
                    *(item[2] for item in block.off_diagonal_ids),
                )
            )
            for block in self.blocks
        )

    def frame_with_updates(self, updates: Mapping[str, float]) -> IndependentValueFrame:
        """Return a compatible source frame with independent public updates."""
        return self.base_frame.with_updates(updates)

    def derive(
        self,
        frame: IndependentValueFrame,
        controlled_ids: tuple[str, ...],
        lower_bounds: tuple[float, ...],
        upper_bounds: tuple[float, ...],
    ) -> FeasibleCoordinates:
        """Compile a role-aware child chart without exposing model internals."""
        return compile_feasible_coordinates(
            self.parameterization,
            frame,
            controlled_ids,
            lower_bounds,
            upper_bounds,
        )

    def decode(  # noqa: C901 - complete feasibility decode and verification
        self,
        solver_vector: Sequence[float],
    ) -> FeasiblePoint:
        if len(solver_vector) != len(self.controlled_ids):
            raise ValueError("Feasible solver vector has the wrong dimension")
        public = {
            param_id: float(value)
            for param_id, value in zip(self.controlled_ids, solver_vector, strict=True)
        }
        rate_excess = {item[0] for item in self.rate_excess_ids}
        floor_rate_ids = {item.rate_id for item in self.rate_floors}
        static_rate_floors = dict(self.static_rate_floors)
        cross_rates = set(self.cross_rate_ids)
        raw_updates = {
            param_id: value
            for param_id, value in public.items()
            if param_id not in rate_excess
            and param_id not in floor_rate_ids
            and param_id not in cross_rates
        }
        frame = self.base_frame.with_updates(raw_updates)
        resolved = self.parameterization.resolve(frame)
        lower_by_id = dict(
            zip(self.controlled_ids, self.public_lower_bounds, strict=True)
        )
        upper_by_id = dict(
            zip(self.controlled_ids, self.public_upper_bounds, strict=True)
        )
        excess_by_rate = {
            rate_id: tuple(
                diagonal_id
                for candidate, diagonal_id in self.rate_excess_ids
                if candidate == rate_id
            )
            for rate_id in self.controlled_ids
        }
        for rate_id in self.controlled_ids:
            diagonal_ids = excess_by_rate.get(rate_id, ())
            if not diagonal_ids:
                continue
            floor = max(
                lower_by_id[rate_id],
                static_rate_floors.get(rate_id, -math.inf),
                *(resolved[rate_id] - resolved[item] for item in diagonal_ids),
            )
            public[rate_id] = floor + public[rate_id]
        frame = self.base_frame.with_updates(
            {key: value for key, value in public.items() if key not in cross_rates}
        )
        resolved = self.parameterization.resolve(frame)
        floors_by_rate = {
            rate_id: tuple(item for item in self.rate_floors if item.rate_id == rate_id)
            for rate_id in self.controlled_ids
        }
        for rate_id in self.controlled_ids:
            specifications = floors_by_rate.get(rate_id, ())
            if not specifications:
                continue
            floor = max(
                lower_by_id[rate_id],
                static_rate_floors.get(rate_id, -math.inf),
                *(_rate_floor(item, resolved) for item in specifications),
            )
            if not math.isfinite(floor):
                raise ScientificFeasibilityError(
                    f"Relaxation rate {rate_id!r} has no finite PSD floor"
                )
            public[rate_id] = floor + public[rate_id]
            frame = frame.with_updates({rate_id: public[rate_id]})
            resolved = self.parameterization.resolve(frame)
        for cross_id in self.cross_rate_ids:
            lower = lower_by_id[cross_id]
            upper = upper_by_id[cross_id]
            values = dict(resolved)
            for block in self.blocks:
                if any(item[2] == cross_id for item in block.off_diagonal_ids):
                    block_lower, block_upper = _conditional_interval(
                        block, cross_id, values
                    )
                    lower = max(lower, block_lower)
                    upper = min(upper, block_upper)
            if lower > upper:
                raise ScientificFeasibilityError(
                    f"Relaxation cross-rate {cross_id!r} has an empty feasible interval"
                )
            rho = public[cross_id]
            public[cross_id] = 0.5 * ((1.0 - rho) * lower + (1.0 + rho) * upper)
            frame = frame.with_updates({cross_id: public[cross_id]})
            resolved = self.parameterization.resolve(frame)
        vector = tuple(public[param_id] for param_id in self.controlled_ids)
        for value, lower, upper in zip(
            vector,
            self.public_lower_bounds,
            self.public_upper_bounds,
            strict=True,
        ):
            if not lower <= value <= upper:
                raise ScientificFeasibilityError(
                    "Feasible-coordinate decode violated a public scalar bound"
                )
        _validate_blocks(self.blocks, resolved)
        return FeasiblePoint(frame, vector)

    def differential(self, solver_vector: Sequence[float]) -> Array:
        """Return d(public coordinates)/d(private solver coordinates)."""
        point = np.asarray(solver_vector, dtype=np.float64)
        dimension = len(point)
        result = np.empty((dimension, dimension), dtype=np.float64)
        eps = _DIFFERENTIAL_STEP_RELATIVE
        lower = np.asarray(self.solver_lower_bounds)
        upper = np.asarray(self.solver_upper_bounds)
        for column in range(dimension):
            step = eps * max(1.0, abs(point[column]))
            left = point.copy()
            right = point.copy()
            if (
                point[column] - step >= lower[column]
                and point[column] + step <= upper[column]
            ):
                left[column] -= step
                right[column] += step
                denominator = 2.0 * step
            elif point[column] + step <= upper[column]:
                right[column] += step
                denominator = step
            else:
                left[column] -= step
                denominator = step
            left_value = np.asarray(self.decode(tuple(left)).vector)
            right_value = np.asarray(self.decode(tuple(right)).vector)
            result[:, column] = (right_value - left_value) / denominator
        return result


def _reference_coefficient(
    expression: object,
    param_id: str,
    constraints: Mapping[str, CompiledConstraint] | None = None,
    seen: frozenset[str] = frozenset(),
) -> float | None:
    if isinstance(expression, ReferenceExpression):
        reference_id = expression.param_id
        if reference_id == param_id:
            return 1.0
        constraint = None if constraints is None else constraints.get(reference_id)
        if constraint is None or reference_id in seen:
            return 0.0
        return _reference_coefficient(
            constraint.expression,
            param_id,
            constraints,
            seen | {reference_id},
        )
    if isinstance(expression, BinaryExpression) and expression.operator in {
        "add",
        "subtract",
    }:
        left = _reference_coefficient(expression.left, param_id, constraints, seen)
        right = _reference_coefficient(expression.right, param_id, constraints, seen)
        if left is None or right is None:
            return None
        return left + right if expression.operator == "add" else left - right
    return None


def _is_additive_reference_expression(expression: object) -> bool:
    if isinstance(expression, ReferenceExpression):
        return True
    return (
        isinstance(expression, BinaryExpression)
        and expression.operator == "add"
        and _is_additive_reference_expression(expression.left)
        and _is_additive_reference_expression(expression.right)
    )


def _derived_diagonal_transforms(
    parameterization: ActiveParameterization,
    controlled: set[str],
    blocks: Sequence[RelaxationPsdBlock],
) -> tuple[tuple[str, str], ...]:
    diagonal_ids = {param_id for block in blocks for param_id in block.diagonal_ids}
    transforms: set[tuple[str, str]] = set()
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    for diagonal_id in diagonal_ids:
        constraint = constraints.get(diagonal_id)
        if constraint is None:
            continue
        expression = constraint.expression
        if _is_additive_reference_expression(expression):
            continue
        dependencies = controlled & set(constraint.dependencies)
        supported = {
            param_id
            for param_id in dependencies
            if _reference_coefficient(expression, param_id, constraints) == 1.0
        }
        if len(supported) == 1:
            transforms.add((supported.pop(), diagonal_id))
        elif dependencies:
            raise FeasibleCoordinateConstructionError(
                "A fitted relaxation diagonal has an unsupported affine controller: "
                f"{diagonal_id!r}"
            )
    return tuple(sorted(transforms))


def _effective_blocks(
    parameterization: ActiveParameterization,
    blocks: Sequence[RelaxationPsdBlock],
) -> tuple[RelaxationPsdBlock, ...]:
    """Map expression-shared cross-rates to their independent controller IDs."""
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }

    def controller(param_id: str) -> str:
        seen: set[str] = set()
        while param_id in constraints and param_id not in seen:
            seen.add(param_id)
            expression = constraints[param_id].expression
            if not isinstance(expression, ReferenceExpression):
                break
            param_id = expression.param_id
        return param_id

    return tuple(
        RelaxationPsdBlock(
            block.domain_id,
            block.state,
            tuple(controller(param_id) for param_id in block.diagonal_ids),
            tuple(
                (row, column, controller(param_id))
                for row, column, param_id in block.off_diagonal_ids
            ),
        )
        for block in blocks
    )


def active_relaxation_blocks(
    parameterization: ActiveParameterization,
) -> tuple[RelaxationPsdBlock, ...]:
    """Return represented PSD blocks after active-scope controller resolution."""
    active = set(parameterization.scope_ids)
    return _effective_blocks(
        parameterization,
        tuple(
            block
            for block in parameterization.program.relaxation_domains.blocks
            if {*block.diagonal_ids, *(item[2] for item in block.off_diagonal_ids)}
            <= active
        ),
    )


def validate_relaxation_state(
    parameterization: ActiveParameterization,
    values: Mapping[str, float],
) -> None:
    """Reject an inadmissible public state before scientific kernel execution."""
    _validate_blocks(active_relaxation_blocks(parameterization), values)


def relaxation_state_is_on_boundary(
    parameterization: ActiveParameterization,
    values: Mapping[str, float],
    *,
    controlled_ids: Sequence[str] | None = None,
) -> bool:
    """Return whether a represented PSD block has a numerically active boundary."""
    controlled = set(controlled_ids or ())
    for block in active_relaxation_blocks(parameterization):
        if controlled and not any(
            _controlled_dependencies(parameterization, param_id, controlled)
            for param_id in (
                *block.diagonal_ids,
                *(item[2] for item in block.off_diagonal_ids),
            )
        ):
            continue
        matrix = _matrix(block, values)
        scale = max(1.0, float(np.max(np.abs(matrix))))
        if float(np.linalg.eigvalsh(matrix)[0]) <= _PSD_EIGENVALUE_RTOL * scale:
            return True
    return False


def _controlled_dependencies(
    parameterization: ActiveParameterization,
    param_id: str,
    controlled: set[str],
    seen: frozenset[str] = frozenset(),
) -> frozenset[str]:
    if param_id in controlled:
        return frozenset((param_id,))
    if param_id in seen:
        return frozenset()
    constraint = next(
        (
            item
            for item in parameterization.program.constraints
            if item.target_id == param_id
        ),
        None,
    )
    if constraint is None:
        return frozenset()
    dependencies: set[str] = set()
    for dependency in constraint.dependencies:
        dependencies.update(
            _controlled_dependencies(
                parameterization,
                dependency,
                controlled,
                seen | {param_id},
            )
        )
    return frozenset(dependencies)


def _is_intrinsic_rate_derivation(
    parameterization: ActiveParameterization,
    param_id: str,
    seen: frozenset[str] = frozenset(),
) -> bool:
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    constraint = constraints.get(param_id)
    if (
        constraint is None
        or constraint.source not in {"baseline", "model"}
        or param_id in seen
    ):
        return False

    def intrinsic(expression: object) -> bool:
        if isinstance(expression, FunctionExpression):
            return True
        if isinstance(expression, LiteralExpression):
            return True
        if isinstance(expression, ReferenceExpression):
            return _is_intrinsic_rate_derivation(
                parameterization,
                expression.param_id,
                seen | {param_id},
            )
        if isinstance(expression, UnaryExpression):
            return intrinsic(expression.operand)
        if isinstance(expression, BinaryExpression):
            return intrinsic(expression.left) and intrinsic(expression.right)
        return False

    return intrinsic(constraint.expression)


def _is_intrinsic_relaxation_block(
    parameterization: ActiveParameterization,
    block: RelaxationPsdBlock,
) -> bool:
    return all(
        _is_intrinsic_rate_derivation(parameterization, param_id)
        for param_id in (
            *block.diagonal_ids,
            *(item[2] for item in block.off_diagonal_ids),
        )
    )


def _derived_pair_floor(
    parameterization: ActiveParameterization,
    rate_id: str,
    other_diagonal_id: str,
    cross_id: str,
) -> RateFloor | None:
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    constraint = constraints.get(other_diagonal_id)
    if constraint is None:
        return None
    expression = constraint.expression
    if _reference_coefficient(expression, rate_id, constraints) == 1.0:
        return RateFloor(
            rate_id,
            RateFloorKind.DERIVED_AFFINE,
            other_diagonal_id,
            cross_id,
        )
    return None


def compile_feasible_coordinates(  # noqa: C901 - complete role-aware chart
    parameterization: ActiveParameterization,
    frame: IndependentValueFrame,
    controlled_ids: tuple[str, ...],
    lower_bounds: tuple[float, ...],
    upper_bounds: tuple[float, ...],
) -> FeasibleCoordinates:
    """Compile an exact supported chart for the active relaxation domains."""
    blocks = active_relaxation_blocks(parameterization)
    chart_blocks = tuple(
        block
        for block in blocks
        if not _is_intrinsic_relaxation_block(parameterization, block)
    )
    controlled = set(controlled_ids)
    rate_excess_ids = _derived_diagonal_transforms(
        parameterization,
        controlled,
        chart_blocks,
    )
    cross_rate_ids: list[str] = []
    for block in chart_blocks:
        fitted = {
            param_id
            for _row, _column, param_id in block.off_diagonal_ids
            if param_id in controlled
        }
        if len(fitted) > 1:
            raise FeasibleCoordinateConstructionError(
                "Jointly fitted cross-rates in one relaxation block require a "
                f"whole-block chart: domain={block.domain_id!r}, ids={sorted(fitted)!r}"
            )
        cross_rate_ids.extend(fitted)
    cross_ids = tuple(
        param_id for param_id in controlled_ids if param_id in cross_rate_ids
    )
    rate_floors: list[RateFloor] = []
    static_rate_floors: dict[str, float] = {}
    resolved_initial = parameterization.resolve(frame)

    def register_rate_floor(specification: RateFloor) -> None:
        dependencies = _controlled_dependencies(
            parameterization,
            specification.other_id,
            controlled,
        ) | _controlled_dependencies(
            parameterization,
            specification.cross_id,
            controlled,
        )
        if dependencies <= {specification.rate_id}:
            floor = _rate_floor(specification, resolved_initial)
            static_rate_floors[specification.rate_id] = max(
                static_rate_floors.get(specification.rate_id, 0.0),
                floor,
            )
        else:
            rate_floors.append(specification)

    for block in chart_blocks:
        if len(block.diagonal_ids) == 3:
            fitted_crosses = {
                item[2] for item in block.off_diagonal_ids if item[2] in controlled
            }
            diagonal_dependencies = set().union(
                *(
                    _controlled_dependencies(parameterization, item, controlled)
                    for item in block.diagonal_ids
                )
            )
            if fitted_crosses and diagonal_dependencies:
                raise FeasibleCoordinateConstructionError(
                    "Jointly fitted longitudinal cross-rates and diagonals require "
                    f"a whole-block chart: domain={block.domain_id!r}"
                )
            if fitted_crosses:
                continue
            resolved = parameterization.resolve(frame)
            if not any(resolved[item[2]] != 0.0 for item in block.off_diagonal_ids):
                continue
            fitted_diagonals = tuple(
                item for item in block.diagonal_ids if item in controlled
            )
            if len(diagonal_dependencies) != 1 or len(fitted_diagonals) != 1:
                if diagonal_dependencies:
                    raise FeasibleCoordinateConstructionError(
                        "Coupled longitudinal diagonal dependencies require a "
                        f"whole-block chart: domain={block.domain_id!r}"
                    )
                continue
            rate_id = fitted_diagonals[0]
            if any(
                rate_id in _controlled_dependencies(parameterization, item, controlled)
                for item in block.diagonal_ids
                if item != rate_id
            ):
                raise FeasibleCoordinateConstructionError(
                    "A fitted diagonal changes several entries of a coupled "
                    f"longitudinal block: domain={block.domain_id!r}"
                )
            floor = _conditional_diagonal_floor(block, rate_id, resolved)
            if not math.isfinite(floor):
                raise ScientificFeasibilityError(
                    f"Longitudinal block {block.domain_id!r} has no finite PSD floor"
                )
            static_rate_floors[rate_id] = max(
                static_rate_floors.get(rate_id, 0.0), floor
            )
            continue
        if len(block.diagonal_ids) != 2 or len(block.off_diagonal_ids) != 1:
            continue
        _row, _column, cross_id = block.off_diagonal_ids[0]
        if cross_id in controlled:
            continue
        resolved = parameterization.resolve(frame)
        if resolved[cross_id] == 0.0:
            continue
        if block.diagonal_ids[0] == block.diagonal_ids[1]:
            rate_id = block.diagonal_ids[0]
            if rate_id in controlled:
                register_rate_floor(
                    RateFloor(
                        rate_id,
                        RateFloorKind.REPEATED_DIAGONAL,
                        rate_id,
                        cross_id,
                    )
                )
            continue
        fitted_diagonals = [item for item in block.diagonal_ids if item in controlled]
        dependencies = {
            item: _controlled_dependencies(parameterization, item, controlled)
            for item in block.diagonal_ids
        }
        if fitted_diagonals:
            if len(fitted_diagonals) > 1:
                raise FeasibleCoordinateConstructionError(
                    "Two independently fitted relaxation diagonals with a held "
                    f"cross-rate need a joint chart: domain={block.domain_id!r}"
                )
            rate_id = next(item for item in controlled_ids if item in fitted_diagonals)
            other_id = next(
                param_id for param_id in block.diagonal_ids if param_id != rate_id
            )
            derived = _derived_pair_floor(
                parameterization,
                rate_id,
                other_id,
                cross_id,
            )
            if derived is not None:
                register_rate_floor(derived)
            elif rate_id not in dependencies[other_id]:
                register_rate_floor(
                    RateFloor(
                        rate_id,
                        RateFloorKind.OTHER_DIAGONAL,
                        other_id,
                        cross_id,
                    )
                )
            else:
                raise FeasibleCoordinateConstructionError(
                    "A fitted relaxation diagonal has an unsupported coupled form: "
                    f"domain={block.domain_id!r}, rate={rate_id!r}"
                )
        elif set().union(*dependencies.values()):
            raise FeasibleCoordinateConstructionError(
                "A controlled parameter can invalidate a derived held relaxation "
                f"diagonal: domain={block.domain_id!r}"
            )
    floor_rate_ids = {item.rate_id for item in rate_floors}
    rate_excess_ids = tuple(
        item for item in rate_excess_ids if item[0] not in floor_rate_ids
    )
    public_start = dict(frame.ordered_items())
    resolved_start = parameterization.resolve(frame)
    _validate_blocks(blocks, resolved_start)
    rate_map = {
        param_id: tuple(
            diagonal_id
            for rate_id, diagonal_id in rate_excess_ids
            if rate_id == param_id
        )
        for param_id in controlled_ids
    }
    floors_by_rate = {
        param_id: tuple(item for item in rate_floors if item.rate_id == param_id)
        for param_id in controlled_ids
    }
    solver_start: list[float] = []
    solver_lower: list[float] = []
    solver_upper: list[float] = []
    lower_by_id = dict(zip(controlled_ids, lower_bounds, strict=True))
    upper_by_id = dict(zip(controlled_ids, upper_bounds, strict=True))
    for param_id in controlled_ids:
        if rate_map[param_id]:
            floor = max(
                lower_by_id[param_id],
                static_rate_floors.get(param_id, -math.inf),
                *(
                    resolved_start[param_id] - resolved_start[diagonal_id]
                    for diagonal_id in rate_map[param_id]
                ),
            )
            excess = public_start[param_id] - floor
            if excess < 0.0:
                raise ScientificFeasibilityError(
                    f"Derived relaxation diagonal for {param_id!r} is negative"
                )
            if upper_by_id[param_id] < np.finfo(np.float64).max:
                raise FeasibleCoordinateConstructionError(
                    f"Finite upper bound for transformed relaxation rate {param_id!r} "
                    "is not yet supported exactly"
                )
            solver_start.append(excess)
            solver_lower.append(0.0)
            solver_upper.append(upper_by_id[param_id])
        elif floors_by_rate[param_id]:
            if upper_by_id[param_id] < np.finfo(np.float64).max:
                raise FeasibleCoordinateConstructionError(
                    f"Finite upper bound for transformed relaxation rate {param_id!r} "
                    "is not yet supported exactly"
                )
            floor = max(
                lower_by_id[param_id],
                static_rate_floors.get(param_id, -math.inf),
                *(
                    _rate_floor(item, resolved_start)
                    for item in floors_by_rate[param_id]
                ),
            )
            if not math.isfinite(floor):
                raise ScientificFeasibilityError(
                    f"Relaxation rate {param_id!r} has no finite PSD floor"
                )
            excess = public_start[param_id] - floor
            if excess < 0.0:
                raise ScientificFeasibilityError(
                    f"Relaxation rate {param_id!r} is below its PSD floor"
                )
            solver_start.append(excess)
            solver_lower.append(0.0)
            solver_upper.append(upper_by_id[param_id])
        elif param_id in cross_ids:
            lower = lower_by_id[param_id]
            upper = upper_by_id[param_id]
            values = dict(resolved_start)
            for block in blocks:
                if any(item[2] == param_id for item in block.off_diagonal_ids):
                    block_lower, block_upper = _conditional_interval(
                        block, param_id, values
                    )
                    lower = max(lower, block_lower)
                    upper = min(upper, block_upper)
            width = upper - lower
            value = public_start[param_id]
            if width < 0.0 or not lower <= value <= upper:
                raise ScientificFeasibilityError(
                    f"Initial cross-rate {param_id!r} is outside its PSD interval"
                )
            rho = 0.0 if width == 0.0 else (2.0 * value - lower - upper) / width
            solver_start.append(rho)
            solver_lower.append(-1.0)
            solver_upper.append(1.0)
        else:
            solver_start.append(public_start[param_id])
            solver_lower.append(
                max(lower_by_id[param_id], static_rate_floors.get(param_id, -math.inf))
            )
            solver_upper.append(upper_by_id[param_id])
    return FeasibleCoordinates(
        parameterization,
        frame,
        controlled_ids,
        lower_bounds,
        upper_bounds,
        rate_excess_ids,
        tuple(rate_floors),
        tuple(sorted(static_rate_floors.items())),
        cross_ids,
        blocks,
        tuple(solver_start),
        tuple(solver_lower),
        tuple(solver_upper),
    )
