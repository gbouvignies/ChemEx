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
from chemex.parameters.relaxation import (
    RelaxationPsdBlock,
    relaxation_blocks_identity,
)
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
class _FeasibilityProjectionProvenance:
    """Private immutable closure proof compiled with one feasibility chart."""

    controlled_domain_groups: tuple[frozenset[str], ...]
    has_root_projection_authority: bool


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
    _projection_provenance: _FeasibilityProjectionProvenance | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    solver_start: tuple[float, ...]
    solver_lower_bounds: tuple[float, ...]
    solver_upper_bounds: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(("unsealed-feasibility-chart",)),
        )

    @property
    def controlled_domain_groups(self) -> tuple[frozenset[str], ...]:
        """Return private root-compiled controlled-domain closure proofs."""
        provenance = self._projection_provenance
        if provenance is None:
            raise FeasibleCoordinateConstructionError(
                "Feasibility chart lacks construction-owned projection provenance"
            )
        return provenance.controlled_domain_groups

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
    ) -> FeasibleCoordinates | None:
        """Compile a role-aware child chart without exposing model internals."""
        _ = self.controlled_domain_groups
        return compile_feasible_coordinates(
            self.parameterization,
            frame,
            controlled_ids,
            lower_bounds,
            upper_bounds,
        )

    def project_component(
        self,
        frame: IndependentValueFrame,
        controlled_ids: tuple[str, ...],
        lower_bounds: tuple[float, ...],
        upper_bounds: tuple[float, ...],
    ) -> FeasibleCoordinates | None:
        """Project an exact root chart onto one closed fit component."""
        provenance = self._projection_provenance
        if provenance is None or not provenance.has_root_projection_authority:
            raise FeasibleCoordinateConstructionError(
                "Component feasibility projection requires the compiled root chart"
            )
        selected = set(controlled_ids)
        root_indices = {
            param_id: index for index, param_id in enumerate(self.controlled_ids)
        }
        if (
            frame != self.base_frame
            or tuple(
                param_id for param_id in self.controlled_ids if param_id in selected
            )
            != controlled_ids
            or len(selected) != len(controlled_ids)
            or len(self.controlled_domain_groups) != len(self.blocks)
            or any(
                dependencies & selected and not dependencies <= selected
                for dependencies in self.controlled_domain_groups
            )
        ):
            raise FeasibleCoordinateConstructionError(
                "Component feasibility projection is not a closed root-chart subset"
            )
        expected_lower = tuple(
            self.public_lower_bounds[root_indices[param_id]]
            for param_id in controlled_ids
        )
        expected_upper = tuple(
            self.public_upper_bounds[root_indices[param_id]]
            for param_id in controlled_ids
        )
        if lower_bounds != expected_lower or upper_bounds != expected_upper:
            raise FeasibleCoordinateConstructionError(
                "Component feasibility projection changed root public bounds"
            )
        projected = _seal_projection_provenance(
            type(self)(
                self.parameterization,
                frame,
                controlled_ids,
                lower_bounds,
                upper_bounds,
                tuple(item for item in self.rate_excess_ids if item[0] in selected),
                tuple(item for item in self.rate_floors if item.rate_id in selected),
                tuple(item for item in self.static_rate_floors if item[0] in selected),
                tuple(
                    param_id for param_id in self.cross_rate_ids if param_id in selected
                ),
                self.blocks,
                tuple(
                    self.solver_start[root_indices[param_id]]
                    for param_id in controlled_ids
                ),
                tuple(
                    self.solver_lower_bounds[root_indices[param_id]]
                    for param_id in controlled_ids
                ),
                tuple(
                    self.solver_upper_bounds[root_indices[param_id]]
                    for param_id in controlled_ids
                ),
            ),
            tuple(
                dependencies & selected
                for dependencies in self.controlled_domain_groups
            ),
            has_root_projection_authority=False,
        )
        return None if projected.is_noop else projected

    def decode(  # noqa: C901 - complete feasibility decode and verification
        self,
        solver_vector: Sequence[float],
    ) -> FeasiblePoint:
        _ = self.controlled_domain_groups
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
        _ = self.controlled_domain_groups
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


@dataclass(frozen=True, slots=True)
class _AffineExpression:
    """Exact scalar constant, linear coefficients, and reference provenance."""

    constant: float
    coefficients: Mapping[str, float]
    reference_ids: frozenset[str]


def _scaled_affine(
    expression: _AffineExpression,
    factor: float,
) -> _AffineExpression | None:
    constant = factor * expression.constant
    coefficients: dict[str, float] = {}
    for param_id, coefficient in expression.coefficients.items():
        scaled = factor * coefficient
        if scaled != 0.0:
            coefficients[param_id] = scaled
    if not math.isfinite(constant) or not all(
        math.isfinite(coefficient) for coefficient in coefficients.values()
    ):
        return None
    return _AffineExpression(constant, coefficients, expression.reference_ids)


def _combined_affine(
    left: _AffineExpression,
    right: _AffineExpression,
    right_factor: float,
) -> _AffineExpression | None:
    constant = left.constant + right_factor * right.constant
    coefficients = dict(left.coefficients)
    for param_id, coefficient in right.coefficients.items():
        combined = coefficients.get(param_id, 0.0) + right_factor * coefficient
        if combined == 0.0:
            coefficients.pop(param_id, None)
        else:
            coefficients[param_id] = combined
    if not math.isfinite(constant) or not all(
        math.isfinite(coefficient) for coefficient in coefficients.values()
    ):
        return None
    return _AffineExpression(
        constant,
        coefficients,
        left.reference_ids | right.reference_ids,
    )


def _binary_affine_expression(
    expression: BinaryExpression,
    constraints: Mapping[str, CompiledConstraint] | None,
    seen: frozenset[str],
) -> _AffineExpression | None:
    left = _affine_expression(expression.left, constraints, seen)
    right = _affine_expression(expression.right, constraints, seen)
    if left is None or right is None:
        return None
    if expression.operator in {"add", "subtract"}:
        factor = 1.0 if expression.operator == "add" else -1.0
        return _combined_affine(left, right, factor)
    if expression.operator == "multiply":
        if not left.reference_ids:
            return _scaled_affine(right, left.constant)
        if not right.reference_ids:
            return _scaled_affine(left, right.constant)
        return None
    if expression.operator == "divide":
        if right.reference_ids or right.constant == 0.0:
            return None
        return _scaled_affine(left, 1.0 / right.constant)
    return None


def _affine_expression(
    expression: object,
    constraints: Mapping[str, CompiledConstraint] | None = None,
    seen: frozenset[str] = frozenset(),
) -> _AffineExpression | None:
    if isinstance(expression, LiteralExpression):
        if not math.isfinite(expression.value):
            return None
        return _AffineExpression(expression.value, {}, frozenset())
    if isinstance(expression, ReferenceExpression):
        reference_id = expression.param_id
        constraint = None if constraints is None else constraints.get(reference_id)
        if constraint is None:
            return _AffineExpression(
                0.0,
                {reference_id: 1.0},
                frozenset((reference_id,)),
            )
        if reference_id in seen:
            return None
        return _affine_expression(
            constraint.expression,
            constraints,
            seen | {reference_id},
        )
    if isinstance(expression, UnaryExpression):
        operand = _affine_expression(expression.operand, constraints, seen)
        if operand is None:
            return None
        factor = 1.0 if expression.operator == "positive" else -1.0
        return _scaled_affine(operand, factor)
    if isinstance(expression, BinaryExpression):
        return _binary_affine_expression(expression, constraints, seen)
    return None


def _reference_coefficient(
    expression: object,
    param_id: str,
    constraints: Mapping[str, CompiledConstraint] | None = None,
) -> float | None:
    affine = _affine_expression(expression, constraints)
    if affine is None:
        return None
    return affine.coefficients.get(param_id, 0.0)


def _is_model_additive_reference_expression(
    expression: object,
    constraints: Mapping[str, CompiledConstraint],
    seen: frozenset[str] = frozenset(),
) -> bool:
    if isinstance(expression, ReferenceExpression):
        reference_id = expression.param_id
        constraint = constraints.get(reference_id)
        if constraint is None:
            return True
        if constraint.source not in {"baseline", "model"} or reference_id in seen:
            return False
        return _is_model_additive_reference_expression(
            constraint.expression,
            constraints,
            seen | {reference_id},
        )
    return (
        isinstance(expression, BinaryExpression)
        and expression.operator == "add"
        and _is_model_additive_reference_expression(
            expression.left,
            constraints,
            seen,
        )
        and _is_model_additive_reference_expression(
            expression.right,
            constraints,
            seen,
        )
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
        if constraint.source in {
            "baseline",
            "model",
        } and _is_model_additive_reference_expression(expression, constraints):
            continue
        dependencies = _controlled_dependencies(
            parameterization,
            diagonal_id,
            controlled,
        )
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


def _seal_projection_provenance(
    chart: FeasibleCoordinates,
    controlled_domain_groups: tuple[frozenset[str], ...],
    *,
    has_root_projection_authority: bool,
) -> FeasibleCoordinates:
    """Seal compiler-owned closure proof and complete scientific chart identity."""
    controlled = frozenset(chart.controlled_ids)
    if len(controlled_domain_groups) != len(chart.blocks) or any(
        not dependencies <= controlled for dependencies in controlled_domain_groups
    ):
        raise FeasibleCoordinateConstructionError(
            "Feasibility dependency provenance differs from its root chart"
        )
    object.__setattr__(
        chart,
        "_projection_provenance",
        _FeasibilityProjectionProvenance(
            controlled_domain_groups,
            has_root_projection_authority,
        ),
    )
    object.__setattr__(
        chart,
        "identity",
        _identity(
            (
                chart.parameterization.evaluator_identity,
                (
                    chart.base_frame.parameterization_identity,
                    chart.base_frame.program_fingerprint,
                    chart.base_frame.occurrence_identity,
                    chart.base_frame.revision,
                    chart.base_frame.ordered_items(),
                ),
                chart.controlled_ids,
                chart.public_lower_bounds,
                chart.public_upper_bounds,
                chart.rate_excess_ids,
                tuple(
                    (item.rate_id, item.kind, item.other_id, item.cross_id)
                    for item in chart.rate_floors
                ),
                chart.static_rate_floors,
                chart.cross_rate_ids,
                relaxation_blocks_identity(chart.blocks),
                tuple(
                    tuple(sorted(dependencies))
                    for dependencies in controlled_domain_groups
                ),
                chart.solver_start,
                chart.solver_lower_bounds,
                chart.solver_upper_bounds,
            )
        ),
    )
    return chart


def _controlled_domain_groups(
    parameterization: ActiveParameterization,
    blocks: Sequence[RelaxationPsdBlock],
    controlled_ids: Sequence[str],
) -> tuple[frozenset[str], ...]:
    """Resolve every PSD domain to immutable root-controlled dependencies."""
    controlled = frozenset(controlled_ids)
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    cache: dict[str, frozenset[str]] = {}

    def resolve(param_id: str, seen: frozenset[str] = frozenset()) -> frozenset[str]:
        if param_id in cache:
            return cache[param_id]
        if param_id in controlled:
            result = frozenset((param_id,))
        elif param_id in seen or (constraint := constraints.get(param_id)) is None:
            result = frozenset()
        else:
            result = frozenset(
                dependency
                for source in constraint.dependencies
                for dependency in resolve(source, seen | {param_id})
            )
        cache[param_id] = result
        return result

    return tuple(
        frozenset(
            dependency
            for param_id in (
                *block.diagonal_ids,
                *(item[2] for item in block.off_diagonal_ids),
            )
            for dependency in resolve(param_id)
        )
        for block in blocks
    )


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
) -> FeasibleCoordinates | None:
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
    chart = _seal_projection_provenance(
        FeasibleCoordinates(
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
        ),
        _controlled_domain_groups(parameterization, blocks, controlled_ids),
        has_root_projection_authority=True,
    )
    return None if chart.is_noop else chart
