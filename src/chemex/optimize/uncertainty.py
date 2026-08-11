"""ChemEx-owned accepted-point uncertainty evidence qualification (#598).

This module is deliberately separate from backend execution evidence.  It
linearizes the complete native residual map again at one accepted external
coordinate vector, constructs full-rank scaled-SVD covariance, and derives
typed marginal, correlation, and constrained-propagation artifacts.

The entry point remains an isolated qualification seam; production fitting and
reporting do not consume these artifacts yet.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, fields, is_dataclass
from enum import StrEnum
from itertools import pairwise
from numbers import Real
from typing import Literal, Protocol, cast
from uuid import uuid4

import numpy as np
from scipy.linalg import svd

from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    BinaryExpression,
    CompiledConstraint,
    FunctionExpression,
    LiteralExpression,
    ParameterRole,
    ReferenceExpression,
    ScalarExpression,
    UnaryExpression,
    scientific_callable_fingerprint,
)
from chemex.typing import Array

_SCHEMA_VERSION = 1
_RESIDUAL_LINEARIZATION_VERSION = "external-real-finite-difference-v1"
_CONSTRAINT_LINEARIZATION_VERSION = "constraint-forward-chain-v1"
_COVARIANCE_VERSION = "full-rank-scaled-svd-cholesky-v2"
_FACTOR_VERSION = "canonical-information-cholesky-gram-v2"
_REDUCTION_VERSION = "fixed-pairwise-binary64-v1"
_MARGINAL_VERSION = "covariance-diagonal-square-root-v1"
_CORRELATION_VERSION = "ordered-double-division-v1"
_PROPAGATION_VERSION = "constraint-gram-factor-v1"
_EPSILON = float(np.finfo(np.float64).eps)
_NOMINAL_STEP_FACTOR = _EPSILON ** (1.0 / 3.0)
_BOUNDARY_THRESHOLD = 3.0


class UncertaintyConstructionError(ValueError):
    """Raised only for malformed policy or artifact construction."""


class FunctionPartialFailure(UncertaintyConstructionError):
    """Typed failure retaining a scientific-function partial trajectory."""

    def __init__(
        self,
        category: str,
        message: str,
        trajectory_fingerprint: str,
        termination_details: tuple[str, ...],
    ) -> None:
        super().__init__(message)
        self.category = category
        self.trajectory_fingerprint = trajectory_fingerprint
        self.termination_details = termination_details


class DerivationTermination(RuntimeError):
    """Internal cooperative cancellation/interruption checkpoint signal."""

    def __init__(
        self,
        terminal: OperationTerminal,
        partial_artifact: object | None = None,
    ) -> None:
        super().__init__(terminal.value)
        self.terminal = terminal
        self.partial_artifact = partial_artifact


class ResidualVarianceScaling(StrEnum):
    """Closed interpretation of supplied observation uncertainties."""

    ABSOLUTE_OBSERVATION_UNCERTAINTIES = "absolute_observation_uncertainties"
    ESTIMATED_COMMON_RESIDUAL_VARIANCE = "estimated_common_residual_variance"


class ClaimState(StrEnum):
    """One closed scientific-claim assessment."""

    SATISFIED = "satisfied"
    VIOLATED = "violated"
    INDETERMINATE = "indeterminate"
    NOT_APPLICABLE = "not_applicable"


class OperationTerminal(StrEnum):
    """Lifecycle result of one evidence derivation operation."""

    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


class ParameterUnit(StrEnum):
    """Closed unit vocabulary carried by uncertainty evidence."""

    DIMENSIONLESS = "dimensionless"
    FRACTION = "fraction"
    RATE_PER_SECOND = "s^-1"
    FREQUENCY_HZ = "Hz"
    CHEMICAL_SHIFT_PPM = "ppm"
    TIME_SECONDS = "s"
    TEMPERATURE_CELSIUS = "degC"
    CONCENTRATION_MOLAR = "mol/L"
    ENERGY_PER_MOLE = "J/mol"
    ENTROPY_PER_MOLE_KELVIN = "J/mol/K"


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _record_value(value: object) -> object:
    """Encode typed immutable evidence with canonical binary64 values."""
    if isinstance(value, StrEnum):
        return value.value
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return {"binary64": _float_token(value)}
    if isinstance(value, tuple):
        return [_record_value(item) for item in value]
    if is_dataclass(value) and not isinstance(value, type):
        return {
            "artifact_type": type(value).__name__,
            **{
                item.name: _record_value(getattr(value, item.name))
                for item in fields(value)
                if item.metadata.get("record", True)
            },
        }
    raise TypeError(f"Unsupported uncertainty record value {type(value).__name__}")


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real binary64 scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return 0.0 if result == 0.0 else result


def _float_token(value: float) -> str:
    return float(value).hex()


def _vector_tokens(values: Sequence[float]) -> tuple[str, ...]:
    return tuple(_float_token(value) for value in values)


def _matrix_tokens(
    matrix: Sequence[Sequence[float]],
) -> tuple[tuple[str, ...], ...]:
    return tuple(_vector_tokens(row) for row in matrix)


def _canonical_matrix(
    values: Sequence[Sequence[float]],
    *,
    rows: int,
    columns: int,
    name: str,
) -> tuple[tuple[float, ...], ...]:
    matrix = tuple(
        tuple(
            _finite(value, name=f"{name}[{row_index},{column_index}]")
            for column_index, value in enumerate(row)
        )
        for row_index, row in enumerate(values)
    )
    if len(matrix) != rows or any(len(row) != columns for row in matrix):
        raise ValueError(f"{name} must have shape ({rows}, {columns})")
    return matrix


def _validated_extent(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise UncertaintyConstructionError(f"{name} must be a non-negative integer")
    return value


def _validated_coordinate_contract(
    raw_scales: Sequence[tuple[str, float]],
    raw_units: Sequence[tuple[str, ParameterUnit]],
) -> tuple[
    tuple[tuple[str, float], ...],
    tuple[tuple[str, ParameterUnit], ...],
]:
    scales = tuple(
        (param_id, _finite(value, name=f"scale {param_id!r}"))
        for param_id, value in raw_scales
    )
    ids = tuple(param_id for param_id, _value in scales)
    if (
        not ids
        or any(not param_id for param_id in ids)
        or len(set(ids)) != len(ids)
        or any(value <= 0.0 for _param_id, value in scales)
    ):
        raise UncertaintyConstructionError(
            "Coordinate scales must be unique positive finite ID/value pairs"
        )
    units = tuple(raw_units)
    if tuple(param_id for param_id, _unit in units) != ids or any(
        not unit for _param_id, unit in units
    ):
        raise UncertaintyConstructionError(
            "Coordinate units must explicitly match coordinate-scale order"
        )
    return scales, units


@dataclass(frozen=True, slots=True)
class FunctionFiniteDifferenceCapability:
    """Explicit numerical partial capability for one scientific function output."""

    function_id: str
    component: str | None
    argument_scales: tuple[float, ...]
    output_scale: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        scales = tuple(
            _finite(value, name=f"function argument scale[{index}]")
            for index, value in enumerate(self.argument_scales)
        )
        output_scale = _finite(self.output_scale, name="function output scale")
        if (
            not self.function_id
            or not scales
            or any(value <= 0.0 for value in scales)
            or output_scale <= 0.0
        ):
            raise UncertaintyConstructionError(
                "Function finite-difference scales must be positive and explicit"
            )
        object.__setattr__(self, "argument_scales", scales)
        object.__setattr__(self, "output_scale", output_scale)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-function-finite-difference-capability",
                (
                    self.function_id,
                    self.component,
                    _vector_tokens(scales),
                    _float_token(output_scale),
                    _CONSTRAINT_LINEARIZATION_VERSION,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class FunctionAnalyticPartialCapability:
    """Versioned analytic partials selected for one scientific function output."""

    function_id: str
    component: str | None
    implementation_identity: str
    partials: tuple[Callable[..., float], ...] = field(compare=False, repr=False)
    implementation_fingerprints: tuple[str, ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            not self.function_id
            or not self.implementation_identity
            or not self.partials
            or any(not callable(partial) for partial in self.partials)
        ):
            raise UncertaintyConstructionError(
                "Analytic scientific-function capability is incomplete"
            )
        fingerprints = tuple(
            scientific_callable_fingerprint(partial) for partial in self.partials
        )
        object.__setattr__(self, "implementation_fingerprints", fingerprints)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-function-analytic-partial-capability",
                (
                    self.function_id,
                    self.component,
                    self.implementation_identity,
                    fingerprints,
                    _CONSTRAINT_LINEARIZATION_VERSION,
                ),
            ),
        )

    def validate_implementations(self) -> None:
        """Reject a mutable/rebound analytic implementation after compilation."""
        current = tuple(
            scientific_callable_fingerprint(partial) for partial in self.partials
        )
        if current != self.implementation_fingerprints:
            raise UncertaintyConstructionError(
                "Analytic partial implementation changed after capability compilation"
            )


FunctionLinearizationCapability = (
    FunctionFiniteDifferenceCapability | FunctionAnalyticPartialCapability
)


@dataclass(frozen=True, slots=True)
class CompiledConstraintLinearizationCapabilities:
    """Capability selection frozen against one constraint program and scope."""

    parameterization_identity: str
    constraint_program_identity: str
    output_scope: tuple[str, ...]
    capabilities: tuple[FunctionLinearizationCapability, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        keys = tuple((item.function_id, item.component) for item in self.capabilities)
        if len(set(keys)) != len(keys):
            raise UncertaintyConstructionError(
                "Compiled function capabilities must be unique"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "compiled-constraint-linearization-capabilities",
                (
                    self.parameterization_identity,
                    self.constraint_program_identity,
                    self.output_scope,
                    tuple(item.identity for item in self.capabilities),
                ),
            ),
        )


def _expression_function_keys(
    expression: ScalarExpression,
) -> tuple[tuple[str, str | None], ...]:
    if isinstance(expression, FunctionExpression):
        nested = tuple(
            key
            for argument in expression.arguments
            for key in _expression_function_keys(argument)
        )
        return ((expression.function_id, expression.component), *nested)
    if isinstance(expression, UnaryExpression):
        return _expression_function_keys(expression.operand)
    if isinstance(expression, BinaryExpression):
        return _expression_function_keys(expression.left) + _expression_function_keys(
            expression.right
        )
    return ()


def _expression_function_arities(
    expression: ScalarExpression,
) -> tuple[tuple[tuple[str, str | None], int], ...]:
    if isinstance(expression, FunctionExpression):
        nested = tuple(
            item
            for argument in expression.arguments
            for item in _expression_function_arities(argument)
        )
        return (
            ((expression.function_id, expression.component), len(expression.arguments)),
            *nested,
        )
    if isinstance(expression, UnaryExpression):
        return _expression_function_arities(expression.operand)
    if isinstance(expression, BinaryExpression):
        return _expression_function_arities(
            expression.left
        ) + _expression_function_arities(expression.right)
    return ()


def _required_function_keys(
    parameterization: ActiveParameterization,
    scope: tuple[str, ...],
) -> set[tuple[str, str | None]]:
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    required: set[tuple[str, str | None]] = set()
    visited: set[str] = set()

    def visit(target_id: str) -> None:
        if target_id in visited or target_id not in constraints:
            return
        visited.add(target_id)
        constraint = constraints[target_id]
        for dependency in constraint.dependencies:
            visit(dependency)
        required.update(_expression_function_keys(constraint.expression))

    for target_id in scope:
        visit(target_id)
    return required


def _required_function_arities(
    parameterization: ActiveParameterization,
    scope: tuple[str, ...],
) -> dict[tuple[str, str | None], int]:
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    required: dict[tuple[str, str | None], int] = {}
    visited: set[str] = set()

    def visit(target_id: str) -> None:
        if target_id in visited or target_id not in constraints:
            return
        visited.add(target_id)
        constraint = constraints[target_id]
        for dependency in constraint.dependencies:
            visit(dependency)
        for key, arity in _expression_function_arities(constraint.expression):
            previous = required.setdefault(key, arity)
            if previous != arity:
                raise UncertaintyConstructionError(
                    f"Scientific function {key!r} has inconsistent compiled arity"
                )

    for target_id in scope:
        visit(target_id)
    return required


def compile_constraint_linearization_capabilities(
    parameterization: ActiveParameterization,
    output_scope: Sequence[str],
    capabilities: Sequence[FunctionLinearizationCapability],
) -> CompiledConstraintLinearizationCapabilities:
    """Freeze scientific-function capability selection before derivation."""
    scope = tuple(output_scope)
    if len(set(scope)) != len(scope) or any(
        param_id not in parameterization.scope_ids for param_id in scope
    ):
        raise UncertaintyConstructionError(
            "Compiled constraint-linearization scope is invalid"
        )
    selected = tuple(capabilities)
    selected_keys = {(item.function_id, item.component) for item in selected}

    required_keys = _required_function_keys(parameterization, scope)
    missing = required_keys - selected_keys
    if missing:
        raise UncertaintyConstructionError(
            "Compiled constraint linearization lacks capabilities for "
            f"{sorted(missing, key=repr)!r}"
        )
    required_arities = _required_function_arities(parameterization, scope)
    for capability in selected:
        key = (capability.function_id, capability.component)
        arity = required_arities.get(key)
        if arity is None:
            continue
        if isinstance(capability, FunctionFiniteDifferenceCapability):
            if len(capability.argument_scales) != arity:
                raise UncertaintyConstructionError(
                    f"Numerical capability {key!r} has incompatible arity"
                )
        else:
            capability.validate_implementations()
            if len(capability.partials) != arity:
                raise UncertaintyConstructionError(
                    f"Analytic capability {key!r} has incompatible arity"
                )
            for partial in capability.partials:
                try:
                    inspect.signature(partial).bind(*((0.0,) * arity))
                except TypeError as error:
                    raise UncertaintyConstructionError(
                        f"Analytic partial for {key!r} has incompatible domain arity"
                    ) from error
    return CompiledConstraintLinearizationCapabilities(
        parameterization.identity,
        parameterization.program.fingerprint,
        scope,
        selected,
    )


@dataclass(frozen=True, slots=True)
class UncertaintyPolicy:
    """Explicit resolved v1 numerical policy; no defaults are inferred."""

    calibration_identity: str
    numerical_compatibility_requirement: str
    coordinate_scales: tuple[tuple[str, float], ...]
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    residual_variance_scaling: ResidualVarianceScaling
    relative_step_tolerance: float
    roundoff_multiplier: float
    smaller_step_extent: int
    larger_step_extent: int
    svd_driver: Literal["gesdd", "gesvd"]
    rank_absolute_tolerance: float
    rank_relative_tolerance: float
    weak_relative_tolerance: float
    singular_value_cluster_relative_tolerance: float
    conditioning_limit: float
    correlation_roundoff_multiplier: float
    affine_feasibility_policy: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.calibration_identity:
            raise UncertaintyConstructionError(
                "Uncertainty policy requires a calibration identity"
            )
        if not self.numerical_compatibility_requirement:
            raise UncertaintyConstructionError(
                "Uncertainty policy requires numerical compatibility semantics"
            )
        scales, units = _validated_coordinate_contract(
            self.coordinate_scales,
            self.coordinate_units,
        )
        relative = _finite(
            self.relative_step_tolerance,
            name="relative step tolerance",
        )
        roundoff = _finite(
            self.roundoff_multiplier,
            name="roundoff multiplier",
        )
        rank_absolute = _finite(
            self.rank_absolute_tolerance,
            name="absolute rank tolerance",
        )
        rank_relative = _finite(
            self.rank_relative_tolerance,
            name="relative rank tolerance",
        )
        weak_relative = _finite(
            self.weak_relative_tolerance,
            name="weak-subspace relative tolerance",
        )
        cluster_relative = _finite(
            self.singular_value_cluster_relative_tolerance,
            name="singular-value cluster relative tolerance",
        )
        conditioning = _finite(
            self.conditioning_limit,
            name="conditioning limit",
        )
        correlation_roundoff = _finite(
            self.correlation_roundoff_multiplier,
            name="correlation roundoff multiplier",
        )
        if relative < 0.0 or roundoff < 0.0:
            raise UncertaintyConstructionError(
                "Finite-difference tolerances must be non-negative"
            )
        if (
            rank_absolute < 0.0
            or rank_relative < 0.0
            or weak_relative < rank_relative
            or cluster_relative < 0.0
        ):
            raise UncertaintyConstructionError("Rank thresholds must be non-negative")
        if conditioning <= 0.0 or correlation_roundoff < 0.0:
            raise UncertaintyConstructionError(
                "Conditioning and correlation policies are invalid"
            )
        _validated_extent(
            self.smaller_step_extent,
            name="smaller step extent",
        )
        _validated_extent(
            self.larger_step_extent,
            name="larger step extent",
        )
        if self.svd_driver not in {"gesdd", "gesvd"}:
            raise UncertaintyConstructionError(
                "SVD driver must be explicitly 'gesdd' or 'gesvd'"
            )
        if not self.affine_feasibility_policy:
            raise UncertaintyConstructionError(
                "Affine-feasibility semantics must be explicit"
            )
        object.__setattr__(self, "coordinate_scales", scales)
        object.__setattr__(self, "coordinate_units", units)
        object.__setattr__(self, "relative_step_tolerance", relative)
        object.__setattr__(self, "roundoff_multiplier", roundoff)
        object.__setattr__(self, "rank_absolute_tolerance", rank_absolute)
        object.__setattr__(self, "rank_relative_tolerance", rank_relative)
        object.__setattr__(self, "weak_relative_tolerance", weak_relative)
        object.__setattr__(
            self,
            "singular_value_cluster_relative_tolerance",
            cluster_relative,
        )
        object.__setattr__(self, "conditioning_limit", conditioning)
        object.__setattr__(
            self,
            "correlation_roundoff_multiplier",
            correlation_roundoff,
        )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-uncertainty-policy",
                (
                    self.calibration_identity,
                    self.numerical_compatibility_requirement,
                    tuple(
                        (param_id, _float_token(value)) for param_id, value in scales
                    ),
                    units,
                    self.residual_variance_scaling.value,
                    _float_token(relative),
                    _float_token(roundoff),
                    self.smaller_step_extent,
                    self.larger_step_extent,
                    self.svd_driver,
                    _float_token(rank_absolute),
                    _float_token(rank_relative),
                    _float_token(weak_relative),
                    _float_token(cluster_relative),
                    _float_token(conditioning),
                    _float_token(correlation_roundoff),
                    self.affine_feasibility_policy,
                    _float_token(_BOUNDARY_THRESHOLD),
                    _RESIDUAL_LINEARIZATION_VERSION,
                    _CONSTRAINT_LINEARIZATION_VERSION,
                    _COVARIANCE_VERSION,
                    _REDUCTION_VERSION,
                    _MARGINAL_VERSION,
                    _CORRELATION_VERSION,
                    _PROPAGATION_VERSION,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ClaimAssessment:
    name: str
    state: ClaimState
    detail: str = ""


@dataclass(frozen=True, slots=True)
class EvidenceFailure:
    stage: str
    category: str
    message: str
    source_identity: str
    trajectory_fingerprint: str | None = None
    termination_details: tuple[str, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-uncertainty-failure",
                (
                    self.stage,
                    self.category,
                    self.message,
                    self.source_identity,
                    self.trajectory_fingerprint,
                    self.termination_details,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class DerivationOperation:
    occurrence_identity: str = field(compare=False)
    resolved_environment_identity: str = field(compare=False)
    stage: str
    request_identity: str
    terminal: OperationTerminal
    artifact_identity: str | None = None
    failure: EvidenceFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.resolved_environment_identity:
            raise ValueError("Derivation operation requires resolved environment")
        if self.terminal is OperationTerminal.COMPLETED:
            if self.artifact_identity is None or self.failure is not None:
                raise ValueError("Completed operation lacks one exact artifact")
        elif self.artifact_identity is not None or self.failure is None:
            raise ValueError("Non-completed operation has inconsistent evidence")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-uncertainty-operation",
                (
                    self.stage,
                    self.request_identity,
                    self.terminal.value,
                    self.artifact_identity,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class LinearizationColumn:
    param_id: str
    orientation: str
    nominal_step: float
    represented_displacements: tuple[float, ...]
    fine_estimate: tuple[float, ...]
    companion_estimate: tuple[float, ...]
    discrepancy: float
    derivative_scale: float
    roundoff_allowance: float
    attempt_count: int
    termination: str
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-linearization-column",
                (
                    self.param_id,
                    self.orientation,
                    _float_token(self.nominal_step),
                    _vector_tokens(self.represented_displacements),
                    _vector_tokens(self.fine_estimate),
                    _vector_tokens(self.companion_estimate),
                    _float_token(self.discrepancy),
                    _float_token(self.derivative_scale),
                    _float_token(self.roundoff_allowance),
                    self.attempt_count,
                    self.termination,
                    self.trajectory_fingerprint,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class FunctionPartialDiagnostic:
    """Separate reliability evidence for one scientific-function partial."""

    capability_identity: str
    function_id: str
    component: str | None
    argument_index: int
    method: str
    estimate: float
    companion_estimate: float | None
    represented_displacements: tuple[float, ...]
    discrepancy: float
    derivative_scale: float
    roundoff_allowance: float
    attempt_count: int
    trajectory_fingerprint: str
    claims: tuple[ClaimAssessment, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-function-partial-diagnostic",
                (
                    self.capability_identity,
                    self.function_id,
                    self.component,
                    self.argument_index,
                    self.method,
                    _float_token(self.estimate),
                    None
                    if self.companion_estimate is None
                    else _float_token(self.companion_estimate),
                    _vector_tokens(self.represented_displacements),
                    _float_token(self.discrepancy),
                    _float_token(self.derivative_scale),
                    _float_token(self.roundoff_allowance),
                    self.attempt_count,
                    self.trajectory_fingerprint,
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                ),
            ),
        )


def _validate_accepted_anchor(
    accepted: AcceptedFitResult,
    semantic_identity: str,
    occurrence_identity: str,
) -> None:
    if (
        semantic_identity != accepted.identity
        or occurrence_identity != accepted.occurrence_identity
    ):
        raise ValueError(
            "Uncertainty artifact differs from its exact accepted occurrence"
        )


def _named_claim(claims: Sequence[ClaimAssessment], name: str) -> ClaimAssessment:
    matches = tuple(item for item in claims if item.name == name)
    if len(matches) != 1:
        raise ValueError(f"Evidence must contain exactly one {name!r} claim")
    return matches[0]


@dataclass(frozen=True, slots=True)
class ResidualJacobianEvidence:
    request_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    accepted_evaluation_identity: str
    problem_identity: str
    source_snapshot_occurrence_identity: str
    source_revision: int
    parameterization_identity: str
    evaluator_parameterization_identity: str
    constraint_program_identity: str
    evaluation_plan_identity: str
    policy_identity: str
    calibration_identity: str
    numerical_compatibility_requirement: str
    controlled_ids: tuple[str, ...]
    accepted_vector: tuple[float, ...]
    coordinate_scales: tuple[float, ...]
    matrix: tuple[tuple[float, ...], ...]
    columns: tuple[LinearizationColumn, ...]
    trajectory_fingerprint: str
    complete_reliable: bool = True
    accepted_anchor: AcceptedFitResult = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_problem: OptimizationProblem = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_parameterization: ActiveParameterization = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_engine: EvaluationEngine = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_policy: UncertaintyPolicy = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_accepted_anchor(
            self.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        residual_count = len(self.matrix)
        coordinate_count = len(self.controlled_ids)
        matrix = _canonical_matrix(
            self.matrix,
            rows=residual_count,
            columns=coordinate_count,
            name="residual Jacobian",
        )
        if (
            residual_count < 1
            or coordinate_count < 1
            or len(self.columns) != coordinate_count
            or len(self.coordinate_scales) != coordinate_count
        ):
            raise ValueError("Residual Jacobian has inconsistent scope")
        if (
            self.accepted_vector != self.accepted_anchor.vector
            or self.accepted_evaluation_identity
            != self.accepted_anchor.evaluation_result.identity
            or self.parameterization_identity
            != self.accepted_anchor.parameterization_identity
            or self.evaluator_parameterization_identity
            != self.accepted_anchor.evaluator_parameterization_identity
            or self.problem_identity != self.accepted_anchor.problem_identity
            or self.source_snapshot_occurrence_identity
            != self.accepted_anchor.source_occurrence_identity
            or self.source_revision != self.accepted_anchor.source_revision
            or self.controlled_ids != self.accepted_anchor.controlled_ids
            or tuple(item.param_id for item in self.columns) != self.controlled_ids
            or any(value <= 0.0 for value in self.coordinate_scales)
            or any(len(item.fine_estimate) != residual_count for item in self.columns)
            or any(
                len(item.companion_estimate) != residual_count for item in self.columns
            )
        ):
            raise ValueError("Residual Jacobian provenance is internally inconsistent")
        derived_matrix = tuple(
            tuple(column.fine_estimate[row] for column in self.columns)
            for row in range(residual_count)
        )
        expected_trajectory = _identity(
            "residual-linearization-trajectory",
            tuple(
                (item.param_id, item.trajectory_fingerprint) for item in self.columns
            ),
        )
        if (
            matrix != derived_matrix
            or not self.complete_reliable
            or self.trajectory_fingerprint != expected_trajectory
        ):
            raise ValueError("Residual Jacobian differs from its accepted columns")
        evaluator = self.source_engine.new_evaluator()
        base = _evaluate_vector(
            self.accepted_anchor.vector,
            problem=self.source_problem,
            parameterization=self.source_parameterization,
            evaluator=evaluator,
        )
        if (
            isinstance(base, EvaluationFailure)
            or base.identity != self.accepted_evaluation_identity
        ):
            raise ValueError("Residual Jacobian baseline cannot be replayed")
        base_residuals = tuple(float(value) for value in base.residuals)
        for column, (param_id, scale, declared) in enumerate(
            zip(
                self.controlled_ids,
                self.coordinate_scales,
                self.columns,
                strict=True,
            )
        ):
            replayed, _trajectory, failure = _linearize_residual_column(
                self.accepted_anchor,
                problem=self.source_problem,
                parameterization=self.source_parameterization,
                evaluator=evaluator,
                policy=self.source_policy,
                base_residuals=base_residuals,
                column=column,
                param_id=param_id,
                scale=scale,
                cancellation_probe=None,
            )
            if failure is not None or replayed != declared:
                raise ValueError(
                    "Residual Jacobian stencil evidence cannot be replayed"
                )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-residual-jacobian-evidence",
                (
                    self.request_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.accepted_evaluation_identity,
                    self.problem_identity,
                    self.source_snapshot_occurrence_identity,
                    self.source_revision,
                    self.parameterization_identity,
                    self.evaluator_parameterization_identity,
                    self.constraint_program_identity,
                    self.evaluation_plan_identity,
                    self.policy_identity,
                    self.calibration_identity,
                    self.numerical_compatibility_requirement,
                    self.controlled_ids,
                    _vector_tokens(self.accepted_vector),
                    _vector_tokens(self.coordinate_scales),
                    _matrix_tokens(matrix),
                    tuple(column.identity for column in self.columns),
                    self.trajectory_fingerprint,
                    self.complete_reliable,
                    _RESIDUAL_LINEARIZATION_VERSION,
                ),
            ),
        )

    @property
    def residual_count(self) -> int:
        return len(self.matrix)

    @property
    def coordinate_count(self) -> int:
        return len(self.controlled_ids)


@dataclass(frozen=True, slots=True)
class SingularSubspaceEvidence:
    indices: tuple[int, ...]
    singular_values: tuple[float, ...]
    classification: Literal[
        "isolated_identifiable",
        "isolated_weak",
        "isolated_null",
        "clustered_identifiable",
        "clustered_weak",
        "clustered_null",
    ]
    projector: tuple[tuple[float, ...], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        dimension = len(self.projector)
        projector = _canonical_matrix(
            self.projector,
            rows=dimension,
            columns=dimension,
            name="singular subspace projector",
        )
        if (
            not self.indices
            or len(self.indices) != len(self.singular_values)
            or tuple(sorted(self.indices)) != self.indices
            or (len(self.indices) == 1) != self.classification.startswith("isolated")
        ):
            raise ValueError("Singular subspace evidence has inconsistent membership")
        object.__setattr__(self, "projector", projector)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-singular-invariant-subspace",
                (
                    self.indices,
                    _vector_tokens(self.singular_values),
                    self.classification,
                    _matrix_tokens(projector),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class RankDiagnostic:
    request_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    source_jacobian_identity: str
    controlled_ids: tuple[str, ...]
    singular_values: tuple[float, ...]
    normalized_singular_values: tuple[float, ...]
    scaled_column_norms: tuple[float, ...]
    threshold: float
    weak_threshold: float
    rank: int
    identifiable_projector: tuple[tuple[float, ...], ...]
    null_projector: tuple[tuple[float, ...], ...]
    weak_projector: tuple[tuple[float, ...], ...]
    subspaces: tuple[SingularSubspaceEvidence, ...]
    source_jacobian: ResidualJacobianEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_policy: UncertaintyPolicy = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - complete SVD integrity chain
        source = self.source_jacobian
        _validate_accepted_anchor(
            source.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        dimension = len(self.controlled_ids)
        identifiable = _canonical_matrix(
            self.identifiable_projector,
            rows=dimension,
            columns=dimension,
            name="identifiable-subspace projector",
        )
        null = _canonical_matrix(
            self.null_projector,
            rows=dimension,
            columns=dimension,
            name="null-subspace projector",
        )
        weak = _canonical_matrix(
            self.weak_projector,
            rows=dimension,
            columns=dimension,
            name="weak-subspace projector",
        )
        if (
            self.source_jacobian_identity != source.identity
            or self.controlled_ids != source.controlled_ids
            or len(self.singular_values) != dimension
            or len(self.normalized_singular_values) != dimension
            or len(self.scaled_column_norms) != dimension
            or not 0 <= self.rank <= dimension
            or any(
                value < 0.0 or not math.isfinite(value)
                for value in self.singular_values
            )
            or any(left < right for left, right in pairwise(self.singular_values))
            or self.rank
            != sum(value > self.threshold for value in self.singular_values)
            or self.threshold
            != self.source_policy.rank_absolute_tolerance
            + self.source_policy.rank_relative_tolerance
            * (self.singular_values[0] if self.singular_values else 0.0)
        ):
            raise ValueError("SVD/rank evidence has inconsistent lineage or spectrum")
        largest = self.singular_values[0] if self.singular_values else 0.0
        expected_normalized = tuple(
            0.0 if largest == 0.0 else value / largest for value in self.singular_values
        )
        scaled = (
            np.asarray(source.matrix)
            * np.asarray(source.coordinate_scales)[np.newaxis, :]
        )
        expected_norms = tuple(
            float(np.linalg.norm(scaled[:, index])) for index in range(dimension)
        )
        information_eigenvalues = np.linalg.eigvalsh(scaled.T @ scaled)[::-1]
        with np.errstate(over="ignore", invalid="ignore"):
            expected_squared = np.square(np.asarray(self.singular_values))
        tolerance = 512.0 * _EPSILON * max(1, dimension)
        if (
            self.normalized_singular_values != expected_normalized
            or not np.allclose(
                self.scaled_column_norms, expected_norms, rtol=tolerance, atol=0.0
            )
            or not np.allclose(
                expected_squared,
                information_eigenvalues,
                rtol=tolerance,
                atol=tolerance,
            )
        ):
            raise ValueError("SVD/rank numerical evidence does not derive from J_z")
        identity_matrix = np.eye(dimension)
        for name, candidate in (
            ("identifiable", identifiable),
            ("null", null),
            ("weak", weak),
        ):
            array = np.asarray(candidate)
            if not (
                np.allclose(array, array.T, rtol=0.0, atol=tolerance)
                and np.allclose(array @ array, array, rtol=0.0, atol=tolerance)
            ):
                raise ValueError(f"{name} subspace evidence is not a projector")
        if not np.allclose(
            np.asarray(identifiable) + np.asarray(null),
            identity_matrix,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError("Identifiable and null projectors are not complementary")
        covered = tuple(index for item in self.subspaces for index in item.indices)
        if covered != tuple(range(dimension)):
            raise ValueError("Singular subspaces do not partition parameter order")
        expected_groups = _singular_cluster_indices(
            self.singular_values,
            rank_threshold=self.threshold,
            weak_threshold=self.weak_threshold,
            cluster_relative_tolerance=(
                self.source_policy.singular_value_cluster_relative_tolerance
            ),
        )
        if tuple(item.indices for item in self.subspaces) != expected_groups:
            raise ValueError(
                "Singular subspaces violate the declared clustering policy"
            )
        subspace_sum = np.zeros((dimension, dimension), dtype=np.float64)
        information = scaled.T @ scaled
        for item in self.subspaces:
            if item.singular_values != tuple(
                self.singular_values[index] for index in item.indices
            ):
                raise ValueError("Singular subspace spectrum is inconsistent")
            expected_classification = (
                "isolated_" if len(item.indices) == 1 else "clustered_"
            ) + _singular_direction_class(
                self.singular_values,
                item.indices[0],
                rank_threshold=self.threshold,
                weak_threshold=self.weak_threshold,
            )
            if item.classification != expected_classification:
                raise ValueError("Singular subspace classification is inconsistent")
            item_projector = np.asarray(item.projector)
            if not np.allclose(
                item_projector @ information,
                information @ item_projector,
                rtol=tolerance,
                atol=tolerance,
            ) or not np.allclose(
                np.trace(item_projector @ information),
                sum(value * value for value in item.singular_values),
                rtol=tolerance,
                atol=tolerance,
            ):
                raise ValueError(
                    "Singular projector does not match its invariant spectrum"
                )
            subspace_sum += item_projector
        if not np.allclose(
            subspace_sum,
            identity_matrix,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError("Singular invariant-subspace projectors are incomplete")
        object.__setattr__(self, "identifiable_projector", identifiable)
        object.__setattr__(self, "null_projector", null)
        object.__setattr__(self, "weak_projector", weak)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-rank-diagnostic",
                (
                    self.request_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_jacobian_identity,
                    self.controlled_ids,
                    _vector_tokens(self.singular_values),
                    _vector_tokens(self.normalized_singular_values),
                    _vector_tokens(self.scaled_column_norms),
                    _float_token(self.threshold),
                    _float_token(self.weak_threshold),
                    self.rank,
                    _matrix_tokens(identifiable),
                    _matrix_tokens(null),
                    _matrix_tokens(weak),
                    tuple(item.identity for item in self.subspaces),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class CovarianceEvidence:
    request_identity: str
    source_jacobian_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    problem_identity: str
    parameterization_identity: str
    evaluation_plan_identity: str
    policy_identity: str
    calibration_identity: str
    numerical_compatibility_requirement: str
    controlled_ids: tuple[str, ...]
    units: tuple[ParameterUnit, ...]
    accepted_vector: tuple[float, ...]
    coordinate_scales: tuple[float, ...]
    residual_variance_scaling: ResidualVarianceScaling
    retained_residual_count: int
    controlled_coordinate_count: int
    profiled_normalization_count: int
    nominal_residual_degrees_of_freedom: int
    chi_square: float
    residual_variance_scale: float
    rank: int
    singular_values: tuple[float, ...]
    rank_threshold: float
    jacobian_condition: float
    information_condition: float
    unscaled_covariance: tuple[tuple[float, ...], ...]
    factor: tuple[tuple[float, ...], ...]
    covariance: tuple[tuple[float, ...], ...]
    claims: tuple[ClaimAssessment, ...]
    rank_diagnostic_identity: str
    source_jacobian: ResidualJacobianEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_rank_diagnostic: RankDiagnostic = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    accepted_anchor: AcceptedFitResult = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_problem: OptimizationProblem = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_policy: UncertaintyPolicy = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_engine: EvaluationEngine = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    factorization: str = _FACTOR_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_accepted_anchor(
            self.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        dimension = len(self.controlled_ids)
        unscaled = _canonical_matrix(
            self.unscaled_covariance,
            rows=dimension,
            columns=dimension,
            name="unscaled covariance",
        )
        factor = _canonical_matrix(
            self.factor,
            rows=dimension,
            columns=dimension,
            name="covariance factor",
        )
        covariance = _canonical_matrix(
            self.covariance,
            rows=dimension,
            columns=dimension,
            name="scaled covariance",
        )
        if self.rank != dimension or self.controlled_coordinate_count != dimension:
            raise ValueError("Covariance requires full controlled-coordinate rank")
        source = self.source_jacobian
        rank_source = self.source_rank_diagnostic
        accepted = self.accepted_anchor
        problem = self.source_problem
        policy = self.source_policy
        engine = self.source_engine
        expected_residual_count = source.residual_count
        expected_normalization_count = sum(
            1 for profile in engine.plan.profiles if profile.is_scaled
        )
        expected_degrees = (
            expected_residual_count - dimension - expected_normalization_count
        )
        expected_scale, scale_failure = _residual_variance_scale(
            accepted,
            degrees_of_freedom=expected_degrees,
            policy=policy,
        )
        if scale_failure is not None or expected_scale is None:
            raise ValueError("Covariance has no valid residual-variance derivation")
        if (
            self.source_jacobian_identity != source.identity
            or self.rank_diagnostic_identity != rank_source.identity
            or rank_source.source_jacobian_identity != source.identity
            or self.request_identity != rank_source.request_identity
            or self.problem_identity != problem.identity
            or self.problem_identity != accepted.problem_identity
            or self.parameterization_identity != accepted.parameterization_identity
            or self.evaluation_plan_identity != accepted.evaluation_result.plan_identity
            or self.policy_identity != policy.identity
            or self.calibration_identity != policy.calibration_identity
            or self.numerical_compatibility_requirement
            != policy.numerical_compatibility_requirement
            or self.controlled_ids != source.controlled_ids
            or self.controlled_ids != accepted.controlled_ids
            or self.accepted_vector != accepted.vector
            or self.coordinate_scales != source.coordinate_scales
            or self.units != tuple(unit for _param_id, unit in policy.coordinate_units)
            or self.residual_variance_scaling is not policy.residual_variance_scaling
            or self.retained_residual_count != expected_residual_count
            or self.profiled_normalization_count != expected_normalization_count
            or self.nominal_residual_degrees_of_freedom != expected_degrees
            or self.chi_square != accepted.chi_square
            or self.residual_variance_scale != expected_scale
            or self.rank != rank_source.rank
            or self.singular_values != rank_source.singular_values
            or self.rank_threshold != rank_source.threshold
            or self.factorization != _FACTOR_VERSION
        ):
            raise ValueError("Covariance derivation lineage is internally inconsistent")
        expected_unscaled, expected_factor, expected_covariance = (
            _canonical_covariance_reduction(
                source.matrix,
                source.coordinate_scales,
                expected_scale,
            )
        )
        if (
            unscaled != expected_unscaled
            or factor != expected_factor
            or covariance != expected_covariance
        ):
            raise ValueError("Covariance arrays do not match the canonical reduction")
        expected_jacobian_condition = self.singular_values[0] / self.singular_values[-1]
        if (
            self.jacobian_condition != expected_jacobian_condition
            or self.information_condition
            != expected_jacobian_condition * expected_jacobian_condition
        ):
            raise ValueError("Covariance conditioning claims do not match the spectrum")
        (
            expected_interior,
            expected_boundary,
            expected_affine,
            expected_controlled_affine,
        ) = _boundary_claims(
            accepted,
            problem,
            covariance,
            expected_scale,
            policy.affine_feasibility_policy,
        )
        expected_full_dimensional = _full_dimensional_feasible_interior_claim(
            accepted,
            problem,
            expected_interior,
        )
        expected_claim_states = {
            "AUTHORITATIVE_LINEAGE": ClaimState.SATISFIED,
            "LOCAL_LINEARIZATION_REGULARITY": ClaimState.SATISFIED,
            "EFFECTIVE_OBSERVATION_SUFFICIENCY": ClaimState.SATISFIED,
            "FULL_COLUMN_RANK": ClaimState.SATISFIED,
            "COVARIANCE_ARITHMETIC_INTEGRITY": ClaimState.SATISFIED,
            "FULL_DIMENSIONAL_FEASIBLE_INTERIOR": expected_full_dimensional.state,
            "INTERIOR_POINT": expected_interior.state,
            "BOUNDARY_SEPARATION": expected_boundary.state,
            "AFFINE_FEASIBILITY": expected_affine.state,
            "CONTROLLED_AFFINE_SEPARATION": expected_controlled_affine.state,
            "CONDITIONING_ADEQUACY": (
                ClaimState.SATISFIED
                if expected_jacobian_condition <= policy.conditioning_limit
                else ClaimState.VIOLATED
            ),
            "PROFILED_NORMALIZATION_REGULARITY": (
                ClaimState.SATISFIED
                if _profiled_normalization_regular(accepted, engine)
                else ClaimState.VIOLATED
            ),
        }
        for name, state in expected_claim_states.items():
            if _named_claim(self.claims, name).state is not state:
                raise ValueError(f"Covariance claim {name!r} is inconsistent")
        usable_inputs = (
            expected_claim_states["CONDITIONING_ADEQUACY"],
            expected_full_dimensional.state,
            expected_interior.state,
            expected_boundary.state,
            expected_controlled_affine.state
            if expected_controlled_affine.state is not ClaimState.NOT_APPLICABLE
            else ClaimState.SATISFIED,
            expected_claim_states["PROFILED_NORMALIZATION_REGULARITY"],
            (
                ClaimState.SATISFIED
                if policy.residual_variance_scaling
                is ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
                or expected_scale > 0.0
                else ClaimState.VIOLATED
            ),
        )
        expected_usable = (
            ClaimState.SATISFIED
            if all(state is ClaimState.SATISFIED for state in usable_inputs)
            else (
                ClaimState.INDETERMINATE
                if ClaimState.INDETERMINATE in usable_inputs
                else ClaimState.VIOLATED
            )
        )
        if (
            _named_claim(self.claims, "USABLE_LOCAL_COVARIANCE").state
            is not expected_usable
        ):
            raise ValueError("Covariance usability claim is inconsistent")
        object.__setattr__(self, "unscaled_covariance", unscaled)
        object.__setattr__(self, "factor", factor)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-local-covariance-evidence",
                (
                    self.request_identity,
                    self.source_jacobian_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.problem_identity,
                    self.parameterization_identity,
                    self.evaluation_plan_identity,
                    self.policy_identity,
                    self.calibration_identity,
                    self.numerical_compatibility_requirement,
                    self.controlled_ids,
                    self.units,
                    _vector_tokens(self.accepted_vector),
                    _vector_tokens(self.coordinate_scales),
                    self.residual_variance_scaling.value,
                    self.retained_residual_count,
                    self.controlled_coordinate_count,
                    self.profiled_normalization_count,
                    self.nominal_residual_degrees_of_freedom,
                    _float_token(self.chi_square),
                    _float_token(self.residual_variance_scale),
                    self.rank,
                    _vector_tokens(self.singular_values),
                    _float_token(self.rank_threshold),
                    _float_token(self.jacobian_condition),
                    _float_token(self.information_condition),
                    _matrix_tokens(unscaled),
                    _matrix_tokens(factor),
                    _matrix_tokens(covariance),
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                    self.rank_diagnostic_identity,
                    self.factorization,
                    _COVARIANCE_VERSION,
                    _REDUCTION_VERSION,
                ),
            ),
        )

    def claim(self, name: str) -> ClaimState:
        return next(item.state for item in self.claims if item.name == name)

    @property
    def usable(self) -> bool:
        return self.claim("USABLE_LOCAL_COVARIANCE") is ClaimState.SATISFIED


@dataclass(frozen=True, slots=True)
class ScalarEvidenceEntry:
    param_id: str
    value: float | None
    outcome: str
    raw_value: float | None = None


@dataclass(frozen=True, slots=True)
class MarginalErrorEvidence:
    source_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    source_family: str
    output_ids: tuple[str, ...]
    units: tuple[ParameterUnit, ...]
    entries: tuple[ScalarEvidenceEntry, ...]
    claims: tuple[ClaimAssessment, ...]
    scope_reportable: bool
    source_artifact: CovarianceEvidence | ConstrainedPropagationEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        source = self.source_artifact
        source_output_ids = (
            source.output_ids
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.controlled_ids
        )
        if (
            self.source_identity != source.identity
            or self.accepted_result_identity != source.accepted_result_identity
            or self.accepted_occurrence_identity != source.accepted_occurrence_identity
            or self.output_ids != source_output_ids
        ):
            raise ValueError("Marginal-error evidence has inconsistent source lineage")
        source_units = (
            source.output_units
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.units
        )
        expected_family = (
            "constrained_propagation"
            if isinstance(source, ConstrainedPropagationEvidence)
            else "local_covariance"
        )
        if (
            self.units != source_units
            or self.source_family != expected_family
            or len(self.entries) != len(self.output_ids)
            or self.claims[: len(source.claims)] != source.claims
        ):
            raise ValueError("Marginal-error scope differs from its covariance source")
        source_covariance = source.covariance
        residual_degenerate = (
            source.source_covariance.residual_variance_scale == 0.0
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.residual_variance_scale == 0.0
        )
        expected_entries: list[ScalarEvidenceEntry] = []
        for index, param_id in enumerate(self.output_ids):
            variance = source_covariance[index][index]
            if residual_degenerate and variance == 0.0:
                expected_entries.append(
                    ScalarEvidenceEntry(
                        param_id, None, "RESIDUAL_VARIANCE_DEGENERACY", 0.0
                    )
                )
            elif variance > 0.0 and math.isfinite(variance):
                expected_entries.append(
                    ScalarEvidenceEntry(param_id, math.sqrt(variance), "AVAILABLE")
                )
            else:
                expected_entries.append(
                    ScalarEvidenceEntry(
                        param_id,
                        None,
                        "INVALID_VARIANCE",
                        variance if math.isfinite(variance) else None,
                    )
                )
        if self.entries != tuple(expected_entries):
            raise ValueError("Marginal errors do not derive from the source covariance")
        complete = all(item.value is not None for item in expected_entries)
        if _named_claim(self.claims, "MARGINAL_VARIANCE_VALIDITY").state is not (
            ClaimState.SATISFIED if complete else ClaimState.VIOLATED
        ):
            raise ValueError("Marginal validity claim is inconsistent")
        report_claim = _named_claim(self.claims, "MARGINAL_SCOPE_REPORTABILITY")
        source_reportable = (
            source.source_covariance.usable
            and source.claim("LOCAL_FIRST_ORDER_DEGENERACY") is ClaimState.SATISFIED
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.usable
        )
        expected_reportable = source_reportable and complete
        if (
            report_claim.state
            is not (
                ClaimState.SATISFIED if expected_reportable else ClaimState.VIOLATED
            )
            or self.scope_reportable is not expected_reportable
        ):
            raise ValueError("Marginal reportability claim is inconsistent")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-marginal-error-evidence",
                (
                    self.source_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_family,
                    self.output_ids,
                    self.units,
                    tuple(
                        (
                            item.param_id,
                            None if item.value is None else _float_token(item.value),
                            item.outcome,
                            None
                            if item.raw_value is None
                            else _float_token(item.raw_value),
                        )
                        for item in self.entries
                    ),
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                    self.scope_reportable,
                    _MARGINAL_VERSION,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class CorrelationEntry:
    value: float | None
    outcome: str
    raw_value: float | None = None
    excess: float | None = None
    tolerance: float | None = None


def _expected_correlation_entries(
    covariance: tuple[tuple[float, ...], ...],
    *,
    residual_variance_degenerate: bool,
    policy: UncertaintyPolicy,
) -> tuple[tuple[CorrelationEntry, ...], ...]:
    dimension = len(covariance)
    entries: list[list[CorrelationEntry]] = [
        [CorrelationEntry(None, "UNAVAILABLE") for _right in range(dimension)]
        for _left in range(dimension)
    ]
    standard_errors = tuple(
        math.sqrt(covariance[index][index]) if covariance[index][index] > 0.0 else None
        for index in range(dimension)
    )
    tolerance = policy.correlation_roundoff_multiplier * _EPSILON
    for left in range(dimension):
        if standard_errors[left] is not None and not residual_variance_degenerate:
            entries[left][left] = CorrelationEntry(1.0, "AVAILABLE", 1.0)
        else:
            entries[left][left] = CorrelationEntry(
                None,
                "RESIDUAL_VARIANCE_DEGENERACY"
                if residual_variance_degenerate
                else "INVALID_VARIANCE",
            )
        for right in range(left + 1, dimension):
            left_error = standard_errors[left]
            right_error = standard_errors[right]
            if residual_variance_degenerate:
                entry = CorrelationEntry(None, "RESIDUAL_VARIANCE_DEGENERACY")
            elif left_error is None or right_error is None:
                entry = CorrelationEntry(None, "INVALID_VARIANCE")
            else:
                raw = (covariance[left][right] / left_error) / right_error
                if not math.isfinite(raw):
                    entry = CorrelationEntry(None, "NON_FINITE_CORRELATION")
                else:
                    excess = max(0.0, abs(raw) - 1.0)
                    if abs(raw) <= 1.0:
                        entry = CorrelationEntry(
                            raw, "AVAILABLE", raw, excess, tolerance
                        )
                    elif excess <= tolerance:
                        entry = CorrelationEntry(
                            math.copysign(1.0, raw),
                            "ENDPOINT_CANONICALIZED",
                            raw,
                            excess,
                            tolerance,
                        )
                    else:
                        entry = CorrelationEntry(
                            None,
                            "CORRELATION_RANGE_VIOLATION",
                            raw,
                            excess,
                            tolerance,
                        )
            entries[left][right] = entry
            entries[right][left] = entry
    return tuple(tuple(row) for row in entries)


@dataclass(frozen=True, slots=True)
class CorrelationEvidence:
    source_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    source_family: str
    output_ids: tuple[str, ...]
    units: tuple[ParameterUnit, ...]
    entries: tuple[tuple[CorrelationEntry, ...], ...]
    claims: tuple[ClaimAssessment, ...]
    scope_reportable: bool
    source_artifact: CovarianceEvidence | ConstrainedPropagationEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_policy: UncertaintyPolicy = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        source = self.source_artifact
        source_output_ids = (
            source.output_ids
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.controlled_ids
        )
        source_units = (
            source.output_units
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.units
        )
        source_covariance = source.covariance
        residual_degenerate = (
            source.source_covariance.residual_variance_scale == 0.0
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.residual_variance_scale == 0.0
        )
        expected_entries = _expected_correlation_entries(
            source_covariance,
            residual_variance_degenerate=residual_degenerate,
            policy=self.source_policy,
        )
        if (
            self.source_identity != source.identity
            or self.accepted_result_identity != source.accepted_result_identity
            or self.accepted_occurrence_identity != source.accepted_occurrence_identity
            or self.output_ids != source_output_ids
            or self.units != source_units
            or self.entries != expected_entries
            or self.source_family
            != (
                "constrained_propagation"
                if isinstance(source, ConstrainedPropagationEvidence)
                else "local_covariance"
            )
            or self.claims[: len(source.claims)] != source.claims
        ):
            raise ValueError("Correlation evidence does not derive from its source")
        complete = all(
            item.value is not None for row in expected_entries for item in row
        )
        if _named_claim(self.claims, "CORRELATION_ENTRY_VALIDITY").state is not (
            ClaimState.SATISFIED if complete else ClaimState.VIOLATED
        ):
            raise ValueError("Correlation validity claim is inconsistent")
        source_reportable = (
            source.source_covariance.usable
            and source.claim("LOCAL_FIRST_ORDER_DEGENERACY") is ClaimState.SATISFIED
            if isinstance(source, ConstrainedPropagationEvidence)
            else source.usable
        )
        expected_reportable = source_reportable and complete
        if (
            _named_claim(self.claims, "CORRELATION_SCOPE_REPORTABILITY").state
            is not (
                ClaimState.SATISFIED if expected_reportable else ClaimState.VIOLATED
            )
            or self.scope_reportable is not expected_reportable
        ):
            raise ValueError("Correlation reportability claim is inconsistent")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-correlation-evidence",
                (
                    self.source_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_family,
                    self.output_ids,
                    self.units,
                    tuple(
                        tuple(
                            (
                                None
                                if item.value is None
                                else _float_token(item.value),
                                item.outcome,
                                None
                                if item.raw_value is None
                                else _float_token(item.raw_value),
                                None
                                if item.excess is None
                                else _float_token(item.excess),
                                None
                                if item.tolerance is None
                                else _float_token(item.tolerance),
                            )
                            for item in row
                        )
                        for row in self.entries
                    ),
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                    self.scope_reportable,
                    _CORRELATION_VERSION,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ConstraintJacobianEvidence:
    request_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    accepted_evaluation_identity: str
    problem_identity: str
    parameterization_identity: str
    constraint_program_identity: str
    capability_selection_identity: str
    controlled_ids: tuple[str, ...]
    output_ids: tuple[str, ...]
    output_units: tuple[ParameterUnit, ...]
    output_scales: tuple[float, ...]
    matrix: tuple[tuple[float, ...], ...]
    structural_dependencies: tuple[tuple[str, ...], ...]
    function_partial_diagnostics: tuple[FunctionPartialDiagnostic, ...]
    claims: tuple[ClaimAssessment, ...]
    accepted_anchor: AcceptedFitResult = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_problem: OptimizationProblem = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_parameterization: ActiveParameterization = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_policy: UncertaintyPolicy = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_capabilities: CompiledConstraintLinearizationCapabilities = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_accepted_anchor(
            self.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        if len(self.output_scales) != len(self.output_ids) or any(
            value <= 0.0 or not math.isfinite(value) for value in self.output_scales
        ):
            raise ValueError("Constraint-output scales must be positive and complete")
        matrix = _canonical_matrix(
            self.matrix,
            rows=len(self.output_ids),
            columns=len(self.controlled_ids),
            name="constraint-output Jacobian",
        )
        if (
            self.accepted_evaluation_identity
            != self.accepted_anchor.evaluation_result.identity
            or self.problem_identity != self.accepted_anchor.problem_identity
            or self.parameterization_identity
            != self.accepted_anchor.parameterization_identity
            or self.controlled_ids != self.accepted_anchor.controlled_ids
            or len(set(self.output_ids)) != len(self.output_ids)
            or len(self.output_units) != len(self.output_ids)
            or len(self.structural_dependencies) != len(self.output_ids)
            or any(
                any(param_id not in self.controlled_ids for param_id in row)
                for row in self.structural_dependencies
            )
        ):
            raise ValueError("Constraint Jacobian provenance or scope is inconsistent")
        expected_claims = (
            ClaimAssessment(
                "CONSTRAINT_OUTPUT_LINEARIZATION_COMPLETE",
                ClaimState.SATISFIED,
            ),
            ClaimAssessment(
                "CONSTRAINT_OUTPUT_LINEARIZATION_RELIABLE",
                ClaimState.SATISFIED,
            ),
        )
        if self.claims != expected_claims:
            raise ValueError("Constraint Jacobian claims are inconsistent")
        rows, dependencies, diagnostics, failure = _constraint_rows(
            self.accepted_anchor,
            self.source_problem,
            self.source_parameterization,
            self.source_policy,
            self.source_capabilities,
            self.output_ids,
            None,
        )
        if (
            failure is not None
            or rows != matrix
            or dependencies != self.structural_dependencies
            or diagnostics is None
            or tuple(item.identity for item in diagnostics)
            != tuple(item.identity for item in self.function_partial_diagnostics)
            or self.problem_identity != self.source_problem.identity
            or self.parameterization_identity != self.source_parameterization.identity
            or self.constraint_program_identity
            != self.source_parameterization.program.fingerprint
            or self.capability_selection_identity != self.source_capabilities.identity
        ):
            raise ValueError(
                "Constraint Jacobian does not match its compiled derivation"
            )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-constraint-jacobian-evidence",
                (
                    self.request_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.accepted_evaluation_identity,
                    self.problem_identity,
                    self.parameterization_identity,
                    self.constraint_program_identity,
                    self.capability_selection_identity,
                    self.controlled_ids,
                    self.output_ids,
                    self.output_units,
                    _vector_tokens(self.output_scales),
                    _matrix_tokens(matrix),
                    self.structural_dependencies,
                    tuple(item.identity for item in self.function_partial_diagnostics),
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                    _CONSTRAINT_LINEARIZATION_VERSION,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ConstrainedPropagationEvidence:
    request_identity: str
    source_covariance_identity: str
    source_constraint_jacobian_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    controlled_ids: tuple[str, ...]
    output_ids: tuple[str, ...]
    output_units: tuple[ParameterUnit, ...]
    output_scales: tuple[float, ...]
    factor: tuple[tuple[float, ...], ...]
    covariance: tuple[tuple[float, ...], ...]
    claims: tuple[ClaimAssessment, ...]
    source_covariance: CovarianceEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    source_constraint_jacobian: ConstraintJacobianEvidence = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    accepted_anchor: AcceptedFitResult = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_accepted_anchor(
            self.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        rows = len(self.output_ids)
        columns = len(self.controlled_ids)
        if len(self.output_scales) != rows or any(
            value <= 0.0 or not math.isfinite(value) for value in self.output_scales
        ):
            raise ValueError("Propagated-output scales must be positive and complete")
        factor = _canonical_matrix(
            self.factor,
            rows=rows,
            columns=columns,
            name="propagated covariance factor",
        )
        covariance = _canonical_matrix(
            self.covariance,
            rows=rows,
            columns=rows,
            name="propagated covariance",
        )
        source_covariance = self.source_covariance
        source_jacobian = self.source_constraint_jacobian
        if (
            self.source_covariance_identity != source_covariance.identity
            or self.source_constraint_jacobian_identity != source_jacobian.identity
            or self.controlled_ids != source_covariance.controlled_ids
            or self.controlled_ids != source_jacobian.controlled_ids
            or self.output_ids != source_jacobian.output_ids
            or self.output_units != source_jacobian.output_units
            or self.output_scales != source_jacobian.output_scales
        ):
            raise ValueError("Constrained propagation lineage is inconsistent")
        expected_factor = tuple(
            tuple(
                _finite(
                    _pairwise_sum(
                        tuple(
                            source_jacobian.matrix[row][inner]
                            * source_covariance.factor[inner][column]
                            for inner in range(columns)
                        )
                    ),
                    name="validated propagated factor",
                )
                for column in range(columns)
            )
            for row in range(rows)
        )
        if factor != expected_factor or covariance != _gram_matrix(expected_factor):
            raise ValueError("Constrained covariance differs from G_S L")
        inherited_claims = tuple(
            ClaimAssessment(f"SOURCE_COVARIANCE::{item.name}", item.state, item.detail)
            for item in source_covariance.claims
        ) + tuple(
            ClaimAssessment(
                f"SOURCE_CONSTRAINT_LINEARIZATION::{item.name}",
                item.state,
                item.detail,
            )
            for item in source_jacobian.claims
        )
        if self.claims[: len(inherited_claims)] != inherited_claims:
            raise ValueError(
                "Constrained propagation inherited claims are inconsistent"
            )
        per_output_degeneracy = tuple(
            all(value == 0.0 for value in source_jacobian.matrix[index])
            and bool(source_jacobian.structural_dependencies[index])
            for index in range(rows)
        )
        expected_states = {
            "EXACT_LINEAGE": ClaimState.SATISFIED,
            "CONSTRAINT_LINEARIZATION_REGULARITY": ClaimState.SATISFIED,
            "GRAM_ARITHMETIC_INTEGRITY": ClaimState.SATISFIED,
            "LOCAL_FIRST_ORDER_DEGENERACY": (
                ClaimState.VIOLATED
                if any(per_output_degeneracy)
                else ClaimState.SATISFIED
            ),
            "RESIDUAL_VARIANCE_NONDEGENERACY": (
                ClaimState.VIOLATED
                if source_covariance.residual_variance_scale == 0.0
                else ClaimState.SATISFIED
            ),
            "SOURCE_COVARIANCE_USABILITY": source_covariance.claim(
                "USABLE_LOCAL_COVARIANCE"
            ),
        }
        expected_states.update(
            {
                f"OUTPUT_FIRST_ORDER_NONDEGENERACY::{param_id}": (
                    ClaimState.VIOLATED
                    if per_output_degeneracy[index]
                    else ClaimState.SATISFIED
                )
                for index, param_id in enumerate(self.output_ids)
            }
        )
        for name, state in expected_states.items():
            if _named_claim(self.claims, name).state is not state:
                raise ValueError(
                    f"Constrained propagation claim {name!r} is inconsistent"
                )
        object.__setattr__(self, "factor", factor)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-constrained-propagation-evidence",
                (
                    self.request_identity,
                    self.source_covariance_identity,
                    self.source_constraint_jacobian_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.controlled_ids,
                    self.output_ids,
                    self.output_units,
                    _vector_tokens(self.output_scales),
                    _matrix_tokens(factor),
                    _matrix_tokens(covariance),
                    tuple(
                        (item.name, item.state.value, item.detail)
                        for item in self.claims
                    ),
                    _PROPAGATION_VERSION,
                    _REDUCTION_VERSION,
                ),
            ),
        )

    def claim(self, name: str) -> ClaimState:
        return next(item.state for item in self.claims if item.name == name)


@dataclass(frozen=True, slots=True)
class UncertaintyEvidence:
    accepted_result_identity: str
    accepted_occurrence_identity: str
    request_identity: str
    policy_identity: str
    resolved_environment_identity: str = field(compare=False)
    operations: tuple[DerivationOperation, ...]
    failures: tuple[EvidenceFailure, ...]
    residual_jacobian: ResidualJacobianEvidence | None = None
    rank_diagnostic: RankDiagnostic | None = None
    covariance: CovarianceEvidence | None = None
    marginal_errors: MarginalErrorEvidence | None = None
    correlations: CorrelationEvidence | None = None
    constraint_jacobian: ConstraintJacobianEvidence | None = None
    constrained_propagation: ConstrainedPropagationEvidence | None = None
    constrained_marginal_errors: MarginalErrorEvidence | None = None
    constrained_correlations: CorrelationEvidence | None = None
    accepted_anchor: AcceptedFitResult = field(
        repr=False,
        compare=False,
        metadata={"record": False},
        kw_only=True,
    )
    requested_output_scope: tuple[str, ...] = field(kw_only=True)
    requested_output_units: tuple[ParameterUnit, ...] = field(kw_only=True)
    requested_output_scales: tuple[float, ...] = field(kw_only=True)
    source_problem: OptimizationProblem = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_parameterization: ActiveParameterization = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_engine: EvaluationEngine = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_policy: UncertaintyPolicy = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    source_capabilities: CompiledConstraintLinearizationCapabilities = field(
        repr=False, compare=False, metadata={"record": False}, kw_only=True
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - complete bundle integrity chain
        _validate_accepted_anchor(
            self.accepted_anchor,
            self.accepted_result_identity,
            self.accepted_occurrence_identity,
        )
        expected_request = _identity(
            "native-uncertainty-request",
            (
                self.accepted_anchor.identity,
                self.accepted_anchor.occurrence_identity,
                self.source_problem.identity,
                self.source_parameterization.identity,
                self.source_engine.plan.identity,
                self.source_policy.identity,
                self.requested_output_scope,
                self.requested_output_units,
                _vector_tokens(self.requested_output_scales),
                self.source_capabilities.identity,
            ),
        )
        if (
            self.request_identity != expected_request
            or self.policy_identity != self.source_policy.identity
            or self.source_capabilities.output_scope != self.requested_output_scope
            or len(self.requested_output_units) != len(self.requested_output_scope)
            or len(self.requested_output_scales) != len(self.requested_output_scope)
        ):
            raise ValueError(
                "Uncertainty request does not match its exact source inputs"
            )
        artifacts = tuple(
            item
            for item in (
                self.residual_jacobian,
                self.rank_diagnostic,
                self.covariance,
                self.marginal_errors,
                self.correlations,
                self.constraint_jacobian,
                self.constrained_propagation,
                self.constrained_marginal_errors,
                self.constrained_correlations,
            )
            if item is not None
        )
        if any(
            item.accepted_result_identity != self.accepted_result_identity
            or item.accepted_occurrence_identity != self.accepted_occurrence_identity
            for item in artifacts
        ):
            raise ValueError("Evidence bundle mixes accepted fit occurrences")
        if self.residual_jacobian is not None and (
            self.rank_diagnostic is not None
            and self.rank_diagnostic.source_jacobian_identity
            != self.residual_jacobian.identity
            or self.covariance is not None
            and self.covariance.source_jacobian_identity
            != self.residual_jacobian.identity
        ):
            raise ValueError("Evidence bundle mixes residual Jacobian sources")
        if self.covariance is not None and (
            self.rank_diagnostic is None
            or self.covariance.rank_diagnostic_identity != self.rank_diagnostic.identity
            or self.marginal_errors is not None
            and self.marginal_errors.source_identity != self.covariance.identity
            or self.correlations is not None
            and self.correlations.source_identity != self.covariance.identity
        ):
            raise ValueError("Evidence bundle has incompatible covariance descendants")
        if self.constrained_propagation is not None and (
            self.covariance is None
            or self.constraint_jacobian is None
            or self.constrained_propagation.source_covariance_identity
            != self.covariance.identity
            or self.constrained_propagation.source_constraint_jacobian_identity
            != self.constraint_jacobian.identity
            or self.constrained_marginal_errors is not None
            and self.constrained_marginal_errors.source_identity
            != self.constrained_propagation.identity
            or self.constrained_correlations is not None
            and self.constrained_correlations.source_identity
            != self.constrained_propagation.identity
        ):
            raise ValueError("Evidence bundle has incompatible constrained descendants")
        stage_artifacts = {
            stage: artifact.identity
            for stage, artifact in (
                ("residual_linearization", self.residual_jacobian),
                ("covariance", self.covariance),
                ("marginal_errors", self.marginal_errors),
                ("correlations", self.correlations),
                ("constraint_linearization", self.constraint_jacobian),
                ("constrained_propagation", self.constrained_propagation),
                ("constrained_marginal_errors", self.constrained_marginal_errors),
                ("constrained_correlations", self.constrained_correlations),
            )
            if artifact is not None
        }
        completed = {
            item.stage: item.artifact_identity
            for item in self.operations
            if item.terminal is OperationTerminal.COMPLETED
        }
        if completed != stage_artifacts:
            raise ValueError("Evidence operations do not bind the assembled artifacts")
        if self.residual_jacobian is not None:
            expected_residual_request = _identity(
                "native-residual-linearization-request",
                (
                    self.request_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_engine.plan.identity,
                    self.source_policy.identity,
                ),
            )
            if self.residual_jacobian.request_identity != expected_residual_request:
                raise ValueError("Residual Jacobian request lineage is inconsistent")
        if self.rank_diagnostic is not None:
            expected_covariance_request = _identity(
                "native-covariance-request",
                (
                    self.request_identity,
                    None
                    if self.residual_jacobian is None
                    else self.residual_jacobian.identity,
                    self.source_policy.identity,
                ),
            )
            if self.rank_diagnostic.request_identity != expected_covariance_request or (
                self.covariance is not None
                and self.covariance.request_identity != expected_covariance_request
            ):
                raise ValueError("Covariance request lineage is inconsistent")
        if self.constraint_jacobian is not None:
            expected_constraint_request = _identity(
                "native-constraint-linearization-request",
                (
                    self.request_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.requested_output_scope,
                    self.requested_output_units,
                    _vector_tokens(self.requested_output_scales),
                    self.source_policy.identity,
                ),
            )
            if (
                self.constraint_jacobian.request_identity != expected_constraint_request
                or self.constraint_jacobian.output_ids != self.requested_output_scope
                or self.constraint_jacobian.output_units != self.requested_output_units
                or self.constraint_jacobian.output_scales
                != self.requested_output_scales
            ):
                raise ValueError("Constraint Jacobian request lineage is inconsistent")
        if self.constrained_propagation is not None:
            expected_propagation_request = _identity(
                "native-constrained-propagation-request",
                (
                    self.request_identity,
                    None if self.covariance is None else self.covariance.identity,
                    None
                    if self.constraint_jacobian is None
                    else self.constraint_jacobian.identity,
                    self.requested_output_scope,
                    self.source_policy.identity,
                ),
            )
            if (
                self.constrained_propagation.request_identity
                != expected_propagation_request
            ):
                raise ValueError(
                    "Constrained propagation request lineage is inconsistent"
                )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-uncertainty-evidence-bundle",
                (
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.request_identity,
                    self.policy_identity,
                    tuple(item.identity for item in self.operations),
                    tuple(item.identity for item in self.failures),
                    None
                    if self.residual_jacobian is None
                    else self.residual_jacobian.identity,
                    None
                    if self.rank_diagnostic is None
                    else self.rank_diagnostic.identity,
                    None if self.covariance is None else self.covariance.identity,
                    None
                    if self.marginal_errors is None
                    else self.marginal_errors.identity,
                    None if self.correlations is None else self.correlations.identity,
                    None
                    if self.constraint_jacobian is None
                    else self.constraint_jacobian.identity,
                    None
                    if self.constrained_propagation is None
                    else self.constrained_propagation.identity,
                    None
                    if self.constrained_marginal_errors is None
                    else self.constrained_marginal_errors.identity,
                    None
                    if self.constrained_correlations is None
                    else self.constrained_correlations.identity,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        """Serialize the closed typed bundle without runtime machinery."""
        return {
            "artifact_type": "native_uncertainty_evidence_bundle",
            "schema_version": _SCHEMA_VERSION,
            "payload": _record_value(self),
        }


@dataclass(frozen=True, slots=True)
class _DerivativeValue:
    value: float
    gradient: tuple[float, ...]
    structural_dependencies: frozenset[str]
    function_partial_diagnostics: tuple[FunctionPartialDiagnostic, ...] = ()


class _IdentityBearing(Protocol):
    @property
    def identity(self) -> str: ...


def _pairwise_sum(values: Sequence[float]) -> float:
    terms = [0.0 if value == 0.0 else value for value in values]
    while len(terms) > 1:
        reduced: list[float] = []
        for index in range(0, len(terms) - 1, 2):
            partial = terms[index] + terms[index + 1]
            if not math.isfinite(partial):
                raise ArithmeticError("fixed pairwise reduction became non-finite")
            reduced.append(0.0 if partial == 0.0 else partial)
        if len(terms) % 2:
            reduced.append(terms[-1])
        terms = reduced
    return terms[0] if terms else 0.0


def _pairwise_sum_feasibility(values: Sequence[float]) -> float:
    """Use the canonical reduction order while retaining exceptional values."""
    terms = list(values)
    while len(terms) > 1:
        reduced = [
            terms[index] + terms[index + 1] for index in range(0, len(terms) - 1, 2)
        ]
        if len(terms) % 2:
            reduced.append(terms[-1])
        terms = reduced
    return terms[0] if terms else 0.0


def _gram_matrix(factor: Sequence[Sequence[float]]) -> tuple[tuple[float, ...], ...]:
    rows = len(factor)
    result = [[0.0] * rows for _index in range(rows)]
    for left in range(rows):
        for right in range(left, rows):
            value = _pairwise_sum(
                tuple(
                    factor[left][column] * factor[right][column]
                    for column in range(len(factor[left]))
                )
            )
            value = _finite(value, name="Gram matrix entry")
            result[left][right] = value
            result[right][left] = value
    return tuple(tuple(row) for row in result)


def _canonical_covariance_reduction(
    jacobian: Sequence[Sequence[float]],
    coordinate_scales: Sequence[float],
    residual_variance_scale: float,
) -> tuple[
    tuple[tuple[float, ...], ...],
    tuple[tuple[float, ...], ...],
    tuple[tuple[float, ...], ...],
]:
    """Reduce the full-rank information matrix through one canonical factor."""
    matrix = np.asarray(jacobian, dtype=np.float64)
    scales = np.asarray(coordinate_scales, dtype=np.float64)
    scaled = matrix * scales[np.newaxis, :]
    information = scaled.T @ scaled
    inverse_scaled = np.linalg.solve(
        information,
        np.eye(information.shape[0], dtype=np.float64),
    )
    scale_matrix = np.diag(scales)
    unscaled_array = scale_matrix @ inverse_scaled @ scale_matrix
    unscaled_array = 0.5 * (unscaled_array + unscaled_array.T)
    unscaled = _canonical_matrix(
        unscaled_array,
        rows=len(scales),
        columns=len(scales),
        name="canonical unscaled covariance",
    )
    unscaled_factor = np.linalg.cholesky(unscaled_array)
    factor_array = math.sqrt(residual_variance_scale) * unscaled_factor
    factor = _canonical_matrix(
        factor_array,
        rows=len(scales),
        columns=len(scales),
        name="canonical covariance factor",
    )
    return unscaled, factor, _gram_matrix(factor)


def _has_positive_scale_zero_variance(
    covariance: Sequence[Sequence[float]],
    residual_variance_scale: float,
) -> bool:
    return residual_variance_scale > 0.0 and any(
        row[index] == 0.0 for index, row in enumerate(covariance)
    )


def _invariant_singular_subspaces(
    right: Array,
    singular_values: tuple[float, ...],
    *,
    rank_threshold: float,
    weak_threshold: float,
    cluster_relative_tolerance: float,
) -> tuple[SingularSubspaceEvidence, ...]:
    """Retain projectors, never basis orientation, for clustered directions."""

    groups = _singular_cluster_indices(
        singular_values,
        rank_threshold=rank_threshold,
        weak_threshold=weak_threshold,
        cluster_relative_tolerance=cluster_relative_tolerance,
    )

    def projector(indices: tuple[int, ...]) -> tuple[tuple[float, ...], ...]:
        factor = tuple(
            tuple(
                _finite(right[row, column], name="right singular subspace")
                for column in indices
            )
            for row in range(right.shape[0])
        )
        return _gram_matrix(factor)

    return tuple(
        SingularSubspaceEvidence(
            indices,
            tuple(singular_values[index] for index in indices),
            ("isolated_" if len(indices) == 1 else "clustered_")
            + _singular_direction_class(
                singular_values,
                indices[0],
                rank_threshold=rank_threshold,
                weak_threshold=weak_threshold,
            ),
            projector(indices),
        )
        for indices in groups
    )


def _singular_cluster_indices(
    singular_values: tuple[float, ...],
    *,
    rank_threshold: float,
    weak_threshold: float,
    cluster_relative_tolerance: float,
) -> tuple[tuple[int, ...], ...]:
    largest = singular_values[0] if singular_values else 0.0

    groups: list[tuple[int, ...]] = []
    current: list[int] = []
    for index, value in enumerate(singular_values):
        if current:
            previous = singular_values[current[-1]]
            if not (
                _singular_direction_class(
                    singular_values,
                    index,
                    rank_threshold=rank_threshold,
                    weak_threshold=weak_threshold,
                )
                == _singular_direction_class(
                    singular_values,
                    current[-1],
                    rank_threshold=rank_threshold,
                    weak_threshold=weak_threshold,
                )
                and abs(previous - value)
                <= cluster_relative_tolerance * max(largest, previous, value)
            ):
                groups.append(tuple(current))
                current = []
        current.append(index)
    if current:
        groups.append(tuple(current))
    return tuple(groups)


def _singular_direction_class(
    singular_values: tuple[float, ...],
    index: int,
    *,
    rank_threshold: float,
    weak_threshold: float,
) -> Literal["identifiable", "weak", "null"]:
    value = singular_values[index]
    if value <= rank_threshold:
        return "null"
    if value <= weak_threshold:
        return "weak"
    return "identifiable"


def _singular_spectrum_error(
    values: Sequence[float],
    expected_count: int,
) -> tuple[str, str] | None:
    if len(values) != expected_count:
        return (
            "incomplete_singular_spectrum",
            "Economical SVD did not produce one value per coordinate",
        )
    if any(value < 0.0 for value in values) or any(
        left < right for left, right in pairwise(values)
    ):
        return (
            "invalid_singular_spectrum",
            "SVD singular spectrum must be non-negative and descending",
        )
    return None


def _max_norm(values: Sequence[float]) -> float:
    return max((abs(value) for value in values), default=0.0)


def _three_point_derivative(
    first: Sequence[float],
    center: Sequence[float],
    second: Sequence[float],
    first_displacement: float,
    second_displacement: float,
) -> tuple[float, ...]:
    first_coefficient = -second_displacement / (
        first_displacement * (first_displacement - second_displacement)
    )
    center_coefficient = -(first_displacement + second_displacement) / (
        first_displacement * second_displacement
    )
    second_coefficient = -first_displacement / (
        second_displacement * (second_displacement - first_displacement)
    )
    result = tuple(
        _finite(
            _pairwise_sum(
                (
                    first_coefficient * left,
                    center_coefficient * middle,
                    second_coefficient * right,
                )
            ),
            name="finite-difference derivative",
        )
        for left, middle, right in zip(first, center, second, strict=True)
    )
    return result


def _step_candidates(
    nominal: float,
    maximum: float,
    policy: UncertaintyPolicy,
) -> tuple[float, ...]:
    candidates: list[tuple[float, float]] = []
    for exponent in range(-policy.smaller_step_extent, policy.larger_step_extent + 1):
        candidate = nominal * math.ldexp(1.0, exponent)
        if math.isfinite(candidate) and 0.0 < candidate <= maximum:
            candidates.append((abs(float(exponent)), candidate))
    if math.isfinite(maximum) and maximum > 0.0:
        distance = abs(math.log2(maximum / nominal)) if nominal > 0.0 else math.inf
        candidates.append((distance, maximum))
    unique: dict[str, tuple[float, float]] = {}
    for distance, candidate in candidates:
        token = _float_token(candidate)
        current = unique.get(token)
        if current is None or (distance, candidate) < current:
            unique[token] = (distance, candidate)
    return tuple(
        candidate
        for _distance, candidate in sorted(unique.values(), key=lambda item: item)
    )


def _lineage_failure(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    policy: UncertaintyPolicy,
) -> str | None:
    if not accepted_occurrence_is_authoritative(accepted):
        return "accepted_occurrence_identity_mismatch"
    scales = tuple(param_id for param_id, _value in policy.coordinate_scales)
    if scales != problem.controlled_ids:
        return "coordinate_scale_scope_mismatch"
    if (
        accepted.problem_identity != problem.identity
        or accepted.parameterization_identity != problem.parameterization_identity
        or accepted.evaluator_parameterization_identity
        != problem.evaluator_parameterization_identity
        or accepted.controlled_ids != problem.controlled_ids
        or accepted.source_occurrence_identity
        != problem.source_snapshot.occurrence_identity
        or accepted.source_revision != problem.source_snapshot.revision
        or accepted.evaluation_result.plan_identity != problem.evaluation_plan_identity
    ):
        return "accepted_result_lineage_mismatch"
    try:
        problem.validate_parameterization(parameterization)
    except ValueError:
        return "parameterization_lineage_mismatch"
    if (
        engine.plan.identity != problem.evaluation_plan_identity
        or engine.plan.parameterization_identity
        != problem.evaluator_parameterization_identity
        or engine.plan.constraint_program_identity
        != problem.constraint_program_identity
    ):
        return "evaluation_engine_lineage_mismatch"
    return None


def _evaluate_vector(
    vector: tuple[float, ...],
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    evaluator: BoundEvaluator,
) -> EvaluationResult | EvaluationFailure:
    lifecycle = problem.lifecycle_frame(vector, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    return evaluator.evaluate(frame)


def _column_orientations(
    value: float,
    lower: float,
    upper: float,
) -> tuple[tuple[str, float], ...]:
    negative_distance = value - lower
    positive_distance = upper - value
    central_distance = min(negative_distance, positive_distance) / 2.0
    candidates: list[tuple[str, float]] = []
    if central_distance > 0.0:
        candidates.append(("centered", central_distance))
    negative_maximum = negative_distance / 4.0
    positive_maximum = positive_distance / 4.0
    one_sided = (
        (
            ("one_sided_positive", positive_maximum),
            ("one_sided_negative", negative_maximum),
        )
        if positive_maximum >= negative_maximum
        else (
            ("one_sided_negative", negative_maximum),
            ("one_sided_positive", positive_maximum),
        )
    )
    candidates.extend(item for item in one_sided if item[1] > 0.0)
    return tuple(candidates)


def _stencil_vectors(
    vector: tuple[float, ...],
    column: int,
    step: float,
    orientation: str,
) -> tuple[tuple[float, ...], ...] | None:
    center = vector[column]
    multipliers = (
        (-2.0, -1.0, 1.0, 2.0)
        if orientation == "centered"
        else (
            (1.0, 2.0, 4.0)
            if orientation == "one_sided_positive"
            else (-1.0, -2.0, -4.0)
        )
    )
    vectors: list[tuple[float, ...]] = []
    represented: set[float] = {center}
    for multiplier in multipliers:
        updated = center + multiplier * step
        if not math.isfinite(updated) or updated in represented:
            return None
        represented.add(updated)
        values = list(vector)
        values[column] = updated
        vectors.append(tuple(values))
    return tuple(vectors)


def _evaluate_stencil(
    vectors: tuple[tuple[float, ...], ...],
    *,
    column: int,
    center: float,
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    evaluator: BoundEvaluator,
) -> tuple[
    list[tuple[float, ...]] | None,
    list[float] | None,
    EvidenceFailure | None,
]:
    results: list[tuple[float, ...]] = []
    displacements: list[float] = []
    for trial in vectors:
        result = _evaluate_vector(
            trial,
            problem=problem,
            parameterization=parameterization,
            evaluator=evaluator,
        )
        if isinstance(result, EvaluationFailure):
            if result.validity == "INVALID_TRIAL":
                return None, None, None
            return (
                None,
                None,
                EvidenceFailure(
                    "residual_linearization",
                    "fatal_stencil_evaluation_failure",
                    result.category,
                    accepted.identity,
                ),
            )
        results.append(tuple(float(item) for item in result.residuals))
        displacements.append(trial[column] - center)
    return results, displacements, None


def _stencil_estimates(
    orientation: str,
    results: Sequence[Sequence[float]],
    base_residuals: tuple[float, ...],
    displacements: Sequence[float],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if orientation == "centered":
        coarse = _three_point_derivative(
            results[0],
            base_residuals,
            results[3],
            displacements[0],
            displacements[3],
        )
        fine = _three_point_derivative(
            results[1],
            base_residuals,
            results[2],
            displacements[1],
            displacements[2],
        )
        return fine, coarse
    fine = _three_point_derivative(
        results[0],
        base_residuals,
        results[1],
        displacements[0],
        displacements[1],
    )
    coarse = _three_point_derivative(
        results[1],
        base_residuals,
        results[2],
        displacements[1],
        displacements[2],
    )
    return fine, coarse


def _assess_stencil(
    fine: tuple[float, ...],
    coarse: tuple[float, ...],
    results: Sequence[Sequence[float]],
    base_residuals: tuple[float, ...],
    displacements: Sequence[float],
    policy: UncertaintyPolicy,
    output_scale: float,
) -> tuple[float, float, float, bool]:
    discrepancy = _max_norm(
        tuple(left - right for left, right in zip(fine, coarse, strict=True))
    )
    derivative_scale = max(_max_norm(fine), _max_norm(coarse))
    function_scale = max(
        output_scale,
        _max_norm(base_residuals),
        *(_max_norm(result) for result in results),
    )
    minimum_displacement = min(abs(item) for item in displacements)
    roundoff = (
        policy.roundoff_multiplier * _EPSILON * function_scale / minimum_displacement
    )
    threshold = policy.relative_step_tolerance * derivative_scale + roundoff
    return discrepancy, derivative_scale, roundoff, discrepancy <= threshold


def _linearize_residual_column(
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    evaluator: BoundEvaluator,
    policy: UncertaintyPolicy,
    base_residuals: tuple[float, ...],
    column: int,
    param_id: str,
    scale: float,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> tuple[LinearizationColumn | None, tuple[object, ...], EvidenceFailure | None]:
    center = accepted.vector[column]
    orientations = _column_orientations(
        center,
        problem.lower_bounds[column],
        problem.upper_bounds[column],
    )
    nominal = _NOMINAL_STEP_FACTOR * scale
    trajectory: list[object] = []
    if not orientations:
        return (
            None,
            tuple(trajectory),
            EvidenceFailure(
                "residual_linearization",
                "exhausted_exact_feasible_distance",
                f"No nonzero feasible displacement exists for {param_id}",
                accepted.identity,
                _identity("residual-column-trajectory", trajectory),
                ("exhausted_exact_feasible_distance",),
            ),
        )
    attempt = 0
    for orientation, maximum in orientations:
        orientation_had_feasible_stencil = False
        for step in _step_candidates(nominal, maximum, policy):
            _raise_if_terminated(cancellation_probe)
            attempt += 1
            vectors = _stencil_vectors(accepted.vector, column, step, orientation)
            if vectors is None:
                trajectory.append((orientation, _float_token(step), "not_distinct"))
                continue
            results, displacements, failure = _evaluate_stencil(
                vectors,
                column=column,
                center=center,
                accepted=accepted,
                problem=problem,
                parameterization=parameterization,
                evaluator=evaluator,
            )
            if failure is not None:
                trajectory.append(
                    (
                        orientation,
                        _float_token(step),
                        "fatal_evaluation_failure",
                        failure.identity,
                    )
                )
                return (
                    None,
                    tuple(trajectory),
                    EvidenceFailure(
                        failure.stage,
                        failure.category,
                        failure.message,
                        failure.source_identity,
                        _identity("residual-column-trajectory", trajectory),
                        ("fatal_stencil_evaluation_termination",),
                    ),
                )
            if results is None or displacements is None:
                trajectory.append((orientation, _float_token(step), "invalid_trial"))
                continue
            orientation_had_feasible_stencil = True
            fine, coarse = _stencil_estimates(
                orientation,
                results,
                base_residuals,
                displacements,
            )
            discrepancy, derivative_scale, roundoff, reliable = _assess_stencil(
                fine,
                coarse,
                results,
                base_residuals,
                displacements,
                policy,
                1.0,
            )
            trajectory.append(
                (
                    orientation,
                    _float_token(step),
                    tuple(_float_token(item) for item in displacements),
                    _float_token(discrepancy),
                    _float_token(derivative_scale),
                    _float_token(roundoff),
                    reliable,
                )
            )
            if reliable:
                fingerprint = _identity("residual-column-trajectory", trajectory)
                return (
                    LinearizationColumn(
                        param_id,
                        orientation,
                        nominal,
                        tuple(displacements),
                        fine,
                        coarse,
                        discrepancy,
                        derivative_scale,
                        roundoff,
                        attempt,
                        "successful_reliable_estimate",
                        fingerprint,
                    ),
                    tuple(trajectory),
                    None,
                )
        if orientation == "centered" and orientation_had_feasible_stencil:
            break
    category = (
        "loss_of_distinct_binary64_representation"
        if trajectory
        and all(
            len(item) == 3 and item[2] == "not_distinct"
            for item in trajectory
            if isinstance(item, tuple)
        )
        else "exhausted_declared_step_search_extent"
    )
    return (
        None,
        tuple(trajectory),
        EvidenceFailure(
            "residual_linearization",
            category,
            f"No reliable feasible stencil for {param_id}",
            accepted.identity,
            _identity("residual-column-trajectory", trajectory),
            (
                "exhausted_declared_smaller_step_extent",
                "exhausted_declared_larger_step_extent",
            ),
        ),
    )


def _linearize_residuals(
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    policy: UncertaintyPolicy,
    request_identity: str,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> tuple[ResidualJacobianEvidence | None, EvidenceFailure | None]:
    _raise_if_terminated(cancellation_probe)
    evaluator = engine.new_evaluator()
    base = _evaluate_vector(
        accepted.vector,
        problem=problem,
        parameterization=parameterization,
        evaluator=evaluator,
    )
    if isinstance(base, EvaluationFailure):
        return None, EvidenceFailure(
            "residual_linearization",
            "accepted_point_evaluation_failed",
            base.category,
            accepted.identity,
        )
    if base.identity != accepted.evaluation_result.identity:
        return None, EvidenceFailure(
            "residual_linearization",
            "fresh_accepted_point_mismatch",
            "Fresh uncertainty baseline differs from accepted materialization",
            accepted.identity,
        )
    _raise_if_terminated(cancellation_probe)
    base_residuals = tuple(float(value) for value in base.residuals)
    scales = tuple(value for _param_id, value in policy.coordinate_scales)
    column_values: list[tuple[float, ...]] = []
    column_evidence: list[LinearizationColumn] = []
    for column, (param_id, scale) in enumerate(policy.coordinate_scales):
        _raise_if_terminated(cancellation_probe)
        accepted_column, _trajectory, failure = _linearize_residual_column(
            accepted,
            problem=problem,
            parameterization=parameterization,
            evaluator=evaluator,
            policy=policy,
            base_residuals=base_residuals,
            column=column,
            param_id=param_id,
            scale=scale,
            cancellation_probe=cancellation_probe,
        )
        if failure is not None:
            return None, failure
        if accepted_column is None:
            return None, EvidenceFailure(
                "residual_linearization",
                "reliable_stencil_unavailable",
                f"No reliable feasible stencil for {param_id}",
                accepted.identity,
            )
        column_evidence.append(accepted_column)
        column_values.append(accepted_column.fine_estimate)
    matrix = tuple(
        tuple(column[row] for column in column_values)
        for row in range(len(base_residuals))
    )
    trajectory_fingerprint = _identity(
        "residual-linearization-trajectory",
        tuple((item.param_id, item.trajectory_fingerprint) for item in column_evidence),
    )
    return (
        ResidualJacobianEvidence(
            request_identity,
            accepted.identity,
            accepted.occurrence_identity,
            accepted.evaluation_result.identity,
            problem.identity,
            problem.source_snapshot.occurrence_identity,
            problem.source_snapshot.revision,
            parameterization.identity,
            parameterization.evaluator_identity,
            parameterization.program.fingerprint,
            engine.plan.identity,
            policy.identity,
            policy.calibration_identity,
            policy.numerical_compatibility_requirement,
            problem.controlled_ids,
            accepted.vector,
            scales,
            matrix,
            tuple(column_evidence),
            trajectory_fingerprint,
            accepted_anchor=accepted,
            source_problem=problem,
            source_parameterization=parameterization,
            source_engine=engine,
            source_policy=policy,
        ),
        None,
    )


def _profiled_normalization_regular(
    accepted: AcceptedFitResult,
    engine: EvaluationEngine,
) -> bool:
    for descriptor, profile in zip(
        engine.plan.profiles,
        accepted.evaluation_result.profiles,
        strict=True,
    ):
        if not descriptor.is_scaled:
            continue
        retained = descriptor.retained_observation_indices
        denominator = _pairwise_sum(
            tuple(
                (
                    float(
                        accepted.evaluation_result.unscaled_calculations[
                            descriptor.observation_offset + index
                        ]
                    )
                    / descriptor.uncertainties[index]
                )
                ** 2
                for index in retained
            )
        )
        if denominator <= 0.0 or not math.isfinite(profile.normalization_factor):
            return False
    return True


def _box_boundary_claims(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_scale: float,
) -> tuple[ClaimAssessment, ClaimAssessment]:
    strict_interior = True
    separation = ClaimState.SATISFIED
    details: list[str] = []
    for index, (value, lower, upper) in enumerate(
        zip(
            accepted.vector,
            problem.lower_bounds,
            problem.upper_bounds,
            strict=True,
        )
    ):
        if not lower <= value <= upper:
            return (
                ClaimAssessment("INTERIOR_POINT", ClaimState.VIOLATED),
                ClaimAssessment("BOUNDARY_SEPARATION", ClaimState.VIOLATED),
            )
        if value in (lower, upper):
            strict_interior = False
        for label, slack in (
            ("lower", value - lower),
            ("upper", upper - value),
        ):
            if math.isinf(slack):
                continue
            variance = covariance[index][index]
            if residual_variance_scale == 0.0 or variance == 0.0:
                separation = ClaimState.INDETERMINATE
                details.append(f"{label}[{index}]: zero scaled variance")
                continue
            if variance < 0.0 or not math.isfinite(variance):
                separation = ClaimState.INDETERMINATE
                details.append(f"{label}[{index}]: invalid directional variance")
                continue
            zeta = float(slack / math.sqrt(variance))
            details.append(f"{label}[{index}] zeta={zeta.hex()}")
            if not math.isfinite(zeta):
                separation = ClaimState.INDETERMINATE
            elif zeta <= _BOUNDARY_THRESHOLD:
                separation = ClaimState.VIOLATED
    return (
        ClaimAssessment(
            "INTERIOR_POINT",
            ClaimState.SATISFIED if strict_interior else ClaimState.VIOLATED,
        ),
        ClaimAssessment("BOUNDARY_SEPARATION", separation, "; ".join(details)),
    )


def _classify_affine_restriction(
    label: str,
    coefficient_array: Array,
    upper_bound: float,
    full_values: Array,
    controlled_frame_indices: tuple[int, ...],
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_scale: float,
) -> tuple[ClaimState, str]:
    linear_value = _pairwise_sum_feasibility(
        tuple(
            float(coefficient) * float(value)
            for coefficient, value in zip(
                coefficient_array,
                full_values,
                strict=True,
            )
        )
    )
    slack = float(upper_bound - linear_value)
    if not math.isfinite(slack):
        return ClaimState.INDETERMINATE, f"{label}: non-finite slack"
    controlled = tuple(
        float(coefficient_array[index]) for index in controlled_frame_indices
    )
    if all(value == 0.0 for value in controlled):
        state = ClaimState.SATISFIED if slack >= 0.0 else ClaimState.VIOLATED
        return state, f"{label}: held-only slack={slack.hex()}"
    if slack < 0.0:
        return ClaimState.VIOLATED, f"{label}: negative slack={slack.hex()}"
    if residual_variance_scale == 0.0:
        return ClaimState.INDETERMINATE, f"{label}: phi=0 covariance degeneracy"
    variance = _pairwise_sum_feasibility(
        tuple(
            left_value * covariance[left][right] * right_value
            for left, left_value in enumerate(controlled)
            for right, right_value in enumerate(controlled)
        )
    )
    if not math.isfinite(variance) or variance < 0.0:
        return ClaimState.INDETERMINATE, f"{label}: invalid directional variance"
    if variance == 0.0:
        return ClaimState.INDETERMINATE, f"{label}: zero directional variance"
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        zeta = float(slack / math.sqrt(variance))
    if not math.isfinite(zeta):
        return ClaimState.INDETERMINATE, f"{label}: non-finite zeta"
    state = ClaimState.SATISFIED if zeta > _BOUNDARY_THRESHOLD else ClaimState.VIOLATED
    return state, f"{label}: zeta={zeta.hex()}"


def _affine_feasibility_claim(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_scale: float,
    affine_feasibility_policy: str,
    *,
    controlled_only: bool = False,
) -> ClaimAssessment:
    claim_name = (
        "CONTROLLED_AFFINE_SEPARATION" if controlled_only else "AFFINE_FEASIBILITY"
    )
    if not problem.affine_half_spaces and not problem.affine_equalities:
        return ClaimAssessment(
            claim_name,
            ClaimState.NOT_APPLICABLE,
            affine_feasibility_policy,
        )
    controlled_ids = set(problem.controlled_ids)
    accepted_controlled = dict(
        zip(problem.controlled_ids, accepted.vector, strict=True)
    )
    full_values = np.asarray(
        tuple(
            accepted_controlled.get(param_id, value)
            for param_id, value in problem.independent_items
        ),
        dtype=np.float64,
    )
    controlled_frame_indices = tuple(
        index
        for index, (param_id, _value) in enumerate(problem.independent_items)
        if param_id in controlled_ids
    )
    restrictions = [
        (
            item.restriction_id,
            np.asarray(item.coefficients, dtype=np.float64),
            item.upper_bound,
        )
        for item in problem.affine_half_spaces
    ]
    restrictions.extend(
        (
            f"{item.restriction_id}:positive",
            np.asarray(item.coefficients, dtype=np.float64),
            item.value,
        )
        for item in problem.affine_equalities
    )
    if controlled_only:
        restrictions = [
            item
            for item in restrictions
            if any(item[1][index] != 0.0 for index in controlled_frame_indices)
        ]
    if not restrictions:
        return ClaimAssessment(
            claim_name,
            ClaimState.NOT_APPLICABLE,
            "No affine restriction has a controlled-coordinate coefficient",
        )
    restrictions.extend(
        (
            f"{item.restriction_id}:negative",
            -np.asarray(item.coefficients, dtype=np.float64),
            -item.value,
        )
        for item in problem.affine_equalities
    )
    assessments = tuple(
        _classify_affine_restriction(
            label,
            coefficients,
            upper_bound,
            full_values,
            controlled_frame_indices,
            covariance,
            residual_variance_scale,
        )
        for label, coefficients, upper_bound in restrictions
    )
    states = tuple(state for state, _detail in assessments)
    if ClaimState.VIOLATED in states:
        affine_state = ClaimState.VIOLATED
    elif ClaimState.INDETERMINATE in states:
        affine_state = ClaimState.INDETERMINATE
    else:
        affine_state = ClaimState.SATISFIED
    return ClaimAssessment(
        claim_name,
        affine_state,
        "; ".join(detail for _state, detail in assessments),
    )


def _boundary_claims(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_scale: float,
    affine_feasibility_policy: str,
) -> tuple[
    ClaimAssessment,
    ClaimAssessment,
    ClaimAssessment,
    ClaimAssessment,
]:
    box = _box_boundary_claims(
        accepted,
        problem,
        covariance,
        residual_variance_scale,
    )
    if box[0].state is ClaimState.VIOLATED and any(
        not lower <= value <= upper
        for value, lower, upper in zip(
            accepted.vector,
            problem.lower_bounds,
            problem.upper_bounds,
            strict=True,
        )
    ):
        affine = ClaimAssessment(
            "AFFINE_FEASIBILITY",
            ClaimState.VIOLATED,
            "Accepted point is outside its box",
        )
    else:
        affine = _affine_feasibility_claim(
            accepted,
            problem,
            covariance,
            residual_variance_scale,
            affine_feasibility_policy,
        )
    controlled_affine = _affine_feasibility_claim(
        accepted,
        problem,
        covariance,
        residual_variance_scale,
        affine_feasibility_policy,
        controlled_only=True,
    )
    return (*box, affine, controlled_affine)


def _residual_variance_scale(
    accepted: AcceptedFitResult,
    *,
    degrees_of_freedom: int,
    policy: UncertaintyPolicy,
) -> tuple[float | None, EvidenceFailure | None]:
    if not math.isfinite(accepted.chi_square):
        return None, EvidenceFailure(
            "covariance",
            "non_finite_chi_square",
            "Accepted chi-square is non-finite",
            accepted.identity,
        )
    if accepted.chi_square < 0.0:
        return None, EvidenceFailure(
            "covariance",
            "negative_chi_square",
            "Accepted chi-square is negative",
            accepted.identity,
        )
    if (
        policy.residual_variance_scaling
        is ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
    ):
        return 1.0, None
    scale = accepted.chi_square / degrees_of_freedom
    if not math.isfinite(scale) or scale < 0.0:
        return None, EvidenceFailure(
            "covariance",
            "invalid_residual_variance_scale",
            "Estimated residual variance is negative or non-finite",
            accepted.identity,
        )
    return float(scale), None


def _full_dimensional_feasible_interior_claim(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    box_interior: ClaimAssessment,
) -> ClaimAssessment:
    controlled_ids = set(problem.controlled_ids)
    controlled_frame_indices = tuple(
        index
        for index, (param_id, _value) in enumerate(problem.independent_items)
        if param_id in controlled_ids
    )
    relevant_equalities = tuple(
        item
        for item in problem.affine_equalities
        if any(item.coefficients[index] != 0.0 for index in controlled_frame_indices)
    )
    if relevant_equalities:
        return ClaimAssessment(
            "FULL_DIMENSIONAL_FEASIBLE_INTERIOR",
            ClaimState.VIOLATED,
            "Affine equality restrictions remove a full-dimensional local interior",
        )
    accepted_controlled = dict(
        zip(problem.controlled_ids, accepted.vector, strict=True)
    )
    values = np.asarray(
        tuple(
            accepted_controlled.get(param_id, value)
            for param_id, value in problem.independent_items
        )
    )
    slacks = tuple(
        float(
            restriction.upper_bound
            - _pairwise_sum_feasibility(
                tuple(
                    coefficient * float(value)
                    for coefficient, value in zip(
                        restriction.coefficients,
                        values,
                        strict=True,
                    )
                )
            )
        )
        for restriction in problem.affine_half_spaces
        if any(
            restriction.coefficients[index] != 0.0 for index in controlled_frame_indices
        )
    )
    state = (
        ClaimState.SATISFIED
        if box_interior.state is ClaimState.SATISFIED
        and all(math.isfinite(slack) and slack > 0.0 for slack in slacks)
        else ClaimState.VIOLATED
    )
    detail = (
        "Accepted point is in the strict box and affine-half-space interior"
        if slacks
        else "Accepted point is in the strict box interior"
    )
    return ClaimAssessment("FULL_DIMENSIONAL_FEASIBLE_INTERIOR", state, detail)


def _covariance_from_jacobian(  # noqa: C901 - ordered scientific derivation gates
    accepted: AcceptedFitResult,
    jacobian: ResidualJacobianEvidence,
    *,
    problem: OptimizationProblem,
    engine: EvaluationEngine,
    policy: UncertaintyPolicy,
    request_identity: str,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> tuple[CovarianceEvidence | None, RankDiagnostic | None, EvidenceFailure | None]:
    _raise_if_terminated(cancellation_probe)
    residual_count = jacobian.residual_count
    coordinate_count = jacobian.coordinate_count
    normalization_count = sum(
        1 for profile in engine.plan.profiles if profile.is_scaled
    )
    degrees_of_freedom = residual_count - coordinate_count - normalization_count
    if residual_count - normalization_count < coordinate_count:
        return (
            None,
            None,
            EvidenceFailure(
                "covariance",
                "insufficient_effective_observations",
                "m - g is smaller than the controlled-coordinate count",
                jacobian.identity,
            ),
        )
    if (
        policy.residual_variance_scaling
        is ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE
        and degrees_of_freedom <= 0
    ):
        return (
            None,
            None,
            EvidenceFailure(
                "covariance",
                "non_positive_nominal_residual_degrees_of_freedom",
                "Estimated common residual variance requires positive nu",
                jacobian.identity,
            ),
        )
    matrix = np.asarray(jacobian.matrix, dtype=np.float64)
    scales = np.asarray(jacobian.coordinate_scales, dtype=np.float64)
    scaled = matrix * scales[np.newaxis, :]
    _raise_if_terminated(cancellation_probe)
    try:
        _left, singular_array, right_transpose = svd(
            scaled,
            full_matrices=False,
            compute_uv=True,
            overwrite_a=False,
            check_finite=True,
            lapack_driver=policy.svd_driver,
        )
    except Exception as error:  # noqa: BLE001 - declared third-party kernel fence
        return (
            None,
            None,
            EvidenceFailure(
                "covariance",
                "svd_failure",
                str(error),
                jacobian.identity,
            ),
        )
    _raise_if_terminated(cancellation_probe)
    singular_values = tuple(
        _finite(value, name=f"singular value[{index}]")
        for index, value in enumerate(singular_array)
    )
    spectrum_error = _singular_spectrum_error(singular_values, coordinate_count)
    if spectrum_error is not None:
        return (
            None,
            None,
            EvidenceFailure(
                "covariance",
                spectrum_error[0],
                spectrum_error[1],
                jacobian.identity,
            ),
        )
    largest = singular_values[0] if singular_values else 0.0
    threshold = (
        policy.rank_absolute_tolerance + policy.rank_relative_tolerance * largest
    )
    threshold = _finite(threshold, name="rank threshold")
    rank = sum(value > threshold for value in singular_values)
    normalized = tuple(
        0.0 if largest == 0.0 else value / largest for value in singular_values
    )
    column_norms = tuple(
        _finite(
            math.sqrt(
                _pairwise_sum(tuple(float(value) ** 2 for value in scaled[:, index]))
            ),
            name=f"scaled column norm[{index}]",
        )
        for index in range(coordinate_count)
    )
    # Individual vectors inside a repeated/clustered singular subspace are not
    # authoritative.  Only projectors are retained below.
    right = np.asarray(right_transpose.T, dtype=np.float64)
    weak_threshold = _finite(
        max(threshold, policy.weak_relative_tolerance * largest),
        name="weak-subspace threshold",
    )

    def projector(indices: Sequence[int]) -> tuple[tuple[float, ...], ...]:
        factor = tuple(
            tuple(
                _finite(right[row, column], name="right singular vector")
                for column in indices
            )
            for row in range(coordinate_count)
        )
        return _gram_matrix(factor)

    identifiable_projector = projector(tuple(range(rank)))
    null_projector = projector(tuple(range(rank, coordinate_count)))
    weak_projector = projector(
        tuple(
            index
            for index, value in enumerate(singular_values)
            if threshold < value <= weak_threshold
        )
    )
    subspaces = _invariant_singular_subspaces(
        right,
        singular_values,
        rank_threshold=threshold,
        weak_threshold=weak_threshold,
        cluster_relative_tolerance=(policy.singular_value_cluster_relative_tolerance),
    )
    rank_diagnostic = RankDiagnostic(
        request_identity,
        accepted.identity,
        accepted.occurrence_identity,
        jacobian.identity,
        jacobian.controlled_ids,
        singular_values,
        normalized,
        column_norms,
        threshold,
        weak_threshold,
        rank,
        identifiable_projector,
        null_projector,
        weak_projector,
        subspaces,
        source_jacobian=jacobian,
        source_policy=policy,
    )
    terminal = _cancellation_terminal(cancellation_probe)
    if terminal is not None:
        raise DerivationTermination(terminal, rank_diagnostic)
    if rank != coordinate_count:
        return (
            None,
            rank_diagnostic,
            EvidenceFailure(
                "covariance",
                "rank_deficient",
                f"Scaled Jacobian rank {rank} is below {coordinate_count}",
                jacobian.identity,
            ),
        )
    smallest = singular_values[-1]
    jacobian_condition = _finite(largest / smallest, name="Jacobian condition")
    information_condition = jacobian_condition * jacobian_condition
    if not math.isfinite(information_condition):
        return (
            None,
            rank_diagnostic,
            EvidenceFailure(
                "covariance",
                "non_finite_information_condition",
                "Squared Jacobian condition is non-finite",
                jacobian.identity,
            ),
        )
    residual_variance_scale, variance_failure = _residual_variance_scale(
        accepted,
        degrees_of_freedom=degrees_of_freedom,
        policy=policy,
    )
    if variance_failure is not None or residual_variance_scale is None:
        return None, rank_diagnostic, variance_failure
    _raise_if_terminated(cancellation_probe)
    chi_square = _finite(accepted.chi_square, name="accepted chi-square")
    unscaled_covariance, factor, covariance = _canonical_covariance_reduction(
        jacobian.matrix,
        jacobian.coordinate_scales,
        residual_variance_scale,
    )
    if _has_positive_scale_zero_variance(covariance, residual_variance_scale):
        return (
            None,
            rank_diagnostic,
            EvidenceFailure(
                "covariance",
                "covariance_factor_underflow",
                "Positive-scale full-rank covariance factor produced zero variance",
                jacobian.identity,
            ),
        )
    interior, boundary, affine, controlled_affine = _boundary_claims(
        accepted,
        problem,
        covariance,
        residual_variance_scale,
        policy.affine_feasibility_policy,
    )
    full_dimensional_interior = _full_dimensional_feasible_interior_claim(
        accepted,
        problem,
        interior,
    )
    conditioning = ClaimAssessment(
        "CONDITIONING_ADEQUACY",
        ClaimState.SATISFIED
        if jacobian_condition <= policy.conditioning_limit
        else ClaimState.VIOLATED,
        f"kappa_J={float(jacobian_condition).hex()}",
    )
    normalization = ClaimAssessment(
        "PROFILED_NORMALIZATION_REGULARITY",
        ClaimState.SATISFIED
        if _profiled_normalization_regular(accepted, engine)
        else ClaimState.VIOLATED,
    )
    variance_nondegeneracy = ClaimAssessment(
        "RESIDUAL_VARIANCE_NONDEGENERACY",
        ClaimState.NOT_APPLICABLE
        if policy.residual_variance_scaling
        is ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
        else (
            ClaimState.SATISFIED
            if residual_variance_scale > 0.0
            else ClaimState.VIOLATED
        ),
    )
    required = (
        conditioning.state,
        ClaimState.SATISFIED,
        full_dimensional_interior.state,
        interior.state,
        boundary.state,
        controlled_affine.state
        if controlled_affine.state is not ClaimState.NOT_APPLICABLE
        else ClaimState.SATISFIED,
        normalization.state,
        variance_nondegeneracy.state
        if variance_nondegeneracy.state is not ClaimState.NOT_APPLICABLE
        else ClaimState.SATISFIED,
    )
    claims = (
        ClaimAssessment("AUTHORITATIVE_LINEAGE", ClaimState.SATISFIED),
        ClaimAssessment("LOCAL_LINEARIZATION_REGULARITY", ClaimState.SATISFIED),
        ClaimAssessment("EFFECTIVE_OBSERVATION_SUFFICIENCY", ClaimState.SATISFIED),
        ClaimAssessment("FULL_COLUMN_RANK", ClaimState.SATISFIED),
        conditioning,
        full_dimensional_interior,
        interior,
        boundary,
        affine,
        controlled_affine,
        variance_nondegeneracy,
        normalization,
        ClaimAssessment("COVARIANCE_ARITHMETIC_INTEGRITY", ClaimState.SATISFIED),
        ClaimAssessment(
            "USABLE_LOCAL_COVARIANCE",
            ClaimState.SATISFIED
            if all(state is ClaimState.SATISFIED for state in required)
            else (
                ClaimState.INDETERMINATE
                if ClaimState.INDETERMINATE in required
                else ClaimState.VIOLATED
            ),
        ),
    )
    covariance_evidence = CovarianceEvidence(
        request_identity,
        jacobian.identity,
        accepted.identity,
        accepted.occurrence_identity,
        problem.identity,
        accepted.parameterization_identity,
        accepted.evaluation_result.plan_identity,
        policy.identity,
        policy.calibration_identity,
        policy.numerical_compatibility_requirement,
        accepted.controlled_ids,
        tuple(unit for _param_id, unit in policy.coordinate_units),
        accepted.vector,
        jacobian.coordinate_scales,
        policy.residual_variance_scaling,
        residual_count,
        coordinate_count,
        normalization_count,
        degrees_of_freedom,
        chi_square,
        residual_variance_scale,
        rank,
        singular_values,
        threshold,
        jacobian_condition,
        information_condition,
        unscaled_covariance,
        factor,
        covariance,
        claims,
        rank_diagnostic.identity,
        source_jacobian=jacobian,
        source_rank_diagnostic=rank_diagnostic,
        accepted_anchor=accepted,
        source_problem=problem,
        source_policy=policy,
        source_engine=engine,
    )
    return covariance_evidence, rank_diagnostic, None


def _marginal_errors(
    *,
    source_identity: str,
    accepted_result_identity: str,
    accepted_occurrence_identity: str,
    source_family: str,
    output_ids: tuple[str, ...],
    units: tuple[ParameterUnit, ...],
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_degenerate: bool,
    source_reportable: bool,
    inherited_claims: tuple[ClaimAssessment, ...],
    source_artifact: CovarianceEvidence | ConstrainedPropagationEvidence,
) -> MarginalErrorEvidence:
    entries: list[ScalarEvidenceEntry] = []
    for index, param_id in enumerate(output_ids):
        variance = covariance[index][index]
        if residual_variance_degenerate and variance == 0.0:
            entries.append(
                ScalarEvidenceEntry(
                    param_id,
                    None,
                    "RESIDUAL_VARIANCE_DEGENERACY",
                    0.0,
                )
            )
        elif variance > 0.0 and math.isfinite(variance):
            entries.append(
                ScalarEvidenceEntry(param_id, math.sqrt(variance), "AVAILABLE")
            )
        else:
            entries.append(
                ScalarEvidenceEntry(
                    param_id,
                    None,
                    "INVALID_VARIANCE",
                    variance if math.isfinite(variance) else None,
                )
            )
    complete = all(item.value is not None for item in entries)
    reportable = source_reportable and complete
    claims = (
        *inherited_claims,
        ClaimAssessment(
            "MARGINAL_VARIANCE_VALIDITY",
            ClaimState.SATISFIED if complete else ClaimState.VIOLATED,
        ),
        ClaimAssessment(
            "MARGINAL_SCOPE_REPORTABILITY",
            ClaimState.SATISFIED if reportable else ClaimState.VIOLATED,
        ),
    )
    return MarginalErrorEvidence(
        source_identity,
        accepted_result_identity,
        accepted_occurrence_identity,
        source_family,
        output_ids,
        units,
        tuple(entries),
        claims,
        reportable,
        source_artifact=source_artifact,
    )


def _correlations(
    *,
    source_identity: str,
    accepted_result_identity: str,
    accepted_occurrence_identity: str,
    source_family: str,
    output_ids: tuple[str, ...],
    units: tuple[ParameterUnit, ...],
    covariance: tuple[tuple[float, ...], ...],
    residual_variance_degenerate: bool,
    source_reportable: bool,
    policy: UncertaintyPolicy,
    inherited_claims: tuple[ClaimAssessment, ...],
    source_artifact: CovarianceEvidence | ConstrainedPropagationEvidence,
) -> CorrelationEvidence:
    entries = _expected_correlation_entries(
        covariance,
        residual_variance_degenerate=residual_variance_degenerate,
        policy=policy,
    )
    complete = all(item.value is not None for row in entries for item in row)
    reportable = source_reportable and complete
    claims = (
        *inherited_claims,
        ClaimAssessment(
            "CORRELATION_ENTRY_VALIDITY",
            ClaimState.SATISFIED if complete else ClaimState.VIOLATED,
        ),
        ClaimAssessment(
            "CORRELATION_SCOPE_REPORTABILITY",
            ClaimState.SATISFIED if reportable else ClaimState.VIOLATED,
        ),
    )
    return CorrelationEvidence(
        source_identity,
        accepted_result_identity,
        accepted_occurrence_identity,
        source_family,
        output_ids,
        units,
        entries,
        claims,
        reportable,
        source_artifact=source_artifact,
        source_policy=policy,
    )


def _combine_gradients(
    left: tuple[float, ...],
    right: tuple[float, ...],
    operation: str,
    left_value: float,
    right_value: float,
) -> tuple[float, ...]:
    if operation == "add":
        return tuple(a + b for a, b in zip(left, right, strict=True))
    if operation == "subtract":
        return tuple(a - b for a, b in zip(left, right, strict=True))
    if operation == "multiply":
        return tuple(
            a * right_value + left_value * b for a, b in zip(left, right, strict=True)
        )
    if right_value == 0.0:
        raise ArithmeticError("constraint derivative divides by zero")
    return tuple(
        (a * right_value - left_value * b) / (right_value * right_value)
        for a, b in zip(left, right, strict=True)
    )


def _scientific_function_output(
    expression: FunctionExpression,
    function: Callable[..., object],
    arguments: Sequence[float],
) -> float:
    result = function(*arguments)
    if expression.component is None:
        return _finite(result, name=f"function {expression.function_id!r} output")
    if not isinstance(result, Mapping) or expression.component not in result:
        raise ValueError(
            f"Scientific function {expression.function_id!r} lacks component "
            f"{expression.component!r}"
        )
    return _finite(
        cast("Mapping[str, object]", result)[expression.component],
        name=(
            f"function {expression.function_id!r} component {expression.component!r}"
        ),
    )


def _scientific_function_partial(
    expression: FunctionExpression,
    function: Callable[..., object],
    arguments: tuple[float, ...],
    argument_index: int,
    capability: FunctionFiniteDifferenceCapability,
    policy: UncertaintyPolicy,
) -> tuple[float, FunctionPartialDiagnostic]:
    center = arguments[argument_index]
    nominal = _NOMINAL_STEP_FACTOR * capability.argument_scales[argument_index]
    base = (_scientific_function_output(expression, function, arguments),)
    trajectory: list[object] = []
    for attempt, step in enumerate(
        _step_candidates(nominal, math.inf, policy),
        start=1,
    ):
        trial_arguments: list[tuple[float, ...]] = []
        displacements: list[float] = []
        represented = {center}
        for multiplier in (-2.0, -1.0, 1.0, 2.0):
            updated = center + multiplier * step
            if not math.isfinite(updated) or updated in represented:
                break
            represented.add(updated)
            trial = list(arguments)
            trial[argument_index] = updated
            trial_arguments.append(tuple(trial))
            displacements.append(updated - center)
        if len(trial_arguments) != 4:
            trajectory.append((_float_token(step), "not_distinct"))
            continue
        try:
            results = tuple(
                (_scientific_function_output(expression, function, trial),)
                for trial in trial_arguments
            )
        except (ArithmeticError, FloatingPointError, TypeError, ValueError):
            trajectory.append((_float_token(step), "domain_failure"))
            continue
        fine, coarse = _stencil_estimates(
            "centered",
            results,
            base,
            displacements,
        )
        discrepancy, derivative_scale, roundoff, reliable = _assess_stencil(
            fine,
            coarse,
            results,
            base,
            displacements,
            policy,
            capability.output_scale,
        )
        trajectory.append(
            (
                _float_token(step),
                tuple(_float_token(item) for item in displacements),
                _float_token(discrepancy),
                _float_token(derivative_scale),
                _float_token(roundoff),
                reliable,
            )
        )
        if reliable:
            fingerprint = _identity(
                "scientific-function-partial-trajectory", trajectory
            )
            return (
                fine[0],
                FunctionPartialDiagnostic(
                    capability.identity,
                    expression.function_id,
                    expression.component,
                    argument_index,
                    "centered_two_scale_numerical",
                    fine[0],
                    coarse[0],
                    tuple(displacements),
                    discrepancy,
                    derivative_scale,
                    roundoff,
                    attempt,
                    fingerprint,
                    (
                        ClaimAssessment(
                            "FUNCTION_PARTIAL_RELIABILITY",
                            ClaimState.SATISFIED,
                        ),
                    ),
                ),
            )
    fingerprint = _identity("scientific-function-partial-trajectory", trajectory)
    if trajectory and all(
        isinstance(item, tuple) and len(item) == 2 and item[1] == "not_distinct"
        for item in trajectory
    ):
        category = "function_partial_representation_loss"
    elif trajectory and all(
        isinstance(item, tuple)
        and len(item) == 2
        and item[1] in {"not_distinct", "domain_failure"}
        for item in trajectory
    ):
        category = "function_partial_active_domain_failure"
    else:
        category = "unreliable_or_nondifferentiable_function_partial"
    raise FunctionPartialFailure(
        category,
        f"No reliable numerical partial for function {expression.function_id!r} "
        f"argument {argument_index}",
        fingerprint,
        (
            "exhausted_declared_smaller_step_extent",
            "exhausted_declared_larger_step_extent",
        ),
    )


def _analytic_function_partials(
    expression: FunctionExpression,
    arguments: tuple[float, ...],
    capability: FunctionAnalyticPartialCapability,
) -> tuple[tuple[float, ...], tuple[FunctionPartialDiagnostic, ...]]:
    capability.validate_implementations()
    if len(capability.partials) != len(arguments):
        raise UncertaintyConstructionError(
            f"Scientific function {expression.function_id!r} capability has the "
            "wrong arity"
        )
    estimates: list[float] = []
    diagnostics: list[FunctionPartialDiagnostic] = []
    for index, partial in enumerate(capability.partials):
        try:
            estimate = _finite(
                partial(*arguments),
                name=f"analytic scientific-function partial[{index}]",
            )
        except (ArithmeticError, FloatingPointError, TypeError, ValueError) as error:
            raise FunctionPartialFailure(
                "analytic_function_partial_domain_failure",
                str(error),
                capability.identity,
                ("versioned_analytic_partial_failed",),
            ) from error
        estimates.append(estimate)
        diagnostics.append(
            FunctionPartialDiagnostic(
                capability.identity,
                expression.function_id,
                expression.component,
                index,
                "versioned_analytic_partial",
                estimate,
                None,
                (),
                0.0,
                abs(estimate),
                0.0,
                1,
                capability.identity,
                (
                    ClaimAssessment(
                        "FUNCTION_PARTIAL_RELIABILITY",
                        ClaimState.SATISFIED,
                        capability.identity,
                    ),
                ),
            )
        )
    return tuple(estimates), tuple(diagnostics)


def _differentiate_function(
    expression: FunctionExpression,
    values: Mapping[str, _DerivativeValue],
    dimension: int,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
) -> _DerivativeValue:
    arguments = tuple(
        _differentiate_expression(
            argument,
            values,
            dimension,
            parameterization,
            policy,
            compiled_capabilities,
        )
        for argument in expression.arguments
    )
    capability = next(
        (
            item
            for item in compiled_capabilities.capabilities
            if (item.function_id, item.component)
            == (expression.function_id, expression.component)
        ),
        None,
    )
    if capability is None:
        raise UncertaintyConstructionError(
            f"Scientific function {expression.function_id!r} lacks an explicitly "
            "declared linearization capability"
        )
    function = parameterization.binder[expression.function_id]
    argument_values = tuple(item.value for item in arguments)
    if isinstance(capability, FunctionFiniteDifferenceCapability):
        if len(capability.argument_scales) != len(arguments):
            raise UncertaintyConstructionError(
                f"Scientific function {expression.function_id!r} capability has the "
                "wrong arity"
            )
        partial_results = tuple(
            _scientific_function_partial(
                expression,
                function,
                argument_values,
                index,
                capability,
                policy,
            )
            for index in range(len(arguments))
        )
        partials = tuple(item[0] for item in partial_results)
        diagnostics = tuple(item[1] for item in partial_results)
    else:
        partials, diagnostics = _analytic_function_partials(
            expression,
            argument_values,
            capability,
        )
    value = _scientific_function_output(expression, function, argument_values)
    gradient = tuple(
        _finite(
            _pairwise_sum(
                tuple(
                    partial * arguments[index].gradient[coordinate]
                    for index, partial in enumerate(partials)
                )
            ),
            name="scientific-function constraint derivative",
        )
        for coordinate in range(dimension)
    )
    structural = frozenset().union(
        *(item.structural_dependencies for item in arguments)
    )
    inherited_diagnostics = tuple(
        diagnostic
        for argument in arguments
        for diagnostic in argument.function_partial_diagnostics
    )
    return _DerivativeValue(
        value,
        gradient,
        structural,
        inherited_diagnostics + diagnostics,
    )


def _differentiate_expression(
    expression: ScalarExpression,
    values: Mapping[str, _DerivativeValue],
    dimension: int,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
) -> _DerivativeValue:
    zero = (0.0,) * dimension
    if isinstance(expression, LiteralExpression):
        return _DerivativeValue(expression.value, zero, frozenset())
    if isinstance(expression, ReferenceExpression):
        return values[expression.param_id]
    if isinstance(expression, UnaryExpression):
        operand = _differentiate_expression(
            expression.operand,
            values,
            dimension,
            parameterization,
            policy,
            compiled_capabilities,
        )
        if expression.operator == "positive":
            return operand
        return _DerivativeValue(
            -operand.value,
            tuple(-item for item in operand.gradient),
            operand.structural_dependencies,
            operand.function_partial_diagnostics,
        )
    if isinstance(expression, BinaryExpression):
        left = _differentiate_expression(
            expression.left,
            values,
            dimension,
            parameterization,
            policy,
            compiled_capabilities,
        )
        right = _differentiate_expression(
            expression.right,
            values,
            dimension,
            parameterization,
            policy,
            compiled_capabilities,
        )
        if expression.operator == "add":
            value = left.value + right.value
        elif expression.operator == "subtract":
            value = left.value - right.value
        elif expression.operator == "multiply":
            value = left.value * right.value
        else:
            value = left.value / right.value
        return _DerivativeValue(
            _finite(value, name="constraint linearization value"),
            tuple(
                _finite(item, name="constraint derivative")
                for item in _combine_gradients(
                    left.gradient,
                    right.gradient,
                    expression.operator,
                    left.value,
                    right.value,
                )
            ),
            left.structural_dependencies | right.structural_dependencies,
            left.function_partial_diagnostics + right.function_partial_diagnostics,
        )
    if isinstance(expression, FunctionExpression):
        return _differentiate_function(
            expression,
            values,
            dimension,
            parameterization,
            policy,
            compiled_capabilities,
        )
    raise TypeError("Unknown constraint expression")


def _constraint_scope_failure(
    accepted: AcceptedFitResult,
    parameterization: ActiveParameterization,
    output_scope: tuple[str, ...],
) -> EvidenceFailure | None:
    if len(set(output_scope)) != len(output_scope):
        return EvidenceFailure(
            "constraint_linearization",
            "duplicate_output_scope",
            "Constraint-output scope must contain unique stable IDs",
            accepted.identity,
        )
    if any(param_id not in parameterization.scope_ids for param_id in output_scope):
        return EvidenceFailure(
            "constraint_linearization",
            "unknown_output_scope",
            "Constraint-output scope contains an unknown stable ID",
            accepted.identity,
        )
    return None


def _initial_derivative_values(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
) -> dict[str, _DerivativeValue]:
    controlled_index = {
        param_id: index for index, param_id in enumerate(problem.controlled_ids)
    }
    dimension = len(problem.controlled_ids)
    zero = (0.0,) * dimension
    resolved = accepted.evaluation_result.resolved_values
    values: dict[str, _DerivativeValue] = {}
    for param_id in parameterization.independent_ids:
        if param_id in controlled_index:
            gradient = tuple(
                1.0 if index == controlled_index[param_id] else 0.0
                for index in range(dimension)
            )
            structural = frozenset((param_id,))
        else:
            gradient = zero
            structural = frozenset()
        values[param_id] = _DerivativeValue(resolved[param_id], gradient, structural)
    return values


def _differentiate_target(
    target_id: str,
    *,
    values: dict[str, _DerivativeValue],
    constraints: Mapping[str, CompiledConstraint],
    resolved: Mapping[str, float],
    dimension: int,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
) -> _DerivativeValue:
    existing = values.get(target_id)
    if existing is not None:
        return existing
    constraint = constraints[target_id]
    dependencies = constraint.dependencies
    for dependency in dependencies:
        _differentiate_target(
            dependency,
            values=values,
            constraints=constraints,
            resolved=resolved,
            dimension=dimension,
            parameterization=parameterization,
            policy=policy,
            compiled_capabilities=compiled_capabilities,
        )
    differentiated = _differentiate_expression(
        constraint.expression,
        values,
        dimension,
        parameterization,
        policy,
        compiled_capabilities,
    )
    if differentiated.value != resolved[target_id]:
        raise UncertaintyConstructionError(
            f"Constraint linearization value differs for {target_id}"
        )
    values[target_id] = differentiated
    return differentiated


def _constraint_rows(
    accepted: AcceptedFitResult,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
    output_scope: tuple[str, ...],
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> tuple[
    tuple[tuple[float, ...], ...] | None,
    tuple[tuple[str, ...], ...] | None,
    tuple[FunctionPartialDiagnostic, ...] | None,
    EvidenceFailure | None,
]:
    values = _initial_derivative_values(accepted, problem, parameterization)
    constraints = {
        constraint.target_id: constraint
        for constraint in parameterization.program.constraints
    }
    resolved = accepted.evaluation_result.resolved_values
    rows: list[tuple[float, ...]] = []
    dependencies: list[tuple[str, ...]] = []
    diagnostics: dict[str, FunctionPartialDiagnostic] = {}
    try:
        for param_id in output_scope:
            _raise_if_terminated(cancellation_probe)
            role = parameterization.role(param_id)
            if (
                param_id in parameterization.independent_ids
                and param_id not in problem.controlled_ids
            ):
                return (
                    None,
                    None,
                    None,
                    EvidenceFailure(
                        "constraint_linearization",
                        "held_independent_outside_uncertainty_basis",
                        f"Held parameter {param_id} cannot be requested",
                        accepted.identity,
                    ),
                )
            value = _differentiate_target(
                param_id,
                values=values,
                constraints=constraints,
                resolved=resolved,
                dimension=len(problem.controlled_ids),
                parameterization=parameterization,
                policy=policy,
                compiled_capabilities=compiled_capabilities,
            )
            _raise_if_terminated(cancellation_probe)
            if role is ParameterRole.DERIVED and not value.structural_dependencies:
                return (
                    None,
                    None,
                    None,
                    EvidenceFailure(
                        "constraint_linearization",
                        "derived_output_outside_controlled_basis",
                        f"Derived parameter {param_id} has no controlled dependency path",
                        accepted.identity,
                    ),
                )
            rows.append(value.gradient)
            dependencies.append(
                tuple(
                    item
                    for item in problem.controlled_ids
                    if item in value.structural_dependencies
                )
            )
            diagnostics.update(
                (item.identity, item) for item in value.function_partial_diagnostics
            )
    except FunctionPartialFailure as error:
        return (
            None,
            None,
            None,
            EvidenceFailure(
                "constraint_linearization",
                error.category,
                str(error),
                accepted.identity,
                error.trajectory_fingerprint,
                error.termination_details,
            ),
        )
    except (ArithmeticError, KeyError, TypeError, ValueError) as error:
        return (
            None,
            None,
            None,
            EvidenceFailure(
                "constraint_linearization",
                "constraint_linearization_failure",
                str(error),
                accepted.identity,
            ),
        )
    return tuple(rows), tuple(dependencies), tuple(diagnostics.values()), None


def _linearize_constraints(
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
    output_scope: tuple[str, ...],
    output_units: tuple[ParameterUnit, ...],
    output_scales: tuple[float, ...],
    request_identity: str,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> tuple[ConstraintJacobianEvidence | None, EvidenceFailure | None]:
    if not output_scope:
        return None, None
    scope_failure = _constraint_scope_failure(
        accepted,
        parameterization,
        output_scope,
    )
    if scope_failure is not None:
        return None, scope_failure
    rows, dependencies, diagnostics, row_failure = _constraint_rows(
        accepted,
        problem,
        parameterization,
        policy,
        compiled_capabilities,
        output_scope,
        cancellation_probe,
    )
    if (
        row_failure is not None
        or rows is None
        or dependencies is None
        or diagnostics is None
    ):
        return None, row_failure
    return (
        ConstraintJacobianEvidence(
            request_identity,
            accepted.identity,
            accepted.occurrence_identity,
            accepted.evaluation_result.identity,
            problem.identity,
            parameterization.identity,
            parameterization.program.fingerprint,
            compiled_capabilities.identity,
            problem.controlled_ids,
            output_scope,
            output_units,
            output_scales,
            rows,
            dependencies,
            diagnostics,
            (
                ClaimAssessment(
                    "CONSTRAINT_OUTPUT_LINEARIZATION_COMPLETE",
                    ClaimState.SATISFIED,
                ),
                ClaimAssessment(
                    "CONSTRAINT_OUTPUT_LINEARIZATION_RELIABLE",
                    ClaimState.SATISFIED,
                ),
            ),
            accepted_anchor=accepted,
            source_problem=problem,
            source_parameterization=parameterization,
            source_policy=policy,
            source_capabilities=compiled_capabilities,
        ),
        None,
    )


def _propagate_constraints(
    accepted: AcceptedFitResult,
    covariance: CovarianceEvidence,
    constraint_jacobian: ConstraintJacobianEvidence,
    *,
    request_identity: str,
    policy: UncertaintyPolicy,
) -> ConstrainedPropagationEvidence:
    gradient = constraint_jacobian.matrix
    factor = tuple(
        tuple(
            _finite(
                _pairwise_sum(
                    tuple(
                        gradient[row][inner] * covariance.factor[inner][column]
                        for inner in range(len(covariance.controlled_ids))
                    )
                ),
                name=f"propagated factor[{row},{column}]",
            )
            for column in range(len(covariance.controlled_ids))
        )
        for row in range(len(constraint_jacobian.output_ids))
    )
    propagated = _gram_matrix(factor)
    per_output_degeneracy = tuple(
        all(value == 0.0 for value in constraint_jacobian.matrix[index])
        and bool(constraint_jacobian.structural_dependencies[index])
        for index in range(len(factor))
    )
    local_degeneracy = any(per_output_degeneracy)
    residual_variance_degenerate = covariance.residual_variance_scale == 0.0
    if not residual_variance_degenerate:
        for index, row in enumerate(factor):
            if any(value != 0.0 for value in gradient[index]) and all(
                value == 0.0 for value in row
            ):
                raise ArithmeticError(
                    "Nonzero constraint-gradient row collapsed to a zero propagated factor"
                )
            if any(value != 0.0 for value in row) and propagated[index][index] == 0.0:
                raise ArithmeticError(
                    "Nonzero propagated factor row underflowed to zero variance"
                )
    scaled_gradient = (
        np.asarray(gradient, dtype=np.float64)
        * np.asarray(covariance.coordinate_scales, dtype=np.float64)[np.newaxis, :]
        / np.asarray(constraint_jacobian.output_scales, dtype=np.float64)[:, np.newaxis]
    )
    output_rank_claim: ClaimAssessment
    try:
        gradient_singular = tuple(
            _finite(value, name="constraint-output singular value")
            for value in svd(
                scaled_gradient,
                full_matrices=False,
                compute_uv=False,
                overwrite_a=False,
                check_finite=True,
                lapack_driver=policy.svd_driver,
            )
        )
        if any(value < 0.0 for value in gradient_singular) or any(
            left < right for left, right in pairwise(gradient_singular)
        ):
            raise ValueError("invalid scaled constraint-output singular spectrum")
        largest = gradient_singular[0] if gradient_singular else 0.0
        rank_threshold = _finite(
            policy.rank_absolute_tolerance + policy.rank_relative_tolerance * largest,
            name="constraint-output rank threshold",
        )
        output_rank = sum(value > rank_threshold for value in gradient_singular)
        output_rank_deficient = output_rank < len(constraint_jacobian.output_ids)
        output_rank_claim = ClaimAssessment(
            "OUTPUT_RANK_DEFICIENCY_EXPECTED",
            ClaimState.SATISFIED
            if output_rank_deficient
            else ClaimState.NOT_APPLICABLE,
            f"scaled-rank={output_rank}; rows={len(constraint_jacobian.output_ids)}; "
            f"threshold={float(rank_threshold).hex()}",
        )
    except Exception as error:  # noqa: BLE001 - diagnostic kernel fence
        output_rank_claim = ClaimAssessment(
            "OUTPUT_RANK_DEFICIENCY_EXPECTED",
            ClaimState.INDETERMINATE,
            f"diagnostic scaled-SVD unavailable: {error}",
        )
    inherited_claims = tuple(
        ClaimAssessment(
            f"SOURCE_COVARIANCE::{item.name}",
            item.state,
            item.detail,
        )
        for item in covariance.claims
    ) + tuple(
        ClaimAssessment(
            f"SOURCE_CONSTRAINT_LINEARIZATION::{item.name}",
            item.state,
            item.detail,
        )
        for item in constraint_jacobian.claims
    )
    output_claims = tuple(
        ClaimAssessment(
            f"OUTPUT_FIRST_ORDER_NONDEGENERACY::{param_id}",
            ClaimState.VIOLATED
            if per_output_degeneracy[index]
            else ClaimState.SATISFIED,
        )
        for index, param_id in enumerate(constraint_jacobian.output_ids)
    )
    claims = (
        *inherited_claims,
        ClaimAssessment("EXACT_LINEAGE", ClaimState.SATISFIED),
        ClaimAssessment(
            "CONSTRAINT_LINEARIZATION_REGULARITY",
            ClaimState.SATISFIED,
        ),
        ClaimAssessment("GRAM_ARITHMETIC_INTEGRITY", ClaimState.SATISFIED),
        output_rank_claim,
        ClaimAssessment(
            "LOCAL_FIRST_ORDER_DEGENERACY",
            ClaimState.VIOLATED if local_degeneracy else ClaimState.SATISFIED,
        ),
        ClaimAssessment(
            "RESIDUAL_VARIANCE_NONDEGENERACY",
            ClaimState.VIOLATED
            if residual_variance_degenerate
            else ClaimState.SATISFIED,
        ),
        ClaimAssessment(
            "SOURCE_COVARIANCE_USABILITY",
            covariance.claim("USABLE_LOCAL_COVARIANCE"),
        ),
        *output_claims,
    )
    return ConstrainedPropagationEvidence(
        request_identity,
        covariance.identity,
        constraint_jacobian.identity,
        accepted.identity,
        accepted.occurrence_identity,
        covariance.controlled_ids,
        constraint_jacobian.output_ids,
        constraint_jacobian.output_units,
        constraint_jacobian.output_scales,
        factor,
        propagated,
        claims,
        source_covariance=covariance,
        source_constraint_jacobian=constraint_jacobian,
        accepted_anchor=accepted,
    )


def _operation(
    stage: str,
    request_identity: str,
    artifact: _IdentityBearing | None,
    failure: EvidenceFailure | None,
    *,
    resolved_environment_identity: str,
    terminal_override: OperationTerminal | None = None,
) -> DerivationOperation:
    if terminal_override in {
        OperationTerminal.CANCELLED,
        OperationTerminal.INTERRUPTED,
    }:
        terminal_failure = EvidenceFailure(
            stage,
            terminal_override.value,
            f"Derivation operation {terminal_override.value}",
            request_identity,
        )
        return DerivationOperation(
            uuid4().hex,
            resolved_environment_identity,
            stage,
            request_identity,
            terminal_override,
            failure=terminal_failure,
        )
    if artifact is not None:
        return DerivationOperation(
            uuid4().hex,
            resolved_environment_identity,
            stage,
            request_identity,
            OperationTerminal.COMPLETED,
            artifact.identity,
        )
    if failure is None:
        failure = EvidenceFailure(
            stage,
            "source_artifact_unavailable",
            "A required source artifact was unavailable",
            request_identity,
        )
    return DerivationOperation(
        uuid4().hex,
        resolved_environment_identity,
        stage,
        request_identity,
        OperationTerminal.FAILED,
        failure=failure,
    )


def _record_operation(
    operations: list[DerivationOperation],
    failures: list[EvidenceFailure],
    operation: DerivationOperation,
) -> None:
    operations.append(operation)
    if operation.failure is not None:
        failures.append(operation.failure)


def _cancellation_terminal(
    probe: Callable[[], OperationTerminal | None] | None,
) -> OperationTerminal | None:
    if probe is None:
        return None
    terminal = probe()
    if terminal not in {
        None,
        OperationTerminal.CANCELLED,
        OperationTerminal.INTERRUPTED,
    }:
        raise UncertaintyConstructionError(
            "Cancellation probe may return only cancelled, interrupted, or None"
        )
    return terminal


def _raise_if_terminated(
    probe: Callable[[], OperationTerminal | None] | None,
) -> None:
    terminal = _cancellation_terminal(probe)
    if terminal is not None:
        raise DerivationTermination(terminal)


@dataclass(frozen=True, slots=True)
class _CovarianceBranch:
    operations: tuple[DerivationOperation, ...]
    failures: tuple[EvidenceFailure, ...]
    residual_jacobian: ResidualJacobianEvidence | None
    rank_diagnostic: RankDiagnostic | None
    covariance: CovarianceEvidence | None
    marginal_errors: MarginalErrorEvidence | None
    correlations: CorrelationEvidence | None


def _derive_covariance_branch(  # noqa: C901 - ordered cancellation phase ledger
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    policy: UncertaintyPolicy,
    request_identity: str,
    resolved_environment_identity: str,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> _CovarianceBranch:
    operations: list[DerivationOperation] = []
    failures: list[EvidenceFailure] = []
    residual_request = _identity(
        "native-residual-linearization-request",
        (
            request_identity,
            accepted.identity,
            accepted.occurrence_identity,
            engine.plan.identity,
            policy.identity,
        ),
    )
    residual_jacobian: ResidualJacobianEvidence | None = None
    residual_failure: EvidenceFailure | None = None
    residual_terminal: OperationTerminal | None = None
    try:
        residual_jacobian, residual_failure = _linearize_residuals(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=policy,
            request_identity=residual_request,
            cancellation_probe=cancellation_probe,
        )
    except DerivationTermination as termination:
        residual_terminal = termination.terminal
    except (ArithmeticError, TypeError, ValueError) as error:
        residual_failure = EvidenceFailure(
            "residual_linearization",
            "invalid_linearization_arithmetic",
            str(error),
            accepted.identity,
        )
    _record_operation(
        operations,
        failures,
        _operation(
            "residual_linearization",
            residual_request,
            residual_jacobian,
            residual_failure,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=residual_terminal,
        ),
    )
    if residual_terminal is not None:
        return _CovarianceBranch(
            tuple(operations),
            tuple(failures),
            None,
            None,
            None,
            None,
            None,
        )
    covariance_request = _identity(
        "native-covariance-request",
        (
            request_identity,
            None if residual_jacobian is None else residual_jacobian.identity,
            policy.identity,
        ),
    )
    covariance: CovarianceEvidence | None = None
    rank_diagnostic: RankDiagnostic | None = None
    covariance_failure: EvidenceFailure | None = None
    covariance_terminal: OperationTerminal | None = None
    if residual_jacobian is not None:
        try:
            covariance, rank_diagnostic, covariance_failure = _covariance_from_jacobian(
                accepted,
                residual_jacobian,
                problem=problem,
                engine=engine,
                policy=policy,
                request_identity=covariance_request,
                cancellation_probe=cancellation_probe,
            )
        except DerivationTermination as termination:
            covariance_terminal = termination.terminal
            if isinstance(termination.partial_artifact, RankDiagnostic):
                rank_diagnostic = termination.partial_artifact
        except (ArithmeticError, TypeError, ValueError) as error:
            covariance_failure = EvidenceFailure(
                "covariance",
                "invalid_covariance_arithmetic",
                str(error),
                residual_jacobian.identity,
            )
    _record_operation(
        operations,
        failures,
        _operation(
            "covariance",
            covariance_request,
            covariance,
            covariance_failure,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=covariance_terminal,
        ),
    )
    if covariance_terminal is not None:
        return _CovarianceBranch(
            tuple(operations),
            tuple(failures),
            residual_jacobian,
            rank_diagnostic,
            None,
            None,
            None,
        )
    marginal_errors: MarginalErrorEvidence | None = None
    correlations: CorrelationEvidence | None = None
    if covariance is not None:
        terminal = _cancellation_terminal(cancellation_probe)
        if terminal is not None:
            operation = _operation(
                "marginal_errors",
                covariance.identity,
                None,
                None,
                resolved_environment_identity=resolved_environment_identity,
                terminal_override=terminal,
            )
            _record_operation(operations, failures, operation)
            return _CovarianceBranch(
                tuple(operations),
                tuple(failures),
                residual_jacobian,
                rank_diagnostic,
                covariance,
                None,
                None,
            )
        variance_degenerate = covariance.residual_variance_scale == 0.0
        marginal_errors = _marginal_errors(
            source_identity=covariance.identity,
            accepted_result_identity=accepted.identity,
            accepted_occurrence_identity=accepted.occurrence_identity,
            source_family="local_covariance",
            output_ids=covariance.controlled_ids,
            units=covariance.units,
            covariance=covariance.covariance,
            residual_variance_degenerate=variance_degenerate,
            source_reportable=covariance.usable,
            inherited_claims=covariance.claims,
            source_artifact=covariance,
        )
        marginal_request = _identity(
            "native-marginal-error-request",
            (covariance.identity, _MARGINAL_VERSION),
        )
        _record_operation(
            operations,
            failures,
            _operation(
                "marginal_errors",
                marginal_request,
                marginal_errors,
                None,
                resolved_environment_identity=resolved_environment_identity,
            ),
        )
        correlation_request = _identity(
            "native-correlation-request",
            (covariance.identity, policy.identity, _CORRELATION_VERSION),
        )
        terminal = _cancellation_terminal(cancellation_probe)
        if terminal is not None:
            _record_operation(
                operations,
                failures,
                _operation(
                    "correlations",
                    correlation_request,
                    None,
                    None,
                    resolved_environment_identity=resolved_environment_identity,
                    terminal_override=terminal,
                ),
            )
            return _CovarianceBranch(
                tuple(operations),
                tuple(failures),
                residual_jacobian,
                rank_diagnostic,
                covariance,
                marginal_errors,
                None,
            )
        correlations = _correlations(
            source_identity=covariance.identity,
            accepted_result_identity=accepted.identity,
            accepted_occurrence_identity=accepted.occurrence_identity,
            source_family="local_covariance",
            output_ids=covariance.controlled_ids,
            units=covariance.units,
            covariance=covariance.covariance,
            residual_variance_degenerate=variance_degenerate,
            source_reportable=covariance.usable,
            policy=policy,
            inherited_claims=covariance.claims,
            source_artifact=covariance,
        )
        _record_operation(
            operations,
            failures,
            _operation(
                "correlations",
                correlation_request,
                correlations,
                None,
                resolved_environment_identity=resolved_environment_identity,
            ),
        )
    return _CovarianceBranch(
        tuple(operations),
        tuple(failures),
        residual_jacobian,
        rank_diagnostic,
        covariance,
        marginal_errors,
        correlations,
    )


@dataclass(frozen=True, slots=True)
class _ConstraintBranch:
    operations: tuple[DerivationOperation, ...]
    failures: tuple[EvidenceFailure, ...]
    jacobian: ConstraintJacobianEvidence | None
    propagation: ConstrainedPropagationEvidence | None
    marginal_errors: MarginalErrorEvidence | None
    correlations: CorrelationEvidence | None


def _derive_constraint_branch(
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    policy: UncertaintyPolicy,
    compiled_capabilities: CompiledConstraintLinearizationCapabilities,
    request_identity: str,
    output_scope: tuple[str, ...],
    output_units: tuple[ParameterUnit, ...],
    output_scales: tuple[float, ...],
    covariance: CovarianceEvidence | None,
    resolved_environment_identity: str,
    cancellation_probe: Callable[[], OperationTerminal | None] | None,
) -> _ConstraintBranch:
    if not output_scope:
        return _ConstraintBranch((), (), None, None, None, None)
    operations: list[DerivationOperation] = []
    failures: list[EvidenceFailure] = []
    constraint_request = _identity(
        "native-constraint-linearization-request",
        (
            request_identity,
            accepted.identity,
            accepted.occurrence_identity,
            output_scope,
            output_units,
            _vector_tokens(output_scales),
            policy.identity,
        ),
    )
    jacobian: ConstraintJacobianEvidence | None = None
    constraint_failure: EvidenceFailure | None = None
    constraint_terminal: OperationTerminal | None = None
    try:
        _raise_if_terminated(cancellation_probe)
        jacobian, constraint_failure = _linearize_constraints(
            accepted,
            problem=problem,
            parameterization=parameterization,
            policy=policy,
            compiled_capabilities=compiled_capabilities,
            output_scope=output_scope,
            output_units=output_units,
            output_scales=output_scales,
            request_identity=constraint_request,
            cancellation_probe=cancellation_probe,
        )
    except DerivationTermination as termination:
        constraint_terminal = termination.terminal
    _record_operation(
        operations,
        failures,
        _operation(
            "constraint_linearization",
            constraint_request,
            jacobian,
            constraint_failure,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=constraint_terminal,
        ),
    )
    if constraint_terminal is not None:
        return _ConstraintBranch(
            tuple(operations), tuple(failures), None, None, None, None
        )
    propagation_request = _identity(
        "native-constrained-propagation-request",
        (
            request_identity,
            None if covariance is None else covariance.identity,
            None if jacobian is None else jacobian.identity,
            output_scope,
            policy.identity,
        ),
    )
    propagation: ConstrainedPropagationEvidence | None = None
    propagation_failure: EvidenceFailure | None = None
    propagation_terminal: OperationTerminal | None = None
    if covariance is not None and jacobian is not None:
        try:
            _raise_if_terminated(cancellation_probe)
            propagation = _propagate_constraints(
                accepted,
                covariance,
                jacobian,
                request_identity=propagation_request,
                policy=policy,
            )
            _raise_if_terminated(cancellation_probe)
        except DerivationTermination as termination:
            propagation_terminal = termination.terminal
        except (ArithmeticError, TypeError, ValueError) as error:
            propagation_failure = EvidenceFailure(
                "constrained_propagation",
                "gram_propagation_failure",
                str(error),
                propagation_request,
            )
    _record_operation(
        operations,
        failures,
        _operation(
            "constrained_propagation",
            propagation_request,
            propagation,
            propagation_failure,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=propagation_terminal,
        ),
    )
    if propagation_terminal is not None:
        return _ConstraintBranch(
            tuple(operations), tuple(failures), jacobian, None, None, None
        )
    marginal: MarginalErrorEvidence | None = None
    correlations: CorrelationEvidence | None = None
    if propagation is not None and covariance is not None:
        variance_degenerate = covariance.residual_variance_scale == 0.0
        source_reportable = (
            covariance.usable
            and propagation.claim("LOCAL_FIRST_ORDER_DEGENERACY")
            is ClaimState.SATISFIED
        )
        marginal = _marginal_errors(
            source_identity=propagation.identity,
            accepted_result_identity=accepted.identity,
            accepted_occurrence_identity=accepted.occurrence_identity,
            source_family="constrained_propagation",
            output_ids=propagation.output_ids,
            units=propagation.output_units,
            covariance=propagation.covariance,
            residual_variance_degenerate=variance_degenerate,
            source_reportable=source_reportable,
            inherited_claims=propagation.claims,
            source_artifact=propagation,
        )
        correlations = _correlations(
            source_identity=propagation.identity,
            accepted_result_identity=accepted.identity,
            accepted_occurrence_identity=accepted.occurrence_identity,
            source_family="constrained_propagation",
            output_ids=propagation.output_ids,
            units=propagation.output_units,
            covariance=propagation.covariance,
            residual_variance_degenerate=variance_degenerate,
            source_reportable=source_reportable,
            policy=policy,
            inherited_claims=propagation.claims,
            source_artifact=propagation,
        )
        marginal_request = _identity(
            "native-constrained-marginal-error-request",
            (propagation.identity, _MARGINAL_VERSION),
        )
        correlation_request = _identity(
            "native-constrained-correlation-request",
            (propagation.identity, policy.identity, _CORRELATION_VERSION),
        )
        _record_operation(
            operations,
            failures,
            _operation(
                "constrained_marginal_errors",
                marginal_request,
                marginal,
                None,
                resolved_environment_identity=resolved_environment_identity,
            ),
        )
        _record_operation(
            operations,
            failures,
            _operation(
                "constrained_correlations",
                correlation_request,
                correlations,
                None,
                resolved_environment_identity=resolved_environment_identity,
            ),
        )
    return _ConstraintBranch(
        tuple(operations),
        tuple(failures),
        jacobian,
        propagation,
        marginal,
        correlations,
    )


def derive_uncertainty_evidence(  # noqa: C901 - fail-closed phase orchestration
    accepted: AcceptedFitResult,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    policy: UncertaintyPolicy,
    constrained_scope: Sequence[str] = (),
    constrained_units: Sequence[tuple[str, ParameterUnit]] = (),
    constrained_scales: Sequence[tuple[str, float]] = (),
    compiled_constraint_linearization: (
        CompiledConstraintLinearizationCapabilities | None
    ) = None,
    cancellation_probe: Callable[[], OperationTerminal | None] | None = None,
    resolved_environment_identity: str,
) -> UncertaintyEvidence:
    """Derive one closed immutable evidence bundle from an exact accepted fit."""
    if not resolved_environment_identity:
        raise UncertaintyConstructionError(
            "Uncertainty evidence requires a resolved environment identity"
        )
    output_scope = tuple(constrained_scope)
    compiled_capabilities = (
        compile_constraint_linearization_capabilities(
            parameterization,
            (),
            (),
        )
        if compiled_constraint_linearization is None and not output_scope
        else compiled_constraint_linearization
    )
    if compiled_capabilities is None or (
        compiled_capabilities.parameterization_identity != parameterization.identity
        or compiled_capabilities.constraint_program_identity
        != parameterization.program.fingerprint
        or compiled_capabilities.output_scope != output_scope
    ):
        raise UncertaintyConstructionError(
            "Constraint linearization requires an exact compiled capability selection"
        )
    unit_pairs = tuple(constrained_units)
    if tuple(param_id for param_id, _unit in unit_pairs) != output_scope or any(
        not unit for _param_id, unit in unit_pairs
    ):
        raise UncertaintyConstructionError(
            "Constrained-output units must explicitly match requested scope order"
        )
    output_units = tuple(unit for _param_id, unit in unit_pairs)
    scale_pairs = tuple(
        (param_id, _finite(value, name=f"constrained scale {param_id!r}"))
        for param_id, value in constrained_scales
    )
    if tuple(param_id for param_id, _value in scale_pairs) != output_scope or any(
        value <= 0.0 for _param_id, value in scale_pairs
    ):
        raise UncertaintyConstructionError(
            "Constrained-output scales must explicitly match requested scope order"
        )
    output_scales = tuple(value for _param_id, value in scale_pairs)
    request_identity = _identity(
        "native-uncertainty-request",
        (
            accepted.identity,
            accepted.occurrence_identity,
            problem.identity,
            parameterization.identity,
            engine.plan.identity,
            policy.identity,
            output_scope,
            output_units,
            _vector_tokens(output_scales),
            compiled_capabilities.identity,
        ),
    )
    initial_terminal = _cancellation_terminal(cancellation_probe)
    if initial_terminal is not None:
        operation = _operation(
            "residual_linearization",
            request_identity,
            None,
            None,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=initial_terminal,
        )
        terminal_failure = cast("EvidenceFailure", operation.failure)
        return UncertaintyEvidence(
            accepted.identity,
            accepted.occurrence_identity,
            request_identity,
            policy.identity,
            resolved_environment_identity,
            (operation,),
            (terminal_failure,),
            accepted_anchor=accepted,
            requested_output_scope=output_scope,
            requested_output_units=output_units,
            requested_output_scales=output_scales,
            source_problem=problem,
            source_parameterization=parameterization,
            source_engine=engine,
            source_policy=policy,
            source_capabilities=compiled_capabilities,
        )
    lineage_category = _lineage_failure(
        accepted,
        problem,
        parameterization,
        engine,
        policy,
    )
    if lineage_category is not None:
        failure = EvidenceFailure(
            "residual_linearization",
            lineage_category,
            "Uncertainty inputs do not have exact accepted-fit lineage",
            accepted.identity,
        )
        operation = _operation(
            "residual_linearization",
            request_identity,
            None,
            failure,
            resolved_environment_identity=resolved_environment_identity,
        )
        return UncertaintyEvidence(
            accepted.identity,
            accepted.occurrence_identity,
            request_identity,
            policy.identity,
            resolved_environment_identity,
            (operation,),
            (failure,),
            accepted_anchor=accepted,
            requested_output_scope=output_scope,
            requested_output_units=output_units,
            requested_output_scales=output_scales,
            source_problem=problem,
            source_parameterization=parameterization,
            source_engine=engine,
            source_policy=policy,
            source_capabilities=compiled_capabilities,
        )

    covariance_branch = _derive_covariance_branch(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        request_identity=request_identity,
        resolved_environment_identity=resolved_environment_identity,
        cancellation_probe=cancellation_probe,
    )
    if covariance_branch.operations[-1].terminal in {
        OperationTerminal.CANCELLED,
        OperationTerminal.INTERRUPTED,
    }:
        return UncertaintyEvidence(
            accepted.identity,
            accepted.occurrence_identity,
            request_identity,
            policy.identity,
            resolved_environment_identity,
            covariance_branch.operations,
            covariance_branch.failures,
            covariance_branch.residual_jacobian,
            covariance_branch.rank_diagnostic,
            covariance_branch.covariance,
            covariance_branch.marginal_errors,
            covariance_branch.correlations,
            accepted_anchor=accepted,
            requested_output_scope=output_scope,
            requested_output_units=output_units,
            requested_output_scales=output_scales,
            source_problem=problem,
            source_parameterization=parameterization,
            source_engine=engine,
            source_policy=policy,
            source_capabilities=compiled_capabilities,
        )
    branch_terminal = (
        _cancellation_terminal(cancellation_probe) if output_scope else None
    )
    if branch_terminal is not None:
        operation = _operation(
            "constraint_linearization",
            request_identity,
            None,
            None,
            resolved_environment_identity=resolved_environment_identity,
            terminal_override=branch_terminal,
        )
        terminal_failure = cast("EvidenceFailure", operation.failure)
        return UncertaintyEvidence(
            accepted.identity,
            accepted.occurrence_identity,
            request_identity,
            policy.identity,
            resolved_environment_identity,
            covariance_branch.operations + (operation,),
            covariance_branch.failures + (terminal_failure,),
            covariance_branch.residual_jacobian,
            covariance_branch.rank_diagnostic,
            covariance_branch.covariance,
            covariance_branch.marginal_errors,
            covariance_branch.correlations,
            accepted_anchor=accepted,
            requested_output_scope=output_scope,
            requested_output_units=output_units,
            requested_output_scales=output_scales,
            source_problem=problem,
            source_parameterization=parameterization,
            source_engine=engine,
            source_policy=policy,
            source_capabilities=compiled_capabilities,
        )
    constraint_branch = _derive_constraint_branch(
        accepted,
        problem=problem,
        parameterization=parameterization,
        policy=policy,
        compiled_capabilities=compiled_capabilities,
        request_identity=request_identity,
        output_scope=output_scope,
        output_units=output_units,
        output_scales=output_scales,
        covariance=covariance_branch.covariance,
        resolved_environment_identity=resolved_environment_identity,
        cancellation_probe=cancellation_probe,
    )

    combined_operations = covariance_branch.operations + constraint_branch.operations
    combined_failures = covariance_branch.failures + constraint_branch.failures
    constraint_terminated = bool(constraint_branch.operations) and (
        constraint_branch.operations[-1].terminal
        in {OperationTerminal.CANCELLED, OperationTerminal.INTERRUPTED}
    )
    if not constraint_terminated:
        final_terminal = _cancellation_terminal(cancellation_probe)
        if final_terminal is not None:
            final_operation = _operation(
                "final_bundle_assembly",
                request_identity,
                None,
                None,
                resolved_environment_identity=resolved_environment_identity,
                terminal_override=final_terminal,
            )
            combined_operations += (final_operation,)
            combined_failures += (cast("EvidenceFailure", final_operation.failure),)

    return UncertaintyEvidence(
        accepted.identity,
        accepted.occurrence_identity,
        request_identity,
        policy.identity,
        resolved_environment_identity,
        combined_operations,
        combined_failures,
        covariance_branch.residual_jacobian,
        covariance_branch.rank_diagnostic,
        covariance_branch.covariance,
        covariance_branch.marginal_errors,
        covariance_branch.correlations,
        constraint_branch.jacobian,
        constraint_branch.propagation,
        constraint_branch.marginal_errors,
        constraint_branch.correlations,
        accepted_anchor=accepted,
        requested_output_scope=output_scope,
        requested_output_units=output_units,
        requested_output_scales=output_scales,
        source_problem=problem,
        source_parameterization=parameterization,
        source_engine=engine,
        source_policy=policy,
        source_capabilities=compiled_capabilities,
    )
