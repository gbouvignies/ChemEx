"""Immutable native parameter roles, constraints, and resolved scalar values.

This module compiles active method roles for authoritative native evaluation
without mutating parameter definitions, configuration, or committed Analysis
Values.
"""

from __future__ import annotations

import ast
import hashlib
import inspect
import json
import math
import operator
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import ClassVar, Literal, cast
from uuid import uuid4

import numpy as np

from chemex.configuration.method_plan import (
    RoleAction as MethodRoleAction,
)
from chemex.configuration.methods import Method
from chemex.nmr.rates import rate_functions
from chemex.parameters.name import ParamName
from chemex.parameters.relaxation import SealedRelaxationDomains
from chemex.parameters.sealed import (
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.userfunctions import user_function_registry
from chemex.parameters.values import AnalysisValuesSnapshot

_PUBLIC_DECIMAL = re.compile(r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")
_PUBLIC_OPERATORS = frozenset({"+", "-", "*", "/", "(", ")"})


class ParameterRole(StrEnum):
    """One method-scoped role for an active parameter."""

    FIT = "fit"
    FIX = "fix"
    DERIVED = "derived"


class ParameterizationError(ValueError):
    """Base class for stable native compilation and resolution failures."""

    code: ClassVar[str] = "parameterization_error"

    def __init__(self, detail: str, **context: object) -> None:
        self.detail = detail
        self.context = MappingProxyType(dict(sorted(context.items())))
        rendered = ", ".join(
            f"{name}={value!r}" for name, value in self.context.items()
        )
        suffix = f" ({rendered})" if rendered else ""
        super().__init__(f"{self.code}: {detail}{suffix}")


class NoParameterMatchError(ParameterizationError):
    code = "no_match"


class AmbiguousParameterReferenceError(ParameterizationError):
    code = "ambiguity"


class ConstraintSelfReferenceError(ParameterizationError):
    code = "self_reference"


class ConstraintCycleError(ParameterizationError):
    code = "cycle"


class ConstraintDomainError(ParameterizationError):
    code = "domain_error"


class ConstraintEvaluationError(ParameterizationError):
    code = "evaluation_error"


class NonFiniteParameterValueError(ParameterizationError):
    code = "non_finite"


class IncompleteParameterDependenciesError(ParameterizationError):
    code = "incomplete_dependencies"


class IncompatibleParameterizationInputError(ParameterizationError):
    code = "incompatible_input"


class ConstraintProgramMismatchError(ParameterizationError):
    code = "program_mismatch"


class UnsupportedConstraintExpressionError(ParameterizationError):
    code = "unsupported_expression"


class ModelDerivationOverrideError(ParameterizationError):
    code = "model_derivation_override"


@dataclass(frozen=True, slots=True)
class ParameterDeclarationContribution:
    """One construction contribution to a parameter's scientific baseline.

    Estimation support is a scientific capability. Requiring an independent
    value and fitting it by default are separate baseline choices.
    """

    param_id: str
    supports_estimation: bool
    model_expression: str
    contributor: str
    model_owned: bool = False
    requires_independent: bool = False
    fits_by_default: bool = False
    report_only: bool = False


@dataclass(frozen=True, slots=True)
class ParameterDeclaration:
    """Sealed scientific inputs used to build each method-local baseline.

    Estimation support permits an explicit FIT override. Requiring an
    independent value chooses baseline FIT/FIX rather than DERIVED, while
    ``fits_by_default`` distinguishes FIT from FIX.
    """

    param_id: str
    supports_estimation: bool
    model_expression: str = ""
    model_owned: bool = False
    requires_independent: bool = False
    fits_by_default: bool = False
    report_only: bool = False


def baseline_parameter_role(declaration: ParameterDeclaration) -> ParameterRole:
    """Return the authoritative method-local role for one sealed declaration."""
    if declaration.model_expression and (
        declaration.model_owned or not declaration.requires_independent
    ):
        return ParameterRole.DERIVED
    if (
        declaration.requires_independent
        and declaration.fits_by_default
        and declaration.supports_estimation
    ):
        return ParameterRole.FIT
    return ParameterRole.FIX


def _fingerprint(kind: str, records: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "schema": 1, "records": records},
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class SealedParameterDeclarations(Mapping[str, ParameterDeclaration]):
    """Canonical immutable baseline roles and model-owned derivations."""

    _items: tuple[ParameterDeclaration, ...]
    _index: Mapping[str, ParameterDeclaration] = field(
        init=False,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        items = tuple(self._items)
        index = MappingProxyType({item.param_id: item for item in items})
        if len(index) != len(items):
            msg = "Sealed parameter declarations contain duplicate IDs"
            raise ValueError(msg)
        object.__setattr__(self, "_items", items)
        object.__setattr__(self, "_index", index)
        object.__setattr__(
            self,
            "identity",
            _fingerprint(
                "parameter-declarations",
                tuple(
                    (
                        item.param_id,
                        item.supports_estimation,
                        item.model_expression,
                        item.model_owned,
                        item.requires_independent,
                        item.fits_by_default,
                        item.report_only,
                    )
                    for item in items
                ),
            ),
        )

    def __iter__(self) -> Iterator[str]:
        return (item.param_id for item in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, key: str) -> ParameterDeclaration:
        return self._index[key]


def seal_parameter_declarations(
    definitions: SealedDefinitions,
    contributions: Mapping[str, Sequence[ParameterDeclarationContribution]],
) -> SealedParameterDeclarations:
    """Seal additive capabilities, baseline independence, and model derivations."""
    items: list[ParameterDeclaration] = []
    for definition in definitions:
        contributed = tuple(contributions.get(definition.param_id, ()))
        if not contributed:
            raise IncompleteParameterDependenciesError(
                "No scientific baseline was contributed for a sealed parameter",
                param_id=definition.param_id,
            )
        expressions = sorted(
            {
                item.model_expression.strip()
                for item in contributed
                if item.model_expression
            }
        )
        if len(expressions) > 1:
            raise IncompatibleParameterizationInputError(
                "Contributors disagree on a model-owned expression",
                param_id=definition.param_id,
                expressions=tuple(expressions),
            )
        supports_estimation = any(item.supports_estimation for item in contributed)
        requires_independent = any(item.requires_independent for item in contributed)
        fits_by_default = any(item.fits_by_default for item in contributed)
        report_only_values = {item.report_only for item in contributed}
        if len(report_only_values) > 1:
            raise IncompatibleParameterizationInputError(
                "Contributors disagree on report-only derivation semantics",
                param_id=definition.param_id,
            )
        report_only = report_only_values.pop()
        expression = expressions[0] if expressions else ""
        model_owned = any(
            item.model_owned and item.model_expression.strip() == expression
            for item in contributed
        )
        items.append(
            ParameterDeclaration(
                definition.param_id,
                supports_estimation,
                expression,
                model_owned,
                requires_independent,
                fits_by_default,
                report_only,
            )
        )
    return SealedParameterDeclarations(tuple(items))


@dataclass(frozen=True, slots=True)
class SealedParameterModel:
    """One immutable aggregate referencing the #582 sealed artifacts."""

    model_name: str
    model_identity: str
    definitions: SealedDefinitions
    configuration: SealedConfiguration
    declarations: SealedParameterDeclarations
    relaxation_domains: SealedRelaxationDomains = field(
        default_factory=SealedRelaxationDomains,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        definition_ids = tuple(item.param_id for item in self.definitions)
        if self.configuration.definitions_identity != self.definitions.identity:
            msg = "Sealed configuration does not belong to the definitions"
            raise ValueError(msg)
        if tuple(item.param_id for item in self.configuration) != definition_ids:
            msg = "Sealed configuration ordering does not match definitions"
            raise ValueError(msg)
        if tuple(self.declarations) != definition_ids:
            msg = "Sealed declarations ordering does not match definitions"
            raise ValueError(msg)
        object.__setattr__(
            self,
            "identity",
            _fingerprint(
                "sealed-parameter-model",
                (
                    self.model_name,
                    self.model_identity,
                    self.definitions.identity,
                    self.configuration.identity,
                    self.declarations.identity,
                    self.relaxation_domains.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class LiteralExpression:
    value: float


@dataclass(frozen=True, slots=True)
class ReferenceExpression:
    param_id: str


@dataclass(frozen=True, slots=True)
class UnaryExpression:
    operator: Literal["positive", "negative"]
    operand: ScalarExpression


@dataclass(frozen=True, slots=True)
class BinaryExpression:
    operator: Literal["add", "subtract", "multiply", "divide"]
    left: ScalarExpression
    right: ScalarExpression


@dataclass(frozen=True, slots=True)
class FunctionExpression:
    function_id: str
    arguments: tuple[ScalarExpression, ...]
    component: str | None = None


type ScalarExpression = (
    LiteralExpression
    | ReferenceExpression
    | UnaryExpression
    | BinaryExpression
    | FunctionExpression
)


@dataclass(frozen=True, slots=True)
class CompiledConstraint:
    target_id: str
    expression: ScalarExpression
    dependencies: tuple[str, ...]
    source: str
    expression_text: str


_SCIENTIFIC_FUNCTION_SEMANTICS_VERSION = 1


def _semantic_value_record(value: object) -> object:
    """Encode deterministic callable state without process-local identities."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, np.generic):
        return _semantic_value_record(value.item())
    if isinstance(value, int):
        return ("integer", str(value))
    if isinstance(value, float):
        return ("float", value.hex())
    if isinstance(value, np.ndarray):
        array = np.ascontiguousarray(value)
        return (
            "ndarray",
            array.dtype.str,
            tuple(array.shape),
            hashlib.sha256(array.tobytes()).hexdigest(),
        )
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise IncompatibleParameterizationInputError(
                "Scientific function state must use string mapping keys",
            )
        return tuple(
            (key, _semantic_value_record(item)) for key, item in sorted(value.items())
        )
    if isinstance(value, (tuple, list)):
        return tuple(_semantic_value_record(item) for item in value)
    if isinstance(value, type):
        return (
            "type",
            value.__module__,
            value.__qualname__,
            *_source_record(value),
        )
    raise IncompatibleParameterizationInputError(
        "Scientific function has unsupported mutable implementation state",
        state_type=type(value).__qualname__,
    )


def _source_record(
    implementation: type[object] | Callable[..., object],
) -> tuple[str, str]:
    try:
        source = inspect.getsource(implementation).encode()
        source_file = inspect.getsourcefile(implementation)
    except (OSError, TypeError) as error:
        raise IncompatibleParameterizationInputError(
            "Scientific function implementation source is unavailable",
            implementation=type(implementation).__qualname__,
        ) from error
    module_digest = ""
    if source_file is not None:
        try:
            module_digest = hashlib.sha256(Path(source_file).read_bytes()).hexdigest()
        except OSError as error:
            raise IncompatibleParameterizationInputError(
                "Scientific function module source is unavailable",
                source_file=source_file,
            ) from error
    return hashlib.sha256(source).hexdigest(), module_digest


def _class_runtime_record(implementation_type: type[object]) -> object:
    """Record runtime class data and Python methods used by callable instances."""
    owners: list[object] = []
    for owner in reversed(implementation_type.__mro__):
        data: list[tuple[str, object]] = []
        methods: list[tuple[str, object]] = []
        for name, value in sorted(vars(owner).items()):
            method = (
                value.__func__
                if isinstance(value, (classmethod, staticmethod))
                else value
            )
            if inspect.isfunction(method):
                methods.append((name, _scientific_function_record(method)))
            elif not name.startswith("__") and not callable(value):
                data.append((name, _semantic_value_record(value)))
        if data or methods:
            owners.append(
                (owner.__module__, owner.__qualname__, tuple(data), tuple(methods))
            )
    return tuple(owners)


def _scientific_function_record(function: Callable[..., object]) -> object:
    if function is max:
        return (
            _SCIENTIFIC_FUNCTION_SEMANTICS_VERSION,
            "chemex-scalar-maximum",
        )
    unwrapped = inspect.unwrap(function)
    if inspect.isfunction(unwrapped):
        closure = tuple(
            _semantic_value_record(cell.cell_contents)
            for cell in (unwrapped.__closure__ or ())
        )
        return (
            _SCIENTIFIC_FUNCTION_SEMANTICS_VERSION,
            "function",
            unwrapped.__module__,
            unwrapped.__qualname__,
            *_source_record(unwrapped),
            _semantic_value_record(unwrapped.__defaults__),
            _semantic_value_record(unwrapped.__kwdefaults__),
            closure,
        )
    if callable(function):
        implementation_type = type(function)
        return (
            _SCIENTIFIC_FUNCTION_SEMANTICS_VERSION,
            "callable-instance",
            implementation_type.__module__,
            implementation_type.__qualname__,
            *_source_record(implementation_type),
            _class_runtime_record(implementation_type),
            _semantic_value_record(vars(function)),
        )
    raise IncompatibleParameterizationInputError(
        "Scientific function registry entry is not callable",
        function_type=type(function).__qualname__,
    )


def scientific_callable_fingerprint(function: Callable[..., object]) -> str:
    """Return ChemEx's canonical scientific implementation fingerprint."""
    return _fingerprint(
        "scientific-callable-implementation",
        _scientific_function_record(function),
    )


@dataclass(frozen=True, slots=True)
class ScientificFunctionBinder:
    """Immutable trusted binding of the existing model-owned scalar functions."""

    model_name: str
    _functions: Mapping[str, Callable[..., object]]
    _implementation_records: Mapping[str, object] = field(
        init=False,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        functions = MappingProxyType(dict(sorted(self._functions.items())))
        object.__setattr__(self, "_functions", functions)
        implementation_records = MappingProxyType(
            {
                name: _scientific_function_record(function)
                for name, function in functions.items()
            }
        )
        object.__setattr__(self, "_implementation_records", implementation_records)
        records = tuple(implementation_records.items())
        object.__setattr__(self, "identity", _fingerprint("function-binder", records))

    def __contains__(self, function_id: object) -> bool:
        return function_id in self._functions

    def __getitem__(self, function_id: str) -> Callable[..., object]:
        return self._functions[function_id]

    def validate_implementations(self) -> None:
        """Reject mutation of a callable after its implementation was bound."""
        for function_id, function in self._functions.items():
            expected = self._implementation_records[function_id]
            try:
                actual = _scientific_function_record(function)
            except ParameterizationError as error:
                raise ConstraintProgramMismatchError(
                    "A bound scientific function no longer has its compiled semantics",
                    function_id=function_id,
                ) from error
            if actual != expected:
                raise ConstraintProgramMismatchError(
                    "A bound scientific function changed after program compilation",
                    function_id=function_id,
                )

    def worker_bindings(self) -> tuple[tuple[str, Callable[..., object]], ...]:
        """Return the validated callable bindings needed by an isolated worker."""
        self.validate_implementations()
        return tuple(self._functions.items())

    @classmethod
    def for_model(cls, model_name: str) -> ScientificFunctionBinder:
        functions: dict[str, Callable[..., object]] = {
            name: cast("Callable[..., object]", function)
            for name, function in (
                rate_functions | user_function_registry.get(model_name)
            ).items()
        }
        functions["max"] = max
        return cls(model_name, functions)


@dataclass(frozen=True, slots=True)
class IndependentValueFrame:
    """One complete ordered independent-value input for a compiled program."""

    parameterization_identity: str
    program_fingerprint: str
    occurrence_identity: str
    revision: int
    _items: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        items = tuple(self._items)
        if len({param_id for param_id, _value in items}) != len(items):
            msg = "Independent-value frame contains duplicate parameter IDs"
            raise ValueError(msg)
        object.__setattr__(self, "_items", items)

    def with_updates(self, updates: Mapping[str, float]) -> IndependentValueFrame:
        unknown = set(updates) - {param_id for param_id, _value in self._items}
        if unknown:
            raise IncompatibleParameterizationInputError(
                "Independent-value updates contain unknown IDs",
                unknown_ids=tuple(sorted(unknown)),
            )
        return type(self)(
            parameterization_identity=self.parameterization_identity,
            program_fingerprint=self.program_fingerprint,
            occurrence_identity=self.occurrence_identity,
            revision=self.revision,
            _items=tuple(
                (param_id, updates.get(param_id, value))
                for param_id, value in self._items
            ),
        )

    def ordered_items(self) -> tuple[tuple[str, float], ...]:
        """Return the complete independent frame in canonical program order."""
        return self._items


@dataclass(frozen=True, slots=True)
class ResolvedParameterValues(Mapping[str, float]):
    """Immutable complete independent and derived scalar values for active scope."""

    parameterization_identity: str
    program_fingerprint: str
    occurrence_identity: str
    revision: int
    _items: tuple[tuple[str, float], ...]
    _index: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        items = tuple(self._items)
        index = MappingProxyType(dict(items))
        if len(index) != len(items):
            msg = "Resolved parameter values contain duplicate IDs"
            raise ValueError(msg)
        object.__setattr__(self, "_items", items)
        object.__setattr__(self, "_index", index)

    def __iter__(self) -> Iterator[str]:
        return (param_id for param_id, _value in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, key: str) -> float:
        return self._index[key]


@dataclass(frozen=True, slots=True)
class ConstraintProgram:
    """Immutable restricted value program over one deterministic active scope."""

    parameter_model_identity: str
    model_identity: str
    definitions_identity: str
    configuration_identity: str
    function_binder_identity: str
    scope_ids: tuple[str, ...]
    independent_ids: tuple[str, ...]
    derived_ids: tuple[str, ...]
    constraints: tuple[CompiledConstraint, ...]
    evaluation_order: tuple[str, ...]
    relaxation_domains: SealedRelaxationDomains = field(
        default_factory=SealedRelaxationDomains,
    )
    fingerprint: str = field(init=False)

    def __post_init__(self) -> None:
        records = (
            self.parameter_model_identity,
            self.model_identity,
            self.definitions_identity,
            self.configuration_identity,
            self.function_binder_identity,
            self.scope_ids,
            self.independent_ids,
            self.derived_ids,
            tuple(
                (
                    item.target_id,
                    _expression_record(item.expression),
                    item.dependencies,
                )
                for item in self.constraints
            ),
            self.evaluation_order,
            self.relaxation_domains.identity,
        )
        object.__setattr__(
            self, "fingerprint", _fingerprint("constraint-program", records)
        )


@dataclass(frozen=True, slots=True)
class ActiveParameterization:
    """Fresh method-scoped roles plus a restricted immutable value program."""

    program: ConstraintProgram
    binder: ScientificFunctionBinder = field(repr=False, compare=False)
    occurrence_identity: str
    source_revision: int
    _roles: tuple[tuple[str, ParameterRole], ...]
    identity: str = field(init=False)
    evaluator_identity: str = field(init=False)
    _role_index: Mapping[str, ParameterRole] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _ordered_constraints: tuple[CompiledConstraint, ...] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _affected_constraint_positions: Mapping[str, tuple[int, ...]] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        # Validate the exact callable references once, when this immutable
        # parameterization occurrence becomes executable. Repeated resolution
        # then uses the already-bound functions without source introspection.
        self.binder.validate_implementations()
        roles = tuple(self._roles)
        role_index = MappingProxyType(dict(roles))
        if tuple(role_index) != self.program.scope_ids:
            msg = "Parameter roles do not match the constraint-program scope"
            raise ValueError(msg)
        object.__setattr__(self, "_roles", roles)
        object.__setattr__(self, "_role_index", role_index)
        constraints = {item.target_id: item for item in self.program.constraints}
        ordered_constraints = tuple(
            constraints[target_id] for target_id in self.program.evaluation_order
        )
        outgoing: dict[str, list[str]] = {}
        for constraint in ordered_constraints:
            for dependency in constraint.dependencies:
                outgoing.setdefault(dependency, []).append(constraint.target_id)
        positions = {
            constraint.target_id: position
            for position, constraint in enumerate(ordered_constraints)
        }
        affected: dict[str, tuple[int, ...]] = {}
        for independent_id in self.program.independent_ids:
            pending = list(outgoing.get(independent_id, ()))
            targets: set[str] = set()
            while pending:
                target_id = pending.pop()
                if target_id in targets:
                    continue
                targets.add(target_id)
                pending.extend(outgoing.get(target_id, ()))
            affected[independent_id] = tuple(
                sorted(positions[target_id] for target_id in targets)
            )
        object.__setattr__(self, "_ordered_constraints", ordered_constraints)
        object.__setattr__(
            self,
            "_affected_constraint_positions",
            MappingProxyType(affected),
        )
        object.__setattr__(
            self,
            "identity",
            _fingerprint(
                "active-parameterization",
                (
                    self.program.fingerprint,
                    self.occurrence_identity,
                    self.source_revision,
                    tuple((param_id, role.value) for param_id, role in roles),
                ),
            ),
        )
        # The occurrence and source revision protect the lifecycle boundary.
        # Evaluation is deliberately downstream of that check: two valid
        # occurrences with the same compiled scientific program must have the
        # same evaluator-facing identity.
        object.__setattr__(
            self,
            "evaluator_identity",
            _fingerprint(
                "evaluator-parameterization",
                (
                    self.program.fingerprint,
                    tuple((param_id, role.value) for param_id, role in roles),
                ),
            ),
        )

    @property
    def independent_ids(self) -> tuple[str, ...]:
        return self.program.independent_ids

    @property
    def derived_ids(self) -> tuple[str, ...]:
        return self.program.derived_ids

    @property
    def scope_ids(self) -> tuple[str, ...]:
        return self.program.scope_ids

    def role(self, param_id: str) -> ParameterRole:
        return self._role_index[param_id]

    def frame_from_snapshot(
        self,
        snapshot: AnalysisValuesSnapshot,
    ) -> IndependentValueFrame:
        expected = (
            self.occurrence_identity,
            self.source_revision,
            self.program.model_identity,
            self.program.definitions_identity,
            self.program.configuration_identity,
        )
        actual = (
            snapshot.occurrence_identity,
            snapshot.revision,
            snapshot.model_identity,
            snapshot.definitions_identity,
            snapshot.configuration_identity,
        )
        if actual != expected:
            raise IncompatibleParameterizationInputError(
                "Analysis Values snapshot is incompatible with the parameterization",
                expected=expected,
                actual=actual,
            )
        try:
            items = tuple(
                (param_id, snapshot[param_id]) for param_id in self.independent_ids
            )
        except KeyError as error:
            raise IncompleteParameterDependenciesError(
                "Analysis Values snapshot lacks an independent parameter",
                param_id=str(error.args[0]),
            ) from error
        return IndependentValueFrame(
            parameterization_identity=self.identity,
            program_fingerprint=self.program.fingerprint,
            occurrence_identity=self.occurrence_identity,
            revision=self.source_revision,
            _items=items,
        )

    def _resolve_values(
        self,
        frame: IndependentValueFrame,
        previous: Mapping[str, float] | None = None,
    ) -> dict[str, float]:
        _validate_frame(self, frame)
        independent_values = {
            param_id: _finite_scalar(value, param_id=param_id)
            for param_id, value in frame._items
        }
        if previous is None:
            values = independent_values
            constraint_positions = range(len(self._ordered_constraints))
        else:
            values = dict(previous)
            changed_ids = {
                param_id
                for param_id, value in independent_values.items()
                if value != previous[param_id]
            }
            values.update(independent_values)
            constraint_positions = sorted(
                {
                    position
                    for param_id in changed_ids
                    for position in self._affected_constraint_positions[param_id]
                }
            )
        for position in constraint_positions:
            constraint = self._ordered_constraints[position]
            values[constraint.target_id] = _evaluate_expression(
                constraint.expression,
                values,
                self.binder,
                target_id=constraint.target_id,
                position=position,
            )
        return values

    def resolve(self, frame: IndependentValueFrame) -> ResolvedParameterValues:
        values = self._resolve_values(frame)
        return ResolvedParameterValues(
            parameterization_identity=self.identity,
            program_fingerprint=self.program.fingerprint,
            occurrence_identity=self.occurrence_identity,
            revision=self.source_revision,
            _items=tuple((param_id, values[param_id]) for param_id in self.scope_ids),
        )


@dataclass(frozen=True, slots=True)
class StaticParameterization:
    """A value-independent parameter program for one Method Step scope."""

    program: ConstraintProgram
    binder: ScientificFunctionBinder = field(repr=False, compare=False)
    roles: tuple[tuple[str, ParameterRole], ...]

    @property
    def fit_ids(self) -> tuple[str, ...]:
        return tuple(
            param_id for param_id, role in self.roles if role is ParameterRole.FIT
        )

    def bind(self, snapshot: AnalysisValuesSnapshot) -> ActiveParameterization:
        expected = (
            self.program.model_identity,
            self.program.definitions_identity,
            self.program.configuration_identity,
        )
        actual = (
            snapshot.model_identity,
            snapshot.definitions_identity,
            snapshot.configuration_identity,
        )
        if actual != expected:
            raise IncompatibleParameterizationInputError(
                "Analysis Values snapshot does not belong to the compiled Method Step",
                expected=expected,
                actual=actual,
            )
        return ActiveParameterization(
            self.program,
            self.binder,
            snapshot.occurrence_identity,
            snapshot.revision,
            self.roles,
        )


@dataclass(frozen=True, slots=True)
class ReportableParameterSet:
    """The finite derived values selected once for uncertainty and output."""

    parameterization: ActiveParameterization
    values: ResolvedParameterValues
    report_only_ids: tuple[str, ...]

    @property
    def report_only_values(self) -> Mapping[str, float]:
        return MappingProxyType(
            {param_id: self.values[param_id] for param_id in self.report_only_ids}
        )


def compatible_reference_context(
    candidate: ParamDefinition,
    target: ParamDefinition,
    selector: ParamName,
) -> tuple[int, frozenset[str], int] | None:
    """Return spin specificity, matching condition fields, and extra context."""
    matched: set[str] = set()
    spin_specificity = 0
    extras = 0
    candidate_spin = SpinSystem.from_name(candidate.spin_system_name)
    target_spin = SpinSystem.from_name(target.spin_system_name)
    if candidate_spin and not selector.spin_system:
        if target_spin and (
            candidate_spin == target_spin
            or candidate_spin.match(target_spin)
            or target_spin.match(candidate_spin)
        ):
            spin_specificity = 2
        elif target_spin and candidate_spin.shares_group_site(target_spin):
            spin_specificity = 1
        elif target_spin:
            return None
        else:
            extras += 1

    selector_conditions = selector.conditions.model_dump()
    target_conditions = dict(target.condition_entries)
    for name, value in candidate.condition_entries:
        if selector_conditions.get(name) is not None:
            continue
        target_value = target_conditions.get(name)
        if target_value is None:
            extras += 1
        elif target_value == value:
            matched.add(name)
        else:
            return None
    return spin_specificity, frozenset(matched), extras


@dataclass(frozen=True, slots=True)
class _ExpressionCompileContext:
    definitions: SealedDefinitions
    binder: ScientificFunctionBinder
    target_id: str


def _compile_literal(node: ast.Constant, target_id: str) -> LiteralExpression:
    if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
        raise UnsupportedConstraintExpressionError(
            "Only finite numeric scalar literals are supported",
            target_id=target_id,
        )
    try:
        value = float(node.value)
    except OverflowError as error:
        raise NonFiniteParameterValueError(
            "Constraint literal exceeds the supported finite scalar range",
            target_id=target_id,
        ) from error
    if not math.isfinite(value):
        raise NonFiniteParameterValueError(
            "Constraint literal is non-finite",
            target_id=target_id,
            value=value,
        )
    return LiteralExpression(0.0 if value == 0.0 else value)


def _compile_reference(
    node: ast.Name,
    context: _ExpressionCompileContext,
) -> ReferenceExpression:
    if not node.id.startswith("__"):
        raise UnsupportedConstraintExpressionError(
            "Bare names are not supported in scalar constraints",
            target_id=context.target_id,
            name=node.id,
        )
    param_id = node.id
    if param_id not in context.definitions:
        raise IncompleteParameterDependenciesError(
            "Model expression references an unknown sealed parameter",
            target_id=context.target_id,
            dependency_id=param_id,
        )
    if param_id == context.target_id:
        raise ConstraintSelfReferenceError(
            "Model expression directly references its target",
            target_id=context.target_id,
        )
    return ReferenceExpression(param_id)


def _compile_unary(
    node: ast.UnaryOp,
    context: _ExpressionCompileContext,
) -> UnaryExpression:
    unary: Literal["positive", "negative"] = (
        "positive" if isinstance(node.op, ast.UAdd) else "negative"
    )
    return UnaryExpression(unary, _compile_ast(node.operand, context))


def _compile_binary(
    node: ast.BinOp,
    context: _ExpressionCompileContext,
) -> BinaryExpression:
    operators: dict[
        type[ast.operator],
        Literal["add", "subtract", "multiply", "divide"],
    ] = {
        ast.Add: "add",
        ast.Sub: "subtract",
        ast.Mult: "multiply",
        ast.Div: "divide",
    }
    return BinaryExpression(
        operators[type(node.op)],
        _compile_ast(node.left, context),
        _compile_ast(node.right, context),
    )


def _compile_function(
    node: ast.Call,
    context: _ExpressionCompileContext,
) -> FunctionExpression:
    if (
        not isinstance(node.func, ast.Name)
        or node.keywords
        or node.func.id not in context.binder
    ):
        function_id = node.func.id if isinstance(node.func, ast.Name) else ""
        raise UnsupportedConstraintExpressionError(
            "Model expression requests an unsupported scientific function",
            target_id=context.target_id,
            function_id=function_id,
        )
    if node.func.id != "max" or len(node.args) != 2:
        raise UnsupportedConstraintExpressionError(
            "Only the two-argument scalar maximum may be used without a component",
            target_id=context.target_id,
            function_id=node.func.id,
        )
    return FunctionExpression(
        node.func.id,
        tuple(_compile_ast(argument, context) for argument in node.args),
    )


def _compile_function_component(
    node: ast.Subscript,
    context: _ExpressionCompileContext,
) -> FunctionExpression:
    if not isinstance(node.value, ast.Call):
        raise UnsupportedConstraintExpressionError(
            "Only model-owned scientific-function components may be selected",
            target_id=context.target_id,
        )
    call = node.value
    if (
        not isinstance(call.func, ast.Name)
        or call.keywords
        or call.func.id == "max"
        or call.func.id not in context.binder
    ):
        function_id = call.func.id if isinstance(call.func, ast.Name) else ""
        raise UnsupportedConstraintExpressionError(
            "Model expression requests an unsupported component function",
            target_id=context.target_id,
            function_id=function_id,
        )
    component_node = node.slice
    if not isinstance(component_node, ast.Constant) or not isinstance(
        component_node.value,
        str,
    ):
        raise UnsupportedConstraintExpressionError(
            "Scientific-function component must be a literal string",
            target_id=context.target_id,
        )
    return FunctionExpression(
        call.func.id,
        tuple(_compile_ast(argument, context) for argument in call.args),
        component_node.value,
    )


def _compile_ast(
    node: ast.AST,
    context: _ExpressionCompileContext,
) -> ScalarExpression:
    if isinstance(node, ast.Constant):
        return _compile_literal(node, context.target_id)
    if isinstance(node, ast.Name):
        return _compile_reference(node, context)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        return _compile_unary(node, context)
    if isinstance(node, ast.BinOp) and isinstance(
        node.op,
        (ast.Add, ast.Sub, ast.Mult, ast.Div),
    ):
        return _compile_binary(node, context)
    if isinstance(node, ast.Call):
        return _compile_function(node, context)
    if isinstance(node, ast.Subscript):
        return _compile_function_component(node, context)
    raise UnsupportedConstraintExpressionError(
        "Expression contains syntax outside ChemEx scalar constraint semantics",
        target_id=context.target_id,
        syntax=type(node).__name__,
    )


def _parse_expression(
    text: str,
    *,
    definitions: SealedDefinitions,
    binder: ScientificFunctionBinder,
    target_id: str,
) -> ScalarExpression:
    try:
        parsed = ast.parse(text, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise UnsupportedConstraintExpressionError(
            "Constraint is not a valid scalar expression",
            target_id=target_id,
            expression=text,
        ) from error
    return _compile_ast(
        parsed.body,
        _ExpressionCompileContext(
            definitions,
            binder,
            target_id,
        ),
    )


def _dependencies(expression: ScalarExpression) -> tuple[str, ...]:
    ordered: dict[str, None] = {}

    def visit(node: ScalarExpression) -> None:
        if isinstance(node, ReferenceExpression):
            ordered.setdefault(node.param_id, None)
        elif isinstance(node, UnaryExpression):
            visit(node.operand)
        elif isinstance(node, BinaryExpression):
            visit(node.left)
            visit(node.right)
        elif isinstance(node, FunctionExpression):
            for argument in node.arguments:
                visit(argument)

    visit(expression)
    return tuple(ordered)


def _expression_record(expression: ScalarExpression) -> object:
    if isinstance(expression, LiteralExpression):
        return ("literal", float(expression.value).hex())
    if isinstance(expression, ReferenceExpression):
        return ("reference", expression.param_id)
    if isinstance(expression, UnaryExpression):
        return ("unary", expression.operator, _expression_record(expression.operand))
    if isinstance(expression, BinaryExpression):
        return (
            "binary",
            expression.operator,
            _expression_record(expression.left),
            _expression_record(expression.right),
        )
    return (
        "function",
        expression.function_id,
        expression.component,
        tuple(_expression_record(item) for item in expression.arguments),
    )


def _topological_order(
    derived_ids: tuple[str, ...],
    constraints: Mapping[str, CompiledConstraint],
    order: Mapping[str, int],
) -> tuple[str, ...]:
    derived = set(derived_ids)
    incoming = {
        param_id: sum(
            dependency in derived for dependency in constraints[param_id].dependencies
        )
        for param_id in derived_ids
    }
    outgoing: dict[str, list[str]] = {param_id: [] for param_id in derived_ids}
    for target_id in derived_ids:
        for dependency in constraints[target_id].dependencies:
            if dependency in derived:
                outgoing[dependency].append(target_id)
    ready = sorted(
        (param_id for param_id, count in incoming.items() if count == 0),
        key=order.__getitem__,
    )
    result: list[str] = []
    while ready:
        param_id = ready.pop(0)
        result.append(param_id)
        for dependent in sorted(outgoing[param_id], key=order.__getitem__):
            incoming[dependent] -= 1
            if incoming[dependent] == 0:
                ready.append(dependent)
                ready.sort(key=order.__getitem__)
    if len(result) != len(derived_ids):
        cycle_ids = _find_constraint_cycle(derived_ids, constraints, order)
        raise ConstraintCycleError(
            "Constraint dependency graph contains a cycle",
            param_ids=cycle_ids,
            constraints=tuple(
                (
                    param_id,
                    constraints[param_id].source,
                    constraints[param_id].expression_text,
                )
                for param_id in cycle_ids
            ),
        )
    return tuple(result)


def _find_constraint_cycle(
    derived_ids: tuple[str, ...],
    constraints: Mapping[str, CompiledConstraint],
    order: Mapping[str, int],
) -> tuple[str, ...]:
    """Return the first exact dependency cycle in stable definition order."""
    derived = set(derived_ids)
    state: dict[str, int] = dict.fromkeys(derived_ids, 0)
    stack: list[str] = []
    positions: dict[str, int] = {}

    def visit(param_id: str) -> tuple[str, ...] | None:
        state[param_id] = 1
        positions[param_id] = len(stack)
        stack.append(param_id)
        dependencies = sorted(
            (
                dependency
                for dependency in constraints[param_id].dependencies
                if dependency in derived
            ),
            key=order.__getitem__,
        )
        for dependency in dependencies:
            if state[dependency] == 0:
                if cycle := visit(dependency):
                    return cycle
            elif state[dependency] == 1:
                return tuple(stack[positions[dependency] :])
        stack.pop()
        positions.pop(param_id)
        state[param_id] = 2
        return None

    for param_id in derived_ids:
        if state[param_id] == 0 and (cycle := visit(param_id)):
            return cycle
    raise AssertionError("Unresolved constraint graph did not contain a cycle")


def _validate_parameterization_inputs(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    required_ids: Sequence[str] | set[str],
) -> set[str]:
    expected = (
        parameter_model.model_identity,
        parameter_model.definitions.identity,
        parameter_model.configuration.identity,
    )
    actual = (
        snapshot.model_identity,
        snapshot.definitions_identity,
        snapshot.configuration_identity,
    )
    if actual != expected:
        raise IncompatibleParameterizationInputError(
            "Analysis Values snapshot does not belong to the sealed parameter model",
            expected=expected,
            actual=actual,
        )
    required = set(required_ids)
    unknown_required = required - set(parameter_model.declarations)
    if unknown_required:
        raise IncompleteParameterDependenciesError(
            "Required scope contains unknown sealed parameter IDs",
            param_ids=tuple(sorted(unknown_required)),
        )
    if not required:
        raise IncompleteParameterDependenciesError("Required parameter scope is empty")
    return required


def _static_parameterization_from_scope(
    parameter_model: SealedParameterModel,
    binder: ScientificFunctionBinder,
    active: set[str],
    roles_by_id: Mapping[str, ParameterRole],
    compiled: Mapping[str, CompiledConstraint],
) -> StaticParameterization:
    definitions = parameter_model.definitions
    definition_order = {
        definition.param_id: position for position, definition in enumerate(definitions)
    }
    scope_ids = tuple(
        definition.param_id
        for definition in definitions
        if definition.param_id in active
    )
    roles = tuple((param_id, roles_by_id[param_id]) for param_id in scope_ids)
    independent_ids = tuple(
        param_id
        for param_id, role in roles
        if role in (ParameterRole.FIT, ParameterRole.FIX)
    )
    derived_ids = tuple(
        param_id for param_id, role in roles if role is ParameterRole.DERIVED
    )
    evaluation_order = _topological_order(derived_ids, compiled, definition_order)
    program = ConstraintProgram(
        parameter_model_identity=parameter_model.identity,
        model_identity=parameter_model.model_identity,
        definitions_identity=definitions.identity,
        configuration_identity=parameter_model.configuration.identity,
        function_binder_identity=binder.identity,
        scope_ids=scope_ids,
        independent_ids=independent_ids,
        derived_ids=derived_ids,
        constraints=tuple(compiled[param_id] for param_id in derived_ids),
        evaluation_order=evaluation_order,
        relaxation_domains=SealedRelaxationDomains(
            tuple(
                block
                for block in parameter_model.relaxation_domains.blocks
                if {*block.diagonal_ids, *(item[2] for item in block.off_diagonal_ids)}
                <= active
            )
        ),
    )
    return StaticParameterization(program, binder, roles)


def compile_static_parameterization(  # noqa: C901 - dependency closure for one resolved scope
    parameter_model: SealedParameterModel,
    roles: Mapping[str, ParameterRole],
    method_constraints: Mapping[str, CompiledConstraint],
    required_ids: Sequence[str] | set[str],
) -> StaticParameterization:
    """Project resolved roles and constraints without reading Analysis Values."""
    required = set(required_ids)
    unknown = required - set(parameter_model.declarations)
    if unknown:
        raise IncompleteParameterDependenciesError(
            "Required scope contains unknown sealed parameter IDs",
            param_ids=tuple(sorted(unknown)),
        )
    if not required:
        raise IncompleteParameterDependenciesError("Required parameter scope is empty")
    binder = ScientificFunctionBinder.for_model(parameter_model.model_name)
    active = set(required)
    compiled: dict[str, CompiledConstraint] = {}
    pending = True
    while pending:
        pending = False
        for definition in parameter_model.definitions:
            param_id = definition.param_id
            if param_id not in active or param_id in compiled:
                continue
            if roles[param_id] is not ParameterRole.DERIVED:
                continue
            constraint = method_constraints.get(param_id)
            if constraint is None:
                expression_text = parameter_model.declarations[
                    param_id
                ].model_expression
                expression = _parse_expression(
                    expression_text,
                    definitions=parameter_model.definitions,
                    binder=binder,
                    target_id=param_id,
                )
                constraint = CompiledConstraint(
                    param_id,
                    expression,
                    _dependencies(expression),
                    "model"
                    if parameter_model.declarations[param_id].model_owned
                    else "baseline",
                    expression_text,
                )
            compiled[param_id] = constraint
            for dependency in constraint.dependencies:
                if dependency not in parameter_model.declarations:
                    raise IncompleteParameterDependenciesError(
                        "Constraint dependency is absent from the sealed model",
                        target_id=param_id,
                        dependency_id=dependency,
                    )
                if dependency not in active:
                    active.add(dependency)
                    pending = True
    return _static_parameterization_from_scope(
        parameter_model, binder, active, roles, compiled
    )


def compile_active_parameterization(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    method: Method,
    required_ids: Sequence[str] | set[str],
) -> ActiveParameterization:
    """Adapt the legacy Python Method through the Method semantic compiler."""
    from chemex.configuration.method_input import normalize_method_plan
    from chemex.configuration.method_plan import MethodFormatError
    from chemex.configuration.method_validation import resolve_method_plan

    _validate_parameterization_inputs(parameter_model, snapshot, required_ids)
    try:
        plan = normalize_method_plan({"DEFAULT": method})
        resolved = resolve_method_plan(plan, parameter_model)[0]
    except MethodFormatError as error:
        context = dict(error.detail_context)
        match error.detail_code:
            case "model_derivation_override":
                raise ModelDerivationOverrideError(error.message, **context) from error
            case "incompatible_input":
                raise IncompatibleParameterizationInputError(
                    error.message, **context
                ) from error
            case "no_match":
                raise NoParameterMatchError(error.message, **context) from error
            case "self_reference":
                raise ConstraintSelfReferenceError(error.message, **context) from error
            case "ambiguity":
                raise AmbiguousParameterReferenceError(
                    error.message, **context
                ) from error
            case "cycle":
                constraints = cast(
                    tuple[tuple[str, str, str], ...],
                    context.get("constraints", ()),
                )
                context["constraints"] = tuple(
                    (
                        param_id,
                        source,
                        method.constraints[int(source.rsplit(":", 1)[1])]
                        .split("=", maxsplit=1)[1]
                        .strip(),
                    )
                    for param_id, source, _text in constraints
                )
                raise ConstraintCycleError(error.message, **context) from error
            case "non_finite":
                raise NonFiniteParameterValueError(error.message) from error
            case _:
                raise UnsupportedConstraintExpressionError(error.message) from error
    static = compile_static_parameterization(
        parameter_model, resolved.roles, resolved.constraints, required_ids
    )
    return static.bind(snapshot)


def compile_active_parameterization_from_actions(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    actions: Sequence[MethodRoleAction],
    required_ids: Sequence[str] | set[str],
) -> ActiveParameterization:
    """Compatibility preview through the authoritative Method resolver."""
    from chemex.configuration.method_plan import FormatOrigin, MethodPlan, StepPlan
    from chemex.configuration.method_validation import resolve_method_plan

    _validate_parameterization_inputs(parameter_model, snapshot, required_ids)
    plan = MethodPlan(
        FormatOrigin.V2,
        (StepPlan("DEFAULT", role_actions=tuple(actions)),),
    )
    resolved = resolve_method_plan(plan, parameter_model)[0]
    return compile_static_parameterization(
        parameter_model, resolved.roles, resolved.constraints, required_ids
    ).bind(snapshot)


def extend_parameterization_for_report_only_outputs(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    parameterization: ActiveParameterization,
    report_only_ids: Sequence[str],
) -> ActiveParameterization:
    """Add finite model-owned report constraints without changing fit roles."""
    requested = tuple(dict.fromkeys(report_only_ids))
    if not requested:
        return ActiveParameterization(
            parameterization.program,
            parameterization.binder,
            snapshot.occurrence_identity,
            snapshot.revision,
            tuple(
                (param_id, parameterization.role(param_id))
                for param_id in parameterization.scope_ids
            ),
        )
    if (
        parameterization.program.parameter_model_identity != parameter_model.identity
        or snapshot.model_identity != parameter_model.model_identity
        or snapshot.definitions_identity != parameter_model.definitions.identity
        or snapshot.configuration_identity != parameter_model.configuration.identity
    ):
        raise IncompatibleParameterizationInputError(
            "Report outputs do not belong to the active parameter model"
        )

    active = set(parameterization.scope_ids)
    constraints = {
        constraint.target_id: constraint
        for constraint in parameterization.program.constraints
    }
    for param_id in requested:
        declaration = parameter_model.declarations[param_id]
        if not declaration.report_only or not declaration.model_expression:
            raise ParameterizationError(
                "Report extension accepts only model-owned report-only derivations",
                param_id=param_id,
            )
        expression = _parse_expression(
            declaration.model_expression,
            definitions=parameter_model.definitions,
            binder=parameterization.binder,
            target_id=param_id,
        )
        dependencies = _dependencies(expression)
        missing = set(dependencies) - active
        if missing:
            raise IncompleteParameterDependenciesError(
                "Report-only derivation depends on the inactive parameter scope",
                target_id=param_id,
                param_ids=tuple(sorted(missing)),
            )
        constraints[param_id] = CompiledConstraint(
            param_id,
            expression,
            dependencies,
            "model",
            declaration.model_expression,
        )
        active.add(param_id)

    definition_order = {
        definition.param_id: position
        for position, definition in enumerate(parameter_model.definitions)
    }
    scope_ids = tuple(
        definition.param_id
        for definition in parameter_model.definitions
        if definition.param_id in active
    )
    roles = tuple(
        (
            param_id,
            ParameterRole.DERIVED
            if param_id in requested
            else parameterization.role(param_id),
        )
        for param_id in scope_ids
    )
    independent_ids = tuple(
        param_id
        for param_id, role in roles
        if role in (ParameterRole.FIT, ParameterRole.FIX)
    )
    derived_ids = tuple(
        param_id for param_id, role in roles if role is ParameterRole.DERIVED
    )
    evaluation_order = _topological_order(
        derived_ids,
        constraints,
        definition_order,
    )
    program = ConstraintProgram(
        parameter_model_identity=parameter_model.identity,
        model_identity=parameter_model.model_identity,
        definitions_identity=snapshot.definitions_identity,
        configuration_identity=snapshot.configuration_identity,
        function_binder_identity=parameterization.binder.identity,
        scope_ids=scope_ids,
        independent_ids=independent_ids,
        derived_ids=derived_ids,
        constraints=tuple(constraints[param_id] for param_id in derived_ids),
        evaluation_order=evaluation_order,
        relaxation_domains=parameterization.program.relaxation_domains,
    )
    return ActiveParameterization(
        program,
        parameterization.binder,
        snapshot.occurrence_identity,
        snapshot.revision,
        roles,
    )


def build_initial_analysis_values(
    parameter_model: SealedParameterModel,
) -> Mapping[str, float]:
    """Natively fill missing model-derived revision-zero configuration values."""
    configuration = parameter_model.configuration
    deferred_ids = deferred_derived_ids(parameter_model)
    missing_ids = tuple(
        config.param_id
        for config in configuration
        if config.effective_value is None and config.param_id not in deferred_ids
    )
    for param_id in missing_ids:
        if not parameter_model.declarations[param_id].model_expression:
            raise IncompleteParameterDependenciesError(
                "A configured parameter lacks both a value and model derivation",
                param_id=param_id,
            )
    bootstrap_snapshot = AnalysisValuesSnapshot(
        occurrence_identity=f"bootstrap:{uuid4().hex}",
        model_identity=parameter_model.model_identity,
        definitions_identity=parameter_model.definitions.identity,
        configuration_identity=configuration.identity,
        revision=0,
        _items=tuple(
            (config.param_id, config.effective_value)
            for config in configuration
            if config.effective_value is not None
            and config.param_id not in deferred_ids
        ),
    )
    roles = {
        param_id: (
            ParameterRole.DERIVED
            if declaration.model_expression
            else baseline_parameter_role(declaration)
        )
        for param_id, declaration in parameter_model.declarations.items()
    }
    parameterization = compile_static_parameterization(
        parameter_model,
        roles,
        {},
        set(parameter_model.declarations) - set(deferred_ids),
    ).bind(bootstrap_snapshot)
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(bootstrap_snapshot)
    )
    return MappingProxyType(
        {
            config.param_id: (
                resolved[config.param_id]
                if config.effective_value is None
                else config.effective_value
            )
            for config in configuration
            if config.param_id not in deferred_ids
        }
    )


def deferred_derived_ids(
    parameter_model: SealedParameterModel,
) -> tuple[str, ...]:
    """Keep public report-only derivations out of central materialization."""
    return report_only_derived_ids(parameter_model)


def report_only_derived_ids(
    parameter_model: SealedParameterModel,
) -> tuple[str, ...]:
    """Return model-owned outputs resolved only when explicitly needed."""
    return tuple(
        declaration.param_id
        for declaration in parameter_model.declarations.values()
        if declaration.report_only
    )


def _validate_frame(
    parameterization: ActiveParameterization,
    frame: IndependentValueFrame,
) -> None:
    if frame.program_fingerprint != parameterization.program.fingerprint:
        raise ConstraintProgramMismatchError(
            "Independent-value frame names a different constraint program",
            expected=parameterization.program.fingerprint,
            actual=frame.program_fingerprint,
        )
    if (
        frame.parameterization_identity != parameterization.identity
        or frame.occurrence_identity != parameterization.occurrence_identity
        or frame.revision != parameterization.source_revision
    ):
        raise IncompatibleParameterizationInputError(
            "Independent-value frame belongs to another parameterization occurrence",
            expected_identity=parameterization.identity,
            actual_identity=frame.parameterization_identity,
            expected_occurrence=parameterization.occurrence_identity,
            actual_occurrence=frame.occurrence_identity,
        )
    frame_ids = tuple(param_id for param_id, _value in frame._items)
    if frame_ids != parameterization.independent_ids:
        raise IncompleteParameterDependenciesError(
            "Independent-value frame is incomplete or incorrectly ordered",
            expected_ids=parameterization.independent_ids,
            actual_ids=frame_ids,
        )


def _finite_scalar(value: object, *, param_id: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ConstraintEvaluationError(
            "Constraint value is not a real scalar",
            param_id=param_id,
            value=value,
        )
    try:
        scalar = float(value)
    except OverflowError as error:
        raise NonFiniteParameterValueError(
            "Constraint value exceeds the supported finite scalar range",
            param_id=param_id,
        ) from error
    if not math.isfinite(scalar):
        raise NonFiniteParameterValueError(
            "Constraint value is non-finite",
            param_id=param_id,
            value=scalar,
        )
    return 0.0 if scalar == 0.0 else scalar


@dataclass(frozen=True, slots=True)
class _EvaluationContext:
    values: Mapping[str, float]
    binder: ScientificFunctionBinder
    target_id: str
    position: int


def _evaluate_reference(
    expression: ReferenceExpression,
    context: _EvaluationContext,
) -> float:
    try:
        return context.values[expression.param_id]
    except KeyError as error:
        raise IncompleteParameterDependenciesError(
            "A dependency was unavailable during constraint evaluation",
            target_id=context.target_id,
            dependency_id=expression.param_id,
            position=context.position,
        ) from error


def _evaluate_binary(
    expression: BinaryExpression,
    context: _EvaluationContext,
) -> float:
    left = cast("float", _evaluate_node(expression.left, context))
    right = cast("float", _evaluate_node(expression.right, context))
    operations: dict[str, Callable[[float, float], float]] = {
        "add": operator.add,
        "subtract": operator.sub,
        "multiply": operator.mul,
        "divide": operator.truediv,
    }
    try:
        return operations[expression.operator](left, right)
    except (ArithmeticError, ValueError) as error:
        raise ConstraintDomainError(
            "Arithmetic constraint operation is outside its domain",
            target_id=context.target_id,
            operator=expression.operator,
            left=left,
            right=right,
            position=context.position,
        ) from error
    except TypeError as error:
        raise ConstraintEvaluationError(
            "Arithmetic constraint operands are not compatible scalars",
            target_id=context.target_id,
            operator=expression.operator,
            left=left,
            right=right,
            position=context.position,
        ) from error


def _evaluate_function(
    expression: FunctionExpression,
    context: _EvaluationContext,
) -> float:
    arguments = tuple(
        _evaluate_node(argument, context) for argument in expression.arguments
    )
    try:
        result = context.binder[expression.function_id](*arguments)
    except (ArithmeticError, FloatingPointError, ValueError) as error:
        raise ConstraintDomainError(
            "Scientific constraint function failed its value domain",
            target_id=context.target_id,
            function_id=expression.function_id,
            component=expression.component,
            arguments=arguments,
            position=context.position,
        ) from error
    except Exception as error:
        raise ConstraintEvaluationError(
            "Scientific constraint function could not be evaluated",
            target_id=context.target_id,
            function_id=expression.function_id,
            component=expression.component,
            arguments=arguments,
            position=context.position,
        ) from error
    if expression.component is None:
        return _finite_scalar(result, param_id=context.target_id)
    if not isinstance(result, Mapping) or expression.component not in result:
        raise ConstraintDomainError(
            "Scientific constraint function did not return its declared component",
            target_id=context.target_id,
            function_id=expression.function_id,
            component=expression.component,
            arguments=arguments,
            position=context.position,
        )
    return _finite_scalar(
        cast("Mapping[str, object]", result)[expression.component],
        param_id=context.target_id,
    )


def _evaluate_node(
    expression: ScalarExpression,
    context: _EvaluationContext,
) -> object:
    if isinstance(expression, LiteralExpression):
        return expression.value
    if isinstance(expression, ReferenceExpression):
        return _evaluate_reference(expression, context)
    if isinstance(expression, UnaryExpression):
        operand = cast("float", _evaluate_node(expression.operand, context))
        try:
            return +operand if expression.operator == "positive" else -operand
        except TypeError as error:
            raise ConstraintEvaluationError(
                "Unary constraint operand is not a compatible scalar",
                target_id=context.target_id,
                operator=expression.operator,
                operand=operand,
                position=context.position,
            ) from error
    if isinstance(expression, BinaryExpression):
        return _evaluate_binary(expression, context)
    return _evaluate_function(expression, context)


def _evaluate_expression(
    expression: ScalarExpression,
    values: Mapping[str, float],
    binder: ScientificFunctionBinder,
    *,
    target_id: str,
    position: int,
) -> float:
    with np.errstate(all="raise"):
        result = _evaluate_node(
            expression,
            _EvaluationContext(values, binder, target_id, position),
        )
    return _finite_scalar(result, param_id=target_id)
