"""Immutable native parameter roles, constraints, and resolved scalar values.

This module is a non-authoritative migration seam.  It compiles the currently
supported method declarations beside the legacy lmfit path and never mutates
parameter definitions, configuration, or committed Analysis Values.
"""

from __future__ import annotations

import ast
import hashlib
import inspect
import io
import json
import math
import operator
import re
import tokenize
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import ClassVar, Literal, cast
from uuid import uuid4

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.nmr.rates import rate_functions
from chemex.parameters.name import ParamName, matches_parameter_index_selector
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
    """One construction contribution to a parameter's scientific baseline."""

    param_id: str
    supports_estimation: bool
    model_expression: str
    contributor: str
    model_owned: bool = False


@dataclass(frozen=True, slots=True)
class ParameterDeclaration:
    """Sealed scientific inputs used to build each method-local baseline."""

    param_id: str
    supports_estimation: bool
    model_expression: str = ""
    model_owned: bool = False


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
    """Seal additive estimation support and model derivations in definition order."""
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

    def __post_init__(self) -> None:
        roles = tuple(self._roles)
        role_index = MappingProxyType(dict(roles))
        if tuple(role_index) != self.program.scope_ids:
            msg = "Parameter roles do not match the constraint-program scope"
            raise ValueError(msg)
        object.__setattr__(self, "_roles", roles)
        object.__setattr__(self, "_role_index", role_index)
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

    def resolve(self, frame: IndependentValueFrame) -> ResolvedParameterValues:
        _validate_frame(self, frame)
        self.binder.validate_implementations()
        values = {
            param_id: _finite_scalar(value, param_id=param_id)
            for param_id, value in frame._items
        }
        constraints = {item.target_id: item for item in self.program.constraints}
        for position, target_id in enumerate(self.program.evaluation_order):
            constraint = constraints[target_id]
            values[target_id] = _evaluate_expression(
                constraint.expression,
                values,
                self.binder,
                target_id=target_id,
                position=position,
            )
        return ResolvedParameterValues(
            parameterization_identity=self.identity,
            program_fingerprint=self.program.fingerprint,
            occurrence_identity=self.occurrence_identity,
            revision=self.source_revision,
            _items=tuple((param_id, values[param_id]) for param_id in self.scope_ids),
        )


@dataclass(frozen=True, slots=True)
class _RoleRule:
    role: ParameterRole
    selector: str
    expression_text: str
    matches: tuple[str, ...]
    ordinal: int


def _definition_name(definition: ParamDefinition) -> ParamName:
    return ParamName(
        definition.name,
        SpinSystem.from_name(definition.spin_system_name),
        Conditions.model_construct(None, **dict(definition.condition_entries)),
    )


def _match_selector(
    selector_text: str,
    definitions: SealedDefinitions,
    *,
    role: ParameterRole,
    ordinal: int,
) -> tuple[str, ...]:
    selector = ParamName.from_section(selector_text)
    matches = tuple(
        definition.param_id
        for definition in definitions
        if matches_parameter_index_selector(selector, _definition_name(definition))
    )
    if not matches:
        raise NoParameterMatchError(
            "No sealed parameter matches a method selector",
            selector=selector_text,
            role=role.value,
            ordinal=ordinal,
        )
    return matches


def _scan_selectors(text: str) -> tuple[tuple[int, int, str], ...]:
    selectors: list[tuple[int, int, str]] = []
    position = 0
    while position < len(text):
        if text[position] != "[":
            position += 1
            continue
        start = position
        depth = 1
        position += 1
        while position < len(text) and depth:
            if text[position] == "[":
                depth += 1
            elif text[position] == "]":
                depth -= 1
            position += 1
        if depth:
            raise UnsupportedConstraintExpressionError(
                "Unclosed parameter selector",
                expression=text,
            )
        selectors.append((start, position, text[start + 1 : position - 1]))
    return tuple(selectors)


def _split_constraint(text: str) -> tuple[str, str]:
    if text.count("=") != 1:
        raise UnsupportedConstraintExpressionError(
            "A method constraint must contain exactly one '='",
            expression=text,
        )
    left, right = text.split("=", maxsplit=1)
    selectors = _scan_selectors(left)
    if (
        len(selectors) != 1
        or left[: selectors[0][0]].strip()
        or left[selectors[0][1] :].strip()
    ):
        raise UnsupportedConstraintExpressionError(
            "Constraint left side must be one bracketed selector",
            expression=text,
        )
    if not right.strip():
        raise UnsupportedConstraintExpressionError(
            "Constraint right side cannot be empty",
            expression=text,
        )
    return selectors[0][2], right.strip()


def _validate_public_numeric_tokens(
    source: str,
    references: Mapping[str, str],
    expression: str,
) -> None:
    try:
        tokens = tokenize.generate_tokens(io.StringIO(source).readline)
        for token in tokens:
            if token.type in {
                tokenize.NEWLINE,
                tokenize.NL,
                tokenize.ENDMARKER,
                tokenize.ENCODING,
            }:
                continue
            if token.type == tokenize.NUMBER and _PUBLIC_DECIMAL.fullmatch(
                token.string
            ):
                continue
            if token.type == tokenize.NAME and token.string in references:
                continue
            if token.type == tokenize.OP and token.string in _PUBLIC_OPERATORS:
                continue
            raise UnsupportedConstraintExpressionError(
                "Method expression contains syntax outside ChemEx numeric grammar",
                expression=expression,
                token=token.string,
            )
    except (IndentationError, tokenize.TokenError) as error:
        raise UnsupportedConstraintExpressionError(
            "Constraint is not a valid scalar expression",
            expression=expression,
        ) from error


def _validate_public_expression_syntax(text: str) -> None:
    source, references = _replace_selectors(text)
    _validate_public_numeric_tokens(source, references, text)
    try:
        parsed = ast.parse(source, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise UnsupportedConstraintExpressionError(
            "Constraint is not a valid scalar expression",
            expression=text,
        ) from error

    def validate(node: ast.AST) -> None:
        if isinstance(node, ast.Constant):
            _compile_literal(node, "method-syntax")
            return
        if isinstance(node, ast.Name) and node.id in references:
            return
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            validate(node.operand)
            return
        if isinstance(node, ast.BinOp) and isinstance(
            node.op,
            (ast.Add, ast.Sub, ast.Mult, ast.Div),
        ):
            validate(node.left)
            validate(node.right)
            return
        raise UnsupportedConstraintExpressionError(
            "Method expression contains unsupported scalar syntax",
            expression=text,
            syntax=type(node).__name__,
        )

    validate(parsed.body)


def _build_rules(
    method: Method, definitions: SealedDefinitions
) -> tuple[_RoleRule, ...]:
    rules: list[_RoleRule] = []
    ordinal = 0
    for text in method.constraints:
        selector, expression = _split_constraint(text)
        _validate_public_expression_syntax(expression)
        rules.append(
            _RoleRule(
                ParameterRole.DERIVED,
                selector,
                expression,
                _match_selector(
                    selector,
                    definitions,
                    role=ParameterRole.DERIVED,
                    ordinal=ordinal,
                ),
                ordinal,
            )
        )
        ordinal += 1
    for role, selectors in (
        (ParameterRole.FIX, method.fix),
        (ParameterRole.FIT, method.fit),
    ):
        for selector in selectors:
            rules.append(
                _RoleRule(
                    role,
                    selector,
                    "",
                    _match_selector(
                        selector,
                        definitions,
                        role=role,
                        ordinal=ordinal,
                    ),
                    ordinal,
                )
            )
            ordinal += 1
    return tuple(rules)


def _role_for(
    param_id: str,
    declaration: ParameterDeclaration,
    rules: Sequence[_RoleRule],
) -> tuple[ParameterRole, str, str]:
    role = (
        ParameterRole.DERIVED
        if declaration.model_expression
        else ParameterRole.FIT
        if declaration.supports_estimation
        else ParameterRole.FIX
    )
    expression = declaration.model_expression
    source = "model" if declaration.model_owned else "baseline"
    for rule in rules:
        if param_id not in rule.matches:
            continue
        role = rule.role
        expression = rule.expression_text
        source = f"method-rule:{rule.ordinal}"
    return role, expression, source


def _validate_model_derivation_authority(
    rules: Sequence[_RoleRule],
    declarations: SealedParameterDeclarations,
) -> None:
    for rule in rules:
        derived_matches = tuple(
            param_id
            for param_id in rule.matches
            if declarations[param_id].model_expression
            and declarations[param_id].model_owned
        )
        if derived_matches:
            raise ModelDerivationOverrideError(
                "Method rule cannot override a model-owned derivation",
                selector=rule.selector,
                role=rule.role.value,
                ordinal=rule.ordinal,
                param_ids=derived_matches,
            )


def _replace_selectors(right: str) -> tuple[str, Mapping[str, str]]:
    selectors = _scan_selectors(right)
    replacements: dict[str, str] = {}
    parts: list[str] = []
    position = 0
    for index, (start, end, selector) in enumerate(selectors):
        placeholder = f"__selector_{index}"
        parts.extend((right[position:start], placeholder))
        replacements[placeholder] = selector
        position = end
    parts.append(right[position:])
    return "".join(parts), MappingProxyType(replacements)


def _compatible_context(
    candidate: ParamDefinition,
    target: ParamDefinition,
    selector: ParamName,
) -> tuple[frozenset[str], int] | None:
    matched: set[str] = set()
    extras = 0
    candidate_spin = SpinSystem.from_name(candidate.spin_system_name)
    target_spin = SpinSystem.from_name(target.spin_system_name)
    if candidate_spin and not selector.spin_system:
        if target_spin and (
            candidate_spin == target_spin
            or candidate_spin.match(target_spin)
            or target_spin.match(candidate_spin)
        ):
            matched.add("spin_system")
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
    return frozenset(matched), extras


def _resolve_reference(
    selector_text: str,
    target_id: str,
    definitions: SealedDefinitions,
) -> str:
    selector = ParamName.from_section(selector_text)
    candidates = tuple(
        definition
        for definition in definitions
        if matches_parameter_index_selector(selector, _definition_name(definition))
    )
    if not candidates:
        raise NoParameterMatchError(
            "No sealed parameter matches a constraint reference",
            selector=selector_text,
            target_id=target_id,
        )
    non_self = tuple(item for item in candidates if item.param_id != target_id)
    if not non_self:
        raise ConstraintSelfReferenceError(
            "Constraint reference resolves only to its target",
            selector=selector_text,
            target_id=target_id,
        )
    target = definitions[target_id]
    ranked = tuple(
        (candidate, context)
        for candidate in non_self
        if (context := _compatible_context(candidate, target, selector)) is not None
    )
    if not ranked:
        raise NoParameterMatchError(
            "No context-compatible parameter matches a constraint reference",
            selector=selector_text,
            target_id=target_id,
        )
    minimum_extras = min(context[1] for _candidate, context in ranked)
    eligible = tuple(
        (candidate, context[0])
        for candidate, context in ranked
        if context[1] == minimum_extras
    )
    maximal = tuple(
        candidate
        for candidate, fields in eligible
        if not any(fields < other_fields for _other, other_fields in eligible)
    )
    if len(maximal) != 1:
        raise AmbiguousParameterReferenceError(
            "Constraint reference has multiple equally specific matches",
            selector=selector_text,
            target_id=target_id,
            candidate_ids=tuple(item.param_id for item in maximal),
        )
    return maximal[0].param_id


@dataclass(frozen=True, slots=True)
class _ExpressionCompileContext:
    references: Mapping[str, str]
    definitions: SealedDefinitions
    binder: ScientificFunctionBinder
    target_id: str
    model_owned: bool


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
    if node.id in context.references:
        param_id = _resolve_reference(
            context.references[node.id],
            context.target_id,
            context.definitions,
        )
    elif context.model_owned and node.id.startswith("__"):
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
    else:
        raise UnsupportedConstraintExpressionError(
            "Bare names are not supported in scalar constraints",
            target_id=context.target_id,
            name=node.id,
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
        not context.model_owned
        or not isinstance(node.func, ast.Name)
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
    if not context.model_owned or not isinstance(node.value, ast.Call):
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
    model_owned: bool,
) -> ScalarExpression:
    source = text
    references: Mapping[str, str] = MappingProxyType({})
    if not model_owned:
        source, references = _replace_selectors(text)
    try:
        parsed = ast.parse(source, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise UnsupportedConstraintExpressionError(
            "Constraint is not a valid scalar expression",
            target_id=target_id,
            expression=text,
        ) from error
    return _compile_ast(
        parsed.body,
        _ExpressionCompileContext(
            references,
            definitions,
            binder,
            target_id,
            model_owned,
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


def _compile_active_scope(
    parameter_model: SealedParameterModel,
    rules: Sequence[_RoleRule],
    binder: ScientificFunctionBinder,
    required: set[str],
) -> tuple[
    set[str],
    dict[str, tuple[ParameterRole, str, str]],
    dict[str, CompiledConstraint],
]:
    definitions = parameter_model.definitions
    active = set(required)
    role_data: dict[str, tuple[ParameterRole, str, str]] = {}
    compiled: dict[str, CompiledConstraint] = {}
    pending = True
    while pending:
        pending = False
        for definition in definitions:
            param_id = definition.param_id
            if param_id not in active or param_id in role_data:
                continue
            declaration = parameter_model.declarations[param_id]
            role, expression_text, source = _role_for(param_id, declaration, rules)
            role_data[param_id] = (role, expression_text, source)
            if role is not ParameterRole.DERIVED:
                continue
            expression = _parse_expression(
                expression_text,
                definitions=definitions,
                binder=binder,
                target_id=param_id,
                model_owned=source in {"model", "baseline"},
            )
            dependencies = _dependencies(expression)
            compiled[param_id] = CompiledConstraint(
                param_id,
                expression,
                dependencies,
                source,
                expression_text,
            )
            for dependency in dependencies:
                if dependency not in parameter_model.declarations:
                    raise IncompleteParameterDependenciesError(
                        "Constraint dependency is absent from the sealed model",
                        target_id=param_id,
                        dependency_id=dependency,
                    )
                if dependency not in active:
                    active.add(dependency)
                    pending = True
    return active, role_data, compiled


def _validate_rules_in_active_scope(
    rules: Sequence[_RoleRule],
    active: set[str],
) -> None:
    inactive_rules = tuple(rule for rule in rules if not (set(rule.matches) & active))
    if not inactive_rules:
        return
    rule = inactive_rules[0]
    raise NoParameterMatchError(
        "Method selector has no match in the active dependency scope",
        selector=rule.selector,
        role=rule.role.value,
        ordinal=rule.ordinal,
    )


def compile_active_parameterization(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    method: Method,
    required_ids: Sequence[str] | set[str],
) -> ActiveParameterization:
    """Compile one fresh method-scoped role and constraint program."""
    required = _validate_parameterization_inputs(
        parameter_model,
        snapshot,
        required_ids,
    )

    definitions = parameter_model.definitions
    rules = _build_rules(method, definitions)
    _validate_model_derivation_authority(rules, parameter_model.declarations)
    binder = ScientificFunctionBinder.for_model(parameter_model.model_name)
    definition_order = {
        definition.param_id: position for position, definition in enumerate(definitions)
    }
    active, role_data, compiled = _compile_active_scope(
        parameter_model,
        rules,
        binder,
        required,
    )
    _validate_rules_in_active_scope(rules, active)

    scope_ids = tuple(
        definition.param_id
        for definition in definitions
        if definition.param_id in active
    )
    roles = tuple((param_id, role_data[param_id][0]) for param_id in scope_ids)
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
        compiled,
        definition_order,
    )
    constraints = tuple(compiled[param_id] for param_id in derived_ids)
    program = ConstraintProgram(
        parameter_model_identity=parameter_model.identity,
        model_identity=parameter_model.model_identity,
        definitions_identity=snapshot.definitions_identity,
        configuration_identity=snapshot.configuration_identity,
        function_binder_identity=binder.identity,
        scope_ids=scope_ids,
        independent_ids=independent_ids,
        derived_ids=derived_ids,
        constraints=constraints,
        evaluation_order=evaluation_order,
    )
    return ActiveParameterization(
        program=program,
        binder=binder,
        occurrence_identity=snapshot.occurrence_identity,
        source_revision=snapshot.revision,
        _roles=roles,
    )


def build_initial_analysis_values(
    parameter_model: SealedParameterModel,
) -> Mapping[str, float]:
    """Natively fill missing model-derived revision-zero configuration values."""
    configuration = parameter_model.configuration
    missing_ids = tuple(
        config.param_id for config in configuration if config.effective_value is None
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
        ),
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        bootstrap_snapshot,
        Method(),
        set(parameter_model.declarations),
    )
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
        }
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
