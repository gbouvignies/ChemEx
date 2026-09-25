from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, cast

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.configuration.method_plan import (
    BinaryExpression,
    ConstrainAction,
    Constraint,
    ConstraintExpression,
    DeSearch,
    FitAction,
    FixAction,
    GridAxis,
    GridSearch,
    GridValues,
    MethodFormatError,
    MethodPlan,
    ParameterSelector,
    SearchScale,
    SelectorExpression,
    SourceRef,
    UnaryExpression,
    render_expression,
)
from chemex.configuration.method_plan import (
    LiteralExpression as MethodLiteralExpression,
)
from chemex.models.kinetic._binding_migration import (
    binding_rate_migration,
    legacy_binding_migration_message,
)
from chemex.parameters.name import ParamName, matches_parameter_index_selector
from chemex.parameters.parameterization import (
    BinaryExpression as ValueBinaryExpression,
)
from chemex.parameters.parameterization import (
    CompiledConstraint,
    ParameterRole,
    ReferenceExpression,
    ScalarExpression,
    ScientificFunctionBinder,
    SealedParameterModel,
    baseline_parameter_role,
    compatible_reference_context,
    compile_model_constraint,
)
from chemex.parameters.parameterization import (
    LiteralExpression as ValueLiteralExpression,
)
from chemex.parameters.parameterization import (
    UnaryExpression as ValueUnaryExpression,
)
from chemex.parameters.sealed import ParamDefinition
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.temperature_shifts import (
    canonical_control_guidance,
)


@dataclass(frozen=True, slots=True)
class _ResolvedConstraint:
    declaration: Constraint
    dependencies: tuple[str, ...]
    compiled: CompiledConstraint


@dataclass(frozen=True, slots=True)
class ResolvedDeCoordinate:
    """One canonical DE coordinate resolved to a stable independent parameter ID."""

    param_id: str
    low: float
    high: float
    scale: SearchScale


@dataclass(frozen=True, slots=True)
class ResolvedGridAxis:
    """One GRID declaration resolved inside the current active FIT scope."""

    param_id: str
    values: tuple[float, ...]
    declaration_ordinal: int


@dataclass(frozen=True, slots=True)
class MatchedGridAxis:
    matches: tuple[str, ...]
    values: tuple[float, ...]
    source: SourceRef
    ordinal: int
    selector_text: str


@dataclass(frozen=True, slots=True)
class ResolvedMethodStep:
    roles: Mapping[str, ParameterRole]
    constraints: Mapping[str, CompiledConstraint]
    grid_axes: tuple[MatchedGridAxis, ...] = ()
    de_coordinates: tuple[ResolvedDeCoordinate, ...] = ()


def _param_name(definition: ParamDefinition) -> ParamName:
    return ParamName(
        definition.name,
        SpinSystem.from_name(definition.spin_system_name),
        Conditions.model_construct(None, **dict(definition.condition_entries)),
    )


def _selector_name(selector: ParameterSelector) -> ParamName:
    return ParamName(
        selector.name,
        SpinSystem.from_name(selector.spin_system or ""),
        Conditions.model_construct(
            None,
            temperature=selector.temperature,
            h_larmor_frq=selector.h_larmor_frq,
            p_total=selector.p_total,
            l_total=selector.l_total,
            d2o=selector.d2o,
        ),
    )


def _source(selector: ParameterSelector, fallback: SourceRef) -> SourceRef:
    return selector.source if selector.source is not None else fallback


def _matches(
    selector: ParameterSelector,
    model: SealedParameterModel,
    fallback: SourceRef,
) -> tuple[str, ...]:
    parsed = _selector_name(selector)
    matches = tuple(
        definition.param_id
        for definition in model.definitions
        if matches_parameter_index_selector(parsed, _param_name(definition))
    )
    if not matches:
        raise MethodFormatError(
            f"No parameter matches selector [{selector.render()}]",
            _source(selector, fallback),
            detail_code="no_match",
            detail_context={"selector": selector.render().lower()},
        )
    return matches


def _resolve_constraint_reference(
    selector: ParameterSelector,
    target_id: str,
    model: SealedParameterModel,
    source: SourceRef,
) -> str:
    parsed = _selector_name(selector)
    candidates = tuple(
        definition
        for definition in model.definitions
        if matches_parameter_index_selector(parsed, _param_name(definition))
    )
    if not candidates:
        raise MethodFormatError(
            f"No parameter matches constraint reference [{selector.render()}]",
            source,
            detail_code="no_match",
            detail_context={
                "selector": selector.render().lower(),
                "target_id": target_id,
            },
        )
    non_self = tuple(item for item in candidates if item.param_id != target_id)
    if not non_self:
        raise MethodFormatError(
            f"Constraint reference [{selector.render()}] resolves only to its target",
            source,
            detail_code="self_reference",
            detail_context={
                "selector": selector.render().lower(),
                "target_id": target_id,
            },
        )
    target = model.definitions[target_id]
    ranked = tuple(
        (candidate, context)
        for candidate in non_self
        if (context := compatible_reference_context(candidate, target, parsed))
        is not None
    )
    if not ranked:
        raise MethodFormatError(
            f"No context-compatible parameter matches [{selector.render()}]",
            source,
            detail_code="no_match",
            detail_context={
                "selector": selector.render().lower(),
                "target_id": target_id,
            },
        )
    maximum_spin_specificity = max(context[0] for _candidate, context in ranked)
    spin_eligible = tuple(
        (candidate, context)
        for candidate, context in ranked
        if context[0] == maximum_spin_specificity
    )
    minimum_extras = min(context[2] for _candidate, context in spin_eligible)
    eligible = tuple(
        (candidate, context[1])
        for candidate, context in spin_eligible
        if context[2] == minimum_extras
    )
    maximal = tuple(
        candidate
        for candidate, fields in eligible
        if not any(fields < other_fields for _other, other_fields in eligible)
    )
    if len(maximal) != 1:
        candidate_ids = tuple(item.param_id for item in maximal)
        raise MethodFormatError(
            f"Constraint reference [{selector.render()}] is ambiguous among "
            f"{candidate_ids}",
            source,
            detail_code="ambiguity",
            detail_context={
                "selector": selector.render().lower(),
                "target_id": target_id,
                "candidate_ids": candidate_ids,
            },
        )
    return maximal[0].param_id


def _references(expression: ConstraintExpression) -> Iterator[ParameterSelector]:
    if isinstance(expression, SelectorExpression):
        yield expression.selector
    elif isinstance(expression, UnaryExpression):
        yield from _references(expression.operand)
    elif isinstance(expression, BinaryExpression):
        yield from _references(expression.left)
        yield from _references(expression.right)


def _action_selectors(
    action: FitAction | FixAction | ConstrainAction,
) -> Iterator[tuple[ParameterSelector, SourceRef]]:
    if isinstance(action, (FitAction, FixAction)):
        for selector in action.selectors:
            yield selector, _source(selector, action.source)
        return
    for constraint in action.constraints:
        yield constraint.target, constraint.source
        for selector in _references(constraint.expression):
            yield selector, _source(selector, constraint.source)


def _method_selectors(
    plan: MethodPlan,
) -> Iterator[tuple[ParameterSelector, SourceRef]]:
    for step in plan.steps:
        for action in step.role_actions:
            yield from _action_selectors(action)
        search = step.search
        if isinstance(search, GridSearch):
            for axis in search.axes:
                yield axis.selector, _source(axis.selector, axis.source)
        elif isinstance(search, DeSearch):
            for coordinate in search.coordinates:
                yield (
                    coordinate.selector,
                    _source(coordinate.selector, coordinate.source),
                )


def _validate_legacy_binding_selectors(
    plan: MethodPlan,
    model: SealedParameterModel,
) -> None:
    migration = binding_rate_migration(model.model_name)
    if migration is None:
        return
    selectors = tuple(_method_selectors(plan))
    grouped: dict[str, tuple[ParamName, set[str], SourceRef]] = {}
    for selector, source in selectors:
        name = selector.name.upper()
        if name not in migration.legacy_names | migration.replacement_names:
            continue
        parsed = _selector_name(selector)
        scope = ParamName("", parsed.spin_system, parsed.conditions)
        _stored_scope, names, _stored_source = grouped.setdefault(
            scope.id_,
            (scope, set(), source),
        )
        names.add(name)
    affected = tuple(
        (scope, names, source)
        for scope, names, source in grouped.values()
        if migration.legacy_names & names
    )
    if not affected:
        return
    messages = tuple(
        legacy_binding_migration_message(migration, names, scope=scope)
        for scope, names, _source in affected
    )
    source = affected[0][2]
    first_legacy_name = sorted(migration.legacy_names & affected[0][1])[0]
    raise MethodFormatError(
        f"{first_legacy_name} is a derived output. " + "\n".join(messages),
        source,
    )


def _check_bounds(
    param_id: str,
    values: tuple[float, ...],
    model: SealedParameterModel,
    source: SourceRef,
) -> None:
    configuration = model.configuration[param_id]
    if any(
        value < configuration.lower_bound or value > configuration.upper_bound
        for value in values
    ):
        raise MethodFormatError(
            f"Search range for {param_id} lies outside physical bounds "
            f"[{configuration.lower_bound}, {configuration.upper_bound}]",
            source,
        )


def _reject_protected(
    matches: tuple[str, ...],
    model: SealedParameterModel,
    source: SourceRef,
    operation: str,
) -> None:
    protected = tuple(
        param_id for param_id in matches if model.declarations[param_id].model_owned
    )
    if protected:
        guidance = tuple(
            dict.fromkeys(
                advice
                for param_id in protected
                if (
                    advice := canonical_control_guidance(
                        model.definitions[param_id].name
                    )
                )
                is not None
            )
        )
        suffix = f"; {'; '.join(guidance)}" if guidance else ""
        raise MethodFormatError(
            f"{operation} cannot override model-owned parameters {protected}{suffix}",
            source,
            detail_code="model_derivation_override",
            detail_context={"param_ids": protected},
        )


def _reject_protected_constants(
    matches: tuple[str, ...],
    model: SealedParameterModel,
    source: SourceRef,
    operation: str,
) -> None:
    protected = tuple(
        param_id
        for param_id in matches
        if model.declarations[param_id].model_owned
        and not model.declarations[param_id].model_expression
    )
    if protected:
        raise MethodFormatError(
            f"{operation} cannot use protected model constants {protected}",
            source,
            detail_code="model_derivation_override",
            detail_context={"param_ids": protected},
        )


def _reject_unestimable(
    matches: tuple[str, ...],
    model: SealedParameterModel,
    source: SourceRef,
) -> None:
    unestimable = tuple(
        param_id
        for param_id in matches
        if not model.declarations[param_id].supports_estimation
    )
    if unestimable:
        raise MethodFormatError(
            f"FIT cannot override parameters that do not support estimation "
            f"{unestimable}",
            source,
            detail_code="incompatible_input",
            detail_context={"param_ids": unestimable},
        )


def _compile_method_expression(
    expression: ConstraintExpression,
    target_id: str,
    model: SealedParameterModel,
    source: SourceRef,
    dependencies: list[str],
) -> ScalarExpression:
    if isinstance(expression, MethodLiteralExpression):
        return ValueLiteralExpression(expression.value)
    if isinstance(expression, SelectorExpression):
        reference = expression.selector
        reference_source = _source(reference, source)
        param_id = _resolve_constraint_reference(
            reference, target_id, model, reference_source
        )
        _reject_protected_constants(
            _matches(reference, model, reference_source),
            model,
            reference_source,
            "Constraint reference",
        )
        dependencies.append(param_id)
        return ReferenceExpression(param_id)
    if isinstance(expression, UnaryExpression):
        return ValueUnaryExpression(
            "positive" if expression.operator == "+" else "negative",
            _compile_method_expression(
                expression.operand, target_id, model, source, dependencies
            ),
        )
    operators = {
        "+": "add",
        "-": "subtract",
        "*": "multiply",
        "/": "divide",
    }
    operator = cast(
        Literal["add", "subtract", "multiply", "divide"],
        operators[expression.operator],
    )
    return ValueBinaryExpression(
        operator,
        _compile_method_expression(
            expression.left, target_id, model, source, dependencies
        ),
        _compile_method_expression(
            expression.right, target_id, model, source, dependencies
        ),
    )


def _apply_actions(
    roles: dict[str, ParameterRole],
    constraints: dict[str, _ResolvedConstraint],
    actions: tuple[FitAction | FixAction | ConstrainAction, ...],
    model: SealedParameterModel,
    ordinal: int = 0,
) -> int:
    for action in actions:
        if isinstance(action, (FitAction, FixAction)):
            role = (
                ParameterRole.FIT
                if isinstance(action, FitAction)
                else ParameterRole.FIX
            )
            for selector in action.selectors:
                matches = _matches(selector, model, action.source)
                _reject_protected(
                    matches, model, _source(selector, action.source), "Method role"
                )
                if isinstance(action, FitAction):
                    _reject_unestimable(
                        matches, model, _source(selector, action.source)
                    )
                roles.update(dict.fromkeys(matches, role))
                for param_id in matches:
                    constraints.pop(param_id, None)
                ordinal += 1
            continue
        for constraint in action.constraints:
            matches = _matches(constraint.target, model, constraint.source)
            _reject_protected(matches, model, constraint.source, "Constraint")
            for param_id in matches:
                dependencies_found: list[str] = []
                expression = _compile_method_expression(
                    constraint.expression,
                    param_id,
                    model,
                    constraint.source,
                    dependencies_found,
                )
                dependencies = tuple(dict.fromkeys(dependencies_found))
                roles[param_id] = ParameterRole.DERIVED
                constraints[param_id] = _ResolvedConstraint(
                    constraint,
                    dependencies,
                    CompiledConstraint(
                        param_id,
                        expression,
                        dependencies,
                        f"method-rule:{ordinal}",
                        render_expression(constraint.expression),
                    ),
                )
            ordinal += 1
    return ordinal


def _find_cycle(
    constraints: Mapping[str, CompiledConstraint], order: tuple[str, ...]
) -> tuple[str, ...] | None:
    state = dict.fromkeys(constraints, 0)
    stack: list[str] = []
    positions: dict[str, int] = {}

    def visit(param_id: str) -> tuple[str, ...] | None:
        state[param_id] = 1
        positions[param_id] = len(stack)
        stack.append(param_id)
        for dependency in constraints[param_id].dependencies:
            if dependency not in constraints:
                continue
            if state[dependency] == 0 and (cycle := visit(dependency)):
                return cycle
            if state[dependency] == 1:
                return tuple(stack[positions[dependency] :])
        stack.pop()
        positions.pop(param_id)
        state[param_id] = 2
        return None

    for param_id in order:
        if (
            param_id in constraints
            and state[param_id] == 0
            and (cycle := visit(param_id))
        ):
            return cycle
    return None


def _validate_constraint_graph(
    constraints: dict[str, _ResolvedConstraint],
    roles: Mapping[str, ParameterRole],
    model: SealedParameterModel,
    binder: ScientificFunctionBinder,
    model_constraints: dict[str, CompiledConstraint],
    step_name: str,
) -> None:
    order = tuple(definition.param_id for definition in model.definitions)
    graph: dict[str, CompiledConstraint] = {}
    for param_id in order:
        if roles[param_id] is not ParameterRole.DERIVED:
            continue
        if param_id in constraints:
            graph[param_id] = constraints[param_id].compiled
            continue
        if param_id not in model_constraints:
            model_constraints[param_id] = compile_model_constraint(
                model, binder, param_id
            )
        graph[param_id] = model_constraints[param_id]
    cycle = _find_cycle(graph, order)
    if cycle is not None:
        method_member = next(
            (param_id for param_id in cycle if param_id in constraints), None
        )
        source = (
            constraints[method_member].declaration.source
            if method_member is not None
            else SourceRef(Path("<sealed-parameter-model>"), step_name, "DERIVATIONS")
        )
        raise MethodFormatError(
            f"Constraint dependency cycle contains {', '.join(cycle)}",
            source,
            detail_code="cycle",
            detail_context={
                "param_ids": cycle,
                "constraints": tuple(
                    (
                        param_id,
                        graph[param_id].source,
                        graph[param_id].expression_text,
                    )
                    for param_id in cycle
                ),
            },
        )


def _validate_grid(
    search: GridSearch,
    roles: dict[str, ParameterRole],
    model: SealedParameterModel,
) -> tuple[MatchedGridAxis, ...]:
    resolved: list[MatchedGridAxis] = []
    for ordinal, axis in enumerate(search.axes):
        matches = _matches(axis.selector, model, axis.source)
        _reject_protected(matches, model, axis.source, "GRID")
        if not any(roles[param_id] is ParameterRole.FIT for param_id in matches):
            raise MethodFormatError(
                "GRID target is not a final independent FIT coordinate",
                axis.source,
            )
        resolved.append(
            MatchedGridAxis(
                matches,
                _grid_values(axis),
                axis.source,
                ordinal,
                axis.selector.render(),
            )
        )
    return tuple(resolved)


def _grid_values(axis: GridAxis) -> tuple[float, ...]:
    spacing = axis.spacing
    if isinstance(spacing, GridValues):
        return spacing.values
    if spacing.scale is SearchScale.LINEAR:
        values = np.linspace(spacing.low, spacing.high, spacing.count)
    else:
        values = np.geomspace(spacing.low, spacing.high, spacing.count)
    return tuple(float(value) for value in values)


def project_grid_axes(
    axes: tuple[MatchedGridAxis, ...],
    model: SealedParameterModel,
    *,
    active_scope_ids: tuple[str, ...],
    final_fit_ids: tuple[str, ...],
) -> tuple[ResolvedGridAxis, ...]:
    """Resolve broad GRID rules against one current active final FIT scope.

    Declarations retain v2's top-to-bottom rule semantics: a later declaration
    replaces an earlier declaration for every concrete active coordinate it
    matches. Parameters outside the active scope are deliberately ignored.
    """
    active_scope = frozenset(active_scope_ids)
    final_fit = frozenset(final_fit_ids)
    concrete: dict[str, ResolvedGridAxis] = {}
    sources: dict[str, SourceRef] = {}
    for axis in axes:
        active_matches = tuple(
            param_id for param_id in axis.matches if param_id in active_scope
        )
        if not active_matches:
            raise MethodFormatError(
                "GRID selector has no applicable coordinate in the current "
                f"active step: [{axis.selector_text}]",
                axis.source,
            )
        matches = tuple(
            param_id for param_id in active_matches if param_id in final_fit
        )
        if not matches:
            raise MethodFormatError(
                "GRID selector has no active final independent FIT coordinate; "
                "active non-FIT matches: " + ", ".join(active_matches),
                axis.source,
            )
        for param_id in matches:
            concrete[param_id] = ResolvedGridAxis(param_id, axis.values, axis.ordinal)
            sources[param_id] = axis.source
    for resolved in concrete.values():
        _check_bounds(
            resolved.param_id,
            resolved.values,
            model,
            sources[resolved.param_id],
        )
    active_order = {param_id: index for index, param_id in enumerate(final_fit_ids)}
    return tuple(
        sorted(
            concrete.values(),
            key=lambda item: (item.declaration_ordinal, active_order[item.param_id]),
        )
    )


def resolve_grid_axes(
    search: GridSearch,
    model: SealedParameterModel,
    *,
    active_scope_ids: tuple[str, ...],
    final_fit_ids: tuple[str, ...],
) -> tuple[ResolvedGridAxis, ...]:
    """Standalone compatibility preview, never an executable-plan fit input.

    This assumes global FIT eligibility; only ``resolve_method_plan`` knows the
    effective Method roles. Fitting must use ``compile_method_plan`` instead.
    """
    roles = dict.fromkeys(model.declarations, ParameterRole.FIT)
    axes = _validate_grid(search, roles, model)
    return project_grid_axes(
        axes,
        model,
        active_scope_ids=active_scope_ids,
        final_fit_ids=final_fit_ids,
    )


def _validate_de(
    search: DeSearch,
    roles: dict[str, ParameterRole],
    model: SealedParameterModel,
) -> tuple[ResolvedDeCoordinate, ...]:
    seen: set[str] = set()
    resolved = resolve_de_coordinates(search, model)
    for coordinate, resolved_coordinate in zip(
        search.coordinates,
        resolved,
        strict=True,
    ):
        param_id = resolved_coordinate.param_id
        if roles[param_id] is not ParameterRole.FIT:
            raise MethodFormatError(
                f"DE target {param_id} is not a final independent FIT coordinate",
                coordinate.source,
            )
        if param_id in seen:
            raise MethodFormatError(
                f"Duplicate DE coordinate {param_id}", coordinate.source
            )
        seen.add(param_id)
        _check_bounds(
            param_id,
            (resolved_coordinate.low, resolved_coordinate.high),
            model,
            coordinate.source,
        )
    return resolved


def resolve_de_coordinates(
    search: DeSearch,
    model: SealedParameterModel,
) -> tuple[ResolvedDeCoordinate, ...]:
    """Resolve validated canonical DE coordinates to stable parameter IDs."""
    resolved: list[ResolvedDeCoordinate] = []
    for coordinate in search.coordinates:
        matches = _matches(coordinate.selector, model, coordinate.source)
        _reject_protected(matches, model, coordinate.source, "DE")
        if len(matches) != 1:
            raise MethodFormatError(
                "Each DE entry must resolve to exactly one final independent "
                f"FIT coordinate; matched {len(matches)}",
                coordinate.source,
            )
        resolved.append(
            ResolvedDeCoordinate(
                matches[0],
                coordinate.range.low,
                coordinate.range.high,
                coordinate.range.scale,
            )
        )
    return tuple(resolved)


def resolve_method_plan(
    plan: MethodPlan, model: SealedParameterModel
) -> tuple[ResolvedMethodStep, ...]:
    _validate_legacy_binding_selectors(plan, model)
    baseline = {
        param_id: baseline_parameter_role(declaration)
        for param_id, declaration in model.declarations.items()
    }
    effective_by_step: dict[str, dict[str, ParameterRole]] = {}
    constraints_by_step: dict[str, dict[str, _ResolvedConstraint]] = {}
    ordinals_by_step: dict[str, int] = {}
    resolved_steps: list[ResolvedMethodStep] = []
    binder = ScientificFunctionBinder.for_model(model.model_name)
    model_constraints: dict[str, CompiledConstraint] = {}
    for step in plan.steps:
        if step.name in effective_by_step:
            raise MethodFormatError(
                "Method step names must be unique",
                SourceRef(Path("<method-plan>"), step.name, "NAME"),
            )
        if step.roles_from is not None and step.roles_from not in effective_by_step:
            raise MethodFormatError(
                "ROLES_FROM must name one unique earlier step",
                SourceRef(Path("<method-plan>"), step.name, "ROLES_FROM"),
            )
        roles = dict(
            baseline if step.roles_from is None else effective_by_step[step.roles_from]
        )
        constraints = dict(
            {} if step.roles_from is None else constraints_by_step[step.roles_from]
        )
        ordinal = _apply_actions(
            roles,
            constraints,
            step.role_actions,
            model,
            0 if step.roles_from is None else ordinals_by_step[step.roles_from],
        )
        _validate_constraint_graph(
            constraints, roles, model, binder, model_constraints, step.name
        )
        effective_by_step[step.name] = roles
        constraints_by_step[step.name] = constraints
        ordinals_by_step[step.name] = ordinal
        grid_axes: tuple[MatchedGridAxis, ...] = ()
        de_coordinates: tuple[ResolvedDeCoordinate, ...] = ()
        if isinstance(step.search, GridSearch):
            grid_axes = _validate_grid(step.search, roles, model)
        elif isinstance(step.search, DeSearch):
            de_coordinates = _validate_de(step.search, roles, model)
        resolved_steps.append(
            ResolvedMethodStep(
                MappingProxyType(roles),
                MappingProxyType(
                    {param_id: item.compiled for param_id, item in constraints.items()}
                ),
                grid_axes,
                de_coordinates,
            )
        )
    return tuple(resolved_steps)


def validate_method_plan(plan: MethodPlan, model: SealedParameterModel) -> None:
    resolve_method_plan(plan, model)
