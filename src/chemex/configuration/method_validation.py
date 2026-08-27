from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

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
)
from chemex.parameters.name import ParamName, matches_parameter_index_selector
from chemex.parameters.parameterization import (
    ParameterRole,
    SealedParameterModel,
    baseline_parameter_role,
    compatible_reference_context,
)
from chemex.parameters.sealed import ParamDefinition
from chemex.parameters.spin_system import SpinSystem


@dataclass(frozen=True, slots=True)
class _ResolvedConstraint:
    declaration: Constraint
    dependencies: tuple[str, ...]


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
            f"No parameter matches constraint reference [{selector.render()}]", source
        )
    non_self = tuple(item for item in candidates if item.param_id != target_id)
    if not non_self:
        raise MethodFormatError(
            f"Constraint reference [{selector.render()}] resolves only to its target",
            source,
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
            f"No context-compatible parameter matches [{selector.render()}]", source
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
        raise MethodFormatError(
            f"{operation} cannot override model-owned parameters {protected}", source
        )


def _reject_unestimable(
    matches: tuple[str, ...],
    model: SealedParameterModel,
    source: SourceRef,
) -> None:
    unestimable = tuple(
        param_id
        for param_id in matches
        if model.declarations[param_id].model_expression
        and not model.declarations[param_id].supports_estimation
    )
    if unestimable:
        raise MethodFormatError(
            f"FIT cannot override parameters that do not support estimation "
            f"{unestimable}",
            source,
        )


def _apply_actions(
    roles: dict[str, ParameterRole],
    constraints: dict[str, _ResolvedConstraint],
    actions: tuple[FitAction | FixAction | ConstrainAction, ...],
    model: SealedParameterModel,
) -> None:
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
            continue
        for constraint in action.constraints:
            matches = _matches(constraint.target, model, constraint.source)
            _reject_protected(matches, model, constraint.source, "Constraint")
            for param_id in matches:
                dependencies = tuple(
                    dict.fromkeys(
                        _resolve_constraint_reference(
                            reference,
                            param_id,
                            model,
                            _source(reference, constraint.source),
                        )
                        for reference in _references(constraint.expression)
                    )
                )
                roles[param_id] = ParameterRole.DERIVED
                constraints[param_id] = _ResolvedConstraint(constraint, dependencies)


def _find_cycle(
    constraints: dict[str, _ResolvedConstraint], order: tuple[str, ...]
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
    constraints: dict[str, _ResolvedConstraint], model: SealedParameterModel
) -> None:
    order = tuple(definition.param_id for definition in model.definitions)
    cycle = _find_cycle(constraints, order)
    if cycle is not None:
        source = constraints[cycle[0]].declaration.source
        raise MethodFormatError(
            f"Constraint dependency cycle contains {', '.join(cycle)}", source
        )


def _validate_grid(
    search: GridSearch,
    roles: dict[str, ParameterRole],
    model: SealedParameterModel,
) -> None:
    for axis in search.axes:
        matches = _matches(axis.selector, model, axis.source)
        if not any(roles[param_id] is ParameterRole.FIT for param_id in matches):
            raise MethodFormatError(
                "GRID target is not a final independent FIT coordinate",
                axis.source,
            )


def _grid_values(axis: GridAxis) -> tuple[float, ...]:
    spacing = axis.spacing
    if isinstance(spacing, GridValues):
        return spacing.values
    if spacing.scale is SearchScale.LINEAR:
        values = np.linspace(spacing.low, spacing.high, spacing.count)
    else:
        values = np.geomspace(spacing.low, spacing.high, spacing.count)
    return tuple(float(value) for value in values)


def resolve_grid_axes(
    search: GridSearch,
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
    for ordinal, axis in enumerate(search.axes):
        active_matches = tuple(
            param_id
            for param_id in _matches(axis.selector, model, axis.source)
            if param_id in active_scope
        )
        if not active_matches:
            raise MethodFormatError(
                "GRID selector has no applicable coordinate in the current "
                f"active step: [{axis.selector.render()}]",
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
        values = _grid_values(axis)
        for param_id in matches:
            concrete[param_id] = ResolvedGridAxis(param_id, values, ordinal)
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


def _validate_de(
    search: DeSearch,
    roles: dict[str, ParameterRole],
    model: SealedParameterModel,
) -> None:
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


def resolve_de_coordinates(
    search: DeSearch,
    model: SealedParameterModel,
) -> tuple[ResolvedDeCoordinate, ...]:
    """Resolve validated canonical DE coordinates to stable parameter IDs."""
    resolved: list[ResolvedDeCoordinate] = []
    for coordinate in search.coordinates:
        matches = _matches(coordinate.selector, model, coordinate.source)
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


def validate_method_plan(plan: MethodPlan, model: SealedParameterModel) -> None:
    baseline = {
        param_id: baseline_parameter_role(declaration)
        for param_id, declaration in model.declarations.items()
    }
    effective_by_step: dict[str, dict[str, ParameterRole]] = {}
    constraints_by_step: dict[str, dict[str, _ResolvedConstraint]] = {}
    for step in plan.steps:
        roles = dict(
            baseline if step.roles_from is None else effective_by_step[step.roles_from]
        )
        constraints = dict(
            {} if step.roles_from is None else constraints_by_step[step.roles_from]
        )
        _apply_actions(roles, constraints, step.role_actions, model)
        _validate_constraint_graph(constraints, model)
        effective_by_step[step.name] = roles
        constraints_by_step[step.name] = constraints
        if isinstance(step.search, GridSearch):
            _validate_grid(step.search, roles, model)
        elif isinstance(step.search, DeSearch):
            _validate_de(step.search, roles, model)
