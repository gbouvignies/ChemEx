from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

from chemex.parameters.spin_system import SpinSystem

if TYPE_CHECKING:
    from chemex.parameters.parameterization import SealedParameterModel


class FormatOrigin(StrEnum):
    V1 = "v1"
    V2 = "v2"


@dataclass(frozen=True, slots=True)
class SourceRef:
    filename: Path
    step: str
    field: str
    index: int | None = None
    start: int | None = None
    end: int | None = None


class MethodFormatError(ValueError):
    def __init__(self, message: str, source: SourceRef) -> None:
        super().__init__(message)
        self.message = message
        self.source = source

    def __str__(self) -> str:
        location = f"{self.source.filename}: [{self.source.step}] {self.source.field}"
        if self.source.index is not None:
            location += f"[{self.source.index}]"
        if self.source.start is not None and self.source.end is not None:
            location += f" characters {self.source.start}:{self.source.end}"
        return f"{location}: {self.message}"


@dataclass(frozen=True, slots=True)
class ParameterSelector:
    name: str
    spin_system: str | None = None
    temperature: float | None = None
    h_larmor_frq: float | None = None
    p_total: float | None = None
    l_total: float | None = None
    d2o: float | None = None
    source: SourceRef | None = field(default=None, compare=False, repr=False)

    def render(self) -> str:
        parts = [self.name]
        if self.spin_system is not None:
            parts.append(f"NUC->{self.spin_system}")
        if self.temperature is not None:
            parts.append(f"T->{self.temperature:.1f}C")
        if self.h_larmor_frq is not None:
            parts.append(f"B0->{self.h_larmor_frq:.1f}MHz")
        if self.p_total is not None:
            parts.append(f"[P]->{render_number(self.p_total)}M")
        if self.l_total is not None:
            parts.append(f"[L]->{render_number(self.l_total)}M")
        if self.d2o is not None:
            parts.append(f"D2O->{render_number(self.d2o)}")
        return ", ".join(parts)


@dataclass(frozen=True, slots=True)
class LiteralExpression:
    value: float


@dataclass(frozen=True, slots=True)
class SelectorExpression:
    selector: ParameterSelector


@dataclass(frozen=True, slots=True)
class UnaryExpression:
    operator: Literal["+", "-"]
    operand: ConstraintExpression


@dataclass(frozen=True, slots=True)
class BinaryExpression:
    operator: Literal["+", "-", "*", "/"]
    left: ConstraintExpression
    right: ConstraintExpression


type ConstraintExpression = (
    LiteralExpression | SelectorExpression | UnaryExpression | BinaryExpression
)


def render_number(value: float) -> str:
    normalized = 0.0 if value == 0.0 else value
    rendered = repr(normalized)
    if "e" not in rendered.lower() and "." not in rendered:
        rendered += ".0"
    return rendered


def render_expression(expression: ConstraintExpression) -> str:
    precedence = {"+": 1, "-": 1, "*": 2, "/": 2}

    def render(node: ConstraintExpression) -> str:
        if isinstance(node, LiteralExpression):
            return render_number(node.value)
        if isinstance(node, SelectorExpression):
            return f"[{node.selector.render()}]"
        if isinstance(node, UnaryExpression):
            operand = render(node.operand)
            if isinstance(node.operand, BinaryExpression):
                operand = f"({operand})"
            return f"{node.operator}{operand}"
        node_precedence = precedence[node.operator]
        left = render(node.left)
        right_text = render(node.right)
        if (
            isinstance(node.left, BinaryExpression)
            and precedence[node.left.operator] < node_precedence
        ):
            left = f"({left})"
        if isinstance(node.right, BinaryExpression) and (
            precedence[node.right.operator] < node_precedence
            or precedence[node.right.operator] == node_precedence
        ):
            right_text = f"({right_text})"
        return f"{left} {node.operator} {right_text}"

    return render(expression)


@dataclass(frozen=True, slots=True)
class Constraint:
    target: ParameterSelector
    expression: ConstraintExpression
    source: SourceRef = field(compare=False, repr=False)

    def render(self) -> str:
        return f"[{self.target.render()}] = {render_expression(self.expression)}"


@dataclass(frozen=True, slots=True)
class FitAction:
    selectors: tuple[ParameterSelector, ...]
    source: SourceRef = field(compare=False, repr=False)


@dataclass(frozen=True, slots=True)
class FixAction:
    selectors: tuple[ParameterSelector, ...]
    source: SourceRef = field(compare=False, repr=False)


@dataclass(frozen=True, slots=True)
class ConstrainAction:
    constraints: tuple[Constraint, ...]
    source: SourceRef = field(compare=False, repr=False)


type RoleAction = FitAction | FixAction | ConstrainAction


@dataclass(frozen=True, slots=True)
class ProfileSelection:
    include: tuple[str, ...] | Literal["*"] | None = None
    exclude: tuple[str, ...] | Literal["*"] | None = None


class SearchScale(StrEnum):
    LINEAR = "lin"
    LOGARITHMIC = "log"


@dataclass(frozen=True, slots=True)
class GridRange:
    scale: SearchScale
    low: float
    high: float
    count: int

    def render(self) -> str:
        return (
            f"{self.scale.value}({render_number(self.low)}, "
            f"{render_number(self.high)}, {self.count})"
        )


@dataclass(frozen=True, slots=True)
class GridValues:
    values: tuple[float, ...]

    def render(self) -> str:
        return f"values({', '.join(render_number(value) for value in self.values)})"


type GridSpacing = GridRange | GridValues


@dataclass(frozen=True, slots=True)
class GridAxis:
    selector: ParameterSelector
    spacing: GridSpacing
    source: SourceRef = field(compare=False, repr=False)

    def render(self) -> str:
        return f"[{self.selector.render()}] = {self.spacing.render()}"


@dataclass(frozen=True, slots=True)
class GridSearch:
    axes: tuple[GridAxis, ...]


@dataclass(frozen=True, slots=True)
class DeRange:
    scale: SearchScale
    low: float
    high: float

    def render(self) -> str:
        return (
            f"{self.scale.value}({render_number(self.low)}, {render_number(self.high)})"
        )


@dataclass(frozen=True, slots=True)
class DeCoordinate:
    selector: ParameterSelector
    range: DeRange
    source: SourceRef = field(compare=False, repr=False)

    def render(self) -> str:
        return f"[{self.selector.render()}] = {self.range.render()}"


@dataclass(frozen=True, slots=True)
class DeSearch:
    seed: int
    coordinates: tuple[DeCoordinate, ...]


@dataclass(frozen=True, slots=True)
class ResamplingRequest:
    replicates: int
    seed: int | None = None


@dataclass(frozen=True, slots=True)
class McmcRequest:
    steps: int
    burn: int | None = None
    seed: int | None = None
    thin: int = 1
    walkers: int | None = None
    workers: int | None = None


@dataclass(frozen=True, slots=True)
class StatisticsPlan:
    mc: ResamplingRequest | None = None
    bs: ResamplingRequest | None = None
    bsn: ResamplingRequest | None = None
    mcmc: McmcRequest | None = None


def profile_selection(
    include: object = None,
    exclude: object = None,
) -> ProfileSelection:
    def normalize(value: object) -> tuple[str, ...] | Literal["*"] | None:
        if isinstance(value, str) and value.lower() in {"*", "all"}:
            return "*"
        if isinstance(value, list):
            if any(
                isinstance(item, str) and item.lower() in {"*", "all"} for item in value
            ):
                return "*"
            if not all(isinstance(item, (int, str, SpinSystem)) for item in value):
                msg = "Profile selections must contain only names or residue numbers"
                raise TypeError(msg)
            normalized: list[str] = []
            for item in value:
                if isinstance(item, SpinSystem):
                    normalized.append(str(item))
                elif isinstance(item, (int, str)):
                    normalized.append(str(SpinSystem.from_name(item)))
            return tuple(normalized)
        return None

    return ProfileSelection(normalize(include), normalize(exclude))


@dataclass(frozen=True, slots=True)
class StepPlan:
    name: str
    selection: ProfileSelection = ProfileSelection()
    roles_from: str | None = None
    role_actions: tuple[RoleAction, ...] = ()
    search: GridSearch | DeSearch | None = None
    statistics: StatisticsPlan | None = None


@dataclass(frozen=True, slots=True)
class MethodPlan:
    format_origin: FormatOrigin
    steps: tuple[StepPlan, ...]

    def effective_role_actions(self) -> Mapping[str, tuple[RoleAction, ...]]:
        """Resolve each step's immutable baseline or inherited action chain."""
        effective: dict[str, tuple[RoleAction, ...]] = {}
        for step in self.steps:
            inherited = () if step.roles_from is None else effective[step.roles_from]
            effective[step.name] = (*inherited, *step.role_actions)
        return MappingProxyType(effective)

    def render(self) -> str:
        lines = ["FORMAT_VERSION = 2"]
        for step in self.steps:
            lines.extend(("", f"[{_toml_key(step.name)}]"))
            lines.extend(_render_selection(step.selection))
            if step.roles_from is not None:
                lines.append(f"ROLES_FROM = {json.dumps(step.roles_from)}")
            if step.role_actions:
                lines.append("ROLES = [")
                for action in step.role_actions:
                    if isinstance(action, FixAction):
                        key = "FIX"
                        values = tuple(item.render() for item in action.selectors)
                    elif isinstance(action, FitAction):
                        key = "FIT"
                        values = tuple(item.render() for item in action.selectors)
                    else:
                        key = "CONSTRAIN"
                        values = tuple(item.render() for item in action.constraints)
                    rendered = ", ".join(json.dumps(value) for value in values)
                    lines.append(f"  {{ {key} = [{rendered}] }},")
                lines.append("]")
            if isinstance(step.search, GridSearch):
                lines.extend(("", f"[{_toml_key(step.name)}.SEARCH.GRID]", "AXES = ["))
                lines.extend(
                    f"  {json.dumps(axis.render())}," for axis in step.search.axes
                )
                lines.append("]")
            elif isinstance(step.search, DeSearch):
                lines.extend(
                    (
                        "",
                        f"[{_toml_key(step.name)}.SEARCH.DE]",
                        f"SEED = {step.search.seed}",
                        "COORDINATES = [",
                    )
                )
                lines.extend(
                    f"  {json.dumps(coordinate.render())},"
                    for coordinate in step.search.coordinates
                )
                lines.append("]")
            lines.extend(_render_statistics(step.name, step.statistics))
        return "\n".join(lines) + "\n"

    def validate(self, parameter_model: SealedParameterModel) -> None:
        from chemex.configuration.method_validation import validate_method_plan

        validate_method_plan(self, parameter_model)


_BARE_TOML_KEY = re.compile(r"[A-Za-z0-9_-]+\Z")


def _toml_key(value: str) -> str:
    return value if _BARE_TOML_KEY.fullmatch(value) else json.dumps(value)


def _render_profile_set(value: tuple[str, ...] | Literal["*"]) -> str:
    if value == "*":
        return json.dumps("ALL")
    return f"[{', '.join(json.dumps(item) for item in value)}]"


def _render_selection(selection: ProfileSelection) -> list[str]:
    lines: list[str] = []
    if selection.include is not None:
        lines.append(f"INCLUDE = {_render_profile_set(selection.include)}")
    if selection.exclude is not None:
        lines.append(f"EXCLUDE = {_render_profile_set(selection.exclude)}")
    return lines


def _render_statistics(name: str, statistics: StatisticsPlan | None) -> list[str]:
    if statistics is None:
        return []

    key = _toml_key(name)
    lines: list[str] = []
    for request_name, request in (
        ("MC", statistics.mc),
        ("BS", statistics.bs),
        ("BSN", statistics.bsn),
    ):
        if request is None:
            continue
        lines.extend(_render_resampling(key, request_name, request))
    if statistics.mcmc is not None:
        lines.extend(_render_mcmc(key, statistics.mcmc))
    return lines


def _render_resampling(
    step_key: str,
    request_name: str,
    request: ResamplingRequest,
) -> list[str]:
    lines = [
        "",
        f"[{step_key}.STATISTICS.{request_name}]",
        f"REPLICATES = {request.replicates}",
    ]
    if request.seed is not None:
        lines.append(f"SEED = {request.seed}")
    return lines


def _render_mcmc(step_key: str, request: McmcRequest) -> list[str]:
    lines = [
        "",
        f"[{step_key}.STATISTICS.MCMC]",
        f"STEPS = {request.steps}",
    ]
    if request.burn is not None:
        lines.append(f"BURN = {request.burn}")
    if request.seed is not None:
        lines.append(f"SEED = {request.seed}")
    if request.thin != 1:
        lines.append(f"# V1_ONLY_THIN = {request.thin}")
    if request.walkers is not None:
        lines.append(f"# V1_ONLY_WALKERS = {request.walkers}")
    if request.workers is not None:
        lines.append(f"# V1_ONLY_WORKERS = {request.workers}")
    return lines
