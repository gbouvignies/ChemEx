from __future__ import annotations

from pathlib import Path
from typing import Any

from chemex.configuration.method_expressions import (
    parse_constraint,
    parse_legacy_grid_axis,
    parse_legacy_selector,
)
from chemex.configuration.method_plan import (
    ConstrainAction,
    Constraint,
    FitAction,
    FixAction,
    FormatOrigin,
    GridSearch,
    McmcRequest,
    MethodPlan,
    ParameterSelector,
    ProfileSelection,
    ResamplingRequest,
    SourceRef,
    StatisticsPlan,
    StepPlan,
    profile_selection,
)


def _selector(text: str, source: SourceRef) -> ParameterSelector:
    return parse_legacy_selector(text, source)


def _constraint(text: str, source: SourceRef) -> Constraint:
    return parse_constraint(text, source, _selector)


def _statistics(value: object) -> StatisticsPlan | None:
    if not isinstance(value, dict):
        return None
    settings = {str(key).lower(): item for key, item in value.items()}

    def resampling(name: str) -> ResamplingRequest | None:
        request = settings.get(name)
        return ResamplingRequest(request, seed=0) if isinstance(request, int) else None

    mcmc_value = settings.get("mcmc")
    mcmc = None
    if isinstance(mcmc_value, int):
        mcmc = McmcRequest(mcmc_value)
    elif isinstance(mcmc_value, dict):
        expanded = {str(key).lower(): item for key, item in mcmc_value.items()}

        def integer(name: str, default: int | None = None) -> int | None:
            item = expanded.get(name, default)
            return (
                item
                if isinstance(item, int) and not isinstance(item, bool)
                else default
            )

        burn_value = expanded.get("burn")
        burn = burn_value if isinstance(burn_value, int) else None
        steps = integer("steps")
        if steps is None:
            return None
        mcmc = McmcRequest(
            steps,
            burn,
            integer("seed"),
            integer("thin", 1) or 1,
            integer("walkers"),
            integer("workers"),
        )
    return StatisticsPlan(resampling("mc"), resampling("bs"), resampling("bsn"), mcmc)


def adapt_v1(raw_steps: list[tuple[Path, str, dict[str, Any]]]) -> MethodPlan:
    steps: list[StepPlan] = []
    previous: str | None = None
    previous_selection = ProfileSelection()
    for filename, name, settings in raw_steps:
        normalized = {str(key).lower(): value for key, value in settings.items()}
        selection = (
            profile_selection(
                normalized.get("include"),
                normalized.get("exclude"),
            )
            if "include" in normalized or "exclude" in normalized
            else previous_selection
        )
        actions = []
        constraints = tuple(
            _constraint(
                text,
                SourceRef(filename, name, "CONSTRAINTS", index),
            )
            for index, text in enumerate(normalized.get("constraints", []))
        )
        if constraints:
            actions.append(
                ConstrainAction(
                    constraints,
                    SourceRef(filename, name, "CONSTRAINTS"),
                )
            )
        for key, action_type in (("fix", FixAction), ("fit", FitAction)):
            texts = normalized.get(key, [])
            if texts:
                actions.append(
                    action_type(
                        tuple(
                            _selector(
                                text,
                                SourceRef(filename, name, key.upper(), index),
                            )
                            for index, text in enumerate(texts)
                        ),
                        SourceRef(filename, name, key.upper()),
                    )
                )
        steps.append(
            StepPlan(
                name=name,
                selection=selection,
                roles_from=previous,
                role_actions=tuple(actions),
                search=(
                    GridSearch(
                        tuple(
                            parse_legacy_grid_axis(
                                text,
                                SourceRef(filename, name, "GRID", index),
                            )
                            for index, text in enumerate(normalized.get("grid", []))
                        )
                    )
                    if normalized.get("grid")
                    else None
                ),
                statistics=_statistics(normalized.get("statistics")),
            )
        )
        previous = name
        previous_selection = selection
    return MethodPlan(FormatOrigin.V1, tuple(steps))
