"""Operational adapters for canonical and legacy method configuration."""

from pathlib import Path

from chemex.configuration.method_plan import (
    GridRange,
    GridSearch,
    GridValues,
    McmcRequest,
    MethodPlan,
    StepPlan,
)
from chemex.configuration.method_v1 import adapt_v1
from chemex.configuration.methods import McmcSettings, Method, Methods, Statistics


def step_names(methods: Methods | MethodPlan) -> tuple[str, ...]:
    """Return method step names in execution order."""
    if isinstance(methods, MethodPlan):
        return tuple(step.name for step in methods.steps)
    return tuple(methods)


def mcmc_settings(request: McmcRequest | None) -> McmcSettings | None:
    """Adapt canonical MCMC settings to the operational v1-shaped model."""
    if request is None:
        return None
    return McmcSettings(
        steps=request.steps,
        burn="auto" if request.burn is None else request.burn,
        thin=request.thin,
        walkers=request.walkers,
        seed=request.seed,
        workers=request.workers,
    )


def _selection_value(value: tuple[str, ...] | str | None) -> list[str] | str | None:
    return list(value) if isinstance(value, tuple) else value


def _grid_entries(step: StepPlan, search: GridSearch) -> list[str]:
    entries: list[str] = []
    for axis in search.axes:
        if isinstance(axis.spacing, GridRange):
            spacing = axis.spacing.render()
        elif isinstance(axis.spacing, GridValues):
            spacing = f"({', '.join(str(value) for value in axis.spacing.values)})"
        else:  # pragma: no cover - closed canonical union
            raise TypeError(f"Unsupported grid spacing in step {step.name!r}")
        entries.append(f"[{axis.selector.render()}] = {spacing}")
    return entries


def operational_method_from_step(step: StepPlan) -> Method:
    """Adapt one canonical step without adding cross-step inheritance."""
    statistics = step.statistics
    settings: dict[str, object] = {
        # Canonical selection is step-local. An omitted INCLUDE therefore
        # restores all profiles that an earlier step may have filtered.
        "include": (
            "*"
            if step.selection.include is None
            else _selection_value(step.selection.include)
        ),
        "grid": (
            _grid_entries(step, step.search)
            if isinstance(step.search, GridSearch)
            else []
        ),
        "statistics": (
            None
            if statistics is None
            else Statistics(
                mc=None if statistics.mc is None else statistics.mc.replicates,
                bs=None if statistics.bs is None else statistics.bs.replicates,
                bsn=None if statistics.bsn is None else statistics.bsn.replicates,
                mcmc=mcmc_settings(statistics.mcmc),
            )
        ),
    }
    if step.selection.exclude is not None:
        settings["exclude"] = _selection_value(step.selection.exclude)
    return Method.model_validate(settings)


def normalize_methods_for_execution(
    methods: Methods | MethodPlan,
) -> tuple[MethodPlan, Methods]:
    """Return canonical semantics plus the existing operational method shape."""
    if isinstance(methods, MethodPlan):
        plan = methods
    else:
        raw_steps = []
        for name, method in methods.items():
            settings: dict[str, object] = {
                "constraints": method.constraints,
                "fix": method.fix,
                "fit": method.fit,
                "grid": method.grid,
            }
            if method.statistics is not None:
                settings["statistics"] = method.statistics.model_dump(exclude_none=True)
            if "include" in method.model_fields_set:
                settings["include"] = method.include
            if "exclude" in method.model_fields_set:
                settings["exclude"] = method.exclude
            raw_steps.append((Path("<runtime-v1-method>"), name, settings))
        plan = adapt_v1(raw_steps)
    operational = {step.name: operational_method_from_step(step) for step in plan.steps}
    return plan, operational
