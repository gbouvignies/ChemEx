"""Compatibility input preparation for canonical Method Plan execution."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from chemex.configuration.method_plan import MethodPlan
from chemex.configuration.method_v1 import adapt_v1
from chemex.configuration.methods import Methods

if TYPE_CHECKING:
    from chemex.parameters.parameterization import SealedParameterModel


def _normalize_method_plan(methods: Methods | MethodPlan) -> MethodPlan:
    if isinstance(methods, MethodPlan):
        return methods

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
    return adapt_v1(raw_steps)


def prepare_method_plan(
    methods: Methods | MethodPlan,
    parameter_model: SealedParameterModel,
) -> MethodPlan:
    """Normalize a supported input and validate its authoritative Method Plan."""
    plan = _normalize_method_plan(methods)
    plan.validate(parameter_model)
    return plan
