"""Test-only bridge for direct numerical tests using compiled Method meaning."""

from __future__ import annotations

from chemex.configuration.method_plan import MethodPlan
from chemex.configuration.method_validation import resolve_method_plan
from chemex.parameters.parameterization import (
    ActiveParameterization,
    compile_static_parameterization,
)
from chemex.runtime import AnalysisSession


def preview_plan_step(
    session: AnalysisSession,
    plan: MethodPlan,
    step_name: str,
    required_ids: set[str],
) -> ActiveParameterization:
    """Bind one test scope through the production global Method resolver."""
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    resolved = resolve_method_plan(plan, model)
    index = next(
        index for index, step in enumerate(plan.steps) if step.name == step_name
    )
    meaning = resolved[index]
    static = compile_static_parameterization(
        model, meaning.roles, meaning.constraints, required_ids
    )
    return static.bind(session.analysis_values.snapshot())
