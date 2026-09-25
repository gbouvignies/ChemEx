"""Test-only bridge for direct numerical tests using Method role actions."""

from __future__ import annotations

from collections.abc import Sequence

from chemex.configuration.method_plan import RoleAction
from chemex.parameters.parameterization import (
    ActiveParameterization,
    compile_active_parameterization_from_actions,
)
from chemex.runtime import AnalysisSession


def preview_actions(
    session: AnalysisSession,
    actions: Sequence[RoleAction],
    required_ids: set[str],
) -> ActiveParameterization:
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    return compile_active_parameterization_from_actions(
        model, session.analysis_values.snapshot(), actions, required_ids
    )
