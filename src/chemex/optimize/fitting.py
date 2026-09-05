"""Run-level fitting helpers and the supported Python compatibility façade."""

import shutil
from pathlib import Path

from chemex.configuration.method_input import prepare_method_plan
from chemex.configuration.method_plan import MethodPlan
from chemex.configuration.methods import Methods
from chemex.containers.experiments import Experiments
from chemex.exceptions import ArtifactPublicationError
from chemex.optimize.method_plan_execution import execute_method_plan
from chemex.run_info import RunInfo
from chemex.runtime import AnalysisSession

_CHEMEX_RESULT_PATHS = (
    "Parameters",
    "Data",
    "Plots",
    "Grid",
    "Groups",
    "All",
    "Components",
    "Statistics",
    "statistics.toml",
)


def invalidate_planned_outputs(plan: MethodPlan, path: Path) -> None:
    """Remove only ChemEx-owned results for every method step planned now."""
    names = tuple(step.name for step in plan.steps)
    if len(names) > 1:
        output_root = path.resolve()
        step_roots = tuple(path / section for section in names)
        if any(
            not step_root.resolve().is_relative_to(output_root)
            for step_root in step_roots
        ):
            msg = "A planned method step root is outside the output directory"
            raise ValueError(msg)
    else:
        step_roots = (path,)
    for step_root in step_roots:
        for name in _CHEMEX_RESULT_PATHS:
            result_path = step_root / name
            try:
                if result_path.is_symlink() or result_path.is_file():
                    result_path.unlink()
                elif result_path.is_dir():
                    shutil.rmtree(result_path)
            except OSError as error:
                raise ArtifactPublicationError(
                    "invalidate planned analysis output",
                    result_path,
                    error,
                ) from error


def run_methods(
    experiments: Experiments,
    methods: Methods | MethodPlan,
    path: Path,
    plot_level: str,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
) -> None:
    parameter_model = session.parameter_factory.sealed_parameter_model
    if parameter_model is None:
        raise RuntimeError("Native parameter model is unavailable")
    plan = prepare_method_plan(methods, parameter_model)
    execute_method_plan(
        experiments,
        plan,
        path,
        plot_level,
        session=session,
        run_info=run_info,
    )
