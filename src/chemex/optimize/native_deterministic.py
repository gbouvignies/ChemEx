"""Production composition for native deterministic method steps."""

from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import cast

import numpy as np

from chemex.configuration.methods import Method
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine, EvaluationResult
from chemex.messages import print_group_name, print_minimizing
from chemex.optimize.direct_trf import AcceptedFitResult, OptimizationProblem
from chemex.optimize.grid_direct_trf import GridDirectTrfInvocation
from chemex.optimize.gridding import (
    GridResult,
    combine_grids,
    make_grids_nd,
    plot_grid_1d,
    plot_grid_2d,
)
from chemex.optimize.grouped_direct_trf import (
    FitDecomposition,
    GroupedDirectTrfInvocation,
)
from chemex.optimize.grouped_grid_direct_trf import GroupedGridDirectTrfOutcome
from chemex.optimize.grouping import Group, create_groups
from chemex.optimize.helper import (
    execute_post_fit,
    execute_post_fit_groups,
    print_header,
    print_values,
)
from chemex.optimize.method_step import (
    EvaluationPurpose,
    MethodStepLifecycle,
    MethodStepStrategy,
    MethodStepWorkflow,
    execute_method_step,
)
from chemex.parameters.parameterization import ActiveParameterization, ParameterRole
from chemex.runtime import AnalysisSession
from chemex.typing import Array

# Match scipy.optimize.least_squares' default maximum when max_nfev is omitted.
_TRF_EVALUATIONS_PER_COORDINATE = 100


@dataclass(frozen=True, slots=True)
class NativeDeterministicFit:
    """Accepted production fit context available to native successor analyses."""

    accepted: AcceptedFitResult
    problem: OptimizationProblem
    parameterization: ActiveParameterization
    engine: EvaluationEngine

    @property
    def objective_request_budget(self) -> int:
        """Return the accepted occurrence's native Direct TRF request budget."""
        return _objective_request_budget(self.problem)


def uses_native_deterministic(method: Method) -> bool:
    """Return whether one whole method occurrence uses native TRF."""
    if method.statistics is not None and method.statistics.mcmc is not None:
        return False
    if method.statistics is not None:
        return not method.grid
    if "fitmethod" not in method.model_fields_set:
        return True
    return method.fitmethod in {"trf", "least_squares"}


def _objective_request_budget(problem: OptimizationProblem) -> int:
    return _TRF_EVALUATIONS_PER_COORDINATE * max(1, len(problem.controlled_ids))


def _has_controlled_parameters(parameterization: ActiveParameterization) -> bool:
    return any(
        parameterization.role(param_id) is ParameterRole.FIT
        for param_id in parameterization.independent_ids
    )


def _materialize_evaluation(
    experiments: Experiments,
    result: EvaluationResult,
) -> None:
    """Copy one authoritative native evaluation into the live output profiles."""
    profiles = tuple(profile for experiment in experiments for profile in experiment)
    if len(profiles) != len(result.profiles):
        raise RuntimeError("Native evaluation profile population changed before output")
    resolved = dict(result.resolved_values.ordered_items())
    for profile, descriptor in zip(profiles, result.profiles, strict=True):
        start = descriptor.observation_offset
        stop = start + descriptor.observation_count
        if descriptor.observation_count != profile.data.size:
            raise RuntimeError("Native evaluation profile shape changed before output")
        profile.data.calc_unscaled = np.array(
            result.unscaled_calculations[start:stop],
            copy=True,
        )
        profile.data.calc = np.array(
            result.normalized_calculations[start:stop],
            copy=True,
        )
        profile.spectrometer.update(
            {
                local_name: resolved[param_id]
                for local_name, param_id in profile.name_map.items()
            }
        )


def _project_residuals(
    subset: Experiments,
    root: Experiments,
    result: EvaluationResult,
) -> Array:
    """Project the root native residual vector onto one established output group."""
    root_profiles = tuple(profile for experiment in root for profile in experiment)
    if len(root_profiles) != len(result.profiles):
        raise RuntimeError("Native evaluation profile population changed before output")
    segments = {
        id(profile): result.residuals[
            descriptor.residual_offset : (
                descriptor.residual_offset + descriptor.residual_count
            )
        ]
        for profile, descriptor in zip(root_profiles, result.profiles, strict=True)
    }
    projected = [
        segments[id(profile)] for experiment in subset for profile in experiment
    ]
    return np.concatenate(projected)


def _build_invocation(
    method: Method,
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    experiments: Experiments,
) -> tuple[
    MethodStepStrategy,
    GroupedDirectTrfInvocation | GridDirectTrfInvocation,
]:
    if method.grid:
        grid = experiments.parameter_store.parse_grid(method.grid)
        invocation = GridDirectTrfInvocation.for_problem(
            problem,
            axes=tuple(
                (param_id, tuple(values.tolist())) for param_id, values in grid.items()
            ),
            objective_request_budget=_objective_request_budget(problem),
        )
        return MethodStepStrategy.GRID_DIRECT_TRF, invocation
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=tuple(
            _objective_request_budget(component.problem)
            for component in decomposition.components
        ),
    )
    return MethodStepStrategy.DIRECT_TRF, invocation


def _write_grid_output(
    outcome: GroupedGridDirectTrfOutcome,
    invocation: GridDirectTrfInvocation,
    path: Path,
    *,
    experiments: Experiments,
    group_scopes: tuple[tuple[Path, frozenset[str]], ...],
) -> None:
    """Render established per-group GRID tables from native seed evidence."""
    grid = {
        axis.param_id: np.asarray(axis.values, dtype=np.float64)
        for axis in invocation.axes
    }
    parameter_store = experiments.parameter_store
    grid_path = path / "Grid"
    grid_path.mkdir(parents=True, exist_ok=True)
    component_results: list[GridResult] = []
    for group_path, controlled_ids in group_scopes:
        matching_indices = tuple(
            index
            for index, component in enumerate(outcome.attempts[0].components)
            if frozenset(component.controlled_ids) == controlled_ids
        )
        if len(matching_indices) != 1:
            raise RuntimeError("Native GRID components do not match output groups")
        component_index = matching_indices[0]
        component_grid = {
            param_id: values
            for param_id, values in grid.items()
            if param_id in controlled_ids
        }
        objective_by_coordinates: dict[tuple[float, ...], float] = {}
        for attempt in outcome.attempts:
            coordinates = dict(attempt.axis_items)
            key = tuple(
                float(cast("int | float", coordinates[param_id]))
                for param_id in component_grid
            )
            candidate = attempt.components[component_index].candidate
            objective = np.inf if candidate is None else candidate.chi_square
            objective_by_coordinates[key] = min(
                objective,
                objective_by_coordinates.get(key, np.inf),
            )
        coordinate_order = tuple(product(*component_grid.values()))
        chisqr = np.asarray(
            [
                objective_by_coordinates[tuple(float(value) for value in coordinates)]
                for coordinates in coordinate_order
            ],
            dtype=np.float64,
        ).reshape(tuple(len(values) for values in component_grid.values()))
        component_results.append(GridResult(component_grid, chisqr))

        basename = group_path if group_path != Path() else Path("grid")
        filename = grid_path / f"{basename}.out"
        filename.parent.mkdir(parents=True, exist_ok=True)
        with filename.open("w", encoding="utf-8") as output:
            output.write(print_header(component_grid, parameter_store=parameter_store))
            for coordinates, objective in zip(
                coordinate_order,
                chisqr.flat,
                strict=True,
            ):
                output.write(print_values(coordinates, float(objective)))

    combined = combine_grids(grid, component_results)
    grids_1d = make_grids_nd(
        grid,
        combined,
        1,
        parameter_store=parameter_store,
    )
    grids_2d = make_grids_nd(
        grid,
        combined,
        2,
        parameter_store=parameter_store,
    )
    plot_grid_1d(grids_1d, grid_path, parameter_store=parameter_store)
    plot_grid_2d(grids_2d, grid_path, parameter_store=parameter_store)


def _write_direct_output(
    experiments: Experiments,
    path: Path,
    plot: str,
    groups: list[Group],
    group_scopes: tuple[tuple[Path, frozenset[str]], ...],
    result: EvaluationResult,
    variable_count: int,
    *,
    aggregate_output: bool,
) -> None:
    """Write established per-group and aggregate output from native evidence."""
    plot_flg = (plot == "normal" and not aggregate_output) or plot == "all"
    controlled_by_path = dict(group_scopes)
    for group in groups:
        if message := group.message:
            print_group_name(message)
        execute_post_fit(
            group.experiments,
            path / group.path,
            plot=plot_flg,
            residuals=_project_residuals(group.experiments, experiments, result),
            nvarys=len(controlled_by_path[group.path]),
        )
    if aggregate_output:
        execute_post_fit_groups(
            experiments,
            path,
            plot,
            residuals=result.residuals,
            nvarys=variable_count,
        )


def run_native_deterministic(
    experiments: Experiments,
    method: Method,
    path: Path,
    plot: str,
    *,
    session: AnalysisSession,
) -> NativeDeterministicFit | None:
    """Execute one complete deterministic method occurrence natively."""
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    if parameter_model is None or configuration is None:
        raise RuntimeError("Native parameter configuration is unavailable")

    normalized_method = method.model_copy(update={"fitmethod": "trf"})
    parameterization = session.compile_current_parameterization(experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    starting_snapshot = session.analysis_values.snapshot()
    if not _has_controlled_parameters(parameterization):
        workflow = MethodStepWorkflow.for_evaluation(
            starting_snapshot=starting_snapshot,
            parameter_model=parameter_model,
            parameterization=parameterization,
            engine=engine,
            method=normalized_method,
            purpose=EvaluationPurpose.NO_OPTIMIZATION_REQUIRED,
        )
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
        )
        if outcome.lifecycle is not MethodStepLifecycle.SUCCESSFUL_NO_STATE_CHANGE:
            msg = f"Native deterministic evaluation failed: {outcome.lifecycle.value}"
            raise RuntimeError(msg)
        result = cast("EvaluationResult", outcome.evaluation_result)
        session.sync_parameter_store_from_analysis_values(
            dict(result.resolved_values.ordered_items())
        )
        _materialize_evaluation(experiments, result)
        execute_post_fit(
            experiments,
            path,
            plot=plot != "nothing",
            residuals=result.residuals,
            nvarys=0,
        )
        return None

    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        starting_snapshot,
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    all_groups = create_groups(experiments)
    groups = [group for group in all_groups if group.experiments]
    aggregate_output = len(all_groups) > 1
    group_scopes = tuple(
        (
            group.path,
            frozenset(
                param_id
                for param_id, parameter in experiments.parameter_store.get_parameters(
                    group.experiments.param_ids
                ).items()
                if parameter.vary and not parameter.expr
            ),
        )
        for group in groups
    )
    strategy, invocation = _build_invocation(
        method,
        problem,
        decomposition,
        experiments,
    )
    workflow = MethodStepWorkflow.for_optimization(
        starting_snapshot=starting_snapshot,
        parameter_model=parameter_model,
        parameterization=parameterization,
        engine=engine,
        method=normalized_method,
        problem=problem,
        decomposition=decomposition,
        strategy=strategy,
        invocation=invocation,
    )
    print_minimizing()
    outcome = execute_method_step(workflow, analysis_values=session.analysis_values)
    if outcome.lifecycle is not MethodStepLifecycle.COMMITTED:
        msg = (
            "Native deterministic fit did not commit: "
            f"{outcome.primary_terminal or outcome.lifecycle.value}"
        )
        raise RuntimeError(msg)

    accepted = outcome.accepted_result
    if accepted is None:
        raise RuntimeError("Committed native fit lacks its accepted result")
    result = accepted.evaluation_result
    fit = NativeDeterministicFit(accepted, problem, parameterization, engine)
    session.sync_parameter_store_from_analysis_values(
        dict(result.resolved_values.ordered_items())
    )
    _materialize_evaluation(experiments, result)
    variable_count = len(problem.controlled_ids)
    if method.grid:
        _write_grid_output(
            cast("GroupedGridDirectTrfOutcome", outcome.primary_execution),
            cast("GridDirectTrfInvocation", invocation),
            path,
            experiments=experiments,
            group_scopes=group_scopes,
        )
        if aggregate_output:
            execute_post_fit_groups(
                experiments,
                path,
                plot,
                residuals=result.residuals,
                nvarys=variable_count,
            )
        else:
            execute_post_fit(
                experiments,
                path,
                plot=plot != "nothing",
                residuals=result.residuals,
                nvarys=variable_count,
            )
        return fit

    _write_direct_output(
        experiments,
        path,
        plot,
        groups,
        group_scopes,
        result,
        variable_count,
        aggregate_output=aggregate_output,
    )
    return fit
