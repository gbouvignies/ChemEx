"""Production composition for native deterministic method steps."""

from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import cast

import numpy as np

from chemex.configuration.methods import Method
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine, EvaluationResult
from chemex.messages import (
    MinimizationProgressReporter,
    console,
    print_group_name,
    print_minimizing,
)
from chemex.native_provenance import ProvenanceEnvironment
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    DirectTrfInvocation,
    OptimizationProblem,
)
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
    GroupedDirectTrfOutcome,
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
    DerivationDisposition,
    EvaluationPurpose,
    MethodStepLifecycle,
    MethodStepStrategy,
    MethodStepWorkflow,
    UncertaintyDerivationRequest,
    execute_method_step,
)
from chemex.optimize.uncertainty import (
    ParameterUnit,
    ResidualVarianceScaling,
    RootAnchoredBlockCovarianceEvidence,
    UncertaintyConstructionError,
    UncertaintyEvidence,
    UncertaintyPolicy,
    compile_constraint_linearization_capabilities,
    derive_root_anchored_block_covariance,
)
from chemex.parameters.parameterization import ActiveParameterization, ParameterRole
from chemex.printers.parameters import (
    ParameterUncertaintyView,
    parameter_uncertainty_view,
)
from chemex.runtime import AnalysisSession
from chemex.typing import Array

# The finalized #573/#664 policy uses 2000 * (nvars + 1), with numerical
# Jacobian requests counted in the same total objective-request ceiling.
_TRF_OBJECTIVE_REQUESTS_PER_DIMENSION = 2000


@dataclass(frozen=True, slots=True)
class NativeDeterministicFit:
    """Accepted production fit context available to native successor analyses."""

    accepted: AcceptedFitResult
    problem: OptimizationProblem
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    uncertainty: UncertaintyEvidence | None = None
    block_uncertainty: RootAnchoredBlockCovarianceEvidence | None = None

    @property
    def objective_request_budget(self) -> int:
        """Return the accepted occurrence's native Direct TRF request budget."""
        return _objective_request_budget(self.problem)


def _objective_request_budget(problem: OptimizationProblem) -> int:
    coordinate_count = max(1, len(problem.controlled_ids))
    return _TRF_OBJECTIVE_REQUESTS_PER_DIMENSION * (coordinate_count + 1)


def _product_x_scale(problem: OptimizationProblem) -> tuple[float, ...]:
    """Scale native product coordinates without changing physical semantics."""
    return tuple(max(1.0, abs(value)) for value in problem.start)


def _product_uncertainty_request(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
) -> tuple[UncertaintyDerivationRequest, tuple[str, ...]]:
    """Resolve the normal absolute-observation-sigma covariance policy."""
    environment = ProvenanceEnvironment.from_current_process()
    policy = UncertaintyPolicy(
        calibration_identity="native-product-local-covariance-numerics-v1",
        numerical_compatibility_requirement=(
            "binary64-scipy-economical-svd-fixed-pairwise-v1"
        ),
        coordinate_scales=tuple(
            zip(problem.controlled_ids, _product_x_scale(problem), strict=True)
        ),
        coordinate_units=tuple(
            (param_id, ParameterUnit.UNSPECIFIED) for param_id in problem.controlled_ids
        ),
        residual_variance_scaling=(
            ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
        ),
        relative_step_tolerance=1.0e-4,
        roundoff_multiplier=64.0,
        smaller_step_extent=8,
        larger_step_extent=8,
        svd_driver="gesdd",
        rank_absolute_tolerance=0.0,
        rank_relative_tolerance=1.0e-12,
        weak_relative_tolerance=1.0e-6,
        singular_value_cluster_relative_tolerance=1.0e-10,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1",
    )
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    controlled = frozenset(problem.controlled_ids)

    def depends_on_controlled(
        param_id: str,
        visiting: frozenset[str] = frozenset(),
    ) -> bool:
        if param_id in controlled:
            return True
        if param_id in visiting or param_id not in constraints:
            return False
        constraint = constraints[param_id]
        return any(
            depends_on_controlled(dependency, visiting | {param_id})
            for dependency in constraint.dependencies
        )

    propagation_candidates = tuple(
        param_id
        for param_id in parameterization.scope_ids
        if parameterization.role(param_id) is ParameterRole.DERIVED
        and depends_on_controlled(param_id)
    )
    supported: list[str] = []
    unsupported: list[str] = []
    for param_id in propagation_candidates:
        try:
            compile_constraint_linearization_capabilities(
                parameterization,
                (param_id,),
                (),
            )
        except UncertaintyConstructionError:
            unsupported.append(param_id)
        else:
            supported.append(param_id)
    constrained_scope = tuple(supported)
    capabilities = compile_constraint_linearization_capabilities(
        parameterization,
        constrained_scope,
        (),
    )
    return (
        UncertaintyDerivationRequest(
            policy,
            constrained_scope=constrained_scope,
            constrained_units=tuple(
                (param_id, ParameterUnit.UNSPECIFIED) for param_id in constrained_scope
            ),
            constrained_scales=tuple((param_id, 1.0) for param_id in constrained_scope),
            compiled_capabilities=capabilities,
            resolved_environment_identity=environment.identity,
        ),
        tuple(unsupported),
    )


def _product_uncertainty_result(
    outcome: object,
    controlled_ids: tuple[str, ...],
    unsupported_constrained_ids: tuple[str, ...],
) -> tuple[
    UncertaintyEvidence | None,
    ParameterUncertaintyView,
    tuple[str, str] | None,
]:
    derivations = getattr(outcome, "derivations", ())
    for derivation in derivations:
        for artifact in derivation.artifacts:
            if isinstance(artifact, UncertaintyEvidence):
                return (
                    artifact,
                    parameter_uncertainty_view(
                        artifact,
                        unsupported_constrained_ids,
                    ),
                    None,
                )
        if derivation.stage == "uncertainty" and derivation.disposition in {
            DerivationDisposition.CANCELLED,
            DerivationDisposition.INTERRUPTED,
            DerivationDisposition.NOT_STARTED_BY_WORKFLOW_STOP,
        }:
            terminal = derivation.disposition.value
            reason = (
                "derivation interrupted"
                if derivation.disposition is DerivationDisposition.INTERRUPTED
                else "derivation cancelled"
            )
            view = ParameterUncertaintyView(
                unavailable_reasons=(
                    *((param_id, reason) for param_id in controlled_ids),
                    *(
                        (param_id, "unsupported constrained derivative")
                        for param_id in unsupported_constrained_ids
                    ),
                )
            )
            return None, view, (terminal, reason)
    raise RuntimeError("Committed native fit lacks requested uncertainty evidence")


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
            x_scale=_product_x_scale(problem),
        )
        return MethodStepStrategy.GRID_DIRECT_TRF, invocation
    invocation = GroupedDirectTrfInvocation(
        decomposition.root_problem_identity,
        decomposition.identity,
        tuple(
            DirectTrfInvocation.for_problem(
                component.problem,
                objective_request_budget=_objective_request_budget(component.problem),
                x_scale=_product_x_scale(component.problem),
            )
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
    uncertainty: ParameterUncertaintyView,
    uncertainty_evidence: UncertaintyEvidence | None,
    uncertainty_status: tuple[str, str] | None,
    block_uncertainty: RootAnchoredBlockCovarianceEvidence | None,
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
            uncertainty=uncertainty,
            uncertainty_evidence=uncertainty_evidence,
            uncertainty_status=uncertainty_status,
            block_uncertainty=block_uncertainty,
        )
    if aggregate_output:
        execute_post_fit_groups(
            experiments,
            path,
            plot,
            residuals=result.residuals,
            nvarys=variable_count,
            uncertainty=uncertainty,
            uncertainty_evidence=uncertainty_evidence,
            uncertainty_status=uncertainty_status,
            block_uncertainty=block_uncertainty,
        )


def run_native_deterministic(
    experiments: Experiments,
    method: Method,
    path: Path,
    plot: str,
    *,
    session: AnalysisSession,
    step_name: str = "",
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
    uncertainty_request, unsupported_constrained_ids = _product_uncertainty_request(
        problem,
        parameterization,
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
        derivations=(uncertainty_request,),
    )
    group_labels = {
        controlled_ids: (
            group_path.name.split("_", maxsplit=1)[-1] if group_path.name else ""
        )
        for group_path, controlled_ids in group_scopes
    }
    progress = MinimizationProgressReporter(
        console,
        interactive=console.is_terminal,
        retained_observation_count=engine.plan.retained_observation_count,
        controlled_parameter_count=len(problem.controlled_ids),
        step_name=step_name,
        group_labels=group_labels,
        grid=bool(method.grid),
    )
    print_minimizing()
    with progress:
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
            progress_observer=progress.observe,
        )
        accepted_result = outcome.accepted_result
        progress.finish(
            final_chi_square=(
                None if accepted_result is None else accepted_result.chi_square
            ),
            terminal_status=outcome.lifecycle.value,
        )
    if outcome.lifecycle is MethodStepLifecycle.INTERRUPTED:
        raise KeyboardInterrupt("Native deterministic fit interrupted")
    if outcome.lifecycle is not MethodStepLifecycle.COMMITTED:
        primary = outcome.primary_execution
        failures = (
            tuple(
                (component.controlled_ids, component.failure)
                for component in primary.components
                if component.failure is not None
            )
            if isinstance(primary, GroupedDirectTrfOutcome)
            else ()
        )
        detail = (
            ""
            if not failures
            else ": "
            + "; ".join(
                f"{controlled_ids!r}: {failure.category}: {failure.message}"
                for controlled_ids, failure in failures
            )
        )
        msg = (
            "Native deterministic fit did not commit: "
            f"{outcome.primary_terminal or outcome.lifecycle.value}{detail}"
        )
        raise RuntimeError(msg)

    accepted = outcome.accepted_result
    if accepted is None:
        raise RuntimeError("Committed native fit lacks its accepted result")
    result = accepted.evaluation_result
    uncertainty_evidence, uncertainty, uncertainty_status = _product_uncertainty_result(
        outcome,
        problem.controlled_ids,
        unsupported_constrained_ids,
    )
    block_uncertainty = (
        None
        if uncertainty_evidence is None
        else derive_root_anchored_block_covariance(
            uncertainty_evidence,
            tuple(component.controlled_ids for component in decomposition.components),
        )
    )
    if block_uncertainty is not None:
        block_errors = tuple(
            (entry.param_id, entry.value)
            for block in block_uncertainty.blocks
            for entry in block.marginal_errors
            if entry.value is not None
        )
        block_unavailable = tuple(
            (param_id, block.unavailable_reason)
            for block in block_uncertainty.blocks
            if block.unavailable_reason
            for param_id in block.controlled_ids
        )
        uncertainty = ParameterUncertaintyView(
            uncertainty.standard_errors + block_errors,
            uncertainty.unavailable_reasons + block_unavailable,
        )
    fit = NativeDeterministicFit(
        accepted,
        problem,
        parameterization,
        engine,
        uncertainty_evidence,
        block_uncertainty,
    )
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
                uncertainty=uncertainty,
                uncertainty_evidence=uncertainty_evidence,
                uncertainty_status=uncertainty_status,
                block_uncertainty=block_uncertainty,
            )
        else:
            execute_post_fit(
                experiments,
                path,
                plot=plot != "nothing",
                residuals=result.residuals,
                nvarys=variable_count,
                uncertainty=uncertainty,
                uncertainty_evidence=uncertainty_evidence,
                uncertainty_status=uncertainty_status,
                block_uncertainty=block_uncertainty,
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
        uncertainty=uncertainty,
        uncertainty_evidence=uncertainty_evidence,
        uncertainty_status=uncertainty_status,
        block_uncertainty=block_uncertainty,
    )
    return fit
