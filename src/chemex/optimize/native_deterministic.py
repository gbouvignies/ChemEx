"""Production composition for native deterministic method steps."""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import reduce
from itertools import product
from operator import and_
from pathlib import Path
from typing import cast

import numpy as np

from chemex.configuration.methods import Method
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationResult,
    ResolvedEvaluationValues,
)
from chemex.messages import (
    MinimizationProgressReporter,
    UncertaintyProgressReporter,
    console,
    print_minimizing,
)
from chemex.native_provenance import ProvenanceEnvironment
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    DirectTrfInvocation,
    FitCommitOperation,
    FitCommitTerminal,
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
from chemex.optimize.helper import (
    execute_post_fit,
    print_header,
    print_values,
)
from chemex.optimize.method_step import (
    DerivationDisposition,
    EvaluationPurpose,
    MethodStepLifecycle,
    MethodStepOutcome,
    MethodStepStrategy,
    MethodStepWorkflow,
    UncertaintyDerivationRequest,
    execute_method_step,
)
from chemex.optimize.uncertainty import (
    MissingFunctionLinearizationCapability,
    ParameterUnit,
    ResidualVarianceScaling,
    RootAnchoredBlockCovarianceEvidence,
    UncertaintyEvidence,
    UncertaintyPolicy,
    UncertaintyUnavailableKind,
    compile_constraint_linearization_capabilities,
    derive_root_anchored_block_covariance,
)
from chemex.parameters.database import ParameterStore
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.printers.parameters import (
    BOUNDARY_WARNING_TEXT,
    ParameterUncertaintyView,
    constrained_boundary_warning_ids,
    parameter_uncertainty_view,
    uncertainty_unavailable_reason,
)
from chemex.runtime import AnalysisSession

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
    parameter_model: SealedParameterModel
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
        calibration_identity="native-product-local-covariance-numerics-v2",
        numerical_compatibility_requirement=(
            "binary64-retained-scipy-2point-gesdd-column-equilibrated-v1"
        ),
        coordinate_scales=tuple((param_id, 1.0) for param_id in problem.controlled_ids),
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
        # Legacy qualification knobs; covariance rank uses sigma_max*max(m,n)*eps.
        rank_absolute_tolerance=0.0,
        rank_relative_tolerance=0.0,
        weak_relative_tolerance=1.0e-6,
        singular_value_cluster_relative_tolerance=1.0e-10,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1",
        residual_jacobian_strategy="retained-backend-or-accepted-2-point",
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
        except MissingFunctionLinearizationCapability:
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
    outcome: MethodStepOutcome,
    controlled_ids: tuple[str, ...],
    supported_constrained_ids: tuple[str, ...],
    unsupported_constrained_ids: tuple[str, ...],
) -> tuple[
    UncertaintyEvidence | None,
    ParameterUncertaintyView,
    tuple[str, str] | None,
]:
    for derivation in outcome.derivations:
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
            reason = uncertainty_unavailable_reason(
                UncertaintyUnavailableKind.DERIVATION_STOPPED
            )
            view = ParameterUncertaintyView(
                unavailable_reasons=(
                    *((param_id, reason) for param_id in controlled_ids),
                    *((param_id, reason) for param_id in supported_constrained_ids),
                    *(
                        (
                            param_id,
                            uncertainty_unavailable_reason(
                                UncertaintyUnavailableKind.UNSUPPORTED_CONSTRAINED_DERIVATIVE
                            ),
                        )
                        for param_id in unsupported_constrained_ids
                    ),
                )
            )
            return None, view, (terminal, reason)
    raise RuntimeError("Committed native fit lacks requested uncertainty evidence")


def _product_block_uncertainty(
    evidence: UncertaintyEvidence | None,
    decomposition: FitDecomposition,
) -> tuple[
    RootAnchoredBlockCovarianceEvidence | None,
    tuple[str, str] | None,
]:
    """Derive proof-bound block evidence without invalidating a committed fit."""
    if evidence is None:
        return None, None
    try:
        return (
            derive_root_anchored_block_covariance(
                evidence,
                decomposition.partition_proof,
            ),
            None,
        )
    except KeyboardInterrupt:
        return None, (
            "interrupted",
            uncertainty_unavailable_reason(
                UncertaintyUnavailableKind.DERIVATION_STOPPED
            ),
        )


def _uncertainty_progress_status(
    evidence: UncertaintyEvidence | None,
    block_evidence: RootAnchoredBlockCovarianceEvidence | None,
    uncertainty: ParameterUncertaintyView,
    controlled_ids: tuple[str, ...],
) -> str:
    """Summarize fitted covariance availability without overstating constraints."""
    root_boundary_warning = (
        evidence is not None
        and evidence.covariance is not None
        and evidence.covariance.boundary_warning
    )
    block_boundary_warning = block_evidence is not None and any(
        block.boundary_warning for block in block_evidence.blocks
    )
    if root_boundary_warning or block_boundary_warning:
        return "covariance available with boundary warnings"
    available_ids = {param_id for param_id, _value in uncertainty.standard_errors}
    if (
        evidence is not None
        and evidence.covariance is not None
        and all(param_id in available_ids for param_id in controlled_ids)
    ):
        return "covariance available"
    if block_evidence is not None and any(
        block.covariance is not None for block in block_evidence.blocks
    ):
        return "covariance partially available"
    reasons = dict(uncertainty.unavailable_reasons)
    reason = next(
        (reasons[param_id] for param_id in controlled_ids if param_id in reasons),
        "insufficient information",
    )
    return f"uncertainty unavailable: {reason}"


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


def _execute_with_phase_progress(
    workflow: MethodStepWorkflow,
    session: AnalysisSession,
    progress: MinimizationProgressReporter,
    uncertainty_progress: UncertaintyProgressReporter,
) -> MethodStepOutcome:
    """Close minimization at commit and time the following uncertainty phase."""
    minimization_finished = False

    def finish_minimization_at_commit(
        accepted: AcceptedFitResult,
        operation: FitCommitOperation,
    ) -> None:
        nonlocal minimization_finished
        progress.finish(
            final_chi_square=accepted.chi_square,
            terminal_status=operation.terminal.value,
        )
        minimization_finished = True
        if operation.terminal is FitCommitTerminal.COMMITTED:
            uncertainty_progress.start()

    print_minimizing()
    with progress:
        outcome = execute_method_step(
            workflow,
            analysis_values=session.analysis_values,
            progress_observer=progress.observe,
            commit_completed_observer=finish_minimization_at_commit,
        )
        if not minimization_finished:
            accepted_result = outcome.accepted_result
            progress.finish(
                final_chi_square=(
                    None if accepted_result is None else accepted_result.chi_square
                ),
                terminal_status=outcome.lifecycle.value,
            )
    return outcome


def _build_invocation(
    method: Method,
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    parameter_store: ParameterStore,
) -> tuple[
    MethodStepStrategy,
    GroupedDirectTrfInvocation | GridDirectTrfInvocation,
]:
    if method.grid:
        grid = parameter_store.parse_grid(method.grid)
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
    parameter_model: SealedParameterModel,
    accepted_values: Mapping[str, float],
) -> None:
    """Render aggregate GRID artifacts from native seed evidence."""
    grid = {
        axis.param_id: np.asarray(axis.values, dtype=np.float64)
        for axis in invocation.axes
    }
    parameter_names = {
        definition.param_id: parameter_name_from_definition(definition)
        for definition in parameter_model.definitions
    }
    grid_path = path / "Grid"
    grid_path.mkdir(parents=True, exist_ok=True)
    component_results: list[GridResult] = []
    for component_index, component in enumerate(outcome.attempts[0].components):
        controlled_ids = frozenset(component.controlled_ids)
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

    with (grid_path / "grid.out").open("w", encoding="utf-8") as output:
        output.write(print_header(grid, parameter_names=parameter_names))
        for attempt in outcome.attempts:
            coordinates = (
                float(cast("int | float", value))
                for _param_id, value in attempt.axis_items
            )
            objective = np.inf if attempt.objective is None else attempt.objective
            output.write(print_values(coordinates, objective))

    combined = combine_grids(grid, component_results)
    grids_1d = make_grids_nd(
        grid,
        combined,
        1,
        parameter_names=parameter_names,
    )
    grids_2d = make_grids_nd(
        grid,
        combined,
        2,
        parameter_names=parameter_names,
    )
    plot_grid_1d(
        grids_1d,
        grid_path,
        parameter_names=parameter_names,
        accepted_values=accepted_values,
    )
    plot_grid_2d(
        grids_2d,
        grid_path,
        parameter_names=parameter_names,
        accepted_values=accepted_values,
    )


def _fit_component_labels(
    decomposition: FitDecomposition,
    parameter_model: SealedParameterModel,
) -> dict[frozenset[str], str]:
    """Return concise optional labels for transient component presentation."""
    try:
        labels: dict[frozenset[str], str] = {}
        for component in decomposition.components:
            param_names = tuple(
                parameter_name_from_definition(parameter_model.definitions[param_id])
                for param_id in component.controlled_ids
            )
            if param_names:
                labels[frozenset(component.controlled_ids)] = reduce(
                    and_, param_names
                ).folder
    except Exception:  # noqa: BLE001 - optional reporting is non-scientific
        return {}
    else:
        return labels


def _commit_resolved_continuity_if_changed(
    analysis_values: AnalysisValues,
    starting_snapshot: AnalysisValuesSnapshot,
    resolved_values: ResolvedEvaluationValues,
) -> None:
    """Commit a successful derived-only transition when its values changed."""
    resolved_items = resolved_values.ordered_items()
    if not any(
        value != starting_snapshot[param_id] for param_id, value in resolved_items
    ):
        return
    analysis_values.commit(
        dict(resolved_items),
        expected=starting_snapshot,
        scope=tuple(resolved_values),
    )


def run_native_deterministic(
    experiments: Experiments,
    method: Method,
    path: Path,
    plot: str,
    *,
    session: AnalysisSession,
    parameterization: ActiveParameterization | None = None,
) -> NativeDeterministicFit | None:
    """Execute one complete deterministic method occurrence natively."""
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    if parameter_model is None or configuration is None:
        raise RuntimeError("Native parameter configuration is unavailable")

    normalized_method = method.model_copy(update={"fitmethod": "trf"})
    if parameterization is None:
        parameterization = session.compile_parameterization(
            method,
            experiments.param_ids,
        )
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
        _commit_resolved_continuity_if_changed(
            session.analysis_values,
            starting_snapshot,
            result.resolved_values,
        )
        _materialize_evaluation(experiments, result)
        execute_post_fit(
            experiments,
            path,
            plot=plot != "nothing",
            residuals=result.residuals,
            nvarys=0,
            parameter_model=parameter_model,
            parameter_values=result.resolved_values,
            parameterization=parameterization,
        )
        return None

    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        starting_snapshot,
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    strategy, invocation = _build_invocation(
        method,
        problem,
        decomposition,
        session.parameters,
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
    component_labels = _fit_component_labels(decomposition, parameter_model)
    progress = MinimizationProgressReporter(
        console,
        interactive=console.is_terminal,
        retained_observation_count=engine.plan.retained_observation_count,
        controlled_parameter_count=len(problem.controlled_ids),
        component_labels=component_labels,
        grid=bool(method.grid),
    )
    uncertainty_progress = UncertaintyProgressReporter(console)
    outcome = _execute_with_phase_progress(
        workflow,
        session,
        progress,
        uncertainty_progress,
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
        uncertainty_request.constrained_scope,
        unsupported_constrained_ids,
    )
    block_uncertainty, block_status = _product_block_uncertainty(
        uncertainty_evidence,
        decomposition,
    )
    uncertainty_status = block_status or uncertainty_status
    if block_uncertainty is not None:
        block_errors = tuple(
            (entry.param_id, entry.value)
            for block in block_uncertainty.blocks
            for entry in block.marginal_errors
            if entry.value is not None
        )
        block_unavailable = tuple(
            (param_id, uncertainty_unavailable_reason(block.unavailable_kind))
            for block in block_uncertainty.blocks
            if block.unavailable_kind is not None
            for param_id in block.controlled_ids
        )
        block_constrained_errors = tuple(
            (entry.param_id, entry.value)
            for block in block_uncertainty.blocks
            for entry in block.constrained_marginal_errors
            if entry.value is not None
        )
        block_controlled_warning_ids = block_uncertainty.simple_bound_warning_ids
        block_constrained_warning_ids = constrained_boundary_warning_ids(
            block_uncertainty.source_bundle,
            block_controlled_warning_ids,
        )
        block_warnings = tuple(
            (entry.param_id, BOUNDARY_WARNING_TEXT)
            for block in block_uncertainty.blocks
            for entry in block.marginal_errors
            if entry.value is not None
            and entry.param_id in block_controlled_warning_ids
        ) + tuple(
            (entry.param_id, BOUNDARY_WARNING_TEXT)
            for block in block_uncertainty.blocks
            for entry in block.constrained_marginal_errors
            if entry.value is not None
            and entry.param_id in block_constrained_warning_ids
        )
        recovered_ids = {
            param_id for param_id, _value in (*block_errors, *block_constrained_errors)
        }
        uncertainty = ParameterUncertaintyView(
            standard_errors=(
                uncertainty.standard_errors + block_errors + block_constrained_errors
            ),
            warnings=uncertainty.warnings + block_warnings,
            unavailable_reasons=(
                tuple(
                    (param_id, reason)
                    for param_id, reason in uncertainty.unavailable_reasons
                    if param_id not in recovered_ids
                )
                + block_unavailable
            ),
        )
    uncertainty_progress.finish(
        _uncertainty_progress_status(
            uncertainty_evidence,
            block_uncertainty,
            uncertainty,
            problem.controlled_ids,
        )
    )
    fit = NativeDeterministicFit(
        accepted,
        problem,
        parameterization,
        engine,
        parameter_model,
        uncertainty_evidence,
        block_uncertainty,
    )
    _materialize_evaluation(experiments, result)
    variable_count = len(problem.controlled_ids)
    if method.grid:
        _write_grid_output(
            cast("GroupedGridDirectTrfOutcome", outcome.primary_execution),
            cast("GridDirectTrfInvocation", invocation),
            path,
            parameter_model=parameter_model,
            accepted_values=result.resolved_values,
        )
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
            parameter_model=parameter_model,
            parameter_values=result.resolved_values,
            parameterization=parameterization,
            fitted_ids=problem.controlled_ids,
        )
        return fit

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
        parameter_model=parameter_model,
        parameter_values=result.resolved_values,
        parameterization=parameterization,
        fitted_ids=problem.controlled_ids,
    )
    return fit
