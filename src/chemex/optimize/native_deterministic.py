"""Production composition for deterministic evaluation and fitting."""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import reduce
from operator import and_
from pathlib import Path
from typing import Never

import numpy as np

from chemex.configuration.method_plan import DeSearch, GridSearch, SearchScale
from chemex.configuration.method_validation import (
    resolve_de_coordinates,
    resolve_grid_axes,
)
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
    ResolvedEvaluationValues,
)
from chemex.exceptions import ArtifactPublicationError, ChemExError
from chemex.messages import (
    MinimizationProgressReporter,
    UncertaintyProgressReporter,
    console,
    print_minimizing,
    print_running_de,
)
from chemex.optimize.de_direct_trf import (
    DeCoordinateSemantics,
    DeSearchExecution,
    DeSearchInvocation,
    DeSearchTerminal,
    execute_de_search,
)
from chemex.optimize.deterministic_uncertainty import (
    AcceptedDeterministicFitFacts,
    ContinuousTrfBasis,
    DerivationDisposition,
    DeterministicUncertainty,
    ProfiledGridBasis,
    derive_deterministic_uncertainty,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    DirectTrfInvocation,
    FitCommitOperation,
    FitCommitTerminal,
    OptimizationProblem,
    TerminalFailure,
    build_bounded_independent_frame,
    execute_fit_commit,
)
from chemex.optimize.grouped_direct_trf import (
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfOutcome,
    execute_grouped_direct_trf,
)
from chemex.optimize.helper import execute_post_fit
from chemex.optimize.profiled_grid import (
    ProfiledGridOutcome,
    execute_profiled_grid,
)
from chemex.parameters.feasible_coordinates import validate_relaxation_state
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
from chemex.printers.grid import write_grid_output
from chemex.printers.parameters import (
    uncertainty_progress_status,
)
from chemex.run_info import RunInfo, mark_failure_stage
from chemex.runtime import AnalysisSession
from chemex.runtime.environment import RuntimeEnvironment

# The finalized #573/#664 policy uses 2000 * (nvars + 1), with numerical
# Jacobian requests counted in the same total objective-request ceiling.
_TRF_OBJECTIVE_REQUESTS_PER_DIMENSION = 2000


type DeterministicTerminalRecord = (
    DeSearchExecution
    | EvaluationFailure
    | FitCommitOperation
    | GroupedDirectTrfOutcome
    | ProfiledGridOutcome
)


class NativeDeterministicAnalysisError(ChemExError, RuntimeError):
    """Recognized deterministic terminal state with no committed fit."""

    def __init__(
        self,
        message: str,
        *,
        reason: str,
        outcome: DeterministicTerminalRecord,
        failures: tuple[TerminalFailure, ...],
    ) -> None:
        super().__init__(message)
        self.terminal = "failed"
        self.reason = reason
        self.outcome = outcome
        self.failures = failures


class NativeDeterministicInternalError(RuntimeError):
    """Unexpected deterministic failure retaining its terminal evidence."""

    def __init__(
        self,
        message: str,
        *,
        outcome: DeterministicTerminalRecord,
        failures: tuple[TerminalFailure, ...] = (),
    ) -> None:
        super().__init__(message)
        self.outcome = outcome
        self.failures = failures


def _is_known_analysis_failure(failure: TerminalFailure) -> bool:
    evaluation_failure = failure.evaluation_failure
    return failure.category in {
        "non_converged",
        "objective_budget_exhausted",
        "objective_scalarization_failure",
    } or (
        evaluation_failure is not None
        and evaluation_failure.validity == "INVALID_TRIAL"
    )


def _known_grouped_analysis_error(
    outcome: GroupedDirectTrfOutcome,
) -> NativeDeterministicAnalysisError | None:
    failures = tuple(
        component.failure
        for component in outcome.components
        if component.failure is not None
    )
    if (
        outcome.terminal.value != "component_failure"
        or not failures
        or not all(_is_known_analysis_failure(failure) for failure in failures)
    ):
        return None
    categories = {failure.category for failure in failures}
    if "objective_budget_exhausted" in categories:
        reason = "The objective-evaluation budget was exhausted."
    elif categories == {"non_converged"}:
        reason = "The optimizer did not converge."
    else:
        reason = "The optimizer encountered an invalid numerical trial."
    return NativeDeterministicAnalysisError(
        "Native deterministic fitting did not produce an accepted fit.",
        reason=reason,
        outcome=outcome,
        failures=failures,
    )


def _propagate_uncertainty_interruption(
    finalization_error: BaseException | None = None,
) -> Never:
    error = KeyboardInterrupt("Native deterministic uncertainty derivation interrupted")
    mark_failure_stage(error, "uncertainty")
    if finalization_error is not None:
        error.add_note("ChemEx could not publish interrupted uncertainty output.")
        raise error from finalization_error
    raise error


def _propagate_output_failure(
    error: BaseException,
    path: Path,
    *,
    interrupted: bool = False,
) -> Never:
    """Add deterministic output context without reclassifying non-I/O failures."""
    if isinstance(error, ArtifactPublicationError):
        if interrupted:
            object.__setattr__(error, "terminal", "interrupted")
        raise error
    if isinstance(error, OSError):
        filename = error.filename
        destination = Path(filename) if isinstance(filename, str) else path
        failure = ArtifactPublicationError(
            "publish deterministic analysis output",
            destination,
            error,
        )
        mark_failure_stage(failure, "output")
        if interrupted:
            object.__setattr__(failure, "terminal", "interrupted")
        raise failure from error
    mark_failure_stage(error, "output")
    raise error


@dataclass(frozen=True, slots=True)
class NativeDeterministicFit:
    """Accepted production fit context available to native successor analyses."""

    accepted: AcceptedFitResult
    problem: OptimizationProblem
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    parameter_model: SealedParameterModel
    deterministic_uncertainty: DeterministicUncertainty

    def __post_init__(self) -> None:
        if self.deterministic_uncertainty.accepted_anchor is not self.accepted:
            raise ValueError(
                "Deterministic uncertainty belongs to a different accepted fit"
            )

    @property
    def objective_request_budget(self) -> int:
        """Return the accepted occurrence's native Direct TRF request budget."""
        return _objective_request_budget(self.problem)


def _objective_request_budget(problem: OptimizationProblem) -> int:
    coordinate_count = max(1, len(problem.controlled_ids))
    return _TRF_OBJECTIVE_REQUESTS_PER_DIMENSION * (coordinate_count + 1)


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


def _execute_and_commit_aggregate(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation | None,
    grid_axes: Mapping[str, tuple[float, ...]] | None,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    analysis_values: AnalysisValues,
    progress: MinimizationProgressReporter,
    uncertainty_progress: UncertaintyProgressReporter,
) -> tuple[
    GroupedDirectTrfOutcome | ProfiledGridOutcome,
    FitCommitOperation | None,
]:
    """Execute components, accept one fresh aggregate, and atomically commit it."""
    token = CancellationToken()
    print_minimizing()
    with progress:
        if invocation is not None:
            outcome: GroupedDirectTrfOutcome | ProfiledGridOutcome = (
                execute_grouped_direct_trf(
                    problem,
                    decomposition,
                    invocation,
                    parameterization,
                    engine,
                    cancellation=token,
                    progress_observer=progress.observe,
                )
            )
        else:
            if grid_axes is None:
                raise RuntimeError("Profiled GRID execution lacks resolved axes")
            outcome = execute_profiled_grid(
                problem,
                grid_axes,
                parameterization,
                engine,
                objective_request_budget=_objective_request_budget(problem),
                cancellation=token,
                progress_observer=progress.observe,
            )
        accepted = outcome.accepted_result
        authority = outcome.commit_authority
        if accepted is None or authority is None:
            progress.finish(
                final_chi_square=None,
                terminal_status=outcome.terminal.value,
            )
            return outcome, None
        if token.is_cancelled:
            progress.finish(
                final_chi_square=accepted.chi_square,
                terminal_status="cancelled",
            )
            raise KeyboardInterrupt("Native deterministic fit cancelled before commit")
        operation = execute_fit_commit(
            accepted,
            authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=analysis_values,
        )
        progress.finish(
            final_chi_square=accepted.chi_square,
            terminal_status=operation.terminal.value,
        )
    if operation.terminal is FitCommitTerminal.COMMITTED and grid_axes is None:
        uncertainty_progress.start()
    return outcome, operation


def _build_invocation(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
) -> GroupedDirectTrfInvocation:
    invocation = GroupedDirectTrfInvocation(
        decomposition.root_problem_identity,
        decomposition.identity,
        tuple(
            DirectTrfInvocation.for_problem(
                component.problem,
                objective_request_budget=_objective_request_budget(component.problem),
            )
            for component in decomposition.components
        ),
    )
    return invocation


def _run_product_de_search(
    search: DeSearch,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    parameter_model: SealedParameterModel,
) -> OptimizationProblem:
    resolved_coordinates = resolve_de_coordinates(search, parameter_model)
    coordinates = tuple(
        (
            coordinate.param_id,
            coordinate.low,
            coordinate.high,
            (
                DeCoordinateSemantics.LINEAR
                if coordinate.scale is SearchScale.LINEAR
                else DeCoordinateSemantics.LOG
            ),
        )
        for coordinate in resolved_coordinates
    )
    invocation = DeSearchInvocation.for_product_problem(
        problem,
        search_coordinates=coordinates,
        root_seed=search.seed,
    )
    print_running_de()
    outcome = execute_de_search(problem, invocation, parameterization, engine)
    if outcome.terminal is DeSearchTerminal.INTERRUPTED:
        raise KeyboardInterrupt("Selected-coordinate DE search interrupted")
    candidate = outcome.best_candidate
    if not outcome.restart_eligible or candidate is None:
        failure = outcome.failure
        if (
            outcome.terminal is DeSearchTerminal.BUDGET_EXHAUSTED
            and failure is not None
        ):
            raise NativeDeterministicAnalysisError(
                "Selected-coordinate DE search produced no eligible candidate.",
                reason="The objective-evaluation budget was exhausted.",
                outcome=outcome,
                failures=(failure,),
            )
        raise NativeDeterministicInternalError(
            "Selected-coordinate DE search produced no eligible candidate",
            outcome=outcome,
            failures=() if failure is None else (failure,),
        )
    return problem.restart_from(candidate.full_vector)


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
) -> AnalysisValuesSnapshot | None:
    """Commit a successful derived-only transition when its values changed."""
    resolved_items = resolved_values.ordered_items()
    if not any(
        value != starting_snapshot[param_id] for param_id, value in resolved_items
    ):
        return None
    return analysis_values.commit(
        dict(resolved_items),
        expected=starting_snapshot,
        scope=tuple(resolved_values),
    )


def run_native_deterministic(  # noqa: C901 - closed Direct/GRID/DE product dispatcher
    experiments: Experiments,
    path: Path,
    plot: str,
    *,
    session: AnalysisSession,
    parameterization: ActiveParameterization,
    search: GridSearch | DeSearch | None = None,
    run_info: RunInfo | None = None,
    step_name: str = "DEFAULT",
) -> NativeDeterministicFit | None:
    """Execute one complete deterministic method occurrence natively."""
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    if parameter_model is None or configuration is None:
        raise RuntimeError("Native parameter configuration is unavailable")

    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    starting_snapshot = session.analysis_values.snapshot()
    if not _has_controlled_parameters(parameterization):
        independent_frame = build_bounded_independent_frame(
            parameterization,
            configuration,
            starting_snapshot,
        )
        validate_relaxation_state(
            parameterization,
            parameterization.resolve(independent_frame),
        )
        frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            independent_frame,
        )
        evaluated = engine.new_evaluator().evaluate(frame)
        if isinstance(evaluated, EvaluationFailure):
            if evaluated.validity == "INVALID_TRIAL":
                failure = TerminalFailure(
                    evaluated.category,
                    evaluated.message,
                    evaluated,
                )
                raise NativeDeterministicAnalysisError(
                    "Native deterministic evaluation did not produce a valid result.",
                    reason="The requested parameter state produced an invalid numerical trial.",
                    outcome=evaluated,
                    failures=(failure,),
                )
            failure = TerminalFailure(
                evaluated.category,
                evaluated.message,
                evaluated,
            )
            raise NativeDeterministicInternalError(
                "Native deterministic evaluation failed",
                outcome=evaluated,
                failures=(failure,),
            )
        result = evaluated
        continuity_snapshot = _commit_resolved_continuity_if_changed(
            session.analysis_values,
            starting_snapshot,
            result.resolved_values,
        )
        if continuity_snapshot is not None and run_info is not None:
            run_info.publish_restart(continuity_snapshot)
        _materialize_evaluation(experiments, result)
        try:
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
        except (Exception, KeyboardInterrupt) as error:  # noqa: BLE001
            _propagate_output_failure(error, path)
        return None

    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        starting_snapshot,
    )
    if isinstance(search, DeSearch):
        if run_info is not None:
            run_info.record_stochastic_operation(step_name, "de", search.seed)
        problem = _run_product_de_search(
            search,
            problem,
            parameterization,
            engine,
            parameter_model,
        )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    if isinstance(search, GridSearch):
        resolved_grid_axes = resolve_grid_axes(
            search,
            parameter_model,
            active_scope_ids=parameterization.scope_ids,
            final_fit_ids=problem.controlled_ids,
        )
        grid_axes: Mapping[str, tuple[float, ...]] | None = {
            axis.param_id: axis.values for axis in resolved_grid_axes
        }
        invocation = None
    else:
        grid_axes = None
        invocation = _build_invocation(problem, decomposition)
    resolved_environment_identity = RuntimeEnvironment.from_current_process().identity
    component_labels = _fit_component_labels(decomposition, parameter_model)
    progress = MinimizationProgressReporter(
        console,
        interactive=console.is_terminal,
        retained_observation_count=engine.plan.retained_observation_count,
        controlled_parameter_count=len(problem.controlled_ids),
        component_labels=component_labels,
        grid=isinstance(search, GridSearch),
    )
    uncertainty_progress = UncertaintyProgressReporter(console)
    outcome, commit = _execute_and_commit_aggregate(
        problem,
        decomposition,
        invocation,
        grid_axes,
        parameterization,
        engine,
        session.analysis_values,
        progress,
        uncertainty_progress,
    )
    if outcome.terminal.value in {"cancelled", "interrupted"}:
        raise KeyboardInterrupt("Native deterministic fit interrupted")
    if commit is None:
        if isinstance(outcome, GroupedDirectTrfOutcome):
            known_error = _known_grouped_analysis_error(outcome)
            if known_error is not None:
                raise known_error
        failures = (
            tuple(
                component.failure
                for component in outcome.components
                if component.failure is not None
            )
            if isinstance(outcome, GroupedDirectTrfOutcome)
            else ()
        )
        if isinstance(outcome, ProfiledGridOutcome) and outcome.failure is not None:
            failures = (outcome.failure,)
        raise NativeDeterministicInternalError(
            f"Native deterministic fit did not commit: {outcome.terminal.value}",
            outcome=outcome,
            failures=failures,
        )
    if commit.terminal is not FitCommitTerminal.COMMITTED:
        failure = commit.failure
        raise NativeDeterministicInternalError(
            "Native deterministic fit commit failed",
            outcome=commit,
        )
    committed_snapshot = commit.committed_snapshot
    if committed_snapshot is None:
        raise RuntimeError("Committed native fit lacks its central-value snapshot")
    if run_info is not None:
        run_info.publish_restart(committed_snapshot)

    accepted = outcome.accepted_result
    if accepted is None:
        raise RuntimeError("Committed native fit lacks its accepted result")
    result = accepted.evaluation_result
    uncertainty_facts = AcceptedDeterministicFitFacts(
        accepted,
        problem,
        parameterization,
        engine,
        (
            ProfiledGridBasis()
            if isinstance(search, GridSearch)
            else ContinuousTrfBasis(decomposition.partition_proof)
        ),
        resolved_environment_identity,
    )
    try:
        deterministic_uncertainty = derive_deterministic_uncertainty(uncertainty_facts)
    except (Exception, KeyboardInterrupt) as error:
        mark_failure_stage(error, "uncertainty")
        raise
    if deterministic_uncertainty.disposition is DerivationDisposition.WITHHELD:
        uncertainty_progress.skip_grid()
    else:
        uncertainty_progress.finish(
            uncertainty_progress_status(deterministic_uncertainty)
        )
    fit = NativeDeterministicFit(
        accepted,
        problem,
        parameterization,
        engine,
        parameter_model,
        deterministic_uncertainty,
    )
    uncertainty_interrupted = deterministic_uncertainty.derivation_interrupted
    try:
        _materialize_evaluation(experiments, result)
    except (Exception, KeyboardInterrupt) as error:
        if uncertainty_interrupted:
            _propagate_uncertainty_interruption(error)
        raise
    variable_count = len(problem.controlled_ids)
    if isinstance(search, GridSearch):
        try:
            if not isinstance(outcome, ProfiledGridOutcome):
                raise TypeError("GRID execution returned a non-GRID outcome")
            write_grid_output(
                outcome,
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
                deterministic_uncertainty=deterministic_uncertainty,
                parameter_model=parameter_model,
                parameter_values=result.resolved_values,
                parameterization=parameterization,
                fitted_ids=problem.controlled_ids,
            )
        except (Exception, KeyboardInterrupt) as error:  # noqa: BLE001
            if uncertainty_interrupted:
                if isinstance(error, OSError):
                    _propagate_output_failure(error, path, interrupted=True)
                _propagate_uncertainty_interruption(error)
            _propagate_output_failure(error, path)
        if uncertainty_interrupted:
            _propagate_uncertainty_interruption()
        return fit

    try:
        execute_post_fit(
            experiments,
            path,
            plot=plot != "nothing",
            residuals=result.residuals,
            nvarys=variable_count,
            deterministic_uncertainty=deterministic_uncertainty,
            parameter_model=parameter_model,
            parameter_values=result.resolved_values,
            parameterization=parameterization,
            fitted_ids=problem.controlled_ids,
        )
    except (Exception, KeyboardInterrupt) as error:  # noqa: BLE001
        if uncertainty_interrupted:
            if isinstance(error, OSError):
                _propagate_output_failure(error, path, interrupted=True)
            _propagate_uncertainty_interruption(error)
        _propagate_output_failure(error, path)
    if uncertainty_interrupted:
        _propagate_uncertainty_interruption()
    return fit
