"""Behavioral qualification tests for the bounded native Direct-TRF slice (#586).

The public seams are the canonical optimization problem/invocation, the typed
terminal outcome, accepted-result materialization, and revision-checked commit.
The production runner uses the native deterministic implementation qualified here.
"""

from __future__ import annotations

import copy
import dataclasses
import math
import os
import pickle
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import SimpleNamespace
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.method_execution import normalize_methods_for_execution
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    AffineEquality,
    AffineHalfSpace,
    CancellationToken,
    CandidateSummary,
    DirectTrfCandidateTerminal,
    DirectTrfConstructionError,
    DirectTrfInterrupted,
    DirectTrfInvocation,
    DirectTrfOutcomeTerminal,
    DirectTrfScalePolicy,
    DirectTrfTerminal,
    LiveFitCommitAuthority,
    MaterializationTerminal,
    ObjectiveScalarizationError,
    OptimizationProblem,
    ResidualJacobianSource,
    RootMaterializationFailure,
    canonical_chi_square,
    commit_accepted_fit,
    execute_direct_trf,
    execute_direct_trf_candidate,
    materialize_root_candidate,
)
from chemex.optimize.grouped_direct_trf import FitDecomposition
from chemex.optimize.progress import ProgressEvent, ProgressPhase
from chemex.parameters.parameterization import (
    ActiveParameterization,
    compile_active_parameterization,
)
from chemex.parameters.sealed import ParamConfig, SealedConfiguration
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot, StaleAnalysisValuesError
from chemex.runtime import AnalysisSession, ExecutionSettings
from chemex.runtime.execution import NATIVE_THREAD_ENV_VARS
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
CPMG_ROOT = ROOT / "examples/Experiments/CPMG_15N_IP"


def _shipped_method(path: Path, section: str) -> Method:
    plan = read_method_plan([path])
    _plan, operational = normalize_methods_for_execution(plan)
    return operational[section]


def _qualification_fit(
    *, budget: int = 80
) -> tuple[
    AnalysisSession,
    Experiments,
    ActiveParameterization,
    EvaluationEngine,
    OptimizationProblem,
    DirectTrfInvocation,
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(
            include=[SpinSystem.from_name("G2N-HN")],
            exclude=None,
        ),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    plan = read_method_plan([METHOD])
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()["DEFAULT"], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=budget,
    )
    return session, experiments, parameterization, engine, problem, invocation


def _full_frame_coefficients(
    problem: OptimizationProblem,
    param_id: str,
    coefficient: float = 1.0,
) -> tuple[float, ...]:
    return tuple(
        coefficient if independent_id == param_id else 0.0
        for independent_id, _value in problem.independent_items
    )


def test_lifecycle_frame_rejects_exact_box_and_affine_infeasibility() -> None:
    _session, _experiments, parameterization, _engine, problem, _invocation = (
        _qualification_fit()
    )
    param_id = problem.controlled_ids[0]
    active_value = problem.start[0]
    constrained = dataclasses.replace(
        problem,
        affine_half_spaces=(
            AffineHalfSpace(
                "active-upper",
                _full_frame_coefficients(problem, param_id),
                active_value,
            ),
        ),
    )
    constrained.lifecycle_frame(constrained.start, parameterization)
    with pytest.raises(DirectTrfConstructionError, match="affine half-space"):
        constrained.lifecycle_frame(
            (float(np.nextafter(active_value, math.inf)),),
            parameterization,
        )
    equality = dataclasses.replace(
        problem,
        affine_equalities=(
            AffineEquality(
                "fixed-coordinate",
                _full_frame_coefficients(problem, param_id),
                active_value,
            ),
        ),
    )
    with pytest.raises(DirectTrfConstructionError, match="affine equality"):
        equality.lifecycle_frame(
            (float(np.nextafter(active_value, -math.inf)),),
            parameterization,
        )
    with pytest.raises(DirectTrfConstructionError, match="box bounds"):
        problem.lifecycle_frame((problem.lower_bounds[0] - 1.0,), parameterization)


def test_problem_construction_rejects_out_of_bounds_held_independent_value() -> None:
    session, _experiments, parameterization, engine, problem, _invocation = (
        _qualification_fit()
    )
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    held_id = next(
        param_id
        for param_id, _value in problem.held_items
        if math.isfinite(configuration[param_id].lower_bound)
    )
    invalid_value = configuration[held_id].lower_bound - 1.0
    invalid_snapshot = dataclasses.replace(
        problem.source_snapshot,
        _items=tuple(
            (param_id, invalid_value if param_id == held_id else value)
            for param_id, value in problem.source_snapshot.items()
        ),
    )

    with pytest.raises(DirectTrfConstructionError, match="Independent parameter"):
        OptimizationProblem.from_native(
            engine.plan,
            parameterization,
            configuration,
            invalid_snapshot,
        )


def test_root_owned_child_derivation_preserves_context_and_rejects_forgery() -> None:
    _session, _experiments, _parameterization, _engine, problem, _invocation = (
        _qualification_fit()
    )
    param_id = problem.controlled_ids[0]
    constrained = dataclasses.replace(
        problem,
        affine_half_spaces=(
            AffineHalfSpace(
                "derived-root-upper",
                _full_frame_coefficients(problem, param_id),
                problem.start[0] + 1.0,
            ),
        ),
        affine_equalities=(
            AffineEquality(
                "derived-root-fixed",
                _full_frame_coefficients(problem, param_id),
                problem.start[0],
            ),
        ),
    )
    child = constrained.derive_profiled_grid_point(
        factor_identity="qualified-profiled-grid-factor",
        point_ordinal=0,
        projected_plan_identity=constrained.evaluation_plan_identity,
        grid_items=(),
        controlled_ids=constrained.controlled_ids,
    )

    assert child.source_snapshot is constrained.source_snapshot
    assert child.evaluation_plan_identity == constrained.evaluation_plan_identity
    assert child.parameterization_identity == constrained.parameterization_identity
    assert (
        child.evaluator_parameterization_identity
        == constrained.evaluator_parameterization_identity
    )
    assert child.constraint_program_identity == constrained.constraint_program_identity
    assert child.configuration_identity == constrained.configuration_identity
    assert child.independent_items == constrained.independent_items
    assert child.held_items == constrained.held_items
    assert child.lower_bounds == constrained.lower_bounds
    assert child.upper_bounds == constrained.upper_bounds
    assert child.commit_scope == constrained.commit_scope
    assert child.scalarization_version == constrained.scalarization_version
    assert child.affine_half_spaces == constrained.affine_half_spaces
    assert child.affine_equalities == constrained.affine_equalities
    constrained.validate_derived_problem(child)

    forged_fields = (
        {"evaluation_plan_identity": "foreign-plan"},
        {"parameterization_identity": "foreign-parameterization"},
        {"evaluator_parameterization_identity": "foreign-evaluator"},
        {"constraint_program_identity": "foreign-constraints"},
        {"configuration_identity": "foreign-configuration"},
        {
            "source_snapshot": dataclasses.replace(
                child.source_snapshot,
                revision=child.source_snapshot.revision + 1,
            ),
        },
        {
            "independent_items": (
                (child.independent_items[0][0], child.independent_items[0][1] + 1.0),
                *child.independent_items[1:],
            ),
        },
        {"lower_bounds": (child.lower_bounds[0] - 1.0,)},
        {"upper_bounds": (child.start[0] + 2.0,)},
        {"commit_scope": tuple(reversed(child.commit_scope))},
        {"scalarization_version": "foreign-scalarization"},
        {"affine_half_spaces": ()},
        {"affine_equalities": ()},
    )
    for replacement in forged_fields:
        with pytest.raises(DirectTrfConstructionError):
            forged = dataclasses.replace(child, **replacement)
            constrained.validate_derived_problem(forged)


def test_pre_cancelled_shared_root_materialization_never_binds_or_promotes() -> None:
    session, _experiments, parameterization, engine, problem, _invocation = (
        _qualification_fit()
    )
    before_state = session.analysis_values.snapshot()
    cancelled_before = CancellationToken()
    cancelled_before.cancel()
    binding_calls = 0
    evaluation_calls = 0

    def reject_binding(_engine: EvaluationEngine) -> BoundEvaluator:
        nonlocal binding_calls
        binding_calls += 1
        raise AssertionError("pre-cancelled materialization must not bind an evaluator")

    def reject_evaluation(
        _evaluator: BoundEvaluator,
        _frame: EvaluationFrame,
    ) -> EvaluationResult:
        nonlocal evaluation_calls
        evaluation_calls += 1
        raise AssertionError("pre-cancelled materialization must not evaluate")

    with (
        patch.object(EvaluationEngine, "new_evaluator", reject_binding),
        patch.object(BoundEvaluator, "evaluate", reject_evaluation),
    ):
        outcome = materialize_root_candidate(
            problem,
            parameterization,
            engine,
            vector=problem.start,
            invocation_identity="root-materialization-invocation",
            execution_identity=lambda _evaluated: "root-materialization-execution",
            cancellation=cancelled_before,
        )

    assert binding_calls == 0
    assert evaluation_calls == 0
    assert isinstance(outcome, RootMaterializationFailure)
    assert outcome.terminal is MaterializationTerminal.CANCELLED
    assert outcome.evaluation_count == 0
    # This exact terminal type carries no materialized candidate, accepted result,
    # or commit authority, so promotion cannot have begun.
    assert session.analysis_values.snapshot() == before_state


def test_shared_root_materialization_preserves_later_cancellation_and_interruption_gates() -> (
    None
):
    _session, _experiments, parameterization, engine, problem, _invocation = (
        _qualification_fit()
    )

    def interrupt_binding(_engine: EvaluationEngine) -> BoundEvaluator:
        raise KeyboardInterrupt

    with patch.object(EvaluationEngine, "new_evaluator", interrupt_binding):
        interrupted_binding = materialize_root_candidate(
            problem,
            parameterization,
            engine,
            vector=problem.start,
            invocation_identity="root-materialization-invocation",
            execution_identity=lambda _evaluated: "root-materialization-execution",
        )

    assert isinstance(interrupted_binding, RootMaterializationFailure)
    assert interrupted_binding.terminal is MaterializationTerminal.INTERRUPTED
    assert interrupted_binding.evaluation_count == 0

    cancelled_during = CancellationToken()
    original_evaluate = BoundEvaluator.evaluate

    def evaluate_then_cancel(
        evaluator: BoundEvaluator,
        frame: EvaluationFrame,
    ) -> EvaluationResult:
        evaluated = original_evaluate(evaluator, frame)
        assert isinstance(evaluated, EvaluationResult)
        cancelled_during.cancel()
        return evaluated

    with patch.object(BoundEvaluator, "evaluate", evaluate_then_cancel):
        during = materialize_root_candidate(
            problem,
            parameterization,
            engine,
            vector=problem.start,
            invocation_identity="root-materialization-invocation",
            execution_identity=lambda _evaluated: "root-materialization-execution",
            cancellation=cancelled_during,
        )

    assert isinstance(during, RootMaterializationFailure)
    assert during.terminal is MaterializationTerminal.CANCELLED
    assert during.evaluation_count == 1

    def interrupt_validation(_evaluated: EvaluationResult) -> str:
        raise KeyboardInterrupt

    interrupted = materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=problem.start,
        invocation_identity="root-materialization-invocation",
        execution_identity=interrupt_validation,
    )

    assert isinstance(interrupted, RootMaterializationFailure)
    assert interrupted.terminal is MaterializationTerminal.INTERRUPTED
    assert interrupted.evaluation_count == 1


def _presentation_values(
    session: AnalysisSession, scope: tuple[str, ...]
) -> dict[str, float | None]:
    return {
        param_id: parameter.value
        for param_id, parameter in session.parameters.get_parameters(set(scope)).items()
    }


def _backend_result(
    candidate: Array,
    residuals: Array,
    *,
    status: int,
    success: bool,
    message: str,
    optimality: float,
) -> object:
    return SimpleNamespace(
        x=candidate,
        fun=residuals,
        status=status,
        success=success,
        message=message,
        nfev=1,
        njev=0,
        cost=0.5 * float(np.dot(residuals, residuals)),
        optimality=optimality,
        active_mask=np.zeros_like(candidate, dtype=np.int64),
        jac=np.zeros((residuals.size, candidate.size), dtype=np.float64),
    )


def test_accepted_fit_retains_exact_final_backend_residual_jacobian() -> None:
    (
        _session,
        _experiments,
        parameterization,
        engine,
        problem,
        invocation,
    ) = _qualification_fit()

    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
    )

    accepted = outcome.accepted_result
    assert accepted is not None
    retained = accepted.final_residual_jacobian
    assert retained is not None
    assert retained.source is ResidualJacobianSource.SCIPY_FINAL_2_POINT
    assert retained.controlled_ids == problem.controlled_ids
    assert retained.final_vector == accepted.vector
    assert retained.final_residuals == tuple(accepted.evaluation_result.residuals)
    assert retained.shape == (
        accepted.evaluation_result.residuals.size,
        len(problem.controlled_ids),
    )
    assert retained.derivative_method == "numerical 2-point"
    assert retained.diff_step_policy == "scipy-default-relative-step"
    assert retained.loss_policy == "linear"
    assert retained.external_coordinate_policy == "physical-external-unscaled"
    assert retained.trust_region_scale_policy == (
        "scipy-adaptive-inverse-jacobian-column-norm-v1"
    )
    assert len(retained.matrix_binary64_sha256) == 64
    assert np.isfinite(np.asarray(retained.matrix)).all()
    assert outcome.execution.backend is not None
    assert outcome.execution.backend.final_residual_jacobian is retained


@pytest.mark.parametrize("malformation", ("wrong-shape", "non-finite"))
def test_malformed_final_backend_jacobian_fails_closed(malformation: str) -> None:
    (
        _session,
        _experiments,
        parameterization,
        engine,
        problem,
        invocation,
    ) = _qualification_fit()

    def malformed_backend(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        result = _backend_result(
            candidate,
            residuals,
            status=1,
            success=True,
            message="converged",
            optimality=0.0,
        )
        result.jac = (
            np.zeros((residuals.size, candidate.size + 1), dtype=np.float64)
            if malformation == "wrong-shape"
            else np.full((residuals.size, candidate.size), np.nan)
        )
        return result

    with patch("chemex.optimize.direct_trf.least_squares", malformed_backend):
        outcome = execute_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is DirectTrfOutcomeTerminal.SOLVER_UNSUCCESSFUL
    assert outcome.execution.terminal is DirectTrfTerminal.IMPLEMENTATION_FAILURE
    assert outcome.execution.failure is not None
    assert outcome.execution.failure.category == "malformed_backend_result"
    assert outcome.accepted_result is None


def _one_request_converged_backend(
    fun: Callable[[Array], Array],
    x0: Array,
    *,
    after_request: Callable[[], None] | None = None,
) -> object:
    candidate = np.asarray(x0, dtype=np.float64) * 0.9
    residuals = fun(candidate)
    if after_request is not None:
        after_request()
    return _backend_result(
        candidate,
        residuals,
        status=1,
        success=True,
        message="gradient tolerance satisfied",
        optimality=0.0,
    )


def test_representative_single_component_fit_materializes_and_commits_atomically() -> (
    None
):
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)

    first = execute_direct_trf(problem, invocation, parameterization, engine)
    second = execute_direct_trf(problem, invocation, parameterization, engine)

    assert problem.controlled_ids == ("__R1A_A_G2N_H_800_0MHZ",)
    assert problem.lower_bounds == (0.0,)
    assert problem.upper_bounds == (float("inf"),)
    assert problem.start == (6.87922079444668,)
    assert first.terminal is DirectTrfOutcomeTerminal.ACCEPTED
    assert first.execution.terminal is DirectTrfTerminal.CONVERGED
    assert first.accepted_result is not None
    assert first.commit_authority is not None
    assert first.materialization is not None
    assert first.materialization.terminal is MaterializationTerminal.SUCCESS
    assert first.materialization.evaluation_count == 1
    assert first.materialization.cache_hits == 0
    # The literals come from the independent legacy least_squares baseline at
    # c6378fb0. A 2e-8 relative tolerance covers the observed native ordered-
    # normalization plus finite-difference/TRF rounding while remaining far
    # below any scientifically meaningful change in this relaxation rate.
    assert first.accepted_result.vector == (
        pytest.approx(2.3474211504, rel=2.0e-8, abs=1.0e-10),
    )
    assert first.accepted_result.chi_square == pytest.approx(
        13.2171307054,
        rel=1.0e-9,
        abs=1.0e-10,
    )
    assert first.execution.counters.solver_requests_received <= 80
    assert first.execution.counters.objective_requests_accepted == (
        first.execution.counters.objective_evaluations_completed
    )
    assert first.execution.identity == second.execution.identity
    assert second.materialization is not None
    assert second.accepted_result is not None
    assert second.commit_authority is not None
    assert first.materialization.identity == second.materialization.identity
    assert first.accepted_result.identity == second.accepted_result.identity
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before

    receipt = commit_accepted_fit(
        first.accepted_result,
        first.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )

    committed = session.analysis_values.snapshot()
    assert receipt.old_revision == 0
    assert receipt.new_revision == 1
    assert receipt.scope == problem.commit_scope
    assert committed.revision == 1
    assert committed[problem.controlled_ids[0]] == pytest.approx(
        2.3474211504,
        rel=2.0e-8,
        abs=1.0e-10,
    )
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            first.accepted_result,
            first.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == committed

    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    plan = read_method_plan([METHOD])
    revision_one_parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()["DEFAULT"],
        experiments.param_ids,
    )
    revision_one_engine = EvaluationEngine.from_experiments(
        experiments,
        revision_one_parameterization,
    )
    revision_one_problem = OptimizationProblem.from_native(
        revision_one_engine.plan,
        revision_one_parameterization,
        configuration,
        committed,
    )
    revision_one_invocation = DirectTrfInvocation.for_problem(
        revision_one_problem,
        objective_request_budget=invocation.objective_request_budget,
        scale_policy=invocation.scale_policy,
        ftol=invocation.ftol,
        xtol=invocation.xtol,
        gtol=invocation.gtol,
    )
    revision_one = execute_direct_trf(
        revision_one_problem,
        revision_one_invocation,
        revision_one_parameterization,
        revision_one_engine,
    )
    assert revision_one.accepted_result is not None
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            revision_one.accepted_result,
            first.commit_authority,
            problem=revision_one_problem,
            parameterization=revision_one_parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == committed
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_solver_requests_use_lean_residuals_and_acceptance_materializes_fresh() -> None:
    _session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    complete_calls = 0
    residual_calls = 0
    original_evaluate = BoundEvaluator.evaluate
    original_evaluate_residuals = BoundEvaluator.evaluate_residuals

    def counted_complete(
        evaluator: BoundEvaluator,
        frame: EvaluationFrame,
    ) -> EvaluationResult | EvaluationFailure:
        nonlocal complete_calls
        complete_calls += 1
        return original_evaluate(evaluator, frame)

    def counted_residuals(
        evaluator: BoundEvaluator,
        frame: EvaluationFrame,
    ) -> Array | EvaluationFailure:
        nonlocal residual_calls
        residual_calls += 1
        return original_evaluate_residuals(evaluator, frame)

    with (
        patch.object(BoundEvaluator, "evaluate", counted_complete),
        patch.object(BoundEvaluator, "evaluate_residuals", counted_residuals),
    ):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.terminal is DirectTrfOutcomeTerminal.ACCEPTED
    assert outcome.materialization is not None
    assert outcome.materialization.evaluation_count == 1
    assert outcome.materialization.cache_hits == 0
    assert outcome.execution.counters.solver_requests_received == 16
    assert complete_calls == 2  # one preflight and one fresh accepted materialization
    assert residual_calls == outcome.execution.counters.objective_evaluations_completed


def test_cpmg_step1_direct_trf_preserves_requests_and_reuses_profile_kernels() -> None:
    method_path = CPMG_ROOT / "Methods/method.toml"
    plan = read_method_plan([method_path])
    method = _shipped_method(method_path, "STEP1")
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        sorted((CPMG_ROOT / "Experiments").glob("*.toml")),
        method.selection,
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([CPMG_ROOT / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()["STEP1"], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    root_problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    decomposition = FitDecomposition.from_root(
        root_problem,
        parameterization,
        engine,
    )
    assert len(decomposition.components) == 1
    component = decomposition.components[0]
    component_engine = engine.project_profiles(component.root_profile_indices)
    invocation = DirectTrfInvocation.for_problem(
        component.problem,
        objective_request_budget=2000 * (len(component.controlled_ids) + 1),
    )
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original_calculate = pulse_type.calculate
    kernel_calls = 0

    def counted_kernel(self: object, spectrometer: object, data: object) -> Array:
        nonlocal kernel_calls
        kernel_calls += 1
        return original_calculate(self, spectrometer, data)

    with patch.object(pulse_type, "calculate", counted_kernel):
        outcome = execute_direct_trf_candidate(
            component.problem,
            invocation,
            parameterization,
            component_engine,
        )

    assert outcome.terminal is DirectTrfCandidateTerminal.SUCCESS
    assert outcome.candidate is not None
    assert outcome.execution.backend is not None
    assert outcome.execution.backend.nfev == 7
    assert outcome.execution.backend.njev == 7
    assert outcome.execution.counters.solver_requests_received == 126
    assert outcome.execution.counters.objective_evaluations_completed == 126
    assert kernel_calls == 360
    assert outcome.candidate.chi_square == pytest.approx(285.8191490381348)


def test_direct_trf_reports_objective_evaluations_without_iteration_terminology() -> (
    None
):
    _session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    events: list[ProgressEvent] = []

    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
        progress_observer=events.append,
    )

    assert outcome.terminal is DirectTrfOutcomeTerminal.ACCEPTED
    assert events[0].phase is ProgressPhase.STARTED
    assert events[-1].phase is ProgressPhase.TERMINATED
    evaluated = [event for event in events if event.phase is ProgressPhase.EVALUATED]
    assert len(evaluated) == outcome.execution.counters.objective_evaluations_completed
    assert [event.objective_evaluations_completed for event in evaluated] == list(
        range(1, len(evaluated) + 1)
    )
    assert all(event.current_chi_square is not None for event in evaluated)
    assert all(event.best_chi_square is not None for event in evaluated)
    assert events[-1].terminal_status == "accepted"
    assert not any("iteration" in name for name in ProgressEvent.__dataclass_fields__)


def test_progress_observation_has_no_numerical_or_kernel_effect() -> None:
    disabled = _qualification_fit()
    enabled = _qualification_fit()
    (
        disabled_parameterization,
        disabled_engine,
        disabled_problem,
        disabled_invocation,
    ) = disabled[2:]
    enabled_parameterization, enabled_engine, enabled_problem, enabled_invocation = (
        enabled[2:]
    )
    pulse_type = type(next(iter(disabled[1])).profiles[0].pulse_sequence)
    original_calculate = cast("Callable[..., Array]", pulse_type.calculate)
    kernel_calls = [0, 0]
    active_run = 0

    def counted_kernel(self: object, spectrometer: object, data: object) -> Array:
        kernel_calls[active_run] += 1
        return original_calculate(self, spectrometer, data)

    events: list[ProgressEvent] = []

    with patch.object(pulse_type, "calculate", counted_kernel):
        without_progress = execute_direct_trf(
            disabled_problem,
            disabled_invocation,
            disabled_parameterization,
            disabled_engine,
        )
        active_run = 1
        with_progress = execute_direct_trf(
            enabled_problem,
            enabled_invocation,
            enabled_parameterization,
            enabled_engine,
            progress_observer=events.append,
        )

    assert len(events) > 2
    assert kernel_calls[0] == kernel_calls[1]
    assert without_progress.terminal is with_progress.terminal
    assert without_progress.execution.counters == with_progress.execution.counters
    assert without_progress.execution.backend == with_progress.execution.backend
    assert without_progress.accepted_result is not None
    assert with_progress.accepted_result is not None
    first = without_progress.accepted_result
    second = with_progress.accepted_result
    assert first.vector == second.vector
    assert first.chi_square == second.chi_square
    assert first.evaluation_result.residuals.tobytes() == (
        second.evaluation_result.residuals.tobytes()
    )
    assert without_progress.materialization is not None
    assert with_progress.materialization is not None
    assert without_progress.materialization.cache_hits == (
        with_progress.materialization.cache_hits
    )
    assert without_progress.materialization.cache_misses == (
        with_progress.materialization.cache_misses
    )


def test_progress_observer_failure_is_suppressed_after_first_event() -> None:
    _session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    observer_calls = 0

    def failing_observer(_event: ProgressEvent) -> None:
        nonlocal observer_calls
        observer_calls += 1
        raise RuntimeError("non-scientific progress failure")

    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
        progress_observer=failing_observer,
    )

    assert outcome.terminal is DirectTrfOutcomeTerminal.ACCEPTED
    assert observer_calls == 1


def test_live_commit_authority_is_atomic_under_concurrent_use() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    assert invocation.scale_policy is (
        DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM
    )
    outcome = execute_direct_trf(problem, invocation, parameterization, engine)
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    barrier = Barrier(2)

    def commit_once() -> object:
        barrier.wait()
        try:
            return commit_accepted_fit(
                outcome.accepted_result,
                outcome.commit_authority,
                problem=problem,
                parameterization=parameterization,
                analysis_values=session.analysis_values,
            )
        except DirectTrfConstructionError as error:
            return error

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = tuple(executor.map(lambda _index: commit_once(), range(2)))

    errors = tuple(
        result for result in results if isinstance(result, DirectTrfConstructionError)
    )
    receipts = tuple(result for result in results if result not in errors)
    assert len(errors) == len(receipts) == 1
    assert "exact live fit commit authority" in str(errors[0])
    assert session.analysis_values.snapshot().revision == before.revision + 1
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_commit_rejects_absent_foreign_or_wrongly_bound_live_authority() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    first = execute_direct_trf(problem, invocation, parameterization, engine)
    second = execute_direct_trf(problem, invocation, parameterization, engine)
    assert first.accepted_result is not None
    assert first.commit_authority is not None
    assert second.accepted_result is not None
    assert second.commit_authority is not None
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)

    with pytest.raises(TypeError, match="minted only"):
        LiveFitCommitAuthority()
    with pytest.raises(TypeError, match="cannot be copied"):
        copy.copy(first.commit_authority)
    with pytest.raises(TypeError, match="cannot be copied"):
        copy.deepcopy(first.commit_authority)
    with pytest.raises(TypeError, match="cannot be serialized"):
        pickle.dumps(first.commit_authority)
    with pytest.raises(AttributeError):
        first.commit_authority.origin_context_identity = "foreign-origin"  # type: ignore[attr-defined]
    with pytest.raises(AttributeError):
        first.commit_authority._problem_identity = "foreign-problem"  # type: ignore[attr-defined]
    equivalent_evidence = dataclasses.replace(first.accepted_result)
    assert equivalent_evidence is not first.accepted_result
    assert equivalent_evidence.identity == first.accepted_result.identity
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            equivalent_evidence,
            first.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            first.accepted_result,
            second.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            first.accepted_result,
            None,  # type: ignore[arg-type] - adversarial runtime absence
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    wrong_snapshot = dataclasses.replace(
        problem.source_snapshot,
        occurrence_identity="foreign-analysis-values",
        revision=problem.source_snapshot.revision + 1,
    )
    wrong_problem = dataclasses.replace(problem, source_snapshot=wrong_snapshot)
    with pytest.raises(
        DirectTrfConstructionError,
        match="commit context",
    ):
        commit_accepted_fit(
            first.accepted_result,
            first.commit_authority,
            problem=wrong_problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_canonical_scalarization_and_candidate_order_are_explicit() -> None:
    assert canonical_chi_square((1.0, 2.0, 3.0)) == 14.0
    with pytest.raises(ObjectiveScalarizationError):
        canonical_chi_square((1.0, float("inf")))

    candidates = (
        CandidateSummary((2.0,), 4.0, 0),
        CandidateSummary((1.0,), 4.0, 2),
        CandidateSummary((1.0,), 4.0, 1),
        CandidateSummary((9.0,), 3.0, 3),
    )
    assert min(candidates, key=CandidateSummary.ordering_key) == candidates[3]
    assert min(candidates[:3], key=CandidateSummary.ordering_key) == candidates[2]


def test_multicoordinate_problem_preserves_canonical_vector_and_bound_alignment() -> (
    None
):
    session, experiments, _parameterization, _engine, _problem, _invocation = (
        _qualification_fit()
    )
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    source = session.analysis_values.snapshot()
    r1a_id = "__R1A_A_G2N_H_800_0MHZ"
    qualification_values = {
        "__PB": (0.17, 0.05, 0.40),
        r1a_id: (7.25, 2.0, 15.0),
    }

    def qualify_config(item: ParamConfig) -> ParamConfig:
        value, lower, upper = qualification_values.get(
            item.param_id,
            (item.effective_value, item.lower_bound, item.upper_bound),
        )
        return ParamConfig(item.param_id, value, lower, upper)

    configs = tuple(qualify_config(item) for item in parameter_model.configuration)
    configuration = SealedConfiguration(
        configs,
        {},
        definitions_identity=parameter_model.definitions.identity,
    )
    qualified_model = dataclasses.replace(
        parameter_model,
        configuration=configuration,
    )
    snapshot = AnalysisValuesSnapshot(
        source.occurrence_identity,
        source.model_identity,
        source.definitions_identity,
        configuration.identity,
        source.revision,
        tuple(
            (
                param_id,
                qualification_values.get(param_id, (value, 0.0, 0.0))[0],
            )
            for param_id, value in source.items()
        ),
    )
    method = Method(fit=["PB"], fix=["KEX_AB"])

    def reconstruct() -> tuple[ActiveParameterization, OptimizationProblem]:
        parameterization = compile_active_parameterization(
            qualified_model,
            snapshot,
            method,
            experiments.param_ids,
        )
        engine = EvaluationEngine.from_experiments(experiments, parameterization)
        problem = OptimizationProblem.from_native(
            engine.plan,
            parameterization,
            configuration,
            snapshot,
        )
        return parameterization, problem

    first_parameterization, first = reconstruct()
    second_parameterization, second = reconstruct()

    expected_ids = ("__PB", r1a_id)
    assert first.controlled_ids == expected_ids
    assert (
        tuple(
            param_id
            for param_id, _value in first.independent_items
            if param_id in expected_ids
        )
        == expected_ids
    )
    assert first.start == (0.17, 7.25)
    assert first.lower_bounds == (0.05, 2.0)
    assert first.upper_bounds == (0.40, 15.0)
    assert (
        first.lifecycle_frame(
            first.start,
            first_parameterization,
        ).ordered_items()
        == first.independent_items
    )
    assert first_parameterization.identity == second_parameterization.identity
    assert first.independent_items == second.independent_items
    assert first.start == second.start
    assert first.lower_bounds == second.lower_bounds
    assert first.upper_bounds == second.upper_bounds
    assert first.identity == second.identity


@pytest.mark.parametrize(
    "budget",
    (True, False, 2.0, 2.5, "2", None, 0, -1),
)
def test_objective_request_budget_requires_a_positive_integer(budget: object) -> None:
    with pytest.raises(
        DirectTrfConstructionError,
        match="Objective-request budget must be a positive integer",
    ):
        DirectTrfInvocation(
            "qualification-problem",
            cast("int", budget),
        )


def test_invocation_rejects_unrecognized_scale_policy() -> None:
    with pytest.raises(
        DirectTrfConstructionError,
        match="Unsupported Direct TRF scale policy",
    ):
        DirectTrfInvocation(
            "qualification-problem",
            10,
            cast("DirectTrfScalePolicy", "start-magnitude"),
        )


def test_non_convergence_keeps_last_iterate_diagnostic_and_commits_nothing() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)

    def non_converged(
        fun: Callable[[Array], Array], x0: Array, **settings: object
    ) -> object:
        bounds = settings.pop("bounds")
        assert isinstance(bounds, tuple)
        np.testing.assert_array_equal(bounds[0], problem.lower_bounds)
        np.testing.assert_array_equal(bounds[1], problem.upper_bounds)
        assert settings.pop("x_scale") == "jac"
        assert settings == {
            "method": "trf",
            "jac": "2-point",
            "diff_step": None,
            "tr_solver": "exact",
            "tr_options": {},
            "loss": "linear",
            "f_scale": 1.0,
            "ftol": 1.0e-8,
            "xtol": 1.0e-8,
            "gtol": 1.0e-8,
            "max_nfev": 80,
            "verbose": 0,
        }
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        return _backend_result(
            candidate,
            residuals,
            status=0,
            success=False,
            message="maximum evaluations reached",
            optimality=1.0,
        )

    with patch("chemex.optimize.direct_trf.least_squares", non_converged):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.terminal is DirectTrfOutcomeTerminal.SOLVER_UNSUCCESSFUL
    assert outcome.execution.terminal is DirectTrfTerminal.NON_CONVERGED
    assert outcome.execution.best_candidate is not None
    assert outcome.materialization is None
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_invocation_execution_settings_drive_the_solver_environment() -> None:
    _session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    invocation = dataclasses.replace(
        invocation,
        execution_settings=ExecutionSettings(workers=1, native_threads=7),
    )
    previous = {name: os.environ.get(name) for name in NATIVE_THREAD_ENV_VARS}

    def inspect_execution_settings(
        fun: Callable[[Array], Array], x0: Array, **settings: object
    ) -> object:
        assert "workers" not in settings
        assert all(os.environ.get(name) == "7" for name in NATIVE_THREAD_ENV_VARS)
        return _one_request_converged_backend(fun, x0)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        inspect_execution_settings,
    ):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.terminal is DirectTrfOutcomeTerminal.ACCEPTED
    assert {name: os.environ.get(name) for name in NATIVE_THREAD_ENV_VARS} == previous


def test_direct_trf_rejects_a_parallel_worker_claim() -> None:
    _session, _experiments, _parameterization, _engine, _problem, invocation = (
        _qualification_fit()
    )

    with pytest.raises(
        DirectTrfConstructionError,
        match="exactly one ChemEx worker",
    ):
        dataclasses.replace(
            invocation,
            execution_settings=ExecutionSettings(workers=2),
        )


def test_budget_exhaustion_has_exact_counters_and_no_accepted_candidate() -> None:
    session, _experiments, parameterization, engine, problem, _invocation = (
        _qualification_fit()
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=2,
    )
    before = session.analysis_values.snapshot()

    def exhaust(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        for _ in range(3):
            fun(np.array(x0, copy=True))
        raise AssertionError("the third request must stop the backend")

    with patch("chemex.optimize.direct_trf.least_squares", exhaust):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.execution.terminal is DirectTrfTerminal.BUDGET_EXHAUSTED
    assert outcome.execution.counters.solver_requests_received == 3
    assert outcome.execution.counters.objective_requests_accepted == 2
    assert outcome.execution.counters.objective_evaluations_completed == 2
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_cancellation_and_interruption_never_materialize_or_commit() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    token = CancellationToken()
    token.cancel()

    cancelled = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
        cancellation=token,
    )

    assert cancelled.terminal is DirectTrfOutcomeTerminal.CANCELLED
    assert cancelled.execution.terminal is DirectTrfTerminal.CANCELLED
    assert cancelled.execution.counters.solver_requests_received == 0
    assert cancelled.materialization is None
    assert cancelled.accepted_result is None
    assert session.analysis_values.snapshot() == before

    def interrupt_after_iterate(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        fun(np.asarray(x0, dtype=np.float64) * 0.9)
        raise KeyboardInterrupt

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            interrupt_after_iterate,
        ),
        pytest.raises(DirectTrfInterrupted) as interrupted,
    ):
        execute_direct_trf(problem, invocation, parameterization, engine)
    assert interrupted.value.execution.terminal is DirectTrfTerminal.INTERRUPTED
    assert interrupted.value.execution.best_candidate is not None
    assert session.analysis_values.snapshot() == before


def test_exact_start_preflight_interruption_freezes_typed_attempt_evidence() -> None:
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)

    with (
        patch.object(pulse_type, "calculate", side_effect=KeyboardInterrupt),
        pytest.raises(DirectTrfInterrupted) as interrupted,
    ):
        execute_direct_trf(problem, invocation, parameterization, engine)

    execution = interrupted.value.execution
    assert execution.terminal is DirectTrfTerminal.INTERRUPTED
    assert execution.counters.solver_requests_received == 0
    assert execution.counters.objective_requests_accepted == 0
    assert execution.counters.objective_evaluations_completed == 0
    assert execution.preflight_evaluation_identity is None
    assert execution.best_candidate is None
    assert execution.final_candidate is None
    assert execution.backend is None
    assert execution.failure is not None
    assert execution.failure.category == "interrupted"
    assert interrupted.value.materialization is None
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_cancellation_after_exact_start_preflight_never_enters_solver() -> None:
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    token = CancellationToken()
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)

    def calculate_then_cancel(*args: object, **kwargs: object) -> Array:
        result = original(*args, **kwargs)
        token.cancel()
        return result

    with (
        patch.object(pulse_type, "calculate", calculate_then_cancel),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=AssertionError("solver must not be entered"),
        ) as backend,
    ):
        outcome = execute_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    backend.assert_not_called()
    assert outcome.terminal is DirectTrfOutcomeTerminal.CANCELLED
    assert outcome.execution.terminal is DirectTrfTerminal.CANCELLED
    assert outcome.execution.counters.solver_requests_received == 0
    assert outcome.execution.counters.objective_requests_accepted == 0
    assert outcome.execution.counters.objective_evaluations_completed == 0
    assert outcome.execution.preflight_evaluation_identity is not None
    assert outcome.execution.best_candidate is None
    assert outcome.execution.final_candidate is None
    assert outcome.execution.backend is None
    assert outcome.execution.failure is not None
    assert outcome.execution.failure.category == "cancelled"
    assert outcome.materialization is None
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_cancellation_after_a_valid_last_iterate_suppresses_materialization() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    token = CancellationToken()

    def cancel_after_iterate(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        return _one_request_converged_backend(
            fun,
            x0,
            after_request=token.cancel,
        )

    with patch("chemex.optimize.direct_trf.least_squares", cancel_after_iterate):
        outcome = execute_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert outcome.terminal is DirectTrfOutcomeTerminal.CANCELLED
    assert outcome.execution.terminal is DirectTrfTerminal.CANCELLED
    assert outcome.execution.best_candidate is not None
    assert outcome.materialization is None
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_native_problem_construction_failure_leaves_legacy_state_untouched() -> None:
    session, _experiments, parameterization, engine, _problem, _invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, parameterization.scope_ids)
    incomplete_plan = dataclasses.replace(engine.plan, resolved_ids=())
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None

    with pytest.raises(DirectTrfConstructionError):
        OptimizationProblem.from_native(
            incomplete_plan,
            parameterization,
            configuration,
            before,
        )

    assert session.analysis_values.snapshot() == before
    assert (
        _presentation_values(session, parameterization.scope_ids) == presentation_before
    )


def test_native_preflight_evaluation_failure_leaves_legacy_state_untouched() -> None:
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)

    with patch.object(
        pulse_type,
        "calculate",
        side_effect=RuntimeError("qualification kernel failure"),
    ):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.terminal is DirectTrfOutcomeTerminal.SOLVER_UNSUCCESSFUL
    assert outcome.execution.terminal is DirectTrfTerminal.PREFLIGHT_INVALID
    assert outcome.materialization is None
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_fresh_materialization_failure_cannot_accept_or_commit_last_iterate() -> None:
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)
    fail_materialization = False

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        if fail_materialization:
            raise RuntimeError("fresh materialization kernel failure")
        return original(*args, **kwargs)

    def arm_failure() -> None:
        nonlocal fail_materialization
        fail_materialization = True

    def converge_then_fail(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        return _one_request_converged_backend(
            fun,
            x0,
            after_request=arm_failure,
        )

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch("chemex.optimize.direct_trf.least_squares", converge_then_fail),
    ):
        outcome = execute_direct_trf(problem, invocation, parameterization, engine)

    assert outcome.terminal is DirectTrfOutcomeTerminal.MATERIALIZATION_FAILURE
    assert outcome.execution.terminal is DirectTrfTerminal.CONVERGED
    assert outcome.materialization is not None
    assert outcome.materialization.terminal is MaterializationTerminal.FAILURE
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_materialization_interruption_freezes_typed_evidence_and_commits_nothing() -> (
    None
):
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    presentation_before = _presentation_values(session, problem.commit_scope)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)
    interrupt_materialization = False

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        if interrupt_materialization:
            raise KeyboardInterrupt
        return original(*args, **kwargs)

    def arm_interruption() -> None:
        nonlocal interrupt_materialization
        interrupt_materialization = True

    def converge_then_interrupt(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        return _one_request_converged_backend(
            fun,
            x0,
            after_request=arm_interruption,
        )

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            converge_then_interrupt,
        ),
        pytest.raises(DirectTrfInterrupted) as interrupted,
    ):
        execute_direct_trf(problem, invocation, parameterization, engine)

    assert interrupted.value.execution.terminal is DirectTrfTerminal.CONVERGED
    assert interrupted.value.materialization is not None
    assert (
        interrupted.value.materialization.terminal
        is MaterializationTerminal.INTERRUPTED
    )
    assert interrupted.value.materialization.evaluation_count == 1
    assert session.analysis_values.snapshot() == before
    assert _presentation_values(session, problem.commit_scope) == presentation_before


def test_cancellation_during_materialization_suppresses_accepted_result() -> None:
    session, experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    before = session.analysis_values.snapshot()
    token = CancellationToken()
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)
    cancel_materialization = False

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        result = original(*args, **kwargs)
        if cancel_materialization:
            token.cancel()
        return result

    def arm_cancellation() -> None:
        nonlocal cancel_materialization
        cancel_materialization = True

    def converge_then_cancel_materialization(
        fun: Callable[[Array], Array], x0: Array, **_settings: object
    ) -> object:
        return _one_request_converged_backend(
            fun,
            x0,
            after_request=arm_cancellation,
        )

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            converge_then_cancel_materialization,
        ),
    ):
        outcome = execute_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert outcome.terminal is DirectTrfOutcomeTerminal.CANCELLED
    assert outcome.execution.terminal is DirectTrfTerminal.CONVERGED
    assert outcome.materialization is not None
    assert outcome.materialization.terminal is MaterializationTerminal.CANCELLED
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_stale_commit_rejects_the_complete_accepted_scope_atomically() -> None:
    session, _experiments, parameterization, engine, problem, invocation = (
        _qualification_fit()
    )
    outcome = execute_direct_trf(problem, invocation, parameterization, engine)
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    start = session.analysis_values.snapshot()
    session.analysis_values.commit({"__PB": 0.01}, expected=start, scope=("__PB",))
    advanced = session.analysis_values.snapshot()

    with pytest.raises(StaleAnalysisValuesError):
        commit_accepted_fit(
            outcome.accepted_result,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            outcome.accepted_result,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert session.analysis_values.snapshot() == advanced
    assert advanced["__PB"] == 0.01
    assert advanced[problem.controlled_ids[0]] == start[problem.controlled_ids[0]]
