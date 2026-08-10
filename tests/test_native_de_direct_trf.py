"""Behavioral qualification for selected-coordinate native DE -> TRF (#597).

The public seams are the immutable two-stage workflow invocation, typed search
and polish outcomes, and the existing revision-checked Direct TRF commit
boundary.  Legacy optimizer dispatch remains outside this qualification.
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest
from scipy.optimize import Bounds

from chemex.configuration.methods import Method, Selection, read_methods
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
from chemex.optimize.de_direct_trf import (
    DeCoordinateSemantics,
    DeDirectTrfInterrupted,
    DeDirectTrfInvocation,
    DeDirectTrfTerminal,
    DeSearchCoordinate,
    DeSearchTerminal,
    commit_de_accepted_fit,
    execute_de_direct_trf,
)
from chemex.optimize.direct_trf import (
    CancellationToken,
    DirectTrfConstructionError,
    OptimizationProblem,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


def _qualification_problem(
    method: Method,
) -> tuple[
    AnalysisSession,
    Experiments,
    ActiveParameterization,
    EvaluationEngine,
    OptimizationProblem,
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    return session, experiments, parameterization, engine, problem


def test_selected_coordinate_plan_is_bounded_seeded_and_non_authoritative() -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    pb_id, r1_id = problem.controlled_ids
    pb_index = problem.controlled_ids.index(pb_id)
    search_lower = problem.lower_bounds[pb_index]
    search_upper = problem.upper_bounds[pb_index]

    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=(
            (pb_id, search_lower, search_upper, DeCoordinateSemantics.LINEAR),
        ),
        root_seed=597,
        de_objective_request_budget=6,
        polish_objective_request_budget=5,
        population_multiplier=4,
        maximum_generations=3,
    )

    assert invocation.root_problem_identity == problem.identity
    assert tuple(item.param_id for item in invocation.search_coordinates) == (pb_id,)
    assert invocation.search_coordinates[0].physical_bounds == (
        search_lower,
        search_upper,
    )
    assert invocation.search_problem.controlled_ids == (pb_id,)
    assert dict(invocation.search_problem.held_items)[r1_id] == problem.start[1]
    assert invocation.search_problem.lower_bounds == (problem.lower_bounds[0],)
    assert invocation.search_problem.upper_bounds == (problem.upper_bounds[0],)
    assert not invocation.search_problem.acceptance_authority
    assert invocation.root_seed == 597
    assert invocation.population.dimension == 1
    assert invocation.population.multiplier == 4
    assert invocation.population.size == 5
    assert invocation.population.maximum_generations == 3
    assert invocation.de_objective_request_budget == 6
    assert invocation.polish_objective_request_budget == 5
    assert (
        invocation.identity
        == DeDirectTrfInvocation.for_problem(
            problem,
            search_coordinates=((pb_id, search_lower, search_upper, "linear"),),
            root_seed=597,
            de_objective_request_budget=6,
            polish_objective_request_budget=5,
            population_multiplier=4,
            maximum_generations=3,
        ).identity
    )


@pytest.mark.parametrize("root_seed", (-1, 2**64, True, None))
def test_de_requires_an_explicit_unsigned_64_bit_root_seed(root_seed: object) -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    param_id = problem.controlled_ids[0]
    index = problem.controlled_ids.index(param_id)

    with pytest.raises(DirectTrfConstructionError, match="unsigned 64-bit root seed"):
        DeDirectTrfInvocation.for_problem(
            problem,
            search_coordinates=(
                (
                    param_id,
                    problem.lower_bounds[index],
                    problem.upper_bounds[index],
                    "linear",
                ),
            ),
            root_seed=root_seed,
            de_objective_request_budget=5,
            polish_objective_request_budget=5,
        )


def test_log_coordinate_map_has_exact_endpoints_and_stable_interior_round_trip() -> (
    None
):
    coordinate = DeSearchCoordinate("log-rate", 0.1, 10.0, "log", 0)
    solver_lower, solver_upper = coordinate.solver_bounds

    assert coordinate.to_physical(solver_lower) == 0.1
    assert coordinate.to_physical(solver_upper) == 10.0
    assert coordinate.to_physical(coordinate.to_solver(2.5)) == pytest.approx(
        2.5,
        rel=2.0e-15,
        abs=0.0,
    )


@pytest.mark.parametrize(
    ("lower", "upper"),
    ((0.0, 1.0), (-1.0, 1.0), (1.0, math.inf), (1.0, 1.0)),
)
def test_log_coordinate_rejects_non_positive_non_finite_or_empty_bounds(
    lower: float,
    upper: float,
) -> None:
    with pytest.raises(DirectTrfConstructionError):
        DeSearchCoordinate("log-rate", lower, upper, "log", 0)


def test_log_coordinate_rejects_a_start_outside_the_declared_search_range() -> None:
    coordinate = DeSearchCoordinate("log-rate", 0.1, 10.0, "log", 0)

    with pytest.raises(DirectTrfConstructionError, match="outside its search range"):
        coordinate.to_solver(10.1)


def _trf_backend_result(candidate: Array, residuals: Array) -> SimpleNamespace:
    return SimpleNamespace(
        x=candidate,
        fun=residuals,
        status=1,
        success=True,
        message="qualified TRF polish",
        nfev=1,
        njev=0,
        cost=0.5 * float(np.dot(residuals, residuals)),
        optimality=0.0,
        active_mask=np.zeros_like(candidate, dtype=np.int64),
    )


def _bounded_de_invocation(
    problem: OptimizationProblem,
    *,
    de_budget: int = 6,
) -> DeDirectTrfInvocation:
    (param_id,) = problem.controlled_ids
    lower = max(problem.lower_bounds[0], problem.start[0] * 0.5)
    upper = min(problem.upper_bounds[0], problem.start[0] * 1.5)
    return DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=((param_id, lower, upper, "linear"),),
        root_seed=597,
        de_objective_request_budget=de_budget,
        polish_objective_request_budget=5,
        population_multiplier=4,
        maximum_generations=2,
    )


def _eligible_search_backend(
    fun: Callable[[Array], float],
    _bounds: Bounds,
    **kwargs: object,
) -> object:
    candidate = np.asarray(kwargs["x0"], dtype=np.float64)
    objective = fun(candidate)
    population = np.repeat(candidate[None, :], 5, axis=0)
    return SimpleNamespace(
        x=candidate,
        fun=objective,
        success=True,
        message="Optimization terminated successfully.",
        nit=1,
        nfev=1,
        population=population,
        population_energies=np.full(5, objective),
    )


def _successful_polish_backend(
    fun: Callable[[Array], Array],
    x0: Array,
    **_kwargs: object,
) -> object:
    candidate = np.asarray(x0, dtype=np.float64)
    return _trf_backend_result(candidate, fun(candidate))


def test_de_selects_canonical_restart_then_runs_one_full_trf_polish_and_commit() -> (
    None
):
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    search_lower = max(problem.lower_bounds[0], problem.start[0] * 0.5)
    search_upper = min(problem.upper_bounds[0], problem.start[0] * 1.5)
    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=((param_id, search_lower, search_upper, "linear"),),
        root_seed=597,
        de_objective_request_budget=6,
        polish_objective_request_budget=5,
        population_multiplier=4,
        maximum_generations=1,
    )
    de_calls = 0
    polish_starts: list[tuple[float, ...]] = []

    def bounded_search_backend(
        fun: Callable[[Array], float],
        bounds: Bounds,
        **kwargs: object,
    ) -> object:
        nonlocal de_calls
        de_calls += 1
        assert kwargs["strategy"] == "best1bin"
        assert kwargs["init"] == "latinhypercube"
        assert kwargs["updating"] == "deferred"
        assert kwargs["polish"] is False
        assert kwargs["workers"] == 1
        assert kwargs["vectorized"] is False
        assert isinstance(kwargs["rng"], np.random.Generator)
        assert kwargs["x0"] is not None
        lower = float(bounds.lb[0])
        upper = float(bounds.ub[0])
        first = np.asarray(kwargs["x0"], dtype=np.float64)
        second = np.asarray(((lower + upper) * 0.5,), dtype=np.float64)
        first_objective = fun(first)
        second_objective = fun(second)
        # Deliberately return the worse backend x: ChemEx owns restart ordering.
        returned = first if first_objective >= second_objective else second
        return SimpleNamespace(
            x=returned,
            fun=max(first_objective, second_objective),
            success=False,
            message="Maximum number of iterations has been exceeded.",
            nit=1,
            nfev=2,
            population=np.asarray((first, second, first, second, first)),
            population_energies=np.asarray(
                (
                    first_objective,
                    second_objective,
                    first_objective,
                    second_objective,
                    first_objective,
                )
            ),
        )

    def polishing_backend(
        fun: Callable[[Array], Array],
        x0: Array,
        **_kwargs: object,
    ) -> object:
        polish_starts.append(tuple(float(value) for value in x0))
        candidate = np.asarray(x0, dtype=np.float64)
        residuals = fun(candidate)
        return _trf_backend_result(candidate, residuals)

    before = session.analysis_values.snapshot()
    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            bounded_search_backend,
        ),
        patch("chemex.optimize.direct_trf.least_squares", polishing_backend),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert de_calls == 1
    assert outcome.terminal is DeDirectTrfTerminal.ACCEPTED
    assert outcome.search.terminal is DeSearchTerminal.GENERATION_LIMIT
    assert outcome.search.restart_eligible
    assert outcome.search.best_candidate is not None
    assert polish_starts == [outcome.search.best_candidate.full_vector]
    assert outcome.polish_problem is not None
    assert outcome.polish_problem.controlled_ids == problem.controlled_ids
    assert outcome.polish_problem.start == outcome.search.best_candidate.full_vector
    assert not outcome.polish_problem.acceptance_authority
    assert outcome.polish is not None
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    assert outcome.accounting.de_budget == 6
    assert outcome.accounting.de_counters.objective_requests_accepted == 2
    assert outcome.accounting.polish_budget == 5
    assert outcome.accounting.polish_counters is not None
    assert outcome.accounting.search_to_polish_transfers == 1
    assert outcome.accounting.root_materializations == 1
    assert session.analysis_values.snapshot() == before

    receipt = commit_de_accepted_fit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )

    assert receipt.old_revision == before.revision
    assert receipt.new_revision == before.revision + 1


def test_de_interruption_freezes_typed_evidence_and_never_starts_polish() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    search_lower = max(problem.lower_bounds[0], problem.start[0] * 0.5)
    search_upper = min(problem.upper_bounds[0], problem.start[0] * 1.5)
    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=((param_id, search_lower, search_upper, "linear"),),
        root_seed=597,
        de_objective_request_budget=6,
        polish_objective_request_budget=5,
        population_multiplier=4,
    )

    def interrupted_backend(
        fun: Callable[[Array], float],
        _bounds: Bounds,
        **kwargs: object,
    ) -> object:
        fun(np.asarray(kwargs["x0"], dtype=np.float64))
        raise KeyboardInterrupt

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            interrupted_backend,
        ),
        patch("chemex.optimize.direct_trf.least_squares") as polish_backend,
        pytest.raises(DeDirectTrfInterrupted) as interrupted,
    ):
        execute_de_direct_trf(problem, invocation, parameterization, engine)

    outcome = interrupted.value.outcome
    assert outcome.terminal is DeDirectTrfTerminal.INTERRUPTED
    assert outcome.search.terminal is DeSearchTerminal.INTERRUPTED
    assert outcome.search.counters.objective_requests_accepted == 1
    assert outcome.polish is None
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None
    assert outcome.accounting.search_to_polish_transfers == 0
    polish_backend.assert_not_called()


def test_de_preflight_interruption_is_typed_before_backend_entry() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)

    with (
        patch.object(BoundEvaluator, "evaluate", side_effect=KeyboardInterrupt),
        patch("chemex.optimize.de_direct_trf.differential_evolution") as search_backend,
        pytest.raises(DeDirectTrfInterrupted) as interrupted,
    ):
        execute_de_direct_trf(problem, invocation, parameterization, engine)

    outcome = interrupted.value.outcome
    assert outcome.terminal is DeDirectTrfTerminal.INTERRUPTED
    assert outcome.search.terminal is DeSearchTerminal.INTERRUPTED
    assert outcome.search.preflight_evaluation_identity is None
    assert outcome.search.counters.objective_requests_accepted == 0
    assert outcome.polish is None
    assert outcome.commit_authority is None
    search_backend.assert_not_called()


def test_initial_de_evaluator_binding_interruption_is_typed() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)

    with (
        patch.object(engine, "new_evaluator", side_effect=KeyboardInterrupt),
        patch("chemex.optimize.de_direct_trf.differential_evolution") as search_backend,
        pytest.raises(DeDirectTrfInterrupted) as interrupted,
    ):
        execute_de_direct_trf(problem, invocation, parameterization, engine)

    outcome = interrupted.value.outcome
    assert outcome.terminal is DeDirectTrfTerminal.INTERRUPTED
    assert outcome.search.terminal is DeSearchTerminal.INTERRUPTED
    assert outcome.search.counters.objective_requests_accepted == 0
    assert outcome.search.preflight_evaluation_identity is None
    assert outcome.polish is None
    search_backend.assert_not_called()


def test_initial_de_evaluator_binding_failure_is_typed_without_backend_entry() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)

    with (
        patch.object(
            engine,
            "new_evaluator",
            side_effect=RuntimeError("binding failed"),
        ),
        patch("chemex.optimize.de_direct_trf.differential_evolution") as search_backend,
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is DeDirectTrfTerminal.EXECUTION_FAILURE
    assert outcome.search.terminal is DeSearchTerminal.IMPLEMENTATION_FAILURE
    assert outcome.search.failure is not None
    assert outcome.search.failure.category == "de_evaluator_binding_failure"
    assert outcome.search.counters.objective_requests_accepted == 0
    assert outcome.polish is None
    search_backend.assert_not_called()


def test_root_materialization_binding_interruption_is_typed() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    original_new_evaluator = engine.new_evaluator
    evaluator_count = 0

    def interrupted_root_binding() -> BoundEvaluator:
        nonlocal evaluator_count
        evaluator_count += 1
        if evaluator_count == 4:
            raise KeyboardInterrupt
        return original_new_evaluator()

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            _eligible_search_backend,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            _successful_polish_backend,
        ),
        patch.object(engine, "new_evaluator", interrupted_root_binding),
        pytest.raises(DeDirectTrfInterrupted) as interrupted,
    ):
        execute_de_direct_trf(problem, invocation, parameterization, engine)

    outcome = interrupted.value.outcome
    assert evaluator_count == 4
    assert outcome.terminal is DeDirectTrfTerminal.INTERRUPTED
    assert outcome.root_materialization is not None
    assert outcome.root_materialization.terminal.value == "interrupted"
    assert outcome.root_materialization.evaluation_count == 0
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None


def test_trf_polish_interruption_has_no_acceptance_or_backend_fallback() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    search_lower = max(problem.lower_bounds[0], problem.start[0] * 0.5)
    search_upper = min(problem.upper_bounds[0], problem.start[0] * 1.5)
    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=((param_id, search_lower, search_upper, "linear"),),
        root_seed=597,
        de_objective_request_budget=6,
        polish_objective_request_budget=5,
        population_multiplier=4,
    )
    polish_calls = 0

    def eligible_search(
        fun: Callable[[Array], float],
        _bounds: Bounds,
        **kwargs: object,
    ) -> object:
        candidate = np.asarray(kwargs["x0"], dtype=np.float64)
        objective = fun(candidate)
        population = np.repeat(candidate[None, :], 5, axis=0)
        return SimpleNamespace(
            x=candidate,
            fun=objective,
            success=True,
            message="Optimization terminated successfully.",
            nit=1,
            nfev=1,
            population=population,
            population_energies=np.full(5, objective),
        )

    def interrupted_polish(
        fun: Callable[[Array], Array],
        x0: Array,
        **_kwargs: object,
    ) -> object:
        nonlocal polish_calls
        polish_calls += 1
        fun(np.asarray(x0, dtype=np.float64))
        raise KeyboardInterrupt

    with (
        patch("chemex.optimize.de_direct_trf.differential_evolution", eligible_search),
        patch("chemex.optimize.direct_trf.least_squares", interrupted_polish),
        pytest.raises(DeDirectTrfInterrupted) as interrupted,
    ):
        execute_de_direct_trf(problem, invocation, parameterization, engine)

    outcome = interrupted.value.outcome
    assert polish_calls == 1
    assert outcome.terminal is DeDirectTrfTerminal.INTERRUPTED
    assert outcome.search.restart_eligible
    assert outcome.polish is not None
    assert outcome.polish.execution.terminal.value == "interrupted"
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None
    assert outcome.accounting.search_to_polish_transfers == 1
    assert outcome.accounting.root_materializations == 0


def test_all_invalid_de_population_has_no_restart_or_polish_authority() -> None:
    _session, experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)
    reject_trials = False

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        calculated = original(*args, **kwargs)
        if reject_trials:
            return np.full_like(calculated, np.nan)
        return calculated

    def rejected_population(
        fun: Callable[[Array], float],
        bounds: Bounds,
        **kwargs: object,
    ) -> object:
        nonlocal reject_trials
        candidate = np.asarray(bounds.ub, dtype=np.float64)
        reject_trials = True
        objectives = tuple(fun(candidate) for _ in range(5))
        assert all(math.isinf(value) for value in objectives)
        population = np.repeat(candidate[None, :], 5, axis=0)
        return SimpleNamespace(
            x=candidate,
            fun=math.inf,
            success=False,
            message="Maximum number of iterations has been exceeded.",
            nit=1,
            nfev=5,
            population=population,
            population_energies=np.full(5, math.inf),
        )

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            rejected_population,
        ),
        patch("chemex.optimize.direct_trf.least_squares") as polish_backend,
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is DeDirectTrfTerminal.SEARCH_UNSUCCESSFUL
    assert outcome.search.terminal is DeSearchTerminal.NO_VALID_CANDIDATE
    assert outcome.search.best_candidate is None
    assert outcome.search.valid_candidate_count == 0
    assert outcome.search.rejected_trial_count == 5
    assert outcome.search.counters.objective_evaluations_completed == 5
    assert outcome.polish is None
    assert outcome.commit_authority is None
    polish_backend.assert_not_called()


def test_corrupt_successful_evaluation_is_fatal_not_a_rejected_trial() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    original_new_evaluator = engine.new_evaluator

    class CorruptingEvaluator:
        def __init__(self, delegate: BoundEvaluator) -> None:
            self._delegate = delegate
            self.compatibility_identity = delegate.compatibility_identity
            self.calls = 0

        @property
        def cache_statistics(self) -> object:
            return self._delegate.cache_statistics

        def evaluate(self, frame: EvaluationFrame) -> object:
            self.calls += 1
            evaluated = self._delegate.evaluate(frame)
            if self.calls == 2 and isinstance(evaluated, EvaluationResult):
                return SimpleNamespace(
                    residuals=evaluated.residuals,
                    resolved_values={},
                    identity="corrupt-successful-evaluation",
                )
            return evaluated

    evaluator_count = 0

    def staged_evaluator() -> BoundEvaluator | CorruptingEvaluator:
        nonlocal evaluator_count
        evaluator_count += 1
        evaluator = original_new_evaluator()
        return CorruptingEvaluator(evaluator) if evaluator_count == 1 else evaluator

    def one_corrupt_request(
        fun: Callable[[Array], float],
        _bounds: Bounds,
        **kwargs: object,
    ) -> object:
        fun(np.asarray(kwargs["x0"], dtype=np.float64))
        raise AssertionError("Corrupt evidence must terminate the DE stage")

    with (
        patch.object(engine, "new_evaluator", staged_evaluator),
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            one_corrupt_request,
        ),
        patch("chemex.optimize.direct_trf.least_squares") as polish_backend,
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is DeDirectTrfTerminal.EXECUTION_FAILURE
    assert outcome.search.terminal is DeSearchTerminal.IMPLEMENTATION_FAILURE
    assert outcome.search.failure is not None
    assert outcome.search.failure.category == "de_candidate_evidence_failure"
    assert outcome.search.rejected_trial_count == 0
    assert outcome.polish is None
    polish_backend.assert_not_called()


def test_de_budget_refusal_is_uncharged_and_best_partial_candidate_is_polished() -> (
    None
):
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem, de_budget=6)

    def exhaust_after_valid_candidates(
        fun: Callable[[Array], float],
        _bounds: Bounds,
        **kwargs: object,
    ) -> object:
        candidate = np.asarray(kwargs["x0"], dtype=np.float64)
        for _ in range(7):
            fun(candidate)
        raise AssertionError("The seventh request must stop the DE backend")

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            exhaust_after_valid_candidates,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            _successful_polish_backend,
        ),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is DeDirectTrfTerminal.ACCEPTED
    assert outcome.search.terminal is DeSearchTerminal.BUDGET_EXHAUSTED
    assert outcome.search.restart_eligible
    assert outcome.search.counters.solver_requests_received == 7
    assert outcome.search.counters.objective_requests_accepted == 6
    assert outcome.search.counters.objective_evaluations_completed == 6
    assert outcome.accounting.search_to_polish_transfers == 1


def test_de_cancellation_stops_before_polish_without_fallback() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    token = CancellationToken()

    def cancel_after_one_candidate(
        fun: Callable[[Array], float],
        _bounds: Bounds,
        **kwargs: object,
    ) -> object:
        candidate = np.asarray(kwargs["x0"], dtype=np.float64)
        fun(candidate)
        token.cancel()
        fun(candidate)
        raise AssertionError("Cancellation must stop the DE backend")

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            cancel_after_one_candidate,
        ),
        patch("chemex.optimize.direct_trf.least_squares") as polish_backend,
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert outcome.terminal is DeDirectTrfTerminal.CANCELLED
    assert outcome.search.terminal is DeSearchTerminal.CANCELLED
    assert outcome.search.best_candidate is not None
    assert not outcome.search.restart_eligible
    assert outcome.polish is None
    assert outcome.accepted_result is None
    polish_backend.assert_not_called()


def test_non_converged_polish_is_terminal_and_never_falls_back_to_de_candidate() -> (
    None
):
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    polish_calls = 0

    def non_converged_polish(
        fun: Callable[[Array], Array],
        x0: Array,
        **_kwargs: object,
    ) -> object:
        nonlocal polish_calls
        polish_calls += 1
        candidate = np.asarray(x0, dtype=np.float64)
        residuals = fun(candidate)
        result = _trf_backend_result(candidate, residuals)
        result.status = 0
        result.success = False
        return result

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            _eligible_search_backend,
        ),
        patch("chemex.optimize.direct_trf.least_squares", non_converged_polish),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert polish_calls == 1
    assert outcome.terminal is DeDirectTrfTerminal.POLISH_UNSUCCESSFUL
    assert outcome.search.restart_eligible
    assert outcome.polish is not None
    assert outcome.polish.execution.terminal.value == "non_converged"
    assert outcome.root_materialization is None
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None
    assert outcome.accounting.root_materializations == 0


def test_root_materialization_failure_never_accepts_or_retries_polish() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    original_new_evaluator = engine.new_evaluator
    evaluator_count = 0
    polish_calls = 0

    class FailingEvaluator:
        def __init__(self, delegate: BoundEvaluator) -> None:
            self.compatibility_identity = delegate.compatibility_identity
            self.cache_statistics = delegate.cache_statistics

        def evaluate(self, _frame: object) -> EvaluationFailure:
            return EvaluationFailure(
                engine.plan.identity,
                parameterization.evaluator_identity,
                "kernel",
                "non_finite_calculation",
                "INVALID_TRIAL",
            )

    def staged_evaluator() -> BoundEvaluator | FailingEvaluator:
        nonlocal evaluator_count
        evaluator_count += 1
        evaluator = original_new_evaluator()
        return FailingEvaluator(evaluator) if evaluator_count == 4 else evaluator

    def counted_polish(
        fun: Callable[[Array], Array],
        x0: Array,
        **kwargs: object,
    ) -> object:
        nonlocal polish_calls
        polish_calls += 1
        return _successful_polish_backend(fun, x0, **kwargs)

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            _eligible_search_backend,
        ),
        patch("chemex.optimize.direct_trf.least_squares", counted_polish),
        patch.object(engine, "new_evaluator", staged_evaluator),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert evaluator_count == 4
    assert polish_calls == 1
    assert outcome.terminal is DeDirectTrfTerminal.MATERIALIZATION_FAILURE
    assert outcome.root_materialization is not None
    assert outcome.root_materialization.terminal.value == "failure"
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None
    assert outcome.accounting.root_materializations == 1


def test_cancellation_after_root_evaluation_suppresses_acceptance() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    token = CancellationToken()
    original_new_evaluator = engine.new_evaluator
    evaluator_count = 0

    class CancellingEvaluator:
        def __init__(self, delegate: BoundEvaluator) -> None:
            self._delegate = delegate
            self.compatibility_identity = delegate.compatibility_identity

        @property
        def cache_statistics(self) -> object:
            return self._delegate.cache_statistics

        def evaluate(
            self,
            frame: EvaluationFrame,
        ) -> EvaluationResult | EvaluationFailure:
            evaluated = self._delegate.evaluate(frame)
            token.cancel()
            return evaluated

    def staged_evaluator() -> BoundEvaluator | CancellingEvaluator:
        nonlocal evaluator_count
        evaluator_count += 1
        evaluator = original_new_evaluator()
        return CancellingEvaluator(evaluator) if evaluator_count == 4 else evaluator

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            _eligible_search_backend,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            _successful_polish_backend,
        ),
        patch.object(engine, "new_evaluator", staged_evaluator),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert evaluator_count == 4
    assert outcome.terminal is DeDirectTrfTerminal.CANCELLED
    assert outcome.root_materialization is not None
    assert outcome.root_materialization.terminal.value == "cancelled"
    assert outcome.root_materialization.evaluation_count == 1
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None


def test_invocation_is_serializable_and_commit_rejects_tampered_lineage() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem)
    record = invocation.to_record()
    restored = json.loads(json.dumps(record, allow_nan=False, sort_keys=True))
    assert restored == record
    assert restored["identity"] == invocation.identity
    search_record = cast("dict[str, object]", restored["search_problem"])
    assert search_record["controlled_ids"] == list(
        invocation.search_problem.controlled_ids
    )
    assert search_record["captured_held_items"] == [
        [param_id, value] for param_id, value in invocation.search_problem.held_items
    ]
    assert search_record["physical_start"] == list(invocation.search_problem.start)
    assert search_record["solver_x0"] == [
        coordinate.to_solver(value)
        for coordinate, value in zip(
            invocation.search_coordinates,
            invocation.search_problem.start,
            strict=True,
        )
    ]

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            _eligible_search_backend,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            _successful_polish_backend,
        ),
    ):
        outcome = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    accepted = outcome.accepted_result
    invalid_provenance = dataclasses.replace(
        accepted.de_polish_provenance,
        search_candidate_identity="foreign-de-candidate",
    )
    invalid = dataclasses.replace(
        accepted,
        de_polish_provenance=invalid_provenance,
    )
    before = session.analysis_values.snapshot()

    with pytest.raises(
        DirectTrfConstructionError,
        match="canonical polish authority",
    ):
        commit_de_accepted_fit(
            invalid,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert session.analysis_values.snapshot() == before


def test_real_scipy_de_is_deterministic_for_the_same_root_seed() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    invocation = _bounded_de_invocation(problem, de_budget=10)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        _successful_polish_backend,
    ):
        first = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )
        second = execute_de_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert first.terminal is second.terminal is DeDirectTrfTerminal.ACCEPTED
    assert first.search.terminal is second.search.terminal
    assert first.search.counters == second.search.counters
    assert first.search.counters.objective_requests_accepted == 10
    assert first.search.best_candidate == second.search.best_candidate
    assert first.search.backend == second.search.backend
    assert first.search.identity == second.search.identity
