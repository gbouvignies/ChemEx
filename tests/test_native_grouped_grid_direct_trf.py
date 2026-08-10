"""Behavioral qualification for exact grouped native GRID -> TRF (#596)."""

from __future__ import annotations

import dataclasses
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    CancellationToken,
    DirectTrfConstructionError,
    OptimizationProblem,
)
from chemex.optimize.grid_direct_trf import (
    GridDirectTrfInvocation,
    GridSeedDisposition,
)
from chemex.optimize.grouped_direct_trf import ComponentDisposition, FitDecomposition
from chemex.optimize.grouped_grid_direct_trf import (
    GroupedGridDirectTrfTerminal,
    GroupedGridSeedOutcome,
    commit_grouped_grid_accepted_fit,
    execute_grouped_grid_direct_trf,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


def _grouped_grid_problem() -> tuple[
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
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameterization = session.compile_parameterization(
        read_methods([METHOD])["DEFAULT"],
        experiments.param_ids,
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
    return session, experiments, parameterization, engine, problem


def _backend_result(
    candidate: Array,
    residuals: Array,
    *,
    status: int = 1,
    success: bool = True,
) -> object:
    return SimpleNamespace(
        x=candidate,
        fun=residuals,
        status=status,
        success=success,
        message="qualified grouped GRID backend result",
        nfev=1,
        njev=0,
        cost=0.5 * float(np.dot(residuals, residuals)),
        optimality=0.0,
        active_mask=np.zeros_like(candidate, dtype=np.int64),
    )


def _component_objective(
    attempt: GroupedGridSeedOutcome,
    component_ordinal: int,
) -> float:
    candidate = attempt.components[component_ordinal].candidate
    assert candidate is not None
    return candidate.chi_square


def test_each_seed_uses_exact_components_and_only_selected_aggregate_commits() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    (axis_id,) = problem.controlled_ids[:1]
    axis_start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (axis_start, axis_start * 1.1)),),
        objective_request_budget=10,
    )
    before = session.analysis_values.snapshot()

    def successful_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate))

    with patch("chemex.optimize.direct_trf.least_squares", successful_backend):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedGridDirectTrfTerminal.ACCEPTED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.ELIGIBLE,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.selection is not None
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    assert session.analysis_values.snapshot() == before

    for seed, attempt in zip(invocation.seeds, outcome.attempts, strict=True):
        assert seed.problem is not None
        assert seed.problem.source_snapshot is problem.source_snapshot
        seed_decomposition = FitDecomposition.from_root(
            seed.problem,
            parameterization,
            engine,
        )
        assert attempt.seed_decomposition is not None
        assert attempt.seed_decomposition.root_problem_identity == seed.problem.identity
        assert attempt.seed_decomposition.identity == seed_decomposition.identity
        assert seed_decomposition.identity == decomposition.identity
        assert tuple(item.component_identity for item in attempt.components) == tuple(
            component.identity for component in decomposition.components
        )
        assert tuple(
            item.candidate.problem_identity
            for item in attempt.components
            if item.candidate is not None
        ) == tuple(
            component.problem.identity for component in seed_decomposition.components
        )
        assert all(not hasattr(item, "accepted_result") for item in attempt.components)
        assert not hasattr(attempt, "accepted_result")
        assert attempt.candidate is not None
        assert attempt.candidate.candidate.problem_identity == problem.identity

    selected = outcome.selection.selected_record
    assert selected.candidate is not None
    assert outcome.accepted_result.vector == selected.candidate.vector
    assert outcome.accepted_result.chi_square == selected.candidate.objective
    assert outcome.accepted_result.problem_identity == problem.identity

    receipt = commit_grouped_grid_accepted_fit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )

    assert receipt.old_revision == before.revision
    assert receipt.new_revision == before.revision + 1


def test_selection_orders_whole_aggregates_and_never_combines_marginal_minima() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    component_count = len(decomposition.components)
    (axis_id,) = problem.controlled_ids[:1]
    axis_start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (axis_start, axis_start * 1.1)),),
        objective_request_budget=10,
    )
    backend_calls = 0

    def crossed_component_backend(
        fun: Callable[[Array], Array],
        x0: Array,
        **kwargs: object,
    ) -> object:
        nonlocal backend_calls
        call_ordinal = backend_calls
        backend_calls += 1
        seed_ordinal, component_ordinal = divmod(call_ordinal, component_count)
        lower, upper = kwargs["bounds"]
        lower_array = np.asarray(lower, dtype=np.float64)
        upper_array = np.asarray(upper, dtype=np.float64)
        start = np.asarray(x0, dtype=np.float64)
        low = np.minimum(np.maximum(start * 0.8, lower_array), upper_array)
        high = np.minimum(np.maximum(start * 1.2, lower_array), upper_array)
        low_residuals = fun(low)
        high_residuals = fun(high)
        low_objective = float(np.dot(low_residuals, low_residuals))
        high_objective = float(np.dot(high_residuals, high_residuals))
        better, worse = (low, high) if low_objective < high_objective else (high, low)
        if component_ordinal == 0:
            candidate = better if seed_ordinal == 0 else worse
        elif component_ordinal == 1:
            candidate = worse if seed_ordinal == 0 else better
        else:
            candidate = better
        residuals = fun(candidate)
        return _backend_result(candidate, residuals)

    with patch("chemex.optimize.direct_trf.least_squares", crossed_component_backend):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedGridDirectTrfTerminal.ACCEPTED
    assert outcome.selection is not None
    assert outcome.selection.ordering_policy == (
        "aggregate-chi-square-vector-seed-ordinal-v1"
    )
    assert all(attempt.candidate is not None for attempt in outcome.attempts)

    aggregate_vectors: list[tuple[float, ...]] = []
    for attempt in outcome.attempts:
        assignments = {
            param_id: value
            for component in attempt.components
            if component.candidate is not None
            for param_id, value in zip(
                component.controlled_ids,
                component.candidate.vector,
                strict=True,
            )
        }
        expected = tuple(assignments[param_id] for param_id in problem.controlled_ids)
        assert attempt.candidate is not None
        assert attempt.candidate.vector == expected
        aggregate_vectors.append(expected)

    marginal_sources: list[int] = []
    marginal_assignments: dict[str, float] = {}
    for component_ordinal, component in enumerate(decomposition.components):
        source_ordinal, source = min(
            enumerate(outcome.attempts),
            key=lambda item: _component_objective(item[1], component_ordinal),
        )
        marginal_sources.append(source_ordinal)
        source_candidate = source.components[component_ordinal].candidate
        assert source_candidate is not None
        marginal_assignments.update(
            zip(component.controlled_ids, source_candidate.vector, strict=True)
        )
    marginal_vector = tuple(
        marginal_assignments[param_id] for param_id in problem.controlled_ids
    )

    assert set(marginal_sources[:2]) == {0, 1}
    assert marginal_vector not in aggregate_vectors
    selected = outcome.selection.selected_record
    assert selected.candidate is not None
    assert selected.candidate.vector in aggregate_vectors
    assert selected.candidate.objective == min(
        attempt.candidate.objective
        for attempt in outcome.attempts
        if attempt.candidate is not None
    )


def test_commit_rejects_tampered_component_lineage_before_generic_authority() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    (axis_id,) = problem.controlled_ids[:1]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (problem.start[0],)),),
        objective_request_budget=10,
    )

    def successful_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate))

    with patch("chemex.optimize.direct_trf.least_squares", successful_backend):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.selection is not None
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    selected = outcome.selection.selected_record
    first_component = selected.components[0]
    assert first_component.candidate is not None
    foreign_candidate = dataclasses.replace(
        first_component.candidate,
        problem_identity="foreign-component-problem",
    )
    foreign_component = dataclasses.replace(
        first_component,
        candidate=foreign_candidate,
    )
    foreign_selected = dataclasses.replace(
        selected,
        components=(foreign_component, *selected.components[1:]),
    )
    selection_type = type(outcome.selection)
    foreign_selection = selection_type.from_candidate_records((foreign_selected,))
    foreign_provenance = dataclasses.replace(
        outcome.accepted_result.grouped_grid_provenance,
        selection_identity=foreign_selection.identity,
        selected_outcome_identity=foreign_selected.identity,
    )
    foreign_accepted = dataclasses.replace(
        outcome.accepted_result,
        grouped_grid_provenance=foreign_provenance,
        selection=foreign_selection,
        selected_outcome=foreign_selected,
        origin_context_identity=foreign_provenance.identity,
    )
    before = session.analysis_values.snapshot()

    with pytest.raises(
        DirectTrfConstructionError,
        match="canonical aggregate authority",
    ):
        commit_grouped_grid_accepted_fit(
            foreign_accepted,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert session.analysis_values.snapshot() == before


def test_cancellation_stops_later_seeds_without_aggregate_authority() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    (axis_id,) = problem.controlled_ids[:1]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (problem.start[0], problem.start[0] * 1.1)),),
        objective_request_budget=10,
    )
    token = CancellationToken()
    backend_calls = 0

    def cancelling_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls
        backend_calls += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        token.cancel()
        return _backend_result(candidate, residuals)

    before = session.analysis_values.snapshot()
    with patch("chemex.optimize.direct_trf.least_squares", cancelling_backend):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert backend_calls == 1
    assert outcome.terminal is GroupedGridDirectTrfTerminal.CANCELLED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.CANCELLED,
        GridSeedDisposition.NOT_STARTED,
    )
    assert outcome.selection is None
    assert outcome.accepted_result is None
    assert outcome.commit_authority is None
    assert session.analysis_values.snapshot() == before


def test_failed_seed_keeps_component_evidence_and_later_seed_can_win() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    (axis_id,) = problem.controlled_ids[:1]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (problem.start[0], problem.start[0] * 1.1)),),
        objective_request_budget=10,
    )
    backend_calls = 0

    def first_component_does_not_converge(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls
        backend_calls += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        if backend_calls == 1:
            return _backend_result(
                candidate,
                residuals,
                status=0,
                success=False,
            )
        return _backend_result(candidate, residuals)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        first_component_does_not_converge,
    ):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert backend_calls == 2 * len(decomposition.components)
    assert outcome.terminal is GroupedGridDirectTrfTerminal.ACCEPTED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.NON_CONVERGED,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.attempts[0].candidate is None
    assert outcome.attempts[0].components[0].disposition is (
        ComponentDisposition.FAILED
    )
    assert outcome.selection is not None
    assert outcome.selection.selected_seed_ordinal == 1


def test_equal_aggregate_objective_and_vector_tie_selects_first_seed() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_grid_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    (axis_id,) = problem.controlled_ids[:1]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, (problem.start[0], problem.start[0] * 1.1)),),
        objective_request_budget=10,
    )
    starts = dict(zip(problem.controlled_ids, problem.start, strict=True))
    backend_calls = 0

    def tied_backend(
        fun: Callable[[Array], Array], _x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls
        component = decomposition.components[
            backend_calls % len(decomposition.components)
        ]
        backend_calls += 1
        candidate = np.asarray(
            tuple(starts[param_id] for param_id in component.controlled_ids),
            dtype=np.float64,
        )
        return _backend_result(candidate, fun(candidate))

    with patch("chemex.optimize.direct_trf.least_squares", tied_backend):
        outcome = execute_grouped_grid_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedGridDirectTrfTerminal.ACCEPTED
    assert outcome.selection is not None
    first, second = outcome.attempts
    assert first.candidate is not None
    assert second.candidate is not None
    assert first.candidate.objective == second.candidate.objective
    assert first.candidate.vector == second.candidate.vector
    assert outcome.selection.selected_seed_ordinal == 0
