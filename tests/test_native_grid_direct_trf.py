"""Behavioral qualification for single-component native GRID -> TRF (#595).

The public seams are the immutable Cartesian workflow invocation, typed seed
outcomes and selection provenance, and the existing revision-checked Direct TRF
commit boundary.  Legacy GRID execution remains outside this qualification.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    AffineHalfSpace,
    CancellationToken,
    DirectTrfConstructionError,
    DirectTrfInvocation,
    OptimizationProblem,
    commit_accepted_fit,
    execute_direct_trf,
)
from chemex.optimize.grid_direct_trf import (
    GridDirectTrfInterrupted,
    GridDirectTrfInvocation,
    GridDirectTrfTerminal,
    GridSeedDisposition,
    GridSelection,
    commit_grid_accepted_fit,
    execute_grid_direct_trf,
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


def test_grid_seed_problems_preserve_root_affine_feasibility() -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    param_id = problem.controlled_ids[0]
    restriction = AffineHalfSpace(
        "grid-root-upper",
        tuple(
            1.0 if independent_id == param_id else 0.0
            for independent_id, _value in problem.independent_items
        ),
        problem.start[0] + 1.0,
    )
    constrained = dataclasses.replace(
        problem,
        affine_half_spaces=(restriction,),
    )
    invocation = GridDirectTrfInvocation.for_problem(
        constrained,
        axes=((param_id, (constrained.start[0],)),),
        objective_request_budget=5,
    )
    (seed,) = invocation.seeds
    assert seed.problem is not None
    assert seed.problem.affine_half_spaces == (restriction,)
    with pytest.raises(DirectTrfConstructionError, match="feasibility"):
        dataclasses.replace(seed.problem, affine_half_spaces=())


def test_cartesian_seed_plan_is_canonical_physical_and_rooted_once() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        Method(fit=["PB"], fix=["KEX_AB"]),
    )
    pb_id, r1_id = problem.controlled_ids
    pb_lower, r1_lower = problem.lower_bounds
    pb_upper, _r1_upper = problem.upper_bounds

    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=(
            (r1_id, (r1_lower, r1_lower)),
            (pb_id, (pb_lower, pb_upper)),
        ),
        objective_request_budget=5,
    )

    assert tuple(axis.param_id for axis in invocation.axes) == (pb_id, r1_id)
    assert tuple(axis.declaration_ordinal for axis in invocation.axes) == (1, 0)
    assert tuple(seed.ordinal for seed in invocation.seeds) == (0, 1, 2, 3)
    assert tuple(seed.axis_items for seed in invocation.seeds) == (
        ((pb_id, pb_lower), (r1_id, r1_lower)),
        ((pb_id, pb_lower), (r1_id, r1_lower)),
        ((pb_id, pb_upper), (r1_id, r1_lower)),
        ((pb_id, pb_upper), (r1_id, r1_lower)),
    )
    assert tuple(seed.start for seed in invocation.seeds) == (
        (pb_lower, r1_lower),
        (pb_lower, r1_lower),
        (pb_upper, r1_lower),
        (pb_upper, r1_lower),
    )
    assert len({seed.identity for seed in invocation.seeds}) == 4
    assert all(seed.problem is not None for seed in invocation.seeds)
    assert all(seed.invocation is not None for seed in invocation.seeds)
    assert all(
        seed.problem is not None
        and seed.problem.source_snapshot is problem.source_snapshot
        and seed.problem.controlled_ids == problem.controlled_ids
        and seed.problem.held_items == problem.held_items
        and seed.problem.lower_bounds == problem.lower_bounds
        and seed.problem.upper_bounds == problem.upper_bounds
        and seed.problem.commit_scope == problem.commit_scope
        and not seed.problem.acceptance_authority
        for seed in invocation.seeds
    )
    assert (
        invocation.identity
        == GridDirectTrfInvocation.for_problem(
            problem,
            axes=(
                (r1_id, (r1_lower, r1_lower)),
                (pb_id, (pb_lower, pb_upper)),
            ),
            objective_request_budget=5,
        ).identity
    )
    first_seed = invocation.seeds[0]
    replacement_seed = invocation.seeds[2]
    assert first_seed.problem is not None
    assert first_seed.invocation is not None
    with pytest.raises(DirectTrfConstructionError, match="no acceptance authority"):
        execute_direct_trf(
            first_seed.problem,
            first_seed.invocation,
            parameterization,
            engine,
        )
    with pytest.raises(DirectTrfConstructionError, match="differs from its seed"):
        dataclasses.replace(
            first_seed,
            problem=replacement_seed.problem,
            invocation=replacement_seed.invocation,
        )


def test_duplicate_axis_declarations_are_rejected_after_canonical_resolution() -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(
            Method(fit=["PB"], fix=["KEX_AB"]),
        )
    )
    param_id = problem.controlled_ids[0]

    with pytest.raises(ValueError, match="Duplicate GRID axis"):
        GridDirectTrfInvocation.for_problem(
            problem,
            axes=((param_id, (0.1,)), (param_id, (0.2,))),
            objective_request_budget=5,
        )


@pytest.mark.parametrize(
    "unusable",
    (np.nan, np.inf, -np.inf, "not-a-real-coordinate"),
    ids=("nan", "positive-infinity", "negative-infinity", "non-real"),
)
def test_unusable_axis_value_is_seed_local_and_later_seeds_run(
    unusable: object,
) -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    axes = ((param_id, (start, unusable, start * 1.1)),)
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=axes,
        objective_request_budget=5,
    )
    backend_calls = 0

    def finite_only_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls
        backend_calls += 1
        assert np.isfinite(x0).all()
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", finite_only_backend):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert backend_calls == 2
    assert tuple(attempt.seed_ordinal for attempt in outcome.attempts) == (0, 1, 2)
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.ELIGIBLE,
        GridSeedDisposition.UNUSABLE_VALUE,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.attempts[1].failure is not None
    assert outcome.attempts[1].failure.category == "grid_seed_unusable_value"
    assert outcome.attempts[2].execution_identity is not None
    assert tuple(row.seed_ordinal for row in outcome.table_rows) == (0, 1, 2)
    assert (
        invocation.identity
        == GridDirectTrfInvocation.for_problem(
            problem,
            axes=axes,
            objective_request_budget=5,
        ).identity
    )


def _backend_result(
    candidate: Array,
    residuals: Array,
    *,
    status: int,
    success: bool,
) -> object:
    return SimpleNamespace(
        x=candidate,
        fun=residuals,
        status=status,
        success=success,
        message="qualified backend result",
        nfev=1,
        njev=0,
        cost=0.5 * float(np.dot(residuals, residuals)),
        optimality=0.0 if success else 1.0,
        active_mask=np.zeros_like(candidate, dtype=np.int64),
    )


def _legacy_values(
    session: AnalysisSession, scope: tuple[str, ...]
) -> dict[str, float]:
    return {
        param_id: parameter.value
        for param_id, parameter in session.parameters.build_lmfit_params(
            set(scope)
        ).items()
    }


def test_failed_seed_continues_and_canonical_tie_commits_only_selected_once() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    (lower,) = problem.lower_bounds
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (lower - 1.0, start * 1.1, start, start)),),
        objective_request_budget=5,
    )
    before = session.analysis_values.snapshot()
    backend_calls = 0

    def mixed_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls
        backend_calls += 1
        if backend_calls == 1:
            candidate = np.asarray(x0, dtype=np.float64)
            return _backend_result(
                candidate,
                fun(candidate),
                status=0,
                success=False,
            )
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(
            candidate,
            fun(candidate),
            status=1,
            success=True,
        )

    with patch("chemex.optimize.direct_trf.least_squares", mixed_backend):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert backend_calls == 3
    assert outcome.terminal is GridDirectTrfTerminal.ACCEPTED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.OUT_OF_BOUNDS,
        GridSeedDisposition.NON_CONVERGED,
        GridSeedDisposition.ELIGIBLE,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.selection is not None
    assert outcome.selection.selected_seed_ordinal == 2
    assert outcome.selection.ordering_policy == "chi-square-vector-seed-ordinal-v1"
    assert len(outcome.selection.eligible_candidate_identities) == 2
    assert outcome.attempts[2].objective == outcome.attempts[3].objective
    assert outcome.attempts[2].candidate is not None
    assert outcome.attempts[3].candidate is not None
    assert (
        outcome.attempts[2].candidate.identity != outcome.attempts[3].candidate.identity
    )
    assert all(not hasattr(attempt, "accepted_result") for attempt in outcome.attempts)
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    assert outcome.accepted_result.problem_identity == problem.identity
    grid_selection = outcome.accepted_result.grid_selection_provenance
    assert grid_selection is not None
    assert grid_selection.workflow_invocation_identity == invocation.identity
    assert grid_selection.selection_identity == outcome.selection.identity
    assert outcome.accepted_result.origin_context_identity == grid_selection.identity
    assert grid_selection.selected_seed_identity == outcome.attempts[2].seed_identity
    assert grid_selection.selected_seed_ordinal == 2
    assert session.analysis_values.snapshot() == before
    legacy_before = _legacy_values(session, problem.commit_scope)

    replaced_provenance = (
        dataclasses.replace(
            grid_selection,
            workflow_invocation_identity="foreign-grid-workflow",
        ),
        dataclasses.replace(
            grid_selection,
            selection_identity="foreign-grid-selection",
        ),
        dataclasses.replace(
            grid_selection,
            selected_seed_identity="foreign-grid-seed",
        ),
        dataclasses.replace(
            grid_selection,
            selected_seed_ordinal=grid_selection.selected_seed_ordinal + 1,
        ),
        dataclasses.replace(
            grid_selection,
            grid_candidate_identity="foreign-grid-candidate",
        ),
        dataclasses.replace(
            grid_selection,
            candidate_execution_identity="foreign-execution",
        ),
        None,
    )
    for invalid_provenance in replaced_provenance:
        invalid_accepted = dataclasses.replace(
            outcome.accepted_result,
            grid_selection_provenance=invalid_provenance,
        )
        with pytest.raises(
            DirectTrfConstructionError,
            match="canonical selection authority",
        ):
            commit_grid_accepted_fit(
                invalid_accepted,
                outcome.commit_authority,
                problem=problem,
                parameterization=parameterization,
                analysis_values=session.analysis_values,
            )
        assert session.analysis_values.snapshot() == before
        assert _legacy_values(session, problem.commit_scope) == legacy_before

    invalid_origin = dataclasses.replace(
        outcome.accepted_result,
        origin_context_identity="foreign-origin-context",
    )
    assert invalid_origin.identity != outcome.accepted_result.identity
    with pytest.raises(
        DirectTrfConstructionError,
        match="canonical selection authority",
    ):
        commit_grid_accepted_fit(
            invalid_origin,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == before
    assert _legacy_values(session, problem.commit_scope) == legacy_before

    receipt = commit_grid_accepted_fit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )

    assert receipt.old_revision == before.revision
    assert receipt.new_revision == before.revision + 1
    committed = session.analysis_values.snapshot()
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live Direct TRF commit authority",
    ):
        commit_grid_accepted_fit(
            outcome.accepted_result,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == committed


def test_plain_base_evidence_cannot_recreate_grid_commit_authority() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start,)),),
        objective_request_budget=5,
    )

    def successful_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", successful_backend):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    accepted = outcome.accepted_result
    plain_evidence = AcceptedFitResult(
        accepted.occurrence_identity,
        accepted.problem_identity,
        accepted.invocation_identity,
        accepted.execution_identity,
        accepted.materialization_identity,
        accepted.parameterization_identity,
        accepted.evaluator_parameterization_identity,
        accepted.source_occurrence_identity,
        accepted.source_revision,
        accepted.controlled_ids,
        accepted.vector,
        accepted.chi_square,
        accepted.evaluation_result,
        accepted.commit_scope,
        accepted.commit_items,
        accepted.origin_context_identity,
        occurrence_witness=None,
    )
    assert plain_evidence.identity == accepted.identity
    replaced_evidence = dataclasses.replace(accepted)
    assert replaced_evidence is not accepted
    assert replaced_evidence.identity == accepted.identity
    before = session.analysis_values.snapshot()
    legacy_before = _legacy_values(session, problem.commit_scope)

    for reconstructed in (plain_evidence, replaced_evidence):
        with pytest.raises(
            DirectTrfConstructionError,
            match="exact live Direct TRF commit authority",
        ):
            commit_accepted_fit(
                reconstructed,
                outcome.commit_authority,
                problem=problem,
                parameterization=parameterization,
                analysis_values=session.analysis_values,
            )
    with pytest.raises(TypeError, match="authority"):
        commit_accepted_fit(
            plain_evidence,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert (
        plain_evidence.identity
        == AcceptedFitResult(
            plain_evidence.occurrence_identity,
            plain_evidence.problem_identity,
            plain_evidence.invocation_identity,
            plain_evidence.execution_identity,
            plain_evidence.materialization_identity,
            plain_evidence.parameterization_identity,
            plain_evidence.evaluator_parameterization_identity,
            plain_evidence.source_occurrence_identity,
            plain_evidence.source_revision,
            plain_evidence.controlled_ids,
            plain_evidence.vector,
            plain_evidence.chi_square,
            plain_evidence.evaluation_result,
            plain_evidence.commit_scope,
            plain_evidence.commit_items,
            plain_evidence.origin_context_identity,
            occurrence_witness=None,
        ).identity
    )
    assert session.analysis_values.snapshot() == before
    assert _legacy_values(session, problem.commit_scope) == legacy_before


def test_grid_commit_rejects_foreign_direct_origin_authority_atomically() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    grid_invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (problem.start[0],)),),
        objective_request_budget=5,
    )
    direct_invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=5,
    )

    def successful_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", successful_backend):
        grid_outcome = execute_grid_direct_trf(
            problem,
            grid_invocation,
            parameterization,
            engine,
        )
        direct_outcome = execute_direct_trf(
            problem,
            direct_invocation,
            parameterization,
            engine,
        )

    assert grid_outcome.accepted_result is not None
    assert grid_outcome.commit_authority is not None
    assert direct_outcome.commit_authority is not None
    before = session.analysis_values.snapshot()
    legacy_before = _legacy_values(session, problem.commit_scope)

    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live Direct TRF commit authority",
    ):
        commit_grid_accepted_fit(
            grid_outcome.accepted_result,
            direct_outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )

    assert session.analysis_values.snapshot() == before
    assert _legacy_values(session, problem.commit_scope) == legacy_before


def test_selection_requires_one_exact_canonical_seed_candidate_record() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start * 1.1)),),
        objective_request_budget=5,
    )

    def tied_backend(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", tied_backend):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.selection is not None
    eligible = tuple(
        attempt
        for attempt in outcome.attempts
        if attempt.disposition is GridSeedDisposition.ELIGIBLE
    )
    winner, other = eligible
    assert winner.candidate is not None

    mismatches = (
        (other.seed_identity, other.seed_ordinal, winner.candidate.identity),
        (winner.seed_identity, other.seed_ordinal, winner.candidate.identity),
        (other.seed_identity, winner.seed_ordinal, winner.candidate.identity),
    )
    for seed_identity, seed_ordinal, candidate_identity in mismatches:
        with pytest.raises(ValueError, match="exact canonical candidate record"):
            GridSelection.from_candidate_records(
                selected_seed_identity=seed_identity,
                selected_seed_ordinal=seed_ordinal,
                selected_candidate_identity=candidate_identity,
                candidate_records=eligible,
            )

    selection = GridSelection.from_candidate_records(
        selected_seed_identity=winner.seed_identity,
        selected_seed_ordinal=winner.seed_ordinal,
        selected_candidate_identity=winner.candidate.identity,
        candidate_records=eligible,
    )
    rendered = dataclasses.replace(outcome, selection=selection)
    selected_rows = tuple(row for row in rendered.table_rows if row.selected)
    selected_points = tuple(
        point for point in rendered.diagnostic_points((param_id,)) if point.selected
    )
    assert len(selected_rows) == len(selected_points) == 1
    assert selected_rows[0].seed_identity == winner.seed_identity
    assert selected_rows[0].candidate_identity == winner.candidate.identity
    assert selected_points[0].seed_identity == winner.seed_identity
    assert selected_points[0].candidate_identity == winner.candidate.identity


def test_invalid_trial_continues_to_later_seed() -> None:
    session, experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start * 1.1)),),
        objective_request_budget=5,
    )
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate
    reject_trial = False
    backend_calls = 0

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        calculated = original(*args, **kwargs)
        if reject_trial:
            return np.full_like(calculated, np.nan)
        return calculated

    def invalid_then_converged(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal backend_calls, reject_trial
        backend_calls += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        if backend_calls == 1:
            reject_trial = True
            try:
                fun(candidate)
            finally:
                reject_trial = False
            raise AssertionError("invalid trial must stop the first backend")
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch("chemex.optimize.direct_trf.least_squares", invalid_then_converged),
    ):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert backend_calls == 2
    assert outcome.terminal is GridDirectTrfTerminal.ACCEPTED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.INVALID_TRIAL,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.attempts[0].candidate is None
    assert outcome.selection is not None
    assert outcome.selection.selected_seed_ordinal == 1
    assert session.analysis_values.snapshot().revision == 0


def test_cancellation_and_interruption_stop_later_seeds() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start, start)),),
        objective_request_budget=5,
    )
    token = CancellationToken()
    cancellation_backend_calls = 0

    def cancel_after_trial(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal cancellation_backend_calls
        cancellation_backend_calls += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        if cancellation_backend_calls == 2:
            token.cancel()
        return _backend_result(candidate, residuals, status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", cancel_after_trial):
        cancelled = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert cancelled.terminal is GridDirectTrfTerminal.CANCELLED
    assert tuple(attempt.disposition for attempt in cancelled.attempts) == (
        GridSeedDisposition.ELIGIBLE,
        GridSeedDisposition.CANCELLED,
        GridSeedDisposition.NOT_STARTED,
    )
    assert cancelled.selection is None
    assert cancelled.accepted_result is None

    interruption_backend_calls = 0

    def interrupt_after_trial(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        nonlocal interruption_backend_calls
        interruption_backend_calls += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        if interruption_backend_calls == 2:
            raise KeyboardInterrupt
        return _backend_result(candidate, residuals, status=1, success=True)

    with (
        patch("chemex.optimize.direct_trf.least_squares", interrupt_after_trial),
        pytest.raises(GridDirectTrfInterrupted) as interrupted,
    ):
        execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert interrupted.value.outcome.terminal is GridDirectTrfTerminal.INTERRUPTED
    assert tuple(
        attempt.disposition for attempt in interrupted.value.outcome.attempts
    ) == (
        GridSeedDisposition.ELIGIBLE,
        GridSeedDisposition.INTERRUPTED,
        GridSeedDisposition.NOT_STARTED,
    )
    assert interrupted.value.outcome.selection is None
    assert interrupted.value.outcome.accepted_result is None
    assert session.analysis_values.snapshot().revision == 0


def test_table_and_two_dimensional_diagnostics_retain_exact_seed_provenance() -> None:
    _session, _experiments, parameterization, engine, problem = _qualification_problem(
        Method(fit=["PB"], fix=["KEX_AB"]),
    )
    pb_id, r1_id = problem.controlled_ids
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((r1_id, (problem.start[1],)), (pb_id, (problem.start[0],))),
        objective_request_budget=5,
    )

    def converge(fun: Callable[[Array], Array], x0: Array, **_kwargs: object) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", converge):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GridDirectTrfTerminal.ACCEPTED
    (row,) = outcome.table_rows
    assert row.seed_identity == invocation.seeds[0].identity
    assert row.seed_ordinal == 0
    assert row.axis_items == (
        (pb_id, problem.start[0]),
        (r1_id, problem.start[1]),
    )
    assert row.objective == outcome.attempts[0].objective
    assert row.disposition is GridSeedDisposition.ELIGIBLE
    assert row.selected
    (point_1d,) = outcome.diagnostic_points((pb_id,))
    (point_2d,) = outcome.diagnostic_points((r1_id, pb_id))
    assert point_1d.coordinates == ((pb_id, problem.start[0]),)
    assert point_2d.coordinates == (
        (r1_id, problem.start[1]),
        (pb_id, problem.start[0]),
    )
    assert point_1d.objective == point_2d.objective == row.objective
    assert point_1d.seed_identity == point_2d.seed_identity == row.seed_identity


def test_no_converged_candidate_has_no_selection_or_commit_authority() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start * 1.1)),),
        objective_request_budget=5,
    )

    def never_converge(
        fun: Callable[[Array], Array], x0: Array, **_kwargs: object
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=0, success=False)

    with patch("chemex.optimize.direct_trf.least_squares", never_converge):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GridDirectTrfTerminal.NO_ELIGIBLE_CANDIDATE
    assert all(
        attempt.disposition is GridSeedDisposition.NON_CONVERGED
        and attempt.candidate is None
        and attempt.objective is None
        for attempt in outcome.attempts
    )
    assert outcome.selection is None
    assert outcome.materialization is None
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot().revision == 0


def test_seed_local_materialization_failure_continues_to_later_seed() -> None:
    session, experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start)),),
        objective_request_budget=5,
    )
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate
    kernel_calls = 0

    def fail_first_local_materialization(*args: object, **kwargs: object) -> Array:
        nonlocal kernel_calls
        kernel_calls += 1
        if kernel_calls == 3:
            raise RuntimeError("qualified seed-local materialization failure")
        return original(*args, **kwargs)

    def converge(fun: Callable[[Array], Array], x0: Array, **_kwargs: object) -> object:
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with (
        patch.object(pulse_type, "calculate", fail_first_local_materialization),
        patch("chemex.optimize.direct_trf.least_squares", converge),
    ):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert kernel_calls == 7
    assert outcome.terminal is GridDirectTrfTerminal.ACCEPTED
    assert tuple(attempt.disposition for attempt in outcome.attempts) == (
        GridSeedDisposition.MATERIALIZATION_FAILURE,
        GridSeedDisposition.ELIGIBLE,
    )
    assert outcome.selection is not None
    assert outcome.selection.selected_seed_ordinal == 1
    assert outcome.accepted_result is not None
    assert session.analysis_values.snapshot().revision == 0


def test_selected_root_materialization_failure_never_falls_back() -> None:
    session, experiments, parameterization, engine, problem = _qualification_problem(
        read_methods([METHOD])["DEFAULT"]
    )
    (param_id,) = problem.controlled_ids
    start = problem.start[0]
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((param_id, (start, start)),),
        objective_request_budget=5,
    )
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate
    kernel_calls = 0

    def fail_seventh_evaluation(*args: object, **kwargs: object) -> Array:
        nonlocal kernel_calls
        kernel_calls += 1
        if kernel_calls == 7:
            raise RuntimeError("qualified final materialization failure")
        return original(*args, **kwargs)

    def converge(fun: Callable[[Array], Array], x0: Array, **_kwargs: object) -> object:
        candidate = np.asarray(problem.start, dtype=np.float64) * 0.9
        return _backend_result(candidate, fun(candidate), status=1, success=True)

    with (
        patch.object(pulse_type, "calculate", fail_seventh_evaluation),
        patch("chemex.optimize.direct_trf.least_squares", converge),
    ):
        outcome = execute_grid_direct_trf(
            problem,
            invocation,
            parameterization,
            engine,
        )

    assert kernel_calls == 7
    assert outcome.terminal is GridDirectTrfTerminal.MATERIALIZATION_FAILURE
    assert outcome.selection is not None
    assert outcome.selection.selected_seed_ordinal == 0
    assert outcome.materialization is not None
    assert outcome.accepted_result is None
    assert all(
        attempt.disposition is GridSeedDisposition.ELIGIBLE
        for attempt in outcome.attempts
    )
    assert session.analysis_values.snapshot().revision == 0
