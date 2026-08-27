"""Scientific qualification for profiled chi-square GRID (#702)."""

from __future__ import annotations

from itertools import product
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.method_execution import normalize_methods_for_execution
from chemex.configuration.method_validation import resolve_grid_axes
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    CancellationToken,
    DirectTrfConstructionError,
    FitCommitFailureCategory,
    FitCommitTerminal,
    OptimizationProblem,
    execute_fit_commit,
)
from chemex.optimize.grouped_direct_trf import _profile_dependencies
from chemex.optimize.native_deterministic import _write_grid_output
from chemex.optimize.profiled_grid import (
    ProfiledGridFactor,
    ProfiledGridFactorResult,
    ProfiledGridPoint,
    ProfiledGridPointStatus,
    ProfiledGridTerminal,
    aggregate_profiled_grids,
    discover_profiled_grid_factors,
    execute_profiled_grid,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
CPMG = ROOT / "examples/Experiments/CPMG_15N_IP"
RELAXATION = ROOT / "examples/Experiments/RELAXATION_HZNZ"


def _synthetic_factor_results() -> tuple[ProfiledGridFactorResult, ...]:
    values = {"k": (0.0, 1.0), "d1": (0.0, 1.0), "d2": (0.0, 1.0)}

    def result(local: str, nuisance: str, ordinal: int) -> ProfiledGridFactorResult:
        factor = ProfiledGridFactor(
            ordinal=ordinal,
            profile_indices=(ordinal,),
            grid_ids=("k", local),
            nuisance_ids=(nuisance,),
        )
        points: list[ProfiledGridPoint] = []
        for point_ordinal, (k, d) in enumerate(product(values["k"], values[local])):
            # Two coherent joint minima exist: (k=0,d1=d2=1) and
            # (k=1,d1=d2=0). Independent 1D first-argmins form (0,0,0),
            # which is not a joint minimum.
            base = 0.0 if d == 1.0 - k else 5.0
            nuisance_value = k + 2.0 * d if nuisance == "n1" else 2.0 * k + d
            points.append(
                ProfiledGridPoint(
                    ordinal=point_ordinal,
                    axis_items=(("k", k), (local, d)),
                    status=ProfiledGridPointStatus.SUCCESS,
                    chi_square=base,
                    nuisance_items=((nuisance, nuisance_value),),
                    objective_evaluations=3,
                )
            )
        return ProfiledGridFactorResult(factor, tuple(points))

    return result("d1", "n1", 0), result("d2", "n2", 1)


def test_factorized_profiles_match_brute_force_and_select_one_coherent_truth() -> None:
    axes = {"k": (0.0, 1.0), "d1": (0.0, 1.0), "d2": (0.0, 1.0)}
    factors = _synthetic_factor_results()

    aggregate = aggregate_profiled_grids(axes, factors)

    brute = np.empty((2, 2, 2), dtype=np.float64)
    for indices, (k, d1, d2) in zip(
        product(range(2), repeat=3),
        product(axes["k"], axes["d1"], axes["d2"]),
        strict=True,
    ):
        brute[indices] = (0.0 if d1 == 1.0 - k else 5.0) + (
            0.0 if d2 == 1.0 - k else 5.0
        )

    assert aggregate.selection.grid_items == (
        ("k", 0.0),
        ("d1", 1.0),
        ("d2", 1.0),
    )
    assert aggregate.selection.nuisance_items == (("n1", 2.0), ("n2", 1.0))
    assert aggregate.selection.chi_square == 0.0
    assert brute[0, 0, 0] == 10.0  # independent marginal first-argmins are incoherent

    for axis_id, profile in aggregate.profiles_1d.items():
        axis = tuple(axes).index(axis_id)
        expected = np.min(brute, axis=tuple(i for i in range(3) if i != axis))
        np.testing.assert_allclose(profile.chi_square, expected, rtol=0.0, atol=0.0)
    for axis_ids, profile in aggregate.profiles_2d.items():
        kept = tuple(tuple(axes).index(axis_id) for axis_id in axis_ids)
        reduced = tuple(i for i in range(3) if i not in kept)
        expected = np.min(brute, axis=reduced)
        if kept != tuple(sorted(kept)):
            expected = expected.T
        np.testing.assert_allclose(profile.chi_square, expected, rtol=0.0, atol=0.0)


def test_unusable_factor_point_is_not_confused_with_high_finite_chi_square() -> None:
    factor = ProfiledGridFactor(0, (0,), ("g",), ())
    result = ProfiledGridFactorResult(
        factor,
        (
            ProfiledGridPoint(
                0,
                (("g", 0.0),),
                ProfiledGridPointStatus.FAILED,
                failure="synthetic failure",
            ),
            ProfiledGridPoint(
                1,
                (("g", 1.0),),
                ProfiledGridPointStatus.SUCCESS,
                chi_square=1.0e300,
                objective_evaluations=1,
            ),
        ),
    )

    aggregate = aggregate_profiled_grids({"g": (0.0, 1.0)}, (result,))

    assert aggregate.selection.grid_items == (("g", 1.0),)
    assert np.isinf(aggregate.profiles_1d["g"].chi_square[0])
    assert aggregate.profiles_1d["g"].chi_square[1] == 1.0e300


def test_raw_factor_output_distinguishes_failure_from_valid_high_value(
    tmp_path: Path,
) -> None:
    session, _parameterization, _engine, problem = _relaxation_problem()
    (axis_id,) = problem.controlled_ids
    factor = ProfiledGridFactor(0, (0,), (axis_id,), ())
    result = ProfiledGridFactorResult(
        factor,
        (
            ProfiledGridPoint(
                0,
                ((axis_id, 1.0),),
                ProfiledGridPointStatus.FAILED,
                failure="synthetic failure",
            ),
            ProfiledGridPoint(
                1,
                ((axis_id, 3.0),),
                ProfiledGridPointStatus.SUCCESS,
                chi_square=1.0e300,
                objective_evaluations=1,
            ),
        ),
    )
    aggregate = aggregate_profiled_grids({axis_id: (1.0, 3.0)}, (result,))
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    outcome = SimpleNamespace(
        aggregate=aggregate,
        accepted_result=SimpleNamespace(chi_square=1.0e300),
        factors=(result,),
    )

    with (
        patch("chemex.optimize.native_deterministic.plot_grid_1d"),
        patch("chemex.optimize.native_deterministic.plot_grid_2d"),
    ):
        _write_grid_output(  # ty: ignore[invalid-argument-type]
            outcome,
            tmp_path,
            parameter_model=parameter_model,
            accepted_values={axis_id: 3.0},
        )

    (factor_output,) = tuple((tmp_path / "Grid" / "Factors").glob("*.tsv"))
    rows = factor_output.read_text(encoding="utf-8").splitlines()
    assert "\t\tfailed\tfalse\t0\tsynthetic failure" in rows[1]
    assert "\t1e+300\tsuccess\ttrue\t1\t" in rows[2]


def test_factor_discovery_removes_shared_grid_coupling_but_keeps_local_ownership() -> (
    None
):
    factors = discover_profiled_grid_factors(
        (
            frozenset(("k", "p", "d1", "n1")),
            frozenset(("k", "p", "d1", "n1")),
            frozenset(("k", "p", "d2", "n2")),
            frozenset(("k", "p", "d2", "n2")),
        ),
        (True, True, True, True),
        grid_ids=("k", "p", "d1", "d2"),
        controlled_ids=("k", "p", "d1", "d2", "n1", "n2"),
    )

    assert tuple(factor.profile_indices for factor in factors) == ((0, 1), (2, 3))
    assert tuple(factor.grid_ids for factor in factors) == (
        ("k", "p", "d1"),
        ("k", "p", "d2"),
    )
    assert tuple(factor.nuisance_ids for factor in factors) == (("n1",), ("n2",))
    axis_sizes = {"k": 10, "p": 10, "d1": 5, "d2": 5}
    assert (
        sum(
            np.prod([axis_sizes[axis] for axis in factor.grid_ids])
            for factor in factors
        )
        == 1000
    )
    assert np.prod([axis_sizes[axis] for axis in ("k", "p", "d1", "d2")]) == 2500


def test_shipped_cpmg_step1_discovers_five_500_point_residue_factors() -> None:
    session = AnalysisSession.create()
    experiments = build_experiments(
        sorted((CPMG / "Experiments").glob("*.toml")),
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([CPMG / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    plan, operational = normalize_methods_for_execution(
        read_method_plan([CPMG / "Methods/method_grid.toml"])
    )
    step = plan.steps[0]
    experiments.select(operational[step.name].selection)
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()[step.name], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None and parameter_model is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    assert step.search is not None
    axes = resolve_grid_axes(
        step.search,  # ty: ignore[invalid-argument-type]
        parameter_model,
        active_scope_ids=parameterization.scope_ids,
        final_fit_ids=problem.controlled_ids,
    )
    factors = discover_profiled_grid_factors(
        _profile_dependencies(problem, parameterization, engine),
        tuple(
            bool(profile.retained_observation_indices)
            for profile in engine.plan.profiles
        ),
        grid_ids=tuple(axis.param_id for axis in axes),
        controlled_ids=problem.controlled_ids,
    )

    assert len(axes) == 7
    assert len(factors) == 5
    axis_sizes = {axis.param_id: len(axis.values) for axis in axes}
    assert [
        tuple(
            parameter_model.definitions[param_id].name for param_id in factor.grid_ids
        )
        for factor in factors
    ] == [("KEX_AB", "PB", "DW_AB")] * 5
    assert [
        int(np.prod([axis_sizes[param_id] for param_id in factor.grid_ids]))
        for factor in factors
    ] == [500] * 5
    assert (
        sum(
            int(np.prod([axis_sizes[param_id] for param_id in factor.grid_ids]))
            for factor in factors
        )
        == 2500
    )
    assert int(np.prod(tuple(axis_sizes.values()))) == 312_500


def test_shipped_cpmg_step3_discovers_residue_local_20_point_factors() -> None:
    session = AnalysisSession.create()
    experiments = build_experiments(
        sorted((CPMG / "Experiments").glob("*.toml")),
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([CPMG / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    plan, operational = normalize_methods_for_execution(
        read_method_plan([CPMG / "Methods/method_grid.toml"])
    )
    step = plan.steps[2]
    experiments.select(operational[step.name].selection)
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()[step.name], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None and parameter_model is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    assert step.search is not None
    axes = resolve_grid_axes(
        step.search,  # ty: ignore[invalid-argument-type]
        parameter_model,
        active_scope_ids=parameterization.scope_ids,
        final_fit_ids=problem.controlled_ids,
    )
    factors = discover_profiled_grid_factors(
        _profile_dependencies(problem, parameterization, engine),
        tuple(
            bool(profile.retained_observation_indices)
            for profile in engine.plan.profiles
        ),
        grid_ids=tuple(axis.param_id for axis in axes),
        controlled_ids=problem.controlled_ids,
    )

    assert len(axes) == 54
    assert len(factors) == 54
    assert all(len(factor.grid_ids) == 1 for factor in factors)
    assert all(len(factor.nuisance_ids) == 2 for factor in factors)
    assert all(len(axes[index].values) == 20 for index in range(len(axes)))
    assert sum(len(axes[index].values) for index in range(len(factors))) == 1_080
    assert 20**54 > 10**70


def _relaxation_problem():
    session = AnalysisSession.create()
    experiments = build_experiments(
        [RELAXATION / "Experiments/800mhz.toml"],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([RELAXATION / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    parameterization = session.compile_parameterization(
        Method(
            fit=["R1A_A"],
            fix=["PB", "KEX_AB", "ETAZ_A", "R1_A"],
        ),
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
    return session, parameterization, engine, problem


def test_evaluation_only_grid_point_never_releases_axis_and_commits_once() -> None:
    session, parameterization, engine, problem = _relaxation_problem()
    (axis_id,) = problem.controlled_ids
    before = session.analysis_values.snapshot()

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        side_effect=AssertionError("evaluation-only GRID must not run TRF"),
    ):
        outcome = execute_profiled_grid(
            problem,
            {axis_id: (1.0, 3.0)},
            parameterization,
            engine,
            objective_request_budget=10,
            x_scale=(1.0,),
        )

    assert outcome.terminal is ProfiledGridTerminal.ACCEPTED
    assert session.analysis_values.snapshot() == before
    assert outcome.aggregate is not None
    assert outcome.aggregate.selection.grid_items == ((axis_id, 3.0),)
    assert outcome.aggregate.selection.nuisance_items == ()
    assert outcome.accepted_result is not None
    assert outcome.accepted_result.vector == (3.0,)
    assert outcome.commit_authority is not None

    operation = execute_fit_commit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )
    assert operation.terminal is FitCommitTerminal.COMMITTED
    assert session.analysis_values.snapshot().revision == before.revision + 1
    with pytest.raises(DirectTrfConstructionError, match="live fit commit authority"):
        execute_fit_commit(
            outcome.accepted_result,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )


def test_profiled_grid_stale_or_cancelled_execution_cannot_commit() -> None:
    session, parameterization, engine, problem = _relaxation_problem()
    (axis_id,) = problem.controlled_ids
    outcome = execute_profiled_grid(
        problem,
        {axis_id: (1.0, 3.0)},
        parameterization,
        engine,
        objective_request_budget=10,
        x_scale=(1.0,),
    )
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    before = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {axis_id: before[axis_id]}, expected=before, scope=(axis_id,)
    )
    stale = execute_fit_commit(
        outcome.accepted_result,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )
    assert stale.terminal is FitCommitTerminal.FAILED
    assert stale.failure is not None
    assert stale.failure.category is FitCommitFailureCategory.STALE_REVISION

    fresh_session, fresh_parameterization, fresh_engine, fresh_problem = (
        _relaxation_problem()
    )
    token = CancellationToken()
    token.cancel()
    cancelled = execute_profiled_grid(
        fresh_problem,
        {fresh_problem.controlled_ids[0]: (1.0, 3.0)},
        fresh_parameterization,
        fresh_engine,
        objective_request_budget=10,
        x_scale=(1.0,),
        cancellation=token,
    )
    assert cancelled.terminal is ProfiledGridTerminal.CANCELLED
    assert cancelled.accepted_result is None
    assert cancelled.commit_authority is None
    assert fresh_session.analysis_values.snapshot().revision == 0


def test_cpmg_grid_points_hold_axes_while_optimizing_factor_nuisance() -> None:
    session = AnalysisSession.create()
    experiments = build_experiments(
        sorted((CPMG / "Experiments").glob("*.toml")),
        Selection(include=[SpinSystem.from_name("15")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([CPMG / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    plan = read_method_plan([CPMG / "Methods/method_grid.toml"])
    step = plan.steps[0]
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()[step.name], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert configuration is not None and parameter_model is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    assert step.search is not None
    resolved = resolve_grid_axes(
        step.search,  # ty: ignore[invalid-argument-type]
        parameter_model,
        active_scope_ids=parameterization.scope_ids,
        final_fit_ids=problem.controlled_ids,
    )
    axes = {axis.param_id: (axis.values[0],) for axis in resolved}
    children: list[OptimizationProblem] = []
    derive = OptimizationProblem.derive_profiled_grid_point

    def record_child(
        self: OptimizationProblem, **kwargs: object
    ) -> OptimizationProblem:
        child = derive(self, **kwargs)  # ty: ignore[invalid-argument-type]
        children.append(child)
        return child

    with patch.object(
        OptimizationProblem,
        "derive_profiled_grid_point",
        record_child,
    ):
        outcome = execute_profiled_grid(
            problem,
            axes,
            parameterization,
            engine,
            objective_request_budget=12_000,
            x_scale=tuple(max(1.0, abs(value)) for value in problem.start),
        )

    assert outcome.terminal is ProfiledGridTerminal.ACCEPTED
    assert children
    grid_ids = frozenset(axes)
    for child in children:
        assert not grid_ids.intersection(child.controlled_ids)
        held = dict(child.held_items)
        assert {param_id: held[param_id] for param_id in axes} == {
            param_id: values[0] for param_id, values in axes.items()
        }
        assert child.controlled_ids
    assert outcome.aggregate is not None
    assert dict(outcome.aggregate.selection.grid_items) == {
        param_id: values[0] for param_id, values in axes.items()
    }
    assert outcome.aggregate.selection.nuisance_items
