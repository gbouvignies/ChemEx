"""Behavioral qualification for the non-authoritative selected-coordinate DE seam."""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import patch

import pytest

from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.de_direct_trf import (
    DeCoordinateSemantics,
    DeSearchCoordinate,
    DeSearchInvocation,
    DeSearchTerminal,
    execute_de_search,
)
from chemex.optimize.direct_trf import (
    DeSearchProblemDerivation,
    DirectTrfConstructionError,
    OptimizationProblem,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"


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
    assert session.try_build_analysis_values()
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


def test_product_search_plan_is_bounded_seeded_and_non_authoritative() -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    pb_id, held_id = problem.controlled_ids
    pb_index = problem.controlled_ids.index(pb_id)

    invocation = DeSearchInvocation.for_product_problem(
        problem,
        search_coordinates=(
            (
                pb_id,
                problem.lower_bounds[pb_index],
                problem.upper_bounds[pb_index],
                DeCoordinateSemantics.LINEAR,
            ),
        ),
        root_seed=597,
    )

    assert isinstance(invocation.search_problem.derivation, DeSearchProblemDerivation)
    assert invocation.search_problem.controlled_ids == (pb_id,)
    assert dict(invocation.search_problem.held_items)[held_id] == problem.start[1]
    assert not invocation.search_problem.acceptance_authority
    assert invocation.root_seed == 597
    assert invocation.population.dimension == 1
    assert invocation.population.multiplier == 15
    assert invocation.population.size == 15
    assert invocation.population.maximum_generations == 1000
    assert invocation.de_objective_request_budget == 15 * 1001


@pytest.mark.parametrize("root_seed", (-1, 2**64, True, None))
def test_de_requires_an_explicit_unsigned_64_bit_root_seed(root_seed: object) -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    param_id = problem.controlled_ids[0]
    index = problem.controlled_ids.index(param_id)

    with pytest.raises(DirectTrfConstructionError, match="unsigned 64-bit root seed"):
        DeSearchInvocation.for_product_problem(
            problem,
            search_coordinates=(
                (
                    param_id,
                    problem.lower_bounds[index],
                    problem.upper_bounds[index],
                    DeCoordinateSemantics.LINEAR,
                ),
            ),
            root_seed=root_seed,
        )


def test_log_coordinate_has_exact_endpoints_and_stable_round_trip() -> None:
    coordinate = DeSearchCoordinate("KEX_AB", 100.0, 5000.0, "log", 0)

    assert coordinate.to_physical(coordinate.solver_bounds[0]) == 100.0
    assert coordinate.to_physical(coordinate.solver_bounds[1]) == 5000.0
    assert coordinate.to_physical(coordinate.to_solver(597.0)) == pytest.approx(597.0)


def test_captured_value_outside_search_range_is_valid_and_not_clipped() -> None:
    _session, _experiments, _parameterization, _engine, problem = (
        _qualification_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    selected_id = problem.controlled_ids[1]
    selected_index = problem.controlled_ids.index(selected_id)
    captured = problem.start[selected_index]
    assert captured > 4.0

    invocation = DeSearchInvocation.for_product_problem(
        problem,
        search_coordinates=((selected_id, 1.0, 4.0, "log"),),
        root_seed=597,
    )

    assert invocation.search_problem.start == (captured,)
    assert invocation.search_coordinates[0].solver_initial == pytest.approx(
        math.log(2.0)
    )


def test_search_interruption_returns_no_candidate_or_fit_authority() -> None:
    session, _experiments, parameterization, engine, problem = _qualification_problem(
        Method(fit=["PB"], fix=["KEX_AB"])
    )
    selected_id = problem.controlled_ids[0]
    index = problem.controlled_ids.index(selected_id)
    invocation = DeSearchInvocation.for_product_problem(
        problem,
        search_coordinates=(
            (
                selected_id,
                problem.lower_bounds[index],
                problem.upper_bounds[index],
                "linear",
            ),
        ),
        root_seed=597,
    )
    before = session.analysis_values.snapshot()

    with patch(
        "chemex.optimize.de_direct_trf._invoke_de_backend",
        side_effect=KeyboardInterrupt,
    ):
        outcome = execute_de_search(problem, invocation, parameterization, engine)

    assert outcome.terminal is DeSearchTerminal.INTERRUPTED
    assert outcome.best_candidate is None
    assert not outcome.restart_eligible
    assert session.analysis_values.snapshot() == before
