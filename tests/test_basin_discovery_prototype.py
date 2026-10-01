"""Small checks for the non-production #712 research runner."""

from types import SimpleNamespace

import pytest
from benchmarks.basin_discovery import (
    branch_vectors,
    build_context,
    measured,
    outer_points,
    refine,
    selected_indices,
)

from chemex.optimize.profiled_grid import execute_profiled_grid


def test_sobol_design_is_reproducible_and_preserves_separate_supplied_start() -> None:
    coordinates = (("k", 10.0, 1000.0, "log"), ("p", 0.01, 0.2, "linear"))
    start = {"k": 80.0, "p": 0.05}
    points = outer_points(coordinates, start, 3, 597)
    assert points == outer_points(coordinates, start, 3, 597)
    assert points != outer_points(coordinates, start, 3, 598)
    assert points[0] == (80.0, 0.05)
    assert len(points) == 9
    assert all(10.0 <= k <= 1000.0 and 0.01 <= p <= 0.2 for k, p in points)


def test_twenty_independent_sign_ambiguities_do_not_form_a_root_product() -> None:
    ids = ("k", *(f"__DW_AB_{index}N" for index in range(20)))
    problem = SimpleNamespace(controlled_ids=ids, start=(100.0, *((2.0,) * 20)))
    alternatives = []
    for param_id in ids[1:]:
        factor = SimpleNamespace(nuisance_ids=(param_id,))
        normal = SimpleNamespace(nuisance_items=((param_id, 2.2),))
        (vector,) = tuple(branch_vectors(problem, factor, normal))
        changed = tuple(
            key for key, a, b in zip(ids, problem.start, vector, strict=True) if a != b
        )
        assert changed == (param_id,)
        alternatives.append(vector)
    assert len(alternatives) == 20  # 20 normal + 20 alternate TRFs, not 2**20.


def test_zero_start_uses_endpoint_magnitude_and_multidw_stays_inside_factor() -> None:
    problem = SimpleNamespace(
        controlled_ids=("k", "__DW_AB_1N", "__DW_AC_1N", "__DW_AB_2N"),
        start=(100.0, 0.0, -3.0, 4.0),
    )
    factor = SimpleNamespace(nuisance_ids=("__DW_AB_1N", "__DW_AC_1N"))
    normal = SimpleNamespace(nuisance_items=(("__DW_AB_1N", 1.2), ("__DW_AC_1N", -3.4)))
    vectors = tuple(branch_vectors(problem, factor, normal))
    assert len(vectors) == 4
    assert {vector[1:3] for vector in vectors} == {
        (-1.2, -3.0),
        (-1.2, 3.0),
        (1.2, -3.0),
        (1.2, 3.0),
    }
    assert all(vector[0] == 100.0 and vector[3] == 4.0 for vector in vectors)
    zero = SimpleNamespace(controlled_ids=("__DW_AB_1N",), start=(0.0,))
    normal = SimpleNamespace(nuisance_items=(("__DW_AB_1N", 0.0),))
    assert not tuple(
        branch_vectors(zero, SimpleNamespace(nuisance_ids=zero.controlled_ids), normal)
    )


def test_diversity_uses_normalized_log_coordinates_and_excludes_failures() -> None:
    coordinates = (("k", 10.0, 1000.0, "log"),)
    records = [
        {"vector": (1.0,), "chi2": 1.0},
        {"vector": (2.0,), "chi2": 2.0},
        {"vector": (3.0,), "chi2": 3.0},
        {"vector": None, "chi2": None},
    ]
    points = [(10.0,), (11.0,), (1000.0,), (100.0,)]
    assert selected_indices(records, points, coordinates, 2, 0.25) == [0, 2]


def test_real_cpmg_factor_branches_reconstruct_and_refine_the_lower_basin() -> None:
    context = build_context("cpmg")
    root = context.problem
    axes = {
        key: (dict(root.independent_items)[key],) for key, *_ in context.coordinates
    }
    outcome, trace = measured(
        context,
        lambda: execute_profiled_grid(
            root,
            axes,
            context.parameterization,
            context.engine,
            objective_request_budget=2000,
        ),
        branches=True,
        endpoint_mirror=True,
    )
    assert outcome.accepted_result is not None
    # Basin-scale tolerances accommodate finite-difference and host variation;
    # they remain much smaller than the 4.9 chi-square sign-branch difference.
    assert outcome.accepted_result.chi_square == pytest.approx(460.1196200, abs=1e-3)
    assert trace["local_refinements"] == 10
    assert trace["model_calculations"] > 0
    assert all(len(call["ids"]) == 3 for call in trace["calls"])
    final = refine(context, outcome.accepted_result.vector)
    assert final["terminal"] == "accepted"
    assert final["chi2"] == pytest.approx(429.6700333, abs=1e-3)
