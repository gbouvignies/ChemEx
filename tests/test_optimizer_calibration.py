"""Executable prospective-calibration contracts for issue #602."""

from __future__ import annotations

from pathlib import Path

import pytest

import tests.qualification.capture_optimizer_calibration as calibration
from tests.qualification.capture_optimizer_calibration import (
    SPECIFICATION,
    SPECIFICATION_ID,
    candidate_ordering_key,
    derive_de_root,
    holdout_decision,
    select_monotone_budget,
    validate_specification_commit,
)


def test_prospective_specification_freezes_every_acquisition_boundary() -> None:
    assert SPECIFICATION_ID == "chemex-optimizer-calibration-v1"
    assert SPECIFICATION == {
        "truth_cases": {
            "numerical_probes": (
                "trf-routine-quadratic-v1",
                "trf-difficult-rosenbrock-v1",
                "grid-27-seed-coverage-v1",
                "grid-candidate-ordering-v1",
                "de-bounded-search-v1",
            ),
            "routine_trf": "RELAXATION_HZNZ:G2N-HN",
            "difficult_trf": "DCEST_15N_3States:K19N-HN:3st_fork",
            "coupled_grid": "RELAXATION_HZNZ:G2N-HN",
            "decomposed_grid": "RELAXATION_HZNZ:all-five-profiles",
            "de_trf": "DCEST_15N_3States:K19N-HN:3st_fork",
        },
        "difficult_starts": {
            "dce_st": {
                "dw_ab": 2.51,
                "dw_ac": -1.97,
                "kex_ab": 10.0,
                "kex_ac": 10.0,
                "pb": 0.03,
                "pc": 0.03,
            }
        },
        "candidates": {
            "routine_trf_budgets": (8, 16, 32),
            "difficult_trf": (
                ("uniform", 2048, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
                ("physical", 256, (5.0, 5.0, 100.0, 100.0, 0.1, 0.1)),
                ("physical", 512, (5.0, 5.0, 100.0, 100.0, 0.1, 0.1)),
            ),
            "coupled_grid_budgets_per_seed": (16, 32, 64, 128),
            "decomposed_grid_budgets_per_component": (16, 32, 64, 128),
            "de": (
                ("compact-30x16", 5, 16, 510, 512),
                ("compact-48x32", 8, 32, 1584, 512),
            ),
        },
        "grid_coordinates": {
            "coupled_pb": (0.0, 0.05, 0.1),
            "decomposed_first_r1_factor": (0.75, 1.0, 1.25),
            "immutable_grid_seed_count": 27,
        },
        "de_coordinates": (
            ("dw_ab", -5.0, 5.0, "linear"),
            ("dw_ac", -5.0, 5.0, "linear"),
            ("kex_ab", 1.0, 1000.0, "log"),
            ("kex_ac", 1.0, 1000.0, "log"),
            ("pb", 0.001, 0.15, "linear"),
            ("pc", 0.001, 0.15, "linear"),
        ),
        "de_roots": {
            "calibration_count": 2,
            "replay_calibration_index": 0,
            "holdout_count": 1,
            "derivation": "u64be(sha256(ascii(chemex-602-v1|de|{phase}|{index:02d}))[0:8])",
        },
        "disqualifiers": (
            "typed_terminal_not_accepted",
            "non_finite_objective_or_vector",
            "objective_request_budget_exceeded",
            "counter_ordering_violation",
            "solver_backend_counter_missing",
            "grid_seed_or_component_incomplete",
            "candidate_ordering_mismatch",
            "transition_count_not_exactly_one",
            "fresh_root_materialization_count_not_exactly_one",
            "de_objective_not_below_90_percent_of_direct_reference",
            "replay_identity_mismatch",
            "holdout_failure",
            "resource_ceiling_exceeded",
        ),
        "selection": {
            "budget": "smallest-passing-with-all-larger-candidates-passing",
            "difficult_trf": "objective-requests,budget,policy-ordinal",
            "de": "maximum-objective-ratio,total-objective-requests,topology-ordinal",
            "tie_break": "frozen-candidate-declaration-order",
            "edge_only": "grid_saturated_unsupported",
            "runner_up_after_holdout": "forbidden",
        },
        "accounting": {
            "authoritative_unit": "chemex-objective-request",
            "nonfungible": (
                "solver_requests_received",
                "objective_requests_accepted",
                "objective_evaluations_completed",
                "backend_nfev",
                "backend_njev",
                "grid_seed_allocations",
                "grouped_component_allocations",
                "de_population_size",
                "de_generation_limit",
                "de_stage_allocation",
                "trf_polish_allocation",
                "search_to_polish_transfers",
                "root_materializations",
                "cache_hits",
                "cache_misses",
                "elapsed_seconds_diagnostic_only",
            ),
            "retries": 0,
            "borrowing": False,
            "redistribution": False,
        },
        "resource_limits": {
            "routine_trf_objective_requests": 88,
            "difficult_trf_objective_requests": 3328,
            "coupled_grid_objective_requests": 1104,
            "decomposed_grid_objective_requests": 5520,
            "de_trf_objective_requests": 10428,
            "wall_clock_scientific_budget": None,
            "workers": 1,
            "native_threads": 1,
        },
        "holdouts": {
            "routine_trf": "RELAXATION_HZNZ:Y3N-HN",
            "difficult_trf": "DCEST_15N_3States:D20N-HN:3st_fork",
            "coupled_grid": "RELAXATION_HZNZ:Y3N-HN",
            "decomposed_grid": "same-case-selected-policy-exact-replay",
            "de_trf": "one-unopened-derived-root-on-selected-K19N-HN-policy",
        },
        "unsupported": {
            "de_beyond_selected_six_coordinates": True,
            "adaptive_grid_refinement": True,
            "adaptive_de_topology_or_window_extension": True,
            "wall_clock_budgets": True,
            "hidden_retries": True,
        },
    }


def test_root_derivation_is_deterministic_and_phase_separated() -> None:
    calibration = tuple(derive_de_root("calibration", index) for index in range(2))
    holdout = derive_de_root("holdout", 0)

    assert calibration == (
        14007474848671902346,
        10563064000758873076,
    )
    assert holdout == 3611316752655662115
    assert len({*calibration, holdout}) == 3


@pytest.mark.parametrize(
    ("passing", "expected"),
    [
        (set(), ("NO_ADEQUATE_CANDIDATE", None)),
        ({128}, ("GRID_SATURATED", None)),
        ({64, 128}, ("SELECTED", 64)),
        ({32, 64, 128}, ("SELECTED", 32)),
        ({16, 32, 64, 128}, ("SELECTED", 16)),
        ({16, 64, 128}, ("NON_MONOTONE_ADEQUACY", None)),
    ],
)
def test_budget_selection_is_monotone_and_edge_fails_closed(
    passing: set[int],
    expected: tuple[str, int | None],
) -> None:
    assert select_monotone_budget((16, 32, 64, 128), passing) == expected


def test_candidate_ordering_has_one_frozen_total_order() -> None:
    records = (
        {"objective": 1.0, "vector": (1.0, 0.0), "ordinal": 0},
        {"objective": 0.0, "vector": (0.0, 0.0), "ordinal": 1},
        {"objective": 1.0, "vector": (0.0, 1.0), "ordinal": 2},
        {"objective": 1.0, "vector": (-1.0, 0.0), "ordinal": 3},
    )

    assert tuple(
        item["ordinal"] for item in sorted(records, key=candidate_ordering_key)
    ) == (1, 3, 2, 0)


def test_holdout_failure_invalidates_scope_without_runner_up() -> None:
    assert holdout_decision("physical-256", (True,)) == (
        "SELECTED",
        "physical-256",
    )
    assert holdout_decision("physical-256", (False,)) == (
        "HOLDOUT_FAILED_UNSUPPORTED",
        None,
    )


def test_authoritative_acquisition_requires_the_frozen_head_commit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "checkout"
    (root / ".git").mkdir(parents=True)
    (root / ".git/HEAD").write_text("a" * 40 + "\n", encoding="ascii")
    monkeypatch.setattr(calibration, "ROOT", root)

    assert validate_specification_commit("a" * 40) == "a" * 40
    with pytest.raises(
        RuntimeError, match="frozen specification commit.*detached HEAD"
    ):
        validate_specification_commit("b" * 40)

    (root / ".git/HEAD").write_text(
        "ref: refs/heads/codex/602-calibrate-optimizer-policies\n",
        encoding="ascii",
    )
    with pytest.raises(RuntimeError, match="detached HEAD"):
        validate_specification_commit("a" * 40)
