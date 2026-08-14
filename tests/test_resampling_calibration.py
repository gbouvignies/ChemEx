"""Prospective checks for the frozen issue #603 calibration specification."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any, cast

import pytest

from chemex.optimize.native_resampling import OptimizationStrategy
from tests.qualification import capture_resampling_calibration as calibration


def test_frozen_roots_match_derivation_and_are_globally_disjoint() -> None:
    roots = []
    for family, phases in calibration.ROOTS.items():
        for phase, literal_roots in phases.items():
            assert len(literal_roots) == (16 if phase == "calibration" else 8)
            assert literal_roots == tuple(
                calibration.derive_root(family, phase, index)
                for index in range(len(literal_roots))
            )
            roots.extend(literal_roots)
    assert len(roots) == 72
    assert 0 not in roots
    assert len(set(roots)) == len(roots)


def test_independently_computed_truth_matches_frozen_references() -> None:
    computed = calibration.truth_estimands()
    assert computed.keys() == calibration.REFERENCE_TRUTH.keys()
    for family, parameter_truth in computed.items():
        for parameter_id, values in parameter_truth.items():
            assert values == pytest.approx(
                calibration.REFERENCE_TRUTH[family][parameter_id],
                rel=0.0,
                abs=1.0e-15,
            )


def test_fixture_and_policy_are_exactly_frozen_without_running_replicates() -> None:
    assert calibration.CANDIDATE_COUNTS == (64, 128, 256)
    assert (calibration.OUTPUT_SCOPE, calibration.THETA) == (("A", "B"), (1.0, 2.0))
    assert calibration.X == (-2, -1, 0, 1, 2, 3)
    assert calibration.CALCULATED == (-3, -1, 1, 3, 5, 7)
    assert calibration.RESIDUAL == (1, -2, 1, -1, 2, -1)
    assert calibration.OBSERVED == (-2, -3, 2, 2, 7, 6)
    assert calibration.ERRORS == (1, 1, 1, 1, 1, 1)
    assert calibration.PROFILE_INDICES == ((0, 1, 2), (3, 4, 5))
    assert calibration.REFERENCES == (True, False, False, True, False, False)
    assert calibration.NUCLEUS_GROUPS == ("N1", "N2", "N3", "N1", "N2", "N3")
    assert calibration.SUMMARY_POLICY.percentile_method == "median_unbiased"
    assert calibration.SUMMARY_POLICY.identity == calibration.SUMMARY_POLICY_IDENTITY
    assert calibration.THRESHOLDS == {
        "bias": 0.50,
        "spread": 0.35,
        "q2.5": 0.90,
        "median": 0.50,
        "q97.5": 0.90,
        "failure_prevalence": 0.0,
    }
    fixture = calibration.native_fixture()
    assert fixture.dataset.observed == tuple(map(float, calibration.OBSERVED))
    assert fixture.dataset.calculated == tuple(map(float, calibration.CALCULATED))
    assert fixture.dataset.standard_errors == tuple(map(float, calibration.ERRORS))
    assert fixture.dataset.mask == (True,) * 6
    for family, scheme in calibration.FAMILY_SCHEMES.items():
        plan = calibration.make_plan(
            fixture,
            family,
            64,
            calibration.ROOTS[family]["calibration"][0],
        )
        assert plan.scheme is scheme
        assert plan.strategy is OptimizationStrategy.DIRECT_TRF
        assert plan.strategy_settings == ()
        assert plan.output_scope == calibration.OUTPUT_SCOPE
        assert plan.minimum_successful_count == plan.replicate_count == 64


@pytest.mark.parametrize(
    ("passing", "expected"),
    (
        (set(), ("NO_ADEQUATE_CANDIDATE", None)),
        ({64}, ("NON_MONOTONE_ADEQUACY", None)),
        ({128}, ("NON_MONOTONE_ADEQUACY", None)),
        ({256}, ("GRID_SATURATED", None)),
        ({64, 128}, ("NON_MONOTONE_ADEQUACY", None)),
        ({64, 256}, ("NON_MONOTONE_ADEQUACY", None)),
        ({128, 256}, ("SELECTED", 128)),
        ({64, 128, 256}, ("SELECTED", 64)),
    ),
)
def test_selection_table_is_exact(
    passing: set[int], expected: tuple[str, int | None]
) -> None:
    assert calibration.select_candidate(passing) == expected


def test_holdout_is_strict_and_cannot_fall_back_to_an_alternate() -> None:
    assert calibration.holdout_decision(64, (True,) * 8) == ("SELECTED", 64)
    result = calibration.holdout_decision(64, (True,) * 7 + (False,))
    assert result == ("HOLDOUT_FAILED_UNSUPPORTED", None)


def test_replay_signature_contains_every_exact_scientific_field() -> None:
    value = lambda text: SimpleNamespace(value=text)  # noqa: E731
    outcome = SimpleNamespace(
        ordinal=0,
        request_identity="request",
        seed=17,
        draw_identity="draw",
        stage=value("terminal"),
        disposition=value("succeeded"),
        success=SimpleNamespace(identity="success"),
        identity="outcome",
    )
    evidence = SimpleNamespace(
        plan=SimpleNamespace(identity="plan"),
        identity="evidence",
        population_identity="population",
        lifecycle=value("completed"),
        completed_count=1,
        successful_count=1,
        failed_count=0,
        outcomes=(outcome,),
    )
    operation = SimpleNamespace(
        terminal=value("completed"), evidence=evidence, unstarted_ordinals=()
    )
    summary = SimpleNamespace(to_record=lambda: {"policy_identity": "policy"})
    summary_outcome = SimpleNamespace(terminal=value("completed"), summary=summary)
    pair = cast("Any", (operation, summary_outcome))
    assert calibration.replay_equal(pair, pair)
    assert (signature := calibration.replay_signature(*pair)) == {
        "operation_terminal": "completed",
        "unstarted_ordinals": (),
        "plan_identity": "plan",
        "evidence_identity": "evidence",
        "population_identity": "population",
        "lifecycle": "completed",
        "counts": (1, 1, 0),
        "outcomes": (
            (0, "request", 17, "draw", "terminal", "succeeded", "success", "outcome"),
        ),
        "summary_terminal": "completed",
        "summary_record": {"policy_identity": "policy"},
    }
    evidence.identity = "foreign-evidence"
    assert calibration.replay_signature(*pair) != signature
