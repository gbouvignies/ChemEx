"""Prospective checks for the frozen issue #603 calibration specification."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from chemex.migration_core import migration_core_authority_selection
from chemex.optimize.native_resampling import OptimizationStrategy
from tests.qualification import capture_resampling_calibration as calibration

ROOT = Path(__file__).parents[1]
EVIDENCE = ROOT / "tests/fixtures/canonical_resampling_calibration_v1.json"


def _canonical_hash(value: object) -> str:
    encoded = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


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


def test_canonical_evidence_reconstructs_frozen_decisions() -> None:
    record = json.loads(EVIDENCE.read_text(encoding="ascii"))
    identity = record.pop("identity")
    assert identity == _canonical_hash(record)
    assert record["source"] == {
        "dependency_lock_sha256": hashlib.sha256(
            (ROOT / "uv.lock").read_bytes()
        ).hexdigest(),
        "qualification_script_sha256": hashlib.sha256(
            Path(calibration.__file__).read_bytes()
        ).hexdigest(),
        "raw_acquisition_sha256": "1ac0c303f1ac0f077848f7b75b5c12987c8ddec40c9fc1c3f855e34813692c71",
        "specification_commit": "61b6328641bc31d1b0d53e4769c39a616658f24c",
    }
    authority = migration_core_authority_selection()
    lane = record["lane"]
    assert (lane["lane_identity"], lane["attestation_identity"]) == (
        authority.lane_identity,
        authority.attestation_identity,
    )
    assert (lane["environment_identity"], lane["image_digest"]) == (
        authority.environment_identity,
        authority.image_digest,
    )
    assert (lane["workers"], lane["native_threads"]) == (1, 1)
    specification = record["specification"]
    assert specification["candidate_counts"] == list(calibration.CANDIDATE_COUNTS)
    assert specification["scope"] == list(calibration.OUTPUT_SCOPE)
    assert specification["thresholds"] == calibration.THRESHOLDS
    assert (
        specification["summary_policy_identity"] == calibration.SUMMARY_POLICY_IDENTITY
    )
    truth = {
        "fields": list(calibration.TRUTH_FIELDS),
        "values": calibration.truth_estimands(),
    }
    assert specification["truth_identity"] == _canonical_hash(truth)
    for family, phases in calibration.ROOTS.items():
        for phase, roots in phases.items():
            assert specification["root_sets"][family][phase] == {
                "count": len(roots),
                "identity": _canonical_hash(list(roots)),
            }
    for family, result in record["families"].items():
        passing = set()
        for count in calibration.CANDIDATE_COUNTS:
            candidate = result["candidates"][str(count)]
            metrics = candidate["aggregate_metrics"]
            assert set(metrics) == set(calibration.THRESHOLDS)
            assert all(math.isfinite(value) for value in metrics.values())
            assert metrics["failure_prevalence"] == 0.0
            assert candidate["root_count"] == 16
            assert not candidate["structural_disqualifiers"]
            if candidate["status"] == "PASS":
                passing.add(count)
                assert candidate["passing_root_count"] == 16
                assert all(
                    metrics[name] <= limit
                    for name, limit in calibration.THRESHOLDS.items()
                )
            else:
                assert candidate["threshold_disqualifiers"]
                assert all(
                    item["root"] in calibration.ROOTS[family]["calibration"]
                    for item in candidate["threshold_disqualifiers"]
                )
        expected = calibration.select_candidate(passing)
        assert result["selection"] == {
            "status": expected[0],
            "selected_count": expected[1],
        }
        if expected[1] is None:
            assert result["replay"] is result["holdout"] is None
        else:
            replay, holdout = result["replay"], result["holdout"]
            assert replay["status"] == "PASS" and replay["exact_identity_match"]
            assert (replay["serial_workers"], replay["parallel_workers"]) == (1, 2)
            assert replay["root"] == calibration.ROOTS[family]["calibration"][0]
            assert holdout["status"] == "PASS"
            assert holdout["passing_root_count"] == holdout["root_count"] == 8
            assert all(
                holdout["aggregate_metrics"][name] <= limit
                for name, limit in calibration.THRESHOLDS.items()
            )
            assert not holdout["structural_disqualifiers"]
            assert not holdout["threshold_disqualifiers"]
