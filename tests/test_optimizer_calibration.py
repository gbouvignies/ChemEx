"""Executable prospective-calibration contracts for issue #602."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import cast

import pytest

import tests.qualification.capture_optimizer_calibration as calibration
from tests.qualification.capture_optimizer_calibration import (
    SPECIFICATION,
    SPECIFICATION_ID,
    candidate_ordering_key,
    derive_de_root,
    holdout_decision,
    select_monotone_budget,
    specification_identity,
    validate_specification_commit,
)

ROOT = Path(__file__).parent.parent
FAILED_ROUND = ROOT / "tests/fixtures/canonical_optimizer_calibration_failure_v1.json"


def test_prospective_specification_freezes_every_acquisition_boundary() -> None:
    assert SPECIFICATION_ID == "chemex-optimizer-calibration-v1"
    assert (
        specification_identity()
        == "7210f451e44ce5311afdea46dc5131afd3328919c728c27128358e2b08ee91d8"
    )
    assert set(SPECIFICATION) == {
        "truth_cases",
        "difficult_starts",
        "candidates",
        "grid_coordinates",
        "de_coordinates",
        "de_roots",
        "disqualifiers",
        "selection",
        "accounting",
        "resource_limits",
        "holdouts",
        "unsupported",
    }
    accounting = cast("dict[str, object]", SPECIFICATION["accounting"])
    selection = cast("dict[str, object]", SPECIFICATION["selection"])
    grid = cast("dict[str, object]", SPECIFICATION["grid_coordinates"])
    assert accounting["authoritative_unit"] == "chemex-objective-request"
    assert accounting["retries"] == 0
    assert selection["runner_up_after_holdout"] == "forbidden"
    assert grid["immutable_grid_seed_count"] == 27


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


def test_failed_authoritative_round_is_content_identified_and_ineligible() -> None:
    record = json.loads(FAILED_ROUND.read_text(encoding="utf-8"))
    identity = record.pop("identity")
    encoded = json.dumps(
        record,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")

    assert identity == hashlib.sha256(encoded).hexdigest()
    assert record["status"] == "ARCHITECTURAL_FAILURE"
    assert record["eligibility"] == {
        "eligible_claims": [],
        "migration_core_coverage_change": 0,
        "reason": "authoritative result record was not content-identified or serialized",
        "status": "INELIGIBLE_NO_RESULT_RECORD",
    }
    assert record["acquisition"]["authoritative_round"] == 1
    assert record["acquisition"]["retry_count"] == 0
    assert record["acquisition"]["retuned"] is False
