"""Pre-acquisition structural contracts for optimizer calibration v2."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import cast

import pytest

import tests.qualification.capture_optimizer_calibration_v2 as calibration

ROOT = Path(__file__).parent.parent
FAILED_ROUND = ROOT / "tests/fixtures/canonical_optimizer_calibration_failure_v2.json"
FAILED_LAUNCH_2 = (
    ROOT / "tests/fixtures/canonical_optimizer_calibration_failure_v2_launch2.json"
)


def _matched_truth() -> dict[str, object]:
    return {
        "qualification": "REFERENCE_MATCHED",
        "artifacts": (
            {"probe_id": "trf-routine-quadratic-v1"},
            {"probe_id": "trf-difficult-rosenbrock-v1"},
            {"probe_id": "grid-27-seed-coverage-v1"},
            {"probe_id": "grid-candidate-ordering-v1"},
        ),
    }


def test_v2_changes_only_instrument_authority_and_keeps_candidates_frozen() -> None:
    assert calibration.SPECIFICATION_ID == "chemex-optimizer-calibration-v2"
    assert calibration.specification_identity() == (
        "c8db7972f1d43262b3935e33c2b711d6f392b5f3fcb9e78339079de31ca0c9c3"
    )
    assert (
        calibration.SPECIFICATION["candidates"]
        == calibration.v1.SPECIFICATION["candidates"]
    )
    assert (
        calibration.SPECIFICATION["grid_coordinates"]
        == (calibration.v1.SPECIFICATION["grid_coordinates"])
    )
    assert (
        calibration.SPECIFICATION["de_coordinates"]
        == calibration.v1.SPECIFICATION["de_coordinates"]
    )
    assert (
        calibration.SPECIFICATION["de_roots"]
        == calibration.v1.SPECIFICATION["de_roots"]
    )
    assert tuple(
        calibration.derive_de_root("calibration", index) for index in range(2)
    ) == tuple(
        calibration.v1.derive_de_root("calibration", index) for index in range(2)
    )
    repairs = cast("dict[str, object]", calibration.SPECIFICATION["v2_repairs"])
    assert "not selected from v1 candidate outcomes" in str(repairs["authority"])


def test_rejected_candidate_survives_final_assembly_and_content_identity() -> None:
    rejected = calibration.rejected_candidate_record(
        ordinal=0,
        budget=8,
        terminal="budget_exhausted",
        reasons=("typed_terminal_not_accepted",),
    )
    record = calibration.assemble_record(
        specification_commit="a" * 40,
        source={"kind": "synthetic-preflight"},
        lane_records={"kind": "synthetic-native-macos-preflight"},
        truth_probes=_matched_truth(),
        strata={"routine_direct_trf": {"candidates": (rejected,)}},
        elapsed_seconds=0.0,
    )
    encoded = json.dumps(record, allow_nan=False, sort_keys=True)
    restored = json.loads(encoded)
    candidate = restored["strata"]["routine_direct_trf"]["candidates"][0]
    assert candidate["objective"] is None
    assert candidate["vector"] is None
    identity = record.pop("identity")
    assert identity == calibration._identity("canonical-optimizer-calibration", record)


def test_resource_ceiling_sums_objective_requests_and_fails_closed() -> None:
    stratum = {
        "candidates": (
            {
                "local_executions": (
                    {
                        "counters": {"objective_requests_accepted": 7},
                        "backend": {"nfev": 99},
                    },
                ),
                "materialization_accounting": {
                    "candidate_materialization_requests": 2,
                    "authoritative_root_materialization_requests": 1,
                },
            },
        ),
        "operational": {"elapsed_seconds_diagnostic_only": 10_000.0},
    }
    exact = calibration.enforce_resource_ceiling(stratum, 10)
    exceeded = calibration.enforce_resource_ceiling(stratum, 9)
    assert exact == {
        "status": "PASS",
        "ceiling": 10,
        "solver_objective_requests": 7,
        "materialization_objective_requests": 3,
        "total_objective_requests": 10,
        "backend_counters_included": False,
        "wall_time_included": False,
    }
    assert exceeded["status"] == "RESOURCE_CEILING_EXCEEDED"


def test_direct_scientific_adequacy_requires_case_specific_truth() -> None:
    unsupported = calibration.direct_scientific_adequacy(
        _matched_truth(), representative_truth_authority=None
    )
    supported = calibration.direct_scientific_adequacy(
        _matched_truth(), representative_truth_authority="synthetic-truth-v1"
    )
    assert unsupported["status"] == (
        "UNSUPPORTED_INSUFFICIENT_PROSPECTIVE_TRUTH_AUTHORITY"
    )
    assert supported == {"status": "ADEQUATE", "authority": "synthetic-truth-v1"}


def test_grouped_holdout_is_independent_and_replay_remains_replay() -> None:
    assert calibration.validate_grouped_holdout_selection() == "RELAXATION_NZ"
    holdouts = cast("dict[str, object]", calibration.SPECIFICATION["holdouts"])
    assert holdouts["decomposed_grid"] == "RELAXATION_NZ:all-five-profiles"
    repairs = cast("dict[str, object]", calibration.SPECIFICATION["v2_repairs"])
    grouped = cast("dict[str, object]", repairs["grouped_holdout_selection"])
    assert grouped["calibration_replay_is_holdout"] is False


def test_source_guard_rejects_dirty_tracked_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "checkout"
    tracked = root / "src/example.py"
    tracked.parent.mkdir(parents=True)
    tracked.write_text("frozen\n", encoding="utf-8")
    (root / ".git").mkdir()
    (root / ".git/HEAD").write_text("a" * 40 + "\n", encoding="ascii")
    manifest = root / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "files": {"src/example.py": hashlib.sha256(b"frozen\n").hexdigest()},
            },
            allow_nan=False,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    manifest_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()
    monkeypatch.setattr(calibration, "ROOT", root)
    monkeypatch.setattr(calibration.v1, "ROOT", root)
    assert (
        calibration.validate_source_checkout(
            "a" * 40,
            manifest_path=manifest,
            expected_manifest_sha256=manifest_hash,
        )["tracked_file_count"]
        == 1
    )
    tracked.write_text("dirty\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="tracked source differs"):
        calibration.validate_source_checkout(
            "a" * 40,
            manifest_path=manifest,
            expected_manifest_sha256=manifest_hash,
        )


def test_grid_qualification_fails_closed_on_ceiling_or_holdout() -> None:
    base: dict[str, object] = {
        "selection": {"status": "SELECTED"},
        "replay": {"status": "PASS"},
        "holdout": {
            "status": "PASS",
            "kind": "independent-untouched-case",
            "case": "RELAXATION_NZ",
        },
        "candidates": ({"counters": {"objective_requests_accepted": 11}},),
    }
    exceeded = calibration._qualify_grid(
        dict(base), ceiling=10, truth_probes=_matched_truth(), grouped=True
    )
    failed_holdout = dict(base)
    failed_holdout["holdout"] = {
        "status": "DISQUALIFIED",
        "kind": "independent-untouched-case",
        "case": "RELAXATION_NZ",
    }
    holdout_result = calibration._qualify_grid(
        failed_holdout, ceiling=20, truth_probes=_matched_truth(), grouped=True
    )
    exceeded_qualification = cast("dict[str, object]", exceeded["qualification"])
    exceeded_resource = cast("dict[str, object]", exceeded["resource_accounting"])
    holdout_qualification = cast("dict[str, object]", holdout_result["qualification"])
    assert exceeded_qualification["status"] == "UNSUPPORTED"
    assert exceeded_resource["status"] == "RESOURCE_CEILING_EXCEEDED"
    assert holdout_qualification["status"] == "UNSUPPORTED"


def test_native_preflight_exercises_acquisition_through_final_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "checkout"
    tracked = root / "src/instrument.py"
    tracked.parent.mkdir(parents=True)
    tracked.write_text("prospective-v2\n", encoding="utf-8")
    (root / ".git").mkdir()
    (root / ".git/HEAD").write_text("b" * 40 + "\n", encoding="ascii")
    manifest = root / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "files": {
                    "src/instrument.py": hashlib.sha256(b"prospective-v2\n").hexdigest()
                },
            },
            allow_nan=False,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    manifest_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()
    synthetic_grid = {
        "selection": {"status": "SELECTED", "budget_per_seed": 16},
        "replay": {"status": "PASS"},
        "holdout": {
            "status": "PASS",
            "kind": "independent-untouched-case",
            "case": "RELAXATION_NZ",
        },
        "candidates": ({"counters": {"objective_requests_accepted": 1}},),
    }
    monkeypatch.setattr(calibration, "ROOT", root)
    monkeypatch.setattr(calibration, "SOURCE_MANIFEST", manifest)
    monkeypatch.setattr(calibration.v1, "ROOT", root)
    monkeypatch.setattr(
        calibration.v1,
        "_attest",
        lambda _digest: (object(), {"lane": "synthetic-native-macos-preflight"}),
    )
    monkeypatch.setattr(
        calibration.v1, "_qualification_probes", lambda _authority: _matched_truth()
    )
    monkeypatch.setattr(
        calibration.v1,
        "_file_hash",
        lambda _path: hashlib.sha256(b"preflight").hexdigest(),
    )
    monkeypatch.setattr(
        calibration, "_calibrate_coupled_grid", lambda: dict(synthetic_grid)
    )
    monkeypatch.setattr(
        calibration, "_calibrate_grouped_grid", lambda: dict(synthetic_grid)
    )

    record = calibration.acquire("sha256:" + "c" * 64, "b" * 40, manifest_hash)

    assert record["identity"] == calibration._identity(
        "canonical-optimizer-calibration",
        {key: value for key, value in record.items() if key != "identity"},
    )
    strata = cast("dict[str, dict[str, object]]", record["strata"])
    routine_qualification = cast(
        "dict[str, object]", strata["routine_direct_trf"]["qualification"]
    )
    assert routine_qualification["status"] == "UNSUPPORTED"
    assert strata["decomposed_grouped_grid_trf"]["qualification"] == {
        "status": "QUALIFIED",
        "truth_authority_matched": True,
        "replay_passed": True,
        "holdout_passed": True,
        "independent_holdout": True,
        "resource_ceiling_passed": True,
    }
    json.dumps(record, allow_nan=False)


def test_v2_import_failure_is_content_identified_and_ineligible() -> None:
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
    assert record["failure"] == {
        "error_type": "ModuleNotFoundError",
        "message": "No module named 'tests'",
        "stage": "qualification-module-import",
    }
    assert record["acquisition"]["scientific_execution_started"] is False
    assert record["scientific_stages"]["grouped_holdout"] == "UNOPENED"
    assert record["eligibility"]["eligible_claims"] == []
    assert record["eligibility"]["migration_core_coverage_change"] == 0


def test_v2_lane_attestation_failure_is_separate_and_pre_scientific() -> None:
    first_launch_bytes = FAILED_ROUND.read_bytes()
    record = json.loads(FAILED_LAUNCH_2.read_text(encoding="utf-8"))
    identity = record.pop("identity")
    encoded = json.dumps(
        record,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")

    assert identity == hashlib.sha256(encoded).hexdigest()
    assert record["failure"] == {
        "error_type": "LaneAuthorityError",
        "message": "build-context manifest does not describe the running filesystem",
        "stage": "canonical-lane-attestation",
    }
    assert record["acquisition"]["in_container_source_guard"] == "PASS"
    assert (
        record["acquisition"]["authoritative_scientific_acquisition_started"] is False
    )
    assert record["scientific_stages"]["truth_probes"] == "NOT_STARTED"
    assert record["scientific_stages"]["grouped_holdout"] == "UNOPENED"
    assert record["eligibility"]["eligible_claims"] == []
    assert first_launch_bytes == FAILED_ROUND.read_bytes()
