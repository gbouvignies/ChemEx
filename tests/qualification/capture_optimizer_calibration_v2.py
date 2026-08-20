"""One replacement prospective optimizer calibration for issue #602.

V2 repairs only defects in the v1 qualification instrument.  Its candidate
budgets, GRID coordinates, DE topology, roots, ordering, ties, and scientific
tolerances are inherited unchanged from v1.  No v1 optimizer observation is an
authority for this specification.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import time
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePosixPath
from typing import Any, cast

import tests.qualification.capture_optimizer_calibration as v1

from chemex.baselines import LegacyObservationImplementation
from chemex.numerical_lanes import RuntimeEnvironment, prospective_lanes

ROOT = Path(__file__).parents[2]
SPECIFICATION_ID = "chemex-optimizer-calibration-v2"
EXPECTED_CANONICAL_LANE_IDENTITY = (
    "953168c14885b9278a71dadf694633dd10cf3740bedfe00c4abb706fc0974329"
)
EXPECTED_SPECIFICATION_IDENTITY = (
    "c8db7972f1d43262b3935e33c2b711d6f392b5f3fcb9e78339079de31ca0c9c3"
)
SOURCE_MANIFEST = ROOT / "tests/fixtures/optimizer_calibration_source_manifest_v2.json"

SPECIFICATION = copy.deepcopy(v1.SPECIFICATION)
cast("dict[str, object]", SPECIFICATION["holdouts"])["decomposed_grid"] = (
    "RELAXATION_NZ:all-five-profiles"
)
SPECIFICATION["v2_repairs"] = {
    "authority": (
        "instrument-and-specification-correctness-only; not selected from v1 "
        "candidate outcomes"
    ),
    "canonical_rejection": (
        "nullable-objective-and-vector-with-explicit-terminal-and-disqualifiers"
    ),
    "direct_scientific_adequacy": {
        "authority": "issue-591-reference-matched-truth-probes",
        "routine_probe": "trf-routine-quadratic-v1",
        "difficult_probe": "trf-difficult-rosenbrock-v1",
        "representative_case_rule": (
            "unsupported-without-pre-existing-case-specific-truth-authority"
        ),
    },
    "round_resource_ceiling": (
        "sum-accepted-chemex-objective-requests-and-materialization-requests-"
        "per-stratum; fail-closed"
    ),
    "grouped_holdout_selection": {
        "rule": (
            "lexicographically-first-shipped-RELAXATION-star-case-with-complete-"
            "experiment-method-parameter-inputs-after-excluding-RELAXATION_HZNZ"
        ),
        "selected": "RELAXATION_NZ",
        "calibration_replay_is_holdout": False,
    },
    "source_guard": (
        "detached-freeze-commit-plus-content-hashed-tracked-source-manifest"
    ),
}

GROUPED_HOLDOUT_DIRECTORY = ROOT / "examples/Experiments/RELAXATION_NZ"
GROUPED_HOLDOUT_EXPERIMENT = GROUPED_HOLDOUT_DIRECTORY / "Experiments/800mhz.toml"
GROUPED_HOLDOUT_PARAMETERS = GROUPED_HOLDOUT_DIRECTORY / "Parameters/parameters.toml"
GROUPED_HOLDOUT_METHOD = GROUPED_HOLDOUT_DIRECTORY / "Methods/method.toml"


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def _identity(kind: str, value: object) -> str:
    return hashlib.sha256(
        _canonical_json({"kind": kind, "schema_version": 1, "record": value})
    ).hexdigest()


def specification_identity() -> str:
    return _identity(SPECIFICATION_ID, SPECIFICATION)


def derive_de_root(phase: str, index: int) -> int:
    """Keep the v1 roots exactly; v1 outcomes have no influence on v2."""
    return v1.derive_de_root(phase, index)


def _attest(image_digest: str) -> tuple[object, dict[str, object]]:
    """Attest only the published #652 canonical Python 3.13 v2 lane."""
    manifest_directory = ROOT / "src/chemex/numerical_lanes/manifests"
    lane = prospective_lanes(manifest_directory)[0]
    if lane.identity != EXPECTED_CANONICAL_LANE_IDENTITY:
        raise RuntimeError("published canonical v2 lane identity changed")
    authority = lane.attest_current_process(image_digest)
    records = {
        "numerical_lane": lane.to_record(),
        "lane_attestation": authority.to_record(),
        "runtime_environment": RuntimeEnvironment(lane.semantics).to_record(),
    }
    return authority, records


def _implementation_wheel_sha256() -> str:
    digest = os.environ.get("CHEMEX_IMPLEMENTATION_WHEEL_SHA256", "")
    if (
        len(digest) != 64
        or digest != digest.lower()
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise RuntimeError("exact implementation wheel authority is unavailable")
    return digest


def _source_authority(
    specification_commit: str, expected_manifest_sha256: str
) -> dict[str, object]:
    source_guard = validate_source_checkout(
        specification_commit,
        manifest_path=SOURCE_MANIFEST,
        expected_manifest_sha256=expected_manifest_sha256,
    )
    return {
        **source_guard,
        "implementation": LegacyObservationImplementation.from_current_package().to_record(),
        "implementation_wheel_sha256": _implementation_wheel_sha256(),
        "qualification_script_sha256": v1._file_hash(Path(__file__)),
        "qualification_test_sha256": v1._file_hash(
            ROOT / "tests/test_optimizer_calibration_v2.py"
        ),
        "dependency_lock_sha256": v1._file_hash(ROOT / "uv.lock"),
    }


def preflight(
    image_digest: str,
    specification_commit: str,
    expected_manifest_sha256: str,
) -> dict[str, object]:
    """Prove all transport and authority boundaries without running science."""
    source = _source_authority(specification_commit, expected_manifest_sha256)
    if specification_identity() != EXPECTED_SPECIFICATION_IDENTITY:
        raise RuntimeError("frozen v2 acquisition specification changed")
    _authority, lane_records = _attest(image_digest)
    lane = cast("Mapping[str, object]", lane_records.get("numerical_lane"))
    if lane.get("identity") != EXPECTED_CANONICAL_LANE_IDENTITY:
        raise RuntimeError("attested lane is not the published canonical v2 lane")
    record: dict[str, object] = {
        "schema_version": 1,
        "record_version": "optimizer-calibration-v2-preflight",
        "status": "READY_FOR_SCIENTIFIC_ACQUISITION",
        "scientific_execution": "NOT_STARTED",
        "specification": {
            "id": SPECIFICATION_ID,
            "identity": EXPECTED_SPECIFICATION_IDENTITY,
        },
        "source": source,
        "canonical_lane": lane_records,
    }
    record["identity"] = _identity("optimizer-calibration-v2-preflight", record)
    return record


def rejected_candidate_record(
    *, ordinal: int, budget: int, terminal: str, reasons: Sequence[str]
) -> dict[str, object]:
    """Represent a rejected candidate without non-finite numeric sentinels."""
    return {
        "ordinal": ordinal,
        "budget": budget,
        "status": "DISQUALIFIED",
        "reasons": tuple(reasons),
        "terminal": terminal,
        "objective": None,
        "vector": None,
        "counters": {
            "solver_requests_received": 0,
            "objective_requests_accepted": 0,
            "objective_evaluations_completed": 0,
        },
        "backend": None,
        "execution_identity": None,
        "accepted_identity": None,
    }


def direct_scientific_adequacy(
    truth_probes: Mapping[str, object], *, representative_truth_authority: str | None
) -> dict[str, object]:
    """Require existing #591 truth plus case-specific representative authority."""
    artifacts = cast(
        "Sequence[Mapping[str, object]]", truth_probes.get("artifacts", ())
    )
    probe_ids = {str(item.get("probe_id")) for item in artifacts}
    required = {
        "trf-routine-quadratic-v1",
        "trf-difficult-rosenbrock-v1",
    }
    truth_matched = (
        truth_probes.get("qualification") == "REFERENCE_MATCHED"
        and required <= probe_ids
    )
    if not truth_matched:
        return {
            "status": "UNSUPPORTED_TRUTH_AUTHORITY_NOT_MATCHED",
            "authority": None,
        }
    if representative_truth_authority is None:
        return {
            "status": "UNSUPPORTED_INSUFFICIENT_PROSPECTIVE_TRUTH_AUTHORITY",
            "authority": None,
        }
    return {"status": "ADEQUATE", "authority": representative_truth_authority}


def validate_grouped_holdout_selection() -> str:
    """Re-derive the frozen independent holdout without observing results."""
    eligible = tuple(
        path.name
        for path in sorted(ROOT.glob("examples/Experiments/RELAXATION_*"))
        if path.name != "RELAXATION_HZNZ"
        and (path / "Experiments/800mhz.toml").is_file()
        and (path / "Methods/method.toml").is_file()
        and (path / "Parameters/parameters.toml").is_file()
    )
    if not eligible or eligible[0] != "RELAXATION_NZ":
        raise RuntimeError("frozen independent grouped holdout selection changed")
    return eligible[0]


def validate_source_checkout(
    specification_commit: str,
    *,
    manifest_path: Path = SOURCE_MANIFEST,
    expected_manifest_sha256: str,
) -> dict[str, object]:
    """Reject detached-checkout tracked bytes that differ from the freeze manifest."""
    v1.validate_specification_commit(specification_commit)
    manifest_bytes = manifest_path.read_bytes()
    observed_manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
    if observed_manifest_sha256 != expected_manifest_sha256:
        raise RuntimeError("authoritative source manifest identity mismatch")
    manifest = json.loads(manifest_bytes)
    if manifest.get("schema_version") != 1:
        raise RuntimeError("authoritative source manifest schema mismatch")
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise RuntimeError("authoritative source manifest is empty")
    for relative_value, expected_value in sorted(files.items()):
        if not isinstance(relative_value, str) or not isinstance(expected_value, str):
            raise TypeError("authoritative source manifest contains an unsafe entry")
        relative = relative_value
        expected = expected_value
        pure = PurePosixPath(relative)
        if pure.is_absolute() or ".." in pure.parts:
            raise RuntimeError("authoritative source manifest contains an unsafe entry")
        path = ROOT.joinpath(*pure.parts)
        if (
            not path.is_file()
            or hashlib.sha256(path.read_bytes()).hexdigest() != expected
        ):
            raise RuntimeError(f"authoritative tracked source differs: {relative}")
    return {
        "manifest_path": str(manifest_path.relative_to(ROOT)),
        "manifest_sha256": observed_manifest_sha256,
        "tracked_file_count": len(files),
    }


def _materialization_accounting(
    record: Mapping[str, object], kind: str
) -> dict[str, int]:
    attempts = cast("Sequence[Mapping[str, object]]", record.get("attempts", ()))
    candidate_key = "candidate_identity" if kind == "coupled" else "aggregate_identity"
    candidate_requests = sum(item.get(candidate_key) is not None for item in attempts)
    root = record.get("root_materialization")
    root_requests = (
        0
        if not isinstance(root, Mapping)
        else int(cast("int", root.get("evaluation_count", 0)))
    )
    return {
        "candidate_materialization_requests": candidate_requests,
        "authoritative_root_materialization_requests": root_requests,
    }


def _attach_materialization_accounting(
    record: dict[str, object], kind: str
) -> dict[str, object]:
    record["materialization_accounting"] = _materialization_accounting(record, kind)
    return record


def _sum_accepted_requests(value: object) -> int:
    if isinstance(value, Mapping):
        total = 0
        for key, child in value.items():
            if key == "objective_requests_accepted":
                total += int(cast("int", child))
            elif key != "materialization_accounting":
                total += _sum_accepted_requests(child)
        return total
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return sum(_sum_accepted_requests(child) for child in value)
    return 0


def enforce_resource_ceiling(
    stratum: Mapping[str, object], ceiling: int
) -> dict[str, object]:
    materializations = 0

    def visit(value: object) -> None:
        nonlocal materializations
        if isinstance(value, Mapping):
            accounting = value.get("materialization_accounting")
            if isinstance(accounting, Mapping):
                materializations += sum(
                    int(cast("int", item)) for item in accounting.values()
                )
            for child in value.values():
                visit(child)
        elif isinstance(value, Sequence) and not isinstance(
            value, (str, bytes, bytearray)
        ):
            for child in value:
                visit(child)

    visit(stratum)
    solver_requests = _sum_accepted_requests(stratum)
    total = solver_requests + materializations
    return {
        "status": "PASS" if total <= ceiling else "RESOURCE_CEILING_EXCEEDED",
        "ceiling": ceiling,
        "solver_objective_requests": solver_requests,
        "materialization_objective_requests": materializations,
        "total_objective_requests": total,
        "backend_counters_included": False,
        "wall_time_included": False,
    }


def _selection_status(record: Mapping[str, object]) -> str | None:
    selection = record.get("selection")
    return (
        None
        if not isinstance(selection, Mapping)
        else cast("str | None", selection.get("status"))
    )


def _qualify_grid(
    record: dict[str, object],
    *,
    ceiling: int,
    truth_probes: Mapping[str, object],
    grouped: bool,
) -> dict[str, object]:
    resource = enforce_resource_ceiling(record, ceiling)
    probe_ids = {
        str(item.get("probe_id"))
        for item in cast(
            "Sequence[Mapping[str, object]]", truth_probes.get("artifacts", ())
        )
    }
    truth_matched = (
        truth_probes.get("qualification") == "REFERENCE_MATCHED"
        and {"grid-27-seed-coverage-v1", "grid-candidate-ordering-v1"} <= probe_ids
    )
    holdout = record.get("holdout")
    holdout_passed = isinstance(holdout, Mapping) and holdout.get("status") == "PASS"
    replay = record.get("replay")
    replay_passed = isinstance(replay, Mapping) and replay.get("status") == "PASS"
    selected = _selection_status(record) == "SELECTED"
    independent = not grouped or (
        isinstance(holdout, Mapping)
        and holdout.get("case") == "RELAXATION_NZ"
        and holdout.get("kind") == "independent-untouched-case"
    )
    passed = (
        truth_matched
        and selected
        and replay_passed
        and holdout_passed
        and independent
        and resource["status"] == "PASS"
    )
    record["resource_accounting"] = resource
    record["qualification"] = {
        "status": "QUALIFIED" if passed else "UNSUPPORTED",
        "truth_authority_matched": truth_matched,
        "replay_passed": replay_passed,
        "holdout_passed": holdout_passed,
        "independent_holdout": independent,
        "resource_ceiling_passed": resource["status"] == "PASS",
    }
    return record


def _build_grouped_holdout() -> tuple[Any, Any, Any, Any, Any]:
    from chemex.configuration.methods import Selection, read_methods
    from chemex.configuration.parameters import read_defaults
    from chemex.evaluation.native import EvaluationEngine
    from chemex.experiments.builder import build_experiments
    from chemex.optimize.direct_trf import OptimizationProblem
    from chemex.runtime import AnalysisSession

    validate_grouped_holdout_selection()
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [GROUPED_HOLDOUT_EXPERIMENT],
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([GROUPED_HOLDOUT_PARAMETERS]))
    if not session.try_build_analysis_values():
        raise RuntimeError("RELAXATION_NZ native parameter construction failed")
    method = read_methods([GROUPED_HOLDOUT_METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    if configuration is None:
        raise RuntimeError("RELAXATION_NZ configuration did not seal")
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, session.analysis_values.snapshot()
    )
    return session, experiments, parameterization, engine, problem


def _coupled_run(
    problem: Any, parameterization: Any, engine: Any, budget: int
) -> dict[str, object]:
    return _attach_materialization_accounting(
        v1._run_coupled_grid(problem, parameterization, engine, budget), "coupled"
    )


def _grouped_run(
    problem: Any,
    parameterization: Any,
    engine: Any,
    decomposition: Any,
    budget: int,
) -> dict[str, object]:
    return _attach_materialization_accounting(
        v1._run_grouped_grid(problem, parameterization, engine, decomposition, budget),
        "grouped",
    )


def _calibrate_coupled_grid() -> dict[str, object]:
    _s, _e, parameterization, engine, problem = v1._build_relaxation("G2N-HN")
    budgets = cast(
        "tuple[int, ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])[
            "coupled_grid_budgets_per_seed"
        ],
    )
    records = tuple(
        _coupled_run(problem, parameterization, engine, budget) for budget in budgets
    )
    passing = {
        cast("int", item["budget_per_seed"])
        for item in records
        if item["status"] == "PASS"
    }
    status, selected = v1.select_monotone_budget(budgets, passing)
    result: dict[str, object] = {
        "candidates": records,
        "selection": {"status": status, "budget_per_seed": selected},
        "replay": None,
        "holdout": None,
    }
    if selected is None:
        return result
    chosen = next(item for item in records if item["budget_per_seed"] == selected)
    replay = _coupled_run(problem, parameterization, engine, selected)
    replay_passed = v1._grid_signature(replay) == v1._grid_signature(chosen)
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {
            "status": "REPLAY_MISMATCH_UNSUPPORTED",
            "budget_per_seed": None,
        }
        return result
    _hs, _he, hp, heng, hproblem = v1._build_relaxation("Y3N-HN")
    holdout = _coupled_run(hproblem, hp, heng, selected)
    holdout_status, final = v1.holdout_decision(
        str(selected), (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "budget_per_seed": None if final is None else int(final),
    }
    return result


def _calibrate_grouped_grid() -> dict[str, object]:
    from chemex.optimize.grouped_direct_trf import FitDecomposition

    _s, _e, parameterization, engine, problem = v1._build_relaxation(None, grouped=True)
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    budgets = cast(
        "tuple[int, ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])[
            "decomposed_grid_budgets_per_component"
        ],
    )
    records = tuple(
        _grouped_run(problem, parameterization, engine, decomposition, budget)
        for budget in budgets
    )
    passing = {
        cast("int", item["budget_per_component"])
        for item in records
        if item["status"] == "PASS"
    }
    status, selected = v1.select_monotone_budget(budgets, passing)
    result: dict[str, object] = {
        "candidates": records,
        "selection": {"status": status, "budget_per_component": selected},
        "replay": None,
        "holdout": None,
    }
    if selected is None:
        return result
    chosen = next(item for item in records if item["budget_per_component"] == selected)
    replay = _grouped_run(problem, parameterization, engine, decomposition, selected)
    replay_passed = v1._grid_signature(replay) == v1._grid_signature(chosen)
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "kind": "calibration-case-exact-replay",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {
            "status": "REPLAY_MISMATCH_UNSUPPORTED",
            "budget_per_component": None,
        }
        return result
    _hs, _he, hp, heng, hproblem = _build_grouped_holdout()
    hdecomposition = FitDecomposition.from_root(hproblem, hp, heng)
    holdout = _grouped_run(hproblem, hp, heng, hdecomposition, selected)
    holdout["kind"] = "independent-untouched-case"
    holdout["case"] = "RELAXATION_NZ"
    holdout_status, final = v1.holdout_decision(
        str(selected), (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "budget_per_component": None if final is None else int(final),
    }
    return result


def _unsupported_direct(
    truth_probes: Mapping[str, object], kind: str, ceiling: int
) -> dict[str, object]:
    record: dict[str, object] = {
        "candidates": (),
        "selection": {"status": "UNSUPPORTED", "policy": None, "budget": None},
        "replay": None,
        "holdout": None,
        "scientific_adequacy": direct_scientific_adequacy(
            truth_probes, representative_truth_authority=None
        ),
        "qualification": {
            "status": "UNSUPPORTED",
            "reason": (
                f"{kind} representative case lacks pre-existing case-specific truth "
                "authority; v1 observations are ineligible"
            ),
        },
    }
    record["resource_accounting"] = enforce_resource_ceiling(record, ceiling)
    return record


def _unsupported_de(ceiling: int) -> dict[str, object]:
    record: dict[str, object] = {
        "candidates": (),
        "selection": {"status": "DIRECT_REFERENCE_UNQUALIFIED", "topology": None},
        "replay": None,
        "holdout": None,
        "qualification": {
            "status": "UNSUPPORTED",
            "reason": "Direct DCEST reference is scientifically unqualified",
        },
    }
    record["resource_accounting"] = enforce_resource_ceiling(record, ceiling)
    return record


def assemble_record(
    *,
    specification_commit: str,
    source: Mapping[str, object],
    lane_records: Mapping[str, object],
    truth_probes: Mapping[str, object],
    strata: Mapping[str, object],
    elapsed_seconds: float,
) -> dict[str, object]:
    record: dict[str, object] = {
        "schema_version": 1,
        "record_version": "canonical-optimizer-calibration-v2",
        "specification": {
            "id": SPECIFICATION_ID,
            "identity": specification_identity(),
            "record": SPECIFICATION,
        },
        "source": dict(source),
        "canonical_lane": dict(lane_records),
        "authoritative_acquisition": {
            "version": 2,
            "round": 1,
            "retry_count": 0,
            "retuned": False,
            "adaptive_extension": False,
            "v1_observations_used": False,
            "specification_commit": specification_commit,
        },
        "truth_probes": dict(truth_probes),
        "strata": dict(strata),
        "operational": {
            "elapsed_seconds_diagnostic_only": elapsed_seconds,
            "workers": 1,
            "native_threads": 1,
        },
    }
    record["identity"] = _identity("canonical-optimizer-calibration", record)
    _canonical_json(record)
    return record


def _guard(name: str, function: Any) -> dict[str, object]:
    try:
        value = function()
    except Exception as error:  # noqa: BLE001 - the one round preserves failure
        return {
            "status": "ARCHITECTURAL_FAILURE",
            "stage": name,
            "error_type": type(error).__name__,
            "message": str(error),
            "qualification": {"status": "UNSUPPORTED"},
        }
    if not isinstance(value, dict):
        raise TypeError(f"{name} did not return a record")
    return value


def acquire(
    image_digest: str,
    specification_commit: str,
    expected_manifest_sha256: str,
) -> dict[str, object]:
    source = _source_authority(specification_commit, expected_manifest_sha256)
    if specification_identity() != EXPECTED_SPECIFICATION_IDENTITY:
        raise RuntimeError("frozen v2 acquisition specification changed")
    authority, lane_records = _attest(image_digest)
    started = time.perf_counter()
    probes = _guard("truth_probes", lambda: v1._qualification_probes(authority))
    limits = cast("Mapping[str, int]", SPECIFICATION["resource_limits"])
    routine = _unsupported_direct(
        probes, "routine Direct TRF", limits["routine_trf_objective_requests"]
    )
    difficult = _unsupported_direct(
        probes,
        "difficult-start Direct TRF",
        limits["difficult_trf_objective_requests"],
    )
    coupled = _guard(
        "coupled_grid",
        lambda: _qualify_grid(
            _calibrate_coupled_grid(),
            ceiling=limits["coupled_grid_objective_requests"],
            truth_probes=probes,
            grouped=False,
        ),
    )
    grouped = _guard(
        "decomposed_grid",
        lambda: _qualify_grid(
            _calibrate_grouped_grid(),
            ceiling=limits["decomposed_grid_objective_requests"],
            truth_probes=probes,
            grouped=True,
        ),
    )
    de = _unsupported_de(limits["de_trf_objective_requests"])
    return assemble_record(
        specification_commit=specification_commit,
        source=source,
        lane_records=cast("Mapping[str, object]", lane_records),
        truth_probes=probes,
        strata={
            "routine_direct_trf": routine,
            "difficult_direct_trf": difficult,
            "coupled_grid_trf": coupled,
            "decomposed_grouped_grid_trf": grouped,
            "de_trf": de,
        },
        elapsed_seconds=time.perf_counter() - started,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--specification-commit", required=True)
    parser.add_argument("--source-manifest-sha256", required=True)
    parser.add_argument("--preflight", action="store_true")
    arguments = parser.parse_args()
    function = preflight if arguments.preflight else acquire
    record = function(
        arguments.image_digest,
        arguments.specification_commit,
        arguments.source_manifest_sha256,
    )
    arguments.output.write_text(
        json.dumps(record, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
