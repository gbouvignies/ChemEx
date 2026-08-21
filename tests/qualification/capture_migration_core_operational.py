"""Build compact operational evidence after canonical lane attestation."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Protocol, cast

from tests.qualification.native_evaluation_process import (
    evaluation_request,
    run_serial_request,
    run_spawned_requests,
    start_interruptible_spawned_request,
)

from chemex.baselines import (
    CanonicalBaselineValue,
    CaseDefinition,
    CaseSourceAuthority,
    ExecutionSpecification,
    InputMember,
    LegacyObservationImplementation,
    Occurrence,
    ResultBundle,
    ResultMember,
)
from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.migration_core_operational import OperationalReplayCapture
from chemex.numerical_lanes import canonical_lanes
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parents[2]
EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
PARAMETERS = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"
FROZEN_LEGACY_IMPLEMENTATION = LegacyObservationImplementation(
    package_version="2026.6.1",
    lmfit_version="1.3.4",
    source_manifest_hash="b16ee859a1e747043650b3d3529843856192bc8fad5d8a4c8d49b99157410ad5",
)


class LaneAuthority(Protocol):
    """Attestation capability exposed by the frozen canonical image."""

    def to_record(self) -> dict[str, object]: ...


def _member(role: str, path: Path) -> InputMember:
    content = path.read_bytes()
    return InputMember(role, hashlib.sha256(content).hexdigest(), len(content))


def _selected_evaluation(
    name: str = "1N", parameter_file: Path = PARAMETERS
) -> tuple[EvaluationEngine, EvaluationFrame, bytes]:
    method = Method()
    selection = Selection(include=[SpinSystem.from_name(name)], exclude=None)
    session = AnalysisSession.create()
    session.set_model("2st_hd")
    experiments = build_experiments([EXPERIMENT], selection, session=session)
    session.parameters.set_defaults(read_defaults([parameter_file]))
    if not session.try_build_analysis_values():
        raise RuntimeError("Operational capture parameter state is unavailable")
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    request = evaluation_request(
        model="2st_hd",
        experiment_files=[EXPERIMENT],
        parameter_files=[parameter_file],
        selection=selection,
        method=method,
        plan=engine.plan,
        frame=frame,
    )
    return engine, frame, request


def _receipt(content: bytes) -> dict[str, object]:
    record = json.loads(content)
    if not isinstance(record, dict):
        raise TypeError("Operational receipt must be a record")
    return record


def _capture_facts() -> dict[str, dict[str, object]]:
    engine, frame, request = _selected_evaluation()
    expected = engine.new_evaluator().evaluate(frame)
    if not isinstance(expected, EvaluationResult):
        raise TypeError("Operational reference evaluation failed")
    restored_plan = EvaluationPlan.from_record(engine.plan.to_record())
    restored_frame = EvaluationFrame.from_record(frame.to_record())
    restored_result = EvaluationResult.from_record(expected.to_record(), restored_plan)

    serial = _receipt(run_serial_request(request))
    spawned, replay = map(_receipt, run_spawned_requests((request, request)))
    first_engine, _first_frame, first_request = _selected_evaluation("1N")
    second_engine, _second_frame, second_request = _selected_evaluation("2N")
    forward = tuple(
        map(_receipt, run_spawned_requests((first_request, second_request)))
    )
    reverse = tuple(
        map(_receipt, run_spawned_requests((second_request, first_request)))
    )
    forward_outcomes = {item["request_identity"]: item["outcome"] for item in forward}
    reverse_outcomes = {item["request_identity"]: item["outcome"] for item in reverse}

    process, barrier = start_interruptible_spawned_request(request)
    if not barrier.poll(60.0):
        process.terminate()
        process.join()
        raise TimeoutError("Operational interruption worker did not bind")
    owner_pid = int(barrier.recv_bytes())
    process.terminate()
    process.join(60.0)
    try:
        interrupted_outcome = barrier.recv_bytes()
    except EOFError:
        interrupted_outcome = None
    barrier.close()

    cached = engine.new_evaluator()
    first = cached.evaluate(frame)
    second = cached.evaluate(frame)
    changed_id = frame._items[0][0]
    changed_frame = EvaluationFrame(
        frame.parameterization_identity,
        tuple(
            (
                param_id,
                value + max(abs(value) * 1.0e-6, 1.0e-9)
                if param_id == changed_id
                else value,
            )
            for param_id, value in frame._items
        ),
    )
    changed = cached.evaluate(changed_frame)
    if not all(isinstance(item, EvaluationResult) for item in (first, second, changed)):
        raise RuntimeError("Operational cache evaluation failed")

    malformed_record = json.loads(request)
    malformed_record["unknown"] = True
    malformed = _receipt(
        run_spawned_requests((json.dumps(malformed_record).encode("utf-8"),))[0]
    )
    with tempfile.TemporaryDirectory() as directory:
        copied = Path(directory) / "parameters.toml"
        copied.write_bytes(PARAMETERS.read_bytes())
        _stale_engine, _stale_frame, stale_request = _selected_evaluation(
            parameter_file=copied
        )
        copied.write_text(
            copied.read_text(encoding="utf-8").replace("KDH = 10.0", "KDH = 11.0"),
            encoding="utf-8",
        )
        stale = _receipt(run_spawned_requests((stale_request,))[0])

    expected_record = expected.to_record()
    spawned_pids = {spawned["worker_pid"], replay["worker_pid"]}
    plan_ids = set()
    for item in forward:
        outcome = item.get("outcome")
        if isinstance(outcome, dict):
            plan_ids.add(outcome.get("plan_identity"))
    return {
        "serialization": {
            "plan_round_trip": restored_plan == engine.plan,
            "frame_round_trip": restored_frame == frame,
            "result_round_trip": restored_result.to_record() == expected.to_record(),
        },
        "multiprocessing": {
            "fresh_worker_ownership": (
                len(spawned_pids) == 2
                and os.getpid() not in spawned_pids
                and spawned["worker_pid"] == spawned["evaluator_owner_pid"]
                and replay["worker_pid"] == replay["evaluator_owner_pid"]
            ),
            "request_order_independent": (
                forward_outcomes == reverse_outcomes
                and plan_ids
                == {first_engine.plan.identity, second_engine.plan.identity}
            ),
            "interruption_suppressed_result": (
                interrupted_outcome is None
                and process.exitcode is not None
                and process.exitcode != 0
                and process.pid == owner_pid
            ),
        },
        "cache": {
            "hits": cached.cache_statistics.hits,
            "misses": cached.cache_statistics.misses,
            "changed_frame_invalidated": changed.to_record() != expected_record,
        },
        "deterministic_replay": {
            "serial_spawn_identity_match": serial.get("outcome") == expected_record,
            "spawn_replay_identity_match": (
                spawned.get("outcome") == replay.get("outcome") == expected_record
            ),
        },
        "fail_closed": {
            "malformed_request_rejected": (
                malformed.get("status") == "rejected"
                and malformed.get("outcome") is None
            ),
            "stale_request_rejected": (
                stale.get("status") == "rejected" and stale.get("outcome") is None
            ),
        },
    }


def capture(
    *,
    source_commit: str,
    lockfile_hash: str,
    authority: LaneAuthority,
) -> OperationalReplayCapture:
    lane = canonical_lanes()[0]
    authority_record = authority.to_record()
    implementation = FROZEN_LEGACY_IMPLEMENTATION
    case = CaseDefinition.create(
        "migration-core-operational-replay",
        CaseSourceAuthority(source_commit, lockfile_hash),
        {"purpose": "serialization-multiprocessing-cache-replay"},
        (
            _member("capture-runner", Path(__file__)),
            _member(
                "attestation-first-runner",
                Path(__file__).with_name("seal_migration_core_operational.py"),
            ),
            _member(
                "process-boundary-helper",
                ROOT / "tests/qualification/native_evaluation_process.py",
            ),
            _member("experiment", EXPERIMENT),
            _member("parameters", PARAMETERS),
        ),
    )
    specification = ExecutionSpecification.create(
        case,
        implementation,
        workflow={"probe": "selected-dcest-evaluation-process-boundary"},
        lane_reference=lane.identity,
        policy={"cache": "frame-identity", "process_start": "spawn"},
        budget={"worker_response_timeout_seconds": 60},
        seed=None,
        execution_settings={
            "environment_identity": authority_record["environment_identity"],
            "workers": authority_record["workers"],
            "native_threads": authority_record["native_threads"],
        },
        artifact_inventory={"owner_records": ["operational-replay-facts"]},
        roles=("qualification:migration-core-operational-replay",),
        claims=("typed-runtime-facts",),
    )
    requested = Occurrence.requested(
        specification,
        case,
        f"migration-core-operational-replay:{source_commit}",
        cast(Any, authority),
    )
    facts = CanonicalBaselineValue.from_value(_capture_facts())
    facts_content = json.dumps(
        facts.to_record_value(),
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    bundle = ResultBundle.create(
        requested.identity,
        specification.identity,
        implementation,
        (
            ResultMember(
                "operational-replay-facts.json",
                hashlib.sha256(facts_content).hexdigest(),
                len(facts_content),
            ),
        ),
    )
    occurrence = requested.succeeded(bundle)
    return OperationalReplayCapture(
        case,
        specification,
        occurrence,
        bundle,
        facts,
    )
