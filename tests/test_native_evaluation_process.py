"""Real spawn-boundary qualification for the native evaluation records (#593)."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict

from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession
from tests.qualification.native_evaluation_process import (
    evaluation_request,
    run_serial_request,
    run_spawned_requests,
    start_interruptible_spawned_request,
)

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
PARAMETERS = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"


class EvaluationReceipt(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    schema_version: Literal[1]
    request_identity: str
    status: Literal["completed", "rejected"]
    worker_pid: int
    evaluator_owner_pid: int | None
    outcome_kind: Literal["result", "failure"] | None
    outcome: dict[str, object] | None
    rejection: str | None


def _receipt(content: bytes) -> EvaluationReceipt:
    return EvaluationReceipt.model_validate_json(content)


def _selected_evaluation(
    name: str = "1N",
    parameter_file: Path = PARAMETERS,
) -> tuple[EvaluationEngine, EvaluationFrame, bytes]:
    method = Method()
    selection = Selection(include=[SpinSystem.from_name(name)], exclude=None)
    session = AnalysisSession.create()
    session.set_model("2st_hd")
    experiments = build_experiments([EXPERIMENT], selection, session=session)
    session.parameters.set_defaults(read_defaults([parameter_file]))
    assert session.try_build_analysis_values()
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


def test_serial_and_fresh_spawn_replay_the_same_selected_evaluation() -> None:
    engine, frame, request = _selected_evaluation()
    expected = engine.new_evaluator().evaluate(frame)
    assert isinstance(expected, EvaluationResult)

    serial = _receipt(run_serial_request(request))
    spawned, replay = map(_receipt, run_spawned_requests((request, request)))

    assert serial.status == spawned.status == replay.status == "completed"
    assert (
        serial.outcome_kind == spawned.outcome_kind == replay.outcome_kind == "result"
    )
    assert serial.outcome == spawned.outcome == replay.outcome == expected.to_record()
    assert spawned.outcome is not None
    assert EvaluationResult.from_record(spawned.outcome, engine.plan).identity == (
        expected.identity
    )
    assert serial.worker_pid == serial.evaluator_owner_pid == os.getpid()
    assert spawned.worker_pid != replay.worker_pid
    assert {spawned.worker_pid, replay.worker_pid}.isdisjoint({os.getpid()})
    assert spawned.worker_pid == spawned.evaluator_owner_pid
    assert replay.worker_pid == replay.evaluator_owner_pid


def test_two_fresh_workers_are_independent_of_request_order() -> None:
    first_engine, _first_frame, first_request = _selected_evaluation("1N")
    second_engine, _second_frame, second_request = _selected_evaluation("2N")

    forward = tuple(
        map(_receipt, run_spawned_requests((first_request, second_request)))
    )
    reverse = tuple(
        map(_receipt, run_spawned_requests((second_request, first_request)))
    )

    assert len({receipt.worker_pid for receipt in forward}) == 2
    assert len({receipt.worker_pid for receipt in reverse}) == 2
    forward_records = {receipt.request_identity: receipt.outcome for receipt in forward}
    reverse_records = {receipt.request_identity: receipt.outcome for receipt in reverse}
    assert forward_records == reverse_records
    assert {
        receipt.outcome["plan_identity"]
        for receipt in forward
        if receipt.outcome is not None
    } == {first_engine.plan.identity, second_engine.plan.identity}


def test_interrupted_fresh_worker_returns_no_result_and_terminates() -> None:
    engine, frame, request = _selected_evaluation()
    parent_evaluator = engine.new_evaluator()
    expected = parent_evaluator.evaluate(frame)
    assert isinstance(expected, EvaluationResult)
    process, boundary = start_interruptible_spawned_request(request)
    started = boundary.poll(60.0)
    if not started:
        process.terminate()
        process.join()
    assert started
    owner_pid = int(boundary.recv_bytes())
    process.terminate()
    process.join(60.0)
    assert not process.is_alive()
    try:
        outcome = boundary.recv_bytes()
    except EOFError:
        outcome = None
    boundary.close()

    assert outcome is None
    assert process.exitcode is not None and process.exitcode != 0
    assert process.pid == owner_pid != parent_evaluator._owner_pid == os.getpid()
    assert parent_evaluator.evaluate(frame).to_record() == expected.to_record()


def test_evaluation_failure_crosses_the_spawn_boundary_unchanged() -> None:
    engine, frame, _request = _selected_evaluation()
    reversed_frame = EvaluationFrame(
        frame.parameterization_identity, tuple(reversed(frame._items))
    )
    expected = engine.new_evaluator().evaluate(reversed_frame)
    assert isinstance(expected, EvaluationFailure)
    request = evaluation_request(
        model="2st_hd",
        experiment_files=[EXPERIMENT],
        parameter_files=[PARAMETERS],
        selection=Selection(include=[SpinSystem.from_name("1N")], exclude=None),
        method=Method(),
        plan=engine.plan,
        frame=reversed_frame,
    )

    receipt = _receipt(run_spawned_requests((request,))[0])

    assert receipt.status == "completed"
    assert receipt.outcome_kind == "failure"
    assert receipt.outcome == expected.to_record()
    assert receipt.outcome is not None
    assert EvaluationFailure.from_record(receipt.outcome, engine.plan) == expected


def test_malformed_and_stale_requests_fail_closed(tmp_path: Path) -> None:
    _engine, _frame, valid_request = _selected_evaluation()
    malformed_record = json.loads(valid_request)
    malformed_record["unknown"] = True
    malformed = json.dumps(malformed_record).encode()

    copied_parameters = tmp_path / "parameters.toml"
    copied_parameters.write_bytes(PARAMETERS.read_bytes())
    _engine, _frame, stale = _selected_evaluation(parameter_file=copied_parameters)
    copied_parameters.write_text(
        copied_parameters.read_text().replace("KDH = 10.0", "KDH = 11.0")
    )

    malformed_receipt, stale_receipt = map(
        _receipt, run_spawned_requests((malformed, stale))
    )

    assert malformed_receipt.status == stale_receipt.status == "rejected"
    assert malformed_receipt.outcome is stale_receipt.outcome is None
    assert malformed_receipt.rejection is not None
    assert stale_receipt.rejection is not None
    assert "Malformed evaluation-process request" in malformed_receipt.rejection
    assert "does not match serialized plan" in stale_receipt.rejection
