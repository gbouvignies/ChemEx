"""One-shot process-boundary qualification helper for issue #593."""

from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
from collections.abc import Mapping, Sequence
from multiprocessing.connection import Connection
from multiprocessing.process import BaseProcess
from pathlib import Path
from typing import Literal, cast

from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

_SCHEMA_VERSION = 1


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, allow_nan=False, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")


def _identity(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _selection_value(value: object) -> list[str] | str | None:
    if value is None or value == "*":
        return cast("str | None", value)
    return [str(item) for item in cast("Sequence[SpinSystem]", value)]


def evaluation_request(
    *,
    model: str,
    experiment_files: Sequence[Path],
    parameter_files: Sequence[Path],
    selection: Selection,
    method: Method,
    plan: EvaluationPlan,
    frame: EvaluationFrame,
) -> bytes:
    payload: dict[str, object] = {
        "schema_version": _SCHEMA_VERSION,
        "model": model,
        "experiment_files": [str(path.resolve()) for path in experiment_files],
        "parameter_files": [str(path.resolve()) for path in parameter_files],
        "selection": {
            "include": _selection_value(selection.include),
            "exclude": _selection_value(selection.exclude),
        },
        "method": method.model_dump(mode="json", exclude_none=True),
        "plan": plan.to_record(),
        "frame": frame.to_record(),
    }
    payload["identity"] = _identity(payload)
    return _canonical_bytes(payload)


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"Evaluation-process {name} must be a record")
    return cast("Mapping[str, object]", value)


def _selection(value: object) -> Selection:
    record = _mapping(value, "selection")
    if set(record) != {"include", "exclude"}:
        raise ValueError("Malformed evaluation-process selection")

    def restore(item: object) -> list[SpinSystem] | Literal["*"] | None:
        if item is None or item == "*":
            return cast('Literal["*"] | None', item)
        if not isinstance(item, list) or not all(
            isinstance(name, str) for name in item
        ):
            raise TypeError("Evaluation-process selection names must be strings")
        names = cast("list[str]", item)
        return [SpinSystem.from_name(name) for name in names]

    return Selection(restore(record["include"]), restore(record["exclude"]))


def _execute_request(
    content: bytes,
    barrier: Connection | None = None,
) -> bytes:
    worker_pid = os.getpid()
    request_identity = ""
    try:
        record = json.loads(content)
        if not isinstance(record, dict) or set(record) != {
            "schema_version",
            "model",
            "experiment_files",
            "parameter_files",
            "selection",
            "method",
            "plan",
            "frame",
            "identity",
        }:
            raise ValueError("Malformed evaluation-process request")
        if not isinstance(record["identity"], str) or not isinstance(
            record["model"], str
        ):
            raise TypeError("Evaluation-process identity and model must be strings")
        request_identity = record["identity"]
        model = record["model"]
        unsigned = {key: value for key, value in record.items() if key != "identity"}
        if record["schema_version"] != _SCHEMA_VERSION or request_identity != _identity(
            unsigned
        ):
            raise ValueError(
                "Evaluation-process request identity does not match payload"
            )
        experiment_files = _absolute_paths(record["experiment_files"], "experiment")
        parameter_files = _absolute_paths(record["parameter_files"], "parameter")
        method = Method.model_validate(_mapping(record["method"], "method"))
        plan = EvaluationPlan.from_record(_mapping(record["plan"], "plan"))
        frame = EvaluationFrame.from_record(_mapping(record["frame"], "frame"))
        if frame.parameterization_identity != plan.parameterization_identity:
            raise ValueError("Evaluation frame belongs to another plan")

        session = AnalysisSession.create()
        session.set_model(model)
        experiments = build_experiments(
            experiment_files, _selection(record["selection"]), session=session
        )
        session.parameters.set_defaults(read_defaults(parameter_files))
        if not session.try_build_analysis_values():
            raise ValueError("Fresh evaluation parameter state is unavailable")
        parameterization = session.compile_parameterization(
            method, experiments.param_ids
        )
        if (
            parameterization.evaluator_identity != plan.parameterization_identity
            or parameterization.program.fingerprint != plan.constraint_program_identity
        ):
            raise ValueError("Fresh evaluation plan does not match serialized plan")
        evaluator = EvaluationEngine.bind(
            plan, parameterization, experiments
        ).new_evaluator()
        if barrier is not None:
            barrier.send_bytes(str(evaluator._owner_pid).encode("ascii"))
            barrier.recv_bytes()
        outcome = evaluator.evaluate(frame)
        kind = "result" if isinstance(outcome, EvaluationResult) else "failure"
        if not isinstance(outcome, (EvaluationResult, EvaluationFailure)):
            raise TypeError("Native evaluation returned an unknown outcome")
        receipt = {
            "schema_version": _SCHEMA_VERSION,
            "request_identity": request_identity,
            "status": "completed",
            "worker_pid": worker_pid,
            "evaluator_owner_pid": evaluator._owner_pid,
            "outcome_kind": kind,
            "outcome": outcome.to_record(),
            "rejection": None,
        }
    except Exception as error:  # noqa: BLE001 - qualification rejection fence
        receipt = {
            "schema_version": _SCHEMA_VERSION,
            "request_identity": request_identity,
            "status": "rejected",
            "worker_pid": worker_pid,
            "evaluator_owner_pid": None,
            "outcome_kind": None,
            "outcome": None,
            "rejection": f"{type(error).__name__}: {error}",
        }
    return _canonical_bytes(receipt)


def _absolute_paths(value: object, name: str) -> list[Path]:
    if not isinstance(value, list) or not value:
        raise TypeError(f"Evaluation-process {name} files must be a non-empty list")
    paths = [Path(item) for item in value if isinstance(item, str)]
    if len(paths) != len(value) or any(not path.is_absolute() for path in paths):
        raise ValueError(f"Evaluation-process {name} files must be absolute paths")
    return paths


def run_serial_request(content: bytes) -> bytes:
    return _execute_request(content)


def _spawn_worker(connection: Connection, content: bytes) -> None:
    connection.send_bytes(_execute_request(content))
    connection.close()


def run_spawned_requests(contents: Sequence[bytes]) -> tuple[bytes, ...]:
    """Start exactly one fresh spawn process per request, with no persistent pool."""
    context = multiprocessing.get_context("spawn")
    workers: list[tuple[BaseProcess, Connection]] = []
    for content in contents:
        receiving, sending = context.Pipe(duplex=False)
        process = context.Process(target=_spawn_worker, args=(sending, content))
        process.start()
        sending.close()
        workers.append((process, receiving))
    receipts: list[bytes] = []
    for process, connection in workers:
        if not connection.poll(60.0):
            process.terminate()
            process.join()
            raise TimeoutError("Spawned evaluation process did not respond")
        receipts.append(connection.recv_bytes())
        connection.close()
        process.join()
        if process.exitcode != 0:
            raise RuntimeError("Spawned evaluation process exited unsuccessfully")
    return tuple(receipts)


def start_interruptible_spawned_request(
    content: bytes,
) -> tuple[BaseProcess, Connection]:
    """Start a worker paused after local binding and before evaluation."""
    context = multiprocessing.get_context("spawn")
    barrier, child_barrier = context.Pipe()
    process = context.Process(
        target=_execute_request,
        args=(content, child_barrier),
    )
    process.start()
    child_barrier.close()
    return process, barrier
