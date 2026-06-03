from __future__ import annotations

import os
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal

NATIVE_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
AUTO_WORKER_LIMIT = 8
ExecutionCount = int | Literal["auto"]


def available_cpu_count() -> int:
    """Return the CPU count visible to this process."""
    return max(1, os.process_cpu_count() or os.cpu_count() or 1)


def auto_worker_count() -> int:
    """Return a conservative default worker count for CLI auto mode."""
    return min(available_cpu_count(), AUTO_WORKER_LIMIT)


def _resolve_count(value: int) -> int:
    if value == 0:
        return available_cpu_count()
    if value < 0:
        msg = "Execution settings must be non-negative"
        raise ValueError(msg)
    return value


def _resolve_workers(value: ExecutionCount | None) -> int:
    if value is None:
        return auto_worker_count()
    if value == "auto":
        return auto_worker_count()
    return _resolve_count(value)


def _resolve_native_threads(value: ExecutionCount | None) -> int | None:
    if value is None or value == "auto":
        return None
    return _resolve_count(value)


def native_thread_env(
    native_threads: int | None,
    *,
    parallel: bool = False,
) -> dict[str, str]:
    if native_threads is None and not parallel:
        return {}
    thread_count = 1 if native_threads is None else native_threads
    return dict.fromkeys(NATIVE_THREAD_ENV_VARS, str(thread_count))


@contextmanager
def native_thread_environment(env: Mapping[str, str]) -> Iterator[None]:
    if not env:
        yield
        return

    previous_values = {name: os.environ.get(name) for name in env}
    os.environ.update(env)
    try:
        yield
    finally:
        for name, value in previous_values.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def _validate_resolved_count(value: int | None) -> None:
    if value is None:
        return
    if value < 1:
        msg = "Execution settings must resolve to at least one worker/thread"
        raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ExecutionSettings:
    """Runtime controls for CPU-level parallel execution."""

    workers: int = 1
    native_threads: int | None = None

    def __post_init__(self) -> None:
        _validate_resolved_count(self.workers)
        _validate_resolved_count(self.native_threads)

    @classmethod
    def from_counts(
        cls,
        *,
        workers: ExecutionCount | None = None,
        native_threads: ExecutionCount | None = None,
    ) -> ExecutionSettings:
        return cls(
            workers=_resolve_workers(workers),
            native_threads=_resolve_native_threads(native_threads),
        )

    @property
    def is_parallel(self) -> bool:
        return self.workers > 1

    def native_thread_env(self, *, parallel: bool = False) -> dict[str, str]:
        return native_thread_env(self.native_threads, parallel=parallel)
