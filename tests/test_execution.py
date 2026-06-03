from __future__ import annotations

import pytest

from chemex.runtime import execution as execution_module
from chemex.runtime.execution import ExecutionSettings


def test_execution_settings_default_to_serial() -> None:
    settings = ExecutionSettings()

    assert settings.workers == 1
    assert settings.native_threads is None
    assert not settings.is_parallel
    assert settings.native_thread_env() == {}


def test_execution_settings_auto_uses_capped_cpu_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(execution_module.os, "process_cpu_count", lambda: 16)

    settings = ExecutionSettings.from_counts(workers="auto", native_threads="auto")

    assert settings.workers == 8
    assert settings.native_threads is None
    assert settings.is_parallel
    assert settings.native_thread_env() == {}
    assert settings.native_thread_env(parallel=True) == dict.fromkeys(
        execution_module.NATIVE_THREAD_ENV_VARS,
        "1",
    )


def test_execution_settings_zero_uses_available_cpu_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(execution_module.os, "process_cpu_count", lambda: 6)

    settings = ExecutionSettings.from_counts(workers=0, native_threads=0)

    assert settings.workers == 6
    assert settings.native_threads == 6


def test_execution_settings_reject_unresolved_zero() -> None:
    with pytest.raises(ValueError, match="at least one"):
        ExecutionSettings(workers=0)


def test_execution_settings_reject_negative_counts() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        ExecutionSettings.from_counts(workers=-1)

    with pytest.raises(ValueError, match="non-negative"):
        ExecutionSettings.from_counts(native_threads=-1)


def test_execution_settings_native_thread_environment() -> None:
    settings = ExecutionSettings(native_threads=2)

    assert settings.native_thread_env() == dict.fromkeys(
        execution_module.NATIVE_THREAD_ENV_VARS,
        "2",
    )


def test_native_thread_environment_restores_previous_values(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
    env = {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }

    with execution_module.native_thread_environment(env):
        assert execution_module.os.environ["OMP_NUM_THREADS"] == "1"
        assert execution_module.os.environ["OPENBLAS_NUM_THREADS"] == "1"

    assert execution_module.os.environ["OMP_NUM_THREADS"] == "4"
    assert "OPENBLAS_NUM_THREADS" not in execution_module.os.environ
