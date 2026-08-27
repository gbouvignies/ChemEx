"""Behavioral tests for session-owned revisioned analysis values."""

from __future__ import annotations

import dataclasses
from argparse import Namespace
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from typing import cast

import pytest

from chemex import chemex as chemex_module
from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.parameters.sealed import ParamConfig, SealedConfiguration
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import (
    AnalysisValues,
    AnalysisValuesCommitError,
    AnalysisValuesSnapshot,
    IncompatibleAnalysisValuesError,
    InvalidAnalysisValuesCommitError,
    StaleAnalysisValuesError,
)
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
SHIPPED_EXPERIMENT = (
    ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
)
SHIPPED_PARAMETERS = (
    ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
)
SHIPPED_METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
CPMG_EXPERIMENT = ROOT / "examples/Experiments/CPMG_15N_IP/Experiments/500mhz.toml"
CPMG_PARAMETERS = ROOT / "examples/Experiments/CPMG_15N_IP/Parameters/parameters.toml"


def _configuration() -> SealedConfiguration:
    configs = (
        ParamConfig("__PB", 0.08, 0.0, 1.0),
        ParamConfig("__KEX_AB", 400.0, 0.0, 10_000.0),
    )
    return SealedConfiguration(
        _configs=configs,
        _index={config.param_id: index for index, config in enumerate(configs)},
        definitions_identity="definitions-v1",
    )


def _build_shipped_session(
    session: AnalysisSession | None = None,
) -> AnalysisSession:
    session = AnalysisSession.create() if session is None else session
    session.set_model("2st")
    build_experiments(
        [SHIPPED_EXPERIMENT],
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([SHIPPED_PARAMETERS]))
    assert session.try_build_analysis_values()
    return session


def test_revision_zero_model_values_and_constraints_resolve_without_lmfit() -> None:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [CPMG_EXPERIMENT],
        Selection(include=[SpinSystem.from_name("15N-HN")], exclude=None),
        session=session,
    )

    session.parameters.set_defaults(read_defaults([CPMG_PARAMETERS]))
    assert session.try_build_analysis_values()
    initial = session.resolve_current_values(experiments.param_ids)

    parameterization = session.compile_parameterization(
        Method(constraints=["[PB] = [KEX_AB] / 10000.0"]),
        experiments.param_ids,
    )
    constrained = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )

    assert initial["__PB"] == pytest.approx(0.07)
    assert initial["__PA"] == pytest.approx(0.93)
    assert initial["__KAB"] == pytest.approx(28.0)
    assert initial["__KBA"] == pytest.approx(372.0)
    assert constrained["__PB"] == pytest.approx(0.04)
    assert constrained["__PA"] == pytest.approx(0.96)
    assert constrained["__KAB"] == pytest.approx(16.0)
    assert constrained["__KBA"] == pytest.approx(384.0)


def _shipped_args(command: str, output: Path) -> Namespace:
    args = Namespace(
        commands=command,
        model="2st",
        include=[SpinSystem.from_name("G2N-HN")],
        exclude=None,
        experiments=[SHIPPED_EXPERIMENT],
        parameters=[SHIPPED_PARAMETERS],
        method=[SHIPPED_METHOD] if command == "fit" else None,
        output=output,
        plot="nothing",
    )
    if command == "fit":
        args.workers = 1
        args.native_threads = 1
    return args


def _authoritative_output_contents(output: Path) -> dict[Path, bytes]:
    return {
        path.relative_to(output): path.read_bytes()
        for path in output.rglob("*")
        if path.is_file()
        and (
            path.relative_to(output).parts[0] in {"Data", "Parameters"}
            or path.relative_to(output) == Path("statistics.toml")
        )
    }


def test_values_initialize_from_configuration_in_deterministic_order() -> None:
    values = AnalysisValues()

    values.initialize("2st", _configuration())
    snapshot = values.snapshot()

    assert snapshot.revision == 0
    assert snapshot.model_identity == "2st"
    assert snapshot.definitions_identity == "definitions-v1"
    assert snapshot.configuration_identity == _configuration().identity
    assert tuple(snapshot.items()) == (("__PB", 0.08), ("__KEX_AB", 400.0))


def test_snapshot_is_immutable_and_has_a_deterministic_json_round_trip() -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    snapshot = values.snapshot()

    encoded = snapshot.to_json()
    restored = AnalysisValuesSnapshot.from_json(encoded)

    assert restored == snapshot
    assert restored.occurrence_identity == snapshot.occurrence_identity
    assert restored.to_json() == encoded
    assert '"occurrence_identity"' in encoded
    assert '"revision":0' in encoded
    assert '"value":"0x1.47ae147ae147bp-4"' in encoded
    with pytest.raises(dataclasses.FrozenInstanceError):
        snapshot.revision = 1  # ty: ignore[invalid-assignment]
    with pytest.raises(TypeError):
        snapshot["__PB"] = 0.2  # ty: ignore[invalid-assignment]


def test_commit_updates_one_complete_scope_and_increments_one_global_revision() -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    initial = values.snapshot()
    restored_initial = AnalysisValuesSnapshot.from_json(initial.to_json())

    committed = values.commit(
        {"__PB": 0.12},
        expected=restored_initial,
        scope=("__PB",),
    )

    assert committed.revision == 1
    assert tuple(committed.items()) == (("__PB", 0.12), ("__KEX_AB", 400.0))
    assert initial.revision == 0
    assert initial["__PB"] == 0.08


def test_stale_and_incompatible_commits_are_atomic() -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    initial = values.snapshot()
    values.commit({"__PB": 0.12}, expected=initial, scope=("__PB",))
    committed = values.snapshot()

    with pytest.raises(StaleAnalysisValuesError):
        values.commit({"__KEX_AB": 500.0}, expected=initial, scope=("__KEX_AB",))

    incompatible = dataclasses.replace(committed, configuration_identity="other")
    with pytest.raises(IncompatibleAnalysisValuesError):
        values.commit({"__KEX_AB": 500.0}, expected=incompatible, scope=("__KEX_AB",))

    assert values.snapshot() == committed


def test_concurrent_disjoint_commits_use_one_conservative_global_revision() -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    initial = values.snapshot()
    barrier = Barrier(2)

    def commit(candidate: dict[str, float]) -> AnalysisValuesSnapshot | Exception:
        barrier.wait()
        try:
            return values.commit(
                candidate,
                expected=initial,
                scope=tuple(candidate),
            )
        except AnalysisValuesCommitError as error:
            return error

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = tuple(executor.map(commit, ({"__PB": 0.12}, {"__KEX_AB": 500.0})))

    assert sum(isinstance(result, AnalysisValuesSnapshot) for result in results) == 1
    assert sum(isinstance(result, StaleAnalysisValuesError) for result in results) == 1
    assert values.snapshot().revision == 1


@pytest.mark.parametrize(
    ("candidate", "scope"),
    [
        ({}, ("__PB",)),
        ({"__PB": 0.12, "__KEX_AB": 500.0}, ("__PB",)),
        ({"__UNKNOWN": 1.0}, ("__UNKNOWN",)),
        ({"__PB": "not-a-number"}, ("__PB",)),
        ({"__PB": float("nan")}, ("__PB",)),
        ({"__PB": float("inf")}, ("__PB",)),
    ],
)
def test_invalid_or_incomplete_commits_leave_values_and_revision_unchanged(
    candidate: Mapping[str, object],
    scope: tuple[str, ...],
) -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    initial = values.snapshot()

    with pytest.raises(InvalidAnalysisValuesCommitError):
        values.commit(
            cast("Mapping[str, float]", candidate),
            expected=initial,
            scope=scope,
        )

    assert values.snapshot() == initial


def test_commit_accepts_finite_resolved_values_outside_configured_bounds() -> None:
    values = AnalysisValues()
    values.initialize("2st", _configuration())
    initial = values.snapshot()

    committed = values.commit(
        {"__PB": 1.01, "__KEX_AB": -1.0},
        expected=initial,
        scope=("__PB", "__KEX_AB"),
    )

    assert committed.revision == 1
    assert committed["__PB"] == 1.01
    assert committed["__KEX_AB"] == -1.0


def test_real_shipped_configuration_initializes_session_values() -> None:
    session = _build_shipped_session()
    snapshot = session.analysis_values.snapshot()
    parameter_model = session.parameter_factory.sealed_parameter_model

    assert snapshot.revision == 0
    assert parameter_model is not None
    assert tuple(snapshot) == tuple(
        config.param_id for config in parameter_model.configuration
    )
    assert all(
        snapshot[config.param_id] == config.effective_value
        for config in parameter_model.configuration
    )


def test_session_reset_rebuilds_values_at_revision_zero() -> None:
    session = _build_shipped_session()
    initial = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {"__PB": 0.12},
        expected=initial,
        scope=("__PB",),
    )
    assert session.analysis_values.snapshot().revision == 1

    session.reset()

    with pytest.raises(RuntimeError, match="not initialized"):
        session.analysis_values.snapshot()

    rebuilt = _build_shipped_session(session).analysis_values.snapshot()
    assert rebuilt.revision == 0
    assert rebuilt.occurrence_identity != initial.occurrence_identity
    assert rebuilt.model_identity == initial.model_identity
    assert rebuilt.definitions_identity == initial.definitions_identity
    assert rebuilt.configuration_identity == initial.configuration_identity
    assert rebuilt["__PB"] == initial["__PB"]


def test_snapshot_from_before_session_reset_is_rejected_after_rebuild() -> None:
    session = _build_shipped_session()
    previous_occurrence = session.analysis_values.snapshot()

    session.reset()
    rebuilt = _build_shipped_session(session).analysis_values.snapshot()

    with pytest.raises(IncompatibleAnalysisValuesError):
        session.analysis_values.commit(
            {"__PB": 0.12},
            expected=previous_occurrence,
            scope=("__PB",),
        )

    assert session.analysis_values.snapshot() == rebuilt


def test_snapshot_from_another_identical_session_is_rejected() -> None:
    session_a = _build_shipped_session()
    session_b = _build_shipped_session()
    snapshot_a = session_a.analysis_values.snapshot()
    snapshot_b = session_b.analysis_values.snapshot()

    assert snapshot_a.occurrence_identity != snapshot_b.occurrence_identity
    assert snapshot_a.model_identity == snapshot_b.model_identity
    assert snapshot_a.definitions_identity == snapshot_b.definitions_identity
    assert snapshot_a.configuration_identity == snapshot_b.configuration_identity

    with pytest.raises(IncompatibleAnalysisValuesError):
        session_b.analysis_values.commit(
            {"__PB": 0.12},
            expected=snapshot_a,
            scope=("__PB",),
        )

    assert session_b.analysis_values.snapshot() == snapshot_b


def test_model_change_is_rejected_after_analysis_values_are_initialized() -> None:
    session = _build_shipped_session()
    snapshot = session.analysis_values.snapshot()

    with pytest.raises(RuntimeError, match="reset the analysis session"):
        session.set_model("3st")

    assert session.model.name == "2st"
    assert session.analysis_values.snapshot() == snapshot


def test_shipped_fit_commits_native_authoritative_values(
    tmp_path: Path,
) -> None:
    candidate_output = tmp_path / "candidate-enabled"
    session = AnalysisSession.create()

    chemex_module.run(_shipped_args("fit", candidate_output), session=session)

    snapshot = session.analysis_values.snapshot()
    assert snapshot.revision == 1
    assert list((candidate_output / "Data").rglob("*.dat"))
    assert list((candidate_output / "Parameters").rglob("*.toml"))


def test_shipped_simulation_fails_closed_on_native_initialization_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidate_output = tmp_path / "candidate-enabled"
    failed_output = tmp_path / "failed"
    session = AnalysisSession.create()

    chemex_module.run(_shipped_args("simulate", candidate_output), session=session)

    snapshot = session.analysis_values.snapshot()
    assert snapshot.revision == 0

    def fail_native_initialization(*_args: object, **_kwargs: object) -> None:
        msg = "native initialization failed"
        raise RuntimeError(msg)

    monkeypatch.setattr(AnalysisValues, "initialize", fail_native_initialization)

    failed_session = AnalysisSession.create()
    with pytest.raises(RuntimeError, match="native initialization failed"):
        chemex_module.run(
            _shipped_args("simulate", failed_output),
            session=failed_session,
        )

    assert not failed_output.exists()
    assert isinstance(
        failed_session.parameter_factory.native_construction_error,
        RuntimeError,
    )
    with pytest.raises(RuntimeError, match="not initialized"):
        failed_session.analysis_values.snapshot()
