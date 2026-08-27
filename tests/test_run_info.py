from __future__ import annotations

import subprocess
import tomllib
from argparse import Namespace
from datetime import UTC, datetime
from pathlib import Path

import pytest

from chemex import __version__
from chemex import run_info as run_info_module
from chemex.configuration.parameters import read_defaults
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ParameterDeclaration,
    SealedParameterDeclarations,
    SealedParameterModel,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
    extract_condition_entries,
)
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot


def _parameter_inputs(
    parameters: tuple[ParamSetting, ...] = (),
) -> dict[str, object]:
    if not parameters:
        parameters = (ParamSetting(ParamName("PB"), value=0.15, min=0.0, max=1.0),)
    ordered = tuple(sorted(parameters, key=lambda parameter: parameter.param_name))
    definitions = SealedDefinitions(
        tuple(
            ParamDefinition(
                parameter.id_,
                parameter.param_name.name,
                str(parameter.param_name.spin_system),
                extract_condition_entries(parameter.param_name.conditions),
                parameter.value,
                parameter.min,
                parameter.max,
            )
            for parameter in ordered
        ),
        {},
    )
    configuration = SealedConfiguration(
        tuple(
            ParamConfig(
                parameter.id_,
                parameter.value,
                parameter.min,
                parameter.max,
            )
            for parameter in ordered
        ),
        {},
        definitions.identity,
    )
    declarations = SealedParameterDeclarations(
        tuple(
            ParameterDeclaration(
                parameter.id_,
                parameter.vary,
                parameter.expr,
                bool(parameter.expr),
                requires_independent=not bool(parameter.expr),
                fits_by_default=parameter.vary,
            )
            for parameter in ordered
        )
    )
    parameter_model = SealedParameterModel(
        "test",
        "test-model",
        definitions,
        configuration,
        declarations,
    )
    starting_values = AnalysisValuesSnapshot(
        "test-occurrence",
        "test-model",
        definitions.identity,
        configuration.identity,
        0,
        tuple((parameter.id_, parameter.value) for parameter in ordered),
    )
    return {
        "parameter_model": parameter_model,
        "starting_values": starting_values,
    }


def _write_run_info(
    args: Namespace,
    *,
    parameters: tuple[ParamSetting, ...] = (),
    **kwargs: object,
) -> None:
    run_info_module.write_run_info(
        args,
        **_parameter_inputs(parameters),  # ty: ignore[invalid-argument-type]
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _write_input(path: Path, content: str = "[input]\nvalue = 1\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_write_run_info_captures_inputs_parameters_and_runtime(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    experiment_file = _write_input(tmp_path / "experiment.toml")
    parameters_file = _write_input(
        tmp_path / "parameters.toml",
        "[GLOBAL]\nPB = 0.1\n",
    )
    method_file = _write_input(tmp_path / "method.toml")
    output = tmp_path / "Output"

    pb = ParamSetting(
        ParamName("PB"),
        value=0.15,
        min=0.0,
        max=1.0,
    )
    dw = ParamSetting(
        ParamName("DW_AB", SpinSystem.from_name("15N")),
        value=2.5,
        min=-10.0,
        max=10.0,
        brute_step=0.5,
    )
    pa = ParamSetting(
        ParamName("PA"),
        value=0.85,
        expr=f"1.0 - {pb.id_}",
    )
    args = Namespace(
        experiments=[Path("experiment.toml")],
        parameters=[Path("parameters.toml")],
        method=[Path("method.toml")],
        output=Path("Output"),
        model="2st",
        include=[SpinSystem.from_name("15N")],
        exclude=None,
        workers=4,
        native_threads="auto",
    )
    argv = [
        "chemex",
        "fit",
        "-e",
        str(experiment_file),
        "-p",
        str(parameters_file),
    ]

    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    _write_run_info(
        args,
        parameters=(pb, dw, pa),
        argv=argv,
        working_directory=tmp_path,
        timestamp=datetime(2026, 6, 24, 12, 30, tzinfo=UTC),
    )

    run_info_path = output / "run_info"
    assert {path.relative_to(run_info_path) for path in run_info_path.rglob("*")} >= {
        Path("run.toml"),
        Path("parameters_used.toml"),
        Path("outcome.toml"),
        Path("inputs/experiments"),
        Path("inputs/parameters"),
        Path("inputs/methods"),
    }
    assert not (run_info_path / "restart.toml").exists()
    run_text = (run_info_path / "run.toml").read_text(encoding="utf-8")
    run = tomllib.loads(run_text)
    assert tomllib.loads(
        (run_info_path / "outcome.toml").read_text(encoding="utf-8")
    ) == {
        "schema_version": 2,
        "status": "running",
        "latest_committed_revision": 0,
        "restart_revision": 0,
    }
    parameters = tomllib.loads(
        (run_info_path / "parameters_used.toml").read_text(encoding="utf-8"),
    )

    assert run_text.startswith("# ChemEx run information.\n")
    assert run["schema_version"] == 2
    assert run["created_at_utc"] == "2026-06-24T12:30:00+00:00"
    assert run["run"]["kind"] == "fit"
    assert run["run"]["working_directory"] == str(tmp_path)
    assert run["run"]["requested_output_path"] == "Output"
    assert run["run"]["resolved_output_path"] == str(output)
    assert run["run"]["model"] == "2st"
    assert run["selection"] == {"include": ["15N"], "exclude": []}
    assert run["execution"] == {"workers": 4, "native_threads": "automatic"}
    assert run["software"]["chemex"] == __version__
    assert run["software"]["python"]
    assert run["software"]["numpy"]
    assert run["software"]["scipy"]
    assert run["software"]["emcee"]
    assert "platform" not in run["software"]
    assert run["command"]["arguments"] == argv
    assert "git" not in run

    copied_experiment = run["inputs"]["experiments"][0]
    assert copied_experiment["provided_path"] == "experiment.toml"
    assert copied_experiment["resolved_path"] == str(experiment_file.resolve())
    assert (run_info_path / copied_experiment["copied_path"]).exists()
    copied_parameter_path = (
        run_info_path / run["inputs"]["parameters"][0]["copied_path"]
    )
    assert copied_parameter_path.read_text(encoding="utf-8") == ("[GLOBAL]\nPB = 0.1\n")
    copied_method = run["inputs"]["methods"][0]
    assert copied_method["resolved_path"] == str(method_file.resolve())
    assert (run_info_path / copied_method["copied_path"]).exists()

    parameters_text = (run_info_path / "parameters_used.toml").read_text(
        encoding="utf-8",
    )
    assert parameters_text.startswith("# Original invocation start used by ChemEx.\n")
    assert "not the latest fitted state" in parameters_text
    assert parameters["GLOBAL"]["PB"] == [0.15, 0.0, 1.0]
    assert parameters["DW_AB"]["15N"] == [2.5, -10.0, 10.0]
    assert "PA" not in parameters["GLOBAL"]


def test_committed_restart_is_one_normal_parameter_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    pb = ParamSetting(ParamName("PB"), value=0.15, min=0.0, max=1.0)
    dw = ParamSetting(
        ParamName("DW_AB", SpinSystem.from_name("15N")),
        value=2.5,
        min=-10.0,
        max=10.0,
    )
    inputs = _parameter_inputs((pb, dw))
    starting = inputs["starting_values"]
    assert isinstance(starting, AnalysisValuesSnapshot)
    output = tmp_path / "Output"
    args = Namespace(
        experiments=[_write_input(tmp_path / "experiment.toml")],
        parameters=[_write_input(tmp_path / "parameters.toml")],
        method=None,
        output=output,
        model="2st",
        include=None,
        exclude=None,
        workers=1,
        native_threads="auto",
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    run = run_info_module.write_run_info(
        args,
        **inputs,  # ty: ignore[invalid-argument-type]
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )
    committed = AnalysisValuesSnapshot(
        starting.occurrence_identity,
        starting.model_identity,
        starting.definitions_identity,
        starting.configuration_identity,
        1,
        ((pb.id_, 0.083417), (dw.id_, 3.041)),
    )

    run.publish_restart(committed)

    restart = output / "run_info" / "restart.toml"
    parsed = tomllib.loads(restart.read_text(encoding="utf-8"))
    assert parsed["GLOBAL"]["PB"] == [0.083417, 0.0, 1.0]
    assert parsed["DW_AB"]["15N"] == [3.041, -10.0, 10.0]
    defaults = {name.id_: setting for name, setting in read_defaults([restart])}
    assert defaults[pb.id_].value == pytest.approx(0.083417)
    assert defaults[dw.id_].value == pytest.approx(3.041)
    assert run.restart_revision == 1
    assert list((output / "run_info").glob("restart.toml")) == [restart]


def test_stochastic_seed_is_a_full_unsigned_decimal_string(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    inputs = _parameter_inputs()
    output = tmp_path / "Output"
    args = Namespace(
        experiments=[_write_input(tmp_path / "experiment.toml")],
        parameters=[_write_input(tmp_path / "parameters.toml")],
        method=None,
        output=output,
        model="2st",
        include=None,
        exclude=None,
        workers=1,
        native_threads="auto",
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)
    run = run_info_module.write_run_info(
        args,
        **inputs,  # ty: ignore[invalid-argument-type]
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )

    run.record_stochastic_operation("SEARCH", "de", (1 << 64) - 1)

    record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert record["stochastic_operations"] == [
        {"step": "SEARCH", "kind": "de", "seed": "18446744073709551615"}
    ]


def test_seed_record_read_failure_is_classified_as_run_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    inputs = _parameter_inputs()
    output = tmp_path / "Output"
    args = Namespace(
        experiments=[_write_input(tmp_path / "experiment.toml")],
        parameters=[_write_input(tmp_path / "parameters.toml")],
        method=None,
        output=output,
        model="2st",
        include=None,
        exclude=None,
        workers=1,
        native_threads="auto",
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)
    run = run_info_module.write_run_info(
        args,
        **inputs,  # ty: ignore[invalid-argument-type]
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )
    run_path = output / "run_info" / "run.toml"
    original_read_text = Path.read_text

    def fail_run_record_read(path: Path, *args, **kwargs) -> str:
        if path == run_path:
            raise OSError("run record read failed")
        return original_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", fail_run_record_read)

    with pytest.raises(OSError, match="run record read failed") as caught:
        run.record_stochastic_operation("SEARCH", "de", 7)

    assert caught.value.failure_stage == "run_info"
    assert run.stochastic_operations == ()


def test_archived_tomls_are_immutable_byte_copies_alongside_outcome(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    original_bytes = {
        "experiments": b'[experiment]\r\nname = "relaxation_hznz"\r\n',
        "parameters": b"[GLOBAL]\r\nPB = 0.1\r\n",
        "methods": b'[DEFAULT]\r\nFITMETHOD = "trf"\r\n',
    }
    sources = {category: tmp_path / f"{category}.toml" for category in original_bytes}
    for category, source in sources.items():
        source.write_bytes(original_bytes[category])
    output = tmp_path / "Output"
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    _write_run_info(
        Namespace(
            experiments=[sources["experiments"]],
            parameters=[sources["parameters"]],
            method=[sources["methods"]],
            output=output,
        ),
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )

    run_info = output / "run_info"
    run_path = run_info / "run.toml"
    run_bytes = run_path.read_bytes()
    run = tomllib.loads(run_bytes.decode())
    archived = {
        category: run_info / run["inputs"][category][0]["copied_path"]
        for category in original_bytes
    }
    for category, path in archived.items():
        assert path.read_bytes() == original_bytes[category]

    for source in sources.values():
        source.write_bytes(b"[changed]\nvalue = true\n")
    run_info_module.write_run_outcome(output, "complete")

    assert run_path.read_bytes() == run_bytes
    assert tomllib.loads((run_info / "outcome.toml").read_text(encoding="utf-8")) == {
        "schema_version": 2,
        "status": "complete",
        "latest_committed_revision": 0,
        "restart_revision": 0,
    }
    for category, path in archived.items():
        assert path.read_bytes() == original_bytes[category]


def test_run_info_uses_chemex_path_semantics_for_literal_tilde(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    experiment_file = _write_input(tmp_path / "~" / "experiment.toml")
    output = tmp_path / "~" / "Output"
    args = Namespace(
        experiments=[Path("~/experiment.toml")],
        parameters=[],
        method=None,
        output=Path("~/Output"),
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    _write_run_info(
        args,
        argv=["chemex", "fit", "--output", "~/Output"],
        working_directory=tmp_path,
        timestamp=datetime(2026, 6, 24, 12, 30, tzinfo=UTC),
    )

    run_info_path = output / "run_info"
    run = tomllib.loads((run_info_path / "run.toml").read_text(encoding="utf-8"))
    copied_experiment = run["inputs"]["experiments"][0]

    assert run["run"]["resolved_output_path"] == str(output.resolve())
    assert copied_experiment["provided_path"] == "~/experiment.toml"
    assert copied_experiment["resolved_path"] == str(experiment_file.resolve())
    assert (run_info_path / copied_experiment["copied_path"]).read_text(
        encoding="utf-8",
    ) == experiment_file.read_text(encoding="utf-8")


def test_input_files_with_same_basename_are_copied_without_collision(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    experiment_a = _write_input(
        tmp_path / "field1" / "experiment.toml",
        '[input]\nfield = "a"\n',
    )
    experiment_b = _write_input(
        tmp_path / "field2" / "experiment.toml",
        '[input]\nfield = "b"\n',
    )
    output = tmp_path / "Output"
    args = Namespace(
        experiments=[
            Path("field1/experiment.toml"),
            Path("field2/experiment.toml"),
        ],
        parameters=[],
        method=None,
        output=output,
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    _write_run_info(
        args,
        argv=["chemex", "fit"],
        working_directory=tmp_path,
        timestamp=datetime(2026, 6, 24, 12, 30, tzinfo=UTC),
    )
    run_info_path = output / "run_info"
    run = tomllib.loads((run_info_path / "run.toml").read_text(encoding="utf-8"))
    copied_inputs = run["inputs"]["experiments"]
    copied_paths = [item["copied_path"] for item in copied_inputs]

    assert len(set(copied_paths)) == 2
    assert all("__" in Path(path).stem for path in copied_paths)
    assert (run_info_path / copied_paths[0]).read_text(encoding="utf-8") == (
        experiment_a.read_text(encoding="utf-8")
    )
    assert (run_info_path / copied_paths[1]).read_text(encoding="utf-8") == (
        experiment_b.read_text(encoding="utf-8")
    )

    _write_run_info(
        args,
        argv=["chemex", "fit"],
        working_directory=tmp_path,
        timestamp=datetime(2026, 6, 24, 12, 30, tzinfo=UTC),
    )
    repeated_run = tomllib.loads(
        (run_info_path / "run.toml").read_text(encoding="utf-8"),
    )
    repeated_paths = [
        item["copied_path"] for item in repeated_run["inputs"]["experiments"]
    ]
    assert repeated_paths == copied_paths


def test_git_metadata_is_optional_when_git_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(run_info_module, "_find_git_root", lambda _start: tmp_path)

    def unavailable(*_args, **_kwargs) -> subprocess.CompletedProcess[str]:
        raise FileNotFoundError

    monkeypatch.setattr(run_info_module.subprocess, "run", unavailable)

    assert run_info_module._git_metadata() is None


def test_git_metadata_is_omitted_for_repository_not_tracking_chemex(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, ...]] = []
    source_directory = Path(run_info_module.__file__).resolve().parent

    monkeypatch.setattr(
        run_info_module,
        "_find_git_root",
        lambda _start: source_directory,
    )

    def run_git(_repository: Path, *arguments: str) -> str:
        calls.append(arguments)
        if arguments[0] == "ls-files":
            raise subprocess.CalledProcessError(1, arguments)
        pytest.fail("Git metadata should stop when the source file is not tracked")

    monkeypatch.setattr(run_info_module, "_run_git", run_git)

    assert run_info_module._git_metadata() is None
    assert calls[0][:2] == ("ls-files", "--error-unmatch")


def test_same_output_rerun_can_use_archived_inputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    original_input = _write_input(tmp_path / "experiment.toml")
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    _write_run_info(
        Namespace(
            experiments=[original_input],
            parameters=[],
            method=None,
            output=output,
        ),
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )
    archived_input = (
        output / "run_info" / "inputs" / "experiments" / original_input.name
    )

    _write_run_info(
        Namespace(
            experiments=[archived_input],
            parameters=[],
            method=None,
            output=output,
        ),
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )

    copied_input = output / "run_info" / "inputs" / "experiments" / original_input.name
    assert copied_input.read_text(encoding="utf-8") == original_input.read_text(
        encoding="utf-8",
    )


def test_failed_staging_preserves_existing_run_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run_info_path = output / "run_info"
    run_info_path.mkdir(parents=True)
    existing_run = run_info_path / "run.toml"
    existing_run.write_text("schema_version = 1\n", encoding="utf-8")
    missing_input = tmp_path / "missing.toml"
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    with pytest.raises(FileNotFoundError):
        _write_run_info(
            Namespace(
                experiments=[missing_input],
                parameters=[],
                method=None,
                output=output,
            ),
            argv=["chemex", "fit"],
            working_directory=tmp_path,
        )

    assert existing_run.read_text(encoding="utf-8") == "schema_version = 1\n"
    assert not list(output.glob(".run_info-*"))


def test_failed_replacement_restores_existing_run_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    run_info_path = tmp_path / "run_info"
    run_info_path.mkdir()
    existing_run = run_info_path / "run.toml"
    existing_run.write_text("schema_version = 1\n", encoding="utf-8")
    staging_path = tmp_path / ".run_info-staging"
    staging_path.mkdir()
    (staging_path / "run.toml").write_text("schema_version = 2\n", encoding="utf-8")
    original_replace = Path.replace
    error_message = "replacement failed"

    def fail_staging_replace(path: Path, target: Path) -> Path:
        if path == staging_path:
            raise OSError(error_message)
        return original_replace(path, target)

    monkeypatch.setattr(Path, "replace", fail_staging_replace)

    with pytest.raises(OSError, match="replacement failed"):
        run_info_module._replace_run_info(
            staging_path,
            run_info_path,
        )

    assert existing_run.read_text(encoding="utf-8") == "schema_version = 1\n"
    assert staging_path.exists()
    assert not list(tmp_path.glob(".run_info-backup-*"))
