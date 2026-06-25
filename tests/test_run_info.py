from __future__ import annotations

import subprocess
import tomllib
from argparse import Namespace
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemex import __version__
from chemex import chemex as chemex_module
from chemex import run_info as run_info_module
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem


class ParameterStore:
    def __init__(self, parameters: dict[str, ParamSetting]) -> None:
        self.parameters = parameters

    def get_parameters(self, param_ids: object) -> dict[str, ParamSetting]:
        return {
            param_id: parameter
            for param_id, parameter in self.parameters.items()
            if param_id in set(param_ids)
        }


class FitExperiments:
    def __init__(self, parameter_store: ParameterStore, param_ids: set[str]) -> None:
        self.parameter_store = parameter_store
        self.param_ids = param_ids
        self.filter_calls = 0

    def filter(self) -> None:
        self.filter_calls += 1


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
    store = ParameterStore({pb.id_: pb, dw.id_: dw, pa.id_: pa})
    experiments = SimpleNamespace(
        param_ids={pb.id_, dw.id_, pa.id_},
        parameter_store=store,
    )
    args = Namespace(
        experiments=[Path("experiment.toml")],
        parameters=[Path("parameters.toml")],
        method=[Path("method.toml")],
        output=output,
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

    run_info_module.write_run_info(
        args,
        experiments,
        argv=argv,
        working_directory=tmp_path,
        timestamp=datetime(2026, 6, 24, 12, 30, tzinfo=UTC),
    )

    run_info_path = output / "run_info"
    run_text = (run_info_path / "run.toml").read_text(encoding="utf-8")
    run = tomllib.loads(run_text)
    parameters = tomllib.loads(
        (run_info_path / "parameters_used.toml").read_text(encoding="utf-8"),
    )

    assert run_text.startswith("# ChemEx run information.\n")
    assert run["schema_version"] == 1
    assert run["created_at_utc"] == "2026-06-24T12:30:00+00:00"
    assert run["run"]["kind"] == "fit"
    assert run["run"]["working_directory"] == str(tmp_path)
    assert run["run"]["output_directory"] == str(output)
    assert run["chemex"]["version"] == __version__
    assert run["python"]["version"]
    assert run["python"]["platform"]
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
    assert parameters_text.startswith(
        "# Starting independent parameters used by ChemEx for this fit.\n"
    )
    assert "intended to be usable as a starting parameter file" in parameters_text
    assert parameters["GLOBAL"]["PB"] == [0.15, 0.0, 1.0]
    assert parameters["DW_AB"]["15N"] == [2.5, -10.0, 10.0, 0.5]
    assert "PA" not in parameters["GLOBAL"]


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
    experiments = SimpleNamespace(
        param_ids=set(),
        parameter_store=ParameterStore({}),
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    run_info_module.write_run_info(
        args,
        experiments,
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

    run_info_module.write_run_info(
        args,
        experiments,
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


def test_run_fit_creates_run_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    experiment_file = _write_input(tmp_path / "experiment.toml")
    parameters_file = _write_input(tmp_path / "parameters.toml")
    parameter = ParamSetting(ParamName("PB"), value=0.2)
    store = ParameterStore({parameter.id_: parameter})
    experiments = FitExperiments(store, {parameter.id_})
    output = tmp_path / "Output"
    args = Namespace(
        experiments=[experiment_file],
        parameters=[parameters_file],
        method=None,
        output=output,
        plot="nothing",
    )
    run_methods_calls: list[Path] = []

    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)
    monkeypatch.setattr(chemex_module, "print_start_fit", lambda: None)
    monkeypatch.setattr(
        chemex_module,
        "run_methods",
        lambda _experiments, _methods, path, _plot, **_kwargs: run_methods_calls.append(
            path
        ),
    )

    chemex_module.run_fit(
        args,
        experiments,
        SimpleNamespace(),
        argv=["chemex", "fit"],
    )

    assert experiments.filter_calls == 1
    assert run_methods_calls == [output]
    assert (output / "run_info" / "run.toml").exists()
    assert (output / "run_info" / "parameters_used.toml").exists()


def test_git_metadata_is_optional_when_git_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(run_info_module, "_find_git_root", lambda _start: tmp_path)

    def unavailable(*_args, **_kwargs) -> subprocess.CompletedProcess[str]:
        raise FileNotFoundError

    monkeypatch.setattr(run_info_module.subprocess, "run", unavailable)

    assert run_info_module._git_metadata() is None  # noqa: SLF001


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

    assert run_info_module._git_metadata() is None  # noqa: SLF001
    assert calls[0][:2] == ("ls-files", "--error-unmatch")


def test_same_output_rerun_can_use_archived_inputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    original_input = _write_input(tmp_path / "experiment.toml")
    experiments = SimpleNamespace(
        param_ids=set(),
        parameter_store=ParameterStore({}),
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    run_info_module.write_run_info(
        Namespace(
            experiments=[original_input],
            parameters=[],
            method=None,
            output=output,
        ),
        experiments,
        argv=["chemex", "fit"],
        working_directory=tmp_path,
    )
    archived_input = (
        output / "run_info" / "inputs" / "experiments" / original_input.name
    )

    run_info_module.write_run_info(
        Namespace(
            experiments=[archived_input],
            parameters=[],
            method=None,
            output=output,
        ),
        experiments,
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
    experiments = SimpleNamespace(
        param_ids=set(),
        parameter_store=ParameterStore({}),
    )
    monkeypatch.setattr(run_info_module, "_git_metadata", lambda: None)

    with pytest.raises(FileNotFoundError):
        run_info_module.write_run_info(
            Namespace(
                experiments=[missing_input],
                parameters=[],
                method=None,
                output=output,
            ),
            experiments,
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
        run_info_module._replace_run_info(  # noqa: SLF001
            staging_path,
            run_info_path,
        )

    assert existing_run.read_text(encoding="utf-8") == "schema_version = 1\n"
    assert staging_path.exists()
    assert not list(tmp_path.glob(".run_info-backup-*"))
