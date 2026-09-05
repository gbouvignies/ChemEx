from __future__ import annotations

from pathlib import Path

import pytest

import chemex.configuration.parameters as parameter_configuration_module
import chemex.toml as toml_module
from chemex.configuration.parameters import ParameterConfigurationError, read_defaults


def test_read_toml_raises_typed_missing_file_error(tmp_path: Path) -> None:
    missing_file = tmp_path / "missing.toml"

    with pytest.raises(toml_module.TomlReadError) as error_info:
        toml_module.read_toml(missing_file)

    assert error_info.value.filename == missing_file
    assert isinstance(error_info.value.__cause__, FileNotFoundError)


def test_read_toml_raises_typed_parse_error(tmp_path: Path) -> None:
    invalid_toml = tmp_path / "invalid.toml"
    invalid_toml.write_text("[section\nvalue = 1", encoding="utf-8")
    with pytest.raises(toml_module.TomlReadError) as error_info:
        toml_module.read_toml(invalid_toml)

    assert error_info.value.filename == invalid_toml
    assert type(error_info.value.__cause__).__name__ == "TOMLDecodeError"


def test_read_toml_raises_typed_invalid_utf8_error(tmp_path: Path) -> None:
    invalid_toml = tmp_path / "invalid-utf8.toml"
    invalid_toml.write_bytes(b"[section]\nvalue = \xff\n")

    with pytest.raises(toml_module.TomlReadError) as error_info:
        toml_module.read_toml(invalid_toml)

    assert error_info.value.filename == invalid_toml
    assert isinstance(error_info.value.__cause__, UnicodeDecodeError)


def test_internal_toml_loader_type_error_is_not_misclassified(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cause = TypeError("loader contract defect")

    def fail_loader(_filename: Path) -> None:
        raise cause

    monkeypatch.setattr(toml_module, "load_toml", fail_loader)

    with pytest.raises(TypeError, match="loader contract defect") as error_info:
        toml_module.read_toml(tmp_path / "parameters.toml")

    assert error_info.value is cause


def test_read_toml_contextualizes_permission_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    filename = tmp_path / "private.toml"
    cause = PermissionError(13, "Permission denied", str(filename))

    def deny_read(_filename: Path) -> None:
        raise cause

    monkeypatch.setattr(toml_module, "load_toml", deny_read)

    with pytest.raises(toml_module.TomlReadError) as error_info:
        toml_module.read_toml(filename)

    assert error_info.value.filename == filename
    assert error_info.value.__cause__ is cause


def test_parameter_configuration_error_retains_source_and_cause(tmp_path: Path) -> None:
    parameter_file = tmp_path / "parameters.toml"
    parameter_file.write_text('"[global]" = "not-a-table"\n', encoding="utf-8")

    with pytest.raises(ParameterConfigurationError) as error_info:
        read_defaults([parameter_file])

    assert error_info.value.filenames == (parameter_file,)
    assert error_info.value.__cause__ is not None


def test_empty_parameter_values_are_a_typed_configuration_failure(
    tmp_path: Path,
) -> None:
    parameter_file = tmp_path / "parameters.toml"
    parameter_file.write_text("[GLOBAL]\nPB = []\n", encoding="utf-8")

    with pytest.raises(ParameterConfigurationError) as error_info:
        read_defaults([parameter_file])

    assert error_info.value.filenames == (parameter_file,)
    assert error_info.value.__cause__ is not None
    assert "at least 1" in error_info.value.explanation


def test_unexpected_parameter_builder_failure_is_not_misclassified(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    parameter_file = tmp_path / "parameters.toml"
    parameter_file.write_text("[GLOBAL]\nPB = 0.1\n", encoding="utf-8")

    def fail_builder(_config: object) -> None:
        raise TypeError("builder defect")

    monkeypatch.setattr(
        parameter_configuration_module,
        "build_default_list",
        fail_builder,
    )

    with pytest.raises(TypeError, match="builder defect") as error_info:
        parameter_configuration_module.read_defaults([parameter_file])

    assert not isinstance(error_info.value, ParameterConfigurationError)
