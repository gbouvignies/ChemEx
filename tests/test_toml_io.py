from __future__ import annotations

from pathlib import Path

import pytest

import chemex.toml as toml_module


def test_read_toml_exits_with_formatted_missing_file_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    missing_file = tmp_path / "missing.toml"
    recorded: dict[str, Path] = {}

    monkeypatch.setattr(
        toml_module,
        "print_file_not_found",
        lambda filename: recorded.update({"filename": filename}),
    )

    with pytest.raises(SystemExit):
        toml_module.read_toml(missing_file)

    assert recorded == {"filename": missing_file}


def test_read_toml_exits_with_formatted_parse_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    invalid_toml = tmp_path / "invalid.toml"
    invalid_toml.write_text("[section\nvalue = 1", encoding="utf-8")
    recorded: dict[str, object] = {}

    monkeypatch.setattr(
        toml_module,
        "print_toml_error",
        lambda filename, error: recorded.update(
            {"filename": filename, "error_type": type(error).__name__}
        ),
    )

    with pytest.raises(SystemExit):
        toml_module.read_toml(invalid_toml)

    assert recorded["filename"] == invalid_toml
    assert recorded["error_type"] == "TOMLDecodeError"
