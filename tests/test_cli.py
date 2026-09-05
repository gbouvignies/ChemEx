from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from chemex.chemex import _read_fit_methods
from chemex.cli import build_parser
from chemex.configuration.method_plan import FormatOrigin
from chemex.toml import TomlReadError


def _fit_args(method: Path) -> Namespace:
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            "experiment.toml",
            "-p",
            "parameters.toml",
            "-m",
            str(method),
        ]
    )


def test_fit_cli_warns_when_accepting_deprecated_v1_and_least_squares(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    method = tmp_path / "v1.toml"
    method.write_text('[STEP]\nFITMETHOD = "least_squares"\n', encoding="utf-8")

    plan = _read_fit_methods(_fit_args(method))

    rendered = " ".join(capsys.readouterr().out.split())
    assert plan.format_origin is FormatOrigin.V1
    assert "Method format v1 is deprecated" in rendered
    assert 'FITMETHOD and its "least_squares" alias' in rendered
    assert "Use FORMAT_VERSION = 2 for new method files" in rendered


def test_fit_cli_does_not_emit_v1_deprecation_warning_for_v2(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    method = tmp_path / "v2.toml"
    method.write_text("FORMAT_VERSION = 2\n\n[STEP]\n", encoding="utf-8")

    plan = _read_fit_methods(_fit_args(method))

    assert plan.format_origin is FormatOrigin.V2
    assert "Method format v1 is deprecated" not in capsys.readouterr().out


def test_fit_method_reader_preserves_typed_toml_error(
    tmp_path: Path,
) -> None:
    method = tmp_path / "invalid.toml"
    method.write_text("[STEP\n", encoding="utf-8")

    with pytest.raises(TomlReadError) as error_info:
        _read_fit_methods(_fit_args(method))

    assert error_info.value.filename == method
    assert error_info.value.__cause__ is not None


def test_fit_parser_defaults_execution_arguments_to_auto() -> None:
    args = build_parser().parse_args(
        [
            "fit",
            "-e",
            "experiment.toml",
            "-p",
            "parameters.toml",
        ],
    )

    assert args.workers == "auto"
    assert args.native_threads == "auto"
    assert args.analysis_command


def test_fit_parser_accepts_execution_arguments() -> None:
    args = build_parser().parse_args(
        [
            "fit",
            "-e",
            "experiment.toml",
            "-p",
            "parameters.toml",
            "--workers",
            "2",
            "--native-threads",
            "3",
        ],
    )

    assert args.workers == 2
    assert args.native_threads == 3
    assert args.analysis_command


def test_fit_parser_accepts_auto_execution_arguments() -> None:
    args = build_parser().parse_args(
        [
            "fit",
            "-e",
            "experiment.toml",
            "-p",
            "parameters.toml",
            "--workers",
            "auto",
            "--native-threads",
            "auto",
        ],
    )

    assert args.workers == "auto"
    assert args.native_threads == "auto"
    assert args.analysis_command


def test_simulate_parser_does_not_accept_execution_arguments() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                "simulate",
                "-e",
                "experiment.toml",
                "-p",
                "parameters.toml",
                "--workers",
                "2",
            ],
        )
