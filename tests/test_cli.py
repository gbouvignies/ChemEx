from __future__ import annotations

import pytest

from chemex.cli import build_parser


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
