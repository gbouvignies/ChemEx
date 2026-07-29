from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
TYPING_CASES = ROOT / "tests" / "typing_cases"
TY = Path(sys.executable).with_name("ty")


def _run_ty(filename: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(  # noqa: S603
        [
            str(TY),
            "check",
            "--output-format",
            "concise",
            str(TYPING_CASES / filename),
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_valid_experiment_type_family_declarations_type_check() -> None:
    result = _run_ty("valid_experiment_type_families.py")

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "filename",
    [
        "invalid_experiment_type_config.py",
        "invalid_experiment_type_family.py",
    ],
)
def test_incompatible_experiment_type_family_is_rejected(filename: str) -> None:
    result = _run_ty(filename)

    assert result.returncode != 0
    assert "invalid-argument-type" in result.stdout + result.stderr
