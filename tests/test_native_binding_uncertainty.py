"""Product regressions for retained Direct-TRF Jacobians in 2st binding fits."""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import patch

import pytest

from chemex.chemex import run
from chemex.cli import build_parser
from chemex.optimize import native_deterministic as native_deterministic_module
from chemex.optimize.deterministic_uncertainty import DeterministicUncertainty
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
BINDING = ROOT / "examples/Combinations/2stBinding"
EXPERIMENTS = tuple(sorted((BINDING / "Experiments").glob("*.toml")))
PARAMETERS = BINDING / "Parameters/params.toml"
METHOD = BINDING / "Methods/method.toml"


def _arguments(
    output: Path,
    method: Path,
    parameters: tuple[Path, ...] = (PARAMETERS,),
):
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            *(str(path) for path in EXPERIMENTS),
            "-p",
            *(str(path) for path in parameters),
            "-m",
            str(method),
            "-d",
            "2st_binding",
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )


def _capture_product_uncertainty(
    output: Path,
    method: Path,
    parameters: tuple[Path, ...] = (PARAMETERS,),
) -> tuple[AnalysisSession, tuple[DeterministicUncertainty, ...]]:
    captured: list[DeterministicUncertainty] = []
    real_derive = native_deterministic_module.derive_deterministic_uncertainty

    def derive(facts):
        uncertainty = real_derive(facts)
        captured.append(uncertainty)
        return uncertainty

    session = AnalysisSession.create()
    with patch.object(
        native_deterministic_module,
        "derive_deterministic_uncertainty",
        derive,
    ):
        run(_arguments(output, method, parameters), session=session)
    return session, tuple(captured)


def test_fitted_kd_uses_backend_jacobian_before_boundary_reporting(
    tmp_path: Path,
) -> None:
    method = tmp_path / "kd-method.toml"
    method.write_text(
        """[STEP]
INCLUDE = [486, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]
FIT = ["KD"]
FIX = ["KOFF", "DW_AB", "R1_A", "R2_A", "CS_A"]
""",
        encoding="utf-8",
    )
    step1_method = tmp_path / "step1-method.toml"
    step1_method.write_text(
        """[STEP1]
INCLUDE = [486, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]
FIX = ["KD"]
""",
        encoding="utf-8",
    )
    step1_output = tmp_path / "Step1"
    step1_session = AnalysisSession.create()
    run(_arguments(step1_output, step1_method), session=step1_session)
    assert step1_session.analysis_values.snapshot().revision == 1
    output = tmp_path / "Output"

    session, (uncertainty,) = _capture_product_uncertainty(
        output,
        method,
        (PARAMETERS, step1_output / "Parameters" / "fitted.toml"),
    )

    assert session.analysis_values.snapshot().revision == 1
    evidence = uncertainty.root_evidence
    assert evidence is not None
    jacobian = evidence.residual_jacobian
    assert jacobian is not None
    assert jacobian.method == "retained-scipy-final-2-point"
    assert jacobian.evaluation_count == 0
    assert evidence.rank_diagnostic is not None
    assert evidence.covariance is not None
    kd_index = jacobian.controlled_ids.index("__KD")
    kd = evidence.accepted_anchor.vector[kd_index]
    assert kd == pytest.approx(1.0156e-6, rel=2.0e-3)
    assert math.isfinite(evidence.covariance.covariance[kd_index][kd_index])
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "KD" in fitted
    assert "Jacobian unavailable" not in fitted


def test_full_binding_step2_reports_boundary_warning_with_retained_jacobian(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"

    session, uncertainty = _capture_product_uncertainty(output, METHOD)

    assert session.analysis_values.snapshot().revision == 2
    assert len(uncertainty) == 2
    step2 = uncertainty[1].root_evidence
    assert step2 is not None
    assert step2.residual_jacobian is not None
    assert step2.residual_jacobian.method == "retained-scipy-final-2-point"
    assert step2.residual_jacobian.evaluation_count == 0
    assert step2.rank_diagnostic is not None
    assert all(failure.stage != "residual_linearization" for failure in step2.failures)
    fitted = (output / "STEP2" / "Parameters" / "fitted.toml").read_text(
        encoding="utf-8"
    )
    fitted_uncertainties = [line for line in fitted.splitlines() if "# ±" in line]
    warned_uncertainties = [
        line
        for line in fitted_uncertainties
        if "boundary may make uncertainty asymmetric" in line
    ]
    unavailable_uncertainties = [
        line for line in fitted.splitlines() if "error unavailable" in line
    ]
    assert step2.covariance is not None
    assert len(warned_uncertainties) == len(step2.covariance.simple_bound_warning_ids)
    assert 0 < len(warned_uncertainties) < len(fitted_uncertainties)
    assert unavailable_uncertainties == []
    assert "error unavailable: boundary limited" not in fitted
    assert "Jacobian unavailable" not in fitted
