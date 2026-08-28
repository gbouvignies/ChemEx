"""Regression qualification for the small-component TRF scale policy (#710)."""

from __future__ import annotations

import tomllib
from argparse import Namespace
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

import chemex.optimize.grouped_direct_trf as grouped_direct_trf_module
import chemex.optimize.native_deterministic as native_deterministic_module
from chemex.chemex import run
from chemex.cli import build_parser
from chemex.optimize.direct_trf import (
    DirectTrfCandidateOutcome,
    DirectTrfCandidateTerminal,
)
from chemex.optimize.grouped_direct_trf import GroupedDirectTrfInvocation
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
CPMG_ROOT = ROOT / "examples/Experiments/CPMG_CH3_1H_SQ"
CEST_ROOT = ROOT / "examples/Experiments/CEST_13C_LABEL_CN"


def _fit_arguments(example: Path, output: Path) -> Namespace:
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            *(str(path) for path in sorted((example / "Experiments").glob("*.toml"))),
            "-p",
            str(example / "Parameters/parameters.toml"),
            "-m",
            str(example / "Methods/method.toml"),
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
            "--native-threads",
            "1",
        ]
    )


def _run_and_capture_components(
    example: Path,
    output: Path,
) -> tuple[
    AnalysisSession,
    tuple[DirectTrfCandidateOutcome, ...],
    tuple[GroupedDirectTrfInvocation, ...],
]:
    original = grouped_direct_trf_module.execute_direct_trf_candidate
    original_build_invocation = native_deterministic_module._build_invocation
    outcomes: list[DirectTrfCandidateOutcome] = []
    invocations: list[GroupedDirectTrfInvocation] = []

    def capture(*args: Any, **kwargs: Any) -> DirectTrfCandidateOutcome:
        outcome = original(*args, **kwargs)
        outcomes.append(outcome)
        return outcome

    def capture_invocation(*args: Any, **kwargs: Any) -> GroupedDirectTrfInvocation:
        invocation = original_build_invocation(*args, **kwargs)
        invocations.append(invocation)
        return invocation

    session = AnalysisSession.create()
    with (
        patch.object(
            grouped_direct_trf_module,
            "execute_direct_trf_candidate",
            side_effect=capture,
        ),
        patch.object(
            native_deterministic_module,
            "_build_invocation",
            side_effect=capture_invocation,
        ),
    ):
        run(_fit_arguments(example, output), session=session)
    return session, tuple(outcomes), tuple(invocations)


def _cpmg_ids(residue: str) -> tuple[str, ...]:
    return (
        f"__DW_AB_{residue}",
        f"__R2A_A_{residue}_600_2MHZ",
        f"__R2A_A_{residue}_800_2MHZ",
        f"__R2_A_{residue}_600_2MHZ",
        f"__R2_A_{residue}_800_2MHZ",
    )


def _component_map(
    outcomes: tuple[DirectTrfCandidateOutcome, ...],
) -> dict[tuple[str, ...], DirectTrfCandidateOutcome]:
    return {
        outcome.execution.backend.final_residual_jacobian.controlled_ids: outcome
        for outcome in outcomes
        if outcome.execution.backend is not None
    }


def _invocation_scale_map(
    outcomes: tuple[DirectTrfCandidateOutcome, ...],
    invocation: GroupedDirectTrfInvocation,
) -> dict[tuple[str, ...], tuple[float, ...]]:
    assert len(outcomes) == len(invocation.component_invocations)
    return {
        outcome.execution.backend.final_residual_jacobian.controlled_ids: item.x_scale
        for outcome, item in zip(
            outcomes,
            invocation.component_invocations,
            strict=True,
        )
        if outcome.execution.backend is not None
    }


def test_complete_cpmg_ch3_1h_sq_workflow_converges_exact_step3_components(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session, outcomes, invocations = _run_and_capture_components(CPMG_ROOT, output)

    assert tuple(len(item.component_invocations) for item in invocations) == (10, 1, 16)
    assert len(outcomes) == 27
    step3 = outcomes[11:]
    expected = tuple(
        _cpmg_ids(residue)
        for residue in (
            "I43HD1",
            "I44HD1",
            "L24HD1",
            "L24HD2",
            "L25HD1",
            "L25HD2",
            "L52HD1",
            "L52HD2",
            "L55HD1",
            "L55HD2",
            "M40HE",
            "M42HE",
            "V30HG1",
            "V30HG2",
            "V67HG1",
            "V67HG2",
        )
    )
    assert (
        tuple(
            outcome.execution.backend.final_residual_jacobian.controlled_ids
            for outcome in step3
            if outcome.execution.backend is not None
        )
        == expected
    )
    assert all(
        outcome.terminal is DirectTrfCandidateTerminal.SUCCESS for outcome in step3
    )

    components = _component_map(step3)
    scales = _invocation_scale_map(step3, invocations[2])
    # These sub-micro chi-square tolerances cover binary64 finite-difference and
    # backend variation while remaining many orders below the 12k-budget gaps.
    qualified = {
        "I43HD1": (
            200,
            2.8012650851,
            (0.0534523022, 43.5239201, 39.0819921, 34.6039587, 39.2064460),
        ),
        "M42HE": (
            200,
            1.9527249983,
            (0.0352002711, 37.6989759, 34.8058070, 31.4773856, 34.5806268),
        ),
        "V67HG1": (
            200,
            1.7300765470,
            (0.0259395306, 39.3128688, 36.2392556, 33.3733778, 36.3260882),
        ),
        "V67HG2": (
            200,
            2.0364471661,
            (0.0450755787, 39.1482135, 34.9783369, 30.8321366, 35.0750046),
        ),
    }
    for residue, (request_limit, chi_square, vector) in qualified.items():
        outcome = components[_cpmg_ids(residue)]
        assert scales[_cpmg_ids(residue)] == (1.0, 5.0, 5.0, 5.0, 5.0)
        assert outcome.execution.counters.objective_requests_accepted < request_limit
        assert outcome.candidate is not None
        assert outcome.candidate.chi_square == pytest.approx(
            chi_square, rel=0.0, abs=5.0e-7
        )
        # Ill-conditioned relaxation pairs permit small coordinate drift; these
        # bounds still reject a distinct basin or a scientifically changed DW.
        assert outcome.candidate.vector[0] == pytest.approx(
            vector[0], rel=0.0, abs=5.0e-5
        )
        assert outcome.candidate.vector[1:] == pytest.approx(
            vector[1:], rel=0.0, abs=5.0e-3
        )
        assert outcome.execution.backend is not None
        assert outcome.execution.backend.active_mask == (0, 0, 0, 0, 0)

    assert session.analysis_values.snapshot().revision == 3
    statistics = tomllib.loads(
        (output / "STEP3/statistics.toml").read_text(encoding="utf-8")
    )
    # statistics.toml prints six significant figures (half-unit rounding here).
    assert statistics["chi-square"] == pytest.approx(72.1332, rel=0.0, abs=5.0e-5)


def test_cest_13c_qualification_preserves_basin_with_bounded_small_scale(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session, outcomes, invocations = _run_and_capture_components(CEST_ROOT, output)

    assert tuple(len(item.component_invocations) for item in invocations) == (1, 15)
    assert len(outcomes) == 16
    step2 = outcomes[1:]
    assert all(
        outcome.terminal is DirectTrfCandidateTerminal.SUCCESS for outcome in step2
    )
    assert (
        sum(outcome.execution.counters.objective_requests_accepted for outcome in step2)
        < 4_000
    )

    l18cd1_ids = (
        "__CS_A_L18CD1",
        "__DW_AB_L18CD1",
        "__R1_A_L18CD1_598_8MHZ",
        "__R2_A_L18CD1_598_8MHZ",
    )
    l18cd1 = _component_map(step2)[l18cd1_ids]
    assert _invocation_scale_map(step2, invocations[1])[l18cd1_ids] == (
        5.0,
        1.0,
        3.5,
        5.0,
    )
    assert l18cd1.execution.counters.objective_requests_accepted < 3_000
    assert l18cd1.candidate is not None
    # This is tight enough to distinguish the qualified lower basin while
    # accommodating binary64 finite-difference variation across supported hosts.
    assert l18cd1.candidate.chi_square == pytest.approx(
        801.9787003, rel=0.0, abs=2.0e-6
    )
    l18cd1_vector = (24.92877435, 0.19320262, 4.84303720, 6.17970065)
    assert l18cd1.candidate.vector[0] == pytest.approx(
        l18cd1_vector[0], rel=0.0, abs=2.0e-4
    )
    assert l18cd1.candidate.vector[1] == pytest.approx(
        l18cd1_vector[1], rel=0.0, abs=2.0e-5
    )
    assert l18cd1.candidate.vector[2:] == pytest.approx(
        l18cd1_vector[2:], rel=0.0, abs=5.0e-3
    )

    assert session.analysis_values.snapshot().revision == 2
    statistics = tomllib.loads(
        (output / "STEP2/statistics.toml").read_text(encoding="utf-8")
    )
    # statistics.toml prints six significant figures (half-unit rounding here).
    assert statistics["chi-square"] == pytest.approx(2849.45, rel=0.0, abs=5.0e-3)
