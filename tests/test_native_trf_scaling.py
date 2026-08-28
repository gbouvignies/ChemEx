"""Regression qualification for adaptive Jacobian TRF scaling (#710)."""

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
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.optimize.direct_trf import (
    DirectTrfCandidateOutcome,
    DirectTrfCandidateTerminal,
    DirectTrfScalePolicy,
    OptimizationProblem,
    canonical_chi_square,
)
from chemex.optimize.grouped_direct_trf import GroupedDirectTrfInvocation
from chemex.parameters.parameterization import ActiveParameterization
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
    dict[
        tuple[str, ...],
        tuple[OptimizationProblem, ActiveParameterization, EvaluationEngine],
    ],
]:
    original = grouped_direct_trf_module.execute_direct_trf_candidate
    original_build_invocation = native_deterministic_module._build_invocation
    outcomes: list[DirectTrfCandidateOutcome] = []
    invocations: list[GroupedDirectTrfInvocation] = []
    contexts: dict[
        tuple[str, ...],
        tuple[OptimizationProblem, ActiveParameterization, EvaluationEngine],
    ] = {}

    def capture(
        problem: OptimizationProblem,
        invocation: Any,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
        *args: Any,
        **kwargs: Any,
    ) -> DirectTrfCandidateOutcome:
        outcome = original(
            problem,
            invocation,
            parameterization,
            engine,
            *args,
            **kwargs,
        )
        outcomes.append(outcome)
        contexts[problem.controlled_ids] = (problem, parameterization, engine)
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
    return session, tuple(outcomes), tuple(invocations), contexts


def _chi_square_at(
    context: tuple[OptimizationProblem, ActiveParameterization, EvaluationEngine],
    vector: tuple[float, ...],
) -> float:
    problem, parameterization, engine = context
    lifecycle = problem.lifecycle_frame(vector, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    residuals = engine.new_evaluator().evaluate_residuals(frame)
    assert not isinstance(residuals, EvaluationFailure)
    return canonical_chi_square(residuals)


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


def _invocation_policy_map(
    outcomes: tuple[DirectTrfCandidateOutcome, ...],
    invocation: GroupedDirectTrfInvocation,
) -> dict[tuple[str, ...], DirectTrfScalePolicy]:
    assert len(outcomes) == len(invocation.component_invocations)
    return {
        outcome.execution.backend.final_residual_jacobian.controlled_ids: (
            item.scale_policy
        )
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
    session, outcomes, invocations, _contexts = _run_and_capture_components(
        CPMG_ROOT, output
    )

    assert tuple(len(item.component_invocations) for item in invocations) == (10, 1, 16)
    assert {
        component.scale_policy
        for invocation in invocations
        for component in invocation.component_invocations
    } == {DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM}
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
    policies = _invocation_policy_map(step3, invocations[2])
    # These basin-scale tolerances cover supported-host finite-difference and
    # linear-algebra variation while remaining many orders below the 12k-budget
    # gaps. The fitted-vector checks below provide the scientific fingerprint.
    qualified = {
        "I43HD1": (
            200,
            2.8012626983,
            (0.0534525168, 43.5239093, 39.0820438, 34.6039736, 39.2063936),
        ),
        "M42HE": (
            200,
            1.9527243308,
            (0.0352001568, 37.6990053, 34.8055903, 31.4773625, 34.5808355),
        ),
        "V67HG1": (
            200,
            1.7300761633,
            (0.0259394799, 39.3129002, 36.2391310, 33.3733510, 36.3262094),
        ),
        "V67HG2": (
            200,
            2.0364462925,
            (0.0450756842, 39.1482012, 34.9781873, 30.8321486, 35.0751506),
        ),
    }
    for residue, (request_limit, chi_square, vector) in qualified.items():
        outcome = components[_cpmg_ids(residue)]
        assert policies[_cpmg_ids(residue)] is (
            DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM
        )
        assert outcome.execution.counters.objective_requests_accepted < request_limit
        assert outcome.candidate is not None
        assert outcome.candidate.chi_square == pytest.approx(
            chi_square, rel=0.0, abs=5.0e-4
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
    assert statistics["chi-square"] == pytest.approx(72.1331, rel=0.0, abs=5.0e-4)


def test_cest_13c_jacobian_scaled_local_refinement_converges_qualified_basin(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session, outcomes, invocations, contexts = _run_and_capture_components(
        CEST_ROOT, output
    )

    assert tuple(len(item.component_invocations) for item in invocations) == (1, 15)
    assert {
        component.scale_policy
        for invocation in invocations
        for component in invocation.component_invocations
    } == {DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM}
    assert len(outcomes) == 16
    step2 = outcomes[1:]
    assert all(
        outcome.terminal is DirectTrfCandidateTerminal.SUCCESS for outcome in step2
    )
    assert (
        sum(outcome.execution.counters.objective_requests_accepted for outcome in step2)
        < 1_000
    )

    l18cd1_ids = (
        "__CS_A_L18CD1",
        "__DW_AB_L18CD1",
        "__R1_A_L18CD1_598_8MHZ",
        "__R2_A_L18CD1_598_8MHZ",
    )
    l18cd1 = _component_map(step2)[l18cd1_ids]
    assert _invocation_policy_map(step2, invocations[1])[l18cd1_ids] is (
        DirectTrfScalePolicy.ADAPTIVE_INVERSE_JACOBIAN_COLUMN_NORM
    )
    assert l18cd1.execution.counters.objective_requests_accepted < 200
    assert l18cd1.candidate is not None
    # This protects the reproducible operational basin from the shipped start;
    # it does not claim that local TRF has discovered the global minimum.
    assert l18cd1.candidate.chi_square == pytest.approx(1093.42227, rel=0.0, abs=2.0e-3)
    l18cd1_vector = (24.97239884, -0.14243669, 4.83981904, 7.18204366)
    assert l18cd1.candidate.vector[0] == pytest.approx(
        l18cd1_vector[0], rel=0.0, abs=2.0e-4
    )
    assert l18cd1.candidate.vector[1] == pytest.approx(
        l18cd1_vector[1], rel=0.0, abs=2.0e-4
    )
    assert l18cd1.candidate.vector[2:] == pytest.approx(
        l18cd1_vector[2:], rel=0.0, abs=5.0e-3
    )

    lower_basin_vector = (24.92877435, 0.19320262, 4.84303720, 6.17970065)
    assert _chi_square_at(contexts[l18cd1_ids], lower_basin_vector) == pytest.approx(
        801.9787003, rel=0.0, abs=1.0e-3
    )
    assert _chi_square_at(contexts[l18cd1_ids], l18cd1.candidate.vector) == (
        pytest.approx(l18cd1.candidate.chi_square, rel=0.0, abs=1.0e-10)
    )

    assert session.analysis_values.snapshot().revision == 2
    statistics = tomllib.loads(
        (output / "STEP2/statistics.toml").read_text(encoding="utf-8")
    )
    # statistics.toml prints six significant figures (half-unit rounding here).
    assert statistics["chi-square"] == pytest.approx(3140.89, rel=0.0, abs=5.0e-3)
