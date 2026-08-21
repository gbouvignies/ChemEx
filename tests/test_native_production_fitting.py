from __future__ import annotations

import json
import math
import subprocess
import sys
import tomllib
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

import chemex.optimize.direct_trf as direct_trf_module
import chemex.optimize.fitting as fitting_module
import chemex.optimize.mcmc as mcmc_module
import chemex.optimize.native_deterministic as native_deterministic_module
import chemex.optimize.native_mcmc as native_mcmc_module
import chemex.optimize.native_resampling as native_resampling_module
import chemex.optimize.resampling as resampling_module
import chemex.optimize.uncertainty as uncertainty_module
import chemex.run_info as run_info_module
from chemex.chemex import run
from chemex.cli import build_parser
from chemex.optimize.mcmc import NativeMcmcIncompleteError
from chemex.optimize.progress import ProgressPhase
from chemex.optimize.resampling import NativeResamplingIncompleteError
from chemex.optimize.uncertainty import ParameterUnit
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
EXAMPLE = ROOT / "examples/Experiments/RELAXATION_HZNZ"
EXPERIMENT = EXAMPLE / "Experiments/800mhz.toml"
PARAMETERS = EXAMPLE / "Parameters/parameters.toml"
METHOD = EXAMPLE / "Methods/method.toml"
DCEST_EXAMPLE = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH"
DCEST_EXPERIMENT = DCEST_EXAMPLE / "Experiments/3hz.toml"
DCEST_PARAMETERS = DCEST_EXAMPLE / "Parameters/parameters.toml"


def _fit_arguments(
    output: Path,
    method: Path = METHOD,
    *,
    parameters: Path = PARAMETERS,
    include: tuple[str, ...] = ("G2N-HN",),
    plot_level: str = "nothing",
    workers: int = 1,
):
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            str(EXPERIMENT),
            "-p",
            str(parameters),
            "-m",
            str(method),
            "-o",
            str(output),
            "--include",
            *include,
            "--plot",
            plot_level,
            "--workers",
            str(workers),
        ]
    )


def _simulation_arguments(output: Path):
    return build_parser().parse_args(
        [
            "simulate",
            "-e",
            str(EXPERIMENT),
            "-p",
            str(PARAMETERS),
            "-o",
            str(output),
            "--include",
            "G2N-HN",
            "--plot",
            "nothing",
        ]
    )


def test_product_trf_uses_legacy_request_ceiling_and_physical_coordinate_scale() -> (
    None
):
    problem = type(
        "Problem",
        (),
        {
            "controlled_ids": ("a", "b", "c", "d"),
            "start": (0.01, -2.0, 0.0, 150.0),
        },
    )()

    assert (
        native_deterministic_module._objective_request_budget(
            problem  # ty: ignore[invalid-argument-type]
        )
        == 10000
    )
    assert native_deterministic_module._product_x_scale(
        problem  # ty: ignore[invalid-argument-type]
    ) == (1.0, 2.0, 1.0, 150.0)


def test_real_simulation_uses_native_values_and_preserves_back_calculation(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    run(_simulation_arguments(output), session=session)

    data_path = output / "Data" / "800mhz.dat"
    calculated = np.loadtxt(data_path, comments="#", skiprows=2, usecols=1)
    legacy_calculated = np.array(
        [
            1.18750357e06,
            9.66031999e05,
            7.85819221e05,
            6.39171094e05,
            5.19835478e05,
            4.22730273e05,
            3.43723059e05,
        ]
    )
    np.testing.assert_allclose(calculated, legacy_calculated, rtol=1.0e-6)
    assert (output / "Parameters" / "fixed.toml").is_file()
    assert (output / "Parameters" / "constrained.toml").is_file()


def _run_real_fit_cli(output: Path, method: Path, parameters: Path) -> None:
    subprocess.run(  # noqa: S603 - fixed local CLI and repository-owned fixtures
        [
            sys.executable,
            "-m",
            "chemex",
            "fit",
            "-e",
            str(EXPERIMENT),
            "-p",
            str(parameters),
            "-m",
            str(method),
            "-o",
            str(output),
            "--include",
            "G2N-HN",
            "--plot",
            "nothing",
            "--workers",
            "1",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )


def _read_outcome(output: Path) -> dict[str, object]:
    return tomllib.loads(
        (output / "run_info" / "outcome.toml").read_text(encoding="utf-8")
    )


def _assert_truthful_resampling_summary(path: Path) -> None:
    summary = tomllib.loads(path.read_text(encoding="utf-8"))
    distribution = next(iter(summary.values()))
    assert {
        "percentile_68_lower",
        "percentile_68_upper",
        "half_percentile_68_width",
    } <= distribution.keys()
    assert distribution["percentile_68_lower"] <= distribution["median"]
    assert distribution["percentile_68_upper"] >= distribution["median"]
    assert distribution["half_percentile_68_width"] >= 0.0
    assert "stderr" not in distribution


def test_real_direct_fit_uses_native_trf_and_commits_product_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    run(_fit_arguments(output, plot_level="normal"), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    fitted_record = next(line for line in fitted.splitlines() if "=" in line)
    fitted_value = float(fitted_record.split("=", 1)[1].split()[0])
    assert fitted_value == pytest.approx(2.34742, rel=5.0e-6)
    assert "# ±" in fitted_record
    fitted_error = float(fitted_record.split("±", 1)[1])
    assert fitted_error == pytest.approx(0.083366086, rel=3.0e-6)
    assert all(
        parameter.stderr is None
        for parameter in session.parameters.get_parameters(
            session.analysis_values.snapshot()
        ).values()
    )
    covariance = output / "Statistics" / "Covariance" / "evidence.json"
    assert covariance.is_file()
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    propagated_record = next(
        line for line in constrained.splitlines() if "G2N-H" in line
    )
    propagated_error = float(propagated_record.split("±", 1)[1].split()[0])
    assert propagated_error == pytest.approx(fitted_error, rel=1.0e-12)
    assert (output / "Statistics" / "Constrained" / "evidence.json").is_file()
    assert (output / "statistics.toml").is_file()
    data = next((output / "Data").glob("*.dat")).read_text(encoding="utf-8")
    first_calculation = float(
        next(line for line in data.splitlines() if line.lstrip()[:1].isdigit()).split()[
            -1
        ]
    )
    # One intensity unit is 1e-4 of the source experimental uncertainty.
    assert first_calculation == pytest.approx(9.36430201e5, abs=1.0)
    fit_curve = (output / "Plots" / "800mhz.fit").read_text(encoding="utf-8")
    curve = [
        tuple(float(value) for value in line.split())
        for line in fit_curve.splitlines()
        if line.lstrip()[:1].isdigit()
    ]
    time_1, intensity_1 = curve[0]
    time_2, intensity_2 = curve[-1]
    fitted_curve_rate = -math.log(intensity_2 / intensity_1) / (time_2 - time_1)
    assert fitted_curve_rate == pytest.approx(fitted_value, rel=1.0e-3)


def test_relaxation_product_covariance_matches_absolute_sigma_reference(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    spin_systems = ("G2N-HN", "H3N-HN", "K4N-HN", "S5N-HN", "L6N-HN")

    run(
        _fit_arguments(output, include=spin_systems),
        session=AnalysisSession.create(),
    )

    evidence = json.loads(
        (output / "All" / "Statistics" / "Covariance" / "evidence.json").read_text(
            encoding="utf-8"
        )
    )
    assert (
        evidence["covariance"]["residual_variance_scaling"]
        == "absolute_observation_uncertainties"
    )
    entries = evidence["marginal_standard_errors"]["entries"]
    observed = {
        entry["param_id"]: float.fromhex(entry["value"]["binary64"])
        for entry in entries
    }
    expected = {
        "G2": 0.083366086,
        "H3": 0.084798698,
        "K4": 0.085572683,
        "S5": 0.086949663,
        "L6": 0.086999315,
    }
    for residue, reference in expected.items():
        param_id = next(
            param_id for param_id in observed if f"_{residue}N_H_" in param_id
        )
        assert observed[param_id] == pytest.approx(reference, rel=3.0e-6)

    historical = {
        "G2": (0.123732, 2.20286),
        "H3": (0.110470, 1.69712),
        "K4": (0.0727646, 0.72305),
        "S5": (0.0607433, 0.48805),
        "L6": (0.101158, 1.35197),
    }
    for residue, reference in expected.items():
        lmfit_error, reduced_chi_square = historical[residue]
        assert reference * math.sqrt(reduced_chi_square) == pytest.approx(
            lmfit_error,
            rel=5.0e-6,
        )


def test_real_direct_fit_reports_objective_evaluation_progress(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    run(
        _fit_arguments(tmp_path / "Output"),
        session=AnalysisSession.create(),
    )

    rendered = capsys.readouterr().out
    assert "eval 0/" in rendered
    assert "final χ²" in rendered
    assert "iteration" not in rendered.lower()


def test_real_direct_fit_propagates_progress_interrupt_after_cleanup(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    original_observe = native_deterministic_module.MinimizationProgressReporter.observe

    def interrupt_after_evaluation(self, context, event):
        original_observe(self, context, event)
        if event.phase is ProgressPhase.EVALUATED:
            raise KeyboardInterrupt

    with (
        patch.object(
            native_deterministic_module.MinimizationProgressReporter,
            "observe",
            interrupt_after_evaluation,
        ),
        pytest.raises(KeyboardInterrupt),
    ):
        run(
            _fit_arguments(tmp_path / "Output"),
            session=AnalysisSession.create(),
        )

    rendered = capsys.readouterr().out
    assert "interrupted" in rendered
    assert "eval 1 total" in rendered


def test_explicit_trf_and_least_squares_alias_have_identical_native_product_output(
    tmp_path: Path,
) -> None:
    outputs: dict[str, Path] = {}
    for spelling in ("trf", "least_squares"):
        method = tmp_path / f"{spelling}.toml"
        method.write_text(
            f'[DEFAULT]\nFITMETHOD = "{spelling}"\nFIX = ["PB", "KEX_AB"]\n',
            encoding="utf-8",
        )
        output = tmp_path / spelling
        run(
            _fit_arguments(output, method),
            session=AnalysisSession.create(),
        )
        outputs[spelling] = output

    for relative in (
        Path("Parameters/fitted.toml"),
        Path("Data/800mhz.dat"),
        Path("statistics.toml"),
    ):
        assert (outputs["least_squares"] / relative).read_bytes() == (
            outputs["trf"] / relative
        ).read_bytes()


@pytest.mark.parametrize("fitmethod", ("leastsq", "differential_evolution"))
def test_invalid_fitmethod_fails_before_defaults_and_output_invalidation(
    tmp_path: Path,
    fitmethod: str,
) -> None:
    output = tmp_path / "Output"
    output.mkdir()
    sentinel = output / "Data"
    sentinel.mkdir()
    preserved = sentinel / "preserved.dat"
    preserved.write_text("existing result\n", encoding="utf-8")
    method = tmp_path / "method.toml"
    method.write_text(
        f'[DEFAULT]\nFITMETHOD = "{fitmethod}"\n',
        encoding="utf-8",
    )

    with (
        patch(
            "chemex.parameters.database.ParameterStore.set_defaults",
            side_effect=AssertionError("parameter defaults were resolved"),
        ),
        pytest.raises(SystemExit) as raised,
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    assert raised.value.code == 1
    assert preserved.read_text(encoding="utf-8") == "existing result\n"
    assert not (output / "run_info").exists()


def test_parameters_used_supports_real_cli_restart_round_trip(tmp_path: Path) -> None:
    first_output = tmp_path / "First"
    restarted_output = tmp_path / "Restarted"

    _run_real_fit_cli(first_output, METHOD, PARAMETERS)
    restart_parameters = first_output / "run_info" / "parameters_used.toml"
    _run_real_fit_cli(restarted_output, METHOD, restart_parameters)

    def fitted_value(output: Path) -> float:
        fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
        record = next(line for line in fitted.splitlines() if "=" in line)
        return float(record.split("=", 1)[1].split()[0])

    first_value = fitted_value(first_output)
    restarted_value = fitted_value(restarted_output)
    # fitted.toml is rendered to 5 decimal places in scientific notation; this
    # allows roughly one final rendered digit of solver/platform variation.
    assert first_value == pytest.approx(2.34742, rel=5.0e-6)
    assert restarted_value == pytest.approx(first_value, rel=5.0e-6)


def test_failed_deterministic_rerun_invalidates_prior_results_and_is_incomplete(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run(
        _fit_arguments(output, plot_level="normal"),
        session=AnalysisSession.create(),
    )
    user_file = output / "notes.txt"
    user_file.write_text("keep me\n", encoding="utf-8")

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=RuntimeError("rerun backend failed"),
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(output), session=AnalysisSession.create())

    outcome = _read_outcome(output)
    assert outcome["schema_version"] == 1
    assert outcome["status"] == "incomplete"
    assert outcome["terminal"] == "failed"
    assert outcome["failure_type"] == "RuntimeError"
    assert not (output / "Parameters").exists()
    assert not (output / "Data").exists()
    assert not (output / "Plots").exists()
    assert not (output / "statistics.toml").exists()
    assert user_file.read_text(encoding="utf-8") == "keep me\n"


def test_data_writer_failure_after_commit_cannot_complete_or_retain_stale_results(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run(_fit_arguments(output), session=AnalysisSession.create())
    stale_parameter = output / "Parameters" / "stale.toml"
    stale_data = output / "Data" / "stale.dat"
    stale_parameter.write_text("old parameter\n", encoding="utf-8")
    stale_data.write_text("old data\n", encoding="utf-8")
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.containers.experiments.Experiments.write",
            side_effect=RuntimeError("data writer failed"),
        ),
        pytest.raises(RuntimeError, match="data writer failed"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert _read_outcome(output) == {
        "schema_version": 1,
        "status": "incomplete",
        "terminal": "failed",
        "failure_type": "RuntimeError",
        "failure_message": "data writer failed",
    }
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert not stale_parameter.exists()
    assert not stale_data.exists()
    assert not (output / "statistics.toml").exists()


def test_final_complete_outcome_write_failure_is_terminal_and_propagates(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    real_atomic_write = run_info_module.write_text_atomic

    def fail_complete_outcome(destination: Path, content: str) -> None:
        if 'status = "complete"' in content:
            raise OSError("complete outcome publication failed")
        real_atomic_write(destination, content)

    session = AnalysisSession.create()
    with (
        patch.object(
            run_info_module,
            "write_text_atomic",
            side_effect=fail_complete_outcome,
        ),
        pytest.raises(OSError, match="complete outcome publication failed"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert _read_outcome(output) == {
        "schema_version": 1,
        "status": "incomplete",
        "terminal": "failed",
        "failure_type": "OSError",
        "failure_message": "complete outcome publication failed",
    }


def test_real_grouped_direct_fit_uses_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    run(
        _fit_arguments(output, include=("G2N-HN", "H3N-HN")),
        session=session,
    )

    assert session.analysis_values.snapshot().revision == 1
    group_outputs = tuple((output / "Groups").glob("*/Parameters/fitted.toml"))
    assert len(group_outputs) == 2
    assert (output / "All" / "Parameters" / "fitted.toml").is_file()


def test_grouped_covariance_isolates_a_rank_deficient_independent_block(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    original_svd = uncertainty_module.svd
    covariance_svd_call = 0

    def one_bad_block_svd(*args, **kwargs):
        nonlocal covariance_svd_call
        covariance_svd_call += 1
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        if covariance_svd_call in {1, 3}:
            singular[-1] = 0.0
        return left, singular, right

    with patch("chemex.optimize.uncertainty.svd", side_effect=one_bad_block_svd):
        run(
            _fit_arguments(output, include=("G2N-HN", "H3N-HN")),
            session=AnalysisSession.create(),
        )

    group_outputs = sorted((output / "Groups").glob("*/Parameters/fitted.toml"))
    assert len(group_outputs) == 2
    assert (output / "All" / "Statistics" / "Covariance" / "evidence.json").is_file()
    assert all(
        not (
            path.parent.parent / "Statistics" / "Covariance" / "evidence.json"
        ).exists()
        for path in group_outputs
    )
    reports = tuple(path.read_text(encoding="utf-8") for path in group_outputs)
    assert sum("# ±" in report for report in reports) == 1
    assert sum("error unavailable: rank deficient" in report for report in reports) == 1
    constrained_reports = tuple(
        (path.parent / "constrained.toml").read_text(encoding="utf-8")
        for path in group_outputs
    )
    assert tuple("# ±" in report for report in reports) == tuple(
        "# ±" in report for report in constrained_reports
    )
    assert (
        sum(
            "error unavailable: constrained propagation unavailable" in report
            for report in constrained_reports
        )
        == 1
    )
    blocks = json.loads(
        (output / "All" / "Statistics" / "Covariance" / "blocks.json").read_text(
            encoding="utf-8"
        )
    )
    assert len(blocks["partition_proof_identity"]) == 64
    assert sorted(block["unavailable_reason"] for block in blocks["blocks"]) == [
        "",
        "rank deficient",
    ]


def test_interrupted_block_fallback_preserves_committed_root_evidence(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    original_svd = uncertainty_module.svd

    def rank_deficient_svd(*args, **kwargs):
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        singular[-1] = 0.0
        return left, singular, right

    with (
        patch("chemex.optimize.uncertainty.svd", side_effect=rank_deficient_svd),
        patch.object(
            native_deterministic_module,
            "derive_root_anchored_block_covariance",
            side_effect=KeyboardInterrupt,
        ),
    ):
        run(
            _fit_arguments(output, include=("G2N-HN", "H3N-HN")),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    covariance_path = output / "All" / "Statistics" / "Covariance"
    assert (covariance_path / "evidence.json").is_file()
    assert not (covariance_path / "blocks.json").exists()
    status = json.loads((covariance_path / "status.json").read_text(encoding="utf-8"))
    assert status["status"] == "incomplete"
    assert status["terminal"] == "interrupted"
    assert status["reason"] == "block derivation interrupted"
    assert (output / "All" / "Parameters" / "fitted.toml").is_file()


def test_shared_parameter_coupling_is_not_split_into_covariance_blocks(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIT = ["PB"]
FIX = ["KEX_AB"]
""",
        encoding="utf-8",
    )
    original_svd = uncertainty_module.svd

    def rank_deficient_svd(*args, **kwargs):
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        singular[-1] = 0.0
        return left, singular, right

    with patch("chemex.optimize.uncertainty.svd", side_effect=rank_deficient_svd):
        run(
            _fit_arguments(output, method, include=("G2N-HN", "H3N-HN")),
            session=AnalysisSession.create(),
        )

    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "# ±" not in fitted
    assert "error unavailable: rank deficient" in fitted
    assert not (output / "Statistics" / "Covariance" / "blocks.json").exists()


def test_unsupported_scientific_constraint_keeps_product_fitted_errors(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIT = ["D2O"]
FIX = ["CS_A", "DW_AB", "KDH", "PHI", "R1_A", "R2_A", "R2_B"]
""",
        encoding="utf-8",
    )
    arguments = build_parser().parse_args(
        [
            "fit",
            "-e",
            str(DCEST_EXPERIMENT),
            "-p",
            str(DCEST_PARAMETERS),
            "-m",
            str(method),
            "-o",
            str(output),
            "-d",
            "2st_hd",
            "--include",
            "1N",
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )

    run(arguments, session=AnalysisSession.create())

    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "# ±" in fitted
    assert "error unavailable: unsupported constrained derivative" in constrained
    constrained_evidence = (
        output / "Statistics" / "Constrained" / "evidence.json"
    ).read_text(encoding="utf-8")
    assert "unsupported constrained derivative" not in constrained_evidence


def test_real_grouped_direct_progress_has_one_bounded_noninteractive_stream(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    run(
        _fit_arguments(
            tmp_path / "Output",
            include=("G2N-HN", "H3N-HN"),
        ),
        session=AnalysisSession.create(),
    )

    rendered = capsys.readouterr().out
    assert "group 1/2" in rendered
    assert rendered.count("eval 0/") == 1
    assert "eval " in rendered
    assert " total" in rendered


def test_boundary_limited_fit_keeps_central_values_and_withholds_symmetric_errors(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIX = ["KEX_AB"]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "error unavailable: boundary limited" in fitted
    assert "# ±" not in fitted
    evidence = (output / "Statistics" / "Covariance" / "evidence.json").read_text(
        encoding="utf-8"
    )
    assert '"name": "BOUNDARY_SEPARATION"' in evidence
    assert '"state": "violated"' in evidence
    assert all(
        parameter.stderr is None
        for parameter in session.parameters.get_parameters(
            session.analysis_values.snapshot()
        ).values()
    )


def test_rank_deficient_covariance_is_nonfatal_and_replaces_stale_valid_errors(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run(_fit_arguments(output), session=AnalysisSession.create())
    assert "# ±" in (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")

    original_svd = uncertainty_module.svd

    def rank_deficient_svd(*args, **kwargs):
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        singular[-1] = 0.0
        return left, singular, right

    session = AnalysisSession.create()
    with patch("chemex.optimize.uncertainty.svd", side_effect=rank_deficient_svd):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "error unavailable: rank deficient" in fitted
    assert "# ±" not in fitted
    evidence = (output / "Statistics" / "Covariance" / "evidence.json").read_text(
        encoding="utf-8"
    )
    assert '"category": "rank_deficient"' in evidence
    assert '"covariance": null' in evidence


def test_interrupted_covariance_keeps_committed_fit_and_invalidates_stale_evidence(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run(_fit_arguments(output), session=AnalysisSession.create())
    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()

    session = AnalysisSession.create()
    with patch(
        "chemex.optimize.method_step.derive_uncertainty_evidence",
        side_effect=KeyboardInterrupt,
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "error unavailable: derivation interrupted" in fitted
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "error unavailable: derivation interrupted" in constrained
    assert "# ±" not in constrained
    covariance_path = output / "Statistics" / "Covariance"
    assert not (covariance_path / "evidence.json").exists()
    status = json.loads((covariance_path / "status.json").read_text(encoding="utf-8"))
    assert status == {
        "artifact_type": "native_covariance_derivation_status",
        "reason": "derivation interrupted",
        "schema_version": 1,
        "status": "incomplete",
        "terminal": "interrupted",
    }


def _grid_method(path: Path) -> Path:
    path.write_text(
        """[GRID]
FIX = ["PB", "KEX_AB"]
GRID = ["[R1A_A] = (1.0, 3.0)"]
""",
        encoding="utf-8",
    )
    return path


def _statistics_method(path: Path, statistics: str) -> Path:
    path.write_text(
        f"""[DEFAULT]
FIX = ["PB", "KEX_AB"]
STATISTICS = {{ {statistics} }}
""",
        encoding="utf-8",
    )
    return path


def _bounded_parameters(path: Path) -> Path:
    path.write_text(
        PARAMETERS.read_text(encoding="utf-8")
        + """

[R1A_A]
G2N-HN = [1.5, 0.1, 5.0]
""",
        encoding="utf-8",
    )
    return path


def _nested_mcmc_method(path: Path, *, workers: int) -> Path:
    path.write_text(
        f"""[DEFAULT]
FIX = ["PB", "KEX_AB"]

[DEFAULT.STATISTICS.MCMC]
STEPS = 5
BURN = 1
THIN = 2
WALKERS = 4
SEED = 612
WORKERS = {workers}
""",
        encoding="utf-8",
    )
    return path


def test_real_compact_mcmc_fit_is_wholly_native_and_writes_products(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MCMC" = 2')
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()

    with patch.object(
        mcmc_module,
        "execute_mcmc_evidence",
        wraps=mcmc_module.execute_mcmc_evidence,
    ) as native_sampler:
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MCMC"
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "samples.tsv").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["engine"] == "native MCMC"
    assert diagnostics["steps"] == 2
    assert diagnostics["walkers"] == 32
    assert "lmfit_version" not in diagnostics
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "# ±" in fitted
    assert "half_credible_interval_68_width" not in fitted
    plan = native_sampler.call_args.args[1]
    assert plan.coordinate_units[0][1] is ParameterUnit.UNSPECIFIED


def test_compact_mcmc_form_runs_through_real_chemex_cli(tmp_path: Path) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MCMC" = 2')
    parameters = _bounded_parameters(tmp_path / "parameters.toml")

    _run_real_fit_cli(output, method, parameters)

    diagnostics = tomllib.loads(
        (output / "Statistics" / "MCMC" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["engine"] == "native MCMC"
    assert diagnostics["steps"] == 2


def test_seeded_nested_mcmc_settings_run_through_real_chemex_cli(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _nested_mcmc_method(tmp_path / "method.toml", workers=1)
    parameters = _bounded_parameters(tmp_path / "parameters.toml")

    _run_real_fit_cli(output, method, parameters)

    diagnostics = tomllib.loads(
        (output / "Statistics" / "MCMC" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["root_seed"] == 612
    assert diagnostics["requested_burn"] == 1
    assert diagnostics["thin"] == 2
    assert diagnostics["walkers"] == 4
    assert diagnostics["workers"] == 1


def test_seeded_nested_native_mcmc_replays_across_worker_counts(
    tmp_path: Path,
) -> None:
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    serial_method = _nested_mcmc_method(tmp_path / "serial.toml", workers=1)
    parallel_method = _nested_mcmc_method(tmp_path / "parallel.toml", workers=2)
    serial = tmp_path / "Serial"
    parallel = tmp_path / "Parallel"

    run(
        _fit_arguments(serial, serial_method, parameters=parameters),
        session=AnalysisSession.create(),
    )
    run(
        _fit_arguments(parallel, parallel_method, parameters=parameters),
        session=AnalysisSession.create(),
    )

    serial_samples = (serial / "Statistics" / "MCMC" / "samples.tsv").read_text(
        encoding="utf-8"
    )
    parallel_samples = (parallel / "Statistics" / "MCMC" / "samples.tsv").read_text(
        encoding="utf-8"
    )
    assert serial_samples == parallel_samples
    diagnostics = tomllib.loads(
        (parallel / "Statistics" / "MCMC" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["root_seed"] == 612
    assert diagnostics["workers"] == 2
    assert diagnostics["requested_burn"] == 1
    assert diagnostics["thin"] == 2
    assert diagnostics["retained_steps"] == 2
    assert diagnostics["retained_samples"] == 8
    summary = tomllib.loads(
        (parallel / "Statistics" / "MCMC" / "summary.toml").read_text(encoding="utf-8")
    )
    posterior = next(iter(summary.values()))
    assert posterior["prior_lower"] == pytest.approx(0.1)
    assert posterior["prior_upper"] == pytest.approx(5.0)
    assert 0.1 < posterior["median"] < 5.0
    assert posterior["standard_deviation"] > 0.0
    assert posterior["credible_interval_68_lower"] <= posterior["median"]
    assert posterior["credible_interval_68_upper"] >= posterior["median"]
    assert posterior["half_credible_interval_68_width"] >= 0.0
    assert "stderr" not in posterior


def test_native_mcmc_preserves_committed_central_fit(tmp_path: Path) -> None:
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    method = _nested_mcmc_method(tmp_path / "method.toml", workers=1)
    central_session = AnalysisSession.create()
    mcmc_session = AnalysisSession.create()

    run(
        _fit_arguments(tmp_path / "Central", parameters=parameters),
        session=central_session,
    )
    run(
        _fit_arguments(tmp_path / "Mcmc", method, parameters=parameters),
        session=mcmc_session,
    )

    central = central_session.analysis_values.snapshot()
    sampled = mcmc_session.analysis_values.snapshot()
    assert central.revision == sampled.revision == 1
    assert central.items() == sampled.items()


def test_unbounded_native_mcmc_fails_closed_after_committing_central_fit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MCMC" = 2')
    session = AnalysisSession.create()

    with pytest.raises(ValueError, match="finite lower and upper bounds"):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["terminal"] == "failed"
    assert "finite lower and upper bounds" in diagnostics["failure_message"]
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


def test_interrupted_native_mcmc_publishes_only_incomplete_diagnostics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _nested_mcmc_method(tmp_path / "method.toml", workers=1)
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()

    with (
        patch.object(
            native_mcmc_module._RecordingStretchMove,
            "propose",
            side_effect=KeyboardInterrupt,
        ),
        pytest.raises(NativeMcmcIncompleteError, match="interrupted"),
    ):
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["terminal"] == "interrupted"
    assert diagnostics["completed_steps"] == 0
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()
    assert _read_outcome(output)["terminal"] == "interrupted"


def test_native_mcmc_postprocessing_failure_replaces_running_diagnostics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _nested_mcmc_method(tmp_path / "method.toml", workers=1)
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()

    with (
        patch.object(
            mcmc_module,
            "_native_result_from_evidence",
            side_effect=RuntimeError("missing acceptance diagnostics"),
        ),
        pytest.raises(RuntimeError, match="acceptance diagnostics"),
    ):
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["terminal"] == "failed"
    assert diagnostics["completed_steps"] == 5
    assert "acceptance diagnostics" in diagnostics["failure_message"]
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


def test_mixed_mc_and_mcmc_use_one_native_central_fit_without_legacy_fallback(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(
        tmp_path / "method.toml",
        '"MC" = 1, "MCMC" = {STEPS = 2, BURN = 0, WALKERS = 4, SEED = 612}',
    )
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()

    with patch.object(
        fitting_module,
        "run_native_deterministic",
        wraps=fitting_module.run_native_deterministic,
    ) as native_central:
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=session,
        )

    native_central.assert_called_once()
    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Statistics" / "MonteCarlo" / "summary.toml").is_file()
    assert (output / "Statistics" / "MCMC" / "summary.toml").is_file()


def test_failed_mixed_statistics_rerun_removes_later_family_outputs_eagerly(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(
        tmp_path / "method.toml",
        '"MC" = 1, "BS" = 1, "MCMC" = {STEPS = 2, BURN = 0, WALKERS = 4, SEED = 612}',
    )
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    run(
        _fit_arguments(output, method, parameters=parameters),
        session=AnalysisSession.create(),
    )
    assert _read_outcome(output) == {"schema_version": 1, "status": "complete"}
    assert (output / "Statistics" / "Bootstrap" / "diagnostics.toml").is_file()
    assert (output / "Statistics" / "MCMC" / "diagnostics.toml").is_file()

    with (
        patch(
            "chemex.optimize.resampling._native_dataset",
            side_effect=RuntimeError("MC setup failed"),
        ),
        pytest.raises(RuntimeError, match="MC setup failed"),
    ):
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=AnalysisSession.create(),
        )

    monte_carlo = output / "Statistics" / "MonteCarlo"
    assert (
        tomllib.loads((monte_carlo / "diagnostics.toml").read_text(encoding="utf-8"))[
            "status"
        ]
        == "incomplete"
    )
    assert not (output / "Statistics" / "Bootstrap").exists()
    assert not (output / "Statistics" / "MCMC").exists()
    assert _read_outcome(output)["status"] == "incomplete"


def test_native_mcmc_without_fitted_parameters_keeps_documented_skip_behavior(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIX = ["R1A_A", "PB", "KEX_AB"]
STATISTICS = {"MCMC" = 2}
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 0
    assert not (output / "Statistics" / "MCMC").exists()


def test_real_mc_fit_is_wholly_native_and_writes_product_statistics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 2')
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MonteCarlo"
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    _assert_truthful_resampling_summary(statistics / "summary.toml")
    samples = (statistics / "samples.tsv").read_text(encoding="utf-8")
    assert len(samples.splitlines()) == 3
    diagnostics = (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    assert "requested_samples = 2" in diagnostics
    assert "completed_samples = 2" in diagnostics
    assert 'engine = "native direct TRF"' in diagnostics
    assert "root_seed = 0" in diagnostics
    assert 'status = "complete"' in diagnostics
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "# ±" in fitted
    assert "half_percentile_68_width" not in fitted


def test_resampling_replicates_do_not_create_progress_streams(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 2')

    run(
        _fit_arguments(tmp_path / "Output", method),
        session=AnalysisSession.create(),
    )

    rendered = capsys.readouterr().out
    assert rendered.count("eval 0/") == 1


def test_native_statistics_filter_resolution_fails_closed_before_filtering(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 1')
    session = AnalysisSession.create()

    with (
        patch.object(
            session,
            "resolve_current_values",
            side_effect=RuntimeError("native values unavailable"),
        ),
        patch(
            "chemex.containers.experiments.Experiments.filter_from_values"
        ) as filter_from_values,
        pytest.raises(RuntimeError, match="native values unavailable"),
    ):
        run(_fit_arguments(output, method), session=session)

    filter_from_values.assert_not_called()


def test_real_bs_fit_uses_native_refits_and_writes_bootstrap_products(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"BS" = 2')
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "Bootstrap"
    assert (
        len((statistics / "samples.tsv").read_text(encoding="utf-8").splitlines()) == 3
    )
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "diagnostics.toml").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    _assert_truthful_resampling_summary(statistics / "summary.toml")


def test_real_bsn_fit_uses_native_nucleus_resampling_products(tmp_path: Path) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"BSN" = 2')
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "BootstrapNS"
    assert (
        len((statistics / "samples.tsv").read_text(encoding="utf-8").splitlines()) == 3
    )
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "diagnostics.toml").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    _assert_truthful_resampling_summary(statistics / "summary.toml")


def test_seeded_native_mc_products_are_ordered_across_worker_counts(
    tmp_path: Path,
) -> None:
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 3')
    serial = tmp_path / "Serial"
    parallel = tmp_path / "Parallel"

    run(_fit_arguments(serial, method, workers=1), session=AnalysisSession.create())
    run(_fit_arguments(parallel, method, workers=2), session=AnalysisSession.create())

    serial_samples = (serial / "Statistics" / "MonteCarlo" / "samples.tsv").read_text(
        encoding="utf-8"
    )
    parallel_samples = (
        parallel / "Statistics" / "MonteCarlo" / "samples.tsv"
    ).read_text(encoding="utf-8")
    assert serial_samples == parallel_samples
    parallel_diagnostics = (
        parallel / "Statistics" / "MonteCarlo" / "diagnostics.toml"
    ).read_text(encoding="utf-8")
    assert "workers = 2" in parallel_diagnostics
    assert "root_seed = 0" in parallel_diagnostics


def test_native_statistics_do_not_mutate_committed_central_values(
    tmp_path: Path,
) -> None:
    central_session = AnalysisSession.create()
    statistics_session = AnalysisSession.create()
    statistics_method = _statistics_method(
        tmp_path / "method.toml",
        '"MC" = 1, "BS" = 1, "BSN" = 1',
    )

    run(_fit_arguments(tmp_path / "Central"), session=central_session)
    run(
        _fit_arguments(tmp_path / "Statistics", statistics_method),
        session=statistics_session,
    )

    central = central_session.analysis_values.snapshot()
    after_statistics = statistics_session.analysis_values.snapshot()
    assert central.revision == after_statistics.revision == 1
    assert central.items() == after_statistics.items()


def test_failed_native_replicate_keeps_central_fit_and_suppresses_complete_products(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 2')
    session = AnalysisSession.create()
    real_least_squares = direct_trf_module.least_squares
    call_count = 0

    def fail_first_replicate(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise RuntimeError("invalid replicate backend")
        return real_least_squares(*args, **kwargs)

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=fail_first_replicate,
        ),
        pytest.raises(NativeResamplingIncompleteError, match="1 of 2"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    fitted_record = next(line for line in fitted.splitlines() if "=" in line)
    fitted_value = float(fitted_record.split("=", 1)[1].split()[0])
    assert fitted_value == pytest.approx(2.34742, rel=5.0e-6)
    statistics = output / "Statistics" / "MonteCarlo"
    assert (
        len((statistics / "samples.tsv").read_text(encoding="utf-8").splitlines()) == 2
    )
    failures = (statistics / "failures.tsv").read_text(encoding="utf-8")
    assert "direct_trf_solver_unsuccessful" in failures
    diagnostics = (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    assert 'status = "incomplete"' in diagnostics
    assert "requested_samples = 2" in diagnostics
    assert "completed_samples = 1" in diagnostics
    assert "failed_samples = 1" in diagnostics
    assert "unstarted_samples = 0" in diagnostics
    assert 'samples_file = "samples.tsv"' in diagnostics
    assert 'failures_file = "failures.tsv"' in diagnostics
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


def test_interrupted_native_statistics_report_truthful_disposition_counts(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 3')
    session = AnalysisSession.create()
    real_generate = native_resampling_module.generate_resampling_draw

    def interrupt_second_draw(dataset, request):
        if request.ordinal == 1:
            raise KeyboardInterrupt
        return real_generate(dataset, request)

    with (
        patch(
            "chemex.optimize.native_resampling.generate_resampling_draw",
            side_effect=interrupt_second_draw,
        ),
        pytest.raises(NativeResamplingIncompleteError, match="1 of 3"),
    ):
        run(_fit_arguments(output, method), session=session)

    statistics = output / "Statistics" / "MonteCarlo"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    disposition_counts = {
        "completed": diagnostics["completed_samples"],
        "failed": diagnostics["failed_samples"],
        "cancelled": diagnostics["cancelled_samples"],
        "interrupted": diagnostics["interrupted_samples"],
        "unstarted": diagnostics["unstarted_samples"],
    }
    assert disposition_counts == {
        "completed": 1,
        "failed": 0,
        "cancelled": 0,
        "interrupted": 1,
        "unstarted": 1,
    }
    assert sum(disposition_counts.values()) == diagnostics["requested_samples"]
    failures = (statistics / "failures.tsv").read_text(encoding="utf-8")
    assert "\tinterrupted\t" in failures
    assert "\tnot_started\t" in failures
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


def test_native_statistics_setup_failure_publishes_incomplete_diagnostics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 2')
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.resampling._native_dataset",
            side_effect=RuntimeError("manifest unavailable"),
        ),
        pytest.raises(RuntimeError, match="manifest unavailable"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MonteCarlo"
    diagnostics = (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    assert 'status = "incomplete"' in diagnostics
    assert 'terminal = "failed"' in diagnostics
    assert "completed_samples = 0" in diagnostics
    assert "failed_samples = 0" in diagnostics
    assert "unstarted_samples = 2" in diagnostics
    assert 'failure_type = "RuntimeError"' in diagnostics
    assert 'failure_message = "manifest unavailable"' in diagnostics
    assert "samples_file" not in diagnostics
    assert "failures_file" not in diagnostics
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


@pytest.mark.parametrize(
    ("family", "directory"),
    (
        ("MC", "MonteCarlo"),
        ("BS", "Bootstrap"),
        ("BSN", "BootstrapNS"),
    ),
)
def test_native_resampling_publication_failure_is_terminal_and_fail_closed(
    tmp_path: Path,
    family: str,
    directory: str,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(
        tmp_path / "method.toml",
        f'"{family}" = 1',
    )

    with (
        patch(
            "chemex.optimize.resampling._write_resampling_correlations",
            side_effect=RuntimeError("correlation publication failed"),
        ),
        pytest.raises(RuntimeError, match="correlation publication failed"),
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    statistics = output / "Statistics" / directory
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["terminal"] == "failed"
    assert diagnostics["completed_samples"] == 1
    assert diagnostics["failed_samples"] == 0
    assert diagnostics["cancelled_samples"] == 0
    assert diagnostics["interrupted_samples"] == 0
    assert diagnostics["unstarted_samples"] == 0
    assert diagnostics["failure_type"] == "RuntimeError"
    assert diagnostics["failure_message"] == "correlation publication failed"
    assert (
        len((statistics / "samples.tsv").read_text(encoding="utf-8").splitlines()) == 2
    )
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()
    assert _read_outcome(output)["status"] == "incomplete"


def test_native_resampling_materialization_failure_replaces_running_diagnostics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 1')

    with (
        patch.object(
            resampling_module,
            "_as_sample_array",
            side_effect=RuntimeError("sample materialization failed"),
        ),
        pytest.raises(RuntimeError, match="sample materialization failed"),
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    statistics = output / "Statistics" / "MonteCarlo"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["completed_samples"] == 1
    assert diagnostics["failure_message"] == "sample materialization failed"
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()
    assert _read_outcome(output)["status"] == "incomplete"


def test_resampling_diagnostics_failure_preserves_original_publication_error(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 1')
    real_atomic_write = resampling_module.write_text_atomic

    def fail_incomplete_diagnostics(destination: Path, content: str) -> None:
        if (
            destination.name == "diagnostics.toml"
            and 'status = "incomplete"' in content
        ):
            raise OSError("incomplete diagnostics publication failed")
        real_atomic_write(destination, content)

    with (
        patch.object(
            resampling_module,
            "_write_resampling_correlations",
            side_effect=RuntimeError("correlation publication failed"),
        ),
        patch.object(
            resampling_module,
            "write_text_atomic",
            side_effect=fail_incomplete_diagnostics,
        ),
        pytest.raises(RuntimeError, match="correlation publication failed"),
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    statistics = output / "Statistics" / "MonteCarlo"
    assert not (statistics / "diagnostics.toml").exists()
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()
    outcome = _read_outcome(output)
    assert outcome["status"] == "incomplete"
    assert outcome["failure_message"] == "correlation publication failed"


def test_real_grid_fit_uses_native_cartesian_trf_and_writes_grid_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Grid" / "grid.out").is_file()
    assert (output / "Parameters" / "fixed.toml").is_file()


def test_real_grid_progress_aggregates_local_seed_polishes(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    method = _grid_method(tmp_path / "method.toml")

    run(
        _fit_arguments(tmp_path / "Output", method),
        session=AnalysisSession.create(),
    )

    rendered = capsys.readouterr().out
    assert rendered.count("GRID seed") == 1
    assert "final χ²" in rendered


def test_real_grouped_grid_fit_uses_one_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")

    run(
        _fit_arguments(
            output,
            method,
            include=("G2N-HN", "H3N-HN"),
        ),
        session=session,
    )

    assert session.analysis_values.snapshot().revision == 1
    group_grid_outputs = tuple((output / "Grid" / "Groups").glob("*.out"))
    assert len(group_grid_outputs) == 2
    assert (output / "All" / "Parameters" / "fixed.toml").is_file()
    assert not (output / "Groups").exists()


def test_ordered_native_steps_carry_forward_committed_state_and_v1_roles(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[STEP1]
FIX = ["PB", "KEX_AB"]

[STEP2]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 2
    for step in ("STEP1", "STEP2"):
        fitted = (output / step / "Parameters" / "fitted.toml").read_text(
            encoding="utf-8"
        )
        assert "# ±" in fitted
        fixed = (output / step / "Parameters" / "fixed.toml").read_text(
            encoding="utf-8"
        )
        assert "PB" in fixed
        assert "KEX_AB" in fixed


def test_failed_second_step_rerun_eagerly_invalidates_every_planned_step(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[STEP1]
FIX = ["PB", "KEX_AB"]

[STEP2]
""",
        encoding="utf-8",
    )
    run(_fit_arguments(output, method), session=AnalysisSession.create())
    stale_step1 = output / "STEP1" / "Parameters" / "stale.toml"
    stale_step2 = output / "STEP2" / "Parameters" / "stale.toml"
    user_file = output / "STEP2" / "notes.txt"
    stale_step1.write_text("old step 1\n", encoding="utf-8")
    stale_step2.write_text("old step 2\n", encoding="utf-8")
    user_file.write_text("keep me\n", encoding="utf-8")
    real_least_squares = direct_trf_module.least_squares
    call_count = 0

    def fail_second_step(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise RuntimeError("step 2 backend failed")
        return real_least_squares(*args, **kwargs)

    session = AnalysisSession.create()
    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=fail_second_step,
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "STEP1" / "Parameters" / "fitted.toml").is_file()
    assert not stale_step1.exists()
    assert not (output / "STEP2" / "Parameters").exists()
    assert not stale_step2.exists()
    assert user_file.read_text(encoding="utf-8") == "keep me\n"
    assert _read_outcome(output)["status"] == "incomplete"


def test_planned_invalidation_rejects_a_step_root_outside_the_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    external_result = tmp_path / "external" / "Parameters" / "keep.toml"
    external_result.parent.mkdir(parents=True)
    external_result.write_text("keep me\n", encoding="utf-8")
    method = tmp_path / "method.toml"
    method.write_text(
        """["../external"]
FIX = ["PB", "KEX_AB"]

[SAFE]
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="outside the output directory"):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    assert external_result.read_text(encoding="utf-8") == "keep me\n"
    assert _read_outcome(output)["status"] == "incomplete"


def test_native_backend_failure_cannot_commit_or_publish_fitted_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    run(_fit_arguments(output), session=AnalysisSession.create())
    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=RuntimeError("backend failed"),
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 0
    assert not (output / "Parameters").exists()
    assert not (output / "Data").exists()
    assert not (output / "Statistics" / "Covariance").exists()


def test_v1_constraint_and_profile_selection_reach_native_product_fit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIX = ["PB", "KEX_AB"]
CONSTRAINTS = ["[R1A_A, NUC->G2N-H] = 2.0"]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(
        _fit_arguments(
            output,
            method,
            include=("G2N-HN", "H3N-HN"),
        ),
        session=session,
    )

    assert session.analysis_values.snapshot().revision == 1
    constrained = tuple((output / "All" / "Parameters").glob("constrained.toml"))
    assert constrained
    constrained_text = constrained[0].read_text(encoding="utf-8")
    assert "G2N-H" in constrained_text
    assert "2.00000e+00" in constrained_text


def test_all_fixed_v1_method_evaluates_without_optimizer_or_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIX = ["R1A_A", "PB", "KEX_AB"]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 0
    assert (output / "Parameters" / "fixed.toml").is_file()
    assert (output / "statistics.toml").is_file()
    data = next((output / "Data").glob("*.dat")).read_text(encoding="utf-8")
    assert "1.00000000e+32" not in data


def test_explicit_trf_with_statistics_uses_native_production_dispatch(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FITMETHOD = "trf"
STATISTICS = { "MC" = 1 }
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Statistics" / "MonteCarlo" / "samples.tsv").is_file()


def test_grid_with_statistics_warns_and_runs_native_grid_only(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
GRID = ["[R1A_A] = (1.0, 3.0)"]
STATISTICS = { "MC" = 1 }
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    with patch(
        "chemex.optimize.fitting.run_native_resampling_statistics",
        side_effect=AssertionError("statistics after GRID were entered"),
    ):
        run(_fit_arguments(output, method), session=session)

    captured = capsys.readouterr()
    assert "GRID" in captured.out
    assert "STATISTICS" in captured.out
    assert (output / "Grid").is_dir()
    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()
    assert not (output / "Statistics" / "MonteCarlo").exists()
