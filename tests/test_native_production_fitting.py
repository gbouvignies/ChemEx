from __future__ import annotations

import json
import math
import subprocess
import sys
import tomllib
from pathlib import Path
from types import SimpleNamespace
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
import chemex.printers.grid as grid_printer_module
import chemex.run_info as run_info_module
from chemex.chemex import run
from chemex.cli import build_parser
from chemex.evaluation.native import BoundEvaluator, EvaluationFailure
from chemex.optimize.mcmc import NativeMcmcIncompleteError
from chemex.optimize.progress import ProgressPhase
from chemex.optimize.resampling import NativeResamplingIncompleteError
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXAMPLE = ROOT / "examples/Experiments/RELAXATION_HZNZ"
EXPERIMENT = EXAMPLE / "Experiments/800mhz.toml"
PARAMETERS = EXAMPLE / "Parameters/parameters.toml"
METHOD = EXAMPLE / "Methods/method.toml"
DCEST_EXAMPLE = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH"
DCEST_EXPERIMENT = DCEST_EXAMPLE / "Experiments/3hz.toml"
DCEST_PARAMETERS = DCEST_EXAMPLE / "Parameters/parameters.toml"
THREE_STATE_DCEST_EXAMPLE = ROOT / "examples/Experiments/DCEST_15N_3States"
THREE_STATE_DCEST_EXPERIMENT = THREE_STATE_DCEST_EXAMPLE / "Experiments/1.25hz.toml"
THREE_STATE_DCEST_PARAMETERS = THREE_STATE_DCEST_EXAMPLE / "Parameters/parameters.toml"
CEST_1HN_IP_AP_EXAMPLE = ROOT / "examples/Experiments/CEST_1HN_IP_AP"
CEST_1HN_IP_AP_EXPERIMENT = CEST_1HN_IP_AP_EXAMPLE / "Experiments/30hz.toml"
CEST_1HN_IP_AP_PARAMETERS = CEST_1HN_IP_AP_EXAMPLE / "Parameters/parameters.toml"
CEST_1HN_IP_AP_METHOD = CEST_1HN_IP_AP_EXAMPLE / "Methods/method.toml"


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


def _three_state_dcest_arguments(output: Path, method: Path):
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            str(THREE_STATE_DCEST_EXPERIMENT),
            "-p",
            str(THREE_STATE_DCEST_PARAMETERS),
            "-m",
            str(method),
            "-o",
            str(output),
            "--model",
            "3st",
            "--include",
            "K19N-HN",
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )


def _cest_1hn_ip_ap_arguments(output: Path):
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            str(CEST_1HN_IP_AP_EXPERIMENT),
            "-p",
            str(CEST_1HN_IP_AP_PARAMETERS),
            "-m",
            str(CEST_1HN_IP_AP_METHOD),
            "-o",
            str(output),
            "--model",
            "2st.rs",
            "--include",
            "4",
            "--plot",
            "nothing",
            "--workers",
            "1",
            "--native-threads",
            "1",
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

    committed = session.analysis_values.snapshot()
    assert committed.revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    fitted_record = next(line for line in fitted.splitlines() if "=" in line)
    fitted_value = float(fitted_record.split("=", 1)[1].split()[0])
    assert fitted_value == pytest.approx(2.34742, rel=5.0e-6)
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    fitted_id = next(
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name == "R1A_A" and definition.spin_system_name == "G2N-H"
    )
    assert fitted_value == pytest.approx(committed[fitted_id], rel=5.0e-6)
    assert "# ±" in fitted_record
    fitted_error = float(fitted_record.split("±", 1)[1])
    assert fitted_error == pytest.approx(0.083366086, rel=3.0e-6)
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
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()
    assert not (output / "Components").exists()
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


def test_cest_1hn_ip_ap_commits_finite_negative_derived_r2a_b(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    committed_snapshots: list[AnalysisValuesSnapshot] = []
    real_publish = run_info_module.RunInfo.publish_restart

    def record_restart(run_info, snapshot):
        real_publish(run_info, snapshot)
        committed_snapshots.append(snapshot)

    with patch.object(
        run_info_module.RunInfo,
        "publish_restart",
        autospec=True,
        side_effect=record_restart,
    ):
        run(_cest_1hn_ip_ap_arguments(output), session=session)

    snapshot = session.analysis_values.snapshot()
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    r2a_b = next(
        definition
        for definition in parameter_model.definitions
        if definition.name == "R2A_B" and definition.spin_system_name == "4H"
    )
    configuration = parameter_model.configuration[r2a_b.param_id]
    assert configuration.lower_bound == 0.0
    assert [item.revision for item in committed_snapshots] == [1, 2]
    # This is the one-profile form of the shipped regression. A 2e-7 relative
    # tolerance covers backend/finite-difference rounding while remaining seven
    # orders of magnitude from the configured zero bound and lmfit's old clipping.
    assert committed_snapshots[0][r2a_b.param_id] == pytest.approx(
        -0.42454801400787906,
        rel=2.0e-7,
        abs=1.0e-9,
    )
    assert committed_snapshots[0][r2a_b.param_id] < configuration.lower_bound
    assert snapshot.revision == 2


def test_product_fit_rejects_a_stale_aggregate_commit_atomically(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    real_commit = direct_trf_module.execute_fit_commit

    def advance_revision_then_commit(*args, **kwargs):
        current = session.analysis_values.snapshot()
        session.analysis_values.commit(
            dict(current),
            expected=current,
            scope=tuple(current),
        )
        return real_commit(*args, **kwargs)

    with (
        patch.object(
            native_deterministic_module,
            "execute_fit_commit",
            side_effect=advance_revision_then_commit,
        ),
        pytest.raises(RuntimeError, match="stale_revision"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert not (output / "run_info" / "restart.toml").exists()
    outcome = _read_outcome(output)
    assert outcome["latest_committed_revision"] == 1
    assert outcome["restart_revision"] == 0
    assert not (output / "Parameters").exists()


def test_ordinary_fit_requires_no_post_seal_parameter_store_value_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    session = AnalysisSession.create()
    build_analysis_values = session.try_build_analysis_values

    def seal_then_reject_store_value_mutation() -> bool:
        initialized = build_analysis_values()

        def reject_store_value_mutation(_values: object) -> None:
            raise AssertionError("post-seal ParameterStore value mutation")

        monkeypatch.setattr(
            session.parameters.database,
            "_set_values",
            reject_store_value_mutation,
        )

        def reject_store_value_read(_param_ids: object) -> object:
            raise AssertionError("post-seal ParameterStore value read")

        monkeypatch.setattr(
            session.parameters,
            "get_parameters",
            reject_store_value_read,
        )
        return initialized

    monkeypatch.setattr(
        session,
        "try_build_analysis_values",
        seal_then_reject_store_value_mutation,
    )

    run(_fit_arguments(tmp_path / "Output"), session=session)

    assert session.analysis_values.snapshot().revision == 1


def test_v2_role_change_keeps_the_prior_committed_value_without_store_roles(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """
FORMAT_VERSION = 2
[FIRST]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]

[SECOND]
ROLES_FROM = "FIRST"
ROLES = [{ FIX = ["R1A_A"] }]
""",
        encoding="utf-8",
    )
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    first_fitted = (output / "FIRST/Parameters/fitted.toml").read_text(encoding="utf-8")
    second_fixed = (output / "SECOND/Parameters/fixed.toml").read_text(encoding="utf-8")
    first_value = float(
        next(line for line in first_fitted.splitlines() if "G2N-H" in line)
        .split("=", 1)[1]
        .split()[0]
    )
    second_value = float(
        next(line for line in second_fixed.splitlines() if "G2N-H" in line)
        .split("=", 1)[1]
        .split()[0]
    )

    assert session.analysis_values.snapshot().revision == 1
    assert second_value == pytest.approx(first_value, rel=1.0e-12)


def test_constrained_value_becomes_the_next_independent_start(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """
FORMAT_VERSION = 2
[CONSTRAINED]
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
  { CONSTRAIN = ["[R1A_A, NUC->G2N-H] = 2.0"] },
]

[INDEPENDENT]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]
""",
        encoding="utf-8",
    )
    starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_start(*args, **kwargs):
        starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    output = tmp_path / "Output"
    session = AnalysisSession.create()
    published: list[tuple[int, str]] = []
    real_publish = run_info_module.RunInfo.publish_restart

    def record_restart(run_info, snapshot):
        real_publish(run_info, snapshot)
        published.append(
            (snapshot.revision, (run_info.path / "restart.toml").read_text())
        )

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=record_start,
        ),
        patch.object(
            run_info_module.RunInfo,
            "publish_restart",
            autospec=True,
            side_effect=record_restart,
        ),
    ):
        run(_fit_arguments(output, method), session=session)

    assert starts == [(2.0,)]
    assert session.analysis_values.snapshot().revision == 2
    assert [revision for revision, _text in published] == [1, 2]
    assert '"G2N-H" = [2.0,' in published[0][1]
    assert _read_outcome(output)["restart_revision"] == 2


def test_independent_constrained_independent_transition_uses_latest_resolved_value(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """
FORMAT_VERSION = 2
[FIRST_INDEPENDENT]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]

[CONSTRAINED]
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
  { CONSTRAIN = ["[R1A_A, NUC->G2N-H] = 3.0"] },
]

[SECOND_INDEPENDENT]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]
""",
        encoding="utf-8",
    )
    starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_start(*args, **kwargs):
        starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    output = tmp_path / "Output"
    session = AnalysisSession.create()
    published: list[tuple[int, str]] = []
    real_publish = run_info_module.RunInfo.publish_restart

    def record_restart(run_info, snapshot):
        real_publish(run_info, snapshot)
        published.append(
            (snapshot.revision, (run_info.path / "restart.toml").read_text())
        )

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=record_start,
        ),
        patch.object(
            run_info_module.RunInfo,
            "publish_restart",
            autospec=True,
            side_effect=record_restart,
        ),
    ):
        run(_fit_arguments(output, method), session=session)

    assert starts[1] == (3.0,)
    assert session.analysis_values.snapshot().revision == 3
    assert [revision for revision, _text in published] == [1, 2, 3]
    assert '"G2N-H" = [3.0,' in published[1][1]
    assert _read_outcome(output)["restart_revision"] == 3


def test_out_of_bounds_derived_continuity_cannot_become_held_independent(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """
FORMAT_VERSION = 2
[DERIVED]
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
  { CONSTRAIN = ["[R1A_A, NUC->G2N-H] = -1.0"] },
]

[HELD]
ROLES = [{ FIX = ["PB", "KEX_AB", "R1A_A"] }]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    with pytest.raises(
        direct_trf_module.DirectTrfConstructionError,
        match="Independent parameter .* outside its effective bounds",
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)

    committed = session.analysis_values.snapshot()
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    r1a_a = next(
        item.param_id
        for item in parameter_model.definitions
        if item.name == "R1A_A" and item.spin_system_name == "G2N-H"
    )
    assert committed.revision == 1
    assert committed[r1a_a] == -1.0


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
        (output / "Statistics" / "Covariance" / "evidence.json").read_text(
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
    assert "Evaluations" in rendered
    assert "Reduced χ²" in rendered
    minimizing_terminal = rendered.index("Evaluations")
    uncertainty_start = rendered.index("Estimating parameter uncertainties")
    uncertainty_terminal = rendered.index("covariance available")
    assert minimizing_terminal < uncertainty_start < uncertainty_terminal
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
    assert "Evaluations" in rendered


def test_cancellation_after_aggregate_acceptance_blocks_commit_and_later_steps(
    tmp_path: Path,
) -> None:
    method = tmp_path / "two-steps.toml"
    method.write_text(
        """[STEP1]
FIX = ["PB", "KEX_AB"]

[STEP2]
FIX = ["PB", "KEX_AB"]
""",
        encoding="utf-8",
    )
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    real_execute = native_deterministic_module.execute_grouped_direct_trf
    executions = 0

    def cancel_after_acceptance(*args, **kwargs):
        nonlocal executions
        executions += 1
        outcome = real_execute(*args, **kwargs)
        assert outcome.accepted_result is not None
        kwargs["cancellation"].cancel()
        return outcome

    with (
        patch.object(
            native_deterministic_module,
            "execute_grouped_direct_trf",
            side_effect=cancel_after_acceptance,
        ),
        patch.object(
            native_deterministic_module,
            "derive_uncertainty_evidence",
            side_effect=AssertionError("uncertainty ran before a committed fit"),
        ) as uncertainty,
        pytest.raises(KeyboardInterrupt, match="cancelled before commit"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert executions == 1
    assert uncertainty.call_count == 0
    assert session.analysis_values.snapshot().revision == 0
    assert not (output / "STEP1" / "Parameters").exists()
    assert not (output / "STEP2").exists()


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


def test_committed_restart_supports_real_cli_round_trip(tmp_path: Path) -> None:
    first_output = tmp_path / "First"
    restarted_output = tmp_path / "Restarted"

    _run_real_fit_cli(first_output, METHOD, PARAMETERS)
    restart_parameters = first_output / "run_info" / "restart.toml"
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


def test_original_start_supports_real_cli_round_trip(tmp_path: Path) -> None:
    first_output = tmp_path / "First"
    reproduced_output = tmp_path / "Reproduced"

    _run_real_fit_cli(first_output, METHOD, PARAMETERS)
    original_start = first_output / "run_info" / "parameters_used.toml"
    _run_real_fit_cli(reproduced_output, METHOD, original_start)

    first = tomllib.loads(
        (first_output / "run_info" / "parameters_used.toml").read_text(encoding="utf-8")
    )
    reproduced = tomllib.loads(
        (reproduced_output / "run_info" / "parameters_used.toml").read_text(
            encoding="utf-8"
        )
    )
    assert reproduced == first


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
    assert outcome["schema_version"] == 2
    assert outcome["status"] == "incomplete"
    assert outcome["latest_committed_revision"] == 0
    assert outcome["restart_revision"] == 0
    assert outcome["terminal"] == "failed"
    assert outcome["failure_stage"] == "deterministic_fit"
    assert outcome["failure_type"] == "RuntimeError"
    assert not (output / "Parameters").exists()
    assert not (output / "Data").exists()
    assert not (output / "Plots").exists()
    assert not (output / "statistics.toml").exists()
    assert user_file.read_text(encoding="utf-8") == "keep me\n"


def test_planned_output_cleanup_failure_is_classified_as_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"

    with (
        patch(
            "chemex.chemex.invalidate_planned_outputs",
            side_effect=OSError("planned output cleanup failed"),
        ),
        pytest.raises(OSError, match="planned output cleanup failed"),
    ):
        run(_fit_arguments(output), session=AnalysisSession.create())

    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 0,
        "restart_revision": 0,
        "terminal": "failed",
        "failure_stage": "output",
        "failure_type": "OSError",
        "failure_message": "planned output cleanup failed",
    }


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
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 1,
        "restart_revision": 1,
        "terminal": "failed",
        "failure_stage": "output",
        "failure_type": "RuntimeError",
        "failure_message": "data writer failed",
    }
    assert (output / "run_info" / "restart.toml").is_file()
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert not stale_parameter.exists()
    assert not stale_data.exists()
    assert not (output / "statistics.toml").exists()


def test_first_restart_publication_failure_preserves_committed_authority(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    real_atomic_write = run_info_module.write_text_atomic

    def fail_restart(destination: Path, content: str) -> None:
        if destination.name == "restart.toml":
            raise OSError("restart publication failed")
        real_atomic_write(destination, content)

    with (
        patch.object(run_info_module, "write_text_atomic", side_effect=fail_restart),
        pytest.raises(OSError, match="restart publication failed"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert not (output / "run_info" / "restart.toml").exists()
    assert not (output / "Parameters").exists()
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 1,
        "restart_revision": 0,
        "terminal": "failed",
        "failure_stage": "restart_publication",
        "failure_type": "OSError",
        "failure_message": "restart publication failed",
    }


def test_first_restart_serialization_failure_preserves_committed_authority(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    real_serialize = run_info_module.serialize_parameter_file

    def fail_restart_serialization(*args, **kwargs) -> str:
        if kwargs.get("state_kind") == "restart":
            raise RuntimeError("restart serialization failed")
        return real_serialize(*args, **kwargs)

    with (
        patch.object(
            run_info_module,
            "serialize_parameter_file",
            side_effect=fail_restart_serialization,
        ),
        pytest.raises(RuntimeError, match="restart serialization failed"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert not (output / "run_info" / "restart.toml").exists()
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 1,
        "restart_revision": 0,
        "terminal": "failed",
        "failure_stage": "restart_publication",
        "failure_type": "RuntimeError",
        "failure_message": "restart serialization failed",
    }


def test_later_restart_publication_failure_preserves_previous_checkpoint(
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
    real_atomic_write = run_info_module.write_text_atomic
    restart_writes = 0
    first_restart = b""

    def fail_second_restart(destination: Path, content: str) -> None:
        nonlocal restart_writes, first_restart
        if destination.name == "restart.toml":
            restart_writes += 1
            if restart_writes == 2:
                raise OSError("second restart publication failed")
        real_atomic_write(destination, content)
        if destination.name == "restart.toml":
            first_restart = destination.read_bytes()

    with (
        patch.object(
            run_info_module,
            "write_text_atomic",
            side_effect=fail_second_restart,
        ),
        pytest.raises(OSError, match="second restart publication failed"),
    ):
        run(_fit_arguments(output, method), session=session)

    restart = output / "run_info" / "restart.toml"
    assert session.analysis_values.snapshot().revision == 2
    assert restart.read_bytes() == first_restart
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 2,
        "restart_revision": 1,
        "terminal": "failed",
        "failure_stage": "restart_publication",
        "failure_type": "OSError",
        "failure_message": "second restart publication failed",
    }


def test_later_restart_serialization_failure_preserves_previous_checkpoint(
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
    real_serialize = run_info_module.serialize_parameter_file
    real_atomic_write = run_info_module.write_text_atomic
    restart_serializations = 0
    first_restart = b""

    def fail_second_restart_serialization(*args, **kwargs) -> str:
        nonlocal restart_serializations
        if kwargs.get("state_kind") == "restart":
            restart_serializations += 1
            if restart_serializations == 2:
                raise RuntimeError("second restart serialization failed")
        return real_serialize(*args, **kwargs)

    def capture_first_restart(destination: Path, content: str) -> None:
        nonlocal first_restart
        real_atomic_write(destination, content)
        if destination.name == "restart.toml":
            first_restart = destination.read_bytes()

    with (
        patch.object(
            run_info_module,
            "serialize_parameter_file",
            side_effect=fail_second_restart_serialization,
        ),
        patch.object(
            run_info_module,
            "write_text_atomic",
            side_effect=capture_first_restart,
        ),
        pytest.raises(RuntimeError, match="second restart serialization failed"),
    ):
        run(_fit_arguments(output, method), session=session)

    restart = output / "run_info" / "restart.toml"
    assert session.analysis_values.snapshot().revision == 2
    assert restart.read_bytes() == first_restart
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 2,
        "restart_revision": 1,
        "terminal": "failed",
        "failure_stage": "restart_publication",
        "failure_type": "RuntimeError",
        "failure_message": "second restart serialization failed",
    }


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
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 1,
        "restart_revision": 1,
        "terminal": "failed",
        "failure_stage": "run_info",
        "failure_type": "OSError",
        "failure_message": "complete outcome publication failed",
    }


def test_real_grouped_direct_fit_uses_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    run(
        _fit_arguments(
            output,
            include=("G2N-HN", "H3N-HN"),
            plot_level="normal",
        ),
        session=session,
    )

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert (output / "Data").is_dir()
    assert (output / "Plots").is_dir()
    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()
    assert (output / "statistics.toml").is_file()
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()
    assert not (output / "Components").exists()


def test_reused_output_removes_legacy_and_development_result_trees(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    for name in ("All", "Groups", "Components"):
        stale = output / name / "stale.txt"
        stale.parent.mkdir(parents=True)
        stale.write_text("obsolete result\n", encoding="utf-8")

    run(
        _fit_arguments(output, include=("G2N-HN", "H3N-HN")),
        session=AnalysisSession.create(),
    )

    assert (output / "Parameters" / "fitted.toml").is_file()
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()
    assert not (output / "Components").exists()
    assert (output / "run_info" / "outcome.toml").is_file()


def test_grouped_covariance_warns_only_near_bound_independent_coordinate(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "Output"
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        PARAMETERS.read_text(encoding="utf-8")
        + """

[R1A_A]
G2N-H = [2.3, 2.18, 5.0]
""",
        encoding="utf-8",
    )
    original_svd = uncertainty_module.svd
    covariance_svd_call = 0

    def one_bad_block_svd(*args, **kwargs):
        nonlocal covariance_svd_call
        covariance_svd_call += 1
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        if covariance_svd_call == 1:
            singular[-1] = 0.0
        return left, singular, right

    with patch("chemex.optimize.uncertainty.svd", side_effect=one_bad_block_svd):
        run(
            _fit_arguments(
                output,
                parameters=parameters,
                include=("G2N-HN", "H3N-HN"),
            ),
            session=AnalysisSession.create(),
        )

    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert fitted.count("# ±") == 2
    assert fitted.count("boundary may make uncertainty asymmetric") == 1
    g2_record = next(
        line for line in fitted.splitlines() if line.lstrip().startswith("G2N-H")
    )
    h3_record = next(
        line for line in fitted.splitlines() if line.lstrip().startswith("H3N-H")
    )
    assert "boundary may make uncertainty asymmetric" in g2_record
    assert "boundary may make uncertainty asymmetric" not in h3_record
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert constrained.count("# ±") == 2
    assert constrained.count("boundary may make uncertainty asymmetric") == 1
    assert "error unavailable" not in constrained
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()
    blocks = json.loads(
        (output / "Statistics" / "Covariance" / "blocks.json").read_text(
            encoding="utf-8"
        )
    )
    assert len(blocks["partition_proof_identity"]) == 64
    assert all(block["scope_reportable"] for block in blocks["blocks"])
    claim_sets = [
        {claim["name"]: claim["state"] for claim in block["claims"]}
        for block in blocks["blocks"]
    ]
    assert sorted(claims["BOUNDARY_SEPARATION"] for claims in claim_sets) == [
        "satisfied",
        "violated",
    ]
    assert all(
        claims["USABLE_LOCAL_COVARIANCE"] == "satisfied" for claims in claim_sets
    )
    rendered = " ".join(capsys.readouterr().out.split())
    assert rendered.count("covariance available with boundary warnings") == 1


def test_multivariate_boundary_warning_product_is_coordinate_specific(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
FIT = ["D2O", "KDH"]
FIX = ["CS_A", "DW_AB", "PHI", "R1_A", "R2_A", "R2_B"]
""",
        encoding="utf-8",
    )
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        DCEST_PARAMETERS.read_text(encoding="utf-8").replace(
            "KDH = 10.0",
            "KDH = [10.0, 9.99, 10.0]",
            1,
        ),
        encoding="utf-8",
    )
    arguments = build_parser().parse_args(
        [
            "fit",
            "-e",
            str(DCEST_EXPERIMENT),
            "-p",
            str(parameters),
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
    assert fitted.count("# ±") == 2
    assert fitted.count("boundary may make uncertainty asymmetric") == 1
    d2o_record = next(line for line in fitted.splitlines() if "D2O" in line)
    kdh_record = next(line for line in fitted.splitlines() if line.startswith("1 ="))
    assert "boundary may make uncertainty asymmetric" not in d2o_record
    assert "boundary may make uncertainty asymmetric" in kdh_record
    rendered = " ".join(capsys.readouterr().out.split())
    assert rendered.count("covariance available with boundary warnings") == 1


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
    covariance_path = output / "Statistics" / "Covariance"
    assert (covariance_path / "evidence.json").is_file()
    assert not (covariance_path / "blocks.json").exists()
    status = json.loads((covariance_path / "status.json").read_text(encoding="utf-8"))
    assert status["status"] == "incomplete"
    assert status["terminal"] == "interrupted"
    assert status["reason"] == "derivation interrupted/cancelled"
    assert (output / "Parameters" / "fitted.toml").is_file()


def test_shared_parameter_coupling_is_not_split_into_covariance_blocks(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
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
    rendered = " ".join(capsys.readouterr().out.split())
    assert "uncertainty unavailable: rank deficient" in rendered


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
    assert "component 1/2" in rendered
    assert rendered.count("eval 0/") == 1
    assert "eval " in rendered
    assert "Evaluations" in rendered
    assert "group " not in rendered.lower()
    assert rendered.count("Writing results in") == 1
    assert "All groups" not in rendered
    assert "Group " not in rendered


def test_boundary_warning_fit_reports_uncertainty_without_changing_central_values(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
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
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "# ±" in fitted
    assert "boundary may make uncertainty asymmetric" in fitted
    assert "error unavailable: boundary limited" not in fitted
    assert "error unavailable: boundary limited" not in constrained
    assert "error unavailable: constrained propagation unavailable" in constrained
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
    rendered = " ".join(capsys.readouterr().out.split())
    assert rendered.count("covariance available with boundary warnings") == 1


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
        "chemex.optimize.native_deterministic.derive_uncertainty_evidence",
        side_effect=KeyboardInterrupt,
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "error unavailable: derivation interrupted/cancelled" in fitted
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "error unavailable: derivation interrupted/cancelled" in constrained
    assert "# ±" not in constrained
    covariance_path = output / "Statistics" / "Covariance"
    assert not (covariance_path / "evidence.json").exists()
    status = json.loads((covariance_path / "status.json").read_text(encoding="utf-8"))
    assert status == {
        "artifact_type": "native_covariance_derivation_status",
        "reason": "derivation interrupted/cancelled",
        "schema_version": 1,
        "status": "incomplete",
        "terminal": "interrupted",
    }


def test_uncertainty_exception_does_not_roll_back_the_committed_fit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    with (
        patch.object(
            native_deterministic_module,
            "derive_uncertainty_evidence",
            side_effect=RuntimeError("uncertainty failed"),
        ),
        pytest.raises(RuntimeError, match="uncertainty failed"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "run_info" / "restart.toml").is_file()
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 1,
        "restart_revision": 1,
        "terminal": "failed",
        "failure_stage": "uncertainty",
        "failure_type": "RuntimeError",
        "failure_message": "uncertainty failed",
    }
    assert not (output / "Parameters").exists()


def test_product_uncertainty_starts_after_commit_and_cannot_change_central_values(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    real_derive = native_deterministic_module.derive_uncertainty_evidence

    def derive_after_commit(*args, **kwargs):
        committed = session.analysis_values.snapshot()
        assert committed.revision == 1
        evidence = real_derive(*args, **kwargs)
        assert session.analysis_values.snapshot() == committed
        return evidence

    with patch.object(
        native_deterministic_module,
        "derive_uncertainty_evidence",
        side_effect=derive_after_commit,
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 1


def _grid_method(path: Path) -> Path:
    path.write_text(
        """[GRID]
FIX = ["PB", "KEX_AB"]
GRID = ["[R1A_A] = (1.0, 3.0)"]
""",
        encoding="utf-8",
    )
    return path


def _de_method(path: Path, *, seed: int = 597) -> Path:
    path.write_text(
        f"""FORMAT_VERSION = 2
[DE]
ROLES = [{{ FIX = ["PB", "KEX_AB"] }}]

[DE.SEARCH.DE]
SEED = {seed}
COORDINATES = ["[R1A_A] = log(1.0, 4.0)"]
""",
        encoding="utf-8",
    )
    return path


def _successful_de_backend_at_x0(objective, _bounds, **kwargs):
    vector = np.asarray(kwargs["x0"], dtype=np.float64)
    value = objective(vector)
    population_size = max(5, int(kwargs["popsize"]) * vector.size)
    return SimpleNamespace(
        success=True,
        message="Optimization terminated successfully.",
        nit=1,
        nfev=1,
        x=vector,
        fun=value,
        population=np.tile(vector, (population_size, 1)),
        population_energies=np.full(population_size, value),
    )


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


def _automatic_mcmc_method(path: Path, *, steps: int = 5) -> Path:
    path.write_text(
        f"""FORMAT_VERSION = 2

[DEFAULT]
ROLES = [{{ FIX = ["PB", "KEX_AB"] }}]

[DEFAULT.STATISTICS.MCMC]
STEPS = {steps}
SEED = 698
""",
        encoding="utf-8",
    )
    return path


def test_valid_automatic_mcmc_burn_publishes_normal_posterior_outputs(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _automatic_mcmc_method(tmp_path / "method.toml")
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()
    committed_snapshots = []

    def valid_autocorrelation_time(_chain: Array, **_kwargs) -> Array:
        committed_snapshots.append(session.analysis_values.snapshot())
        return np.array([1.0])

    with patch.object(
        mcmc_module.emcee_autocorr,
        "integrated_time",
        side_effect=valid_autocorrelation_time,
    ):
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=session,
        )

    assert len(committed_snapshots) == 1
    assert session.analysis_values.snapshot() == committed_snapshots[0]
    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["discarded_steps"] == 2
    assert diagnostics["retained_steps"] == 3
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "samples.tsv").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    assert not (statistics / "raw_chain.tsv").exists()


@pytest.mark.parametrize(
    ("autocorrelation_outcome", "failure_reason", "autocorrelation_status"),
    [
        (
            ValueError("autocorrelation calculation failed"),
            "calculation failed",
            "unavailable",
        ),
        (
            RuntimeError("unexpected autocorrelation failure\nwith context"),
            "RuntimeError",
            "unavailable",
        ),
        (object(), "invalid", "unavailable"),
        (
            mcmc_module.emcee_autocorr.AutocorrError(np.array([1.0])),
            "unreliable",
            "unreliable_short_chain",
        ),
        (np.array([1.0, 1.0]), "invalid shape", "invalid"),
        (np.array([np.nan]), "invalid", "unavailable"),
        (np.array([0.0]), "invalid", "unavailable"),
        (np.array([-1.0]), "invalid", "unavailable"),
        (np.array([3.0]), "does not leave a retained chain", "reliable"),
    ],
    ids=(
        "calculation-failure",
        "unexpected-calculation-failure",
        "malformed",
        "unreliable",
        "wrong-shape",
        "non-finite",
        "zero",
        "negative",
        "no-retained-window",
    ),
)
def test_automatic_mcmc_burn_failure_preserves_only_incomplete_raw_evidence(
    tmp_path: Path,
    autocorrelation_outcome: object,
    failure_reason: str,
    autocorrelation_status: str,
) -> None:
    central_output = tmp_path / "Central"
    failed_output = tmp_path / "Failed"
    method = _automatic_mcmc_method(tmp_path / "method.toml")
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    central_session = AnalysisSession.create()
    failed_session = AnalysisSession.create()
    stale_statistics = failed_output / "Statistics" / "MCMC"
    stale_statistics.mkdir(parents=True)
    for name in (
        "summary.toml",
        "samples.tsv",
        "correlations.tsv",
        "plots.pdf",
        "raw_chain.tsv",
    ):
        (stale_statistics / name).write_text(
            "stale posterior artifact\n", encoding="utf-8"
        )

    run(
        _fit_arguments(central_output, parameters=parameters),
        session=central_session,
    )
    patch_kwargs = (
        {"side_effect": autocorrelation_outcome}
        if isinstance(autocorrelation_outcome, BaseException)
        else {"return_value": autocorrelation_outcome}
    )
    with (
        patch.object(
            mcmc_module.emcee_autocorr,
            "integrated_time",
            **patch_kwargs,
        ),
        pytest.raises(
            NativeMcmcIncompleteError,
            match="authoritative MCMC posterior summarization was withheld",
        ),
    ):
        run(
            _fit_arguments(failed_output, method, parameters=parameters),
            session=failed_session,
        )

    central = central_session.analysis_values.snapshot()
    failed = failed_session.analysis_values.snapshot()
    assert central.revision == failed.revision == 1
    assert central.items() == failed.items()
    statistics = failed_output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["terminal"] == "failed"
    assert diagnostics["completed_steps"] == 5
    assert diagnostics["sampling_terminal"] == "completed"
    assert diagnostics["posterior_retention"] == "failed"
    assert diagnostics["posterior_summary"] == "withheld"
    assert diagnostics["raw_evidence"] == "diagnostic_chain"
    assert diagnostics["raw_chain_file"] == "raw_chain.tsv"
    assert diagnostics["raw_steps"] == 5
    assert diagnostics["raw_samples"] == 160
    assert diagnostics["autocorrelation_status"] == autocorrelation_status
    assert diagnostics["acceptance_fraction_min"] >= 0.0
    assert diagnostics["acceptance_fraction_max"] <= 1.0
    assert (
        diagnostics["acceptance_fraction_min"]
        <= diagnostics["acceptance_fraction_mean"]
        <= diagnostics["acceptance_fraction_max"]
    )
    assert diagnostics["sampling_seconds"] >= 0.0
    assert diagnostics["result_processing_seconds"] >= 0.0
    assert (
        "authoritative MCMC posterior summarization was withheld"
        in diagnostics["failure_message"]
    )
    assert failure_reason in diagnostics["failure_message"]
    raw_chain = (statistics / "raw_chain.tsv").read_text(encoding="utf-8")
    raw_lines = raw_chain.splitlines()
    header = raw_lines[0]
    assert header.startswith("step\twalker\t[R1A_A")
    assert header.endswith("\tlnprob")
    assert len(raw_lines) == 161
    observed_indices = [
        (int(fields[0]), int(fields[1]))
        for fields in (line.split("\t") for line in raw_lines[1:])
    ]
    assert observed_indices == [
        (step, walker) for step in range(1, 6) for walker in range(32)
    ]
    assert all(
        np.all(np.isfinite([float(value) for value in line.split("\t")[2:]]))
        for line in raw_lines[1:]
    )
    assert "stale posterior artifact" not in raw_chain
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()
    outcome = _read_outcome(failed_output)
    assert outcome["latest_committed_revision"] == 1
    assert outcome["restart_revision"] == 1
    assert outcome["failure_stage"] == "statistics"


def test_explicit_mcmc_burn_survives_autocorrelation_calculation_failure(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _nested_mcmc_method(tmp_path / "method.toml", workers=1)
    parameters = _bounded_parameters(tmp_path / "parameters.toml")

    with patch.object(
        mcmc_module.emcee_autocorr,
        "integrated_time",
        side_effect=RuntimeError("autocorrelation backend failed"),
    ):
        run(
            _fit_arguments(output, method, parameters=parameters),
            session=AnalysisSession.create(),
        )

    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["requested_burn"] == 1
    assert diagnostics["autocorrelation_status"] == "unavailable"
    assert "RuntimeError" in diagnostics["autocorrelation_warning"]
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "samples.tsv").is_file()
    assert not (statistics / "raw_chain.tsv").exists()


def test_real_compact_mcmc_fit_is_wholly_native_and_writes_products(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MCMC" = 2')
    parameters = _bounded_parameters(tmp_path / "parameters.toml")
    session = AnalysisSession.create()

    generated_seed = 0x1234_5678_90AB_CDEF
    with (
        patch.object(mcmc_module.secrets, "randbits", return_value=generated_seed),
        patch.object(
            mcmc_module.emcee_autocorr,
            "integrated_time",
            return_value=np.array([0.25]),
        ),
        patch.object(
            mcmc_module,
            "execute_mcmc_evidence",
            wraps=mcmc_module.execute_mcmc_evidence,
        ) as native_sampler,
    ):
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
    assert diagnostics["root_seed"] == generated_seed
    assert "lmfit_version" not in diagnostics
    run_record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert run_record["stochastic_operations"] == [
        {"step": "DEFAULT", "kind": "mcmc", "seed": str(generated_seed)}
    ]
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "# ±" in fitted
    assert "half_credible_interval_68_width" not in fitted
    plan = native_sampler.call_args.args[1]
    assert plan.coordinate_units[0][1] is ParameterUnit.UNSPECIFIED


def test_short_compact_mcmc_form_fails_closed_through_real_chemex_cli(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MCMC" = 2')
    parameters = _bounded_parameters(tmp_path / "parameters.toml")

    with pytest.raises(subprocess.CalledProcessError):
        _run_real_fit_cli(output, method, parameters)

    statistics = output / "Statistics" / "MCMC"
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "incomplete"
    assert diagnostics["engine"] == "native MCMC"
    assert diagnostics["steps"] == 2
    assert diagnostics["raw_evidence"] == "diagnostic_chain"
    assert (statistics / "raw_chain.tsv").is_file()
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "samples.tsv").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


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
    run_record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert run_record["stochastic_operations"] == [
        {"step": "DEFAULT", "kind": "mcmc", "seed": "612"}
    ]
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
    outcome = _read_outcome(output)
    assert outcome["latest_committed_revision"] == 1
    assert outcome["restart_revision"] == 1
    assert outcome["failure_stage"] == "statistics"


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
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "complete",
        "latest_committed_revision": 1,
        "restart_revision": 1,
    }
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
    assert not (output / "run_info" / "restart.toml").exists()
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
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["root_seed"] == 0


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
    diagnostics = tomllib.loads(
        (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    )
    assert diagnostics["root_seed"] == 0


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


def test_real_grid_fit_writes_profiled_surfaces_and_withholds_covariance(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")

    run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    (factor_output,) = tuple((output / "Grid" / "Factors").glob("*.tsv"))
    (profile_output,) = tuple((output / "Grid" / "Profiles" / "1D").glob("*.tsv"))
    factor_rows = np.genfromtxt(
        factor_output, delimiter="\t", names=True, dtype=None, encoding="utf-8"
    )
    profile_rows = np.genfromtxt(
        profile_output, delimiter="\t", names=True, dtype=None, encoding="utf-8"
    )
    summary = tomllib.loads(
        (output / "Grid" / "summary.toml").read_text(encoding="utf-8")
    )
    assert factor_rows["selected"].dtype.kind == "b"
    assert profile_rows["selected"].dtype.kind == "b"
    factor_selected = factor_rows["selected"]
    profile_selected = profile_rows["selected"]
    assert np.count_nonzero(factor_selected) == 1
    assert np.count_nonzero(profile_selected) == 1
    selected_factor = factor_rows[factor_selected]
    selected_profile = profile_rows[profile_selected]
    factor_axis = factor_rows.dtype.names[2]
    profile_axis = profile_rows.dtype.names[0]
    assert selected_factor[factor_axis][0] == pytest.approx(
        selected_profile[profile_axis][0]
    )
    assert summary["selected_axes"][0]["value"] == pytest.approx(
        selected_factor[factor_axis][0]
    )
    assert (output / "Grid" / "grid_1d.pdf").is_file()
    assert not (output / "Grid" / "grid.out").exists()
    assert not (output / "Statistics" / "Covariance" / "evidence.json").exists()
    covariance_status = json.loads(
        (output / "Statistics" / "Covariance" / "status.json").read_text(
            encoding="utf-8"
        )
    )
    assert "discrete profiled surface" in covariance_status["reason"]
    assert (output / "Parameters" / "fixed.toml").is_file()

    stale = output / "Grid" / "Factors" / "stale.tsv"
    stale.write_text("obsolete\n", encoding="utf-8")
    run(_fit_arguments(output, method), session=AnalysisSession.create())
    assert not stale.exists()


def test_evaluation_only_grid_does_not_manufacture_zero_variable_trf(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    method = _grid_method(tmp_path / "method.toml")

    run(
        _fit_arguments(tmp_path / "Output", method),
        session=AnalysisSession.create(),
    )

    rendered = capsys.readouterr().out
    assert "GRID seed" not in rendered
    assert "Evaluations" in rendered


def test_real_grouped_grid_fit_uses_one_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")
    plot_1d = grid_printer_module.plot_grid_1d
    plot_2d = grid_printer_module.plot_grid_2d

    with (
        patch.object(
            grid_printer_module,
            "plot_grid_1d",
            wraps=plot_1d,
        ) as plotted_1d,
        patch.object(
            grid_printer_module,
            "plot_grid_2d",
            wraps=plot_2d,
        ) as plotted_2d,
    ):
        run(
            _fit_arguments(
                output,
                method,
                include=("G2N-HN", "H3N-HN"),
            ),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    factor_outputs = sorted((output / "Grid" / "Factors").glob("*.tsv"))
    assert len(factor_outputs) == 2
    for factor_output in factor_outputs:
        rows = np.genfromtxt(
            factor_output, delimiter="\t", names=True, dtype=None, encoding="utf-8"
        )
        np.testing.assert_allclose(rows[rows.dtype.names[2]], (1.0, 3.0))
        assert tuple(rows["status"]) == ("success", "success")
        assert all(rows["profiles"])
    summary = tomllib.loads(
        (output / "Grid" / "summary.toml").read_text(encoding="utf-8")
    )
    assert summary["factor_count"] == 2
    profiles_1d = tuple((output / "Grid" / "Profiles" / "1D").glob("*.tsv"))
    profiles_2d = tuple((output / "Grid" / "Profiles" / "2D").glob("*.tsv"))
    assert len(profiles_1d) == 2
    assert len(profiles_2d) == 1
    plotted_1d_grids = plotted_1d.call_args.args[0]
    plotted_2d_grids = plotted_2d.call_args.args[0]

    def numerical_surfaces(paths: tuple[Path, ...]) -> list[tuple[float, ...]]:
        surfaces = []
        for profile_path in paths:
            rows = np.genfromtxt(
                profile_path,
                delimiter="\t",
                names=True,
                dtype=None,
                encoding="utf-8",
            )
            surfaces.append(tuple(np.atleast_1d(rows["chi_square"])))
        return sorted(surfaces)

    assert numerical_surfaces(profiles_1d) == sorted(
        tuple(grid.chisqr.ravel()) for grid in plotted_1d_grids
    )
    assert numerical_surfaces(profiles_2d) == sorted(
        tuple(grid.chisqr.ravel()) for grid in plotted_2d_grids
    )
    assert (output / "Parameters" / "fixed.toml").is_file()
    assert not (output / "Grid" / "Groups").exists()
    assert not (output / "All").exists()
    assert not (output / "Groups").exists()
    assert not (output / "Components").exists()


def test_following_direct_step_starts_from_the_committed_grid_solution(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[GRID]
ROLES = [
  { FIX = ["PB", "KEX_AB", "ETAZ_A", "R1_A"] },
  { FIT = ["R1A_A"] },
]

[GRID.SEARCH.GRID]
AXES = ["[R1A_A] = values(1.0, 3.0)"]

[DIRECT]
ROLES_FROM = "GRID"
ROLES = [{ FIT = ["R1A_A"] }]
""",
        encoding="utf-8",
    )
    starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_start(*args, **kwargs):
        starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    session = AnalysisSession.create()
    with patch(
        "chemex.optimize.direct_trf.least_squares",
        side_effect=record_start,
    ):
        run(_fit_arguments(output, method), session=session)

    assert starts[0] == pytest.approx((3.0,))
    assert session.analysis_values.snapshot().revision == 2
    assert (output / "DIRECT" / "Statistics" / "Covariance" / "evidence.json").is_file()


def test_real_v2_de_reaches_normal_trf_product_path_from_out_of_range_start(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "Output"
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()
    trf_starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_full_trf_start(*args, **kwargs):
        trf_starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        side_effect=record_full_trf_start,
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert len(trf_starts) == 1
    starting_parameters = (output / "run_info" / "parameters_used.toml").read_text(
        encoding="utf-8"
    )
    starting_record = next(
        line for line in starting_parameters.splitlines() if '"G2N-H"' in line
    )
    starting_value = float(starting_record.split("[", 1)[1].split(",", 1)[0])
    assert starting_value > 4.0
    assert 1.0 <= trf_starts[0][0] <= 4.0
    assert (output / "Parameters" / "fitted.toml").is_file()
    assert (output / "Data" / "800mhz.dat").is_file()
    assert (output / "Statistics" / "Covariance" / "evidence.json").is_file()
    assert not (output / "Components").exists()
    assert "Running selected-coordinate DE search" in capsys.readouterr().out


def test_real_cli_v2_de_reaches_normal_product_output(tmp_path: Path) -> None:
    output = tmp_path / "Output"
    method = _de_method(tmp_path / "method-v2-de.toml")

    _run_real_fit_cli(output, method, PARAMETERS)

    assert (output / "Parameters" / "fitted.toml").is_file()
    assert (output / "Data" / "800mhz.dat").is_file()
    assert (output / "run_info" / "outcome.toml").is_file()


def test_v2_de_maps_lin_and_log_coordinates_then_releases_every_fit_coordinate(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2-de.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DE]
ROLES = [{ FIT = ["PB", "KEX_AB", "R1A_A"] }]

[DE.SEARCH.DE]
SEED = 597
COORDINATES = [
  "[PB] = log(0.001, 0.1)",
  "[KEX_AB] = lin(100.0, 500.0)",
]
""",
        encoding="utf-8",
    )
    trf_starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_full_trf_start(*args, **kwargs):
        trf_starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=_successful_de_backend_at_x0,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=record_full_trf_start,
        ),
    ):
        run(
            _fit_arguments(tmp_path / "Output", method),
            session=AnalysisSession.create(),
        )

    assert len(trf_starts) == 1
    assert len(trf_starts[0]) == 3
    assert trf_starts[0][:2] == pytest.approx((300.0, 0.01))
    assert trf_starts[0][2] > 4.0


def test_grouped_v2_de_holds_unselected_coordinate_then_releases_all_components(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2-de.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DE]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]

[DE.SEARCH.DE]
SEED = 597
COORDINATES = [
  "[R1A_A, NUC->G2N-H, B0->800MHz] = log(1.0, 4.0)",
]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()
    trf_starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_component_start(*args, **kwargs):
        trf_starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=_successful_de_backend_at_x0,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=record_component_start,
        ),
    ):
        run(
            _fit_arguments(
                tmp_path / "Output",
                method,
                include=("G2N-HN", "H3N-HN"),
            ),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    assert len(trf_starts) == 2
    assert sorted(start[0] for start in trf_starts) == pytest.approx(
        (2.0, 6.87922079444668)
    )


def test_v2_de_failure_has_no_direct_trf_fallback_or_commit(tmp_path: Path) -> None:
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=RuntimeError("DE backend failed"),
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=AssertionError("direct TRF fallback ran"),
        ),
        pytest.raises(RuntimeError, match="no eligible candidate"),
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)

    assert session.analysis_values.snapshot().revision == 0


def test_v2_de_all_invalid_candidates_has_no_trf_fallback_or_commit(
    tmp_path: Path,
) -> None:
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()
    backend_calls = 0

    def invalid_evaluation(evaluator: BoundEvaluator, frame) -> EvaluationFailure:
        return EvaluationFailure(
            evaluator.plan.identity,
            frame.parameterization_identity,
            "kernel",
            "non_finite_calculation",
            "INVALID_TRIAL",
            message="scientifically invalid DE trial",
        )

    def all_invalid_backend(objective, bounds, **kwargs):
        nonlocal backend_calls
        backend_calls += 1
        x0 = np.asarray(kwargs["x0"], dtype=np.float64)
        lower = np.asarray(bounds.lb, dtype=np.float64)
        assert math.isinf(objective(x0))
        assert math.isinf(objective(lower))
        population_size = max(5, int(kwargs["popsize"]) * x0.size)
        return SimpleNamespace(
            success=True,
            message="Optimization terminated successfully.",
            nit=1,
            nfev=2,
            x=x0,
            fun=math.inf,
            population=np.tile(x0, (population_size, 1)),
            population_energies=np.full(population_size, math.inf),
        )

    with (
        patch.object(
            BoundEvaluator,
            "evaluate",
            autospec=True,
            side_effect=invalid_evaluation,
        ),
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=all_invalid_backend,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=AssertionError("direct TRF fallback ran"),
        ),
        pytest.raises(RuntimeError, match="no eligible candidate"),
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)

    assert backend_calls == 1
    assert session.analysis_values.snapshot().revision == 0


def test_v2_de_interruption_leaves_analysis_values_unchanged(tmp_path: Path) -> None:
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=KeyboardInterrupt,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=AssertionError("TRF ran after DE interruption"),
        ),
        pytest.raises(KeyboardInterrupt, match="DE search interrupted"),
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)

    assert session.analysis_values.snapshot().revision == 0


def test_only_final_trf_determines_de_product_values_and_uncertainty(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()

    with patch(
        "chemex.optimize.de_direct_trf.differential_evolution",
        side_effect=_successful_de_backend_at_x0,
    ):
        run(_fit_arguments(output, method), session=session)

    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    fitted_record = next(line for line in fitted.splitlines() if "G2N-H" in line)
    fitted_value = float(fitted_record.split("=", 1)[1].split()[0])
    fitted_error = float(fitted_record.split("±", 1)[1])
    assert fitted_value == pytest.approx(2.34742, rel=5.0e-6)
    assert fitted_value != pytest.approx(2.0)
    assert fitted_error == pytest.approx(0.083366086, rel=3.0e-6)
    assert session.analysis_values.snapshot().revision == 1


def test_v2_de_candidate_cannot_commit_when_final_trf_fails(tmp_path: Path) -> None:
    method = _de_method(tmp_path / "method-v2-de.toml")
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=_successful_de_backend_at_x0,
        ),
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=RuntimeError("final TRF failed"),
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)

    assert session.analysis_values.snapshot().revision == 0


def test_v2_de_seed_repeats_and_changes_the_stochastic_trajectory(
    tmp_path: Path,
) -> None:
    trajectories: list[tuple[float, ...]] = []

    def seeded_backend(objective, bounds, **kwargs):
        lower = np.asarray(bounds.lb, dtype=np.float64)
        upper = np.asarray(bounds.ub, dtype=np.float64)
        vector = lower + kwargs["rng"].random(lower.size) * (upper - lower)
        trajectories.append(tuple(float(value) for value in vector))
        value = objective(vector)
        population_size = max(5, int(kwargs["popsize"]) * vector.size)
        return SimpleNamespace(
            success=True,
            message="Optimization terminated successfully.",
            nit=1,
            nfev=1,
            x=vector,
            fun=value,
            population=np.tile(vector, (population_size, 1)),
            population_energies=np.full(population_size, value),
        )

    with patch(
        "chemex.optimize.de_direct_trf.differential_evolution",
        side_effect=seeded_backend,
    ):
        for ordinal, seed in enumerate((597, 597, 598)):
            method = _de_method(tmp_path / f"method-{ordinal}.toml", seed=seed)
            run(
                _fit_arguments(tmp_path / f"Output-{ordinal}", method),
                session=AnalysisSession.create(),
            )

    assert trajectories[0] == trajectories[1]
    assert trajectories[0] != trajectories[2]
    for ordinal, seed in enumerate((597, 597, 598)):
        run_record = tomllib.loads(
            (tmp_path / f"Output-{ordinal}" / "run_info" / "run.toml").read_text(
                encoding="utf-8"
            )
        )
        assert run_record["stochastic_operations"] == [
            {"step": "DE", "kind": "de", "seed": str(seed)}
        ]


def test_seed_publication_failure_prevents_de_from_starting(tmp_path: Path) -> None:
    output = tmp_path / "Output"
    method = _de_method(tmp_path / "method-v2-de.toml", seed=597)
    real_atomic_write = run_info_module.write_text_atomic

    def fail_run_update(destination: Path, content: str) -> None:
        if destination.name == "run.toml":
            raise OSError("seed publication failed")
        real_atomic_write(destination, content)

    with (
        patch.object(run_info_module, "write_text_atomic", side_effect=fail_run_update),
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=AssertionError("DE consumed RNG after seed publication failed"),
        ) as backend,
        pytest.raises(OSError, match="seed publication failed"),
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    backend.assert_not_called()
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 0,
        "restart_revision": 0,
        "terminal": "failed",
        "failure_stage": "run_info",
        "failure_type": "OSError",
        "failure_message": "seed publication failed",
    }


def test_seed_record_serialization_failure_prevents_de_from_starting(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _de_method(tmp_path / "method-v2-de.toml", seed=597)
    real_json_dumps = run_info_module.json.dumps

    def fail_seed_record_serialization(value, **kwargs) -> str:
        if value == "DE":
            raise RuntimeError("seed record serialization failed")
        return real_json_dumps(value, **kwargs)

    with (
        patch.object(
            run_info_module.json,
            "dumps",
            side_effect=fail_seed_record_serialization,
        ),
        patch(
            "chemex.optimize.de_direct_trf.differential_evolution",
            side_effect=AssertionError("DE started after seed recording failed"),
        ) as backend,
        pytest.raises(RuntimeError, match="seed record serialization failed"),
    ):
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    backend.assert_not_called()
    assert _read_outcome(output) == {
        "schema_version": 2,
        "status": "incomplete",
        "latest_committed_revision": 0,
        "restart_revision": 0,
        "terminal": "failed",
        "failure_stage": "run_info",
        "failure_type": "RuntimeError",
        "failure_message": "seed record serialization failed",
    }


def test_v2_de_statistics_start_only_after_the_final_trf_commit(tmp_path: Path) -> None:
    method = tmp_path / "method-v2-de.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DE]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]
[DE.SEARCH.DE]
SEED = 597
COORDINATES = ["[R1A_A] = log(1.0, 4.0)"]

[DE.STATISTICS.MC]
REPLICATES = 1
SEED = 7
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(tmp_path / "Output", method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    diagnostics = tomllib.loads(
        (
            tmp_path / "Output" / "Statistics" / "MonteCarlo" / "diagnostics.toml"
        ).read_text(encoding="utf-8")
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["root_seed"] == 7


def test_v2_step_after_de_starts_from_committed_final_trf_values(
    tmp_path: Path,
) -> None:
    method = tmp_path / "method-v2-de.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DE]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]
[DE.SEARCH.DE]
SEED = 597
COORDINATES = ["[R1A_A] = log(1.0, 4.0)"]

[FOLLOWUP]
ROLES = [{ FIX = ["PB", "KEX_AB"] }]
""",
        encoding="utf-8",
    )
    starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_start(*args, **kwargs):
        starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        side_effect=record_start,
    ):
        run(
            _fit_arguments(tmp_path / "Output", method),
            session=AnalysisSession.create(),
        )

    first_fitted = (
        tmp_path / "Output" / "DE" / "Parameters" / "fitted.toml"
    ).read_text(encoding="utf-8")
    first_value = float(
        next(line for line in first_fitted.splitlines() if "G2N-H" in line)
        .split("=", 1)[1]
        .split()[0]
    )
    assert len(starts) == 2
    assert starts[1] == pytest.approx((first_value,), rel=5.0e-6)


def test_selected_coordinate_de_reaches_a_better_three_state_dcest_basin(
    tmp_path: Path,
) -> None:
    direct_method = tmp_path / "direct.toml"
    direct_method.write_text(
        """FORMAT_VERSION = 2
[DIRECT]
ROLES = [
  { FIX = ["KEX_BC"] },
  { FIT = ["PB", "PC", "KEX_AB", "KEX_AC", "DW_AB", "DW_AC"] },
]
""",
        encoding="utf-8",
    )
    de_method = tmp_path / "de.toml"
    de_method.write_text(
        """FORMAT_VERSION = 2
[DE]
ROLES = [
  { FIX = ["KEX_BC"] },
  { FIT = ["PB", "PC", "KEX_AB", "KEX_AC", "DW_AB", "DW_AC"] },
]
[DE.SEARCH.DE]
SEED = 597
COORDINATES = [
  "[PB] = log(0.001, 0.2)",
  "[PC] = log(0.001, 0.2)",
  "[KEX_AB] = log(10.0, 5000.0)",
  "[KEX_AC] = log(10.0, 5000.0)",
  "[DW_AB, NUC->K19N] = lin(-15.0, 15.0)",
  "[DW_AC, NUC->K19N] = lin(-15.0, 15.0)",
]
""",
        encoding="utf-8",
    )
    direct_output = tmp_path / "Direct"
    de_output = tmp_path / "DE"

    run(
        _three_state_dcest_arguments(direct_output, direct_method),
        session=AnalysisSession.create(),
    )
    run(
        _three_state_dcest_arguments(de_output, de_method),
        session=AnalysisSession.create(),
    )

    direct_chi_square = tomllib.loads(
        (direct_output / "statistics.toml").read_text(encoding="utf-8")
    )["chi-square"]
    de_chi_square = tomllib.loads(
        (de_output / "statistics.toml").read_text(encoding="utf-8")
    )["chi-square"]
    assert direct_chi_square == pytest.approx(271.894, rel=1.0e-4)
    assert de_chi_square < 0.8 * direct_chi_square


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
    starts: list[tuple[float, ...]] = []
    real_least_squares = direct_trf_module.least_squares

    def record_start(*args, **kwargs):
        starts.append(tuple(float(value) for value in args[1]))
        return real_least_squares(*args, **kwargs)

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        side_effect=record_start,
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 2
    restart = output / "run_info" / "restart.toml"
    assert restart.is_file()
    assert not list(output.glob("*/run_info/restart.toml"))
    assert _read_outcome(output)["restart_revision"] == 2
    first_fitted = (output / "STEP1" / "Parameters" / "fitted.toml").read_text(
        encoding="utf-8"
    )
    first_value = float(
        next(line for line in first_fitted.splitlines() if "G2N-H" in line)
        .split("=", 1)[1]
        .split()[0]
    )
    assert starts[1] == pytest.approx((first_value,), rel=5.0e-6)
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
    assert (output / "run_info" / "restart.toml").is_file()
    outcome = _read_outcome(output)
    assert outcome["status"] == "incomplete"
    assert outcome["latest_committed_revision"] == 1
    assert outcome["restart_revision"] == 1


def test_interrupted_second_step_preserves_first_committed_restart(
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
    real_least_squares = direct_trf_module.least_squares
    call_count = 0

    def interrupt_second_step(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise KeyboardInterrupt
        return real_least_squares(*args, **kwargs)

    session = AnalysisSession.create()
    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=interrupt_second_step,
        ),
        pytest.raises(KeyboardInterrupt),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "run_info" / "restart.toml").is_file()
    outcome = _read_outcome(output)
    assert outcome["terminal"] == "interrupted"
    assert outcome["latest_committed_revision"] == 1
    assert outcome["restart_revision"] == 1


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
    constrained = output / "Parameters" / "constrained.toml"
    assert constrained.is_file()
    constrained_text = constrained.read_text(encoding="utf-8")
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
    assert not (output / "run_info" / "restart.toml").exists()
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


def test_v1_grid_runs_mc_from_the_accepted_aggregate_grid_fit(
    tmp_path: Path,
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

    run(_fit_arguments(output, method), session=session)

    assert (output / "Grid").is_dir()
    assert not (output / "Statistics" / "Covariance" / "evidence.json").exists()
    assert (output / "Statistics" / "Covariance" / "status.json").is_file()
    diagnostics = tomllib.loads(
        (output / "Statistics" / "MonteCarlo" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["status"] == "complete"
    assert diagnostics["root_seed"] == 0
    run_record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert run_record["stochastic_operations"] == [
        {"step": "DEFAULT", "kind": "mc", "seed": "0"}
    ]


def test_v2_grid_runs_requested_statistics_from_the_accepted_grid_fit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DEFAULT]

[DEFAULT.SEARCH.GRID]
AXES = ["[R1A_A] = values(1.0, 3.0)"]

[DEFAULT.STATISTICS.MC]
REPLICATES = 1
SEED = 7
""",
        encoding="utf-8",
    )

    run(_fit_arguments(output, method), session=AnalysisSession.create())

    assert (output / "Grid").is_dir()
    assert (output / "Statistics" / "MonteCarlo" / "samples.tsv").is_file()
    diagnostics = tomllib.loads(
        (output / "Statistics" / "MonteCarlo" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["root_seed"] == 7
    run_record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert run_record["stochastic_operations"] == [
        {"step": "DEFAULT", "kind": "mc", "seed": "7"}
    ]


def test_v2_omitted_resampling_seed_is_generated_once_and_recorded(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method-v2.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[DEFAULT]

[DEFAULT.STATISTICS.MC]
REPLICATES = 1
""",
        encoding="utf-8",
    )
    generated_seed = 0xFEDC_BA09_8765_4321

    with patch(
        "chemex.optimize.fitting.secrets.randbits", return_value=generated_seed
    ) as generate_seed:
        run(_fit_arguments(output, method), session=AnalysisSession.create())

    diagnostics = tomllib.loads(
        (output / "Statistics" / "MonteCarlo" / "diagnostics.toml").read_text(
            encoding="utf-8"
        )
    )
    assert diagnostics["root_seed"] == generated_seed
    assert generate_seed.call_count == 1
    run_record = tomllib.loads(
        (output / "run_info" / "run.toml").read_text(encoding="utf-8")
    )
    assert run_record["stochastic_operations"] == [
        {"step": "DEFAULT", "kind": "mc", "seed": str(generated_seed)}
    ]
