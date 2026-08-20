from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import patch

import pytest

import chemex.optimize.direct_trf as direct_trf_module
from chemex.chemex import run
from chemex.cli import build_parser
from chemex.optimize.resampling import NativeResamplingIncompleteError
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
EXAMPLE = ROOT / "examples/Experiments/RELAXATION_HZNZ"
EXPERIMENT = EXAMPLE / "Experiments/800mhz.toml"
PARAMETERS = EXAMPLE / "Parameters/parameters.toml"
METHOD = EXAMPLE / "Methods/method.toml"


def _fit_arguments(
    output: Path,
    method: Path = METHOD,
    *,
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
            str(PARAMETERS),
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


def test_real_direct_fit_uses_native_trf_and_commits_product_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    with (
        patch(
            "lmfit.Minimizer.minimize",
            side_effect=AssertionError("legacy lmfit optimizer was called"),
        ),
        patch(
            "chemex.containers.experiments.Experiments.residuals",
            side_effect=AssertionError("legacy evaluation was called"),
        ),
    ):
        run(_fit_arguments(output, plot_level="normal"), session=session)

    assert session.analysis_values.snapshot().revision == 1
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    fitted_record = next(line for line in fitted.splitlines() if "=" in line)
    fitted_value = float(fitted_record.split("=", 1)[1].split()[0])
    assert fitted_value == pytest.approx(2.34742, rel=5.0e-6)
    assert "(error not calculated)" in fitted_record
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


def test_real_grouped_direct_fit_uses_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
    ):
        run(
            _fit_arguments(output, include=("G2N-HN", "H3N-HN")),
            session=session,
        )

    assert session.analysis_values.snapshot().revision == 1
    group_outputs = tuple((output / "Groups").glob("*/Parameters/fitted.toml"))
    assert len(group_outputs) == 2
    assert (output / "All" / "Parameters" / "fitted.toml").is_file()


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


def test_real_mc_fit_is_wholly_native_and_writes_product_statistics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"MC" = 2')
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.parameters.database.ParameterStore.build_lmfit_params",
            side_effect=AssertionError("legacy lmfit Parameters were constructed"),
        ),
        patch(
            "lmfit.Minimizer.minimize",
            side_effect=AssertionError("legacy lmfit optimizer was called"),
        ),
        patch(
            "chemex.containers.experiments.Experiments.residuals",
            side_effect=AssertionError("legacy evaluation was called"),
        ),
        patch(
            "chemex.optimize.fitting.run_resampling_statistics",
            side_effect=AssertionError("legacy resampling was called"),
        ),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    statistics = output / "Statistics" / "MonteCarlo"
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0
    samples = (statistics / "samples.tsv").read_text(encoding="utf-8")
    assert len(samples.splitlines()) == 3
    diagnostics = (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    assert "requested_samples = 2" in diagnostics
    assert "completed_samples = 2" in diagnostics
    assert 'engine = "native direct TRF"' in diagnostics
    assert "root_seed = 0" in diagnostics
    assert 'status = "complete"' in diagnostics
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "(error not calculated)" in fitted


def test_native_statistics_filter_resolution_fails_closed_before_lmfit(
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
        patch("chemex.containers.experiments.Experiments.filter") as legacy_filter,
        pytest.raises(RuntimeError, match="native values unavailable"),
    ):
        run(_fit_arguments(output, method), session=session)

    legacy_filter.assert_not_called()


def test_real_bs_fit_uses_native_refits_and_writes_bootstrap_products(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"BS" = 2')
    session = AnalysisSession.create()

    with (
        patch(
            "lmfit.Minimizer.minimize",
            side_effect=AssertionError("legacy lmfit optimizer was called"),
        ),
        patch(
            "chemex.containers.experiments.Experiments.residuals",
            side_effect=AssertionError("legacy evaluation was called"),
        ),
        patch(
            "chemex.optimize.fitting.run_resampling_statistics",
            side_effect=AssertionError("legacy resampling was called"),
        ),
    ):
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


def test_real_bsn_fit_uses_native_nucleus_resampling_products(tmp_path: Path) -> None:
    output = tmp_path / "Output"
    method = _statistics_method(tmp_path / "method.toml", '"BSN" = 2')
    session = AnalysisSession.create()

    with (
        patch(
            "lmfit.Minimizer.minimize",
            side_effect=AssertionError("legacy lmfit optimizer was called"),
        ),
        patch(
            "chemex.containers.experiments.Experiments.residuals",
            side_effect=AssertionError("legacy evaluation was called"),
        ),
        patch(
            "chemex.optimize.fitting.run_resampling_statistics",
            side_effect=AssertionError("legacy resampling was called"),
        ),
    ):
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
    assert 'failure_type = "RuntimeError"' in diagnostics
    assert 'failure_message = "manifest unavailable"' in diagnostics
    assert not (statistics / "summary.toml").exists()
    assert not (statistics / "correlations.tsv").exists()
    assert not (statistics / "plots.pdf").exists()


def test_real_grid_fit_uses_native_cartesian_trf_and_writes_grid_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Grid" / "grid.out").is_file()
    assert (output / "Parameters" / "fixed.toml").is_file()


def test_real_grouped_grid_fit_uses_one_native_aggregate_commit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()
    method = _grid_method(tmp_path / "method.toml")

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
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

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 2
    for step in ("STEP1", "STEP2"):
        fitted = (output / step / "Parameters" / "fitted.toml").read_text(
            encoding="utf-8"
        )
        assert "(error not calculated)" in fitted
        fixed = (output / step / "Parameters" / "fixed.toml").read_text(
            encoding="utf-8"
        )
        assert "PB" in fixed
        assert "KEX_AB" in fixed


def test_native_step_clears_generic_error_from_preceding_legacy_step(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[LEGACY]
FITMETHOD = "leastsq"
FIX = ["PB", "KEX_AB"]

[NATIVE]
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    run(_fit_arguments(output, method), session=session)

    legacy = (output / "LEGACY" / "Parameters" / "fitted.toml").read_text(
        encoding="utf-8"
    )
    native = (output / "NATIVE" / "Parameters" / "fitted.toml").read_text(
        encoding="utf-8"
    )
    assert "±" in legacy
    assert "(error not calculated)" in native


def test_native_backend_failure_cannot_commit_or_publish_fitted_output(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    session = AnalysisSession.create()

    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=RuntimeError("backend failed"),
        ),
        patch(
            "lmfit.Minimizer.minimize",
            side_effect=AssertionError("legacy lmfit optimizer was called"),
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(output), session=session)

    assert session.analysis_values.snapshot().revision == 0
    assert not (output / "Parameters").exists()
    assert not (output / "Data").exists()


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

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
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

    with patch(
        "lmfit.Minimizer.minimize",
        side_effect=AssertionError("legacy lmfit optimizer was called"),
    ):
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

    with patch(
        "chemex.optimize.fitting._fit_groups",
        side_effect=AssertionError("legacy grouped fit was called"),
    ):
        run(_fit_arguments(output, method), session=session)

    assert session.analysis_values.snapshot().revision == 1
    assert (output / "Statistics" / "MonteCarlo" / "samples.tsv").is_file()


def test_grid_with_statistics_remains_wholly_on_legacy_dispatch(
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

    with (
        patch(
            "chemex.optimize.fitting.run_native_deterministic",
            side_effect=AssertionError("native deterministic dispatch was called"),
        ),
        patch("chemex.optimize.fitting.run_grid") as run_grid,
    ):
        run(_fit_arguments(output, method), session=session)

    run_grid.assert_called_once()
