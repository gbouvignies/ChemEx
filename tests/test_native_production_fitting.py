from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import patch

import pytest

from chemex.chemex import run
from chemex.cli import build_parser
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
            "1",
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


def test_statistics_fallback_maps_explicit_trf_to_legacy_least_squares(
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

    with patch("chemex.optimize.fitting._fit_groups") as fit_groups:
        run(_fit_arguments(output, method), session=session)

    fit_groups.assert_called_once()
    assert fit_groups.call_args.args[3] == "least_squares"


def test_statistics_fallback_preserves_omitted_legacy_fitmethod(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = tmp_path / "method.toml"
    method.write_text(
        """[DEFAULT]
STATISTICS = { "MC" = 1 }
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()

    with patch("chemex.optimize.fitting._fit_groups") as fit_groups:
        run(_fit_arguments(output, method), session=session)

    fit_groups.assert_called_once()
    assert fit_groups.call_args.args[3] == "leastsq"


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
