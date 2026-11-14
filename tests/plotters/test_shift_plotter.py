from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure

from chemex.containers.data import Data
from chemex.plotters import shift as shift_module
from chemex.plotters.shift import ShiftPlotter, plot_shift
from chemex.printers.data import ShiftUnit


@pytest.mark.parametrize(
    ("unit", "expected_suffix"),
    [
        ("PPB", "(ppb)"),
        ("PPM", "(ppm)"),
    ],
)
def test_plot_shift_labels_axes_with_configured_unit(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    unit: ShiftUnit,
    expected_suffix: str,
) -> None:
    captured_labels: dict[str, str] = {}

    def fake_savefig(
        self: Figure, fname: str, *_args: object, **_kwargs: object
    ) -> None:
        ax = self.axes[0]
        captured_labels["xlabel"] = ax.get_xlabel()
        captured_labels["ylabel"] = ax.get_ylabel()
        captured_labels["fname"] = fname

    monkeypatch.setattr(Figure, "savefig", fake_savefig)

    data = Data(exp=np.array([1.0, 2.0]), err=np.array([0.1, 0.1]))
    data.calc = np.array([1.05, 1.95])
    pdf_path = tmp_path / "shifts.pdf"
    plot_shift(pdf_path, data, unit=unit)

    assert captured_labels["fname"] == str(pdf_path)
    assert captured_labels["xlabel"] == rf"$Δδ_\mathregular{{fit}}$ {expected_suffix}"
    assert captured_labels["ylabel"] == rf"$Δδ_\mathregular{{exp}}$ {expected_suffix}"


def test_shift_plotter_forwards_unit_for_plot_and_simulation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        shift_module,
        "print_plot_filename",
        lambda *_args, **_kwargs: None,
    )
    recorded_units: list[ShiftUnit] = []
    real_plot_shift = shift_module.plot_shift

    def spy_plot_shift(
        name_pdf: Path,
        data: Data,
        *,
        unit: ShiftUnit = "PPB",
    ) -> None:
        recorded_units.append(unit)
        real_plot_shift(name_pdf, data, unit=unit)

    monkeypatch.setattr(shift_module, "plot_shift", spy_plot_shift)

    profile_data = Data(exp=np.array([15.0]), err=np.array([0.001]))
    profile_data.calc = np.array([14.998])
    profiles = [SimpleNamespace(data=profile_data)]

    plotter_ppm = ShiftPlotter(Path("c13_sq.toml"), unit="PPM")
    stem_fit = tmp_path / "c13_fit"
    stem_sim = tmp_path / "c13_sim"
    plotter_ppm.plot(stem_fit, profiles)  # type: ignore[arg-type]
    plotter_ppm.plot_simulation(stem_sim, profiles)  # type: ignore[arg-type]

    assert recorded_units == ["PPM", "PPM"]
    assert stem_fit.with_suffix(".pdf").stat().st_size > 0
    assert stem_sim.with_suffix(".pdf").stat().st_size > 0


def test_shift_plotter_rejects_invalid_unit(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Unsupported shift unit"):
        ShiftPlotter(Path("shifts.toml"), unit="HZ")  # type: ignore[arg-type]

    data = Data(exp=np.array([1.0]), err=np.array([0.1]))
    data.calc = np.array([1.0])
    with pytest.raises(ValueError, match="Unsupported shift unit"):
        plot_shift(tmp_path / "shifts.pdf", data, unit="HZ")  # type: ignore[arg-type]
