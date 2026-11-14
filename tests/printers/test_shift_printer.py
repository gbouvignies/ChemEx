from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.chemex import run
from chemex.cli import build_parser
from chemex.configuration.methods import Selection
from chemex.containers.data import Data
from chemex.experiments import experiment_types
from chemex.experiments.loader import register_experiments
from chemex.plotters.shift import ShiftPlotter
from chemex.printers.data import ShiftPrinter
from chemex.runtime.session import AnalysisSession

ROOT = Path(__file__).resolve().parents[2]
ALL_PROFILES = Selection(include=None, exclude=None)


def test_shift_printer_units_precision_and_simulation_toggle() -> None:
    data = Data(
        exp=np.array([14.844195]),
        err=np.array([0.000500]),
    )
    data.calc = np.array([14.844212])

    ppb_printer = ShiftPrinter()
    assert "# Unit: PPB\n" in ppb_printer.header
    assert "SHIFT (EXP)" in ppb_printer.header
    assert "ERROR (EXP)" in ppb_printer.header
    assert "SHIFT (CALC)" in ppb_printer.header
    assert (
        ppb_printer.print("G23N-HN", data)
        == f"{'G23N-HN':>18s} {14.844195: 17.3f} {0.000500: 17.3f} {14.844212: 17.3f}\n"
    )

    ppm_printer = ShiftPrinter(unit="PPM")
    assert "# Unit: PPM\n" in ppm_printer.header
    assert (
        ppm_printer.print("A5CB", data)
        == f"{'A5CB':>18s} {14.844195: 17.6f} {0.000500: 17.6f} {14.844212: 17.6f}\n"
    )

    # Mutating `simulation = True` after construction matches Profile.prepare_for_simulation()
    ppm_printer.simulation = True
    assert "# Unit: PPM\n" in ppm_printer.header
    assert "SHIFT (EXP)" not in ppm_printer.header
    assert "ERROR (EXP)" not in ppm_printer.header
    assert "SHIFT (CALC)" in ppm_printer.header
    sim_row = ppm_printer.print("A5CB", data)
    assert sim_row == f"{'A5CB':>18s} {14.844212: 17.6f}\n"
    assert sim_row.split() == ["A5CB", "14.844212"]

    sim_ppb = ShiftPrinter(simulation=True, unit="PPB")
    assert "# Unit: PPB\n" in sim_ppb.header
    assert "SHIFT (EXP)" not in sim_ppb.header
    assert sim_ppb.print("G23N-HN", data).split() == ["G23N-HN", "14.844"]


def test_shift_printer_rejects_invalid_unit_or_precision() -> None:
    with pytest.raises(ValueError, match="Unsupported shift unit"):
        ShiftPrinter(unit="HZ")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="`precision` must be non-negative"):
        ShiftPrinter(precision=-1)


@pytest.mark.parametrize(
    ("exp_name", "expected_unit"),
    [
        ("shift_15n_sqmq", "PPB"),
        ("shift_15n_sq", "PPM"),
        ("shift_13c_sq", "PPM"),
    ],
)
def test_shift_experiment_types_configure_expected_printer_units(
    tmp_path: Path,
    exp_name: str,
    expected_unit: str,
) -> None:
    register_experiments()
    data_file = tmp_path / "shifts.txt"
    data_file.write_text("G23N-HN  1.0  0.1\nA24N-HN  2.0  0.1\n", encoding="utf-8")
    exp_toml = tmp_path / f"{exp_name}.toml"
    exp_toml.write_text(
        "\n".join(
            [
                "[experiment]",
                f'name = "{exp_name}"',
                "",
                "[conditions]",
                "h_larmor_frq = 600.0",
                "",
                "[data]",
                'path = "shifts.txt"',
                "",
            ]
        ),
        encoding="utf-8",
    )

    session = AnalysisSession.create()
    session.set_model("2st")
    source = experiment_types.open(exp_toml)
    built = experiment_types.build(
        source,
        selection=ALL_PROFILES,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )
    assert isinstance(built.experiment.printer, ShiftPrinter)
    assert built.experiment.printer.unit == expected_unit
    assert isinstance(built.experiment.plotter, ShiftPlotter)
    assert built.experiment.plotter.unit == expected_unit


def test_shift_simulate_and_fit_end_to_end_cli(tmp_path: Path) -> None:
    # 1. Simulate `shift_13c_sq` and verify `# Unit: PPM`, `.6f`, and 2-column rows
    c13_data = tmp_path / "c13_input.txt"
    c13_data.write_text("A5CB   0.000000  0.000500\nL8CD1  0.000000  0.000500\n")
    c13_exp = tmp_path / "c13_sq.toml"
    c13_exp.write_text(
        "\n".join(
            [
                "[experiment]",
                'name = "shift_13c_sq"',
                "",
                "[conditions]",
                "h_larmor_frq = 600.0",
                "",
                "[data]",
                'path = "c13_input.txt"',
                "",
            ]
        )
    )
    c13_pars = tmp_path / "c13_pars.toml"
    c13_pars.write_text(
        "\n".join(
            [
                "[GLOBAL]",
                "PB = 0.10",
                "KEX_AB = 5000.0",
                "R2_A = 10.0",
                "",
                "[CS_A]",
                "A5CB = 0.050000",
                "L8CD1 = 15.000000",
                "",
                "[DW_AB]",
                "A5CB = -2.000000",
                "L8CD1 = -2.000000",
                "",
            ]
        )
    )
    sim_out = tmp_path / "out_sim"
    parser = build_parser()
    sim_args = parser.parse_args(
        [
            "simulate",
            "-e",
            str(c13_exp),
            "-p",
            str(c13_pars),
            "-d",
            "2st",
            "-o",
            str(sim_out),
            "--plot",
            "nothing",
        ]
    )
    session = AnalysisSession.create()
    run(sim_args, session=session)

    sim_dat = (sim_out / "Data" / "c13_sq.dat").read_text(encoding="utf-8")
    assert "# Unit: PPM\n" in sim_dat
    assert "SHIFT (EXP)" not in sim_dat
    assert "ERROR (EXP)" not in sim_dat
    assert "SHIFT (CALC)" in sim_dat
    sim_rows = [
        line.split()
        for line in sim_dat.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert len(sim_rows) == 2
    assert all(len(row) == 2 for row in sim_rows)
    sim_map = {name: float(val_str) for name, val_str in sim_rows}
    # Near 0.05 ppm with DW_AB = -2.0 ppm, shift must move negative (< 0), no mirror bounce!
    assert sim_map["A5CB"] < 0.0
    assert (sim_map["A5CB"] - 0.05) == pytest.approx(sim_map["L8CD1"] - 15.0, abs=2e-6)

    # 2. Fit shipped `examples/Combinations/Shifts` (500 + 800 MHz SQMQ) WITHOUT CS_A!
    sqmq_500 = ROOT / "examples/Combinations/Shifts/Experiments/sqmq_500.toml"
    sqmq_800 = ROOT / "examples/Combinations/Shifts/Experiments/sqmq_800.toml"
    sqmq_pars = tmp_path / "sqmq_no_csa_pars.toml"
    sqmq_pars.write_text(
        "\n".join(
            [
                "[GLOBAL]",
                "PB = 0.023",
                "KEX_AB = 4000.0",
                "DW_AB = 0.0",
                "",
                "[DW_AB]",
                "29N = -4.0",
                "29HN = 0.70",
                "117N = 1.0",
                "117HN = 0.70",
                "",
            ]
        )
    )
    fit_out = tmp_path / "out_fit"
    fit_args = parser.parse_args(
        [
            "fit",
            "-e",
            str(sqmq_500),
            str(sqmq_800),
            "-p",
            str(sqmq_pars),
            "-d",
            "2st",
            "--include",
            "29",
            "117",
            "-o",
            str(fit_out),
            "--plot",
            "nothing",
        ]
    )
    fit_session = AnalysisSession.create()
    run(fit_args, session=fit_session)

    fit_dat = (fit_out / "Data" / "sqmq_500.dat").read_text(encoding="utf-8")
    assert "# Unit: PPB\n" in fit_dat
    assert "SHIFT (EXP)" in fit_dat
    assert "ERROR (EXP)" in fit_dat
    assert "SHIFT (CALC)" in fit_dat
    fit_rows = {
        parts[0]: [float(x) for x in parts[1:]]
        for line in fit_dat.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
        for parts in [line.split()]
    }
    assert len(fit_rows["29N-H"]) == 3
    # Non-zero SQMQ calculated shift even with CS_A omitted, matching experimental signs
    assert fit_rows["29N-H"][2] < -10.0
    assert fit_rows["117N-H"][2] > 5.0
