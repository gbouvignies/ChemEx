from pathlib import Path
from types import SimpleNamespace

import numpy as np

from chemex.containers.data import Data
from chemex.plotters import cest as cest_module
from chemex.plotters import cpmg as cpmg_module


def _make_calc_data() -> Data:
    data = Data(
        exp=np.array([0.0]),
        err=np.array([0.0]),
        metadata=np.array([1.0]),
    )
    data.calc = np.array([0.5])
    return data


def test_cest_plot_simulation_accepts_empty_experimental_data(
    monkeypatch,
    tmp_path: Path,
) -> None:
    recorded: dict[str, Data] = {}
    plotter = cest_module.CestPlotter(Path("13hz.toml"), SimpleNamespace())

    monkeypatch.setattr(cest_module, "print_plot_filename", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cest_module, "create_plot_data_calc", lambda _profile: _make_calc_data())

    def fake_plot_profile(_pdf, _profile, data_exp: Data, data_calc: Data) -> None:
        recorded["data_exp"] = data_exp
        recorded["data_calc"] = data_calc

    monkeypatch.setattr(plotter, "_plot_profile", fake_plot_profile)

    plotter.plot_simulation(tmp_path, [SimpleNamespace(spin_system="G23N-H")])

    assert recorded["data_exp"].size == 0
    assert recorded["data_exp"].mask.size == 0
    assert recorded["data_calc"].calc.tolist() == [0.5]
    assert "[G23N-H]" in (tmp_path / "13hz.sim").read_text(encoding="utf-8")


def test_cpmg_plot_simulation_accepts_empty_experimental_data(
    monkeypatch,
    tmp_path: Path,
) -> None:
    recorded: dict[str, Data] = {}
    plotter = cpmg_module.CpmgPlotter(Path("500mhz.toml"), SimpleNamespace())

    monkeypatch.setattr(cpmg_module, "print_plot_filename", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        cpmg_module,
        "create_plot_data_calc",
        lambda _profile, _config: _make_calc_data(),
    )

    def fake_plot_cpmg(_pdf, _name: str, data_exp: Data, data_calc: Data) -> None:
        recorded["data_exp"] = data_exp
        recorded["data_calc"] = data_calc

    monkeypatch.setattr(cpmg_module, "plot_cpmg", fake_plot_cpmg)

    plotter.plot_simulation(tmp_path, [SimpleNamespace(spin_system="G23N-H")])

    assert recorded["data_exp"].size == 0
    assert recorded["data_exp"].mask.size == 0
    assert recorded["data_calc"].calc.tolist() == [0.5]
    assert "[G23N-H]" in (tmp_path / "500mhz.sim").read_text(encoding="utf-8")
