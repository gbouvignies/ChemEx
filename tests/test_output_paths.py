from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from chemex.containers.experiments import Experiments


@dataclass
class RecordingExperiment:
    filename: Path
    written_stems: list[Path] = field(default_factory=list)
    plotted_stems: list[Path] = field(default_factory=list)
    simulated_stems: list[Path] = field(default_factory=list)

    def write(self, output_stem: Path) -> None:
        self.written_stems.append(output_stem)
        output_stem.with_suffix(".dat").write_text("data", encoding="utf-8")

    def plot(self, output_stem: Path) -> None:
        self.plotted_stems.append(output_stem)
        output_stem.with_suffix(".pdf").write_text("plot", encoding="utf-8")

    def plot_simulation(self, output_stem: Path) -> None:
        self.simulated_stems.append(output_stem)
        output_stem.with_suffix(".sim").write_text("simulation", encoding="utf-8")

    def __len__(self) -> int:
        return 1

    def __bool__(self) -> bool:
        return True


def test_experiments_write_and_plot_use_unique_output_stems(tmp_path: Path) -> None:
    experiments = Experiments(parameter_store=object())  # type: ignore[arg-type]
    experiment_a = RecordingExperiment(Path("set_a/experiment.toml"))
    experiment_b = RecordingExperiment(Path("set_b/experiment.toml"))
    experiments.add(experiment_a)
    experiments.add(experiment_b)

    experiments.write(tmp_path)
    experiments.plot(tmp_path / "Plots")
    experiments.plot_simulation(tmp_path / "SimPlots")

    assert experiment_a.written_stems == [tmp_path / "Data" / "set_a" / "experiment"]
    assert experiment_b.written_stems == [tmp_path / "Data" / "set_b" / "experiment"]
    assert experiment_a.plotted_stems == [tmp_path / "Plots" / "set_a" / "experiment"]
    assert experiment_b.plotted_stems == [tmp_path / "Plots" / "set_b" / "experiment"]
    assert experiment_a.simulated_stems == [
        tmp_path / "SimPlots" / "set_a" / "experiment"
    ]
    assert experiment_b.simulated_stems == [
        tmp_path / "SimPlots" / "set_b" / "experiment"
    ]

    assert (tmp_path / "Data" / "set_a" / "experiment.dat").exists()
    assert (tmp_path / "Data" / "set_b" / "experiment.dat").exists()
    assert (tmp_path / "Plots" / "set_a" / "experiment.pdf").exists()
    assert (tmp_path / "Plots" / "set_b" / "experiment.pdf").exists()
    assert (tmp_path / "SimPlots" / "set_a" / "experiment.sim").exists()
    assert (tmp_path / "SimPlots" / "set_b" / "experiment.sim").exists()


def test_experiments_keep_simple_stems_when_basenames_are_unique(tmp_path: Path) -> None:
    experiments = Experiments(parameter_store=object())  # type: ignore[arg-type]
    experiment_a = RecordingExperiment(Path("set_a/alpha.toml"))
    experiment_b = RecordingExperiment(Path("set_b/beta.toml"))
    experiments.add(experiment_a)
    experiments.add(experiment_b)

    experiments.write(tmp_path)

    assert experiment_a.written_stems == [tmp_path / "Data" / "alpha"]
    assert experiment_b.written_stems == [tmp_path / "Data" / "beta"]


def test_experiments_fall_back_to_hashed_stems_for_path_sentinel_collisions(
    tmp_path: Path,
) -> None:
    experiments = Experiments(parameter_store=object())  # type: ignore[arg-type]
    experiment_a = RecordingExperiment(Path("__absolute__/set/experiment.toml"))
    experiment_b = RecordingExperiment(Path("/set/experiment.toml"))
    experiments.add(experiment_a)
    experiments.add(experiment_b)

    experiments.write(tmp_path)

    stem_a = experiment_a.written_stems[0]
    stem_b = experiment_b.written_stems[0]

    assert stem_a != stem_b
    assert stem_a.name.startswith("experiment__")
    assert stem_b.name.startswith("experiment__")
