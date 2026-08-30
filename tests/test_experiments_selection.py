from dataclasses import dataclass, field
from pathlib import Path

import pytest

from chemex.configuration.method_plan import ProfileSelection
from chemex.configuration.methods import Selection
from chemex.containers.experiments import Experiments


@dataclass
class SelectableExperiment:
    filename: Path
    profile_count: int
    canonical_selections: list[ProfileSelection] = field(default_factory=list)

    def select(self, _selection: Selection) -> None:
        self.profile_count = 0

    def __len__(self) -> int:
        return self.profile_count

    def __bool__(self) -> bool:
        return self.profile_count > 0

    def select_profiles(self, selection: ProfileSelection) -> None:
        self.canonical_selections.append(selection)


def test_experiments_become_false_when_selection_removes_all_profiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    experiments = Experiments()
    experiments.add(SelectableExperiment(Path("a.toml"), 2))
    experiments.add(SelectableExperiment(Path("b.toml"), 1))

    monkeypatch.setattr(
        "chemex.containers.experiments.print_selecting_profiles",
        lambda _selected_nb: None,
    )

    assert experiments
    assert len(experiments) == 3

    experiments.select(Selection(include=[], exclude=None))

    assert len(experiments) == 0
    assert not experiments


def test_canonical_empty_selection_is_applied_to_restore_step_local_baseline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    experiments = Experiments()
    experiment = SelectableExperiment(Path("a.toml"), 2)
    experiments.add(experiment)
    monkeypatch.setattr(
        "chemex.containers.experiments.print_selecting_profiles",
        lambda _selected_nb: None,
    )

    experiments.select_profiles(ProfileSelection())

    assert experiment.canonical_selections == [ProfileSelection()]
