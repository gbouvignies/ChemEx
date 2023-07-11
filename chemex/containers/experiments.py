from __future__ import annotations

from itertools import chain
from random import choices
from typing import TYPE_CHECKING, Literal

import numpy as np

from chemex.messages import print_selecting_profiles
from chemex.parameters import database
from chemex.parameters.spin_system import Group, SpinSystem

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from lmfit.parameter import Parameters

    from chemex.configuration.methods import Selection
    from chemex.containers.experiment import Experiment
    from chemex.typing import ArrayFloat


# Type definitions
SelectionType = list[SpinSystem] | Literal["*", "all"], None


class Experiments:
    def __init__(self) -> None:
        self._experiments: dict[Path, Experiment] = {}

    @property
    def groups(self) -> set[Group]:
        groups: set[Group] = set()
        return groups.union(*(experiment.groups for experiment in self))

    def add(self, experiment: Experiment):
        self._experiments[experiment.filename] = experiment

    def residuals(self, params: Parameters) -> ArrayFloat:
        return np.asarray(
            list(
                chain.from_iterable(experiment.residuals(params) for experiment in self)
            )
        )

    def back_calculate(self) -> None:
        params_lf = database.build_lmfit_params(self.param_ids)
        self.residuals(params_lf)

    def prepare_for_simulation(self) -> None:
        self.back_calculate()
        for experiment in self:
            experiment.prepare_for_simulation()

    def write(self, path: Path):
        path_dat = path / "Data"
        path_dat.mkdir(parents=True, exist_ok=True)
        for experiment in self:
            experiment.write(path_dat)

    def plot(self, path: Path):
        for experiment in self:
            experiment.plot(path)

    def plot_simulation(self, path: Path):
        for experiment in self:
            experiment.plot_simulation(path)

    def select(self, selection: Selection):
        if selection.include is None and selection.exclude is None:
            return
        for experiment in self:
            experiment.select(selection)
        print_selecting_profiles(len(self))

    @property
    def param_id_sets(self) -> list[set[str]]:
        return list(
            chain.from_iterable(experiment.param_id_sets for experiment in self)
        )

    @property
    def param_ids(self) -> set[str]:
        result: set[str] = set()
        return result.union(*(self.param_id_sets))

    def filter(self) -> None:
        params = database.build_lmfit_params(self.param_ids)
        for experiment in self:
            experiment.filter(params)

    def get_relevant_subset(self, param_ids: set[str]) -> Experiments:
        relevant_subset = Experiments()
        for experiment in self:
            if subset := experiment.get_relevant_subset(param_ids):
                relevant_subset.add(subset)
        return relevant_subset

    def __iter__(self) -> Iterator[Experiment]:
        yield from self._experiments.values()

    def __len__(self) -> int:
        return sum(len(experiment) for experiment in self)

    def __bool__(self) -> bool:
        return bool(len(self))


def generate_monte_carlo_experiments(experiments: Experiments) -> Experiments:
    experiments_mc = Experiments()
    for experiment in experiments:
        experiments_mc.add(experiment.monte_carlo())
    return experiments_mc


def generate_bootstrap_experiments(experiments: Experiments) -> Experiments:
    experiments_mc = Experiments()
    for experiment in experiments:
        experiments_mc.add(experiment.bootstrap())
    return experiments_mc


def generate_bootstrap_ns_experiments(experiments: Experiments) -> Experiments:
    groups = experiments.groups
    groups_bs = choices(tuple(groups), k=len(groups))
    experiments_bs = Experiments()
    for experiment in experiments:
        experiments_bs.add(experiment.bootstrap_ns(groups_bs))
    return experiments_bs


def generate_exp_for_statistics(
    experiments: Experiments, statistic_name: str
) -> Experiments:
    generators = {
        "mc": generate_monte_carlo_experiments,
        "bs": generate_bootstrap_experiments,
        "bsn": generate_bootstrap_ns_experiments,
    }
    return generators[statistic_name](experiments)
