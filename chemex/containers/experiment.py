from __future__ import annotations

from dataclasses import dataclass, field
from itertools import chain
from typing import TYPE_CHECKING

import numpy as np

from chemex.messages import (
    print_no_duplicate_warning,
    print_not_implemented_noise_method_warning,
)
from chemex.parameters import database
from chemex.uncertainty import estimate_noise_variance

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from lmfit import Parameters as ParametersLF

    from chemex.configuration.methods import Selection
    from chemex.containers.profile import Profile
    from chemex.parameters.spin_system import Group
    from chemex.plotters.plotter import Plotter
    from chemex.printers.data import Printer


@dataclass
class Experiment:
    filename: Path
    name: str
    profiles: list[Profile]
    filtered_profiles: list[Profile] = field(init=False, default_factory=list)
    printer: Printer
    plotter: Plotter

    def residuals(self, params: ParametersLF) -> list[float]:
        return list(
            chain.from_iterable(profile.residuals(params) for profile in self.profiles)
        )

    def plot(self, path: Path):
        self.plotter.plot(path, self.profiles)

    def plot_simulation(self, path: Path):
        self.plotter.plot_simulation(path, self.profiles)

    def write(self, path: Path):
        filename = (path / self.filename.name).with_suffix(".dat")
        with filename.open("w") as file_dat:
            file_dat.write(self.printer.header)
            for profile in sorted(self.profiles):
                file_dat.write(str(profile))

    def select(self, selection: Selection):
        include = selection.include
        exclude = selection.exclude
        profiles_all = [*self.profiles, *self.filtered_profiles]
        profiles: list[Profile] = []
        filtered: list[Profile] = []
        for profile in profiles_all:
            included = include is None or profile.name.part_of(include)
            excluded = exclude is not None and profile.name.part_of(exclude)
            if included and not excluded:
                profiles.append(profile)
            else:
                filtered.append(profile)
        self.profiles = profiles
        self.filtered_profiles = filtered

    def filter(self, params: ParametersLF):
        for profile in self.profiles:
            profile.filter(params)

    def _any_duplicate(self):
        return any(profile.any_duplicate() for profile in self.profiles)

    def estimate_noise(self, kind: str):
        # TODO: Validation should be moved to the configuration file module
        implemented = ("file", "scatter", "duplicates")
        if kind not in implemented:
            print_not_implemented_noise_method_warning(self.filename, kind, implemented)
            kind = "file"
        if kind == "duplicates" and not self._any_duplicate():
            print_no_duplicate_warning(self.filename)
            kind = "file"
        if kind == "file" or not self.profiles:
            return
        noise_variance_values = [
            estimate_noise_variance[kind](profile.data) for profile in self.profiles
        ]
        noise_mean = np.sqrt(np.mean(noise_variance_values))
        for profile in self.profiles:
            profile.set_noise(noise_mean)

    def prepare_for_simulation(self) -> None:
        for profile in self.profiles:
            profile.prepare_for_simulation()

    def monte_carlo(self) -> Experiment:
        profiles = [profile.monte_carlo() for profile in self.profiles]
        return Experiment(
            self.filename, self.name, profiles, self.printer, self.plotter
        )

    def bootstrap(self) -> Experiment:
        profiles = [profile.bootstrap() for profile in self.profiles]
        return Experiment(
            self.filename, self.name, profiles, self.printer, self.plotter
        )

    def bootstrap_ns(self, groups: list[Group]) -> Experiment:
        """Residue-specific bootstrap."""
        profiles: dict[Group, list[Profile]] = {}
        for profile in self.profiles:
            profiles.setdefault(profile.name.groups["i"], []).append(profile)
        profiles_bs_ns: list[Profile] = []
        for group in groups:
            profiles_bs_ns.extend(profiles.get(group, []))
        return Experiment(
            self.filename, self.name, profiles_bs_ns, self.printer, self.plotter
        )

    @property
    def groups(self) -> set[Group]:
        return {profile.name.groups["i"] for profile in self.profiles}

    @property
    def param_id_sets(self) -> list[set[str]]:
        return [
            set(database.get_parameters(profile.param_ids)) for profile in self.profiles
        ]

    def get_relevant_subset(self, param_ids: set[str]) -> Experiment:
        profiles = [
            profile
            for profile in self.profiles
            if set(database.get_parameters(profile.param_ids)) & param_ids
        ]

        return Experiment(
            self.filename, self.name, profiles, self.printer, self.plotter
        )

    def __iter__(self) -> Iterator[Profile]:
        yield from self.profiles

    def __len__(self) -> int:
        return len(self.profiles)

    def __bool__(self) -> bool:
        return bool(len(self.profiles))
