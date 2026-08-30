from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Self

import numpy as np

from chemex.configuration.method_plan import ProfileSelection
from chemex.configuration.methods import Selection
from chemex.containers.profile import Profile
from chemex.parameters.spin_system import Group, SpinSystem
from chemex.plotters.plotter import Plotter
from chemex.printers.data import Printer
from chemex.uncertainty import estimate_noise_variance


@dataclass(frozen=True)
class UnsupportedNoiseMethodNotice:
    """Report that noise estimation fell back to errors from the data file."""

    filename: Path
    requested: str
    implemented: tuple[str, ...]


@dataclass(frozen=True)
class NoDuplicateNoiseNotice:
    """Report that duplicate-based noise estimation had no duplicate data."""

    filename: Path


type NoiseEstimateNotice = UnsupportedNoiseMethodNotice | NoDuplicateNoiseNotice


@dataclass
class Experiment:
    filename: Path
    name: str
    profiles: list[Profile]
    filtered_profiles: list[Profile] = field(init=False, default_factory=list)
    printer: Printer
    plotter: Plotter
    _noise_notices: tuple[NoiseEstimateNotice, ...] = field(
        init=False,
        repr=False,
        compare=False,
        default=(),
    )

    def back_calculate_from_values(
        self,
        parameter_values: Mapping[str, float],
    ) -> None:
        """Publish calculated profile values from a native resolved mapping."""
        for profile in self.profiles:
            profile.calculate_from_values(parameter_values)

    def plot(self, output_stem: Path) -> None:
        self.plotter.plot(output_stem, self.profiles)

    def plot_simulation(self, output_stem: Path) -> None:
        self.plotter.plot_simulation(output_stem, self.profiles)

    def write(self, output_stem: Path) -> None:
        filename = output_stem.with_suffix(".dat")
        with filename.open("w", encoding="utf-8") as file_dat:
            file_dat.write(self.printer.header)
            for profile in sorted(self.profiles):
                file_dat.write(str(profile))

    def _select_profiles(
        self,
        include: list[SpinSystem] | tuple[SpinSystem, ...] | str | None,
        exclude: list[SpinSystem] | tuple[SpinSystem, ...] | str | None,
    ) -> None:
        profiles_all = [*self.profiles, *self.filtered_profiles]
        profiles: list[Profile] = []
        filtered: list[Profile] = []
        for profile in profiles_all:
            included = include is None or profile.spin_system.part_of(include)
            excluded = exclude is not None and profile.spin_system.part_of(exclude)
            if included and not excluded:
                profiles.append(profile)
            else:
                filtered.append(profile)
        self.profiles = profiles
        self.filtered_profiles = filtered

    def select(self, selection: Selection) -> None:
        """Apply the legacy standalone selection interface."""
        self._select_profiles(selection.include, selection.exclude)

    def select_profiles(self, selection: ProfileSelection) -> None:
        """Apply canonical step-local profile selection semantics."""

        def resolve(
            value: tuple[str, ...] | str | None,
        ) -> tuple[SpinSystem, ...] | str | None:
            if value is None or isinstance(value, str):
                return value
            return tuple(SpinSystem.from_name(name) for name in value)

        self._select_profiles(
            resolve(selection.include),
            resolve(selection.exclude),
        )

    def filter_from_values(self, parameter_values: Mapping[str, float]) -> None:
        """Apply profile filters from resolved native parameter values."""
        for profile in self.profiles:
            profile.filter_from_values(parameter_values)

    def _any_duplicate(self) -> bool:
        return any(profile.any_duplicate() for profile in self.profiles)

    def estimate_noise(
        self,
        kind: str,
        *,
        global_error: bool = True,
    ) -> None:
        """Estimate data noise while preserving the historical ``None`` return."""
        notices: list[NoiseEstimateNotice] = []
        # TODO: Validation should be moved to the configuration file module
        implemented = ("file", "scatter", "duplicates")
        if kind not in implemented:
            notices.append(
                UnsupportedNoiseMethodNotice(self.filename, kind, implemented),
            )
            kind = "file"
        if kind == "duplicates" and not self._any_duplicate():
            notices.append(NoDuplicateNoiseNotice(self.filename))
            kind = "file"
        if kind == "file" or not self.profiles:
            self._noise_notices = tuple(notices)
            return
        if global_error:
            noise_variance_values = [
                estimate_noise_variance[kind](profile.data) for profile in self.profiles
            ]
            noise_mean = np.sqrt(np.mean(noise_variance_values))
            for profile in self.profiles:
                profile.set_noise(noise_mean)
        else:
            for profile in self.profiles:
                profile.set_noise(np.sqrt(estimate_noise_variance[kind](profile.data)))
        self._noise_notices = tuple(notices)

    @property
    def noise_notices(self) -> tuple[NoiseEstimateNotice, ...]:
        """Nonfatal notices produced by the latest noise estimate."""
        return self._noise_notices

    def prepare_for_simulation(self) -> None:
        for profile in self.profiles:
            profile.prepare_for_simulation()

    def monte_carlo(self) -> Self:
        profiles = [profile.monte_carlo() for profile in self.profiles]
        return type(self)(
            self.filename,
            self.name,
            profiles,
            self.printer,
            self.plotter,
        )

    def bootstrap(self) -> Self:
        profiles = [profile.bootstrap() for profile in self.profiles]
        return type(self)(
            self.filename,
            self.name,
            profiles,
            self.printer,
            self.plotter,
        )

    def bootstrap_ns(self, groups: list[Group]) -> Self:
        """Residue-specific bootstrap."""
        profiles: dict[Group, list[Profile]] = {}
        for profile in self.profiles:
            profiles.setdefault(profile.spin_system.groups["i"], []).append(profile)
        profiles_bs_ns: list[Profile] = []
        for group in groups:
            profiles_bs_ns.extend(profiles.get(group, []))
        return type(self)(
            self.filename,
            self.name,
            profiles_bs_ns,
            self.printer,
            self.plotter,
        )

    @property
    def groups(self) -> set[Group]:
        return {profile.spin_system.groups["i"] for profile in self.profiles}

    @property
    def _required_param_id_sets(self) -> list[set[str]]:
        """Return direct profile requirements without expanding constraints."""
        return [set(profile.param_ids) for profile in self.profiles]

    def __iter__(self) -> Iterator[Profile]:
        yield from self.profiles

    def __len__(self) -> int:
        return len(self.profiles)

    def __bool__(self) -> bool:
        return bool(len(self.profiles))
