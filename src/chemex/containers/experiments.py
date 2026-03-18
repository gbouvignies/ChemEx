"""Module for managing and manipulating collections of NMR experiments.

This module contains the Experiments class and functions for generating different
types of simulated experiment datasets for statistical analysis.
"""

from collections.abc import Iterator
from hashlib import blake2b
from itertools import chain
from pathlib import Path
from random import choices
from typing import Literal, Self

import numpy as np
from lmfit.parameter import Parameters

from chemex.configuration.methods import Selection
from chemex.containers.experiment import Experiment
from chemex.messages import print_selecting_profiles
from chemex.parameters.database import ParameterStore
from chemex.parameters.spin_system import Group, SpinSystem
from chemex.typing import Array

# Type definitions
SelectionType = list[SpinSystem] | Literal["*", "all"], None


def _sanitize_output_part(part: str) -> str:
    replacements = {
        "": "__empty__",
        ".": "__current__",
        "..": "__parent__",
    }
    return replacements.get(part, part)


def _filename_parts_for_output(filename: Path) -> tuple[str, ...]:
    output_path = filename.with_suffix("")
    parts = output_path.parts

    if output_path.anchor:
        parts = parts[1:]
        parts = ("__absolute__", *parts)

    return tuple(_sanitize_output_part(part) for part in parts)


def _suffix_path(parts: tuple[str, ...], depth: int) -> Path:
    return Path(*parts[-depth:])


def _append_hash_suffix(path: Path, filename: Path) -> Path:
    digest = blake2b(str(filename).encode("utf-8"), digest_size=4).hexdigest()
    return path.parent / f"{path.name}__{digest}"


def _build_output_stems(filenames: list[Path]) -> dict[Path, Path]:
    output_parts = {filename: _filename_parts_for_output(filename) for filename in filenames}
    depth_by_filename = dict.fromkeys(filenames, 1)

    while True:
        grouped: dict[Path, list[Path]] = {}
        for filename, parts in output_parts.items():
            suffix = _suffix_path(parts, depth_by_filename[filename])
            grouped.setdefault(suffix, []).append(filename)

        duplicates = [
            filenames_for_suffix
            for filenames_for_suffix in grouped.values()
            if len(filenames_for_suffix) > 1
        ]
        if not duplicates:
            break

        progress = False
        for filenames_for_suffix in duplicates:
            for filename in filenames_for_suffix:
                max_depth = len(output_parts[filename])
                if depth_by_filename[filename] < max_depth:
                    depth_by_filename[filename] += 1
                    progress = True

        if not progress:
            colliding_filenames = {
                filename
                for filenames_for_suffix in duplicates
                for filename in filenames_for_suffix
            }
            return {
                filename: (
                    _append_hash_suffix(
                        _suffix_path(output_parts[filename], depth_by_filename[filename]),
                        filename,
                    )
                    if filename in colliding_filenames
                    else _suffix_path(output_parts[filename], depth_by_filename[filename])
                )
                for filename in filenames
            }

    return {
        filename: _suffix_path(output_parts[filename], depth_by_filename[filename])
        for filename in filenames
    }


class Experiments:
    """A collection of NMR experiments.

    This class manages a set of NMR experiments, providing functionalities for
    adding experiments, calculating residuals, simulation preparation, writing
    data, plotting, and more.
    """

    def __init__(self, *, parameter_store: ParameterStore) -> None:
        """Initialize an empty Experiments collection."""
        self._experiments: dict[Path, Experiment] = {}
        self.parameter_store = parameter_store

    @property
    def groups(self) -> set[Group]:
        """Get the set of unique groups across all experiments.

        Returns:
            set[Group]: A set of unique groups.

        """
        return {group for experiment in self for group in experiment.groups}

    def add(self, experiment: Experiment) -> None:
        """Add an experiment to the collection.

        Args:
            experiment (Experiment): The experiment to be added.

        """
        if experiment.filename in self._experiments:
            msg = "Experiment already exists."
            raise ValueError(msg)
        self._experiments[experiment.filename] = experiment

    def residuals(self, params: Parameters) -> Array:
        """Calculate the residuals for all experiments in the collection.

        Args:
            params (Parameters): Parameters for residual calculation.

        Returns:
            Array: Residuals as a NumPy array.

        """
        return np.concatenate([experiment.residuals(params) for experiment in self])

    def back_calculate(self) -> None:
        """Back calculate experiments using the instance parameter store."""
        params_lf = self.parameter_store.build_lmfit_params(self.param_ids)
        self.residuals(params_lf)

    def prepare_for_simulation(self) -> None:
        """Prepare each experiment in the collection for simulation."""
        self.back_calculate()
        for experiment in self:
            experiment.prepare_for_simulation()

    def write(self, path: Path) -> None:
        """Write experiment data to a specified path.

        Args:
            path (Path): Directory path to write data.

        """
        path_dat = path / "Data"
        path_dat.mkdir(parents=True, exist_ok=True)
        output_stems = _build_output_stems(list(self._experiments))
        for experiment in self:
            output_stem = path_dat / output_stems[experiment.filename]
            output_stem.parent.mkdir(parents=True, exist_ok=True)
            experiment.write(output_stem)

    def plot(self, path: Path) -> None:
        """Plot each experiment in the collection.

        Args:
            path (Path): Directory path for saving plots.

        """
        output_stems = _build_output_stems(list(self._experiments))
        for experiment in self:
            output_stem = path / output_stems[experiment.filename]
            output_stem.parent.mkdir(parents=True, exist_ok=True)
            experiment.plot(output_stem)

    def plot_simulation(self, path: Path) -> None:
        """Plot simulations for each experiment.

        Args:
            path (Path): Directory path for saving simulation plots.

        """
        output_stems = _build_output_stems(list(self._experiments))
        for experiment in self:
            output_stem = path / output_stems[experiment.filename]
            output_stem.parent.mkdir(parents=True, exist_ok=True)
            experiment.plot_simulation(output_stem)

    def select(self, selection: Selection) -> None:
        """Select experiments based on the given selection criteria.

        Args:
            selection (Selection): Criteria for selecting experiments.

        """
        if selection.include is None and selection.exclude is None:
            return
        for experiment in self:
            experiment.select(selection)
        print_selecting_profiles(len(self))

    @property
    def param_id_sets(self) -> list[set[str]]:
        """Get a list of parameter ID sets for all experiments.

        Returns:
            list[set[str]]: List of sets containing parameter IDs.

        """
        return list(
            chain.from_iterable(experiment.param_id_sets for experiment in self),
        )

    @property
    def param_ids(self) -> set[str]:
        """Get a set of all unique parameter IDs across experiments.

        Returns:
            set[str]: Set of unique parameter IDs.

        """
        result: set[str] = set()
        return result.union(*(self.param_id_sets))

    def filter(self) -> None:
        """Get a set of all unique parameter IDs across experiments.

        Returns:
            set[str]: Set of unique parameter IDs.

        """
        params = self.parameter_store.build_lmfit_params(self.param_ids)
        for experiment in self:
            experiment.filter(params)

    def get_relevant_subset(self, param_ids: set[str]) -> Self:
        """Get a subset of experiments relevant to specified parameter IDs.

        Args:
            param_ids (set[str]): Parameter IDs to filter experiments.

        Returns:
            Self: Subset of Experiments relevant to given parameter IDs.

        """
        relevant_subset = type(self)(parameter_store=self.parameter_store)
        for experiment in self:
            if subset := experiment.get_relevant_subset(param_ids):
                relevant_subset.add(subset)
        return relevant_subset

    def __iter__(self) -> Iterator[Experiment]:
        """Iterate over the experiments in the collection.

        Returns:
            Iterator[Experiment]: An iterator over experiments.

        """
        return iter(self._experiments.values())

    def __len__(self) -> int:
        """Get the total number of data points in the collection.

        Returns:
            int: Total number of data points.

        """
        return sum(len(experiment) for experiment in self)

    def __bool__(self) -> bool:
        """Check if the collection contains any active profiles.

        Returns:
            bool: True if at least one profile remains selected, False otherwise.

        """
        return len(self) > 0


def generate_monte_carlo_experiments(experiments: Experiments) -> Experiments:
    """Generate a new Experiments collection for Monte Carlo simulation.

    Args:
        experiments (Experiments): Original experiments.

    Returns:
        Experiments: Collection with Monte Carlo simulated data.

    """
    experiments_mc = Experiments(parameter_store=experiments.parameter_store)
    for experiment in experiments:
        experiments_mc.add(experiment.monte_carlo())
    return experiments_mc


def generate_bootstrap_experiments(experiments: Experiments) -> Experiments:
    """Generate a new Experiments collection for Bootstrap simulation.

    Args:
        experiments (Experiments): Original experiments.

    Returns:
        Experiments: Collection with Bootstrap simulated data.

    """
    experiments_mc = Experiments(parameter_store=experiments.parameter_store)
    for experiment in experiments:
        experiments_mc.add(experiment.bootstrap())
    return experiments_mc


def generate_bootstrap_ns_experiments(experiments: Experiments) -> Experiments:
    """Generate a new Experiments collection for Bootstrap NS simulation.

    Args:
        experiments (Experiments): Original experiments.

    Returns:
        Experiments: Collection with Bootstrap NS simulated data.

    """
    groups = experiments.groups
    groups_bs = choices(tuple(groups), k=len(groups))
    experiments_bs = Experiments(parameter_store=experiments.parameter_store)
    for experiment in experiments:
        experiments_bs.add(experiment.bootstrap_ns(groups_bs))
    return experiments_bs


def generate_exp_for_statistics(
    experiments: Experiments,
    statistic_name: str,
) -> Experiments:
    """Generate experiments for a specific statistical method.

    Args:
        experiments (Experiments): Original experiments.
        statistic_name (str): Name of the statistical method.

    Returns:
        Experiments: Collection with simulated data for statistics.

    """
    generators = {
        "mc": generate_monte_carlo_experiments,
        "bs": generate_bootstrap_experiments,
        "bsn": generate_bootstrap_ns_experiments,
    }
    return generators[statistic_name](experiments)
