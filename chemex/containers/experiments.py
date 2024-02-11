"""Module for managing and manipulating collections of NMR experiments.

This module contains the Experiments class and functions for generating different
types of simulated experiment datasets for statistical analysis.
"""
from __future__ import annotations

from collections.abc import Iterator
from itertools import chain
from pathlib import Path
from random import choices
from typing import Literal, Self

import numpy as np
from lmfit.parameter import Parameters

from chemex.configuration.methods import Selection
from chemex.containers.experiment import Experiment
from chemex.messages import print_selecting_profiles
from chemex.parameters import database
from chemex.parameters.spin_system import Group, SpinSystem
from chemex.typing import ArrayFloat

# Type definitions
SelectionType = list[SpinSystem] | Literal["*", "all"], None


class Experiments:
    """A collection of NMR experiments.

    This class manages a set of NMR experiments, providing functionalities for
    adding experiments, calculating residuals, simulation preparation, writing
    data, plotting, and more.
    """

    def __init__(self) -> None:
        """Initialize an empty Experiments collection."""
        self._experiments: dict[Path, Experiment] = {}

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

    def residuals(self, params: Parameters) -> ArrayFloat:
        """Calculate the residuals for all experiments in the collection.

        Args:
            params (Parameters): Parameters for residual calculation.

        Returns:
            ArrayFloat: Residuals as a NumPy array.
        """
        return np.asarray(
            list(
                chain.from_iterable(
                    experiment.residuals(params) for experiment in self
                ),
            ),
        )

    def back_calculate(self) -> None:
        """Back calculate experiments' parameters from the database."""
        params_lf = database.build_lmfit_params(self.param_ids)
        self.residuals(params_lf)

    def prepare_for_simulation(self) -> None:
        """Prepare each experiment in the collection for simulation."""
        self.back_calculate()
        for experiment in self:
            experiment.prepare_for_simulation()

    def write(self, path: Path):
        """Write experiment data to a specified path.

        Args:
            path (Path): Directory path to write data.
        """
        path_dat = path / "Data"
        path_dat.mkdir(parents=True, exist_ok=True)
        for experiment in self:
            experiment.write(path_dat)

    def plot(self, path: Path):
        """Plot each experiment in the collection.

        Args:
            path (Path): Directory path for saving plots.
        """
        for experiment in self:
            experiment.plot(path)

    def plot_simulation(self, path: Path):
        """Plot simulations for each experiment.

        Args:
            path (Path): Directory path for saving simulation plots.
        """
        for experiment in self:
            experiment.plot_simulation(path)

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
        params = database.build_lmfit_params(self.param_ids)
        for experiment in self:
            experiment.filter(params)

    def get_relevant_subset(self, param_ids: set[str]) -> Self:
        """Get a subset of experiments relevant to specified parameter IDs.

        Args:
            param_ids (set[str]): Parameter IDs to filter experiments.

        Returns:
            Self: Subset of Experiments relevant to given parameter IDs.
        """
        relevant_subset = type(self)()
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
        """Check if the collection contains any experiments.

        Returns:
            bool: True if collection is non-empty, False otherwise.
        """
        return bool(self._experiments)


def generate_monte_carlo_experiments(experiments: Experiments) -> Experiments:
    """Generate a new Experiments collection for Monte Carlo simulation.

    Args:
        experiments (Experiments): Original experiments.

    Returns:
        Experiments: Collection with Monte Carlo simulated data.
    """
    experiments_mc = Experiments()
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
    experiments_mc = Experiments()
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
    experiments_bs = Experiments()
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
