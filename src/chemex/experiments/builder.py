from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from pydantic import ValidationError
from rich.live import Live

from chemex.configuration.base import BaseSettings, ExperimentConfiguration
from chemex.configuration.methods import Selection
from chemex.containers.dataset import Dataset
from chemex.containers.experiment import Experiment
from chemex.containers.experiments import Experiments
from chemex.containers.profile import Profile
from chemex.experiments.factories import Creators, factories
from chemex.messages import (
    console,
    get_reading_exp_text,
    print_experiment_name_error,
    print_file_not_found_error,
    print_loading_experiments,
    print_pydantic_parsing_error,
)
from chemex.models.model import ModelSpec
from chemex.parameters.factory import create_parameters
from chemex.runtime import AnalysisSession
from chemex.toml import read_toml

GenericConfig = ExperimentConfiguration[Any, Any, Any]


class NameConfig(BaseSettings):
    name: str


class ExperimentConfig(BaseSettings):
    experiment: NameConfig


def _get_experiment_name(config: Mapping[str, object], filename: Path) -> str:
    try:
        model = ExperimentConfig.model_validate(config)
    except ValidationError:
        print_experiment_name_error(filename)
        sys.exit()

    return model.experiment.name


def _apply_selection(dataset: Dataset, selection: Selection) -> Dataset:
    if (include := selection.include) is not None:
        return [
            (spin_system, data)
            for spin_system, data in dataset
            if spin_system.part_of(include)
        ]

    if (exclude := selection.exclude) is not None:
        return [
            (spin_system, data)
            for spin_system, data in dataset
            if not spin_system.part_of(exclude)
        ]

    return dataset


def _create_dataset(
    filename: Path,
    live: Live,
    factory: Creators,
    config: GenericConfig,
) -> Dataset:
    try:
        dataset = factory.create_dataset(filename.parent, config)
    except FileNotFoundError as e:
        live.stop()
        print_file_not_found_error(e)
        sys.exit(1)
    return dataset


def _create_config(
    filename: Path,
    live: Live,
    factory: Creators,
    config_dict: Mapping[str, object],
    *,
    model_spec: ModelSpec,
) -> GenericConfig:
    try:
        config = factory.create_config(config_dict, model=model_spec)
    except ValidationError as e:
        live.stop()
        print_pydantic_parsing_error(filename, e)
        sys.exit()
    return config


def build_experiment(
    filename: Path,
    selection: Selection,
    *,
    session: AnalysisSession,
) -> Experiment:
    config_dict: Mapping[str, object] = read_toml(filename)

    experiment_name = _get_experiment_name(config_dict, filename)
    parameter_store = session.parameters
    parameter_factory = session.parameter_factory
    model_spec = session.model.spec

    with Live(
        get_reading_exp_text(filename, experiment_name),
        console=console,
    ) as live:
        factory = factories.get(experiment_name)

        config = _create_config(
            filename,
            live,
            factory,
            config_dict,
            model_spec=model_spec,
        )

        printer = factory.create_printer()
        plotter = factory.create_plotter(filename, config)

        dataset = _create_dataset(filename, live, factory, config)
        dataset = _apply_selection(dataset, selection)

        profiles: list[Profile] = []

        for spin_system, data in dataset:
            spectrometer = factory.create_spectrometer(config, spin_system)
            sequence = factory.create_sequence(config)
            filterer = factory.create_filterer(config, spectrometer)
            name_map = create_parameters(
                config,
                spectrometer.liouvillian,
                parameter_factory=parameter_factory,
            )
            profiles.append(
                Profile(
                    data,
                    spectrometer,
                    sequence,
                    name_map,
                    printer,
                    filterer,
                    config.data.scaled,
                ),
            )

        live.update(get_reading_exp_text(filename, experiment_name, len(profiles)))

    experiment = Experiment(
        filename,
        experiment_name,
        profiles,
        printer,
        plotter,
        parameter_store=parameter_store,
    )
    experiment.estimate_noise(config.data.error, global_error=config.data.global_error)

    return experiment


def build_experiments(
    filenames: list[Path] | None,
    selection: Selection,
    *,
    session: AnalysisSession,
) -> Experiments:
    parameter_store = session.parameters

    if not filenames:
        return Experiments(parameter_store=parameter_store)

    print_loading_experiments()

    experiments = Experiments(parameter_store=parameter_store)

    for filename in filenames:
        experiment = build_experiment(filename, selection, session=session)
        experiments.add(experiment)

    parameter_store.sort_parameters()

    return experiments
