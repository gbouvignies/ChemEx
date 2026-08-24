"""Outer orchestration for loading configured Experiments."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NoReturn

from rich.live import Live

from chemex.configuration.methods import Selection
from chemex.containers.dataset import DatasetLoadError
from chemex.containers.experiment import (
    NoDuplicateNoiseNotice,
    UnsupportedNoiseMethodNotice,
)
from chemex.containers.experiments import Experiments
from chemex.experiments import experiment_types
from chemex.experiments.experiment_types import (
    ExperimentBuildError,
    ExperimentConfigurationError,
    ExperimentDataError,
    ExperimentFileError,
    ExperimentNameError,
    ExperimentNotice,
    ExperimentTomlError,
    InvalidExperimentSourceError,
    UnknownExperimentTypeError,
)
from chemex.messages import (
    console,
    get_reading_exp_text,
    print_dataset_error,
    print_experiment_name_error,
    print_experiment_source_error,
    print_experiment_type_error,
    print_file_not_found,
    print_file_not_found_error,
    print_loading_experiments,
    print_no_duplicate_warning,
    print_not_implemented_noise_method_warning,
    print_pydantic_parsing_error,
    print_toml_error,
)
from chemex.runtime import AnalysisSession


def _exit_after_build_error(error: ExperimentBuildError) -> NoReturn:
    if isinstance(error, ExperimentFileError):
        print_file_not_found(error.filename)
    elif isinstance(error, ExperimentTomlError):
        print_toml_error(error.filename, error.error)
    elif isinstance(error, ExperimentNameError):
        print_experiment_name_error(error.filename)
    elif isinstance(error, UnknownExperimentTypeError):
        print_experiment_type_error(error.filename, error.experiment_type_name)
    elif isinstance(error, InvalidExperimentSourceError):
        print_experiment_source_error(error.filename, error.explanation)
    elif isinstance(error, ExperimentConfigurationError):
        print_pydantic_parsing_error(error.filename, error.error)
    elif isinstance(error, ExperimentDataError):
        if isinstance(error.error, DatasetLoadError):
            print_dataset_error(error.error.filename, error.error.explanation)
        else:
            print_file_not_found_error(error.error)
    else:
        raise error
    sys.exit(1)


def _print_notice(notice: ExperimentNotice) -> None:
    if isinstance(notice, UnsupportedNoiseMethodNotice):
        print_not_implemented_noise_method_warning(
            notice.filename,
            notice.requested,
            notice.implemented,
        )
    elif isinstance(notice, NoDuplicateNoiseNotice):
        print_no_duplicate_warning(notice.filename)


def build_experiments(
    filenames: list[Path] | None,
    selection: Selection,
    *,
    session: AnalysisSession,
) -> Experiments:
    """Build all configured Experiments using one session's parameter state."""
    if not filenames:
        session.parameter_factory.try_seal_definitions()
        return Experiments()

    print_loading_experiments()
    experiments = Experiments()

    for filename in filenames:
        try:
            source = experiment_types.open(filename)
        except ExperimentBuildError as error:
            _exit_after_build_error(error)

        with Live(
            get_reading_exp_text(filename, source.experiment_type_name),
            console=console,
        ) as live:
            try:
                result = experiment_types.build(
                    source,
                    selection=selection,
                    model=session.model.spec,
                    parameters=session.parameter_factory,
                )
            except ExperimentBuildError as error:
                live.stop()
                _exit_after_build_error(error)

            live.update(
                get_reading_exp_text(
                    filename,
                    source.experiment_type_name,
                    len(result.experiment),
                ),
            )

        for notice in result.notices:
            _print_notice(notice)
        experiments.add(result.experiment)

    session.parameters.sort()
    session.parameter_factory.try_seal_definitions()
    return experiments
