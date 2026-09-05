"""Outer orchestration for loading configured Experiments."""

from __future__ import annotations

from pathlib import Path

from rich.live import Live

from chemex.configuration.methods import Selection
from chemex.containers.experiment import (
    NoDuplicateNoiseNotice,
    UnsupportedNoiseMethodNotice,
)
from chemex.containers.experiments import Experiments
from chemex.experiments import experiment_types
from chemex.experiments.experiment_types import ExperimentNotice
from chemex.messages import (
    console,
    get_reading_exp_text,
    print_loading_experiments,
    print_no_duplicate_warning,
    print_not_implemented_noise_method_warning,
)
from chemex.runtime import AnalysisSession


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
        source = experiment_types.open(filename)

        with Live(
            get_reading_exp_text(filename, source.experiment_type_name),
            console=console,
        ) as live:
            result = experiment_types.build(
                source,
                selection=selection,
                model=session.model.spec,
                parameters=session.parameter_factory,
            )

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
