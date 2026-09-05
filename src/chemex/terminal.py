"""Presentation-only normalization for trusted terminal failures."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from chemex.exceptions import ChemExError


def _fields(headline: str, *fields: tuple[str, str | None]) -> str:
    lines = [headline]
    lines.extend(f"{label}: {value}" for label, value in fields if value)
    return "\n".join(lines)


def _paths(error: BaseException) -> tuple[Path, ...]:
    paths: list[Path] = []
    for name in (
        "diagnostics_path",
        "evidence_path",
        "samples_path",
        "failures_path",
        "restart_path",
        "outcome_path",
    ):
        value = getattr(error, name, None)
        if isinstance(value, Path) and value not in paths:
            paths.append(value)
    return tuple(paths)


def _preserved_fields(error: BaseException) -> tuple[tuple[str, str | None], ...]:
    paths = _paths(error)
    return (("Preserved", ", ".join(str(path) for path in paths)),) if paths else ()


def _analysis_fields(error: BaseException) -> tuple[tuple[str, str | None], ...]:
    stage = getattr(error, "failure_stage", None)
    step = getattr(error, "method_step", None)
    return (
        ("During", stage.replace("_", " ") if isinstance(stage, str) else None),
        ("Method Step", step if isinstance(step, str) else None),
    )


def format_interruption(error: BaseException) -> str:
    """Format a user interruption from verified orchestration facts only."""
    stage = getattr(error, "failure_stage", None)
    step = getattr(error, "method_step", None)
    context = None
    if isinstance(stage, str):
        context = stage.replace("_", " ")
    if isinstance(step, str):
        context = f"{context or 'analysis'} (Method Step {step})"
    lines = ["Analysis interrupted by user."]
    if context:
        lines.append(f"During: {context}")
    paths = _paths(error)
    if paths:
        lines.append("Preserved:")
        lines.extend(f"  {path}" for path in paths)
    return "\n".join(lines)


def _validation_reason(error: object) -> str:
    errors = getattr(error, "errors", None)
    if not callable(errors):
        return "The configuration does not match the required schema."
    items = errors()
    if not items:
        return "The configuration does not match the required schema."
    first = items[0]
    location = " -> ".join(str(item) for item in first.get("loc", ()))
    message = str(first.get("msg", "Invalid value"))
    return f"{location}: {message}" if location else message


type _Formatter = Callable[[ChemExError], str | None]


def _format_path_error(error: ChemExError) -> str | None:
    from chemex.exceptions import ArtifactPublicationError, InputFileReadError

    if isinstance(error, InputFileReadError):
        return _fields(
            "Input file could not be read.",
            ("Input", str(error.path)),
            ("Reason", str(error.error)),
            ("Next", "Check the input path and read permissions."),
        )
    if isinstance(error, ArtifactPublicationError):
        return _fields(
            "ChemEx could not publish an output artifact.",
            ("Operation", error.operation),
            ("Path", str(error.path)),
            ("Reason", str(error.error)),
            *_analysis_fields(error),
            *_preserved_fields(error),
            ("Next", "Check available space and write permissions."),
        )
    return None


def _format_toml(error: ChemExError) -> str | None:
    from chemex.toml import TomlReadError

    if not isinstance(error, TomlReadError):
        return None
    if isinstance(error.error, OSError):
        return _fields(
            "Input file could not be read.",
            ("Input", str(error.filename)),
            ("Reason", str(error.error)),
            ("Next", "Check the input path and read permissions."),
        )
    return _fields(
        "TOML input is invalid.",
        ("Input", str(error.filename)),
        ("Reason", str(error.error)),
        ("Next", "Check the file against the TOML format."),
    )


def _format_method(error: ChemExError) -> str | None:
    from chemex.configuration.method_plan import MethodFormatError

    if not isinstance(error, MethodFormatError):
        return None
    return _fields(
        "Method Plan is invalid.",
        ("Input", str(error.source.filename)),
        ("Method Step", error.source.step),
        ("Field", error.source.field),
        ("Reason", error.message),
    )


def _format_parameters(error: ChemExError) -> str | None:
    from chemex.configuration.parameters import ParameterConfigurationError

    if not isinstance(error, ParameterConfigurationError):
        return None
    return _fields(
        "Parameter configuration is invalid.",
        ("Input", ", ".join(str(path) for path in error.filenames)),
        ("Reason", error.explanation),
    )


def _format_model(error: ChemExError) -> str | None:
    from chemex.models.model import ModelSelectionError

    if not isinstance(error, ModelSelectionError):
        return None
    return _fields(
        "Exchange model selection is invalid.",
        ("Model", error.name),
        ("Reason", error.explanation),
        ("Available", ", ".join(error.available)),
    )


def _format_experiment(error: ChemExError) -> str | None:
    from chemex.containers.dataset import DatasetLoadError
    from chemex.experiments.experiment_types import (
        ExperimentConfigurationError,
        ExperimentDataError,
        ExperimentFileError,
        ExperimentNameError,
        ExperimentTomlError,
        UnknownExperimentTypeError,
    )

    if isinstance(error, ExperimentFileError):
        return _fields(
            "Experiment input could not be read.",
            ("Input", str(error.filename)),
            ("Reason", str(error.error)),
            ("Next", "Check the input path and read permissions."),
        )
    if isinstance(error, ExperimentTomlError):
        return _fields(
            "Experiment TOML is invalid.",
            ("Input", str(error.filename)),
            ("Reason", str(error.error)),
            ("Next", "Check the file against the TOML format."),
        )
    if isinstance(error, ExperimentNameError):
        return _fields(
            "Experiment name is missing or invalid.",
            ("Input", str(error.filename)),
            ("Reason", "The [experiment] table requires a name field."),
            ("Next", "Run 'chemex info' to list available Experiment Types."),
        )
    if isinstance(error, UnknownExperimentTypeError):
        return _fields(
            "Experiment Type is not available.",
            ("Input", str(error.filename)),
            ("Experiment Type", error.experiment_type_name),
            ("Next", "Run 'chemex info' to list available Experiment Types."),
        )
    if isinstance(error, ExperimentConfigurationError):
        return _fields(
            "Experiment configuration is invalid.",
            ("Input", str(error.filename)),
            ("Reason", _validation_reason(error.error)),
        )
    if isinstance(error, ExperimentDataError):
        if isinstance(error.error, DatasetLoadError):
            return _fields(
                "Experimental data is invalid.",
                ("Input", str(error.error.filename)),
                ("Experiment", str(error.filename)),
                ("Reason", error.error.explanation),
            )
        return _fields(
            "Experimental data could not be read.",
            ("Experiment", str(error.filename)),
            ("Reason", str(error.error)),
            ("Next", "Check the data path and read permissions."),
        )
    return None


def _format_scientific(error: ChemExError) -> str:
    return _fields(
        "The parameter state is scientifically invalid.",
        ("Reason", str(error)),
        *_analysis_fields(error),
        *_preserved_fields(error),
        ("Next", "Check the reported relaxation parameters and constraints."),
    )


def _format_parameter_expression(error: ChemExError) -> str | None:
    from chemex.parameters.database import (
        ConstraintExpressionError,
        GridExpressionError,
    )

    if isinstance(error, ConstraintExpressionError):
        return _fields(
            "Parameter constraint is invalid.",
            ("Constraint", error.expression),
            ("Reason", error.detail),
        )
    if isinstance(error, GridExpressionError):
        return _fields(
            "Parameter grid definition is invalid.",
            ("Grid", error.entry),
            ("Reason", error.detail),
        )
    return None


def _format_deterministic(error: ChemExError) -> str | None:
    from chemex.optimize.native_deterministic import NativeDeterministicAnalysisError

    if not isinstance(error, NativeDeterministicAnalysisError):
        return None
    details = "; ".join(
        failure.message for failure in error.failures if failure.message
    )
    return _fields(
        "Deterministic analysis did not produce an accepted fit.",
        ("Reason", error.reason),
        ("Details", details),
        *_analysis_fields(error),
        *_preserved_fields(error),
        ("Next", "Review the fit settings and starting parameter values."),
    )


def _format_mcmc(error: ChemExError) -> str | None:
    from chemex.optimize.mcmc import McmcConfigurationError, NativeMcmcIncompleteError

    if isinstance(error, McmcConfigurationError):
        return _fields(
            "MCMC configuration is invalid.",
            ("Reason", error.explanation),
            *_analysis_fields(error),
            *_preserved_fields(error),
            ("Next", "Adjust the MCMC settings or fitted parameter bounds."),
        )

    if not isinstance(error, NativeMcmcIncompleteError):
        return None
    diagnostics = getattr(error, "diagnostics_path", None)
    return _fields(
        "MCMC analysis is incomplete; posterior products were withheld.",
        ("Reason", str(error)),
        *_analysis_fields(error),
        *_preserved_fields(error),
        (
            "Next",
            "Inspect the MCMC diagnostics before changing the analysis."
            if isinstance(diagnostics, Path)
            else None,
        ),
    )


def _format_resampling(error: ChemExError) -> str | None:
    from chemex.optimize.resampling import NativeResamplingIncompleteError

    if not isinstance(error, NativeResamplingIncompleteError):
        return None
    diagnostics = getattr(error, "diagnostics_path", None)
    failures = getattr(error, "failures_path", None)
    if isinstance(diagnostics, Path) and isinstance(failures, Path):
        next_action = "Inspect the resampling diagnostics and failed replicates."
    elif isinstance(diagnostics, Path):
        next_action = "Inspect the resampling diagnostics."
    elif isinstance(failures, Path):
        next_action = "Inspect the failed resampling replicates."
    else:
        next_action = None
    return _fields(
        "Resampling analysis is incomplete; summary products were withheld.",
        ("Reason", str(error)),
        *_analysis_fields(error),
        *_preserved_fields(error),
        ("Next", next_action),
    )


_FORMATTERS: dict[str, _Formatter] = {
    "chemex.configuration.method_plan": _format_method,
    "chemex.configuration.parameters": _format_parameters,
    "chemex.exceptions": _format_path_error,
    "chemex.experiments.experiment_types": _format_experiment,
    "chemex.models.model": _format_model,
    "chemex.optimize.mcmc": _format_mcmc,
    "chemex.optimize.native_deterministic": _format_deterministic,
    "chemex.optimize.resampling": _format_resampling,
    "chemex.parameters.database": _format_parameter_expression,
    "chemex.parameters.feasible_coordinates": _format_scientific,
    "chemex.toml": _format_toml,
}


def format_known_failure(error: ChemExError) -> str:
    """Format one trusted ChemEx failure without importing Rich into its source."""
    formatter = _FORMATTERS.get(type(error).__module__)
    if formatter is not None and (message := formatter(error)):
        return message
    return str(error) or "ChemEx could not complete the requested operation."
