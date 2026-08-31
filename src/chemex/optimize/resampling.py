"""Output handling for resampling-based uncertainty analyses."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np

from chemex.atomic import write_text_atomic
from chemex.configuration.method_plan import ResamplingKind, ResamplingRequest
from chemex.containers.experiments import Experiments
from chemex.messages import print_running_statistics
from chemex.optimize.direct_trf import AcceptedFitResult
from chemex.optimize.native_deterministic import NativeDeterministicFit
from chemex.optimize.native_resampling import (
    OperationTerminal,
    ReplicateDisposition,
    ReplicateOutcome,
    ResamplingAnalysisResult,
    ResamplingAnalysisSample,
    ResamplingAnalysisStatus,
    ResamplingAnalysisSummary,
    ResamplingDatasetManifest,
    ResamplingOperation,
    ResamplingPlan,
    ResamplingScheme,
    derive_resampling_analysis_result,
    derive_resampling_diagnostic_samples,
    execute_resampling_evidence,
)
from chemex.optimize.statistics_plot import write_resampling_plots
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.runtime import ExecutionSettings


@dataclass(frozen=True, slots=True)
class _ResamplingMethod:
    statistic_name: str
    message: str
    directory: str
    iterations: int


class NativeResamplingIncompleteError(RuntimeError):
    """Raised after truthful incomplete native statistics have been published."""


def _resampling_result_and_samples(
    operation: ResamplingOperation,
    accepted: AcceptedFitResult,
    method_message: str,
) -> tuple[ResamplingAnalysisResult | None, Sequence[ResamplingAnalysisSample]]:
    if operation.terminal is OperationTerminal.INTERRUPTED:
        return None, derive_resampling_diagnostic_samples(operation, accepted)
    result = derive_resampling_analysis_result(operation, accepted)
    if result is None:
        raise NativeResamplingIncompleteError(
            f"Native {method_message} produced no meaningful Analysis Result"
        )
    return result, result.samples


def _resampling_method(
    kind: ResamplingKind,
    request: ResamplingRequest,
) -> _ResamplingMethod:
    descriptions = {
        "mc": ("Monte Carlo", "MonteCarlo"),
        "bs": ("bootstrap", "Bootstrap"),
        "bsn": ("nucleus-based bootstrap", "BootstrapNS"),
    }
    message, directory = descriptions[kind]
    return _ResamplingMethod(kind, message, directory, request.replicates)


def _quote_toml_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _format_toml_string_list(values: list[str]) -> str:
    if not values:
        return "[]"
    return "[" + ", ".join(_quote_toml_string(value) for value in values) + "]"


def _format_toml_float(value: float) -> str:
    return f"{value:.5e}"


def _format_tsv_float(value: float) -> str:
    return "nan" if not np.isfinite(value) else f"{value:.5e}"


def _format_parameter_names(
    parameter_ids: list[str],
    parameter_model: SealedParameterModel,
) -> tuple[str, ...]:
    return tuple(
        str(parameter_name_from_definition(parameter_model.definitions[param_id]))
        for param_id in parameter_ids
    )


def _write_resampling_summary(
    path: Path,
    *,
    parameter_names: tuple[str, ...],
    summary: ResamplingAnalysisSummary,
) -> None:
    lines: list[str] = []
    for parameter_name, distribution in zip(
        parameter_names,
        summary.distributions,
        strict=True,
    ):
        lines.extend(
            [
                f"[{_quote_toml_string(parameter_name.strip('[]'))}]",
                'interval = "95% percentile"',
                f"sample_count = {distribution.sample_count}",
            ],
        )
        lines.extend(
            [
                f"mean = {_format_toml_float(distribution.mean)}",
                (
                    "standard_deviation = "
                    f"{_format_toml_float(distribution.standard_deviation)}"
                ),
                f"median = {_format_toml_float(distribution.median)}",
                (
                    "percentile_95_lower = "
                    f"{_format_toml_float(distribution.percentile_95_lower)}"
                ),
                (
                    "percentile_95_upper = "
                    f"{_format_toml_float(distribution.percentile_95_upper)}"
                ),
                (
                    "percentile_68_lower = "
                    f"{_format_toml_float(distribution.percentile_68_lower)}"
                ),
                (
                    "percentile_68_upper = "
                    f"{_format_toml_float(distribution.percentile_68_upper)}"
                ),
                (
                    "half_percentile_68_width = "
                    f"{_format_toml_float(distribution.half_percentile_68_width)}"
                ),
            ],
        )
        lines.append("")
    (path / "summary.toml").write_text("\n".join(lines), encoding="utf-8")


def _write_resampling_correlations(
    path: Path,
    *,
    parameter_names: tuple[str, ...],
    summary: ResamplingAnalysisSummary,
) -> np.ndarray:
    width = len(summary.parameter_ids)
    correlations = np.asarray(
        tuple(
            np.nan if item.value is None else item.value
            for item in summary.correlations
        ),
        dtype=np.float64,
    ).reshape((width, width))
    if not parameter_names:
        (path / "correlations.tsv").write_text("", encoding="utf-8")
        return correlations

    with (path / "correlations.tsv").open("w", encoding="utf-8") as fileout:
        fileout.write("parameter\t" + "\t".join(parameter_names) + "\n")
        for name, values in zip(parameter_names, correlations, strict=True):
            row = "\t".join(_format_tsv_float(value) for value in values)
            fileout.write(f"{name}\t{row}\n")
    return correlations


def _write_resampling_diagnostics(
    path: Path,
    *,
    method: str,
    fitmethod: str,
    requested_samples: int,
    completed_samples: int,
    workers: int,
    parameter_ids: list[str],
    engine: str | None = None,
    root_seed: int | None = None,
    status: str | None = None,
) -> None:
    lines = [
        f"method = {_quote_toml_string(method)}",
        f"fitmethod = {_quote_toml_string(fitmethod)}",
        f"requested_samples = {requested_samples}",
        f"completed_samples = {completed_samples}",
        f"workers = {workers}",
        f"parameters = {_format_toml_string_list(parameter_ids)}",
        'samples_file = "samples.tsv"',
        'summary_file = "summary.toml"',
        'correlations_file = "correlations.tsv"',
        'plots_file = "plots.pdf"',
    ]
    if engine is not None:
        lines.append(f"engine = {_quote_toml_string(engine)}")
    if root_seed is not None:
        lines.append(f"root_seed = {root_seed}")
    if status is not None:
        lines.append(f"status = {_quote_toml_string(status)}")
    write_text_atomic(path / "diagnostics.toml", "\n".join(lines) + "\n")


def _native_dataset(
    experiments: Experiments,
    fit: NativeDeterministicFit,
) -> ResamplingDatasetManifest:
    profiles = tuple(profile for experiment in experiments for profile in experiment)
    if len(profiles) != len(fit.engine.plan.profiles):
        raise NativeResamplingIncompleteError(
            "Native resampling profile population differs from the accepted fit"
        )
    references: list[bool] = []
    nucleus_groups: list[str] = []
    descriptors: list[str] = []
    for profile, profile_plan in zip(profiles, fit.engine.plan.profiles, strict=True):
        if profile.data.size != profile_plan.observation_count:
            raise NativeResamplingIncompleteError(
                "Native resampling profile shape differs from the accepted fit"
            )
        references.extend(bool(value) for value in profile.data.refs)
        nucleus_group = str(profile.spin_system.groups["i"])
        nucleus_groups.extend([nucleus_group] * profile.data.size)
        descriptors.extend(
            f"{profile_plan.identity}:{index}"
            for index in range(profile_plan.observation_count)
        )
    return ResamplingDatasetManifest(
        fit.engine.plan,
        tuple(
            float(value)
            for value in fit.accepted.evaluation_result.normalized_calculations
        ),
        tuple(references),
        tuple(nucleus_groups),
        tuple(descriptors),
    )


def _write_native_failures(path: Path, outcomes: Sequence[ReplicateOutcome]) -> None:
    failed = tuple(outcome for outcome in outcomes if outcome.failure is not None)
    if not failed:
        return
    with (path / "failures.tsv").open("w", encoding="utf-8") as output:
        output.write("ordinal\tdisposition\tcategory\tmessage\n")
        for outcome in failed:
            failure = outcome.failure
            if failure is None:
                continue
            message = failure.message.replace("\t", " ").replace("\n", " ")
            output.write(
                f"{outcome.ordinal}\t{outcome.disposition.value}\t"
                f"{failure.category}\t{message}\n"
            )


def _write_native_state_diagnostics(
    path: Path,
    *,
    method: _ResamplingMethod,
    workers: int,
    parameter_ids: tuple[str, ...],
    root_seed: int,
    status: str,
    successful_samples: int | None = None,
    failed_samples: int | None = None,
    cancelled_samples: int | None = None,
    interrupted_samples: int | None = None,
    unstarted_samples: int | None = None,
    terminal: str | None = None,
    failure: BaseException | None = None,
) -> None:
    lines = [
        f"method = {_quote_toml_string(method.message)}",
        'fitmethod = "trf"',
        'engine = "native direct TRF"',
        f"status = {_quote_toml_string(status)}",
        f"requested_samples = {method.iterations}",
        f"workers = {workers}",
        f"parameters = {_format_toml_string_list(list(parameter_ids))}",
        f"root_seed = {root_seed}",
    ]
    if terminal is not None:
        lines.append(f"terminal = {_quote_toml_string(terminal)}")
    sample_counts = (
        ("completed_samples", successful_samples),
        ("failed_samples", failed_samples),
        ("cancelled_samples", cancelled_samples),
        ("interrupted_samples", interrupted_samples),
        ("unstarted_samples", unstarted_samples),
    )
    lines.extend(
        f"{name} = {count}" for name, count in sample_counts if count is not None
    )
    if status == "incomplete":
        if (path / "samples.tsv").is_file():
            lines.append('samples_file = "samples.tsv"')
        if (path / "failures.tsv").is_file():
            lines.append('failures_file = "failures.tsv"')
    if failure is not None:
        message = str(failure).replace("\n", " ")
        lines.extend(
            (
                f"failure_type = {_quote_toml_string(type(failure).__name__)}",
                f"failure_message = {_quote_toml_string(message)}",
            )
        )
    write_text_atomic(path / "diagnostics.toml", "\n".join(lines) + "\n")


def _clear_native_statistics_artifacts(path: Path) -> None:
    for name in (
        "summary.toml",
        "samples.tsv",
        "correlations.tsv",
        "diagnostics.toml",
        "plots.pdf",
        "failures.tsv",
    ):
        (path / name).unlink(missing_ok=True)


def _remove_failed_publication_artifacts(
    path: Path,
    names: Sequence[str],
    failure: BaseException,
) -> None:
    for name in names:
        try:
            (path / name).unlink(missing_ok=True)
        except (Exception, KeyboardInterrupt) as cleanup_error:  # noqa: BLE001
            failure.add_note(f"ChemEx could not remove {name}: {cleanup_error}")


def _run_native_resampling_method(
    experiments: Experiments,
    path: Path,
    method: _ResamplingMethod,
    fit: NativeDeterministicFit,
    *,
    execution: ExecutionSettings,
    root_seed: int,
) -> None:
    statistic_path = path / "Statistics" / method.directory
    statistic_path.mkdir(parents=True, exist_ok=True)
    _clear_native_statistics_artifacts(statistic_path)
    parameter_ids = fit.accepted.controlled_ids
    parameter_names = _format_parameter_names(list(parameter_ids), fit.parameter_model)
    worker_count = min(execution.workers, method.iterations)
    _write_native_state_diagnostics(
        statistic_path,
        method=method,
        workers=worker_count,
        parameter_ids=parameter_ids,
        root_seed=root_seed,
        status="running",
    )
    try:
        dataset = _native_dataset(experiments, fit)
        width = max(6, len(str(method.iterations)))
        plan = ResamplingPlan.for_accepted(
            fit.accepted,
            dataset=dataset,
            source_problem=fit.problem,
            parameterization=fit.parameterization,
            source_engine=fit.engine,
            scheme=ResamplingScheme(method.statistic_name),
            replicate_count=method.iterations,
            replicate_structural_identities=tuple(
                f"production-replicate-{index:0{width}d}"
                for index in range(method.iterations)
            ),
            replicate_component_identities=tuple(
                (f"production-direct-trf-{index:0{width}d}",)
                for index in range(method.iterations)
            ),
            root_seed=root_seed,
            output_scope=fit.problem.commit_scope,
            output_units=("native",) * len(fit.problem.commit_scope),
            minimum_successful_count=method.iterations,
            strategy_settings=(
                (
                    "objective_request_budget",
                    str(fit.objective_request_budget),
                ),
            ),
        )
        operation = execute_resampling_evidence(
            fit.accepted,
            plan,
            execution=execution,
        )
    except (Exception, KeyboardInterrupt) as error:
        terminal = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        _write_native_state_diagnostics(
            statistic_path,
            method=method,
            workers=worker_count,
            parameter_ids=parameter_ids,
            root_seed=root_seed,
            status="incomplete",
            successful_samples=0,
            failed_samples=0,
            cancelled_samples=0,
            interrupted_samples=0,
            unstarted_samples=method.iterations,
            terminal=terminal,
            failure=error,
        )
        raise
    evidence = operation.evidence
    if evidence is None:
        error = NativeResamplingIncompleteError(
            f"Native {method.message} produced no replicate evidence"
        )
        _write_native_state_diagnostics(
            statistic_path,
            method=method,
            workers=worker_count,
            parameter_ids=parameter_ids,
            root_seed=root_seed,
            status="incomplete",
            successful_samples=0,
            failed_samples=0,
            cancelled_samples=0,
            interrupted_samples=0,
            unstarted_samples=method.iterations,
            terminal=operation.terminal.value,
            failure=error,
        )
        raise error
    samples_published = False
    failures_published = False
    result: ResamplingAnalysisResult | None = None
    operational_interruption = operation.terminal is OperationTerminal.INTERRUPTED
    try:
        result, analysis_samples = _resampling_result_and_samples(
            operation,
            fit.accepted,
            method.message,
        )
        complete = (
            result is not None and result.status is ResamplingAnalysisStatus.COMPLETE
        )
        sample_rows = [sample.values for sample in analysis_samples]
        chisqr_rows = [sample.chi_square for sample in analysis_samples]
        samples = np.asarray(sample_rows, dtype=float).reshape(
            (len(sample_rows), len(parameter_ids))
        )
        with (statistic_path / "samples.tsv").open("w", encoding="utf-8") as output:
            output.write("\t".join((*parameter_names, "chisqr")) + "\n")
            for values, chisqr in zip(sample_rows, chisqr_rows, strict=True):
                output.write(
                    "\t".join(
                        _format_tsv_float(value) for value in (*values, float(chisqr))
                    )
                    + "\n"
                )
        samples_published = True
        if not complete:
            _write_native_failures(statistic_path, evidence.outcomes)
            failures_published = True
            disposition_source = evidence if result is None else result
            _write_native_state_diagnostics(
                statistic_path,
                method=method,
                workers=worker_count,
                parameter_ids=parameter_ids,
                root_seed=root_seed,
                status="incomplete",
                successful_samples=disposition_source.disposition_count(
                    ReplicateDisposition.SUCCEEDED
                ),
                failed_samples=disposition_source.disposition_count(
                    ReplicateDisposition.FAILED
                ),
                cancelled_samples=disposition_source.disposition_count(
                    ReplicateDisposition.CANCELLED
                ),
                interrupted_samples=disposition_source.disposition_count(
                    ReplicateDisposition.INTERRUPTED
                ),
                unstarted_samples=disposition_source.disposition_count(
                    ReplicateDisposition.NOT_STARTED
                ),
                terminal=operation.terminal.value,
            )
        else:
            summary = cast("ResamplingAnalysisSummary", result.summary)
            _write_resampling_summary(
                statistic_path,
                parameter_names=parameter_names,
                summary=summary,
            )
            correlations = _write_resampling_correlations(
                statistic_path,
                parameter_names=parameter_names,
                summary=summary,
            )
            write_resampling_plots(
                statistic_path,
                method=method.message,
                fitmethod="trf",
                parameter_names=parameter_names,
                samples=samples,
                chisqr_values=np.asarray(chisqr_rows, dtype=float),
                correlations=correlations,
                summary=summary,
                requested_samples=method.iterations,
                completed_samples=result.disposition_count(
                    ReplicateDisposition.SUCCEEDED
                ),
            )
            _write_resampling_diagnostics(
                statistic_path,
                method=method.message,
                fitmethod="trf",
                requested_samples=method.iterations,
                completed_samples=result.disposition_count(
                    ReplicateDisposition.SUCCEEDED
                ),
                workers=worker_count,
                parameter_ids=list(parameter_ids),
                engine="native direct TRF",
                root_seed=root_seed,
                status="complete",
            )
    except (Exception, KeyboardInterrupt) as error:
        artifacts_to_remove = ["summary.toml", "correlations.tsv", "plots.pdf"]
        if not samples_published:
            artifacts_to_remove.append("samples.tsv")
        if not failures_published:
            artifacts_to_remove.append("failures.tsv")
        _remove_failed_publication_artifacts(
            statistic_path,
            artifacts_to_remove,
            error,
        )
        terminal = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        disposition_source = evidence if result is None else result
        try:
            _write_native_state_diagnostics(
                statistic_path,
                method=method,
                workers=worker_count,
                parameter_ids=parameter_ids,
                root_seed=root_seed,
                status="incomplete",
                successful_samples=disposition_source.disposition_count(
                    ReplicateDisposition.SUCCEEDED
                ),
                failed_samples=disposition_source.disposition_count(
                    ReplicateDisposition.FAILED
                ),
                cancelled_samples=disposition_source.disposition_count(
                    ReplicateDisposition.CANCELLED
                ),
                interrupted_samples=disposition_source.disposition_count(
                    ReplicateDisposition.INTERRUPTED
                ),
                unstarted_samples=disposition_source.disposition_count(
                    ReplicateDisposition.NOT_STARTED
                ),
                terminal=terminal,
                failure=error,
            )
        except (Exception, KeyboardInterrupt) as diagnostics_error:  # noqa: BLE001
            error.add_note(
                "ChemEx could not publish incomplete resampling diagnostics: "
                f"{diagnostics_error}"
            )
            _remove_failed_publication_artifacts(
                statistic_path,
                ("diagnostics.toml",),
                error,
            )
        raise
    if operational_interruption or (
        result is not None and result.status is ResamplingAnalysisStatus.INCOMPLETE
    ):
        raise NativeResamplingIncompleteError(
            f"Native {method.message} completed {evidence.successful_count} of "
            f"{method.iterations} requested replicates"
        )


def run_native_resampling(
    experiments: Experiments,
    path: Path,
    kind: ResamplingKind,
    request: ResamplingRequest,
    fit: NativeDeterministicFit,
    *,
    execution: ExecutionSettings | None = None,
    root_seed: int = 0,
) -> None:
    """Run one canonical resampling request from a committed native fit."""
    settings = ExecutionSettings() if execution is None else execution
    method = _resampling_method(kind, request)
    print_running_statistics(method.message)
    _run_native_resampling_method(
        experiments,
        path,
        method,
        fit,
        execution=settings,
        root_seed=root_seed,
    )
