"""Output handling for resampling-based uncertainty analyses."""

from __future__ import annotations

import multiprocessing
from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rich.progress import track

from chemex.atomic import write_text_atomic
from chemex.configuration.methods import Statistics
from chemex.containers.experiments import Experiments, generate_exp_for_statistics
from chemex.messages import (
    print_calculation_stopped_error,
    print_running_statistics,
    print_value_error,
)
from chemex.optimize.helper import (
    calculate_statistics,
)
from chemex.optimize.minimizer import minimize
from chemex.optimize.native_deterministic import NativeDeterministicFit
from chemex.optimize.native_resampling import (
    OperationTerminal,
    ReplicateDisposition,
    ReplicateOutcome,
    ResamplingDatasetManifest,
    ResamplingPlan,
    ResamplingScheme,
    execute_resampling_evidence,
)
from chemex.optimize.statistics_plot import write_resampling_plots
from chemex.runtime import ExecutionSettings
from chemex.runtime.execution import native_thread_environment


@dataclass(frozen=True, slots=True)
class _ResamplingMethod:
    statistic_name: str
    message: str
    directory: str
    iterations: int


@dataclass(frozen=True, slots=True)
class _ResamplingContext:
    experiments: Experiments
    statistic_name: str
    fitmethod: str
    parameter_ids: tuple[str, ...]


@dataclass(slots=True)
class _WorkerState:
    context: _ResamplingContext | None = None

    def initialize(self, context: _ResamplingContext) -> None:
        self.context = context

    def clear(self) -> None:
        self.context = None


_WORKER_STATE = _WorkerState()


class NativeResamplingIncompleteError(RuntimeError):
    """Raised after truthful incomplete native statistics have been published."""


def _resampling_methods(statistics: Statistics) -> tuple[_ResamplingMethod, ...]:
    methods: list[_ResamplingMethod] = []
    if statistics.mc is not None:
        methods.append(
            _ResamplingMethod("mc", "Monte Carlo", "MonteCarlo", statistics.mc)
        )
    if statistics.bs is not None:
        methods.append(_ResamplingMethod("bs", "bootstrap", "Bootstrap", statistics.bs))
    if statistics.bsn is not None:
        methods.append(
            _ResamplingMethod(
                "bsn",
                "nucleus-based bootstrap",
                "BootstrapNS",
                statistics.bsn,
            )
        )
    return tuple(methods)


def _initialize_worker(context: _ResamplingContext) -> None:
    _WORKER_STATE.initialize(context)


def _clear_worker() -> None:
    _WORKER_STATE.clear()


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
    parameter_store: Any,
) -> tuple[str, ...]:
    if not hasattr(parameter_store, "get_parameters"):
        return tuple(parameter_ids)
    parameters = parameter_store.get_parameters(parameter_ids)
    return tuple(
        str(parameters[param_id].param_name) if param_id in parameters else param_id
        for param_id in parameter_ids
    )


def _sample_parameter_values(
    params_lf: Any,
    parameter_ids: tuple[str, ...],
) -> list[float]:
    return [
        float(params_lf[param_id].value) if param_id in params_lf else np.nan
        for param_id in parameter_ids
    ]


def _calculate_resampling_sample(
    context: _ResamplingContext,
) -> tuple[list[float], float]:
    parameter_store = context.experiments.parameter_store
    exp_stat = generate_exp_for_statistics(context.experiments, context.statistic_name)
    params_lf = parameter_store.build_lmfit_params(exp_stat.param_ids)
    params_fit = minimize(exp_stat, params_lf, context.fitmethod)
    stats = calculate_statistics(exp_stat, params_fit)
    chisqr = stats.get("chisqr", 1e32)
    return _sample_parameter_values(params_fit, context.parameter_ids), float(chisqr)


def _calculate_worker_sample(_index: int) -> tuple[list[float], float]:
    context = _WORKER_STATE.context
    if context is None:
        msg = "Resampling worker context is not initialized"
        raise RuntimeError(msg)
    return _calculate_resampling_sample(context)


def _iter_resampling_samples(
    context: _ResamplingContext,
    iterations: int,
    *,
    execution: ExecutionSettings,
) -> Iterator[tuple[list[float], float]]:
    if execution.workers == 1 or iterations == 1:
        for _index in range(iterations):
            yield _calculate_resampling_sample(context)
        return

    worker_count = min(execution.workers, iterations)
    try:
        with (
            native_thread_environment(execution.native_thread_env(parallel=True)),
            multiprocessing.Pool(
                processes=worker_count,
                initializer=_initialize_worker,
                initargs=(context,),
            ) as pool,
        ):
            yield from pool.imap(_calculate_worker_sample, range(iterations))
    finally:
        _clear_worker()


def _as_sample_array(sample_rows: list[list[float]], width: int) -> np.ndarray:
    if not sample_rows:
        return np.empty((0, width), dtype=float)
    return np.asarray(sample_rows, dtype=float)


def _write_resampling_summary(
    path: Path,
    *,
    parameter_names: tuple[str, ...],
    samples: np.ndarray,
) -> None:
    lines: list[str] = []
    for index, parameter_name in enumerate(parameter_names):
        values = samples[:, index] if samples.size else np.empty(0, dtype=float)
        values = values[np.isfinite(values)]
        lines.extend(
            [
                f"[{_quote_toml_string(parameter_name.strip('[]'))}]",
                'interval = "95% percentile"',
                f"sample_count = {len(values)}",
            ],
        )
        if len(values) > 0:
            lower_95, lower, median, upper, upper_95 = np.percentile(
                values,
                [2.5, 15.87, 50.0, 84.13, 97.5],
            )
            standard_deviation = (
                float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            )
            lines.extend(
                [
                    f"mean = {_format_toml_float(float(np.mean(values)))}",
                    f"standard_deviation = {_format_toml_float(standard_deviation)}",
                    f"median = {_format_toml_float(float(median))}",
                    f"percentile_95_lower = {_format_toml_float(float(lower_95))}",
                    f"percentile_95_upper = {_format_toml_float(float(upper_95))}",
                    f"lower_1sigma = {_format_toml_float(float(lower))}",
                    f"upper_1sigma = {_format_toml_float(float(upper))}",
                    f"stderr = {_format_toml_float(0.5 * float(upper - lower))}",
                ],
            )
        lines.append("")
    (path / "summary.toml").write_text("\n".join(lines), encoding="utf-8")


def _pairwise_correlation(samples: np.ndarray, index_a: int, index_b: int) -> float:
    values_a = samples[:, index_a]
    values_b = samples[:, index_b]
    valid = np.isfinite(values_a) & np.isfinite(values_b)
    if int(np.sum(valid)) < 2:
        return np.nan
    return float(np.corrcoef(values_a[valid], values_b[valid])[0, 1])


def _write_resampling_correlations(
    path: Path,
    *,
    parameter_names: tuple[str, ...],
    samples: np.ndarray,
) -> np.ndarray:
    if not parameter_names:
        (path / "correlations.tsv").write_text("", encoding="utf-8")
        return np.empty((0, 0), dtype=float)
    correlations = np.empty((len(parameter_names), len(parameter_names)), dtype=float)
    for index_a in range(len(parameter_names)):
        for index_b in range(len(parameter_names)):
            correlations[index_a, index_b] = (
                1.0
                if index_a == index_b
                else _pairwise_correlation(samples, index_a, index_b)
            )

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


def _run_resampling_method(
    experiments: Experiments,
    path: Path,
    fitmethod: str,
    method: _ResamplingMethod,
    parameter_ids: tuple[str, ...],
    *,
    execution: ExecutionSettings,
) -> None:
    parameter_store = experiments.parameter_store
    statistic_path = path / "Statistics" / method.directory
    statistic_path.mkdir(parents=True, exist_ok=True)
    samples_tsv = statistic_path / "samples.tsv"
    sample_rows: list[list[float]] = []
    chisqr_rows: list[float] = []
    parameter_names = _format_parameter_names(list(parameter_ids), parameter_store)
    completed_samples = 0
    context = _ResamplingContext(
        experiments=experiments,
        statistic_name=method.statistic_name,
        fitmethod=fitmethod,
        parameter_ids=parameter_ids,
    )
    worker_count = min(execution.workers, method.iterations)

    with samples_tsv.open(mode="w", encoding="utf-8") as file_tsv:
        file_tsv.write("\t".join((*parameter_names, "chisqr")) + "\n")

        try:
            samples = _iter_resampling_samples(
                context,
                method.iterations,
                execution=execution,
            )
            for sample_values, chisqr in track(
                samples,
                total=method.iterations,
                description="   ",
            ):
                file_tsv.write(
                    "\t".join(
                        _format_tsv_float(value)
                        for value in (*sample_values, float(chisqr))
                    )
                    + "\n",
                )
                sample_rows.append(sample_values)
                chisqr_rows.append(float(chisqr))
                completed_samples += 1
        except KeyboardInterrupt:
            print_calculation_stopped_error()
        except ValueError:
            print_value_error()
        finally:
            file_tsv.flush()
            samples = _as_sample_array(sample_rows, len(parameter_ids))
            _write_resampling_summary(
                statistic_path,
                parameter_names=parameter_names,
                samples=samples,
            )
            correlations = _write_resampling_correlations(
                statistic_path,
                parameter_names=parameter_names,
                samples=samples,
            )
            write_resampling_plots(
                statistic_path,
                method=method.message,
                fitmethod=fitmethod,
                parameter_names=parameter_names,
                samples=samples,
                chisqr_values=np.asarray(chisqr_rows, dtype=float),
                correlations=correlations,
                requested_samples=method.iterations,
                completed_samples=completed_samples,
            )
            _write_resampling_diagnostics(
                statistic_path,
                method=method.message,
                fitmethod=fitmethod,
                requested_samples=method.iterations,
                completed_samples=completed_samples,
                workers=worker_count,
                parameter_ids=list(parameter_ids),
            )


def run_resampling_statistics(
    experiments: Experiments,
    path: Path,
    fitmethod: str,
    statistics: Statistics,
    *,
    execution: ExecutionSettings | None = None,
) -> None:
    if statistics.mc is None and statistics.bs is None and statistics.bsn is None:
        return

    execution = ExecutionSettings() if execution is None else execution
    methods = _resampling_methods(statistics)

    parameter_store = experiments.parameter_store
    params_lf = parameter_store.build_lmfit_params(experiments.param_ids)
    parameter_ids = tuple(param.name for param in params_lf.values() if param.vary)

    for method in methods:
        print_running_statistics(method.message)
        _run_resampling_method(
            experiments,
            path,
            fitmethod,
            method,
            parameter_ids,
            execution=execution,
        )


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
    parameter_names = _format_parameter_names(
        list(parameter_ids), experiments.parameter_store
    )
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
    disposition_counts = Counter(outcome.disposition for outcome in evidence.outcomes)
    samples_published = False
    failures_published = False
    try:
        complete = (
            operation.terminal is OperationTerminal.COMPLETED
            and evidence.successful_count == method.iterations
        )
        successful = tuple(
            outcome.success
            for outcome in evidence.outcomes
            if outcome.success is not None
        )
        sample_rows = [
            [dict(success.resolved_items)[param_id] for param_id in parameter_ids]
            for success in successful
        ]
        chisqr_rows = [success.chi_square for success in successful]
        samples = _as_sample_array(sample_rows, len(parameter_ids))
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
            _write_native_state_diagnostics(
                statistic_path,
                method=method,
                workers=worker_count,
                parameter_ids=parameter_ids,
                root_seed=root_seed,
                status="incomplete",
                successful_samples=disposition_counts[ReplicateDisposition.SUCCEEDED],
                failed_samples=disposition_counts[ReplicateDisposition.FAILED],
                cancelled_samples=disposition_counts[ReplicateDisposition.CANCELLED],
                interrupted_samples=disposition_counts[
                    ReplicateDisposition.INTERRUPTED
                ],
                unstarted_samples=disposition_counts[ReplicateDisposition.NOT_STARTED],
                terminal=operation.terminal.value,
            )
        else:
            _write_resampling_summary(
                statistic_path,
                parameter_names=parameter_names,
                samples=samples,
            )
            correlations = _write_resampling_correlations(
                statistic_path,
                parameter_names=parameter_names,
                samples=samples,
            )
            write_resampling_plots(
                statistic_path,
                method=method.message,
                fitmethod="trf",
                parameter_names=parameter_names,
                samples=samples,
                chisqr_values=np.asarray(chisqr_rows, dtype=float),
                correlations=correlations,
                requested_samples=method.iterations,
                completed_samples=evidence.successful_count,
            )
            _write_resampling_diagnostics(
                statistic_path,
                method=method.message,
                fitmethod="trf",
                requested_samples=method.iterations,
                completed_samples=evidence.successful_count,
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
        try:
            _write_native_state_diagnostics(
                statistic_path,
                method=method,
                workers=worker_count,
                parameter_ids=parameter_ids,
                root_seed=root_seed,
                status="incomplete",
                successful_samples=disposition_counts[ReplicateDisposition.SUCCEEDED],
                failed_samples=disposition_counts[ReplicateDisposition.FAILED],
                cancelled_samples=disposition_counts[ReplicateDisposition.CANCELLED],
                interrupted_samples=disposition_counts[
                    ReplicateDisposition.INTERRUPTED
                ],
                unstarted_samples=disposition_counts[ReplicateDisposition.NOT_STARTED],
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
    if not complete:
        raise NativeResamplingIncompleteError(
            f"Native {method.message} completed {evidence.successful_count} of "
            f"{method.iterations} requested replicates"
        )


def run_native_resampling_statistics(
    experiments: Experiments,
    path: Path,
    statistics: Statistics,
    fit: NativeDeterministicFit,
    *,
    execution: ExecutionSettings | None = None,
    root_seed: int = 0,
) -> None:
    """Run existing MC/BS/BSN product analyses from one committed native fit."""
    settings = ExecutionSettings() if execution is None else execution
    for method in _resampling_methods(statistics):
        print_running_statistics(method.message)
        _run_native_resampling_method(
            experiments,
            path,
            method,
            fit,
            execution=settings,
            root_seed=root_seed,
        )
