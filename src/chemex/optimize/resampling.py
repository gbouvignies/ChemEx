"""Output handling for resampling-based uncertainty analyses."""

from __future__ import annotations

import multiprocessing
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rich.progress import track

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
    (path / "diagnostics.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    methods: list[_ResamplingMethod] = []
    if statistics.mc is not None:
        methods.append(
            _ResamplingMethod(
                statistic_name="mc",
                message="Monte Carlo",
                directory="MonteCarlo",
                iterations=statistics.mc,
            ),
        )
    if statistics.bs is not None:
        methods.append(
            _ResamplingMethod(
                statistic_name="bs",
                message="bootstrap",
                directory="Bootstrap",
                iterations=statistics.bs,
            ),
        )
    if statistics.bsn is not None:
        methods.append(
            _ResamplingMethod(
                statistic_name="bsn",
                message="nucleus-based bootstrap",
                directory="BootstrapNS",
                iterations=statistics.bsn,
            ),
        )

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
