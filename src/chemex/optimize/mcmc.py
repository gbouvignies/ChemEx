from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from math import ceil
from pathlib import Path
from time import perf_counter
from types import TracebackType
from typing import Any, Self, cast

import emcee.autocorr as emcee_autocorr
import emcee.ensemble as emcee_ensemble
import numpy as np
from lmfit.parameter import Parameters
from rich.progress import Progress, TaskID

from chemex.configuration.methods import McmcBurnSetting, McmcSettings
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_mcmc_no_vary_warning,
    print_mcmc_unbounded_warning,
)
from chemex.optimize.mcmc_engine import (
    EmceeSamplerResult,
    McmcProblem,
    run_emcee_sampler,
)
from chemex.optimize.statistics_plot import write_mcmc_plots
from chemex.parameters.database import ParameterStore
from chemex.runtime import ExecutionSettings
from chemex.runtime.execution import native_thread_env, native_thread_environment
from chemex.typing import Array


class _RichEmceeProgressBar:
    def __init__(self, total: int) -> None:
        self._total = total
        self._progress = Progress()
        self._task_id: TaskID | None = None

    def __enter__(self) -> Self:
        self._progress.start()
        self._task_id = self._progress.add_task("   ", total=self._total)
        return self

    def __exit__(
        self,
        _exc_type: type[BaseException] | None,
        _exc: BaseException | None,
        _traceback: TracebackType | None,
    ) -> None:
        self._progress.stop()

    def update(self, count: int) -> None:
        if self._task_id is not None:
            self._progress.update(self._task_id, advance=count)


@contextmanager
def _use_rich_emcee_progress() -> Iterator[None]:
    original_get_progress_bar = emcee_ensemble.get_progress_bar

    def get_progress_bar(
        display: bool | str,
        total: int,
        **_kwargs: object,
    ) -> object:
        if display:
            return _RichEmceeProgressBar(total)
        return original_get_progress_bar(display, total, **_kwargs)

    emcee_ensemble_dynamic = cast("Any", emcee_ensemble)
    emcee_ensemble_dynamic.get_progress_bar = get_progress_bar
    try:
        yield
    finally:
        emcee_ensemble_dynamic.get_progress_bar = original_get_progress_bar


@dataclass(frozen=True)
class EffectiveMcmcSettings:
    steps: int
    burn: McmcBurnSetting
    thin: int
    walkers: int
    seed: int | None
    workers: int
    native_threads: int | None
    update_parameters: bool


@dataclass(frozen=True)
class McmcSummary:
    parameter_id: str
    mean: float
    standard_deviation: float
    median: float
    eti_95_lower: float
    eti_95_upper: float
    lower_1sigma: float
    upper_1sigma: float
    stderr: float
    effective_sample_size: float | None = None
    mcse_mean: float | None = None


@dataclass(frozen=True)
class McmcResult:
    var_names: tuple[str, ...]
    chain: Array
    lnprob: Array
    summary: tuple[McmcSummary, ...]
    correlations: Array
    acceptance_fraction: Array
    autocorrelation_time: Array | None
    discarded_steps: int
    burn_in_warning: str | None
    tentative_autocorrelation_time: Array | None = None
    autocorrelation_warning: str | None = None
    raw_chain: Array | None = None
    raw_lnprob: Array | None = None

    @property
    def samples(self) -> Array:
        """Return the flattened retained chain."""
        return self.chain.reshape((-1, len(self.var_names)))

    @property
    def log_probabilities(self) -> Array:
        """Return the flattened retained log probabilities."""
        return self.lnprob.reshape(-1)


_TIMING_KEYS = (
    "sampling_seconds",
    "result_processing_seconds",
    "output_summary_seconds",
    "output_samples_seconds",
    "output_correlations_seconds",
    "output_plots_seconds",
    "output_total_seconds",
    "total_seconds",
)


def resolve_mcmc_settings(
    settings: McmcSettings,
    *,
    nvarys: int,
    execution: ExecutionSettings | None = None,
) -> EffectiveMcmcSettings:
    walkers = settings.walkers or max(32, 2 * nvarys)
    min_walkers = 2 * nvarys
    if walkers < min_walkers:
        msg = (
            f"MCMC requires at least {min_walkers} walkers for {nvarys} fitted "
            f"parameters"
        )
        raise ValueError(msg)

    workers = settings.workers or (execution.workers if execution is not None else 1)
    native_threads = execution.native_threads if execution is not None else None

    return EffectiveMcmcSettings(
        steps=settings.steps,
        burn=settings.burn,
        thin=settings.thin,
        walkers=walkers,
        seed=settings.seed,
        workers=workers,
        native_threads=native_threads,
        update_parameters=settings.update_parameters,
    )


def _varying_parameter_ids(params: Parameters) -> tuple[str, ...]:
    return tuple(
        name
        for name, parameter in params.items()
        if parameter.vary and not parameter.expr
    )


def _format_parameter_ids(
    parameter_ids: tuple[str, ...],
    parameter_store: ParameterStore,
) -> tuple[str, ...]:
    parameters = parameter_store.get_parameters(parameter_ids)
    return tuple(str(parameters[param_id].param_name) for param_id in parameter_ids)


def _find_unbounded_parameter_ids(
    params: Parameters,
    var_names: tuple[str, ...],
) -> tuple[str, ...]:
    return tuple(
        name
        for name in var_names
        if not np.isfinite(params[name].min) or not np.isfinite(params[name].max)
    )


def _validate_parameter_bounds(params: Parameters, var_names: tuple[str, ...]) -> None:
    for name in var_names:
        parameter = params[name]
        if (
            np.isfinite(parameter.min)
            and np.isfinite(parameter.max)
            and parameter.min >= parameter.max
        ):
            msg = f"MCMC parameter {name!r} has an empty bounded interval"
            raise ValueError(msg)


def _parameter_bound(value: float | None, fallback: float) -> float:
    if value is None:
        return fallback
    bound = float(value)
    if np.isnan(bound):
        return fallback
    return bound


def _parameter_bounds(params: Parameters, var_names: tuple[str, ...]) -> Array:
    return np.asarray(
        [
            (
                _parameter_bound(params[name].min, -np.inf),
                _parameter_bound(params[name].max, np.inf),
            )
            for name in var_names
        ],
        dtype=float,
    )


def _jitter_scale(
    value: float,
    lower: float,
    upper: float,
    stderr: float | None,
) -> float:
    scales = [abs(value) * 1.0e-4, 1.0e-8]
    if np.isfinite(lower) and np.isfinite(upper):
        scales.append((upper - lower) * 1.0e-4)
    if stderr is not None and np.isfinite(stderr) and stderr > 0.0:
        scales.append(stderr * 1.0e-2)
    return max(scales)


def _clip_to_bounds(values: Array, lower: float, upper: float) -> Array:
    if np.isfinite(lower) and np.isfinite(upper):
        width = upper - lower
        eps = max(width * 1.0e-10, np.finfo(float).eps)
        return np.clip(values, lower + eps, upper - eps)
    if np.isfinite(lower):
        return np.maximum(values, lower + np.finfo(float).eps)
    if np.isfinite(upper):
        return np.minimum(values, upper - np.finfo(float).eps)
    return values


def _initial_positions(
    params: Parameters,
    var_names: tuple[str, ...],
    *,
    nwalkers: int,
    seed: int | None,
) -> Array:
    rng = np.random.default_rng(seed)
    positions = np.empty((nwalkers, len(var_names)), dtype=float)

    for index, name in enumerate(var_names):
        parameter = params[name]
        value = float(parameter.value)
        lower = float(parameter.min)
        upper = float(parameter.max)
        scale = _jitter_scale(value, lower, upper, parameter.stderr)
        values = value + scale * rng.standard_normal(nwalkers)
        positions[:, index] = _clip_to_bounds(values, lower, upper)

    return positions


def _summarize_chain(
    var_names: tuple[str, ...],
    samples: Array,
    autocorrelation_time: Array | None,
    thin: int,
) -> tuple[McmcSummary, ...]:
    quantiles = np.percentile(samples, [2.5, 15.87, 50.0, 84.13, 97.5], axis=0)
    summaries: list[McmcSummary] = []
    for index, name in enumerate(var_names):
        lower_95, lower, median, upper, upper_95 = quantiles[:, index]
        values = samples[:, index]
        standard_deviation = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
        effective_sample_size = None
        mcse_mean = None
        if autocorrelation_time is not None:
            tau = float(autocorrelation_time[index])
            if np.isfinite(tau) and tau > 0.0:
                tau_in_retained_steps = max(tau / thin, 1.0)
                effective_sample_size = len(values) / tau_in_retained_steps
                mcse_mean = standard_deviation / np.sqrt(effective_sample_size)
        summaries.append(
            McmcSummary(
                parameter_id=name,
                mean=float(np.mean(values)),
                standard_deviation=standard_deviation,
                median=float(median),
                eti_95_lower=float(lower_95),
                eti_95_upper=float(upper_95),
                lower_1sigma=float(lower),
                upper_1sigma=float(upper),
                stderr=0.5 * float(upper - lower),
                effective_sample_size=effective_sample_size,
                mcse_mean=mcse_mean,
            ),
        )
    return tuple(summaries)


def _correlation_matrix(samples: Array) -> Array:
    if samples.shape[1] == 1:
        return np.ones((1, 1), dtype=float)
    return np.corrcoef(samples, rowvar=False)


def _valid_autocorrelation_time(autocorrelation_time: object) -> Array | None:
    autocorrelation_time = np.asarray(
        cast("Array", autocorrelation_time),
        dtype=float,
    )
    if (
        autocorrelation_time.ndim != 1
        or not np.all(np.isfinite(autocorrelation_time))
        or np.any(autocorrelation_time <= 0.0)
    ):
        return None
    return autocorrelation_time


def _estimate_autocorrelation_time(
    chain: Array,
) -> tuple[Array | None, Array | None, str | None]:
    try:
        autocorrelation_time = emcee_autocorr.integrated_time(chain, quiet=False)
    except emcee_autocorr.AutocorrError as error:
        tentative_autocorrelation_time = _valid_autocorrelation_time(error.tau)
        if tentative_autocorrelation_time is None:
            return None, None, "autocorrelation time unavailable"
        return (
            None,
            tentative_autocorrelation_time,
            "chain shorter than 50 times the autocorrelation time; "
            "tentative estimate reported",
        )
    except (FloatingPointError, ValueError):
        return None, None, "autocorrelation time unavailable"

    autocorrelation_time = _valid_autocorrelation_time(autocorrelation_time)
    if autocorrelation_time is None:
        return None, None, "autocorrelation time invalid"
    return autocorrelation_time, None, None


def _resolve_discarded_steps(
    burn: McmcBurnSetting,
    *,
    nsteps: int,
    autocorrelation_time: Array | None,
    autocorrelation_time_reliable: bool = True,
) -> tuple[int, str | None]:
    if burn != "auto":
        return int(burn), None
    if autocorrelation_time is None:
        return 0, "autocorrelation time unavailable; automatic burn-in was not applied"
    max_tau = float(np.max(autocorrelation_time))
    if not np.isfinite(max_tau) or max_tau <= 0.0:
        return 0, "autocorrelation time invalid; automatic burn-in was not applied"
    discarded_steps = ceil(2.0 * max_tau)
    if discarded_steps >= nsteps:
        return (
            0,
            "estimated automatic burn-in is longer than the chain; "
            "automatic burn-in was not applied",
        )
    if not autocorrelation_time_reliable:
        return (
            discarded_steps,
            "autocorrelation time estimate is unreliable; "
            "tentative automatic burn-in was applied",
        )
    return discarded_steps, None


def _apply_sample_window(
    chain: Array,
    lnprob: Array,
    *,
    burn: McmcBurnSetting,
    thin: int,
    autocorrelation_time: Array | None,
    autocorrelation_time_reliable: bool = True,
) -> tuple[Array, Array, int, str | None]:
    discarded_steps, warning = _resolve_discarded_steps(
        burn,
        nsteps=chain.shape[0],
        autocorrelation_time=autocorrelation_time,
        autocorrelation_time_reliable=autocorrelation_time_reliable,
    )
    retained_chain = chain[discarded_steps::thin]
    retained_lnprob = lnprob[discarded_steps::thin]
    if retained_chain.shape[0] < 1:
        msg = "MCMC settings did not retain any samples"
        raise ValueError(msg)
    return retained_chain, retained_lnprob, discarded_steps, warning


def _result_from_emcee(
    result: EmceeSamplerResult,
    settings: EffectiveMcmcSettings,
) -> McmcResult:
    var_names = result.var_names
    chain = result.chain
    lnprob = result.lnprob
    autocorrelation_time, tentative_autocorrelation_time, autocorrelation_warning = (
        _estimate_autocorrelation_time(chain)
    )
    burn_in_autocorrelation_time = (
        autocorrelation_time
        if autocorrelation_time is not None
        else tentative_autocorrelation_time
    )
    retained_chain, retained_lnprob, discarded_steps, burn_in_warning = (
        _apply_sample_window(
            chain,
            lnprob,
            burn=settings.burn,
            thin=settings.thin,
            autocorrelation_time=burn_in_autocorrelation_time,
            autocorrelation_time_reliable=autocorrelation_time is not None,
        )
    )
    samples = retained_chain.reshape((-1, len(var_names)))
    return McmcResult(
        var_names=var_names,
        chain=retained_chain,
        lnprob=retained_lnprob,
        summary=_summarize_chain(
            var_names,
            samples,
            autocorrelation_time,
            settings.thin,
        ),
        correlations=_correlation_matrix(samples),
        acceptance_fraction=result.acceptance_fraction,
        autocorrelation_time=autocorrelation_time,
        discarded_steps=discarded_steps,
        burn_in_warning=burn_in_warning,
        tentative_autocorrelation_time=tentative_autocorrelation_time,
        autocorrelation_warning=autocorrelation_warning,
        raw_chain=chain,
        raw_lnprob=lnprob,
    )


def _quote_toml_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _format_toml_string_list(values: list[str]) -> str:
    if not values:
        return "[]"
    return "[" + ", ".join(_quote_toml_string(value) for value in values) + "]"


def _format_toml_float(value: float) -> str:
    return f"{value:.5e}"


def _format_toml_float_list(values: list[float]) -> str:
    if not values:
        return "[]"
    return "[" + ", ".join(_format_toml_float(value) for value in values) + "]"


def _package_version(package_name: str) -> str:
    try:
        return version(package_name)
    except PackageNotFoundError:
        return "unknown"


def _autocorrelation_status(result: McmcResult) -> str:
    if result.autocorrelation_time is not None:
        return "reliable"
    if result.tentative_autocorrelation_time is not None:
        return "unreliable_short_chain"
    return "unavailable"


def _autocorrelation_steps(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
) -> int:
    if result.raw_chain is not None:
        return int(result.raw_chain.shape[0])
    return int(result.discarded_steps + result.chain.shape[0] * settings.thin)


def _autocorrelation_time_for_reporting(
    result: McmcResult,
) -> tuple[Array | None, bool]:
    if result.autocorrelation_time is not None:
        return result.autocorrelation_time, True
    if result.tentative_autocorrelation_time is not None:
        return result.tentative_autocorrelation_time, False
    return None, False


def _extend_autocorrelation_diagnostics(
    lines: list[str],
    result: McmcResult,
    settings: EffectiveMcmcSettings,
) -> None:
    autocorrelation_time, is_reliable = _autocorrelation_time_for_reporting(result)
    if autocorrelation_time is None:
        warning = result.autocorrelation_warning or "not available"
        lines.append(f"autocorrelation_warning = {_quote_toml_string(warning)}")
        return

    max_tau = float(np.max(autocorrelation_time))
    suffix = "" if is_reliable else "_tentative"
    tau_key = f"autocorrelation_time{suffix}"
    max_tau_key = f"max_autocorrelation_time{suffix}"
    lines.extend(
        [
            f"{tau_key} = {_format_toml_float_list(autocorrelation_time.tolist())}",
            f"{max_tau_key} = {_format_toml_float(max_tau)}",
        ],
    )
    if max_tau <= 0.0:
        lines.append('autocorrelation_warning = "autocorrelation time is not positive"')
        return

    autocorrelation_steps = _autocorrelation_steps(result, settings)
    lines.extend(
        [
            "steps_over_max_autocorrelation_time = "
            f"{_format_toml_float(autocorrelation_steps / max_tau)}",
            f"recommended_min_steps_50tau = {ceil(50.0 * max_tau)}",
            f"recommended_min_steps_100tau = {ceil(100.0 * max_tau)}",
        ],
    )
    if not is_reliable:
        warning = (
            result.autocorrelation_warning
            or "chain shorter than 50 times the autocorrelation time"
        )
        lines.extend(
            [
                f"autocorrelation_warning = {_quote_toml_string(warning)}",
                'effective_sample_size_warning = "not reported: autocorrelation time '
                'estimate is unreliable"',
            ],
        )
        return

    retained_steps_over_tau = result.chain.shape[0] * settings.thin / max_tau
    min_effective_sample_size = min(
        (
            summary.effective_sample_size
            for summary in result.summary
            if summary.effective_sample_size is not None
        ),
        default=None,
    )
    lines.append(
        "retained_steps_over_max_autocorrelation_time = "
        f"{_format_toml_float(retained_steps_over_tau)}",
    )
    if min_effective_sample_size is not None:
        lines.append(
            "min_effective_sample_size = "
            f"{_format_toml_float(min_effective_sample_size)}",
        )
    if retained_steps_over_tau < 50.0:
        lines.append(
            'autocorrelation_warning = "Retained chain length is shorter '
            'than 50 times the maximum autocorrelation time."',
        )


def _extend_timing_diagnostics(
    lines: list[str],
    timings: Mapping[str, float],
) -> None:
    lines.extend(
        f"{key} = {_format_toml_float(timings[key])}"
        for key in _TIMING_KEYS
        if key in timings
    )


def _write_summary(
    result: McmcResult,
    path: Path,
    parameter_store: ParameterStore,
) -> None:
    parameters = parameter_store.get_parameters(result.var_names)
    lines: list[str] = []
    for summary in result.summary:
        parameter = parameters[summary.parameter_id]
        parameter_name = str(parameter.param_name).strip("[]")
        lines.extend(
            [
                f"[{_quote_toml_string(parameter_name)}]",
                'prior = "uniform"',
                f"prior_lower = {_format_toml_float(float(parameter.min))}",
                f"prior_upper = {_format_toml_float(float(parameter.max))}",
                'credible_interval = "95% equal-tailed"',
                f"mean = {_format_toml_float(summary.mean)}",
                f"standard_deviation = {_format_toml_float(summary.standard_deviation)}",
                f"median = {_format_toml_float(summary.median)}",
                f"eti_95_lower = {_format_toml_float(summary.eti_95_lower)}",
                f"eti_95_upper = {_format_toml_float(summary.eti_95_upper)}",
                f"lower_1sigma = {_format_toml_float(summary.lower_1sigma)}",
                f"upper_1sigma = {_format_toml_float(summary.upper_1sigma)}",
                f"stderr = {_format_toml_float(summary.stderr)}",
            ],
        )
        if summary.effective_sample_size is not None:
            lines.append(
                "effective_sample_size = "
                f"{_format_toml_float(summary.effective_sample_size)}",
            )
        if summary.mcse_mean is not None:
            lines.append(f"mcse_mean = {_format_toml_float(summary.mcse_mean)}")
        lines.append("")
    (path / "summary.toml").write_text("\n".join(lines), encoding="utf-8")


def _write_samples(
    result: McmcResult,
    path: Path,
    parameter_store: ParameterStore,
) -> None:
    parameter_names = _format_parameter_ids(result.var_names, parameter_store)
    values = np.column_stack((result.samples, result.log_probabilities))

    with (path / "samples.tsv").open("w", encoding="utf-8") as fileout:
        fileout.write("\t".join((*parameter_names, "lnprob")) + "\n")
        np.savetxt(fileout, values, fmt="%.5e", delimiter="\t")


def _write_correlations(
    result: McmcResult,
    path: Path,
    parameter_store: ParameterStore,
) -> None:
    parameter_names = _format_parameter_ids(result.var_names, parameter_store)

    with (path / "correlations.tsv").open("w", encoding="utf-8") as fileout:
        fileout.write("parameter\t" + "\t".join(parameter_names) + "\n")
        for name, values in zip(
            parameter_names,
            result.correlations,
            strict=True,
        ):
            row = "\t".join(f"{value:.5e}" for value in values)
            fileout.write(f"{name}\t{row}\n")


def _write_diagnostics(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    path: Path,
    parameter_store: ParameterStore,
    unbounded_parameter_ids: tuple[str, ...],
    timings: Mapping[str, float],
) -> None:
    acceptance = result.acceptance_fraction
    unbounded_parameters = list(
        _format_parameter_ids(unbounded_parameter_ids, parameter_store),
    )
    requested_burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        'sampler = "emcee via ChemEx direct EnsembleSampler"',
        f"lmfit_version = {_quote_toml_string(_package_version('lmfit'))}",
        f"emcee_version = {_quote_toml_string(_package_version('emcee'))}",
        'credible_interval = "95% equal-tailed"',
        'convergence_diagnostic = "integrated_autocorrelation_time"',
        f"autocorrelation_status = {_quote_toml_string(_autocorrelation_status(result))}",
        'rhat = "not computed: emcee ensemble walkers are not independent chains"',
        f"steps = {settings.steps}",
        f"requested_burn = {requested_burn}",
        f"discarded_steps = {result.discarded_steps}",
        f"thin = {settings.thin}",
        f"walkers = {settings.walkers}",
        f"workers = {settings.workers}",
        f"retained_steps = {result.chain.shape[0]}",
        f"retained_samples = {len(result.samples)}",
        'samples_file = "samples.tsv"',
        'summary_file = "summary.toml"',
        'correlations_file = "correlations.tsv"',
        'plots_file = "plots.pdf"',
        f"acceptance_fraction_mean = {_format_toml_float(float(np.mean(acceptance)))}",
        f"acceptance_fraction_min = {_format_toml_float(float(np.min(acceptance)))}",
        f"acceptance_fraction_max = {_format_toml_float(float(np.max(acceptance)))}",
        f"unbounded_parameters = {_format_toml_string_list(unbounded_parameters)}",
    ]
    if result.burn_in_warning is not None:
        lines.append(f"burn_in_warning = {_quote_toml_string(result.burn_in_warning)}")
    _extend_autocorrelation_diagnostics(lines, result, settings)
    _extend_timing_diagnostics(lines, timings)
    if unbounded_parameters:
        lines.append(
            'warning = "Finite lower and upper bounds are recommended for MCMC."',
        )

    (path / "diagnostics.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_mcmc_outputs(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    path: Path,
    parameter_store: ParameterStore,
    unbounded_parameter_ids: tuple[str, ...] = (),
    timings: dict[str, float] | None = None,
) -> None:
    timings = {} if timings is None else timings
    output_start = perf_counter()
    path_mcmc = path / "Statistics" / "MCMC"
    path_mcmc.mkdir(parents=True, exist_ok=True)

    phase_start = perf_counter()
    _write_summary(result, path_mcmc, parameter_store)
    timings["output_summary_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    _write_samples(result, path_mcmc, parameter_store)
    timings["output_samples_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    _write_correlations(result, path_mcmc, parameter_store)
    timings["output_correlations_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    write_mcmc_plots(
        result,
        settings,
        path_mcmc,
        parameter_names=_format_parameter_ids(result.var_names, parameter_store),
    )
    timings["output_plots_seconds"] = perf_counter() - phase_start
    timings["output_total_seconds"] = perf_counter() - output_start
    if "sampling_seconds" in timings and "result_processing_seconds" in timings:
        timings["total_seconds"] = (
            timings["sampling_seconds"]
            + timings["result_processing_seconds"]
            + timings["output_total_seconds"]
        )
    _write_diagnostics(
        result,
        settings,
        path_mcmc,
        parameter_store,
        unbounded_parameter_ids,
        timings,
    )


def run_mcmc(
    experiments: Experiments,
    params: Parameters,
    settings: McmcSettings,
    path: Path,
    *,
    execution: ExecutionSettings | None = None,
) -> McmcResult | None:
    var_names = _varying_parameter_ids(params)
    if not var_names:
        print_mcmc_no_vary_warning()
        return None

    effective_settings = resolve_mcmc_settings(
        settings,
        nvarys=len(var_names),
        execution=execution,
    )
    _validate_parameter_bounds(params, var_names)

    parameter_store = experiments.parameter_store
    unbounded_parameter_ids = _find_unbounded_parameter_ids(params, var_names)
    if unbounded_parameter_ids:
        print_mcmc_unbounded_warning(
            list(_format_parameter_ids(unbounded_parameter_ids, parameter_store)),
        )

    timings: dict[str, float] = {}
    phase_start = perf_counter()
    with _use_rich_emcee_progress():
        env = native_thread_env(
            effective_settings.native_threads,
            parallel=effective_settings.workers > 1,
        )
        with native_thread_environment(env):
            result_emcee = run_emcee_sampler(
                McmcProblem(
                    experiments=experiments,
                    params=params,
                    var_names=var_names,
                    bounds=_parameter_bounds(params, var_names),
                ),
                _initial_positions(
                    params,
                    var_names,
                    nwalkers=effective_settings.walkers,
                    seed=effective_settings.seed,
                ),
                steps=effective_settings.steps,
                workers=effective_settings.workers,
                seed=effective_settings.seed,
                progress=True,
            )
    timings["sampling_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    result = _result_from_emcee(result_emcee, effective_settings)
    timings["result_processing_seconds"] = perf_counter() - phase_start
    write_mcmc_outputs(
        result,
        effective_settings,
        path,
        parameter_store,
        unbounded_parameter_ids,
        timings=timings,
    )
    if effective_settings.update_parameters:
        updated_params = params.copy()
        for summary in result.summary:
            updated_params[summary.parameter_id].value = summary.median
            updated_params[summary.parameter_id].stderr = summary.stderr
        updated_params.update_constraints()
        parameter_store.update_from_parameters(updated_params)
    return result
