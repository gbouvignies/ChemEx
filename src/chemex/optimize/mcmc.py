from __future__ import annotations

import secrets
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from math import ceil
from pathlib import Path
from time import perf_counter
from typing import cast

import emcee.autocorr as emcee_autocorr
import numpy as np

from chemex.configuration.methods import McmcBurnSetting, McmcSettings
from chemex.containers.experiments import Experiments
from chemex.optimize.native_deterministic import NativeDeterministicFit
from chemex.optimize.native_mcmc import (
    McmcDiagnosticStatus,
    McmcEvidence,
    McmcOperationTerminal,
    McmcPlan,
    derive_mcmc_diagnostics,
    execute_mcmc_evidence,
    resolve_product_mcmc_policy,
)
from chemex.optimize.statistics_plot import write_mcmc_plots
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.runtime import ExecutionSettings
from chemex.typing import Array


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
    credible_interval_68_lower: float
    credible_interval_68_upper: float
    half_credible_interval_68_width: float
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


class NativeMcmcIncompleteError(RuntimeError):
    """Raised when native product MCMC cannot publish a complete chain."""

    def __init__(
        self,
        message: str,
        *,
        terminal: str = "failed",
        completed_steps: int = 0,
    ) -> None:
        super().__init__(message)
        self.terminal = terminal
        self.completed_steps = completed_steps


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


def _format_parameter_ids(
    parameter_ids: tuple[str, ...],
    parameter_model: SealedParameterModel,
) -> tuple[str, ...]:
    return tuple(
        str(parameter_name_from_definition(parameter_model.definitions[param_id]))
        for param_id in parameter_ids
    )


def _summarize_chain(
    var_names: tuple[str, ...],
    samples: Array,
    autocorrelation_time: Array | None,
    thin: int,
) -> tuple[McmcSummary, ...]:
    quantiles = np.percentile(samples, [2.5, 15.87, 50.0, 84.13, 97.5], axis=0)
    summaries: list[McmcSummary] = []
    for index, name in enumerate(var_names):
        lower_95, lower_68, median, upper_68, upper_95 = quantiles[:, index]
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
                credible_interval_68_lower=float(lower_68),
                credible_interval_68_upper=float(upper_68),
                half_credible_interval_68_width=0.5 * float(upper_68 - lower_68),
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
            (
                "chain shorter than 50 times the autocorrelation time; "
                "tentative estimate reported"
            ),
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
            (
                "estimated automatic burn-in is longer than the chain; "
                "automatic burn-in was not applied"
            ),
        )
    if not autocorrelation_time_reliable:
        return (
            discarded_steps,
            (
                "autocorrelation time estimate is unreliable; "
                "tentative automatic burn-in was applied"
            ),
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
            (
                "steps_over_max_autocorrelation_time = "
                f"{_format_toml_float(autocorrelation_steps / max_tau)}"
            ),
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
                (
                    'effective_sample_size_warning = "not reported: autocorrelation time '
                    'estimate is unreliable"'
                ),
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
    parameter_model: SealedParameterModel,
) -> None:
    lines: list[str] = []
    for summary in result.summary:
        param_id = summary.parameter_id
        parameter_name = str(
            parameter_name_from_definition(parameter_model.definitions[param_id])
        ).strip("[]")
        configuration = parameter_model.configuration[param_id]
        lines.extend(
            [
                f"[{_quote_toml_string(parameter_name)}]",
                'prior = "uniform"',
                f"prior_lower = {_format_toml_float(configuration.lower_bound)}",
                f"prior_upper = {_format_toml_float(configuration.upper_bound)}",
                'credible_interval = "95% equal-tailed"',
                f"mean = {_format_toml_float(summary.mean)}",
                f"standard_deviation = {_format_toml_float(summary.standard_deviation)}",
                f"median = {_format_toml_float(summary.median)}",
                f"eti_95_lower = {_format_toml_float(summary.eti_95_lower)}",
                f"eti_95_upper = {_format_toml_float(summary.eti_95_upper)}",
                (
                    "credible_interval_68_lower = "
                    f"{_format_toml_float(summary.credible_interval_68_lower)}"
                ),
                (
                    "credible_interval_68_upper = "
                    f"{_format_toml_float(summary.credible_interval_68_upper)}"
                ),
                (
                    "half_credible_interval_68_width = "
                    f"{_format_toml_float(summary.half_credible_interval_68_width)}"
                ),
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
    parameter_model: SealedParameterModel,
) -> None:
    parameter_names = _format_parameter_ids(result.var_names, parameter_model)
    values = np.column_stack((result.samples, result.log_probabilities))

    with (path / "samples.tsv").open("w", encoding="utf-8") as fileout:
        fileout.write("\t".join((*parameter_names, "lnprob")) + "\n")
        np.savetxt(fileout, values, fmt="%.5e", delimiter="\t")


def _write_correlations(
    result: McmcResult,
    path: Path,
    parameter_model: SealedParameterModel,
) -> None:
    parameter_names = _format_parameter_ids(result.var_names, parameter_model)

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
    timings: Mapping[str, float],
    *,
    engine: str,
    root_seed: int | None = None,
) -> None:
    acceptance = result.acceptance_fraction
    requested_burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        'status = "complete"',
        f"engine = {_quote_toml_string(engine)}",
        'sampler = "emcee via ChemEx direct EnsembleSampler"',
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
        "unbounded_parameters = []",
    ]
    if root_seed is not None:
        lines.append(f"root_seed = {root_seed}")
    if result.burn_in_warning is not None:
        lines.append(f"burn_in_warning = {_quote_toml_string(result.burn_in_warning)}")
    _extend_autocorrelation_diagnostics(lines, result, settings)
    _extend_timing_diagnostics(lines, timings)
    (path / "diagnostics.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_mcmc_outputs(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    path: Path,
    parameter_model: SealedParameterModel,
    timings: dict[str, float] | None = None,
    *,
    engine: str = "native MCMC",
    root_seed: int | None = None,
) -> None:
    timings = {} if timings is None else timings
    output_start = perf_counter()
    path_mcmc = path / "Statistics" / "MCMC"
    path_mcmc.mkdir(parents=True, exist_ok=True)

    phase_start = perf_counter()
    _write_summary(result, path_mcmc, parameter_model)
    timings["output_summary_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    _write_samples(result, path_mcmc, parameter_model)
    timings["output_samples_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    _write_correlations(result, path_mcmc, parameter_model)
    timings["output_correlations_seconds"] = perf_counter() - phase_start

    phase_start = perf_counter()
    write_mcmc_plots(
        result,
        settings,
        path_mcmc,
        parameter_names=_format_parameter_ids(result.var_names, parameter_model),
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
        timings,
        engine=engine,
        root_seed=root_seed,
    )


def _clear_mcmc_artifacts(path: Path) -> None:
    for name in (
        "summary.toml",
        "samples.tsv",
        "correlations.tsv",
        "diagnostics.toml",
        "plots.pdf",
    ):
        (path / name).unlink(missing_ok=True)


def _write_native_mcmc_state_diagnostics(
    path: Path,
    *,
    settings: EffectiveMcmcSettings,
    root_seed: int,
    parameter_ids: tuple[str, ...],
    status: str,
    terminal: str | None = None,
    completed_steps: int | None = None,
    failure: BaseException | None = None,
) -> None:
    requested_burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        f"status = {_quote_toml_string(status)}",
        'engine = "native MCMC"',
        'sampler = "emcee via ChemEx direct EnsembleSampler"',
        f"steps = {settings.steps}",
        f"requested_burn = {requested_burn}",
        f"thin = {settings.thin}",
        f"walkers = {settings.walkers}",
        f"workers = {settings.workers}",
        f"root_seed = {root_seed}",
        f"parameters = {_format_toml_string_list(list(parameter_ids))}",
    ]
    if terminal is not None:
        lines.append(f"terminal = {_quote_toml_string(terminal)}")
    if completed_steps is not None:
        lines.append(f"completed_steps = {completed_steps}")
    if failure is not None:
        message = str(failure).replace("\n", " ")
        lines.extend(
            (
                f"failure_type = {_quote_toml_string(type(failure).__name__)}",
                f"failure_message = {_quote_toml_string(message)}",
            )
        )
    (path / "diagnostics.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _native_result_from_evidence(
    evidence: McmcEvidence,
    settings: EffectiveMcmcSettings,
) -> McmcResult:
    raw_chain = np.asarray(
        [state.positions for state in evidence.states[1:]],
        dtype=float,
    )
    raw_lnprob = np.asarray(
        [state.log_densities for state in evidence.states[1:]],
        dtype=float,
    )
    diagnostics = derive_mcmc_diagnostics(evidence)
    if (
        diagnostics.status is not McmcDiagnosticStatus.AVAILABLE
        or diagnostics.acceptance_fractions is None
    ):
        msg = "Native MCMC completed without authoritative acceptance diagnostics"
        raise NativeMcmcIncompleteError(
            msg,
            completed_steps=evidence.completed_transition_count,
        )
    autocorrelation_time, tentative_autocorrelation_time, autocorrelation_warning = (
        _estimate_autocorrelation_time(raw_chain)
    )
    burn_autocorrelation_time = (
        autocorrelation_time
        if autocorrelation_time is not None
        else tentative_autocorrelation_time
    )
    retained_chain, retained_lnprob, discarded_steps, burn_in_warning = (
        _apply_sample_window(
            raw_chain,
            raw_lnprob,
            burn=settings.burn,
            thin=settings.thin,
            autocorrelation_time=burn_autocorrelation_time,
            autocorrelation_time_reliable=autocorrelation_time is not None,
        )
    )
    samples = retained_chain.reshape((-1, len(evidence.coordinate_ids)))
    return McmcResult(
        var_names=evidence.coordinate_ids,
        chain=retained_chain,
        lnprob=retained_lnprob,
        summary=_summarize_chain(
            evidence.coordinate_ids,
            samples,
            autocorrelation_time,
            settings.thin,
        ),
        correlations=_correlation_matrix(samples),
        acceptance_fraction=np.asarray(
            diagnostics.acceptance_fractions,
            dtype=float,
        ),
        autocorrelation_time=autocorrelation_time,
        discarded_steps=discarded_steps,
        burn_in_warning=burn_in_warning,
        tentative_autocorrelation_time=tentative_autocorrelation_time,
        autocorrelation_warning=autocorrelation_warning,
        raw_chain=raw_chain,
        raw_lnprob=raw_lnprob,
    )


def run_native_mcmc(
    experiments: Experiments,
    fit: NativeDeterministicFit,
    settings: McmcSettings,
    path: Path,
    *,
    execution: ExecutionSettings | None = None,
) -> McmcResult:
    """Run product MCMC directly from one committed native deterministic fit."""
    effective = resolve_mcmc_settings(
        settings,
        nvarys=len(fit.problem.controlled_ids),
        execution=execution,
    )
    root_seed = secrets.randbits(64) if effective.seed is None else effective.seed
    statistic_path = path / "Statistics" / "MCMC"
    statistic_path.mkdir(parents=True, exist_ok=True)
    _clear_mcmc_artifacts(statistic_path)
    _write_native_mcmc_state_diagnostics(
        statistic_path,
        settings=effective,
        root_seed=root_seed,
        parameter_ids=fit.problem.controlled_ids,
        status="running",
    )
    lower = np.asarray(fit.problem.lower_bounds, dtype=float)
    upper = np.asarray(fit.problem.upper_bounds, dtype=float)
    invalid_bounds = ~np.isfinite(lower) | ~np.isfinite(upper) | (lower >= upper)
    if np.any(invalid_bounds):
        parameter_names = _format_parameter_ids(
            tuple(
                param_id
                for param_id, invalid in zip(
                    fit.problem.controlled_ids,
                    invalid_bounds,
                    strict=True,
                )
                if invalid
            ),
            fit.parameter_model,
        )
        error = ValueError(
            "Native MCMC requires finite lower and upper bounds with lower < upper "
            f"for every fitted parameter; affected: {', '.join(parameter_names)}"
        )
        _write_native_mcmc_state_diagnostics(
            statistic_path,
            settings=effective,
            root_seed=root_seed,
            parameter_ids=fit.problem.controlled_ids,
            status="incomplete",
            terminal="failed",
            completed_steps=0,
            failure=error,
        )
        raise error

    timings: dict[str, float] = {}
    completed_steps = 0
    try:
        policy = resolve_product_mcmc_policy(
            dimension=len(fit.problem.controlled_ids),
            walkers=effective.walkers,
            steps=effective.steps,
            root_seed=root_seed,
        )
        plan = McmcPlan.for_accepted(
            fit.accepted,
            source_problem=fit.problem,
            parameterization=fit.parameterization,
            source_engine=fit.engine,
            policy=policy,
            coordinate_units=tuple(
                (param_id, ParameterUnit.UNSPECIFIED)
                for param_id in fit.problem.controlled_ids
            ),
        )
        phase_start = perf_counter()
        operation = execute_mcmc_evidence(
            fit.accepted,
            plan,
            execution=ExecutionSettings(
                workers=effective.workers,
                native_threads=effective.native_threads,
            ),
        )
        timings["sampling_seconds"] = perf_counter() - phase_start
        completed_steps = (
            0
            if operation.evidence is None
            else operation.evidence.completed_transition_count
        )
        if (
            operation.terminal is not McmcOperationTerminal.COMPLETED
            or operation.evidence is None
        ):
            error = NativeMcmcIncompleteError(
                f"Native MCMC {operation.terminal.value}: {operation.failure_message}",
                terminal=operation.terminal.value,
                completed_steps=completed_steps,
            )
            raise error
        phase_start = perf_counter()
        result = _native_result_from_evidence(operation.evidence, effective)
        timings["result_processing_seconds"] = perf_counter() - phase_start
        write_mcmc_outputs(
            result,
            effective,
            path,
            fit.parameter_model,
            timings=timings,
            engine="native MCMC",
            root_seed=root_seed,
        )
    except (Exception, KeyboardInterrupt) as error:
        _clear_mcmc_artifacts(statistic_path)
        terminal = (
            error.terminal
            if isinstance(error, NativeMcmcIncompleteError)
            else "interrupted"
            if isinstance(error, KeyboardInterrupt)
            else "failed"
        )
        failure_completed_steps = (
            max(error.completed_steps, completed_steps)
            if isinstance(error, NativeMcmcIncompleteError)
            else completed_steps
        )
        _write_native_mcmc_state_diagnostics(
            statistic_path,
            settings=effective,
            root_seed=root_seed,
            parameter_ids=fit.problem.controlled_ids,
            status="incomplete",
            terminal=terminal,
            completed_steps=failure_completed_steps,
            failure=error,
        )
        raise
    else:
        return result
