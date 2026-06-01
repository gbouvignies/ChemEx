from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from math import ceil
from pathlib import Path
from typing import cast

import lmfit
import numpy as np
from lmfit.minimizer import MinimizerResult
from lmfit.parameter import Parameters

from chemex.configuration.methods import McmcBurnSetting, McmcSettings
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_mcmc_no_vary_warning,
    print_mcmc_unbounded_warning,
)
from chemex.parameters.database import ParameterStore
from chemex.typing import Array


@dataclass(frozen=True)
class EffectiveMcmcSettings:
    steps: int
    burn: McmcBurnSetting
    thin: int
    walkers: int
    seed: int | None
    workers: int
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


def resolve_mcmc_settings(
    settings: McmcSettings,
    *,
    nvarys: int,
) -> EffectiveMcmcSettings:
    walkers = settings.walkers or max(32, 2 * nvarys)
    min_walkers = 2 * nvarys
    if walkers < min_walkers:
        msg = (
            f"MCMC requires at least {min_walkers} walkers for {nvarys} fitted "
            f"parameters"
        )
        raise ValueError(msg)

    return EffectiveMcmcSettings(
        steps=settings.steps,
        burn=settings.burn,
        thin=settings.thin,
        walkers=walkers,
        seed=settings.seed,
        workers=settings.workers,
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


def _result_attr(result: MinimizerResult, name: str) -> object:
    return getattr(result, name)


def _autocorrelation_time(result: MinimizerResult) -> Array | None:
    if not hasattr(result, "acor"):
        return None
    autocorrelation_time = np.asarray(
        cast("Array", _result_attr(result, "acor")),
        dtype=float,
    )
    if autocorrelation_time.ndim != 1 or not np.all(np.isfinite(autocorrelation_time)):
        return None
    return autocorrelation_time


def _resolve_discarded_steps(
    burn: McmcBurnSetting,
    *,
    nsteps: int,
    autocorrelation_time: Array | None,
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
    return discarded_steps, None


def _apply_sample_window(
    chain: Array,
    lnprob: Array,
    *,
    burn: McmcBurnSetting,
    thin: int,
    autocorrelation_time: Array | None,
) -> tuple[Array, Array, int, str | None]:
    discarded_steps, warning = _resolve_discarded_steps(
        burn,
        nsteps=chain.shape[0],
        autocorrelation_time=autocorrelation_time,
    )
    retained_chain = chain[discarded_steps::thin]
    retained_lnprob = lnprob[discarded_steps::thin]
    if retained_chain.shape[0] < 1:
        msg = "MCMC settings did not retain any samples"
        raise ValueError(msg)
    return retained_chain, retained_lnprob, discarded_steps, warning


def _result_from_lmfit(
    result: MinimizerResult,
    settings: EffectiveMcmcSettings,
) -> McmcResult:
    var_names = tuple(cast("Sequence[str]", _result_attr(result, "var_names")))
    chain = np.asarray(cast("Array", _result_attr(result, "chain")), dtype=float)
    lnprob = np.asarray(cast("Array", _result_attr(result, "lnprob")), dtype=float)
    autocorrelation_time = _autocorrelation_time(result)
    retained_chain, retained_lnprob, discarded_steps, burn_in_warning = (
        _apply_sample_window(
            chain,
            lnprob,
            burn=settings.burn,
            thin=settings.thin,
            autocorrelation_time=autocorrelation_time,
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
        acceptance_fraction=np.asarray(
            cast("Array", _result_attr(result, "acceptance_fraction")),
            dtype=float,
        ),
        autocorrelation_time=autocorrelation_time,
        discarded_steps=discarded_steps,
        burn_in_warning=burn_in_warning,
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
) -> None:
    acceptance = result.acceptance_fraction
    unbounded_parameters = list(
        _format_parameter_ids(unbounded_parameter_ids, parameter_store),
    )
    requested_burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        'sampler = "emcee via lmfit.Minimizer.emcee"',
        f"lmfit_version = {_quote_toml_string(_package_version('lmfit'))}",
        f"emcee_version = {_quote_toml_string(_package_version('emcee'))}",
        'credible_interval = "95% equal-tailed"',
        'convergence_diagnostic = "integrated_autocorrelation_time"',
        'rhat = "not computed: emcee ensemble walkers are not independent chains"',
        f"steps = {settings.steps}",
        f"requested_burn = {requested_burn}",
        f"discarded_steps = {result.discarded_steps}",
        f"thin = {settings.thin}",
        f"walkers = {settings.walkers}",
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
    if result.autocorrelation_time is None:
        lines.append('autocorrelation_warning = "not available"')
    else:
        max_tau = float(np.max(result.autocorrelation_time))
        lines.append(
            "autocorrelation_time = "
            f"{_format_toml_float_list(result.autocorrelation_time.tolist())}",
        )
        lines.append(f"max_autocorrelation_time = {_format_toml_float(max_tau)}")
        if max_tau <= 0.0:
            lines.append(
                'autocorrelation_warning = "autocorrelation time is not positive"',
            )
        else:
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
) -> None:
    path_mcmc = path / "Statistics" / "MCMC"
    path_mcmc.mkdir(parents=True, exist_ok=True)
    _write_summary(result, path_mcmc, parameter_store)
    _write_samples(result, path_mcmc, parameter_store)
    _write_correlations(result, path_mcmc, parameter_store)
    from chemex.optimize.mcmc_plot import write_mcmc_plots

    write_mcmc_plots(
        result,
        settings,
        path_mcmc,
        parameter_names=_format_parameter_ids(result.var_names, parameter_store),
    )
    _write_diagnostics(
        result,
        settings,
        path_mcmc,
        parameter_store,
        unbounded_parameter_ids,
    )


def run_mcmc(
    experiments: Experiments,
    params: Parameters,
    settings: McmcSettings,
    path: Path,
) -> McmcResult | None:
    var_names = _varying_parameter_ids(params)
    if not var_names:
        print_mcmc_no_vary_warning()
        return None

    effective_settings = resolve_mcmc_settings(settings, nvarys=len(var_names))
    _validate_parameter_bounds(params, var_names)

    parameter_store = experiments.parameter_store
    unbounded_parameter_ids = _find_unbounded_parameter_ids(params, var_names)
    if unbounded_parameter_ids:
        print_mcmc_unbounded_warning(
            list(_format_parameter_ids(unbounded_parameter_ids, parameter_store)),
        )

    minimizer = lmfit.Minimizer(experiments.residuals, params)
    result_lf = minimizer.emcee(
        params=params,
        steps=effective_settings.steps,
        nwalkers=effective_settings.walkers,
        burn=0,
        thin=1,
        pos=_initial_positions(
            params,
            var_names,
            nwalkers=effective_settings.walkers,
            seed=effective_settings.seed,
        ),
        workers=effective_settings.workers,
        is_weighted=True,
        seed=effective_settings.seed,
        progress=False,
    )
    result = _result_from_lmfit(result_lf, effective_settings)
    write_mcmc_outputs(
        result,
        effective_settings,
        path,
        parameter_store,
        unbounded_parameter_ids,
    )
    if effective_settings.update_parameters:
        for summary in result.summary:
            result_lf.params[summary.parameter_id].value = summary.median
            result_lf.params[summary.parameter_id].stderr = summary.stderr
        result_lf.params.update_constraints()
        parameter_store.update_from_parameters(result_lf.params)
    return result
