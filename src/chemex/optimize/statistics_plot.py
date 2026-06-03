from __future__ import annotations

from itertools import combinations
from math import ceil
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from emcee.autocorr import AutocorrError, integrated_time
from matplotlib.backends.backend_pdf import PdfPages

from chemex.typing import Array

if TYPE_CHECKING:
    from chemex.optimize.mcmc import EffectiveMcmcSettings, McmcResult


_CORRELATION_THRESHOLD = 0.5
_PAGE_SIZE = (8.27, 11.69)


def _format_float(value: float | None) -> str:
    if value is None:
        return "n/a"
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.5g}"


def _finite_values(values: Array) -> Array:
    return values[np.isfinite(values)]


def _mcmc_autocorrelation_time(result: McmcResult) -> Array | None:
    if result.autocorrelation_time is not None:
        return result.autocorrelation_time
    return result.tentative_autocorrelation_time


def _mcmc_autocorrelation_status(result: McmcResult) -> str:
    if result.autocorrelation_time is not None:
        return "reliable"
    if result.tentative_autocorrelation_time is not None:
        return "unreliable short chain"
    return "unavailable"


def _correlated_pairs(
    correlations: Array,
    *,
    threshold: float,
) -> list[tuple[int, int, float]]:
    pairs: list[tuple[int, int, float]] = []
    for index_x, index_y in combinations(range(correlations.shape[0]), 2):
        correlation = float(correlations[index_x, index_y])
        if np.isfinite(correlation) and abs(correlation) >= threshold:
            pairs.append((index_x, index_y, correlation))
    pairs.sort(key=lambda pair: abs(pair[2]), reverse=True)
    return pairs


def _save_text_page(pdf: PdfPages, title: str, lines: list[str]) -> None:
    fig = plt.figure(figsize=_PAGE_SIZE)
    fig.text(0.08, 0.95, title, fontsize=16, weight="bold", va="top")
    fig.text(
        0.08,
        0.90,
        "\n".join(lines),
        family="monospace",
        fontsize=9,
        va="top",
    )
    pdf.savefig(fig)
    plt.close(fig)


def _save_histogram_page(
    pdf: PdfPages,
    *,
    values: Array,
    title: str,
    xlabel: str,
    ylabel: str,
    color: str,
    density: bool,
    interval: tuple[float, float] | None = None,
    interval_label: str | None = None,
    median: float | None = None,
    empty_message: str = "No finite samples",
) -> None:
    values = _finite_values(values)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    if len(values):
        ax.hist(values, bins="auto", density=density, color=color, alpha=0.85)
        if interval is not None:
            lower, upper = interval
            ax.axvspan(
                lower,
                upper,
                color=color,
                alpha=0.15,
                label=interval_label,
            )
        if median is not None:
            ax.axvline(
                median,
                color="#C44E52",
                linestyle="--",
                linewidth=1.2,
                label="median",
            )
        handles, _labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
    else:
        ax.text(0.5, 0.5, empty_message, ha="center", va="center")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _save_pair_distribution_page(
    pdf: PdfPages,
    *,
    samples: Array,
    parameter_names: tuple[str, ...],
    pair: tuple[int, int, float],
    title_prefix: str,
    cmap: str,
) -> None:
    index_x, index_y, correlation = pair
    values_x = samples[:, index_x]
    values_y = samples[:, index_y]
    finite = np.isfinite(values_x) & np.isfinite(values_y)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    if np.any(finite):
        image = ax.hist2d(
            values_x[finite],
            values_y[finite],
            bins=40,
            cmap=cmap,
        )
        cbar = fig.colorbar(image[3], ax=ax)
        cbar.set_label("Samples")
    else:
        ax.text(0.5, 0.5, "No finite samples", ha="center", va="center")
    ax.set_title(
        f"{title_prefix}: {parameter_names[index_x]} vs "
        f"{parameter_names[index_y]} (r = {correlation:.3f})",
    )
    ax.set_xlabel(parameter_names[index_x])
    ax.set_ylabel(parameter_names[index_y])
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _save_mcmc_summary_page(
    pdf: PdfPages,
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    parameter_names: tuple[str, ...],
    correlated_pairs: list[tuple[int, int, float]],
) -> None:
    acceptance = result.acceptance_fraction
    burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        "Sampler: emcee via ChemEx direct EnsembleSampler",
        f"Parameters: {len(parameter_names)}",
        f"Steps: {settings.steps}",
        f"Walkers: {settings.walkers}",
        f"Requested burn: {burn}",
        f"Discarded steps: {result.discarded_steps}",
        f"Thin: {settings.thin}",
        f"Retained steps: {result.chain.shape[0]}",
        f"Retained samples: {len(result.samples)}",
        (
            "Acceptance fraction: "
            f"mean={_format_float(float(np.mean(acceptance)))}, "
            f"min={_format_float(float(np.min(acceptance)))}, "
            f"max={_format_float(float(np.max(acceptance)))}"
        ),
    ]
    autocorrelation_time = _mcmc_autocorrelation_time(result)
    lines.append(f"Autocorrelation status: {_mcmc_autocorrelation_status(result)}")
    if autocorrelation_time is None:
        lines.append("Autocorrelation time: unavailable")
    else:
        max_tau = float(np.max(autocorrelation_time))
        lines.append(f"Max autocorrelation time: {_format_float(max_tau)}")
        lines.append(f"Recommended steps (50 tau): {ceil(50.0 * max_tau)}")
        lines.append(f"Recommended steps (100 tau): {ceil(100.0 * max_tau)}")
    if result.autocorrelation_warning is not None:
        lines.append(f"Autocorrelation warning: {result.autocorrelation_warning}")
    if result.burn_in_warning is not None:
        lines.append(f"Burn-in warning: {result.burn_in_warning}")

    lines.extend(["", "Posterior summaries:"])
    for name, summary in zip(parameter_names, result.summary, strict=True):
        lines.append(
            f"{name}: median={_format_float(summary.median)}, "
            f"95% ETI=[{_format_float(summary.eti_95_lower)}, "
            f"{_format_float(summary.eti_95_upper)}], "
            f"ESS={_format_float(summary.effective_sample_size)}",
        )

    lines.extend(["", f"Correlation threshold: |r| >= {_CORRELATION_THRESHOLD:g}"])
    if correlated_pairs:
        for index_x, index_y, correlation in correlated_pairs:
            lines.append(
                f"{parameter_names[index_x]} vs {parameter_names[index_y]}: "
                f"r={correlation:.3f}",
            )
    else:
        lines.append("No parameter pairs exceeded the correlation threshold.")

    _save_text_page(pdf, "MCMC Diagnostics", lines)


def _save_mcmc_distribution_pages(
    pdf: PdfPages,
    result: McmcResult,
    parameter_names: tuple[str, ...],
) -> None:
    samples = result.samples
    for index, (name, summary) in enumerate(
        zip(parameter_names, result.summary, strict=True),
    ):
        _save_histogram_page(
            pdf,
            values=samples[:, index],
            title=f"Posterior distribution: {name}",
            xlabel=name,
            ylabel="Density",
            color="#4C72B0",
            density=True,
            interval=(summary.eti_95_lower, summary.eti_95_upper),
            interval_label="95% ETI",
            median=summary.median,
        )


def _trace_chain(result: McmcResult) -> tuple[Array, bool]:
    if result.raw_chain is not None:
        return result.raw_chain, True
    return result.chain, False


def _trace_lnprob(result: McmcResult) -> tuple[Array, bool]:
    if result.raw_lnprob is not None:
        return result.raw_lnprob, True
    return result.lnprob, False


def _trace_steps(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    length: int,
    *,
    raw: bool,
) -> Array:
    if raw:
        return np.arange(length)
    return result.discarded_steps + np.arange(length) * settings.thin


def _mark_burn_in(ax: plt.Axes, result: McmcResult, *, raw: bool) -> None:
    if raw and result.discarded_steps > 0:
        ax.axvspan(
            0,
            result.discarded_steps,
            color="0.9",
            label="discarded burn-in",
            zorder=-1,
        )
        ax.axvline(result.discarded_steps, color="0.5", linestyle="--", linewidth=1.0)


def _save_trace_pages(
    pdf: PdfPages,
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    parameter_names: tuple[str, ...],
) -> None:
    chain, raw = _trace_chain(result)
    steps = _trace_steps(result, settings, chain.shape[0], raw=raw)
    for index, name in enumerate(parameter_names):
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        _mark_burn_in(ax, result, raw=raw)
        for walker in range(chain.shape[1]):
            ax.plot(steps, chain[:, walker, index], linewidth=0.5, alpha=0.35)
        ax.set_title(f"MCMC trace: {name}")
        ax.set_xlabel("Step")
        ax.set_ylabel(name)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles[:1], labels[:1], loc="best")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def _save_log_probability_page(
    pdf: PdfPages,
    result: McmcResult,
    settings: EffectiveMcmcSettings,
) -> None:
    lnprob, raw = _trace_lnprob(result)
    steps = _trace_steps(result, settings, lnprob.shape[0], raw=raw)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    _mark_burn_in(ax, result, raw=raw)
    for walker in range(lnprob.shape[1]):
        ax.plot(steps, lnprob[:, walker], linewidth=0.5, alpha=0.35)
    ax.set_title("MCMC log-probability trace")
    ax.set_xlabel("Step")
    ax.set_ylabel("lnprob")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[:1], labels[:1], loc="best")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _autocorrelation_monitor_values(chain: Array) -> tuple[Array, Array, Array]:
    nsteps = chain.shape[0]
    if nsteps < 4:
        empty = np.empty(0, dtype=float)
        return empty, empty, empty

    start = max(4, min(100, nsteps // 10))
    lengths = np.unique(np.geomspace(start, nsteps, num=12).astype(int))
    mean_tau: list[float] = []
    max_tau: list[float] = []
    used_lengths: list[int] = []
    for length in lengths:
        try:
            tau = integrated_time(chain[:length], quiet=False)
        except AutocorrError as error:
            tau = error.tau
        except (FloatingPointError, ValueError):
            continue
        tau = np.asarray(tau, dtype=float)
        if tau.ndim != 1 or not np.all(np.isfinite(tau)) or np.any(tau <= 0.0):
            continue
        used_lengths.append(int(length))
        mean_tau.append(float(np.mean(tau)))
        max_tau.append(float(np.max(tau)))
    return (
        np.asarray(used_lengths, dtype=float),
        np.asarray(mean_tau, dtype=float),
        np.asarray(max_tau, dtype=float),
    )


def _save_autocorrelation_monitor_page(pdf: PdfPages, result: McmcResult) -> None:
    chain, _raw = _trace_chain(result)
    lengths, mean_tau, max_tau = _autocorrelation_monitor_values(chain)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    if len(lengths):
        ax.plot(lengths, mean_tau, "o-", label="mean tau")
        ax.plot(lengths, max_tau, "o-", label="max tau")
        ax.plot(lengths, lengths / 50.0, "--", color="0.4", label="N / 50")
        ax.plot(lengths, lengths / 100.0, ":", color="0.4", label="N / 100")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "Autocorrelation unavailable", ha="center", va="center")
    ax.set_title("MCMC autocorrelation monitor")
    ax.set_xlabel("Steps")
    ax.set_ylabel("Integrated autocorrelation time")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _save_mcmc_correlation_pages(
    pdf: PdfPages,
    result: McmcResult,
    parameter_names: tuple[str, ...],
    correlated_pairs: list[tuple[int, int, float]],
) -> None:
    samples = result.samples
    for pair in correlated_pairs:
        _save_pair_distribution_page(
            pdf,
            samples=samples,
            parameter_names=parameter_names,
            pair=pair,
            title_prefix="Posterior pair",
            cmap="Blues",
        )


def write_mcmc_plots(
    result: McmcResult,
    settings: EffectiveMcmcSettings,
    path: Path,
    *,
    parameter_names: tuple[str, ...],
) -> None:
    correlated_pairs = _correlated_pairs(
        result.correlations,
        threshold=_CORRELATION_THRESHOLD,
    )
    with PdfPages(str(path / "plots.pdf")) as pdf:
        _save_mcmc_summary_page(
            pdf,
            result,
            settings,
            parameter_names,
            correlated_pairs,
        )
        _save_mcmc_distribution_pages(pdf, result, parameter_names)
        _save_trace_pages(pdf, result, settings, parameter_names)
        _save_log_probability_page(pdf, result, settings)
        _save_autocorrelation_monitor_page(pdf, result)
        _save_mcmc_correlation_pages(pdf, result, parameter_names, correlated_pairs)


def _save_resampling_summary_page(
    pdf: PdfPages,
    *,
    method: str,
    fitmethod: str,
    parameter_names: tuple[str, ...],
    correlations: Array,
    requested_samples: int,
    completed_samples: int,
    chisqr_values: Array,
) -> list[tuple[int, int, float]]:
    correlated_pairs = _correlated_pairs(
        correlations,
        threshold=_CORRELATION_THRESHOLD,
    )
    finite_chisqr = _finite_values(chisqr_values)
    lines = [
        f"Method: {method}",
        f"Fit method: {fitmethod}",
        f"Parameters: {len(parameter_names)}",
        f"Requested samples: {requested_samples}",
        f"Completed samples: {completed_samples}",
    ]
    if len(finite_chisqr):
        lines.extend(
            [
                f"chisqr mean: {_format_float(float(np.mean(finite_chisqr)))}",
                f"chisqr median: {_format_float(float(np.median(finite_chisqr)))}",
                f"chisqr min: {_format_float(float(np.min(finite_chisqr)))}",
                f"chisqr max: {_format_float(float(np.max(finite_chisqr)))}",
            ],
        )
    else:
        lines.append("chisqr: no finite values")

    lines.extend(["", f"Correlation threshold: |r| >= {_CORRELATION_THRESHOLD:g}"])
    if correlated_pairs:
        for index_x, index_y, correlation in correlated_pairs:
            lines.append(
                f"{parameter_names[index_x]} vs {parameter_names[index_y]}: "
                f"r={correlation:.3f}",
            )
    else:
        lines.append("No parameter pairs exceeded the correlation threshold.")

    _save_text_page(pdf, f"{method} Diagnostics", lines)
    return correlated_pairs


def _save_resampling_distribution_pages(
    pdf: PdfPages,
    *,
    samples: Array,
    parameter_names: tuple[str, ...],
) -> None:
    for index, name in enumerate(parameter_names):
        values = samples[:, index] if samples.size else np.empty(0, dtype=float)
        finite_values = _finite_values(values)
        interval = None
        median = None
        if len(finite_values):
            lower_95, median, upper_95 = np.percentile(
                finite_values,
                [2.5, 50.0, 97.5],
            )
            interval = (float(lower_95), float(upper_95))
            median = float(median)
        _save_histogram_page(
            pdf,
            values=values,
            title=f"Sample distribution: {name}",
            xlabel=name,
            ylabel="Density",
            color="#55A868",
            density=True,
            interval=interval,
            interval_label="95% percentile",
            median=median,
        )


def _save_chisqr_distribution_page(pdf: PdfPages, chisqr_values: Array) -> None:
    values = _finite_values(chisqr_values)
    median = float(np.median(values)) if len(values) else None
    _save_histogram_page(
        pdf,
        values=values,
        title="chisqr distribution",
        xlabel="chisqr",
        ylabel="Samples",
        color="#8172B3",
        density=False,
        median=median,
        empty_message="No finite chisqr values",
    )


def _save_resampling_correlation_pages(
    pdf: PdfPages,
    *,
    samples: Array,
    parameter_names: tuple[str, ...],
    correlated_pairs: list[tuple[int, int, float]],
) -> None:
    for pair in correlated_pairs:
        _save_pair_distribution_page(
            pdf,
            samples=samples,
            parameter_names=parameter_names,
            pair=pair,
            title_prefix="Sample pair",
            cmap="Greens",
        )


def write_resampling_plots(
    path: Path,
    *,
    method: str,
    fitmethod: str,
    parameter_names: tuple[str, ...],
    samples: Array,
    chisqr_values: Array,
    correlations: Array,
    requested_samples: int,
    completed_samples: int,
) -> None:
    with PdfPages(str(path / "plots.pdf")) as pdf:
        correlated_pairs = _save_resampling_summary_page(
            pdf,
            method=method,
            fitmethod=fitmethod,
            parameter_names=parameter_names,
            correlations=correlations,
            requested_samples=requested_samples,
            completed_samples=completed_samples,
            chisqr_values=chisqr_values,
        )
        _save_resampling_distribution_pages(
            pdf,
            samples=samples,
            parameter_names=parameter_names,
        )
        _save_chisqr_distribution_page(pdf, chisqr_values)
        _save_resampling_correlation_pages(
            pdf,
            samples=samples,
            parameter_names=parameter_names,
            correlated_pairs=correlated_pairs,
        )
