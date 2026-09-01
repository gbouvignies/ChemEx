from __future__ import annotations

import secrets
from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Literal, cast

import numpy as np

from chemex.atomic import open_text_atomic, remove_paths_best_effort, write_text_atomic
from chemex.configuration.method_plan import McmcRequest
from chemex.configuration.methods import McmcSettings
from chemex.containers.experiments import Experiments
from chemex.optimize.native_deterministic import NativeDeterministicFit
from chemex.optimize.native_mcmc import (
    EnsembleState,
    McmcAcceptanceSummary,
    McmcAnalysisFailureCategory,
    McmcAnalysisResult,
    McmcAnalysisStatus,
    McmcAutocorrelationReport,
    McmcEvidence,
    McmcOperation,
    McmcOperationTerminal,
    McmcPlan,
    RawMcmcCapture,
    derive_mcmc_analysis_result,
    execute_mcmc_evidence,
    product_mcmc_invalid_bound_ids,
    resolve_product_mcmc_policy,
)
from chemex.optimize.statistics_plot import write_mcmc_plots
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.runtime import ExecutionSettings
from chemex.typing import Array

type McmcBurnSetting = int | Literal["auto"]


@dataclass(frozen=True)
class EffectiveMcmcSettings:
    steps: int
    burn: int | Literal["auto"]
    thin: int
    walkers: int
    seed: int | None
    workers: int
    native_threads: int | None
    update_parameters: bool


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

_MCMC_ARTIFACT_NAMES = (
    "summary.toml",
    "samples.tsv",
    "correlations.tsv",
    "diagnostics.toml",
    "plots.pdf",
    "raw_chain.tsv",
    "raw_capture.tsv",
)


class NativeMcmcIncompleteError(RuntimeError):
    """Raised when native product MCMC cannot publish a complete chain."""

    def __init__(
        self,
        message: str,
        *,
        terminal: str = "failed",
        completed_steps: int = 0,
        preserve_raw_evidence: bool = False,
        autocorrelation_time: Array | None = None,
        autocorrelation_time_reliable: bool = False,
        autocorrelation_warning: str | None = None,
    ) -> None:
        super().__init__(message)
        self.terminal = terminal
        self.completed_steps = completed_steps
        self.preserve_raw_evidence = preserve_raw_evidence
        self.autocorrelation_time = autocorrelation_time
        self.autocorrelation_time_reliable = autocorrelation_time_reliable
        self.autocorrelation_warning = autocorrelation_warning


def _compatibility_error_from_result(
    result: McmcAnalysisResult,
) -> NativeMcmcIncompleteError:
    failure = result.failure
    if failure is None:
        raise RuntimeError("Incomplete MCMC Analysis Result has no typed failure")
    report = result.autocorrelation_report
    autocorrelation_time = (
        None
        if report is None or report.values is None
        else np.asarray(report.values, dtype=float)
    )
    autocorrelation_warning = (
        None
        if failure.category
        is McmcAnalysisFailureCategory.ACCEPTANCE_DIAGNOSTICS_UNAVAILABLE
        or report is None
        else report.warning
    )
    completed_steps = result.evidence.completed_transition_count
    return NativeMcmcIncompleteError(
        failure.message,
        completed_steps=completed_steps,
        preserve_raw_evidence=failure.preserve_raw_evidence,
        autocorrelation_time=autocorrelation_time,
        autocorrelation_time_reliable=report is not None and report.reliable,
        autocorrelation_warning=autocorrelation_warning,
    )


def resolve_mcmc_request(
    request: McmcRequest,
    *,
    nvarys: int,
    execution: ExecutionSettings | None = None,
) -> EffectiveMcmcSettings:
    walkers = request.walkers or max(32, 2 * nvarys)
    min_walkers = 2 * nvarys
    if walkers < min_walkers:
        msg = (
            f"MCMC requires at least {min_walkers} walkers for {nvarys} fitted "
            f"parameters"
        )
        raise ValueError(msg)

    workers = request.workers or (execution.workers if execution is not None else 1)
    native_threads = execution.native_threads if execution is not None else None

    return EffectiveMcmcSettings(
        steps=request.steps,
        burn="auto" if request.burn is None else request.burn,
        thin=request.thin,
        walkers=walkers,
        seed=request.seed,
        workers=workers,
        native_threads=native_threads,
        update_parameters=False,
    )


def resolve_mcmc_settings(
    settings: McmcSettings,
    *,
    nvarys: int,
    execution: ExecutionSettings | None = None,
) -> EffectiveMcmcSettings:
    """Preserve the standalone legacy settings resolver compatibility surface."""
    effective = resolve_mcmc_request(
        McmcRequest(
            steps=settings.steps,
            burn=None if settings.burn == "auto" else settings.burn,
            seed=settings.seed,
            thin=settings.thin,
            walkers=settings.walkers,
            workers=settings.workers,
        ),
        nvarys=nvarys,
        execution=execution,
    )
    return replace(effective, update_parameters=settings.update_parameters)


def _format_parameter_ids(
    parameter_ids: tuple[str, ...],
    parameter_model: SealedParameterModel,
) -> tuple[str, ...]:
    return tuple(
        str(parameter_name_from_definition(parameter_model.definitions[param_id]))
        for param_id in parameter_ids
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


def _extend_autocorrelation_diagnostics(
    lines: list[str],
    report: McmcAutocorrelationReport,
) -> None:
    if report.values is None:
        warning = report.warning or "not available"
        lines.append(f"autocorrelation_warning = {_quote_toml_string(warning)}")
        return

    suffix = "" if report.reliable else "_tentative"
    maximum = cast("float", report.maximum)
    sampled_steps_over_maximum = cast(
        "float",
        report.sampled_steps_over_maximum,
    )
    tau_key = f"autocorrelation_time{suffix}"
    max_tau_key = f"max_autocorrelation_time{suffix}"
    lines.extend(
        [
            f"{tau_key} = {_format_toml_float_list(list(report.values))}",
            f"{max_tau_key} = {_format_toml_float(maximum)}",
            (
                "steps_over_max_autocorrelation_time = "
                f"{_format_toml_float(sampled_steps_over_maximum)}"
            ),
            f"recommended_min_steps_50tau = {report.recommended_min_steps_50tau}",
            f"recommended_min_steps_100tau = {report.recommended_min_steps_100tau}",
        ],
    )
    if report.retained_steps_over_maximum is not None:
        lines.append(
            "retained_steps_over_max_autocorrelation_time = "
            f"{_format_toml_float(report.retained_steps_over_maximum)}"
        )
    if report.minimum_effective_sample_size is not None:
        lines.append(
            "min_effective_sample_size = "
            f"{_format_toml_float(report.minimum_effective_sample_size)}",
        )
    if report.warning is not None:
        lines.append(f"autocorrelation_warning = {_quote_toml_string(report.warning)}")
    if report.effective_sample_size_warning is not None:
        lines.append(
            "effective_sample_size_warning = "
            f"{_quote_toml_string(report.effective_sample_size_warning)}"
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
    result: McmcAnalysisResult,
    path: Path,
    parameter_model: SealedParameterModel,
) -> None:
    lines: list[str] = []
    for summary in result.summary:
        param_id = summary.parameter_id
        parameter_name = str(
            parameter_name_from_definition(parameter_model.definitions[param_id])
        ).strip("[]")
        if summary.prior_lower is None or summary.prior_upper is None:
            raise RuntimeError(
                "Authoritative MCMC Summary has no resolved prior bounds"
            )
        lines.extend(
            [
                f"[{_quote_toml_string(parameter_name)}]",
                'prior = "uniform"',
                f"prior_lower = {_format_toml_float(summary.prior_lower)}",
                f"prior_upper = {_format_toml_float(summary.prior_upper)}",
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
    result: McmcAnalysisResult,
    path: Path,
    parameter_model: SealedParameterModel,
) -> None:
    parameter_names = _format_parameter_ids(result.var_names, parameter_model)
    values = np.column_stack((result.samples, result.log_probabilities))

    with (path / "samples.tsv").open("w", encoding="utf-8") as fileout:
        fileout.write("\t".join((*parameter_names, "lnprob")) + "\n")
        np.savetxt(fileout, values, fmt="%.5e", delimiter="\t")


def _write_raw_chain_evidence(
    evidence: McmcEvidence,
    path: Path,
    parameter_model: SealedParameterModel,
) -> None:
    parameter_names = _format_parameter_ids(evidence.coordinate_ids, parameter_model)
    _write_raw_states(path / "raw_chain.tsv", evidence.states[1:], parameter_names)


def _write_raw_states(
    destination: Path,
    states: tuple[EnsembleState, ...],
    parameter_names: tuple[str, ...],
) -> None:
    """Write already-selected complete ensemble states atomically."""
    with open_text_atomic(destination) as fileout:
        fileout.write("\t".join(("step", "walker", *parameter_names, "lnprob")) + "\n")
        for state in states:
            for walker, (position, lnprob) in enumerate(
                zip(state.positions, state.log_densities, strict=True),
            ):
                values = "\t".join(_format_toml_float(value) for value in position)
                fileout.write(
                    f"{state.ordinal}\t{walker}\t{values}\t"
                    f"{_format_toml_float(lnprob)}\n"
                )


def _write_raw_capture(
    capture: RawMcmcCapture,
    path: Path,
    parameter_model: SealedParameterModel,
    parameter_ids: tuple[str, ...],
) -> None:
    parameter_names = _format_parameter_ids(parameter_ids, parameter_model)
    _write_raw_states(path / "raw_capture.tsv", capture.states, parameter_names)


def _extend_acceptance_diagnostics(
    lines: list[str],
    acceptance: McmcAcceptanceSummary,
) -> None:
    lines.extend(
        (
            (f"acceptance_fraction_mean = {_format_toml_float(acceptance.mean)}"),
            (f"acceptance_fraction_min = {_format_toml_float(acceptance.minimum)}"),
            (f"acceptance_fraction_max = {_format_toml_float(acceptance.maximum)}"),
        )
    )


def _write_correlations(
    result: McmcAnalysisResult,
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
    result: McmcAnalysisResult,
    settings: EffectiveMcmcSettings,
    path: Path,
    timings: Mapping[str, float],
    *,
    engine: str,
    root_seed: int | None = None,
) -> None:
    acceptance = result.acceptance_summary
    autocorrelation = result.autocorrelation_report
    if acceptance is None or autocorrelation is None:
        raise RuntimeError(
            "Authoritative MCMC Analysis Result has no publication diagnostics"
        )
    requested_burn = '"auto"' if settings.burn == "auto" else str(settings.burn)
    lines = [
        'status = "complete"',
        f"engine = {_quote_toml_string(engine)}",
        'sampler = "emcee via ChemEx direct EnsembleSampler"',
        f"emcee_version = {_quote_toml_string(_package_version('emcee'))}",
        'credible_interval = "95% equal-tailed"',
        'convergence_diagnostic = "integrated_autocorrelation_time"',
        f"autocorrelation_status = {_quote_toml_string(autocorrelation.status)}",
        'rhat = "not computed: emcee ensemble walkers are not independent chains"',
        f"steps = {settings.steps}",
        f"requested_burn = {requested_burn}",
        f"discarded_steps = {result.discarded_steps}",
        f"thin = {settings.thin}",
        f"walkers = {settings.walkers}",
        f"workers = {settings.workers}",
        f"retained_steps = {result.retained_step_count}",
        f"retained_samples = {result.retained_sample_count}",
        'samples_file = "samples.tsv"',
        'summary_file = "summary.toml"',
        'correlations_file = "correlations.tsv"',
        'plots_file = "plots.pdf"',
    ]
    _extend_acceptance_diagnostics(lines, acceptance)
    lines.append("unbounded_parameters = []")
    if root_seed is not None:
        lines.append(f"root_seed = {root_seed}")
    _extend_autocorrelation_diagnostics(lines, autocorrelation)
    _extend_timing_diagnostics(lines, timings)
    write_text_atomic(path / "diagnostics.toml", "\n".join(lines) + "\n")


def write_mcmc_outputs(
    result: McmcAnalysisResult,
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
    for name in _MCMC_ARTIFACT_NAMES:
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
    raw_evidence: McmcEvidence | None = None,
    raw_capture: RawMcmcCapture | None = None,
    analysis_result: McmcAnalysisResult | None = None,
    timings: Mapping[str, float] | None = None,
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
    if raw_evidence is not None:
        raw_steps = (
            raw_evidence.completed_transition_count
            if analysis_result is None
            else analysis_result.raw_step_count
        )
        raw_samples = (
            raw_steps * raw_evidence.plan.policy.walkers
            if analysis_result is None
            else analysis_result.raw_sample_count
        )
        lines.extend(
            (
                (
                    "sampling_terminal = "
                    f"{_quote_toml_string(raw_evidence.terminal.value)}"
                ),
                'posterior_summary = "withheld"',
                'raw_evidence = "qualified_chain"',
                'raw_chain_file = "raw_chain.tsv"',
                f"raw_steps = {raw_steps}",
                f"raw_samples = {raw_samples}",
                f"emcee_version = {_quote_toml_string(_package_version('emcee'))}",
            )
        )
        if terminal == "interrupted":
            lines.extend(
                (
                    'posterior_retention = "withheld"',
                    'burn_validity = "withheld"',
                    'convergence = "withheld"',
                    'correlations = "withheld"',
                )
            )
        else:
            lines.extend(
                (
                    'posterior_retention = "failed"',
                    'convergence_diagnostic = "integrated_autocorrelation_time"',
                    'rhat = "not computed: emcee ensemble walkers are not independent chains"',
                )
            )
        if (
            terminal != "interrupted"
            and analysis_result is not None
            and analysis_result.acceptance_summary is not None
        ):
            _extend_acceptance_diagnostics(lines, analysis_result.acceptance_summary)
        if (
            terminal != "interrupted"
            and analysis_result is not None
            and analysis_result.autocorrelation_report is not None
        ):
            lines.append(
                "autocorrelation_status = "
                f"{_quote_toml_string(analysis_result.autocorrelation_report.status)}"
            )
            _extend_autocorrelation_diagnostics(
                lines,
                analysis_result.autocorrelation_report,
            )
    elif raw_capture is not None:
        lines.extend(
            (
                (
                    "sampling_terminal = "
                    f"{_quote_toml_string(raw_capture.terminal.value)}"
                ),
                'posterior_retention = "withheld"',
                'posterior_summary = "withheld"',
                'burn_validity = "withheld"',
                'convergence = "withheld"',
                'correlations = "withheld"',
                'acceptance_fraction = "withheld"',
                'raw_evidence = "execution_capture_unqualified"',
                'raw_capture_file = "raw_capture.tsv"',
                f"captured_states = {raw_capture.complete_state_count}",
                (
                    "completed_transitions = "
                    f"{max(0, raw_capture.complete_state_count - 1)}"
                ),
                f"emcee_version = {_quote_toml_string(_package_version('emcee'))}",
            )
        )
    if timings is not None:
        _extend_timing_diagnostics(lines, timings)
    write_text_atomic(path / "diagnostics.toml", "\n".join(lines) + "\n")


def run_native_mcmc(  # noqa: C901 - closed execution/publication lifecycle
    experiments: Experiments,
    fit: NativeDeterministicFit,
    request: McmcRequest,
    path: Path,
    *,
    execution: ExecutionSettings | None = None,
    seed_recorder: Callable[[int], None] | None = None,
) -> McmcAnalysisResult:
    """Run product MCMC directly from one committed native deterministic fit."""
    effective = resolve_mcmc_request(
        request,
        nvarys=len(fit.problem.controlled_ids),
        execution=execution,
    )
    root_seed = secrets.randbits(64) if effective.seed is None else effective.seed
    if seed_recorder is not None:
        seed_recorder(root_seed)
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
    invalid_bound_ids = product_mcmc_invalid_bound_ids(fit.problem)
    if invalid_bound_ids:
        parameter_names = _format_parameter_ids(
            invalid_bound_ids,
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
    mcmc_operation: McmcOperation | None = None
    mcmc_evidence: McmcEvidence | None = None
    analysis_result: McmcAnalysisResult | None = None
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
        mcmc_operation = execute_mcmc_evidence(
            fit.accepted,
            plan,
            execution=ExecutionSettings(
                workers=effective.workers,
                native_threads=effective.native_threads,
            ),
        )
        mcmc_evidence = mcmc_operation.evidence
        timings["sampling_seconds"] = perf_counter() - phase_start
        completed_steps = max(0, mcmc_operation.complete_state_count - 1)
        if (
            mcmc_operation.terminal is not McmcOperationTerminal.COMPLETED
            or mcmc_operation.evidence is None
        ):
            error = NativeMcmcIncompleteError(
                "Native MCMC "
                f"{mcmc_operation.terminal.value}: "
                f"{mcmc_operation.failure_message}",
                terminal=mcmc_operation.terminal.value,
                completed_steps=completed_steps,
                preserve_raw_evidence=(
                    mcmc_operation.terminal is McmcOperationTerminal.INTERRUPTED
                    and mcmc_operation.evidence is not None
                ),
            )
            raise error
        phase_start = perf_counter()
        try:
            analysis_result = derive_mcmc_analysis_result(
                mcmc_operation.evidence,
                request,
                fit.parameter_model,
            )
        finally:
            timings["result_processing_seconds"] = perf_counter() - phase_start
        if analysis_result.status is McmcAnalysisStatus.INCOMPLETE:
            raise _compatibility_error_from_result(analysis_result)
        write_mcmc_outputs(
            analysis_result,
            effective,
            path,
            fit.parameter_model,
            timings=timings,
            engine="native MCMC",
            root_seed=root_seed,
        )
    except (Exception, KeyboardInterrupt) as error:
        raw_evidence = (
            analysis_result.evidence
            if analysis_result is not None
            and analysis_result.status is McmcAnalysisStatus.INCOMPLETE
            and analysis_result.failure is not None
            and analysis_result.failure.preserve_raw_evidence
            else mcmc_evidence
            if (
                isinstance(error, KeyboardInterrupt)
                or (
                    isinstance(error, NativeMcmcIncompleteError)
                    and error.preserve_raw_evidence
                )
            )
            else None
        )
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
        raw_capture = (
            mcmc_operation.raw_capture
            if terminal == "interrupted"
            and raw_evidence is None
            and mcmc_operation is not None
            and mcmc_operation.raw_capture is not None
            and mcmc_operation.raw_capture.complete_state_count > 0
            else None
        )
        if terminal == "interrupted":
            remove_paths_best_effort(
                (statistic_path / name for name in _MCMC_ARTIFACT_NAMES),
                error,
            )
            published_evidence: McmcEvidence | None = None
            published_capture: RawMcmcCapture | None = None
            if raw_evidence is not None:
                try:
                    _write_raw_chain_evidence(
                        raw_evidence,
                        statistic_path,
                        fit.parameter_model,
                    )
                    published_evidence = raw_evidence
                except (Exception, KeyboardInterrupt):  # noqa: BLE001
                    error.add_note(
                        "ChemEx could not publish an interrupted qualified MCMC chain."
                    )
            elif raw_capture is not None:
                try:
                    _write_raw_capture(
                        raw_capture,
                        statistic_path,
                        fit.parameter_model,
                        fit.problem.controlled_ids,
                    )
                    published_capture = raw_capture
                except (Exception, KeyboardInterrupt):  # noqa: BLE001
                    error.add_note(
                        "ChemEx could not publish interrupted MCMC execution capture."
                    )
            try:
                _write_native_mcmc_state_diagnostics(
                    statistic_path,
                    settings=effective,
                    root_seed=root_seed,
                    parameter_ids=fit.problem.controlled_ids,
                    status="incomplete",
                    terminal=terminal,
                    completed_steps=failure_completed_steps,
                    failure=error,
                    raw_evidence=published_evidence,
                    raw_capture=published_capture,
                    analysis_result=analysis_result,
                    timings=timings,
                )
            except (Exception, KeyboardInterrupt):  # noqa: BLE001
                error.add_note("ChemEx could not publish interrupted MCMC diagnostics.")
            raise
        _clear_mcmc_artifacts(statistic_path)
        if raw_evidence is not None:
            _write_raw_chain_evidence(
                raw_evidence,
                statistic_path,
                fit.parameter_model,
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
            raw_evidence=raw_evidence,
            raw_capture=None,
            analysis_result=analysis_result,
            timings=timings,
        )
        raise
    else:
        return analysis_result
