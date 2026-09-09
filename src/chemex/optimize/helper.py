from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np
from scipy import stats

from chemex.containers.experiments import Experiments
from chemex.exceptions import ArtifactPublicationError
from chemex.messages import (
    print_making_plots,
    print_writing_results,
)
from chemex.optimize.deterministic_uncertainty import (
    DerivationDisposition,
    DeterministicUncertainty,
    InterpretationCompleteness,
)
from chemex.optimize.statistics import FitStatisticsCounts
from chemex.optimize.uncertainty import UncertaintyUnavailableKind
from chemex.parameters.feasible_coordinates import validate_relaxation_state
from chemex.parameters.parameterization import (
    ActiveParameterization,
    SealedParameterModel,
)
from chemex.printers.native_reporting import (
    write_block_recovery,
    write_block_uncertainty,
    write_json,
    write_uncertainty,
)
from chemex.printers.parameters import (
    GRID_UNCERTAINTY_WITHHELD_TEXT,
    uncertainty_unavailable_reason,
    write_parameters,
)
from chemex.typing import Array


def calculate_statistics_from_residuals(
    residuals: Array,
    controlled_coordinate_count: int,
    profiled_normalization_count: int = 0,
) -> dict[str, int | float]:
    """Calculate established fit statistics from an authoritative residual vector."""
    counts = FitStatisticsCounts(
        len(residuals),
        controlled_coordinate_count,
        profiled_normalization_count,
    )
    chisqr = sum(residuals**2)
    degrees_of_freedom = counts.residual_degrees_of_freedom
    estimated_parameter_count = counts.estimated_parameter_count
    redchi = counts.reduced_chi_square(chisqr)
    aic = chisqr + 2 * estimated_parameter_count
    bic = chisqr + np.log(counts.residual_count) * estimated_parameter_count
    _, ks_p_value = stats.kstest(residuals, "norm")
    ks_p_value = float(ks_p_value)
    positive_degrees_of_freedom = counts.positive_residual_degrees_of_freedom
    pvalue = (
        float("nan")
        if positive_degrees_of_freedom is None
        else float(1.0 - stats.chi2.cdf(chisqr, positive_degrees_of_freedom))
    )
    return {
        "ndata": counts.residual_count,
        "nvarys": counts.controlled_coordinate_count,
        "dof": degrees_of_freedom,
        "chisqr": chisqr,
        "redchi": redchi,
        "pvalue": pvalue,
        "ks_pvalue": ks_p_value,
        "aic": aic,
        "bic": bic,
    }


def _write_statistics(
    experiments: Experiments,
    path: Path,
    *,
    residuals: Array,
    controlled_coordinate_count: int,
    profiled_normalization_count: int,
) -> None:
    """Write fitting statistics to a file."""
    stats = calculate_statistics_from_residuals(
        residuals,
        controlled_coordinate_count,
        profiled_normalization_count,
    )
    filename = path / "statistics.toml"
    with filename.open("w", encoding="utf-8") as f:
        f.write(f'"number of data points"                = {stats["ndata"]}\n')
        f.write(f'"number of variables"                  = {stats["nvarys"]}\n')
        f.write(f'"chi-square"                           = {stats["chisqr"]: .5e}\n')
        f.write(f'"reduced-chi-square"                   = {stats["redchi"]: .5e}\n')
        f.write(f'"chi-squared test"                     = {stats["pvalue"]: .5e}\n')
        f.write(
            f'"Kolmogorov-Smirnov test"              = {stats["ks_pvalue"]: .5e}\n',
        )
        f.write(f'"Akaike Information Criterion (AIC)"   = {stats["aic"]: .5e}\n')
        f.write(f'"Bayesian Information Criterion (BIC)" = {stats["bic"]: .5e}\n')


def _write_files(
    experiments: Experiments,
    path: Path,
    *,
    residuals: Array,
    controlled_coordinate_count: int,
    profiled_normalization_count: int,
    deterministic_uncertainty: DeterministicUncertainty | None = None,
    parameter_model: SealedParameterModel,
    parameter_values: Mapping[str, float],
    parameterization: ActiveParameterization,
    fitted_ids: tuple[str, ...] = (),
    report_only_values: Mapping[str, float] | None = None,
) -> None:
    """Write the results of the fit to output files."""
    print_writing_results(path)
    path.mkdir(parents=True, exist_ok=True)
    write_parameters(
        path,
        parameter_model=parameter_model,
        parameter_values=parameter_values,
        parameterization=parameterization,
        fitted_ids=fitted_ids,
        deterministic_uncertainty=deterministic_uncertainty,
        report_only_values=report_only_values,
    )
    experiments.write(path)
    uncertainty_interrupted = (
        deterministic_uncertainty is not None
        and deterministic_uncertainty.derivation_interrupted
    )
    if not uncertainty_interrupted:
        _write_statistics(
            experiments,
            path=path,
            residuals=residuals,
            controlled_coordinate_count=controlled_coordinate_count,
            profiled_normalization_count=profiled_normalization_count,
        )
    if (
        deterministic_uncertainty is not None
        and deterministic_uncertainty.root_evidence is not None
    ):
        statistics_path = path / "Statistics"
        statistics_path.mkdir(exist_ok=True)
        write_uncertainty(statistics_path, deterministic_uncertainty.root_evidence)
        write_block_uncertainty(
            statistics_path,
            deterministic_uncertainty.block_evidence,
        )
        write_block_recovery(
            statistics_path,
            deterministic_uncertainty.block_operation,
        )
    if deterministic_uncertainty is not None and (
        deterministic_uncertainty.disposition is DerivationDisposition.WITHHELD
        or deterministic_uncertainty.completeness
        is InterpretationCompleteness.INCOMPLETE
    ):
        if deterministic_uncertainty.disposition is DerivationDisposition.WITHHELD:
            terminal = "withheld"
            reason = GRID_UNCERTAINTY_WITHHELD_TEXT
        else:
            incomplete_terminal = deterministic_uncertainty.incomplete_terminal
            if incomplete_terminal is None:
                raise ValueError("Incomplete uncertainty lacks its terminal")
            terminal = incomplete_terminal.value
            reason = uncertainty_unavailable_reason(
                UncertaintyUnavailableKind.DERIVATION_STOPPED
            )
        covariance_path = path / "Statistics" / "Covariance"
        covariance_path.mkdir(parents=True, exist_ok=True)
        write_json(
            covariance_path / "status.json",
            {
                "artifact_type": "native_covariance_derivation_status",
                "schema_version": 1,
                "status": "incomplete",
                "terminal": terminal,
                "reason": reason,
            },
        )


def _write_simulation_files(
    experiments: Experiments,
    path: Path,
    *,
    parameter_model: SealedParameterModel,
    parameter_values: Mapping[str, float],
    parameterization: ActiveParameterization,
    report_only_values: Mapping[str, float] | None = None,
) -> None:
    """Write the results of the simulation to output files."""
    print_writing_results(path)
    path.mkdir(parents=True, exist_ok=True)
    write_parameters(
        path,
        parameter_model=parameter_model,
        parameter_values=parameter_values,
        parameterization=parameterization,
        report_only_values=report_only_values,
    )
    experiments.write(path)


def _write_plots(experiments: Experiments, path: Path) -> None:
    """Plot the experimental and fitted data."""
    print_making_plots()

    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    experiments.plot(path=path_)


def _write_simulation_plots(experiments: Experiments, path: Path) -> None:
    """Plot the experimental and fitted data."""
    print_making_plots()

    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    experiments.plot_simulation(path=path_)


def execute_post_fit(
    experiments: Experiments,
    path: Path,
    *,
    plot: bool = False,
    residuals: Array,
    controlled_coordinate_count: int,
    profiled_normalization_count: int,
    deterministic_uncertainty: DeterministicUncertainty | None = None,
    parameter_model: SealedParameterModel,
    parameter_values: Mapping[str, float],
    parameterization: ActiveParameterization,
    fitted_ids: tuple[str, ...] = (),
    report_only_values: Mapping[str, float] | None = None,
) -> None:
    _write_files(
        experiments,
        path,
        residuals=residuals,
        controlled_coordinate_count=controlled_coordinate_count,
        profiled_normalization_count=profiled_normalization_count,
        deterministic_uncertainty=deterministic_uncertainty,
        parameter_model=parameter_model,
        parameter_values=parameter_values,
        parameterization=parameterization,
        fitted_ids=fitted_ids,
        report_only_values=report_only_values,
    )
    uncertainty_interrupted = (
        deterministic_uncertainty is not None
        and deterministic_uncertainty.derivation_interrupted
    )
    if plot and not uncertainty_interrupted:
        _write_plots(experiments, path)


def execute_simulation(
    experiments: Experiments,
    path: Path,
    *,
    parameter_values: Mapping[str, float],
    parameter_model: SealedParameterModel,
    parameterization: ActiveParameterization,
    plot: bool = False,
    report_only_values: Mapping[str, float] | None = None,
) -> None:
    validate_relaxation_state(parameterization, parameter_values)
    experiments.prepare_for_simulation(parameter_values)
    try:
        _write_simulation_files(
            experiments,
            path,
            parameter_model=parameter_model,
            parameter_values=parameter_values,
            parameterization=parameterization,
            report_only_values=report_only_values,
        )
        if plot:
            _write_simulation_plots(experiments, path)
    except OSError as error:
        filename = error.filename
        destination = Path(filename) if isinstance(filename, str) else path
        failure = ArtifactPublicationError(
            "publish simulation output",
            destination,
            error,
        )
        object.__setattr__(failure, "failure_stage", "output")
        raise failure from error
    except (Exception, KeyboardInterrupt) as error:
        if getattr(error, "failure_stage", None) is None:
            object.__setattr__(error, "failure_stage", "output")
        raise


def print_header(
    grid: Iterable[str],
    *,
    parameter_names: Mapping[str, object],
) -> str:
    header_pnames = " ".join(f"{parameter_names[param_id]}" for param_id in grid)
    return f"# {header_pnames} [χ²]\n"


def print_values(values: Iterable[float], chisqr: float) -> str:
    body_values = " ".join(f"{value:.5e}" for value in values)
    return f"  {body_values} {chisqr:.5e}\n"
