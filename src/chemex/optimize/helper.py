from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np
from scipy import stats

from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_chi2,
    print_group_name,
    print_making_plots,
    print_plotting_canceled,
    print_writing_results,
)
from chemex.parameters.database import ParameterStore
from chemex.printers.parameters import write_parameters
from chemex.typing import Array


def calculate_statistics_from_residuals(
    residuals: Array,
    nvarys: int,
) -> dict[str, int | float]:
    """Calculate established fit statistics from an authoritative residual vector."""
    ndata = len(residuals)
    chisqr = sum(residuals**2)
    redchi = chisqr / max(1, ndata - nvarys)
    aic = chisqr + 2 * nvarys
    bic = chisqr + np.log(ndata) * nvarys
    _, ks_p_value = stats.kstest(residuals, "norm")
    ks_p_value = float(ks_p_value)
    pvalue: float = 1.0 - stats.chi2.cdf(chisqr, ndata - nvarys)
    return {
        "ndata": ndata,
        "nvarys": nvarys,
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
    nvarys: int,
) -> None:
    """Write fitting statistics to a file."""
    stats = calculate_statistics_from_residuals(residuals, nvarys)
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
    nvarys: int,
) -> None:
    """Write the results of the fit to output files."""
    print_writing_results(path)
    path.mkdir(parents=True, exist_ok=True)
    write_parameters(experiments, path, parameter_store=experiments.parameter_store)
    experiments.write(path)
    _write_statistics(
        experiments,
        path=path,
        residuals=residuals,
        nvarys=nvarys,
    )


def _write_simulation_files(
    experiments: Experiments,
    path: Path,
) -> None:
    """Write the results of the simulation to output files."""
    print_writing_results(path)
    path.mkdir(parents=True, exist_ok=True)
    write_parameters(experiments, path, parameter_store=experiments.parameter_store)
    experiments.write(path)


def _write_plots(experiments: Experiments, path: Path) -> None:
    """Plot the experimental and fitted data."""
    print_making_plots()

    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_)
    except KeyboardInterrupt:
        print_plotting_canceled()
        raise


def _write_simulation_plots(experiments: Experiments, path: Path) -> None:
    """Plot the experimental and fitted data."""
    print_making_plots()

    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot_simulation(path=path_)
    except KeyboardInterrupt:
        print_plotting_canceled()


def execute_post_fit(
    experiments: Experiments,
    path: Path,
    *,
    plot: bool = False,
    residuals: Array,
    nvarys: int,
) -> None:
    _write_files(experiments, path, residuals=residuals, nvarys=nvarys)
    if plot:
        _write_plots(experiments, path)


def execute_simulation(
    experiments: Experiments,
    path: Path,
    *,
    parameter_values: Mapping[str, float],
    plot: bool = False,
) -> None:
    experiments.prepare_for_simulation(parameter_values)
    _write_simulation_files(experiments, path)
    if plot:
        _write_simulation_plots(experiments, path)


def execute_post_fit_groups(
    experiments: Experiments,
    path: Path,
    plot: str,
    *,
    residuals: Array,
    nvarys: int,
) -> None:
    print_group_name("All groups")
    statistics = calculate_statistics_from_residuals(residuals, nvarys)
    print_chi2(statistics["chisqr"], statistics["redchi"])
    execute_post_fit(
        experiments,
        path / "All",
        plot=(plot != "nothing"),
        residuals=residuals,
        nvarys=nvarys,
    )


def print_header(
    grid: Iterable[str],
    *,
    parameter_store: ParameterStore,
) -> str:
    parameters = parameter_store.get_parameters(grid)
    header_pnames = " ".join(f"{parameters[param_id].param_name}" for param_id in grid)
    return f"# {header_pnames} [χ²]\n"


def print_values(values: Iterable[float], chisqr: float) -> str:
    body_values = " ".join(f"{value:.5e}" for value in values)
    return f"  {body_values} {chisqr:.5e}\n"
