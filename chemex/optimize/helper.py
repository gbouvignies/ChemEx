from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
from lmfit import Parameters as ParametersLF
from scipy import stats

from chemex.containers.experiments import Experiments
from chemex.messages import print_chi2
from chemex.messages import print_group_name
from chemex.messages import print_making_plots
from chemex.messages import print_writing_results
from chemex.parameters import database
from chemex.printers.parameters import write_parameters


def calculate_statistics(
    experiments: Experiments, params_lf: ParametersLF
) -> dict[str, int | float]:
    residuals = experiments.residuals(params_lf)
    ndata = len(residuals)
    nvarys = len(
        [param for param in params_lf.values() if param.vary and not param.expr]
    )
    chisqr = sum(residuals**2)
    redchi = chisqr / max(1, ndata - nvarys)
    _neg2_log_likel = ndata * np.log(chisqr / ndata)
    aic = _neg2_log_likel + 2 * nvarys
    bic = _neg2_log_likel + np.log(ndata) * nvarys
    _, ks_p_value = stats.kstest(residuals, "norm")
    pvalue: float = 1.0 - stats.chi2.cdf(chisqr, ndata - nvarys)  # type: ignore
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


def _write_statistics(experiments: Experiments, path: Path):
    """Write fitting statistics to a file."""
    params_lf = database.build_lmfit_params(experiments.param_ids)
    statistics = calculate_statistics(experiments, params_lf)
    filename = path / "statistics.toml"
    with open(filename, "w") as f:
        f.write(f"\"number of data points\"                = {statistics['ndata']}\n")
        f.write(f"\"number of variables\"                  = {statistics['nvarys']}\n")
        f.write(
            f"\"chi-square\"                           = {statistics['chisqr']: .5e}\n"
        )
        f.write(
            f"\"reduced-chi-square\"                   = {statistics['redchi']: .5e}\n"
        )
        f.write(
            f"\"chi-squared test\"                     = {statistics['pvalue']: .5e}\n"
        )
        f.write(
            f"\"Kolmogorov-Smirnov test\"              = {statistics['ks_pvalue']: .5e}\n"
        )
        f.write(
            f"\"Akaike Information Criterion (AIC)\"   = {statistics['aic']: .5e}\n"
        )
        f.write(
            f"\"Bayesian Information Criterion (BIC)\" = {statistics['bic']: .5e}\n"
        )


def _write_files(experiments: Experiments, path: Path):
    """Write the results of the fit to output files."""
    print_writing_results(path)
    path.mkdir(parents=True, exist_ok=True)
    write_parameters(experiments, path)
    experiments.write(path)
    _write_statistics(experiments, path=path)


def _write_plots(experiments: Experiments, path: Path):
    """Plot the experimental and fitted data."""

    print_making_plots()

    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_)
    except KeyboardInterrupt:
        print("  - Plotting cancelled\n")
    print("")


def execute_post_fit(
    experiments: Experiments,
    path: Path,
    plot: bool = False,
):
    _write_files(experiments, path)
    if plot:
        _write_plots(experiments, path)


def execute_post_fit_groups(experiments: Experiments, path: Path, plot: str) -> None:
    print_group_name("All groups")
    params_lf = database.build_lmfit_params(experiments.param_ids)
    statistics = calculate_statistics(experiments, params_lf)
    print_chi2(statistics["chisqr"], statistics["redchi"])
    execute_post_fit(experiments, path / "All", plot != "nothing")


def print_header(grid: Iterable[str]) -> str:
    parameters = database.get_parameters(grid)
    header_pnames = " ".join(f"{parameters[param_id].param_name}" for param_id in grid)
    return f"# {header_pnames} [χ²]\n"


def print_values(values: Iterable[float], chisqr: float) -> str:
    body_values = " ".join(f"{value:.5e}" for value in values)
    return f"  {body_values} {chisqr:.5e}\n"


def print_values_stat(
    params_lf: ParametersLF, fnames: Iterable[str], chisqr: float
) -> str:
    body_values_list = []
    for fname in fnames:
        if fname in params_lf:
            body_values_list.append(f"{params_lf[fname].value:12.5e}")
        else:
            body_values_list.append(f"{'--':^12s}")
    body_values = " ".join(body_values_list)
    return f"  {body_values} {chisqr:.5e}\n"
