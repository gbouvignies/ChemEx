from __future__ import annotations

import math
import tomllib
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from chemex.containers.experiments import Experiments
from chemex.optimize.helper import (
    _write_statistics,
    calculate_statistics_from_residuals,
)


@pytest.mark.parametrize(
    ("observations", "controlled", "normalizations", "expected_dof"),
    (
        (10, 2, 0, 8),
        (10, 2, 1, 7),
        (10, 2, 3, 5),
    ),
)
def test_fit_statistics_count_all_estimated_quantities(
    observations: int,
    controlled: int,
    normalizations: int,
    expected_dof: int,
) -> None:
    residuals = np.ones(observations)

    result = calculate_statistics_from_residuals(
        residuals,
        controlled_coordinate_count=controlled,
        profiled_normalization_count=normalizations,
    )

    estimated = controlled + normalizations
    assert result["ndata"] == observations
    assert result["nvarys"] == controlled
    assert result["dof"] == expected_dof
    assert result["redchi"] == observations / expected_dof
    assert result["pvalue"] == 1.0 - stats.chi2.cdf(observations, expected_dof)
    assert result["aic"] == observations + 2 * estimated
    assert result["bic"] == observations + math.log(observations) * estimated


@pytest.mark.parametrize(("normalizations", "expected_dof"), ((1, 0), (2, -1)))
def test_nonpositive_residual_dof_makes_dof_dependent_statistics_unavailable(
    normalizations: int,
    expected_dof: int,
) -> None:
    residuals = np.ones(2)

    result = calculate_statistics_from_residuals(
        residuals,
        controlled_coordinate_count=1,
        profiled_normalization_count=normalizations,
    )

    estimated = 1 + normalizations
    assert result["dof"] == expected_dof
    assert result["chisqr"] == 2.0
    assert math.isnan(result["redchi"])
    assert math.isnan(result["pvalue"])
    assert result["aic"] == 2.0 + 2 * estimated
    assert result["bic"] == 2.0 + math.log(2.0) * estimated
    assert result["ks_pvalue"] == stats.kstest(residuals, "norm").pvalue


def test_statistics_toml_serializes_nonpositive_dof_as_nan(tmp_path: Path) -> None:
    _write_statistics(
        Experiments(),
        tmp_path,
        residuals=np.ones(2),
        controlled_coordinate_count=1,
        profiled_normalization_count=1,
    )

    text = (tmp_path / "statistics.toml").read_text(encoding="utf-8")
    result = tomllib.loads(text)

    assert '"reduced-chi-square"                   =  nan' in text
    assert '"chi-squared test"                     =  nan' in text
    assert math.isnan(result["reduced-chi-square"])
    assert math.isnan(result["chi-squared test"])
    assert result["chi-square"] == 2.0
    assert result["Akaike Information Criterion (AIC)"] == 6.0
    assert result["Bayesian Information Criterion (BIC)"] == pytest.approx(
        2.0 + 2.0 * math.log(2.0),
        rel=0.0,
        abs=5.0e-5,
    )


def test_statistics_toml_preserves_variable_count_and_uses_effective_parameters(
    tmp_path: Path,
) -> None:
    _write_statistics(
        Experiments(),
        tmp_path,
        residuals=np.ones(7),
        controlled_coordinate_count=1,
        profiled_normalization_count=1,
    )

    result = tomllib.loads((tmp_path / "statistics.toml").read_text(encoding="utf-8"))

    assert result["number of data points"] == 7
    assert result["number of variables"] == 1
    assert result["chi-square"] == 7.0
    assert result["reduced-chi-square"] == 7.0 / 5.0
    assert result["chi-squared test"] == pytest.approx(
        1.0 - stats.chi2.cdf(7.0, 5),
        rel=0.0,
        abs=5.0e-7,
    )
    assert result["Kolmogorov-Smirnov test"] == pytest.approx(
        stats.kstest(np.ones(7), "norm").pvalue,
        rel=0.0,
        abs=5.0e-7,
    )
    assert result["Akaike Information Criterion (AIC)"] == 11.0
    assert result["Bayesian Information Criterion (BIC)"] == pytest.approx(
        7.0 + 2.0 * math.log(7.0),
        rel=0.0,
        abs=5.0e-5,
    )
