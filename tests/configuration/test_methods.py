from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemex.configuration.methods import Method, Statistics


def test_statistics_parse_mcmc_short_form() -> None:
    statistics = Statistics.model_validate({"MCMC": 5000})

    assert statistics.mcmc is not None
    assert statistics.mcmc.steps == 5000
    assert statistics.mcmc.burn == "auto"
    assert statistics.mcmc.thin == 1
    assert statistics.mcmc.workers is None


def test_statistics_parse_mcmc_expanded_form() -> None:
    method = Method.model_validate(
        {
            "STATISTICS": {
                "MCMC": {
                    "STEPS": 5000,
                    "BURN": 1000,
                    "THIN": 10,
                    "WALKERS": 64,
                    "SEED": 1234,
                    "WORKERS": 2,
                },
            },
        },
    )

    assert method.statistics is not None
    assert method.statistics.mcmc is not None
    assert method.statistics.mcmc.steps == 5000
    assert method.statistics.mcmc.burn == 1000
    assert method.statistics.mcmc.thin == 10
    assert method.statistics.mcmc.walkers == 64
    assert method.statistics.mcmc.seed == 1234
    assert method.statistics.mcmc.workers == 2


def test_statistics_rejects_mcmc_burn_greater_than_steps() -> None:
    with pytest.raises(ValidationError, match="burn must be smaller than steps"):
        Statistics.model_validate({"MCMC": {"STEPS": 100, "BURN": 100}})


def test_statistics_parse_mcmc_auto_burn_case_insensitive() -> None:
    statistics = Statistics.model_validate(
        {"MCMC": {"STEPS": 100, "BURN": "AUTO"}},
    )

    assert statistics.mcmc is not None
    assert statistics.mcmc.burn == "auto"
