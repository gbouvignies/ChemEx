from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemex.configuration.methods import Method, Statistics


def test_fitmethod_defaults_to_canonical_trf() -> None:
    method = Method()

    assert method.fitmethod == "trf"


@pytest.mark.parametrize("fitmethod", ("trf", "TRF", "least_squares"))
def test_fitmethod_accepts_only_trf_surface_and_canonicalizes_alias(
    fitmethod: str,
) -> None:
    method = Method.model_validate({"FITMETHOD": fitmethod})

    assert method.fitmethod == "trf"


@pytest.mark.parametrize(
    "fitmethod",
    ("leastsq", "differential_evolution", "nelder", "arbitrary"),
)
def test_fitmethod_rejects_legacy_and_arbitrary_spellings(fitmethod: str) -> None:
    with pytest.raises(ValidationError, match="FITMETHOD supports only 'trf'"):
        Method.model_validate({"FITMETHOD": fitmethod})


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


def test_statistics_rejects_mcmc_update_parameters_state_mutation() -> None:
    with pytest.raises(
        ValidationError, match="cannot mutate the committed central fit"
    ):
        Statistics.model_validate({"MCMC": {"STEPS": 100, "UPDATE_PARAMETERS": True}})


@pytest.mark.parametrize("seed", (0, (1 << 64) - 1))
def test_statistics_accepts_unsigned_64_bit_mcmc_seed_boundaries(seed: int) -> None:
    statistics = Statistics.model_validate({"MCMC": {"STEPS": 100, "SEED": seed}})

    assert statistics.mcmc is not None
    assert statistics.mcmc.seed == seed


@pytest.mark.parametrize("seed", (-1, 1 << 64))
def test_statistics_rejects_out_of_range_mcmc_seed(seed: int) -> None:
    with pytest.raises(ValidationError, match="unsigned 64-bit integer"):
        Statistics.model_validate({"MCMC": {"STEPS": 100, "SEED": seed}})
