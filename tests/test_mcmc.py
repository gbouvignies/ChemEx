from __future__ import annotations

import pytest

from chemex.configuration.method_plan import McmcRequest
from chemex.configuration.methods import McmcSettings
from chemex.optimize.mcmc import resolve_mcmc_request, resolve_mcmc_settings
from chemex.runtime import ExecutionSettings


def test_resolve_mcmc_settings_defaults_walkers_from_varying_parameters() -> None:
    settings = resolve_mcmc_settings(McmcSettings(steps=100), nvarys=3)

    assert settings.walkers == 32
    assert settings.burn == "auto"
    assert settings.workers == 1


def test_resolve_mcmc_settings_inherits_execution_workers() -> None:
    settings = resolve_mcmc_settings(
        McmcSettings(steps=100),
        nvarys=3,
        execution=ExecutionSettings(workers=4),
    )

    assert settings.workers == 4


def test_resolve_mcmc_settings_method_workers_override_execution() -> None:
    settings = resolve_mcmc_settings(
        McmcSettings(steps=100, workers=2),
        nvarys=3,
        execution=ExecutionSettings(workers=4),
    )

    assert settings.workers == 2


def test_resolve_mcmc_settings_rejects_too_few_walkers() -> None:
    with pytest.raises(ValueError, match="at least 6 walkers"):
        resolve_mcmc_settings(McmcSettings(steps=100, walkers=5), nvarys=3)


def test_resolve_mcmc_request_preserves_canonical_compatibility_values() -> None:
    settings = resolve_mcmc_request(
        McmcRequest(
            steps=100,
            burn=10,
            thin=3,
            walkers=8,
            workers=2,
        ),
        nvarys=3,
        execution=ExecutionSettings(workers=4),
    )

    assert settings.burn == 10
    assert settings.thin == 3
    assert settings.walkers == 8
    assert settings.workers == 2
    assert not settings.update_parameters
