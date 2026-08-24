from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.methods import McmcSettings
from chemex.optimize.mcmc import (
    EffectiveMcmcSettings,
    McmcResult,
    McmcSummary,
    _apply_sample_window,
    resolve_mcmc_settings,
    write_mcmc_outputs,
)
from chemex.parameters.parameterization import (
    ParameterDeclaration,
    SealedParameterDeclarations,
    SealedParameterModel,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)
from chemex.runtime import ExecutionSettings


def _parameter_model() -> SealedParameterModel:
    definitions = SealedDefinitions(
        (
            ParamDefinition("__PB", "PB", "", (), 0.1, 0.0, 1.0),
            ParamDefinition(
                "__KEX_AB",
                "KEX_AB",
                "",
                (),
                200.0,
                1.0,
                5000.0,
            ),
        ),
        {},
    )
    configuration = SealedConfiguration(
        (
            ParamConfig("__PB", 0.1, 0.0, 1.0),
            ParamConfig("__KEX_AB", 200.0, 1.0, 5000.0),
        ),
        {},
        definitions.identity,
    )
    declarations = SealedParameterDeclarations(
        (
            ParameterDeclaration("__PB", True, "", False),
            ParameterDeclaration("__KEX_AB", True, "", False),
        )
    )
    return SealedParameterModel(
        "test",
        "test-model",
        definitions,
        configuration,
        declarations,
    )


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


def test_apply_sample_window_uses_auto_burn_from_autocorrelation_time() -> None:
    chain = np.arange(20.0).reshape(10, 2, 1)
    lnprob = np.arange(20.0).reshape(10, 2)

    retained_chain, retained_lnprob, discarded_steps, warning = _apply_sample_window(
        chain,
        lnprob,
        burn="auto",
        thin=2,
        autocorrelation_time=np.array([1.6]),
    )

    assert discarded_steps == 4
    assert warning is None
    assert retained_chain.shape == (3, 2, 1)
    assert np.array_equal(retained_chain[:, :, 0], chain[4::2, :, 0])
    assert np.array_equal(retained_lnprob, lnprob[4::2])


def test_apply_sample_window_keeps_samples_when_auto_burn_unavailable() -> None:
    chain = np.arange(12.0).reshape(3, 2, 2)
    lnprob = np.arange(6.0).reshape(3, 2)

    retained_chain, retained_lnprob, discarded_steps, warning = _apply_sample_window(
        chain,
        lnprob,
        burn="auto",
        thin=1,
        autocorrelation_time=None,
    )

    assert discarded_steps == 0
    assert "autocorrelation time unavailable" in str(warning)
    assert np.array_equal(retained_chain, chain)
    assert np.array_equal(retained_lnprob, lnprob)


def test_write_mcmc_outputs(tmp_path: Path) -> None:
    parameter_model = _parameter_model()
    settings = EffectiveMcmcSettings(
        steps=4,
        burn="auto",
        thin=1,
        walkers=2,
        seed=1234,
        workers=1,
        native_threads=None,
        update_parameters=False,
    )
    result = McmcResult(
        var_names=("__PB", "__KEX_AB"),
        chain=np.array(
            [
                [[0.10, 200.0], [0.20, 250.0]],
                [[0.30, 300.0], [0.40, 350.0]],
            ],
        ),
        lnprob=np.array([[-1.0, -2.0], [-3.0, -4.0]]),
        summary=(
            McmcSummary(
                parameter_id="__PB",
                mean=0.25,
                standard_deviation=0.13,
                median=0.25,
                eti_95_lower=0.11,
                eti_95_upper=0.39,
                credible_interval_68_lower=0.12,
                credible_interval_68_upper=0.38,
                half_credible_interval_68_width=0.13,
                effective_sample_size=2.0,
                mcse_mean=0.09,
            ),
            McmcSummary(
                parameter_id="__KEX_AB",
                mean=275.0,
                standard_deviation=65.0,
                median=275.0,
                eti_95_lower=205.0,
                eti_95_upper=345.0,
                credible_interval_68_lower=210.0,
                credible_interval_68_upper=340.0,
                half_credible_interval_68_width=65.0,
                effective_sample_size=1.3333333333,
                mcse_mean=56.2916512459,
            ),
        ),
        correlations=np.array([[1.0, 0.5], [0.5, 1.0]]),
        acceptance_fraction=np.array([0.25, 0.50]),
        autocorrelation_time=np.array([2.0, 3.0]),
        discarded_steps=1,
        burn_in_warning=None,
    )

    write_mcmc_outputs(
        result,
        settings,
        tmp_path,
        parameter_model,
        timings={"sampling_seconds": 1.25, "result_processing_seconds": 0.5},
    )

    path_mcmc = tmp_path / "Statistics" / "MCMC"
    summary = (path_mcmc / "summary.toml").read_text(encoding="utf-8")
    samples = (path_mcmc / "samples.tsv").read_text(encoding="utf-8")
    correlations = (path_mcmc / "correlations.tsv").read_text(
        encoding="utf-8",
    )
    diagnostics = (path_mcmc / "diagnostics.toml").read_text(
        encoding="utf-8",
    )

    assert '["PB"]' in summary
    assert 'prior = "uniform"' in summary
    assert "prior_lower = 0.00000e+00" in summary
    assert 'credible_interval = "95% equal-tailed"' in summary
    assert "median = 2.50000e-01" in summary
    assert "eti_95_lower = 1.10000e-01" in summary
    assert "credible_interval_68_lower = 1.20000e-01" in summary
    assert "credible_interval_68_upper = 3.80000e-01" in summary
    assert "half_credible_interval_68_width = 1.30000e-01" in summary
    assert "effective_sample_size = 2.00000e+00" in summary
    assert "mcse_mean = 9.00000e-02" in summary
    assert "stderr" not in summary
    assert "[PB]\t[KEX_AB]\tlnprob" in samples
    assert "2.00000e-01\t2.50000e+02\t-2.00000e+00" in samples
    assert "[PB]" in correlations
    assert "5.00000e-01" in correlations
    assert 'sampler = "emcee via ChemEx direct EnsembleSampler"' in diagnostics
    assert 'samples_file = "samples.tsv"' in diagnostics
    assert 'correlations_file = "correlations.tsv"' in diagnostics
    assert 'requested_burn = "auto"' in diagnostics
    assert "discarded_steps = 1" in diagnostics
    assert "retained_steps = 2" in diagnostics
    assert "retained_samples = 4" in diagnostics
    assert "workers = 1" in diagnostics
    assert "sampling_seconds = 1.25000e+00" in diagnostics
    assert "result_processing_seconds = 5.00000e-01" in diagnostics
    assert "output_summary_seconds" in diagnostics
    assert "output_samples_seconds" in diagnostics
    assert "output_correlations_seconds" in diagnostics
    assert "output_plots_seconds" in diagnostics
    assert "output_total_seconds" in diagnostics
    assert "total_seconds" in diagnostics
    assert "min_effective_sample_size = 1.33333e+00" in diagnostics
    assert "unbounded_parameters = []" in diagnostics
    assert "autocorrelation_time = [2.00000e+00, 3.00000e+00]" in diagnostics
    assert 'autocorrelation_status = "reliable"' in diagnostics
    assert "recommended_min_steps_50tau = 150" in diagnostics
    assert 'plots_file = "plots.pdf"' in diagnostics
    assert (path_mcmc / "plots.pdf").stat().st_size > 0


def test_write_mcmc_outputs_reports_tentative_autocorrelation_time(
    tmp_path: Path,
) -> None:
    parameter_model = _parameter_model()
    settings = EffectiveMcmcSettings(
        steps=10,
        burn="auto",
        thin=1,
        walkers=2,
        seed=None,
        workers=1,
        native_threads=None,
        update_parameters=False,
    )
    result = McmcResult(
        var_names=("__PB",),
        chain=np.array([[[0.30], [0.40]], [[0.50], [0.60]]]),
        lnprob=np.array([[-3.0, -4.0], [-5.0, -6.0]]),
        summary=(
            McmcSummary(
                parameter_id="__PB",
                mean=0.45,
                standard_deviation=0.13,
                median=0.45,
                eti_95_lower=0.31,
                eti_95_upper=0.59,
                credible_interval_68_lower=0.32,
                credible_interval_68_upper=0.58,
                half_credible_interval_68_width=0.13,
            ),
        ),
        correlations=np.array([[1.0]]),
        acceptance_fraction=np.array([0.25, 0.50]),
        autocorrelation_time=None,
        discarded_steps=4,
        burn_in_warning=(
            "autocorrelation time estimate is unreliable; "
            "tentative automatic burn-in was applied"
        ),
        tentative_autocorrelation_time=np.array([1.6]),
        autocorrelation_warning=(
            "chain shorter than 50 times the autocorrelation time; "
            "tentative estimate reported"
        ),
        raw_chain=np.arange(20.0).reshape(10, 2, 1),
        raw_lnprob=np.arange(20.0).reshape(10, 2),
    )

    write_mcmc_outputs(result, settings, tmp_path, parameter_model)

    diagnostics = (tmp_path / "Statistics" / "MCMC" / "diagnostics.toml").read_text(
        encoding="utf-8"
    )
    summary = (tmp_path / "Statistics" / "MCMC" / "summary.toml").read_text(
        encoding="utf-8",
    )

    assert "effective_sample_size" not in summary
    assert 'autocorrelation_status = "unreliable_short_chain"' in diagnostics
    assert "autocorrelation_time_tentative = [1.60000e+00]" in diagnostics
    assert "max_autocorrelation_time_tentative = 1.60000e+00" in diagnostics
    assert "steps_over_max_autocorrelation_time = 6.25000e+00" in diagnostics
    assert "recommended_min_steps_50tau = 80" in diagnostics
    assert "recommended_min_steps_100tau = 160" in diagnostics
    assert "effective_sample_size_warning" in diagnostics
