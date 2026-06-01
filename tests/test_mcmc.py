from __future__ import annotations

from pathlib import Path

import lmfit
import numpy as np
import pytest

from chemex.configuration.methods import McmcSettings
from chemex.optimize.mcmc import (
    EffectiveMcmcSettings,
    McmcResult,
    McmcSummary,
    _apply_sample_window,
    resolve_mcmc_settings,
    run_mcmc,
    write_mcmc_outputs,
)
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting


class DummyParameterStore:
    def __init__(self) -> None:
        self.parameters = {
            "__PB": ParamSetting(ParamName("PB"), value=0.1, min=0.0, max=1.0),
            "__KEX_AB": ParamSetting(
                ParamName("KEX_AB"),
                value=200.0,
                min=1.0,
                max=5000.0,
            ),
        }
        self.updated = False

    def get_parameters(self, param_ids: tuple[str, ...]) -> dict[str, ParamSetting]:
        return {
            param_id: parameter
            for param_id, parameter in self.parameters.items()
            if param_id in param_ids
        }

    def update_from_parameters(self, _params: lmfit.Parameters) -> None:
        self.updated = True


class DummyExperiments:
    def __init__(self) -> None:
        self.parameter_store = DummyParameterStore()

    def residuals(self, _params: lmfit.Parameters) -> np.ndarray:
        msg = "Sampler should not run"
        raise AssertionError(msg)


def test_resolve_mcmc_settings_defaults_walkers_from_varying_parameters() -> None:
    settings = resolve_mcmc_settings(McmcSettings(steps=100), nvarys=3)

    assert settings.walkers == 32
    assert settings.burn == "auto"


def test_resolve_mcmc_settings_rejects_too_few_walkers() -> None:
    with pytest.raises(ValueError, match="at least 6 walkers"):
        resolve_mcmc_settings(McmcSettings(steps=100, walkers=5), nvarys=3)


def test_run_mcmc_skips_when_no_parameters_vary() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, vary=False)

    result = run_mcmc(
        DummyExperiments(),
        params,
        McmcSettings(steps=100),
        Path("unused"),
    )

    assert result is None


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
    store = DummyParameterStore()
    settings = EffectiveMcmcSettings(
        steps=4,
        burn="auto",
        thin=1,
        walkers=2,
        seed=1234,
        workers=1,
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
                lower_1sigma=0.12,
                upper_1sigma=0.38,
                stderr=0.13,
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
                lower_1sigma=210.0,
                upper_1sigma=340.0,
                stderr=65.0,
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
        store,
        unbounded_parameter_ids=("__KEX_AB",),
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
    assert "effective_sample_size = 2.00000e+00" in summary
    assert "[PB]\t[KEX_AB]\tlnprob" in samples
    assert "2.00000e-01\t2.50000e+02\t-2.00000e+00" in samples
    assert "[PB]" in correlations
    assert "5.00000e-01" in correlations
    assert 'sampler = "emcee via lmfit.Minimizer.emcee"' in diagnostics
    assert 'samples_file = "samples.tsv"' in diagnostics
    assert 'correlations_file = "correlations.tsv"' in diagnostics
    assert 'requested_burn = "auto"' in diagnostics
    assert "discarded_steps = 1" in diagnostics
    assert "retained_steps = 2" in diagnostics
    assert "retained_samples = 4" in diagnostics
    assert "min_effective_sample_size = 1.33333e+00" in diagnostics
    assert 'unbounded_parameters = ["[KEX_AB]"]' in diagnostics
    assert "autocorrelation_time = [2.00000e+00, 3.00000e+00]" in diagnostics
