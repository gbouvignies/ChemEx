from __future__ import annotations

from pathlib import Path
from typing import Self

import lmfit
import numpy as np
import pytest
from emcee.autocorr import AutocorrError

from chemex.configuration.methods import McmcSettings
from chemex.optimize import mcmc as mcmc_module
from chemex.optimize import mcmc_engine as mcmc_engine_module
from chemex.optimize.mcmc import (
    EffectiveMcmcSettings,
    McmcResult,
    McmcSummary,
    _apply_sample_window,
    _initial_positions,
    _result_from_emcee,
    resolve_mcmc_settings,
    run_mcmc,
    write_mcmc_outputs,
)
from chemex.optimize.mcmc_engine import (
    EmceeSamplerResult,
    McmcEvaluator,
    McmcProblem,
    evaluate_log_probability,
)
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.runtime import ExecutionSettings


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


class ResidualExperiments:
    def __init__(self, residuals: np.ndarray) -> None:
        self.parameter_store = DummyParameterStore()
        self.residual_values = residuals
        self.calls = 0

    def residuals(self, _params: lmfit.Parameters) -> np.ndarray:
        self.calls += 1
        return self.residual_values


class ConstraintExperiments:
    def __init__(self) -> None:
        self.parameter_store = DummyParameterStore()
        self.constrained_value: float | None = None

    def residuals(self, params: lmfit.Parameters) -> np.ndarray:
        self.constrained_value = float(params["twice_pb"].value)
        return np.array([self.constrained_value, 0.0])


def _single_parameter_problem(
    experiments: object,
    params: lmfit.Parameters,
) -> McmcProblem:
    return McmcProblem(
        experiments=experiments,
        params=params,
        var_names=("__PB",),
        bounds=np.array([[0.0, 1.0]]),
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


def test_log_probability_uses_weighted_vector_residuals() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    experiments = ResidualExperiments(np.array([1.0, 2.0]))
    evaluator = McmcEvaluator(
        experiments=experiments,
        params=params.copy(),
        var_names=("__PB",),
        bounds=np.array([[0.0, 1.0]]),
    )

    log_probability = evaluate_log_probability(np.array([0.2]), evaluator)

    assert log_probability == -2.5
    assert experiments.calls == 1


def test_log_probability_rejects_out_of_bounds_without_residual_call() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    experiments = ResidualExperiments(np.array([1.0, 2.0]))
    evaluator = McmcEvaluator(
        experiments=experiments,
        params=params.copy(),
        var_names=("__PB",),
        bounds=np.array([[0.0, 1.0]]),
    )

    log_probability = evaluate_log_probability(np.array([1.2]), evaluator)

    assert log_probability == -np.inf
    assert experiments.calls == 0


def test_log_probability_updates_constraints_before_residual_call() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    params.add("twice_pb", expr="2.0 * __PB")
    experiments = ConstraintExperiments()
    evaluator = McmcEvaluator(
        experiments=experiments,
        params=params.copy(),
        var_names=("__PB",),
        bounds=np.array([[0.0, 1.0]]),
    )

    log_probability = evaluate_log_probability(np.array([0.3]), evaluator)

    assert experiments.constrained_value == pytest.approx(0.6)
    assert log_probability == pytest.approx(-0.18)


def test_log_probability_raises_for_invalid_residuals() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    experiments = ResidualExperiments(np.array([np.nan]))
    evaluator = McmcEvaluator(
        experiments=experiments,
        params=params.copy(),
        var_names=("__PB",),
        bounds=np.array([[0.0, 1.0]]),
    )

    with pytest.raises(ValueError, match="NaN values detected"):
        evaluate_log_probability(np.array([0.2]), evaluator)


def test_initial_positions_are_seeded_reproducibly() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)

    first = _initial_positions(params, ("__PB",), nwalkers=4, seed=1234)
    second = _initial_positions(params, ("__PB",), nwalkers=4, seed=1234)

    assert np.array_equal(first, second)


def test_run_emcee_sampler_uses_serial_path_without_pool(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class FakeSampler:
        def __init__(
            self,
            nwalkers: int,
            ndim: int,
            log_prob_fn: object,
            *,
            pool: object | None,
        ) -> None:
            captured["sampler_pool"] = pool
            self.acceptance_fraction = np.full(nwalkers, 0.5)
            self._chain = np.zeros((1, nwalkers, ndim))
            self._lnprob = np.zeros((1, nwalkers))

        def run_mcmc(
            self,
            initial_positions: np.ndarray,
            steps: int,
            *,
            progress: bool,
        ) -> None:
            captured["run_mcmc"] = (initial_positions.copy(), steps, progress)

        def get_chain(self) -> np.ndarray:
            return self._chain

        def get_log_prob(self) -> np.ndarray:
            return self._lnprob

    def fail_pool(*_args: object, **_kwargs: object) -> object:
        msg = "Pool should not be created for a serial MCMC run"
        raise AssertionError(msg)

    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    problem = _single_parameter_problem(ResidualExperiments(np.array([0.0])), params)
    initial_positions = np.array([[0.1], [0.2]])

    monkeypatch.setattr(mcmc_engine_module.emcee, "EnsembleSampler", FakeSampler)
    monkeypatch.setattr(mcmc_engine_module.multiprocessing, "Pool", fail_pool)

    result = mcmc_engine_module.run_emcee_sampler(
        problem,
        initial_positions,
        steps=3,
        workers=1,
        seed=1234,
        progress=False,
    )

    assert captured["sampler_pool"] is None
    assert isinstance(captured["run_mcmc"], tuple)
    assert result.chain.shape == (1, 2, 1)


def test_run_emcee_sampler_runs_direct_emcee_serial() -> None:
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    problem = _single_parameter_problem(
        ResidualExperiments(np.array([0.0, 0.0])), params
    )
    initial_positions = np.array([[0.1], [0.2], [0.3], [0.4]])

    result = mcmc_engine_module.run_emcee_sampler(
        problem,
        initial_positions,
        steps=2,
        workers=1,
        seed=1234,
        progress=False,
    )

    assert result.var_names == ("__PB",)
    assert result.chain.shape == (2, 4, 1)
    assert result.lnprob.shape == (2, 4)
    assert result.acceptance_fraction.shape == (4,)


def test_run_emcee_sampler_initializes_multiprocessing_pool(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class FakePool:
        def __init__(
            self,
            *,
            processes: int,
            initializer: object,
            initargs: tuple[object, ...],
        ) -> None:
            captured["pool_processes"] = processes
            captured["pool_initializer"] = initializer
            captured["pool_initargs"] = initargs
            initializer_callable = initializer
            assert callable(initializer_callable)
            initializer_callable(*initargs)

        def __enter__(self) -> Self:
            return self

        def __exit__(
            self,
            _exc_type: object,
            _exc: object,
            _traceback: object,
        ) -> None:
            return None

        def map(self, func: object, iterable: object) -> list[object]:
            assert callable(func)
            return [func(item) for item in iterable]

    class FakeSampler:
        def __init__(
            self,
            nwalkers: int,
            ndim: int,
            log_prob_fn: object,
            *,
            pool: object | None,
        ) -> None:
            captured["sampler_pool"] = pool
            self.acceptance_fraction = np.full(nwalkers, 0.5)
            self._chain = np.zeros((1, nwalkers, ndim))
            self._lnprob = np.zeros((1, nwalkers))

        def run_mcmc(
            self,
            initial_positions: np.ndarray,
            steps: int,
            *,
            progress: bool,
        ) -> None:
            captured["run_mcmc"] = (initial_positions.copy(), steps, progress)

        def get_chain(self) -> np.ndarray:
            return self._chain

        def get_log_prob(self) -> np.ndarray:
            return self._lnprob

    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)
    problem = _single_parameter_problem(ResidualExperiments(np.array([0.0])), params)
    initial_positions = np.array([[0.1], [0.2]])

    monkeypatch.setattr(mcmc_engine_module.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(mcmc_engine_module.emcee, "EnsembleSampler", FakeSampler)

    result = mcmc_engine_module.run_emcee_sampler(
        problem,
        initial_positions,
        steps=3,
        workers=2,
        seed=1234,
        progress=False,
    )

    assert captured["pool_processes"] == 2
    assert captured["pool_initargs"] == (problem,)
    assert captured["sampler_pool"] is not None
    assert result.chain.shape == (1, 2, 1)


def test_run_mcmc_passes_execution_workers_to_direct_sampler(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}
    params = lmfit.Parameters()
    params.add("__PB", value=0.1, min=0.0, max=1.0, vary=True)

    result = McmcResult(
        var_names=("__PB",),
        chain=np.array([[[0.1], [0.2]]]),
        lnprob=np.array([[-1.0, -2.0]]),
        summary=(
            McmcSummary(
                parameter_id="__PB",
                mean=0.15,
                standard_deviation=0.05,
                median=0.15,
                eti_95_lower=0.10,
                eti_95_upper=0.20,
                lower_1sigma=0.10,
                upper_1sigma=0.20,
                stderr=0.05,
            ),
        ),
        correlations=np.array([[1.0]]),
        acceptance_fraction=np.array([0.5, 0.5]),
        autocorrelation_time=None,
        discarded_steps=0,
        burn_in_warning=None,
    )

    def fake_run_emcee_sampler(
        problem: McmcProblem,
        initial_positions: np.ndarray,
        *,
        steps: int,
        workers: int,
        seed: int | None,
        progress: bool,
    ) -> EmceeSamplerResult:
        captured["problem"] = problem
        captured["initial_positions"] = initial_positions
        captured["sampler_kwargs"] = {
            "steps": steps,
            "workers": workers,
            "seed": seed,
            "progress": progress,
        }
        return EmceeSamplerResult(
            var_names=result.var_names,
            chain=result.chain,
            lnprob=result.lnprob,
            acceptance_fraction=result.acceptance_fraction,
        )

    monkeypatch.setattr(mcmc_module, "run_emcee_sampler", fake_run_emcee_sampler)

    def fake_write_mcmc_outputs(*_args: object, **kwargs: object) -> None:
        captured["output_kwargs"] = kwargs
        captured["outputs_written"] = True

    monkeypatch.setattr(mcmc_module, "write_mcmc_outputs", fake_write_mcmc_outputs)

    run_mcmc(
        DummyExperiments(),
        params,
        McmcSettings(steps=10),
        tmp_path,
        execution=ExecutionSettings(workers=7),
    )

    sampler_kwargs = captured["sampler_kwargs"]
    assert isinstance(sampler_kwargs, dict)
    assert sampler_kwargs["workers"] == 7
    assert sampler_kwargs["steps"] == 10
    assert sampler_kwargs["progress"] is True
    assert isinstance(captured["problem"], McmcProblem)
    initial_positions = captured["initial_positions"]
    assert isinstance(initial_positions, np.ndarray)
    assert initial_positions.shape == (32, 1)
    assert captured["outputs_written"]
    output_kwargs = captured["output_kwargs"]
    assert isinstance(output_kwargs, dict)
    timings = output_kwargs["timings"]
    assert isinstance(timings, dict)
    assert timings["sampling_seconds"] >= 0.0
    assert timings["result_processing_seconds"] >= 0.0


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


def test_result_from_emcee_uses_tentative_tau_for_auto_burn(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    chain = np.arange(20.0).reshape(10, 2, 1)
    lnprob = np.arange(20.0).reshape(10, 2)
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

    def raise_autocorr_error(_chain: np.ndarray, *, quiet: bool) -> np.ndarray:
        _ = quiet
        raise AutocorrError(np.array([1.6]), "short chain")

    monkeypatch.setattr(
        mcmc_module.emcee_autocorr,
        "integrated_time",
        raise_autocorr_error,
    )

    result = _result_from_emcee(
        EmceeSamplerResult(
            var_names=("__PB",),
            chain=chain,
            lnprob=lnprob,
            acceptance_fraction=np.array([0.25, 0.5]),
        ),
        settings,
    )

    assert result.autocorrelation_time is None
    assert result.tentative_autocorrelation_time is not None
    assert np.array_equal(result.tentative_autocorrelation_time, np.array([1.6]))
    assert result.discarded_steps == 4
    assert "tentative automatic burn-in was applied" in str(result.burn_in_warning)
    assert "shorter than 50 times" in str(result.autocorrelation_warning)
    assert result.summary[0].effective_sample_size is None


def test_write_mcmc_outputs(tmp_path: Path) -> None:
    store = DummyParameterStore()
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
    assert "effective_sample_size = 2.00000e+00" in summary
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
    assert 'unbounded_parameters = ["[KEX_AB]"]' in diagnostics
    assert "autocorrelation_time = [2.00000e+00, 3.00000e+00]" in diagnostics
    assert 'autocorrelation_status = "reliable"' in diagnostics
    assert "recommended_min_steps_50tau = 150" in diagnostics
    assert 'plots_file = "plots.pdf"' in diagnostics
    assert (path_mcmc / "plots.pdf").stat().st_size > 0


def test_write_mcmc_outputs_reports_tentative_autocorrelation_time(
    tmp_path: Path,
) -> None:
    store = DummyParameterStore()
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
                lower_1sigma=0.32,
                upper_1sigma=0.58,
                stderr=0.13,
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

    write_mcmc_outputs(result, settings, tmp_path, store)

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
