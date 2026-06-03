"""Direct emcee execution for ChemEx MCMC sampling."""

from __future__ import annotations

import multiprocessing
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Protocol, cast

import emcee
import numpy as np
from lmfit.minimizer import coerce_float64
from lmfit.parameter import Parameters

from chemex.containers.experiments import Experiments
from chemex.typing import Array


class _PoolLike(Protocol):
    def map(
        self,
        func: Callable[[Array], float],
        iterable: Iterable[Array],
    ) -> Sequence[float]: ...


@dataclass(frozen=True, slots=True)
class McmcProblem:
    experiments: Experiments
    params: Parameters
    var_names: tuple[str, ...]
    bounds: Array


@dataclass(frozen=True, slots=True)
class EmceeSamplerResult:
    var_names: tuple[str, ...]
    chain: Array
    lnprob: Array
    acceptance_fraction: Array


@dataclass(slots=True)
class McmcEvaluator:
    experiments: Experiments
    params: Parameters
    var_names: tuple[str, ...]
    bounds: Array

    @classmethod
    def from_problem(cls, problem: McmcProblem) -> McmcEvaluator:
        return cls(
            experiments=problem.experiments,
            params=problem.params.copy(),
            var_names=problem.var_names,
            bounds=problem.bounds,
        )


@dataclass(slots=True)
class _WorkerState:
    evaluator: McmcEvaluator | None = None

    def initialize(self, problem: McmcProblem) -> None:
        self.evaluator = McmcEvaluator.from_problem(problem)

    def clear(self) -> None:
        self.evaluator = None


_WORKER_STATE = _WorkerState()


def _initialize_worker(problem: McmcProblem) -> None:
    _WORKER_STATE.initialize(problem)


def _clear_worker() -> None:
    _WORKER_STATE.clear()


def _parameter_within_bounds(theta: Array, bounds: Array) -> bool:
    return bool(np.all(theta >= bounds[:, 0]) and np.all(theta <= bounds[:, 1]))


def evaluate_log_probability(
    theta: Array,
    evaluator: McmcEvaluator,
) -> float:
    theta = np.asarray(theta, dtype=float)
    if not _parameter_within_bounds(theta, evaluator.bounds):
        return -np.inf

    for name, value in zip(evaluator.var_names, theta, strict=True):
        evaluator.params[name].value = float(value)
    evaluator.params.update_constraints()

    residuals = evaluator.experiments.residuals(evaluator.params)
    log_probability = coerce_float64(residuals, nan_policy="raise")
    if len(log_probability) == 0:
        return -1.0e100
    if log_probability.size > 1:
        return float(-0.5 * np.sum(log_probability * log_probability))
    return float(log_probability[0])


def _log_probability(theta: Array) -> float:
    evaluator = _WORKER_STATE.evaluator
    if evaluator is None:
        msg = "MCMC worker context is not initialized"
        raise RuntimeError(msg)
    return evaluate_log_probability(theta, evaluator)


def _seed_sampler(sampler: emcee.EnsembleSampler, seed: int | None) -> None:
    if seed is None:
        return
    rng = np.random.RandomState(seed)
    sampler.random_state = rng.get_state()


def _run_sampler(
    problem: McmcProblem,
    initial_positions: Array,
    *,
    steps: int,
    seed: int | None,
    progress: bool,
    pool: _PoolLike | None,
) -> EmceeSamplerResult:
    nwalkers, ndim = initial_positions.shape
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        _log_probability,
        pool=pool,
    )
    _seed_sampler(sampler, seed)
    sampler.run_mcmc(initial_positions, steps, progress=progress)
    return EmceeSamplerResult(
        var_names=problem.var_names,
        chain=np.asarray(sampler.get_chain(), dtype=float),
        lnprob=np.asarray(sampler.get_log_prob(), dtype=float),
        acceptance_fraction=np.asarray(sampler.acceptance_fraction, dtype=float),
    )


def run_emcee_sampler(
    problem: McmcProblem,
    initial_positions: Array,
    *,
    steps: int,
    workers: int,
    seed: int | None,
    progress: bool = True,
) -> EmceeSamplerResult:
    if workers < 1:
        msg = "MCMC requires at least one worker"
        raise ValueError(msg)

    if workers == 1:
        _initialize_worker(problem)
        try:
            return _run_sampler(
                problem,
                initial_positions,
                steps=steps,
                seed=seed,
                progress=progress,
                pool=None,
            )
        finally:
            _clear_worker()

    try:
        with multiprocessing.Pool(
            processes=workers,
            initializer=_initialize_worker,
            initargs=(problem,),
        ) as pool:
            return _run_sampler(
                problem,
                initial_positions,
                steps=steps,
                seed=seed,
                progress=progress,
                pool=cast("_PoolLike", pool),
            )
    finally:
        _clear_worker()
