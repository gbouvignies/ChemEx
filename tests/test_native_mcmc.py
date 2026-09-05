"""Behavioral qualification tests for native fixed-topology MCMC evidence (#600)."""

from __future__ import annotations

import copy
import dataclasses
import multiprocessing
import os
import pickle
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pydantic import BaseModel
from scipy.stats import truncnorm

import chemex.optimize.native_mcmc as native_mcmc
from chemex.configuration.method_plan import McmcRequest
from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.nmr.spectrometer import Spectrometer
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    OptimizationProblem,
)
from chemex.optimize.mcmc import EffectiveMcmcSettings, write_mcmc_outputs
from chemex.optimize.native_mcmc import (
    EnsembleState,
    ExpertMcmcPolicy,
    InitializationKind,
    McmcAnalysisFailureCategory,
    McmcAnalysisStatus,
    McmcAutocorrelationStatus,
    McmcConstructionError,
    McmcDiagnosticReason,
    McmcDiagnosticStatus,
    McmcEvidence,
    McmcEvidenceLifecycle,
    McmcExecutionStage,
    McmcInitializationOutcome,
    McmcOperationTerminal,
    McmcPlan,
    McmcPolicyKind,
    PosteriorSampleDisposition,
    ProposalKind,
    ResolvedMcmcPolicy,
    build_accepted_point_ensemble,
    build_bounded_latin_hypercube,
    derive_mcmc_analysis_result,
    derive_mcmc_diagnostics,
    derive_mcmc_operation_diagnostics,
    derive_posterior_sample_evidence,
    derive_posterior_summary,
    derive_retained_sample_view,
    execute_mcmc_evidence,
    resolve_product_mcmc_policy,
    validate_raw_mcmc_capture,
)
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ConstraintProgram,
    ParameterDeclaration,
    ParameterRole,
    ScientificFunctionBinder,
    SealedParameterDeclarations,
    SealedParameterModel,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.printers.data import Printer
from chemex.runtime import ExecutionSettings
from chemex.runtime.execution import NATIVE_THREAD_ENV_VARS
from chemex.typing import Array


class _KernelSettings(BaseModel):
    kind: str = "linear-test-kernel"


class _LinearSpectrometer:
    def __init__(self) -> None:
        self.spin_system = SpinSystem.from_name("1N")
        self.values = {"a": 0.0, "b": 0.0}

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _LinearSpectrometer:
        return deepcopy(self)

    def native_kernel_descriptor(self) -> dict[str, object]:
        return {"kind": "linear-test-spectrometer"}


class _LinearPulseSequence:
    settings = _KernelSettings()

    def __init__(
        self,
        evaluation_observer: Any | None = None,
        evaluation_barrier: Any | None = None,
        fail_outside_pid: int | None = None,
        *,
        source_pid: int | None = None,
    ) -> None:
        self._evaluation_observer = evaluation_observer
        self._evaluation_barrier = evaluation_barrier
        self._fail_outside_pid = fail_outside_pid
        self._source_pid = os.getpid() if source_pid is None else source_pid
        self._barrier_used = False

    def __deepcopy__(self, _memo: dict[int, object]) -> _LinearPulseSequence:
        return type(self)(
            self._evaluation_observer,
            self._evaluation_barrier,
            self._fail_outside_pid,
            source_pid=self._source_pid,
        )

    def calculate(self, spectrometer: _LinearSpectrometer, data: Data) -> Array:
        if (
            self._evaluation_barrier is not None
            and os.getpid() != self._source_pid
            and not self._barrier_used
        ):
            self._barrier_used = True
            self._evaluation_barrier.wait(timeout=10.0)
        if self._evaluation_observer is not None:
            self._evaluation_observer.append(
                (
                    os.getpid(),
                    id(self),
                    tuple(os.environ.get(name) for name in NATIVE_THREAD_ENV_VARS),
                )
            )
        if self._fail_outside_pid is not None and os.getpid() != self._fail_outside_pid:
            raise RuntimeError("worker kernel failure")
        return spectrometer.values["a"] + spectrometer.values["b"] * np.asarray(
            data.metadata,
            dtype=np.float64,
        )

    def is_reference(self, metadata: Array) -> Array:
        return metadata < 0.0


def test_linear_pulse_sequence_deepcopy_preserves_source_process_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_pid = os.getpid()
    pulse_sequence = _LinearPulseSequence()
    pulse_sequence._barrier_used = True
    monkeypatch.setattr(os, "getpid", lambda: source_pid + 1)

    worker_copy = deepcopy(pulse_sequence)

    assert worker_copy._source_pid == source_pid
    assert not worker_copy._barrier_used


def _native_context(
    *,
    fit_b: bool = True,
    evaluation_observer: Any | None = None,
    evaluation_barrier: Any | None = None,
    fail_outside_pid: int | None = None,
) -> tuple[
    AcceptedFitResult,
    OptimizationProblem,
    ActiveParameterization,
    EvaluationEngine,
]:
    binder = ScientificFunctionBinder("qualification", {})
    program = ConstraintProgram(
        "parameter-model",
        "model",
        "definitions",
        "configuration",
        binder.identity,
        ("A", "B"),
        ("A", "B"),
        (),
        (),
        (),
    )
    parameterization = ActiveParameterization(
        program,
        binder,
        "source-occurrence",
        4,
        (
            ("A", ParameterRole.FIT),
            ("B", ParameterRole.FIT if fit_b else ParameterRole.FIX),
        ),
    )
    snapshot = AnalysisValuesSnapshot(
        "source-occurrence",
        "model",
        "definitions",
        "configuration",
        4,
        (("A", 0.5), ("B", 1.5)),
    )
    data = Data(
        exp=np.asarray((1.0, 3.0, 5.0), dtype=np.float64),
        err=np.ones(3, dtype=np.float64),
        metadata=np.asarray((0.0, 1.0, 2.0), dtype=np.float64),
    )
    profile = Profile(
        data,
        cast("Spectrometer", _LinearSpectrometer()),
        cast(
            "PulseSequence",
            _LinearPulseSequence(
                evaluation_observer,
                evaluation_barrier,
                fail_outside_pid,
            ),
        ),
        {"a": "A", "b": "B"},
        cast("Printer", None),
        is_scaled=False,
    )
    experiments = cast("Any", (SimpleNamespace(profiles=(profile,)),))
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    controlled_ids = ("A", "B") if fit_b else ("A",)
    held_items = () if fit_b else (("B", 1.5),)
    start = (0.5, 1.5) if fit_b else (0.5,)
    lower = (0.0, 1.0) if fit_b else (0.0,)
    upper = (2.0, 3.0) if fit_b else (2.0,)
    accepted_vector = (1.0, 2.0) if fit_b else (1.0,)
    problem = OptimizationProblem(
        engine.plan.identity,
        parameterization.identity,
        parameterization.evaluator_identity,
        program.fingerprint,
        "configuration",
        snapshot,
        (("A", 0.5), ("B", 1.5)),
        controlled_ids,
        held_items,
        start,
        lower,
        upper,
        ("A", "B"),
    )
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame(accepted_vector, parameterization),
    )
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
    accepted = AcceptedFitResult.for_qualification(
        occurrence_identity="accepted-occurrence",
        problem_identity=problem.identity,
        invocation_identity="accepted-invocation",
        execution_identity="accepted-execution",
        materialization_identity="accepted-materialization",
        parameterization_identity=parameterization.identity,
        evaluator_parameterization_identity=parameterization.evaluator_identity,
        source_occurrence_identity=snapshot.occurrence_identity,
        source_revision=snapshot.revision,
        controlled_ids=problem.controlled_ids,
        vector=accepted_vector,
        chi_square=float(np.sum(result.residuals**2)),
        evaluation_result=result,
        commit_scope=problem.commit_scope,
        commit_items=result.resolved_values.ordered_items(),
        origin_context_identity="accepted-origin",
    )
    return accepted, problem, parameterization, engine


def _resolved_policy() -> ResolvedMcmcPolicy:
    return ExpertMcmcPolicy(
        burn_steps=2,
        retained_steps=3,
        walkers=8,
        expert_provenance="qualification-test-expert",
    ).resolve(dimension=2, root_seed=1234)


def _product_parameter_model() -> SealedParameterModel:
    definitions = SealedDefinitions(
        (
            ParamDefinition("A", "A", "", (), 0.5, 0.0, 2.0),
            ParamDefinition("B", "B", "", (), 1.5, 1.0, 3.0),
        ),
        {},
    )
    configuration = SealedConfiguration(
        (
            ParamConfig("A", 0.5, 0.0, 2.0),
            ParamConfig("B", 1.5, 1.0, 3.0),
        ),
        {},
        definitions.identity,
    )
    declarations = SealedParameterDeclarations(
        (
            ParameterDeclaration(
                "A", True, requires_independent=True, fits_by_default=True
            ),
            ParameterDeclaration(
                "B", True, requires_independent=True, fits_by_default=True
            ),
        )
    )
    return SealedParameterModel(
        "test",
        "test-model",
        definitions,
        configuration,
        declarations,
    )


def _same_semantics_foreign_occurrence(
    accepted: AcceptedFitResult,
) -> AcceptedFitResult:
    return AcceptedFitResult.for_qualification(
        occurrence_identity="foreign-accepted-occurrence",
        problem_identity=accepted.problem_identity,
        invocation_identity=accepted.invocation_identity,
        execution_identity=accepted.execution_identity,
        materialization_identity=accepted.materialization_identity,
        parameterization_identity=accepted.parameterization_identity,
        evaluator_parameterization_identity=(
            accepted.evaluator_parameterization_identity
        ),
        source_occurrence_identity=accepted.source_occurrence_identity,
        source_revision=accepted.source_revision,
        controlled_ids=accepted.controlled_ids,
        vector=accepted.vector,
        chi_square=accepted.chi_square,
        evaluation_result=accepted.evaluation_result,
        commit_scope=accepted.commit_scope,
        commit_items=accepted.commit_items,
        origin_context_identity=accepted.origin_context_identity,
    )


def test_seeded_native_mcmc_workers_execute_in_processes_without_changing_chain() -> (
    None
):
    context = multiprocessing.get_context("spawn")
    with context.Manager() as manager:
        two_worker_observations = manager.list()
        accepted, plan = _plan_context(
            evaluation_observer=two_worker_observations,
            evaluation_barrier=manager.Barrier(2),
        )
        two_worker_observations[:] = []

        parallel = execute_mcmc_evidence(
            accepted,
            plan,
            execution=ExecutionSettings(workers=2),
        )
        parallel_observations = set(two_worker_observations)

        four_worker_observations_proxy = manager.list()
        four_worker_accepted, four_worker_plan = _plan_context(
            evaluation_observer=four_worker_observations_proxy,
            evaluation_barrier=manager.Barrier(4),
        )
        four_worker_observations_proxy[:] = []
        four_worker = execute_mcmc_evidence(
            four_worker_accepted,
            four_worker_plan,
            execution=ExecutionSettings(workers=4),
        )
        four_worker_observations = set(four_worker_observations_proxy)

        serial_observations_proxy = manager.list()
        serial_accepted, serial_plan = _plan_context(
            evaluation_observer=serial_observations_proxy,
        )
        serial_observations_proxy[:] = []
        serial = execute_mcmc_evidence(
            serial_accepted,
            serial_plan,
            execution=ExecutionSettings(workers=1),
        )
        serial_observations = set(serial_observations_proxy)

    parallel_pids = {pid for pid, _evaluator, _environment in parallel_observations}
    serial_pids = {pid for pid, _evaluator, _environment in serial_observations}
    four_worker_pids = {
        pid for pid, _evaluator, _environment in four_worker_observations
    }

    assert parallel.terminal is McmcOperationTerminal.COMPLETED
    assert four_worker.terminal is McmcOperationTerminal.COMPLETED
    assert serial.terminal is McmcOperationTerminal.COMPLETED
    assert parallel.evidence is not None
    assert four_worker.evidence is not None
    assert serial.evidence is not None
    assert len(parallel_pids) == 2
    assert len(four_worker_pids) == 4
    assert os.getpid() not in parallel_pids
    assert os.getpid() not in four_worker_pids
    assert serial_pids == {os.getpid()}
    assert (parallel_pids | four_worker_pids).isdisjoint(
        child.pid for child in multiprocessing.active_children()
    )
    assert {environment for _pid, _evaluator, environment in parallel_observations} == {
        tuple("1" for _name in NATIVE_THREAD_ENV_VARS)
    }
    assert all(
        len(
            {
                evaluator
                for observed_pid, evaluator, _environment in parallel_observations
                if observed_pid == pid
            }
        )
        == 1
        for pid in parallel_pids
    )
    assert parallel.evidence.states == serial.evidence.states
    assert four_worker.evidence.states == serial.evidence.states
    assert parallel.evidence.identity == serial.evidence.identity
    assert four_worker.evidence.identity == serial.evidence.identity
    assert parallel.backend_transition_evidence is not None
    assert four_worker.backend_transition_evidence is not None
    assert serial.backend_transition_evidence is not None
    serial_masks = tuple(
        transition.accepted
        for transition in serial.backend_transition_evidence.transitions
    )
    assert (
        tuple(
            transition.accepted
            for transition in parallel.backend_transition_evidence.transitions
        )
        == serial_masks
    )
    assert (
        tuple(
            transition.accepted
            for transition in four_worker.backend_transition_evidence.transitions
        )
        == serial_masks
    )
    assert parallel.raw_capture is not None
    assert four_worker.raw_capture is not None
    assert serial.raw_capture is not None
    assert parallel.raw_capture.objective_request_count == (
        serial.raw_capture.objective_request_count
    )
    assert parallel.raw_capture.evaluation_request_count == (
        serial.raw_capture.evaluation_request_count
    )
    assert four_worker.raw_capture.objective_request_count == (
        serial.raw_capture.objective_request_count
    )
    assert four_worker.raw_capture.evaluation_request_count == (
        serial.raw_capture.evaluation_request_count
    )
    parallel_diagnostics = derive_mcmc_diagnostics(parallel.evidence)
    four_worker_diagnostics = derive_mcmc_diagnostics(four_worker.evidence)
    serial_diagnostics = derive_mcmc_diagnostics(serial.evidence)
    assert parallel_diagnostics.accepted_counts == serial_diagnostics.accepted_counts
    assert four_worker_diagnostics.accepted_counts == (
        serial_diagnostics.accepted_counts
    )
    assert parallel_diagnostics.acceptance_fractions == (
        serial_diagnostics.acceptance_fractions
    )
    assert four_worker_diagnostics.acceptance_fractions == (
        serial_diagnostics.acceptance_fractions
    )
    assert parallel_diagnostics.mean_acceptance_fraction == (
        serial_diagnostics.mean_acceptance_fraction
    )
    assert four_worker_diagnostics.mean_acceptance_fraction == (
        serial_diagnostics.mean_acceptance_fraction
    )


def test_parallel_worker_failure_enters_typed_failed_lifecycle() -> None:
    context = multiprocessing.get_context("spawn")
    with context.Manager() as manager:
        evaluation_observations = manager.list()
        accepted, plan = _plan_context(
            evaluation_observer=evaluation_observations,
            fail_outside_pid=os.getpid(),
        )

        operation = execute_mcmc_evidence(
            accepted,
            plan,
            execution=ExecutionSettings(workers=2),
        )
        worker_pids = {pid for pid, _evaluator, _environment in evaluation_observations}

    assert operation.terminal is McmcOperationTerminal.FAILED
    assert operation.failure_category == "McmcExecutionError"
    assert "MCMC worker failed: McmcExecutionError" in (operation.failure_message)
    assert "kernel_exception: worker kernel failure" in operation.failure_message
    assert operation.raw_capture is not None
    assert operation.raw_capture.complete_state_count == 0
    assert operation.raw_capture.objective_request_count == plan.policy.walkers
    assert operation.raw_capture.evaluation_request_count == plan.policy.walkers
    assert operation.evidence is None
    assert worker_pids
    assert worker_pids.isdisjoint(
        child.pid for child in multiprocessing.active_children()
    )


def test_parallel_worker_initialization_failure_is_typed_and_does_not_hang(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan_context()
    original = native_mcmc._McmcWorkerContext.from_plan

    def invalid_worker_context(source_plan: McmcPlan) -> Any:
        context = original(source_plan)
        return dataclasses.replace(
            context,
            parameterization=dataclasses.replace(
                context.parameterization,
                identity="foreign-parameterization",
            ),
        )

    monkeypatch.setattr(
        native_mcmc._McmcWorkerContext,
        "from_plan",
        staticmethod(invalid_worker_context),
    )
    children_before = {child.pid for child in multiprocessing.active_children()}

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        execution=ExecutionSettings(workers=2),
    )

    assert operation.terminal is McmcOperationTerminal.FAILED
    assert operation.failure_category == "McmcExecutionError"
    assert "MCMC worker failed: McmcExecutionError" in operation.failure_message
    assert "foreign parameterization" in operation.failure_message
    assert operation.raw_capture is not None
    assert operation.raw_capture.complete_state_count == 0
    assert operation.raw_capture.objective_request_count == plan.policy.walkers
    assert operation.raw_capture.evaluation_request_count == 0
    assert {child.pid for child in multiprocessing.active_children()} == children_before


def test_parallel_keyboard_interrupt_terminates_and_joins_workers() -> None:
    context = multiprocessing.get_context("spawn")
    with context.Manager() as manager:
        evaluation_observations = manager.list()
        accepted, plan = _plan_context(
            evaluation_observer=evaluation_observations,
        )

        def interrupt_after_first_transition(state: EnsembleState) -> None:
            if state.ordinal == 1:
                raise KeyboardInterrupt

        operation = execute_mcmc_evidence(
            accepted,
            plan,
            state_observer=interrupt_after_first_transition,
            execution=ExecutionSettings(workers=2),
        )
        worker_pids = {pid for pid, _evaluator, _environment in evaluation_observations}

    assert operation.terminal is McmcOperationTerminal.INTERRUPTED
    assert operation.complete_state_count == 2
    assert worker_pids
    assert worker_pids.isdisjoint(
        child.pid for child in multiprocessing.active_children()
    )


def test_expert_policy_only_overrides_predeclared_topology() -> None:
    policy = ExpertMcmcPolicy(
        burn_steps=7,
        retained_steps=13,
        walkers=8,
        expert_provenance="user-declared-expert-settings",
    )

    resolved = policy.resolve(dimension=2, root_seed=99)

    assert resolved.kind is McmcPolicyKind.EXPERT
    assert resolved.policy_version == "expert-v1"
    assert resolved.provenance_identity == "user-declared-expert-settings"
    assert resolved.burn_steps == 7
    assert resolved.retained_steps == 13
    assert resolved.walkers == 8
    assert resolved.initialization is InitializationKind.BOUNDED_LATIN_HYPERCUBE
    assert resolved.proposal is ProposalKind.STRETCH
    assert resolved.thin == 1
    assert resolved.objective_request_budget == 8 * 21


def test_initial_ensemble_is_a_seeded_bounded_latin_hypercube() -> None:
    lower = (-2.0, 10.0)
    upper = (2.0, 22.0)

    first = build_bounded_latin_hypercube(
        lower,
        upper,
        walkers=8,
        seed=0xCAFE,
    )
    second = build_bounded_latin_hypercube(
        lower,
        upper,
        walkers=8,
        seed=0xCAFE,
    )
    positions = np.asarray(first)

    assert first == second
    assert positions.shape == (8, 2)
    assert np.all(positions > np.asarray(lower))
    assert np.all(positions < np.asarray(upper))
    strata = np.floor(
        (positions - np.asarray(lower)) / (np.asarray(upper) - np.asarray(lower)) * 8
    ).astype(int)
    for dimension in range(2):
        assert sorted(strata[:, dimension]) == list(range(8))


@pytest.mark.parametrize(
    ("lower", "upper"),
    (
        ((-np.inf,), (1.0,)),
        ((0.0,), (np.inf,)),
        ((1.0,), (1.0,)),
        ((2.0,), (1.0,)),
    ),
)
def test_initial_ensemble_rejects_nonfinite_or_empty_bounds(
    lower: tuple[float, ...],
    upper: tuple[float, ...],
) -> None:
    with pytest.raises(McmcConstructionError, match="finite strictly ordered bounds"):
        build_bounded_latin_hypercube(lower, upper, walkers=4, seed=1)


@pytest.mark.parametrize("seed", (True, -1, 1 << 64, 1.5, np.nan))
def test_initial_ensemble_rejects_noncanonical_seed(seed: object) -> None:
    with pytest.raises(McmcConstructionError, match="unsigned 64-bit seed"):
        build_bounded_latin_hypercube((0.0,), (1.0,), walkers=4, seed=cast("int", seed))


def test_initial_ensemble_rejects_transposed_or_unrepresentably_narrow_bounds() -> None:
    with pytest.raises(McmcConstructionError, match="finite strictly ordered bounds"):
        build_bounded_latin_hypercube(
            cast("tuple[float, ...]", ((0.0, 1.0),)),
            cast("tuple[float, ...]", ((2.0, 3.0),)),
            walkers=4,
            seed=1,
        )
    lower = 1.0
    upper = float(np.nextafter(lower, np.inf))
    with pytest.raises(McmcConstructionError, match="closed bound"):
        build_bounded_latin_hypercube((lower,), (upper,), walkers=4, seed=1)


def test_product_initial_ensemble_is_seeded_accepted_point_jitter() -> None:
    first = build_accepted_point_ensemble(
        (1.0, 2.0),
        (0.0, 1.0),
        (2.0, 3.0),
        walkers=8,
        seed=612,
    )
    second = build_accepted_point_ensemble(
        (1.0, 2.0),
        (0.0, 1.0),
        (2.0, 3.0),
        walkers=8,
        seed=612,
    )
    positions = np.asarray(first)

    assert first == second
    assert np.all(positions > np.asarray((0.0, 1.0)))
    assert np.all(positions < np.asarray((2.0, 3.0)))
    np.testing.assert_allclose(
        positions.mean(axis=0),
        (1.0, 2.0),
        atol=2.0e-4,
    )


def test_product_policy_initializes_from_exact_accepted_fit() -> None:
    accepted, problem, parameterization, engine = _native_context()
    policy = resolve_product_mcmc_policy(
        dimension=2,
        walkers=8,
        steps=3,
        root_seed=612,
    )
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    assert policy.kind is McmcPolicyKind.PRODUCT
    assert policy.initialization is InitializationKind.ACCEPTED_POINT_JITTER
    np.testing.assert_allclose(
        np.asarray(plan.initial_ensemble).mean(axis=0),
        accepted.vector,
        atol=2.0e-4,
    )


def test_product_analysis_result_preserves_explicit_burn_and_summary_policy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, problem, parameterization, engine = _native_context()
    request = McmcRequest(
        steps=5,
        burn=1,
        thin=2,
        walkers=8,
        seed=1234,
    )
    policy = resolve_product_mcmc_policy(
        dimension=2,
        walkers=8,
        steps=request.steps,
        root_seed=1234,
    )
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    monkeypatch.setattr(
        native_mcmc.emcee.autocorr,
        "integrated_time",
        lambda _chain, **_kwargs: np.array([1.25, 1.5]),
    )

    parameter_model = _product_parameter_model()
    result = derive_mcmc_analysis_result(
        operation.evidence,
        request,
        parameter_model,
    )

    assert result.status is McmcAnalysisStatus.COMPLETE
    assert result.failure is None
    assert result.discarded_steps == 1
    assert result.chain.shape == (2, 8, 2)
    expected_raw = np.asarray(
        [state.positions for state in operation.evidence.states[1:]],
        dtype=float,
    )
    np.testing.assert_array_equal(result.chain, expected_raw[1::2])
    assert result.summary[0].parameter_id == "A"
    expected_quantiles = np.percentile(
        result.samples[:, 0],
        [2.5, 15.87, 50.0, 84.13, 97.5],
    )
    np.testing.assert_array_equal(
        np.array(
            [
                result.summary[0].eti_95_lower,
                result.summary[0].credible_interval_68_lower,
                result.summary[0].median,
                result.summary[0].credible_interval_68_upper,
                result.summary[0].eti_95_upper,
            ]
        ),
        expected_quantiles,
    )


def test_product_analysis_result_classifies_automatic_burn_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, problem, parameterization, engine = _native_context()
    request = McmcRequest(
        steps=5,
        walkers=8,
        seed=1234,
    )
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=resolve_product_mcmc_policy(
            dimension=2,
            walkers=8,
            steps=request.steps,
            root_seed=1234,
        ),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    monkeypatch.setattr(
        native_mcmc.emcee.autocorr,
        "integrated_time",
        lambda _chain, **_kwargs: object(),
    )

    parameter_model = _product_parameter_model()
    result = derive_mcmc_analysis_result(
        operation.evidence,
        request,
        parameter_model,
    )

    assert result.status is McmcAnalysisStatus.INCOMPLETE
    assert result.failure is not None
    assert (
        result.failure.category
        is McmcAnalysisFailureCategory.AUTOMATIC_BURN_UNAVAILABLE
    )
    assert "autocorrelation time invalid" in result.failure.message
    assert result.failure.preserve_raw_evidence
    assert result.summary == ()
    assert result.raw_chain is not None
    assert result.raw_chain.shape == (5, 8, 2)


def test_product_auto_burn_uses_realistic_tentative_autocorrelation_time() -> None:
    policy = native_mcmc._ProductMcmcInterpretationPolicy(
        requested_steps=10_000,
        burn="auto",
        thin=1,
    )
    tentative_tau = np.array(
        [
            212.787,
            193.166,
            215.584,
            224.983,
            195.036,
            201.503,
            207.812,
            230.222,
            210.754,
            247.529,
            222.710,
            224.156,
            201.567,
            228.543,
            220.933,
            215.228,
            239.719,
        ]
    )

    discarded_steps, failure = native_mcmc._resolve_product_mcmc_burn(
        policy,
        raw_step_count=10_000,
        coordinate_count=17,
        autocorrelation_time=tentative_tau,
    )

    assert discarded_steps == 496
    assert failure is None


def test_product_auto_burn_publishes_tentative_posterior_without_tau_metrics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, problem, parameterization, engine = _native_context()
    request = McmcRequest(steps=5, walkers=8, seed=1234)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=resolve_product_mcmc_policy(
            dimension=2,
            walkers=8,
            steps=request.steps,
            root_seed=1234,
        ),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    tentative_tau = np.array([0.12, 0.124])
    monkeypatch.setattr(
        native_mcmc.emcee.autocorr,
        "integrated_time",
        lambda _chain, **_kwargs: (_ for _ in ()).throw(
            native_mcmc.emcee.autocorr.AutocorrError(tentative_tau)
        ),
    )

    result = derive_mcmc_analysis_result(
        operation.evidence,
        request,
        _product_parameter_model(),
    )

    assert result.status is McmcAnalysisStatus.COMPLETE
    assert result.failure is None
    assert result.discarded_steps == 1
    assert result.retained_step_count == 4
    assert result.retained_sample_count == 32
    assert result.burn_in_warning == (
        "autocorrelation time estimate is unreliable; tentative automatic "
        "burn-in was applied"
    )
    assert result.autocorrelation_report is not None
    assert (
        result.autocorrelation_report.status
        is McmcAutocorrelationStatus.UNRELIABLE_SHORT_CHAIN
    )
    assert result.autocorrelation_report.values == pytest.approx(tentative_tau)
    assert result.autocorrelation_report.sampled_steps_over_maximum == pytest.approx(
        5.0 / 0.124
    )
    assert result.autocorrelation_report.minimum_effective_sample_size is None
    assert all(
        summary.effective_sample_size is None and summary.mcse_mean is None
        for summary in result.summary
    )


@pytest.mark.parametrize(
    ("autocorrelation_outcome", "expected_status"),
    (
        (np.array([1.25, 1.5]), McmcAutocorrelationStatus.RELIABLE),
        (
            native_mcmc.emcee.autocorr.AutocorrError(np.array([1.6, 1.8])),
            McmcAutocorrelationStatus.UNRELIABLE_SHORT_CHAIN,
        ),
    ),
    ids=("reliable", "tentative"),
)
def test_publication_consumes_predetermined_authoritative_mcmc_result(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    autocorrelation_outcome: object,
    expected_status: McmcAutocorrelationStatus,
) -> None:
    accepted, problem, parameterization, engine = _native_context()
    request = McmcRequest(
        steps=5,
        burn=1,
        thin=2,
        walkers=8,
        seed=1234,
    )
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=resolve_product_mcmc_policy(
            dimension=2,
            walkers=8,
            steps=request.steps,
            root_seed=1234,
        ),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None

    def integrated_time(_chain: Array, **_kwargs: object) -> object:
        if isinstance(autocorrelation_outcome, BaseException):
            raise autocorrelation_outcome
        return autocorrelation_outcome

    monkeypatch.setattr(
        native_mcmc.emcee.autocorr,
        "integrated_time",
        integrated_time,
    )
    parameter_model = _product_parameter_model()
    result = derive_mcmc_analysis_result(
        operation.evidence,
        request,
        parameter_model,
    )
    assert result.autocorrelation_report is not None
    assert result.autocorrelation_report.status is expected_status

    write_mcmc_outputs(
        result,
        EffectiveMcmcSettings(
            steps=5,
            burn=1,
            thin=2,
            walkers=8,
            seed=1234,
            workers=1,
            native_threads=None,
            update_parameters=False,
        ),
        tmp_path,
        parameter_model,
    )

    statistics = tmp_path / "Statistics" / "MCMC"
    diagnostics = (statistics / "diagnostics.toml").read_text(encoding="utf-8")
    assert f'autocorrelation_status = "{expected_status.value}"' in diagnostics
    assert "discarded_steps = 1" in diagnostics
    assert (statistics / "summary.toml").is_file()
    assert (statistics / "samples.tsv").is_file()
    assert (statistics / "correlations.tsv").is_file()
    assert (statistics / "plots.pdf").stat().st_size > 0


def test_mcmc_rejects_private_relaxation_coordinates_without_measure() -> None:
    accepted, problem, parameterization, engine = _native_context()
    transformed_problem = dataclasses.replace(
        problem,
        feasible_coordinates=SimpleNamespace(
            is_noop=False,
            identity="private-relaxation-chart",
            supports_box_only_algorithms=False,
            solver_bounds=(problem.lower_bounds, problem.upper_bounds),
        ),
    )
    compatible = AcceptedFitResult.for_qualification(
        occurrence_identity="relaxation-chart-accepted",
        problem_identity=transformed_problem.identity,
        invocation_identity=accepted.invocation_identity,
        execution_identity=accepted.execution_identity,
        materialization_identity=accepted.materialization_identity,
        parameterization_identity=accepted.parameterization_identity,
        evaluator_parameterization_identity=(
            accepted.evaluator_parameterization_identity
        ),
        source_occurrence_identity=accepted.source_occurrence_identity,
        source_revision=accepted.source_revision,
        controlled_ids=accepted.controlled_ids,
        vector=accepted.vector,
        chi_square=accepted.chi_square,
        evaluation_result=accepted.evaluation_result,
        commit_scope=accepted.commit_scope,
        commit_items=accepted.commit_items,
        origin_context_identity=accepted.origin_context_identity,
    )

    with pytest.raises(McmcConstructionError, match="Jacobian measure"):
        McmcPlan.for_accepted(
            compatible,
            source_problem=transformed_problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=_resolved_policy(),
            coordinate_units=(
                ("A", ParameterUnit.DIMENSIONLESS),
                ("B", ParameterUnit.DIMENSIONLESS),
            ),
        )


@pytest.mark.parametrize("walkers", (2, 3))
def test_one_dimensional_product_policy_preserves_two_walker_minimum(
    walkers: int,
) -> None:
    policy = resolve_product_mcmc_policy(
        dimension=1,
        walkers=walkers,
        steps=2,
        root_seed=612,
    )

    assert policy.walkers == walkers


def test_plan_binds_initial_ensemble_to_exact_accepted_native_lineage() -> None:
    accepted, problem, parameterization, engine = _native_context()

    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    assert plan.accepted_result_identity == accepted.identity
    assert plan.anchor.accepted_occurrence_identity == accepted.occurrence_identity
    assert (
        plan.anchor.source_occurrence_identity
        == problem.source_snapshot.occurrence_identity
    )
    assert plan.anchor.source_revision == problem.source_snapshot.revision
    assert plan.anchor.evaluation_plan_identity == engine.plan.identity
    assert plan.anchor.parameterization_identity == parameterization.identity
    assert plan.anchor.held_items == problem.held_items
    assert (
        plan.anchor.affine_feasibility_identity == problem.affine_feasibility_identity
    )
    assert accepted.vector != problem.start
    assert plan.problem_identity == problem.identity
    assert plan.coordinate_ids == ("A", "B")
    assert plan.lower_bounds == (0.0, 1.0)
    assert plan.upper_bounds == (2.0, 3.0)
    assert len(plan.initial_ensemble) == plan.policy.walkers
    assert all(
        lower < value < upper
        for row in plan.initial_ensemble
        for value, lower, upper in zip(
            row,
            plan.lower_bounds,
            plan.upper_bounds,
            strict=True,
        )
    )
    assert plan.identity


def test_plan_rejects_foreign_problem_parameterization_and_coordinate_order() -> None:
    accepted, problem, parameterization, engine = _native_context()
    policy = _resolved_policy()
    units = (
        ("A", ParameterUnit.DIMENSIONLESS),
        ("B", ParameterUnit.DIMENSIONLESS),
    )
    foreign_problem = dataclasses.replace(problem, start=(0.75, 1.5))
    foreign_parameterization = dataclasses.replace(parameterization, source_revision=5)

    with pytest.raises(McmcConstructionError, match="exact accepted native lineage"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=foreign_problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=units,
        )
    with pytest.raises(McmcConstructionError, match="exact accepted native lineage"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=problem,
            parameterization=foreign_parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=units,
        )
    with pytest.raises(McmcConstructionError, match="coordinate units"):
        McmcPlan.for_accepted(
            accepted,
            source_problem=problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=policy,
            coordinate_units=tuple(reversed(units)),
        )


def test_complete_chain_is_primary_and_flat_samples_are_derived_views() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    source_vector = accepted.vector
    source_items = tuple(problem.source_snapshot.items())

    operation = execute_mcmc_evidence(accepted, plan)
    replay = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    evidence = operation.evidence
    assert replay.evidence is not None
    assert replay.evidence.identity == evidence.identity
    assert evidence.lifecycle is McmcEvidenceLifecycle.COMPLETED
    assert evidence.plan_identity == plan.identity
    assert evidence.accepted_result_identity == accepted.identity
    assert evidence.completed_transition_count == plan.policy.total_steps
    assert tuple(state.ordinal for state in evidence.states) == tuple(
        range(plan.policy.total_steps + 1)
    )
    assert evidence.states[0].positions == plan.initial_ensemble
    assert all(len(state.positions) == plan.policy.walkers for state in evidence.states)
    assert all(
        len(position) == plan.policy.dimension
        for state in evidence.states
        for position in state.positions
    )
    assert all(
        len(state.log_densities) == plan.policy.walkers for state in evidence.states
    )
    assert evidence.objective_request_count == plan.policy.objective_request_budget
    assert accepted.vector == source_vector
    assert tuple(problem.source_snapshot.items()) == source_items

    initial = np.asarray(evidence.states[0].positions)
    metadata = np.asarray((0.0, 1.0, 2.0))
    observed = np.asarray((1.0, 3.0, 5.0))
    expected_initial_log_density = -0.5 * np.sum(
        (initial[:, 0, np.newaxis] + initial[:, 1, np.newaxis] * metadata - observed)
        ** 2,
        axis=1,
    )
    # These small, well-conditioned binary64 calculations use only affine
    # arithmetic and a three-term sum. 1e-14 is roughly 45 machine eps at unit
    # scale: enough for operation-order differences while still detecting a
    # changed seeded trajectory or log-density convention.
    np.testing.assert_allclose(
        evidence.states[0].log_densities,
        expected_initial_log_density,
        rtol=0.0,
        atol=1.0e-14,
    )
    # A local seeded replay is ordinary deterministic evidence, not an
    # independently frozen product reference.
    assert replay.evidence.states == evidence.states
    retained = derive_retained_sample_view(evidence)
    assert retained.source_evidence_identity == evidence.identity
    assert retained.coordinate_ids == plan.coordinate_ids
    assert retained.selected_state_ordinals == (3, 4, 5)
    assert retained.sample_indices[0] == (3, 0)
    assert retained.sample_indices[-1] == (5, 7)
    assert len(retained.samples) == plan.policy.walkers * plan.policy.retained_steps
    assert len(retained.log_densities) == len(retained.samples)
    assert retained.is_complete

    diagnostics = derive_mcmc_diagnostics(evidence)
    assert operation.backend_transition_evidence is not None
    assert diagnostics.source_identity == operation.backend_transition_evidence.identity
    assert diagnostics.state_ordinals == (1, 2, 3, 4, 5)
    assert diagnostics.walker_ordinals == tuple(range(plan.policy.walkers))
    assert diagnostics.acceptance_fractions is not None
    assert diagnostics.mean_acceptance_fraction is not None
    assert len(diagnostics.acceptance_fractions) == plan.policy.walkers
    assert 0.0 <= diagnostics.mean_acceptance_fraction <= 1.0
    assert diagnostics.accepted_counts == (0, 5, 2, 4, 3, 2, 4, 5)
    assert diagnostics.acceptance_fractions == (
        0.0,
        1.0,
        0.4,
        0.8,
        0.6,
        0.4,
        0.8,
        1.0,
    )
    assert diagnostics.mean_acceptance_fraction == 0.625


def test_seeded_native_chain_recovers_known_truncated_normal_posterior() -> None:
    accepted, problem, parameterization, engine = _native_context(fit_b=False)
    policy = resolve_product_mcmc_policy(
        dimension=1,
        walkers=8,
        steps=6000,
        root_seed=657,
    )
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(("A", ParameterUnit.DIMENSIONLESS),),
    )

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    samples = np.asarray(
        [
            position[0]
            for state in operation.evidence.states[1001:]
            for position in state.positions
        ]
    )
    # With B fixed at 1.5, the three unit-error residuals are
    # A - (1.0, 1.5, 2.0), hence N(1.5, 1/3) truncated to [0, 2].
    mean = 1.5
    standard_deviation = 1.0 / np.sqrt(3.0)
    distribution = truncnorm(
        (0.0 - mean) / standard_deviation,
        (2.0 - mean) / standard_deviation,
        loc=mean,
        scale=standard_deviation,
    )

    assert samples.mean() == pytest.approx(distribution.mean(), abs=0.04)
    assert samples.std(ddof=1) == pytest.approx(distribution.std(), abs=0.025)


def test_interruption_preserves_unqualified_complete_state_prefix_without_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    def interrupt_after_second_transition(state: EnsembleState) -> None:
        if state.ordinal == 2:
            raise KeyboardInterrupt

    def forbid_capture_qualification(*_args: object, **_kwargs: object) -> None:
        pytest.fail("interrupted execution must not qualify an incomplete capture")

    monkeypatch.setattr(
        native_mcmc,
        "validate_raw_mcmc_capture",
        forbid_capture_qualification,
    )

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        state_observer=interrupt_after_second_transition,
    )

    assert operation.terminal is McmcOperationTerminal.INTERRUPTED
    assert operation.failure_category == "interrupted"
    assert operation.evidence is None
    assert operation.validation is None
    assert operation.raw_capture is not None
    assert tuple(state.ordinal for state in operation.raw_capture.states) == (0, 1, 2)
    assert operation.raw_capture.objective_request_count == plan.policy.walkers * 3


def test_interruption_after_qualification_preserves_complete_chain_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan_context()
    original = native_mcmc.validate_raw_mcmc_capture
    validation_calls = 0

    def count_validation(*args: object, **kwargs: object):
        nonlocal validation_calls
        validation_calls += 1
        return original(*args, **kwargs)

    def interrupt_before_final_assembly(
        stage: McmcExecutionStage,
        _count: int,
    ) -> None:
        if stage is McmcExecutionStage.BEFORE_FINAL_ASSEMBLY:
            raise KeyboardInterrupt

    monkeypatch.setattr(
        native_mcmc,
        "validate_raw_mcmc_capture",
        count_validation,
    )
    operation = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=interrupt_before_final_assembly,
    )

    assert validation_calls == 1
    assert operation.terminal is McmcOperationTerminal.INTERRUPTED
    assert operation.evidence is not None
    assert operation.evidence.lifecycle is McmcEvidenceLifecycle.COMPLETED
    assert operation.validation is not None
    assert operation.validation.is_complete
    assert operation.raw_capture is not None
    assert operation.raw_capture.terminal is McmcOperationTerminal.COMPLETED


def test_primary_chain_evidence_round_trips_and_rejects_tampering() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    evidence = operation.evidence

    record = evidence.to_record()
    assert operation.raw_capture is not None
    restored = McmcEvidence.from_record(
        record,
        plan=plan,
        raw_capture=operation.raw_capture,
    )

    assert restored.identity == evidence.identity
    assert restored.states == evidence.states
    assert restored.objective_request_count == evidence.objective_request_count

    tampered = deepcopy(record)
    states = cast("list[dict[str, object]]", tampered["states"])
    positions = cast("list[list[str]]", states[1]["positions"])
    positions[0][0] = float(float.fromhex(positions[0][0]) + 0.1).hex()
    with pytest.raises(McmcConstructionError, match="identity"):
        McmcEvidence.from_record(
            tampered,
            plan=plan,
            raw_capture=operation.raw_capture,
        )


def test_execution_rejects_same_semantics_from_foreign_accepted_occurrence() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    foreign = _same_semantics_foreign_occurrence(accepted)
    assert foreign.identity == accepted.identity
    assert foreign.occurrence_identity != accepted.occurrence_identity

    with pytest.raises(McmcConstructionError, match="accepted occurrence"):
        execute_mcmc_evidence(foreign, plan)


def test_capture_integrity_rejects_tampered_captured_log_density() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    first = operation.raw_capture.states[0]
    tampered_first = dataclasses.replace(
        first,
        log_densities=(first.log_densities[0] + 1.0, *first.log_densities[1:]),
    )
    with pytest.raises(McmcConstructionError, match="capture"):
        dataclasses.replace(
            operation.raw_capture,
            states=(tampered_first, *operation.raw_capture.states[1:]),
        )


def _copy_raw_capture(
    source: native_mcmc.RawMcmcCapture,
    **changes: Any,
) -> native_mcmc.RawMcmcCapture:
    values = {
        "plan_identity": source.plan_identity,
        "policy_identity": source.policy_identity,
        "root_seed": source.root_seed,
        "walkers": source.walkers,
        "dimension": source.dimension,
        "total_steps": source.total_steps,
        "backend_execution_occurrence_identity": (
            source.backend_execution_occurrence_identity
        ),
        "terminal": source.terminal,
        "states": source.states,
        "objective_request_count": source.objective_request_count,
        "evaluation_request_count": source.evaluation_request_count,
        "stage": source.stage,
        "initialization_outcome": source.initialization_outcome,
        "failure_category": source.failure_category,
    }
    values.update(changes)
    return native_mcmc.RawMcmcCapture(
        **values,
        _occurrence_witness=native_mcmc._mint_mcmc_evidence_witness("raw-capture"),
    )


@pytest.mark.parametrize("field", ["walkers", "dimension", "total_steps"])
def test_capture_qualification_rejects_frozen_topology_metadata_mutation(
    field: str,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    source = operation.raw_capture
    corrupted = _copy_raw_capture(
        source,
        **{field: cast("int", getattr(source, field)) + 1},
    )

    validation = validate_raw_mcmc_capture(plan, corrupted)

    assert not validation.is_complete
    assert validation.primary_evidence is None
    assert validation.failures[0].category == "capture_topology_mismatch"


@pytest.mark.parametrize("accounting", ["objective", "evaluation"])
def test_capture_qualification_rejects_request_accounting_mutation(
    accounting: str,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    source = operation.raw_capture
    changes = (
        {"objective_request_count": source.objective_request_count - 1}
        if accounting == "objective"
        else {"evaluation_request_count": source.objective_request_count + 1}
    )
    corrupted = _copy_raw_capture(source, **changes)

    validation = validate_raw_mcmc_capture(plan, corrupted)

    assert not validation.is_complete
    assert validation.primary_evidence is None
    assert validation.failures[0].category == "request_accounting_mismatch"


@pytest.mark.parametrize("corruption", ["bounds", "dimension"])
def test_capture_qualification_rejects_validly_encoded_state_corruption(
    corruption: str,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    source = operation.raw_capture
    state = source.states[1]
    positions = [list(row) for row in state.positions]
    if corruption == "bounds":
        positions[0][0] = plan.upper_bounds[0]
        expected_category = "state_outside_frozen_bounds"
    else:
        positions = [row[:-1] for row in positions]
        expected_category = "state_topology_mismatch"
    corrupted_state = EnsembleState(
        state.ordinal,
        tuple(tuple(row) for row in positions),
        state.log_densities,
    )
    corrupted = _copy_raw_capture(
        source,
        states=(source.states[0], corrupted_state, *source.states[2:]),
    )

    validation = validate_raw_mcmc_capture(plan, corrupted)

    assert not validation.is_complete
    assert validation.primary_evidence is None
    assert validation.failures[0].category == expected_category


def test_capture_qualification_rejects_missing_completed_state() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    corrupted = _copy_raw_capture(
        operation.raw_capture,
        states=operation.raw_capture.states[:-1],
    )

    validation = validate_raw_mcmc_capture(plan, corrupted)

    assert not validation.is_complete
    assert validation.primary_evidence is None
    assert validation.failures[0].category == "incomplete_completed_capture"


@pytest.mark.parametrize("corruption", ["reordered", "invalid_ordinal"])
def test_raw_capture_construction_rejects_noncanonical_state_order(
    corruption: str,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    states = list(operation.raw_capture.states)
    if corruption == "reordered":
        states[1], states[2] = states[2], states[1]
    else:
        state = states[1]
        states[1] = EnsembleState(
            state.ordinal + 1,
            state.positions,
            state.log_densities,
        )

    with pytest.raises(McmcConstructionError, match="contiguous"):
        _copy_raw_capture(operation.raw_capture, states=tuple(states))


def test_execution_qualification_does_not_repeat_chain_likelihoods(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan_context()
    original_calculate = _LinearPulseSequence.calculate
    calculation_count = 0

    def count_calculation(
        pulse_sequence: _LinearPulseSequence,
        spectrometer: _LinearSpectrometer,
        data: Data,
    ) -> Array:
        nonlocal calculation_count
        calculation_count += 1
        return original_calculate(pulse_sequence, spectrometer, data)

    monkeypatch.setattr(_LinearPulseSequence, "calculate", count_calculation)

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.raw_capture is not None
    assert calculation_count == operation.raw_capture.evaluation_request_count
    assert operation.evidence is not None
    assert operation.validation is not None
    assert operation.validation.is_complete


def _plan_context(
    *,
    evaluation_observer: Any | None = None,
    evaluation_barrier: Any | None = None,
    fail_outside_pid: int | None = None,
) -> tuple[AcceptedFitResult, McmcPlan]:
    accepted, problem, parameterization, engine = _native_context(
        evaluation_observer=evaluation_observer,
        evaluation_barrier=evaluation_barrier,
        fail_outside_pid=fail_outside_pid,
    )
    return accepted, McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )


def test_cancellation_before_initialization_freezes_zero_state_operation() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()
    cancellation.cancel()

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
    )

    assert operation.terminal is McmcOperationTerminal.CANCELLED
    assert operation.complete_state_count == 0
    assert operation.stage is McmcExecutionStage.BEFORE_INITIALIZATION
    assert operation.evidence is None
    assert operation.raw_capture is not None
    assert operation.raw_capture.policy_identity == plan.policy.identity
    assert operation.raw_capture.root_seed == plan.policy.root_seed
    assert (
        operation.raw_capture.initialization_outcome
        is McmcInitializationOutcome.CANCELLED
    )
    diagnostics = derive_mcmc_operation_diagnostics(operation)
    assert diagnostics.status is McmcDiagnosticStatus.UNAVAILABLE
    assert diagnostics.reason is McmcDiagnosticReason.NO_TRANSITIONS


@pytest.mark.parametrize(
    ("target_stage", "terminal", "category"),
    (
        (
            McmcExecutionStage.INITIALIZING,
            McmcOperationTerminal.INTERRUPTED,
            "interrupted",
        ),
        (
            McmcExecutionStage.BEFORE_TRANSITION,
            McmcOperationTerminal.CANCELLED,
            "cancelled",
        ),
        (
            McmcExecutionStage.DURING_TRANSITION,
            McmcOperationTerminal.INTERRUPTED,
            "interrupted",
        ),
    ),
)
def test_zero_or_initial_state_terminal_checkpoints_are_truthful(
    target_stage: McmcExecutionStage,
    terminal: McmcOperationTerminal,
    category: str,
) -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def stop(stage: McmcExecutionStage, _count: int) -> None:
        if stage is not target_stage:
            return
        if terminal is McmcOperationTerminal.CANCELLED:
            cancellation.cancel()
        else:
            raise KeyboardInterrupt

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=stop,
    )

    assert operation.terminal is terminal
    assert operation.failure_category == category
    assert operation.stage is target_stage
    expected_count = 0 if target_stage is McmcExecutionStage.INITIALIZING else 1
    assert operation.complete_state_count == expected_count


def test_initialization_failure_and_later_failure_preserve_only_complete_prefix() -> (
    None
):
    accepted, plan = _plan_context()

    def fail_initialization(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.INITIALIZING:
            raise RuntimeError("initialization failure")

    initial_failure = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=fail_initialization,
    )
    assert initial_failure.terminal is McmcOperationTerminal.FAILED
    assert initial_failure.complete_state_count == 0
    assert initial_failure.evidence is None
    assert initial_failure.raw_capture is not None
    assert (
        initial_failure.raw_capture.initialization_outcome
        is McmcInitializationOutcome.FAILED
    )

    def fail_after_two_transitions(stage: McmcExecutionStage, count: int) -> None:
        if stage is McmcExecutionStage.AFTER_COMPLETE_STATE and count == 3:
            raise RuntimeError("later failure")

    later_failure = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=fail_after_two_transitions,
    )
    assert later_failure.terminal is McmcOperationTerminal.FAILED
    assert later_failure.complete_state_count == 3
    assert later_failure.evidence is not None
    assert tuple(state.ordinal for state in later_failure.evidence.states) == (0, 1, 2)


def test_cancellation_after_complete_state_and_before_final_assembly() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def cancel_after_one_transition(stage: McmcExecutionStage, count: int) -> None:
        if stage is McmcExecutionStage.AFTER_COMPLETE_STATE and count == 2:
            cancellation.cancel()

    partial = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=cancel_after_one_transition,
    )
    assert partial.terminal is McmcOperationTerminal.CANCELLED
    assert partial.complete_state_count == 2
    assert partial.evidence is not None
    assert tuple(state.ordinal for state in partial.evidence.states) == (0, 1)

    final_cancellation = CancellationToken()

    def cancel_final(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.BEFORE_FINAL_ASSEMBLY:
            final_cancellation.cancel()

    final = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=final_cancellation,
        checkpoint_observer=cancel_final,
    )
    assert final.terminal is McmcOperationTerminal.CANCELLED
    assert final.complete_state_count == plan.policy.total_steps + 1
    assert final.evidence is not None
    assert final.evidence.lifecycle is McmcEvidenceLifecycle.PARTIAL


def test_zero_transition_acceptance_is_typed_unavailable() -> None:
    accepted, plan = _plan_context()
    cancellation = CancellationToken()

    def cancel_before_transition(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.BEFORE_TRANSITION:
            cancellation.cancel()

    operation = execute_mcmc_evidence(
        accepted,
        plan,
        cancellation=cancellation,
        checkpoint_observer=cancel_before_transition,
    )
    assert operation.evidence is not None

    diagnostics = derive_mcmc_diagnostics(operation.evidence)

    assert diagnostics.status is McmcDiagnosticStatus.UNAVAILABLE
    assert diagnostics.reason is McmcDiagnosticReason.NO_TRANSITIONS
    assert diagnostics.observed_transition_count == 0
    assert diagnostics.accepted_counts == (0,) * plan.policy.walkers
    assert diagnostics.acceptance_fractions is None
    assert diagnostics.mean_acceptance_fraction is None

    source = operation.backend_transition_evidence
    assert source is not None
    historical = native_mcmc.BackendTransitionEvidence.from_record(
        source.to_record(),
        source=source,
    )
    historical_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        historical
    )
    assert historical_diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert (
        historical_diagnostic.reason
        is McmcDiagnosticReason.HISTORICAL_OBSERVATION_UNVERIFIED
    )
    assert historical_diagnostic.acceptance_fractions is None
    assert historical_diagnostic.mean_acceptance_fraction is None


@pytest.mark.parametrize(
    ("fractions", "mean"),
    (((2.5,), -4.0), ((float("nan"),), float("nan"))),
)
def test_acceptance_diagnostics_reject_caller_supplied_values(
    fractions: tuple[float, ...],
    mean: float,
) -> None:
    constructor = cast("Any", native_mcmc.McmcDiagnostics)

    with pytest.raises(TypeError):
        constructor(
            "forged-source",
            (1,),
            (0,),
            McmcDiagnosticStatus.AVAILABLE,
            None,
            fractions,
            mean,
        )


def test_acceptance_diagnostics_reject_replacement_and_serialized_tampering() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(source)
    replace = cast("Any", dataclasses.replace)

    with pytest.raises(TypeError, match="init=False"):
        replace(diagnostic, acceptance_fractions=(2.5,))

    object.__setattr__(diagnostic, "mean_acceptance_fraction", -4.0)
    with pytest.raises(McmcConstructionError, match="content"):
        diagnostic.validate_integrity()

    canonical = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(source)
    record = canonical.to_record()
    for field_name, forged_value in (
        ("source_identity", "forged"),
        ("walker_ordinals", list(reversed(canonical.walker_ordinals))),
        ("observed_transition_count", 999),
        ("accepted_counts", [5] * plan.policy.walkers),
        ("acceptance_fractions", [float("nan").hex()]),
        ("mean_acceptance_fraction", (-4.0).hex()),
        ("status", McmcDiagnosticStatus.UNAVAILABLE.value),
        ("identity", "recomputed-outer-hash"),
    ):
        tampered = dict(record)
        tampered[field_name] = forged_value
        with pytest.raises(McmcConstructionError, match="canonical backend evidence"):
            native_mcmc.AcceptanceDiagnostics.from_record(tampered, source=source)


def test_acceptance_diagnostics_reject_forged_source_lineage() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(source)

    object.__setattr__(diagnostic, "source_identity", "forged")

    with pytest.raises(McmcConstructionError, match="source identity"):
        diagnostic.validate_integrity()


@pytest.mark.parametrize(
    ("field_name", "forged_value"),
    (
        ("state_ordinals", (5, 4, 3, 2, 1)),
        ("walker_ordinals", (7, 6, 5, 4, 3, 2, 1, 0)),
        ("observed_transition_count", 999),
        ("accepted_counts", (5,) * 8),
        ("acceptance_fractions", (1.0,) * 8),
        ("mean_acceptance_fraction", 1.0),
        ("status", McmcDiagnosticStatus.UNAVAILABLE),
    ),
)
def test_acceptance_diagnostics_reject_recursive_field_mutation(
    field_name: str,
    forged_value: object,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(source)

    object.__setattr__(diagnostic, field_name, forged_value)
    object.__setattr__(diagnostic, "identity", diagnostic._content_identity())

    with pytest.raises(McmcConstructionError, match="content"):
        diagnostic.validate_integrity()


def test_acceptance_diagnostics_reject_replaced_source_object() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(source)
    hypothetical = source.with_hypothetical_masks(
        tuple(
            (True,) * plan.policy.walkers
            for _transition in range(plan.policy.total_steps)
        ),
        observation_provenance="emcee-stretch-backend-v1",
    )

    object.__setattr__(diagnostic, "source", hypothetical)
    object.__setattr__(diagnostic, "identity", diagnostic._content_identity())

    with pytest.raises(McmcConstructionError, match="source identity"):
        diagnostic.validate_integrity()


@pytest.mark.parametrize(
    "mask_kind",
    ("all_true", "all_false", "alternating"),
)
def test_backend_mask_variants_do_not_change_primary_scientific_evidence(
    mask_kind: str,
) -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.raw_capture is not None
    assert operation.evidence is not None
    backend = operation.backend_transition_evidence
    assert backend is not None
    masks = tuple(
        tuple(
            True
            if mask_kind == "all_true"
            else False
            if mask_kind == "all_false"
            else (transition + walker) % 2 == 0
            for walker in range(plan.policy.walkers)
        )
        for transition in range(plan.policy.total_steps)
    )
    variant = backend.with_hypothetical_masks(
        masks,
        observation_provenance=backend.observation_provenance,
    )

    validation = validate_raw_mcmc_capture(
        plan,
        operation.raw_capture,
    )

    assert validation.primary_evidence is not None
    assert validation.primary_evidence.identity == operation.evidence.identity
    assert validation.primary_evidence.states == operation.evidence.states
    object.__setattr__(
        validation.primary_evidence,
        "backend_transition_evidence",
        backend,
    )
    validation.primary_evidence.validate_integrity()
    original_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        backend
    )
    variant_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        variant
    )
    assert backend.kind is native_mcmc.BackendTransitionEvidenceKind.OBSERVED_EXECUTION
    assert variant.kind is native_mcmc.BackendTransitionEvidenceKind.HYPOTHETICAL
    assert original_diagnostic.status is McmcDiagnosticStatus.AVAILABLE
    assert variant_diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert (
        variant_diagnostic.reason
        is McmcDiagnosticReason.UNVALIDATED_BACKEND_TRANSITIONS
    )
    assert variant_diagnostic.acceptance_fractions is None
    assert variant_diagnostic.mean_acceptance_fraction is None
    assert variant_diagnostic.source_identity != original_diagnostic.source_identity
    assert variant_diagnostic.identity != original_diagnostic.identity

    original_selection = derive_retained_sample_view(operation.evidence)
    variant_selection = derive_retained_sample_view(validation.primary_evidence)
    assert variant_selection.samples == original_selection.samples
    assert variant_selection.log_densities == original_selection.log_densities
    output_units = (
        ("A", ParameterUnit.DIMENSIONLESS),
        ("B", ParameterUnit.DIMENSIONLESS),
    )
    original_posterior = derive_posterior_sample_evidence(
        original_selection,
        output_units,
    )
    variant_posterior = derive_posterior_sample_evidence(
        variant_selection, output_units
    )
    assert variant_posterior.outcomes == original_posterior.outcomes
    original_summary = derive_posterior_summary(original_posterior)
    variant_summary = derive_posterior_summary(variant_posterior)
    assert variant_summary.parameter_summaries == original_summary.parameter_summaries
    assert variant_summary.covariance == original_summary.covariance
    assert variant_summary.correlations == original_summary.correlations


def test_observed_backend_evidence_is_bound_to_its_execution_occurrence() -> None:
    accepted, plan = _plan_context()
    first = execute_mcmc_evidence(accepted, plan)
    second = execute_mcmc_evidence(accepted, plan)
    first_source = first.backend_transition_evidence
    second_source = second.backend_transition_evidence
    assert first_source is not None
    assert second_source is not None
    assert first_source.execution_occurrence_identity != (
        second_source.execution_occurrence_identity
    )

    record = first_source.to_record()
    with pytest.raises(McmcConstructionError, match="execution occurrence"):
        native_mcmc.BackendTransitionEvidence.from_record(
            record,
            source=second_source,
        )

    restored = native_mcmc.BackendTransitionEvidence.from_record(
        record,
        source=first_source,
    )
    assert (
        restored.kind
        is native_mcmc.BackendTransitionEvidenceKind.HISTORICAL_OBSERVATION
    )
    restored_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        restored
    )
    assert restored.transitions == first_source.transitions
    assert restored.observed_mask_payload_identity == (
        first_source.observed_mask_payload_identity
    )
    assert restored_diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert (
        restored_diagnostic.reason
        is McmcDiagnosticReason.HISTORICAL_OBSERVATION_UNVERIFIED
    )
    assert restored_diagnostic.acceptance_fractions is None
    assert restored_diagnostic.mean_acceptance_fraction is None

    assert second.raw_capture is not None
    foreign_masks = tuple(item.accepted for item in first_source.transitions)
    foreign = native_mcmc.BackendTransitionEvidence.from_capture(
        second.raw_capture,
        foreign_masks,
        observation_provenance="emcee-stretch-backend-v1",
    )
    foreign_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        foreign
    )
    assert foreign.kind is native_mcmc.BackendTransitionEvidenceKind.HYPOTHETICAL
    assert foreign_diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert foreign_diagnostic.acceptance_fractions is None

    for field_name, forged_value in (
        ("source_capture_identity", "forged"),
        ("execution_occurrence_identity", "forged"),
        ("observed_mask_payload_identity", "forged"),
        ("masks", [[True] * plan.policy.walkers] * plan.policy.total_steps),
        ("identity", "forged"),
    ):
        tampered = dict(record)
        tampered[field_name] = forged_value
        with pytest.raises(McmcConstructionError):
            native_mcmc.BackendTransitionEvidence.from_record(
                tampered,
                source=first_source,
            )


def test_fabricated_historical_masks_cannot_authorize_acceptance_diagnostics() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    all_true = source.with_hypothetical_masks(
        tuple(
            (True,) * plan.policy.walkers
            for _transition in range(plan.policy.total_steps)
        ),
        observation_provenance=source.observation_provenance,
    )

    fabricated_history = native_mcmc._initialize_backend_transition_evidence(
        native_mcmc.BackendTransitionEvidence,
        source.source_capture,
        native_mcmc.BackendTransitionEvidenceKind.HISTORICAL_OBSERVATION,
        source.observation_provenance,
        source.execution_occurrence_identity,
        all_true.transitions,
        None,
        native_mcmc._mint_mcmc_evidence_witness("backend-transition-evidence"),
    )
    diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        fabricated_history
    )

    assert all(
        all(transition.accepted) for transition in fabricated_history.transitions
    )
    assert diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert diagnostic.reason is McmcDiagnosticReason.HISTORICAL_OBSERVATION_UNVERIFIED
    assert diagnostic.acceptance_fractions is None
    assert diagnostic.mean_acceptance_fraction is None


def test_direct_constructor_cannot_forge_observed_backend_evidence() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    source = operation.backend_transition_evidence
    assert source is not None
    all_true = source.with_hypothetical_masks(
        tuple(
            (True,) * plan.policy.walkers
            for _transition in range(plan.policy.total_steps)
        ),
        observation_provenance=source.observation_provenance,
    )

    with pytest.raises(TypeError, match="canonical factories"):
        native_mcmc.BackendTransitionEvidence(
            source.source_capture,
            native_mcmc.BackendTransitionEvidenceKind.OBSERVED_EXECUTION,
            source.observation_provenance,
            source.execution_occurrence_identity,
            all_true.transitions,
            native_mcmc._mint_mcmc_evidence_witness("backend-transition-evidence"),
        )
    with pytest.raises(McmcConstructionError, match="authority registry"):
        native_mcmc._initialize_backend_transition_evidence(
            native_mcmc.BackendTransitionEvidence,
            source.source_capture,
            native_mcmc.BackendTransitionEvidenceKind.OBSERVED_EXECUTION,
            source.observation_provenance,
            source.execution_occurrence_identity,
            all_true.transitions,
            source._live_observation_authority,
            native_mcmc._mint_mcmc_evidence_witness("backend-transition-evidence"),
        )

    hypothetical_diagnostic = native_mcmc.AcceptanceDiagnostics.from_backend_evidence(
        all_true
    )
    assert all_true.kind is native_mcmc.BackendTransitionEvidenceKind.HYPOTHETICAL
    assert hypothetical_diagnostic.status is McmcDiagnosticStatus.UNAVAILABLE
    assert hypothetical_diagnostic.acceptance_fractions is None


def test_primary_rejects_deterministic_replay_backend_substitution() -> None:
    accepted, plan = _plan_context()
    replay_a = execute_mcmc_evidence(accepted, plan)
    replay_b = execute_mcmc_evidence(accepted, plan)
    evidence_a = replay_a.evidence
    evidence_b = replay_b.evidence
    backend_a = replay_a.backend_transition_evidence
    backend_b = replay_b.backend_transition_evidence
    assert evidence_a is not None
    assert evidence_b is not None
    assert backend_a is not None
    assert backend_b is not None
    assert evidence_a.identity == evidence_b.identity
    assert replay_a.raw_capture is not None
    assert replay_b.raw_capture is not None
    assert replay_a.raw_capture.identity == replay_b.raw_capture.identity
    assert backend_a.execution_occurrence_identity != (
        backend_b.execution_occurrence_identity
    )

    object.__setattr__(evidence_a, "backend_transition_evidence", backend_b)

    with pytest.raises(McmcConstructionError, match="execution occurrence"):
        evidence_a.validate_integrity()
    with pytest.raises(McmcConstructionError, match="execution occurrence"):
        derive_mcmc_diagnostics(evidence_a)


def test_backend_observation_authority_is_exact_nontransferable_capability(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    accepted, plan = _plan_context()
    original_factory = cast(
        "Any",
        native_mcmc._mint_observed_backend_transition_evidence,
    )
    captured: list[Any] = []

    def capture_authority(
        authority: Any,
        raw_capture: Any,
    ) -> Any:
        captured.append(authority)
        return original_factory(authority, raw_capture)

    monkeypatch.setattr(
        native_mcmc,
        "_mint_observed_backend_transition_evidence",
        capture_authority,
    )
    replay_a = execute_mcmc_evidence(accepted, plan)
    replay_b = execute_mcmc_evidence(accepted, plan)
    assert replay_a.raw_capture is not None
    assert replay_b.raw_capture is not None
    authority_a, authority_b = captured

    with pytest.raises(TypeError, match="cannot be copied"):
        copy.copy(authority_a)
    with pytest.raises(TypeError, match="cannot be copied"):
        copy.deepcopy(authority_a)
    with pytest.raises(TypeError, match="cannot be serialized"):
        pickle.dumps(authority_a)
    with pytest.raises(McmcConstructionError, match="foreign execution occurrence"):
        original_factory(authority_b, replay_a.raw_capture)
    with pytest.raises(McmcConstructionError, match="foreign execution occurrence"):
        original_factory(authority_a, replay_b.raw_capture)

    raw_a = replay_a.raw_capture
    copied_raw_a = native_mcmc.RawMcmcCapture(
        raw_a.plan_identity,
        raw_a.policy_identity,
        raw_a.root_seed,
        raw_a.walkers,
        raw_a.dimension,
        raw_a.total_steps,
        raw_a.backend_execution_occurrence_identity,
        raw_a.terminal,
        raw_a.states,
        raw_a.objective_request_count,
        raw_a.evaluation_request_count,
        raw_a.stage,
        raw_a.initialization_outcome,
        raw_a.failure_category,
        native_mcmc._mint_mcmc_evidence_witness("raw-capture"),
    )
    assert copied_raw_a.identity == raw_a.identity
    with pytest.raises(McmcConstructionError, match="foreign execution occurrence"):
        original_factory(authority_a, copied_raw_a)

    absent = object.__new__(native_mcmc._BackendExecutionObservation)
    with pytest.raises(McmcConstructionError, match="foreign execution occurrence"):
        original_factory(absent, replay_a.raw_capture)


def test_short_chain_autocorrelation_and_dependents_are_unreliable() -> None:
    accepted, problem, parameterization, engine = _native_context()
    retained_steps = 16
    policy = ExpertMcmcPolicy(
        burn_steps=0,
        retained_steps=retained_steps,
        walkers=8,
        expert_provenance="short-autocorrelation-qualification",
    ).resolve(dimension=2, root_seed=1234)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)
    posterior = derive_posterior_sample_evidence(selection, plan.coordinate_units)
    summary_policy = native_mcmc.PosteriorSummaryPolicy()
    changed_window_policy = dataclasses.replace(
        summary_policy,
        autocorrelation_window_parameter=4.0,
    )
    changed_tolerance_policy = dataclasses.replace(
        summary_policy,
        autocorrelation_adequacy_tolerance=25.0,
    )
    assert changed_window_policy.identity != summary_policy.identity
    assert changed_tolerance_policy.identity != summary_policy.identity
    summary = derive_posterior_summary(posterior, summary_policy)

    assert summary.policy_identity == summary_policy.identity
    for item in summary.parameter_summaries:
        for estimate in (
            item.autocorrelation_time,
            item.effective_sample_size,
            item.monte_carlo_standard_error,
        ):
            assert estimate.status is McmcDiagnosticStatus.UNAVAILABLE
            assert estimate.reason is McmcDiagnosticReason.UNRELIABLE_AUTOCORRELATION
            assert estimate.value is None


def test_long_chain_autocorrelation_and_dependents_are_reliable() -> None:
    accepted, problem, parameterization, engine = _native_context()
    retained_steps = 2048
    walkers = 8
    policy = ExpertMcmcPolicy(
        burn_steps=256,
        retained_steps=retained_steps,
        walkers=walkers,
        expert_provenance="long-autocorrelation-qualification",
    ).resolve(dimension=2, root_seed=5678)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)
    posterior = derive_posterior_sample_evidence(selection, plan.coordinate_units)
    summary_policy = native_mcmc.PosteriorSummaryPolicy()
    summary = derive_posterior_summary(posterior, summary_policy)

    assert summary.source_identity == posterior.identity
    assert summary.policy_identity == summary_policy.identity
    assert summary.included_labels == selection.sample_indices
    sample_count = retained_steps * walkers
    for item in summary.parameter_summaries:
        tau = item.autocorrelation_time
        effective = item.effective_sample_size
        mcse = item.monte_carlo_standard_error
        assert tau.status is McmcDiagnosticStatus.AVAILABLE
        assert tau.reason is None
        assert tau.value is not None
        assert retained_steps >= (
            summary_policy.autocorrelation_adequacy_tolerance * tau.value
        )
        assert effective.status is McmcDiagnosticStatus.AVAILABLE
        assert effective.value == pytest.approx(sample_count / tau.value)
        assert mcse.status is McmcDiagnosticStatus.AVAILABLE
        assert mcse.value == pytest.approx(
            item.posterior_standard_deviation / np.sqrt(effective.value)
        )


def test_recomputed_outer_hash_cannot_repair_altered_child_diagnostics() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)
    posterior = derive_posterior_sample_evidence(selection, plan.coordinate_units)
    summary = derive_posterior_summary(posterior)

    object.__setattr__(summary.acceptance, "mean_acceptance_fraction", -4.0)
    object.__setattr__(
        summary.acceptance, "identity", summary.acceptance._content_identity()
    )
    object.__setattr__(summary, "identity", summary._content_identity())

    with pytest.raises(McmcConstructionError, match="Acceptance diagnostic"):
        summary.validate_integrity()


def test_complete_scope_posterior_outcomes_and_summary_are_canonical() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)

    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    summary = derive_posterior_summary(posterior)

    assert (
        tuple((item.state_ordinal, item.walker_ordinal) for item in posterior.outcomes)
        == selection.sample_indices
    )
    assert all(
        item.disposition is PosteriorSampleDisposition.SUCCESS
        and item.resolved_items is not None
        and tuple(key for key, _value in item.resolved_items) == ("A", "B")
        for item in posterior.outcomes
    )
    assert summary.included_labels == selection.sample_indices
    assert not summary.excluded_labels
    assert len(summary.parameter_summaries) == 2
    assert all(
        item.credible_interval_name == "equal_tailed_95_percent"
        and item.posterior_standard_deviation >= 0.0
        and item.autocorrelation_time.status is McmcDiagnosticStatus.UNAVAILABLE
        for item in summary.parameter_summaries
    )
    assert len(summary.covariance) == 4
    assert len(summary.correlations) == 4

    with pytest.raises(McmcConstructionError, match="witness"):
        dataclasses.replace(summary, included_labels=())


def test_posterior_samples_combine_exact_held_values_before_scope_projection() -> None:
    accepted, problem, parameterization, engine = _native_context(fit_b=False)
    policy = ExpertMcmcPolicy(
        burn_steps=1,
        retained_steps=2,
        walkers=4,
        expert_provenance="held-value-qualification",
    ).resolve(dimension=1, root_seed=9876)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=policy,
        coordinate_units=(("A", ParameterUnit.DIMENSIONLESS),),
    )
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    selection = derive_retained_sample_view(operation.evidence)

    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )

    assert all(
        outcome.resolved_items is not None and dict(outcome.resolved_items)["B"] == 1.5
        for outcome in posterior.outcomes
    )


def test_primary_tampering_cannot_reach_any_derived_posterior_artifact() -> None:
    accepted, plan = _plan_context()
    operation = execute_mcmc_evidence(accepted, plan)
    assert operation.evidence is not None
    evidence = operation.evidence
    with pytest.raises(McmcConstructionError, match="contiguous"):
        dataclasses.replace(
            evidence,
            states=(evidence.states[0], *reversed(evidence.states[1:])),
        )


def test_all_operation_and_posterior_paths_preserve_source_values_bitwise() -> None:
    accepted, problem, parameterization, engine = _native_context()
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=problem,
        parameterization=parameterization,
        source_engine=engine,
        policy=_resolved_policy(),
        coordinate_units=(
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    before = problem.source_snapshot.to_json()
    completed = execute_mcmc_evidence(accepted, plan)
    cancellation = CancellationToken()
    cancellation.cancel()
    cancelled = execute_mcmc_evidence(accepted, plan, cancellation=cancellation)

    def interrupt(stage: McmcExecutionStage, _count: int) -> None:
        if stage is McmcExecutionStage.INITIALIZING:
            raise KeyboardInterrupt

    interrupted = execute_mcmc_evidence(
        accepted,
        plan,
        checkpoint_observer=interrupt,
    )
    assert completed.evidence is not None
    selection = derive_retained_sample_view(completed.evidence)
    posterior = derive_posterior_sample_evidence(
        selection,
        (
            ("A", ParameterUnit.DIMENSIONLESS),
            ("B", ParameterUnit.DIMENSIONLESS),
        ),
    )
    derive_posterior_summary(posterior)

    assert cancelled.terminal is McmcOperationTerminal.CANCELLED
    assert interrupted.terminal is McmcOperationTerminal.INTERRUPTED
    assert problem.source_snapshot.to_json() == before
    assert accepted.vector == (1.0, 2.0)
