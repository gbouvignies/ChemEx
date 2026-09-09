"""Workflow qualification for Eyring uncertainty, MCMC, and resampling."""

from __future__ import annotations

import math
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pydantic import BaseModel

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultSetting
from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.models.factory import model_factory
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.optimize.deterministic_uncertainty import (
    AcceptedDeterministicFitFacts,
    ContinuousTrfBasis,
    DeterministicUncertainty,
    derive_deterministic_uncertainty,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    DirectTrfInvocation,
    OptimizationProblem,
    canonical_chi_square,
    execute_direct_trf,
)
from chemex.optimize.grouped_direct_trf import FitDecomposition
from chemex.optimize.native_mcmc import (
    McmcOperationTerminal,
    McmcPlan,
    execute_mcmc_evidence,
    resolve_product_mcmc_policy,
)
from chemex.optimize.native_resampling import (
    OperationTerminal,
    OptimizationStrategy,
    ResamplingDatasetManifest,
    ResamplingPlan,
    ResamplingScheme,
    execute_resampling_evidence,
)
from chemex.optimize.uncertainty import ParameterUnit, UncertaintyUnavailableKind
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ActiveParameterization,
    SealedParameterModel,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.data import Printer
from chemex.printers.parameters import write_parameters
from chemex.runtime import AnalysisSession


class _EyringKernelSettings(BaseModel):
    kind: str = "eyring-workflow-qualification"


class _EyringSpectrometer:
    def __init__(self, spin_system: SpinSystem) -> None:
        self.spin_system = spin_system
        self.values = {"kab": 0.0, "pb": 0.0}

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _EyringSpectrometer:
        return deepcopy(self)

    def native_kernel_descriptor(self) -> dict[str, str]:
        return {"kind": "eyring-workflow-spectrometer"}


class _EyringPulseSequence:
    settings = _EyringKernelSettings()

    def calculate(self, spectrometer: _EyringSpectrometer, data: Data) -> np.ndarray:
        metadata = np.asarray(data.metadata, dtype=np.float64)
        return spectrometer.values["kab"] + spectrometer.values["pb"] * metadata

    def is_reference(self, metadata: np.ndarray) -> np.ndarray:
        return np.zeros(metadata.shape, dtype=np.bool_)


def _decimal_eyring_rate(
    temperature: float,
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
) -> float:
    """Evaluate a directional rate independently from exact SI literals."""
    with localcontext() as context:
        context.prec = 100
        gas_constant = Decimal("8.31446261815324")
        frequency_factor = Decimal("1.380649e-23") / Decimal("6.62607015e-34")
        kelvin = Decimal(str(temperature)) + Decimal("273.15")
        activation_enthalpy = Decimal(str(transition_enthalpy)) - Decimal(
            str(initial_enthalpy)
        )
        activation_entropy = Decimal(str(transition_entropy)) - Decimal(
            str(initial_entropy)
        )
        exponent = activation_entropy / gas_constant - activation_enthalpy / (
            gas_constant * kelvin
        )
        return float(frequency_factor * kelvin * exponent.exp())


def _decimal_two_state_population(
    temperature: float,
    state_enthalpy: float,
    state_entropy: float,
) -> float:
    """Evaluate P(B) independently from exact SI Boltzmann weights."""
    with localcontext() as context:
        context.prec = 100
        gas_constant = Decimal("8.31446261815324")
        kelvin = Decimal(str(temperature)) + Decimal("273.15")
        log_weight = Decimal(str(state_entropy)) / gas_constant - Decimal(
            str(state_enthalpy)
        ) / (gas_constant * kelvin)
        weight = log_weight.exp()
        return float(weight / (Decimal(1) + weight))


@dataclass(frozen=True, slots=True)
class _EyringWorkflow:
    session: AnalysisSession
    local_ids: dict[str, str]
    parameter_model: SealedParameterModel
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    problem: OptimizationProblem
    accepted: AcceptedFitResult
    profile: Profile


@dataclass(frozen=True, slots=True)
class _MultiTemperatureEyringWorkflow:
    session: AnalysisSession
    temperatures: tuple[float, ...]
    local_ids: dict[float, dict[str, str]]
    parameter_model: SealedParameterModel
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    problem: OptimizationProblem
    accepted: AcceptedFitResult
    profiles: tuple[Profile, ...]


def _build_eyring_workflow(
    *,
    temperature: float = 25.0,
    true_dh_b: float = 8_000.0,
    true_ds_b: float = 10.0,
    true_dh_ab: float = 65_000.0,
    true_ds_ab: float = 20.0,
    fitted_names: tuple[str, ...] = ("DH_AB", "DH_B"),
    execute_fit: bool = True,
    error: float = 0.02,
) -> _EyringWorkflow:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=temperature,
        p_total=1.0e-3,
        l_total=2.0e-3,
    )
    session = AnalysisSession.create()
    session.set_model("2st_eyring")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    settings = model_factory.create("2st_eyring", conditions)
    local_ids = {
        name: setting.name_setting.get_param_name(spin_system, conditions).id_
        for name, setting in settings.items()
    }
    config = SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=spin_system,
    )
    initial_values = {
        "dh_b": true_dh_b if not execute_fit else true_dh_b - 500.0,
        "ds_b": true_ds_b,
        "dh_ab": true_dh_ab if not execute_fit else true_dh_ab - 500.0,
        "ds_ab": true_ds_ab,
    }
    session.parameters.set_defaults(
        [
            (ParamName.from_section(name), DefaultSetting(value))
            for name, value in initial_values.items()
        ]
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )

    metadata = np.asarray((0.0, 40.0, 90.0, 150.0, 220.0, 300.0))
    true_rate = _decimal_eyring_rate(
        temperature,
        0.0,
        0.0,
        true_dh_ab,
        true_ds_ab,
    )
    true_population = _decimal_two_state_population(
        temperature,
        true_dh_b,
        true_ds_b,
    )
    noise = np.asarray((0.010, -0.015, 0.006, 0.018, -0.011, 0.004))
    observed = true_rate + true_population * metadata
    if execute_fit:
        observed = observed + noise
    data = Data(
        exp=np.asarray(observed, dtype=np.float64),
        err=np.full(metadata.shape, error, dtype=np.float64),
        metadata=metadata,
    )
    profile = Profile(
        data,
        cast("Spectrometer", _EyringSpectrometer(spin_system)),
        cast("PulseSequence", _EyringPulseSequence()),
        {"kab": local_ids["kab"], "pb": local_ids["pb"]},
        cast("Printer", None),
        is_scaled=False,
    )
    experiments = cast("Any", (SimpleNamespace(profiles=(profile,)),))
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    requested_ids = {local_ids[name] for name in ("kab", "kba", "pa", "pb")}
    required_ids = profile.param_ids | requested_ids
    independent_names = {"DH_B", "DS_B", "DH_AB", "DS_AB"}
    parameterization = session.compile_parameterization(
        Method(
            fit=list(fitted_names),
            fix=sorted(independent_names.difference(fitted_names)),
        ),
        required_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    if execute_fit:
        outcome = execute_direct_trf(
            problem,
            DirectTrfInvocation.for_problem(problem, objective_request_budget=100),
            parameterization,
            engine,
        )
        accepted = outcome.accepted_result
        assert accepted is not None
    else:
        lifecycle = problem.lifecycle_frame(problem.start, parameterization)
        frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
        evaluation = engine.new_evaluator().evaluate(frame)
        assert isinstance(evaluation, EvaluationResult)
        accepted = AcceptedFitResult.for_qualification(
            occurrence_identity="eyring-workflow-accepted-occurrence",
            problem_identity=problem.identity,
            invocation_identity="eyring-workflow-invocation",
            execution_identity="eyring-workflow-execution",
            materialization_identity="eyring-workflow-materialization",
            parameterization_identity=parameterization.identity,
            evaluator_parameterization_identity=parameterization.evaluator_identity,
            source_occurrence_identity=problem.source_snapshot.occurrence_identity,
            source_revision=problem.source_snapshot.revision,
            controlled_ids=problem.controlled_ids,
            vector=problem.start,
            chi_square=canonical_chi_square(evaluation.residuals),
            evaluation_result=evaluation,
            commit_scope=problem.commit_scope,
            commit_items=evaluation.resolved_values.ordered_items(),
            origin_context_identity="eyring-workflow-qualification",
        )
    return _EyringWorkflow(
        session,
        local_ids,
        parameter_model,
        parameterization,
        engine,
        problem,
        accepted,
        profile,
    )


def _build_multitemperature_eyring_workflow() -> _MultiTemperatureEyringWorkflow:
    temperatures = (-20.0, 25.0, 70.0)
    true_coordinates = {
        "dh_b": 8_000.0,
        "ds_b": 10.0,
        "dh_ab": 65_000.0,
        "ds_ab": 20.0,
    }
    session = AnalysisSession.create()
    session.set_model("2st_eyring")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    local_ids: dict[float, dict[str, str]] = {}

    for temperature in temperatures:
        conditions = Conditions(
            h_larmor_frq=600.0,
            temperature=temperature,
            p_total=1.0e-3,
            l_total=2.0e-3,
        )
        settings = model_factory.create("2st_eyring", conditions)
        local_ids[temperature] = {
            name: setting.name_setting.get_param_name(spin_system, conditions).id_
            for name, setting in settings.items()
        }
        config = SimpleNamespace(
            conditions=conditions,
            to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
        )
        session.parameter_factory.create_parameters(
            config,  # ty: ignore[invalid-argument-type]
            basis=basis,
            spin_system=spin_system,
        )

    initial_values = {
        "dh_b": true_coordinates["dh_b"] - 300.0,
        "ds_b": true_coordinates["ds_b"] - 1.0,
        "dh_ab": true_coordinates["dh_ab"] - 300.0,
        "ds_ab": true_coordinates["ds_ab"] - 1.0,
    }
    session.parameters.set_defaults(
        [
            (ParamName.from_section(name), DefaultSetting(value))
            for name, value in initial_values.items()
        ]
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )

    metadata = np.asarray((0.0, 40.0, 90.0, 150.0, 220.0, 300.0))
    noise = np.asarray((0.010, -0.015, 0.006, 0.018, -0.011, 0.004))
    profiles: list[Profile] = []
    for temperature in temperatures:
        rate = _decimal_eyring_rate(
            temperature,
            0.0,
            0.0,
            true_coordinates["dh_ab"],
            true_coordinates["ds_ab"],
        )
        population = _decimal_two_state_population(
            temperature,
            true_coordinates["dh_b"],
            true_coordinates["ds_b"],
        )
        data = Data(
            exp=np.asarray(rate + population * metadata + noise, dtype=np.float64),
            err=np.full(metadata.shape, 0.02, dtype=np.float64),
            metadata=metadata,
        )
        profile = Profile(
            data,
            cast("Spectrometer", _EyringSpectrometer(spin_system)),
            cast("PulseSequence", _EyringPulseSequence()),
            {
                "kab": local_ids[temperature]["kab"],
                "pb": local_ids[temperature]["pb"],
            },
            cast("Printer", None),
            is_scaled=False,
        )
        profiles.append(profile)

    profile_tuple = tuple(profiles)
    experiments = cast("Any", (SimpleNamespace(profiles=profile_tuple),))
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    requested_ids = {
        local_ids[temperature][name]
        for temperature in temperatures
        for name in ("kab", "kba", "pa", "pb")
    }
    required_ids = requested_ids | {
        param_id for profile in profile_tuple for param_id in profile.param_ids
    }
    parameterization = session.compile_parameterization(
        Method(fit=["DH_B", "DS_B", "DH_AB", "DS_AB"]),
        required_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    outcome = execute_direct_trf(
        problem,
        DirectTrfInvocation.for_problem(problem, objective_request_budget=300),
        parameterization,
        engine,
    )
    accepted = outcome.accepted_result
    assert accepted is not None
    return _MultiTemperatureEyringWorkflow(
        session,
        temperatures,
        local_ids,
        parameter_model,
        parameterization,
        engine,
        problem,
        accepted,
        profile_tuple,
    )


def _derive_uncertainty(
    workflow: _EyringWorkflow | _MultiTemperatureEyringWorkflow,
) -> DeterministicUncertainty:
    decomposition = FitDecomposition.from_root(
        workflow.problem,
        workflow.parameterization,
        workflow.engine,
    )
    return derive_deterministic_uncertainty(
        AcceptedDeterministicFitFacts(
            workflow.accepted,
            workflow.problem,
            workflow.parameterization,
            workflow.engine,
            ContinuousTrfBasis(decomposition.partition_proof),
            "eyring-workflow-uncertainty",
        )
    )


def test_real_fit_reports_directional_rate_and_population_uncertainty(
    tmp_path: Path,
) -> None:
    workflow = _build_eyring_workflow()
    uncertainty = _derive_uncertainty(workflow)
    output_names = ("kab", "kba", "pa", "pb")

    assert all(
        (conclusion := uncertainty.parameter(workflow.local_ids[name])) is not None
        and conclusion.reportable
        and conclusion.standard_error is not None
        and math.isfinite(conclusion.standard_error)
        for name in output_names
    )
    resolved = dict(workflow.accepted.commit_items)
    write_parameters(
        tmp_path,
        parameter_model=workflow.parameter_model,
        parameter_values=resolved,
        parameterization=workflow.parameterization,
        fitted_ids=workflow.problem.controlled_ids,
        deterministic_uncertainty=uncertainty,
    )
    constrained = (tmp_path / "Parameters" / "constrained.toml").read_text()
    for name in ("KAB", "KBA", "PA", "PB"):
        assert any(name in line and "# ±" in line for line in constrained.splitlines())


def test_one_temperature_is_rank_deficient_but_multiple_temperatures_add_rank() -> None:
    gas_constant = 8.31446261815324

    def row(kelvin: float) -> tuple[float, float]:
        return (-1.0 / (gas_constant * kelvin), 1.0 / gas_constant)

    one_temperature = np.asarray((row(298.15), row(298.15)))
    multiple_temperatures = np.asarray(
        tuple(row(value) for value in (298.15, 303.15, 308.15))
    )

    assert np.linalg.matrix_rank(one_temperature) == 1
    assert np.linalg.matrix_rank(multiple_temperatures) == 2
    covariance = np.linalg.inv(multiple_temperatures.T @ multiple_temperatures)
    correlation = covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])
    assert correlation > 0.999


def test_multitemperature_fit_qualifies_shared_h_s_and_derived_output(
    tmp_path: Path,
) -> None:
    workflow = _build_multitemperature_eyring_workflow()
    uncertainty = _derive_uncertainty(workflow)
    first_ids = workflow.local_ids[workflow.temperatures[0]]

    for name in ("dh_b", "ds_b", "dh_ab", "ds_ab"):
        assert {
            workflow.local_ids[temperature][name]
            for temperature in workflow.temperatures
        } == {first_ids[name]}
    for name in ("kab", "kba", "pa", "pb"):
        assert len(
            {
                workflow.local_ids[temperature][name]
                for temperature in workflow.temperatures
            }
        ) == len(workflow.temperatures)

    root = uncertainty.root_evidence
    assert root is not None
    assert root.rank_diagnostic is not None
    assert root.rank_diagnostic.rank == len(workflow.problem.controlled_ids) == 4
    assert root.covariance is not None
    assert math.isfinite(root.covariance.jacobian_condition)
    assert root.covariance.usable
    assert root.correlations is not None
    dh_index = root.correlations.output_ids.index(first_ids["dh_ab"])
    ds_index = root.correlations.output_ids.index(first_ids["ds_ab"])
    correlation = root.correlations.entries[dh_index][ds_index].value
    assert correlation is not None
    assert 0.98 < abs(correlation) < 1.0

    for name in ("dh_b", "ds_b", "dh_ab", "ds_ab"):
        conclusion = uncertainty.parameter(first_ids[name])
        assert conclusion is not None
        assert conclusion.reportable
        assert conclusion.standard_error is not None
        assert math.isfinite(conclusion.standard_error)
    for temperature in workflow.temperatures:
        for name in ("kab", "kba", "pa", "pb"):
            conclusion = uncertainty.parameter(workflow.local_ids[temperature][name])
            assert conclusion is not None
            assert conclusion.reportable
            assert conclusion.standard_error is not None
            assert math.isfinite(conclusion.standard_error)

    resolved = dict(workflow.accepted.commit_items)
    for temperature in workflow.temperatures:
        expected_rate = _decimal_eyring_rate(
            temperature,
            0.0,
            0.0,
            resolved[first_ids["dh_ab"]],
            resolved[first_ids["ds_ab"]],
        )
        actual_rate = resolved[workflow.local_ids[temperature]["kab"]]
        assert actual_rate == pytest.approx(expected_rate, rel=3.0e-14)

    write_parameters(
        tmp_path,
        parameter_model=workflow.parameter_model,
        parameter_values=resolved,
        parameterization=workflow.parameterization,
        fitted_ids=workflow.problem.controlled_ids,
        deterministic_uncertainty=uncertainty,
    )
    fitted = (tmp_path / "Parameters" / "fitted.toml").read_text()
    constrained = (tmp_path / "Parameters" / "constrained.toml").read_text()
    assert all(
        any(name in line and "# ±" in line for line in fitted.splitlines())
        for name in ("DH_B", "DS_B", "DH_AB", "DS_AB")
    )
    for temperature in workflow.temperatures:
        qualifier = f"{temperature:.1f}C"
        assert qualifier in constrained
    assert all(
        any(name in line and "# ±" in line for line in constrained.splitlines())
        for name in ("KAB", "KBA", "PA", "PB")
    )


def test_one_temperature_h_s_fit_reports_rank_deficient_uncertainty() -> None:
    workflow = _build_eyring_workflow(fitted_names=("DH_AB", "DS_AB"))
    uncertainty = _derive_uncertainty(workflow)

    assert all(
        (conclusion := uncertainty.parameter(param_id)) is not None
        and not conclusion.reportable
        and conclusion.unavailable_kind is UncertaintyUnavailableKind.RANK_DEFICIENT
        for param_id in workflow.problem.controlled_ids
    )


def test_mcmc_retries_an_unavailable_eyring_activation_proposal() -> None:
    workflow = _build_eyring_workflow(
        temperature=-273.14,
        true_dh_b=0.0,
        true_ds_b=0.0,
        true_dh_ab=0.0,
        true_ds_ab=0.0,
        fitted_names=("DH_AB",),
        execute_fit=False,
        error=1.0e300,
    )
    plan = McmcPlan.for_accepted(
        workflow.accepted,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        policy=resolve_product_mcmc_policy(
            dimension=1,
            walkers=8,
            steps=2,
            root_seed=770,
        ),
        coordinate_units=(
            (workflow.problem.controlled_ids[0], ParameterUnit.ENERGY_PER_MOLE),
        ),
    )
    evaluator = workflow.engine.new_evaluator()
    initial_evaluations = tuple(
        evaluator.evaluate(
            EvaluationFrame.from_lifecycle_frame(
                workflow.parameterization,
                workflow.problem.lifecycle_frame(vector, workflow.parameterization),
            )
        )
        for vector in plan.initial_ensemble
    )
    assert any(
        isinstance(evaluation, EvaluationFailure)
        and evaluation.stage == "resolution"
        and evaluation.category == "domain_error"
        and "function_id='eyring_rate'" in evaluation.message
        for evaluation in initial_evaluations
    )

    operation = execute_mcmc_evidence(workflow.accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    initial_state = operation.evidence.states[0]
    assert initial_state.positions != plan.initial_ensemble
    assert all(math.isfinite(value) for value in initial_state.log_densities)


def test_mcmc_retries_an_unrepresentable_eyring_population_proposal() -> None:
    workflow = _build_eyring_workflow(
        temperature=-273.14,
        true_dh_b=61.72329710480702,
        true_ds_b=0.0,
        true_dh_ab=61.72329710480702,
        true_ds_ab=0.0,
        fitted_names=("DH_B",),
        execute_fit=False,
        error=1.0e300,
    )
    plan = McmcPlan.for_accepted(
        workflow.accepted,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        policy=resolve_product_mcmc_policy(
            dimension=1,
            walkers=8,
            steps=2,
            root_seed=771,
        ),
        coordinate_units=(
            (workflow.problem.controlled_ids[0], ParameterUnit.ENERGY_PER_MOLE),
        ),
    )
    evaluator = workflow.engine.new_evaluator()
    initial_evaluations = tuple(
        evaluator.evaluate(
            EvaluationFrame.from_lifecycle_frame(
                workflow.parameterization,
                workflow.problem.lifecycle_frame(vector, workflow.parameterization),
            )
        )
        for vector in plan.initial_ensemble
    )
    assert any(
        isinstance(evaluation, EvaluationFailure)
        and evaluation.stage == "resolution"
        and evaluation.category == "domain_error"
        and "function_id='pop_2st_eyring'" in evaluation.message
        for evaluation in initial_evaluations
    )

    operation = execute_mcmc_evidence(workflow.accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    assert all(
        math.isfinite(value) for value in operation.evidence.states[0].log_densities
    )


def test_monte_carlo_recomputes_eyring_rates_from_resampled_coordinates() -> None:
    workflow = _build_eyring_workflow()
    accepted = workflow.accepted
    size = workflow.engine.plan.observation_count
    dataset = ResamplingDatasetManifest(
        workflow.engine.plan,
        tuple(
            float(value) for value in accepted.evaluation_result.normalized_calculations
        ),
        tuple(False for _ in range(size)),
        tuple("G23" for _ in range(size)),
        tuple(f"eyring-observation-{index}" for index in range(size)),
    )
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=dataset,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        scheme=ResamplingScheme.MONTE_CARLO,
        replicate_count=3,
        replicate_structural_identities=tuple(
            f"eyring-monte-carlo-{index}" for index in range(3)
        ),
        replicate_component_identities=tuple(
            (f"eyring-direct-trf-{index}",) for index in range(3)
        ),
        root_seed=770,
        output_scope=workflow.problem.commit_scope,
        output_units=("native",) * len(workflow.problem.commit_scope),
        minimum_successful_count=3,
        strategy=OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "100"),),
    )

    operation = execute_resampling_evidence(accepted, plan)

    assert operation.terminal is OperationTerminal.COMPLETED
    assert operation.evidence is not None
    assert operation.evidence.successful_count == 3
    for outcome in operation.evidence.outcomes:
        assert outcome.success is not None
        resolved = dict(outcome.success.resolved_items)
        expected_rate = _decimal_eyring_rate(
            25.0,
            0.0,
            0.0,
            resolved[workflow.local_ids["dh_ab"]],
            resolved[workflow.local_ids["ds_ab"]],
        )
        assert resolved[workflow.local_ids["kab"]] == pytest.approx(
            expected_rate,
            rel=3.0e-14,
        )
