"""Focused Direct-TRF, MCMC, resampling, and uncertainty N-state workflows."""

from __future__ import annotations

import math
import tomllib
from argparse import Namespace
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pytest
from pydantic import BaseModel

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.configuration.parameters import (
    DefaultListType,
    DefaultSetting,
    read_defaults,
)
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
from chemex.optimize import de_direct_trf as de_module
from chemex.optimize import native_mcmc as native_mcmc_module
from chemex.optimize import uncertainty as uncertainty_module
from chemex.optimize.de_direct_trf import (
    DeSearchInvocation,
    DeSearchTerminal,
    execute_de_search,
)
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
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import ActiveParameterization, ParameterRole
from chemex.parameters.spin_system import SpinSystem
from chemex.run_info import write_run_info
from chemex.runtime import AnalysisSession


class _NStateKernelSettings(BaseModel):
    kind: str = "nstate-workflow-qualification"


class _NStateSpectrometer:
    def __init__(self, spin_system: SpinSystem) -> None:
        self.spin_system = spin_system
        self.values = {
            "pb": 0.0,
            "pc": 0.0,
            "pd": 0.0,
            "pe": 0.0,
            "pf": 0.0,
        }

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _NStateSpectrometer:
        return deepcopy(self)

    def native_kernel_descriptor(self) -> dict[str, str]:
        return {"kind": "nstate-workflow-spectrometer"}


class _NStatePulseSequence:
    settings = _NStateKernelSettings()

    def calculate(
        self,
        spectrometer: _NStateSpectrometer,
        data: Data,
    ) -> np.ndarray:
        metadata = np.asarray(data.metadata, dtype=np.float64)
        return sum(
            spectrometer.values[name] * metadata**power
            for power, name in enumerate(("pb", "pc", "pd", "pe", "pf"), 1)
            if name in spectrometer.values
        )

    def is_reference(self, metadata: np.ndarray) -> np.ndarray:
        return np.zeros(metadata.shape, dtype=np.bool_)


@dataclass(frozen=True, slots=True)
class _NStateWorkflow:
    session: AnalysisSession
    local_ids: dict[str, str]
    component_local_ids: tuple[dict[str, str], ...]
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    problem: OptimizationProblem


def _fit_workflow(workflow: _NStateWorkflow) -> AcceptedFitResult:
    outcome = execute_direct_trf(
        workflow.problem,
        DirectTrfInvocation.for_problem(
            workflow.problem,
            objective_request_budget=100,
        ),
        workflow.parameterization,
        workflow.engine,
    )
    accepted = outcome.accepted_result
    assert accepted is not None
    return accepted


def _qualified_accepted_at(
    workflow: _NStateWorkflow,
    template: AcceptedFitResult,
    vector: tuple[float, ...],
) -> AcceptedFitResult:
    """Materialize an exact valid vector for isolated production qualification."""
    lifecycle = workflow.problem.lifecycle_frame(vector, workflow.parameterization)
    evaluated = workflow.engine.new_evaluator().evaluate(
        EvaluationFrame.from_lifecycle_frame(workflow.parameterization, lifecycle)
    )
    assert isinstance(evaluated, EvaluationResult)
    return AcceptedFitResult.for_qualification(
        occurrence_identity=f"nstate-exact-{template.occurrence_identity}",
        problem_identity=workflow.problem.identity,
        invocation_identity=template.invocation_identity,
        execution_identity=template.execution_identity,
        materialization_identity=f"exact-{template.materialization_identity}",
        parameterization_identity=template.parameterization_identity,
        evaluator_parameterization_identity=(
            template.evaluator_parameterization_identity
        ),
        source_occurrence_identity=template.source_occurrence_identity,
        source_revision=template.source_revision,
        controlled_ids=template.controlled_ids,
        vector=vector,
        chi_square=canonical_chi_square(evaluated.residuals),
        evaluation_result=evaluated,
        commit_scope=template.commit_scope,
        commit_items=tuple(
            (param_id, evaluated.resolved_values[param_id])
            for param_id in template.commit_scope
        ),
        origin_context_identity=template.origin_context_identity,
    )


def _build_workflow(
    *,
    start: dict[str, float],
    target: dict[str, float],
    fitted_name: str | tuple[str, ...],
    include_rates: bool = False,
    temperatures: tuple[float, ...] = (25.0,),
    model_name: str = "4st_fork",
    defaults: DefaultListType | None = None,
) -> _NStateWorkflow:
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    population_names = tuple(f"p{state}" for state in session.model.states[1:])
    component_local_ids: list[dict[str, str]] = []
    for temperature in temperatures:
        conditions = Conditions(
            h_larmor_frq=600.0,
            temperature=temperature,
            p_total=1.0e-3,
            l_total=2.0e-3,
        )
        settings = model_factory.create(model_name, conditions)
        component_local_ids.append(
            {
                name: setting.name_setting.get_param_name(spin_system, conditions).id_
                for name, setting in settings.items()
            }
        )
        config = SimpleNamespace(
            conditions=conditions,
            to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
        )
        session.parameter_factory.create_parameters(
            config,  # ty: ignore[invalid-argument-type]
            basis=basis,
            spin_system=spin_system,
        )
    session.parameters.set_defaults(
        [
            (ParamName.from_section(name), DefaultSetting(value))
            for name, value in start.items()
        ]
        if defaults is None
        else defaults
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )

    metadata = np.asarray((0.25, 0.5, 0.8, 1.1, 1.5, 2.0))
    observed = sum(
        target[name] * metadata**power for power, name in enumerate(population_names, 1)
    )
    data = Data(
        exp=np.asarray(observed, dtype=np.float64),
        err=np.full(metadata.shape, 1.0e-4),
        metadata=metadata,
    )
    profiles = tuple(
        Profile(
            deepcopy(data),
            cast("Spectrometer", _NStateSpectrometer(spin_system)),
            cast("PulseSequence", _NStatePulseSequence()),
            {name: local_ids[name] for name in population_names},
            cast("Any", None),
            is_scaled=False,
        )
        for local_ids in component_local_ids
    )
    experiments = cast("Any", (SimpleNamespace(profiles=profiles),))
    required_ids = set().union(*(profile.param_ids for profile in profiles))
    required_ids.update(local_ids["pa"] for local_ids in component_local_ids)
    if include_rates:
        for local_ids in component_local_ids:
            required_ids.update((local_ids["kab"], local_ids["kba"]))
    fitted_names = (fitted_name,) if isinstance(fitted_name, str) else fitted_name
    held_names = sorted({name.upper() for name in population_names} - set(fitted_names))
    if include_rates:
        held_names.append("KEX_AB")
    parameterization = session.compile_parameterization(
        Method(fit=list(fitted_names), fix=held_names),
        required_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    return _NStateWorkflow(
        session,
        component_local_ids[0],
        tuple(component_local_ids),
        parameterization,
        engine,
        problem,
    )


def _derive_uncertainty(
    workflow: _NStateWorkflow,
    accepted: AcceptedFitResult,
) -> DeterministicUncertainty:
    decomposition = FitDecomposition.from_root(
        workflow.problem,
        workflow.parameterization,
        workflow.engine,
    )
    return derive_deterministic_uncertainty(
        AcceptedDeterministicFitFacts(
            accepted,
            workflow.problem,
            workflow.parameterization,
            workflow.engine,
            ContinuousTrfBasis(decomposition.partition_proof),
            "nstate-workflow-uncertainty",
        )
    )


@pytest.mark.parametrize(
    ("fitted_name", "start", "target", "boundary_name"),
    (
        (
            "PB",
            {"pb": 0.1, "pc": 0.25, "pd": 0.35},
            {"pb": 0.0, "pc": 0.25, "pd": 0.35},
            "pb",
        ),
        (
            "PD",
            {"pb": 0.2, "pc": 0.3, "pd": 0.4},
            {"pb": 0.2, "pc": 0.3, "pd": 0.5},
            "pa",
        ),
        (
            "PD",
            {"pb": 0.2, "pc": 0.3, "pd": 0.1},
            {"pb": 0.2, "pc": 0.3, "pd": 0.0},
            "pd",
        ),
    ),
)
def test_direct_trf_reaches_closed_simplex_boundaries(
    fitted_name: str,
    start: dict[str, float],
    target: dict[str, float],
    boundary_name: str,
) -> None:
    workflow = _build_workflow(start=start, target=target, fitted_name=fitted_name)

    accepted = _fit_workflow(workflow)
    resolved = accepted.evaluation_result.resolved_values
    assert resolved[workflow.local_ids[boundary_name]] == pytest.approx(
        0.0,
        abs=1.0e-7,
    )
    assert all(
        resolved[workflow.local_ids[name]] >= 0.0 for name in ("pa", "pb", "pc", "pd")
    )
    assert sum(
        resolved[workflow.local_ids[name]] for name in ("pa", "pb", "pc", "pd")
    ) == pytest.approx(1.0)


def test_selected_coordinate_de_uses_public_nstate_simplex_values() -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.2},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.45},
        fitted_name="PD",
    )
    selected_id = workflow.problem.controlled_ids[0]
    invocation = DeSearchInvocation.for_product_problem(
        workflow.problem,
        search_coordinates=((selected_id, 0.0, 1.0, "linear"),),
        root_seed=773,
    )
    chart = invocation.search_problem.feasible_coordinates
    assert chart is not None
    assert chart.population_simplexes
    assert not chart.uses_private_relaxation_coordinates
    assert invocation.search_problem.controlled_ids == (selected_id,)
    assert tuple(item.param_id for item in invocation.search_coordinates) == (
        selected_id,
    )

    def valid_and_exterior_backend(live, _invocation, _solver_start):
        valid = np.asarray((0.45,), dtype=np.float64)
        exterior = np.asarray((0.8,), dtype=np.float64)
        valid_objective = live.objective(valid)
        assert math.isfinite(valid_objective)
        assert live.objective(exterior) == math.inf
        return SimpleNamespace(
            success=True,
            message="Optimization terminated successfully.",
            nit=1,
            nfev=2,
            x=valid,
            fun=valid_objective,
            population=np.tile(valid, (invocation.population.size, 1)),
            population_energies=np.full(
                invocation.population.size,
                valid_objective,
            ),
        )

    with patch.object(
        de_module,
        "_invoke_de_backend",
        side_effect=valid_and_exterior_backend,
    ):
        outcome = execute_de_search(
            workflow.problem,
            invocation,
            workflow.parameterization,
            workflow.engine,
        )

    assert outcome.terminal is DeSearchTerminal.POPULATION_CONVERGED
    assert outcome.valid_candidate_count == 1
    assert outcome.rejected_trial_count == 1
    assert outcome.best_candidate is not None
    assert outcome.best_candidate.selected_vector == (0.45,)


def test_mcmc_accepts_exact_nstate_simplex_boundary_and_retries_exterior() -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.4},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.5},
        fitted_name="PD",
    )
    template = _fit_workflow(workflow)
    accepted = _qualified_accepted_at(workflow, template, (0.5,))
    assert accepted.evaluation_result.resolved_values[workflow.local_ids["pa"]] == 0.0
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        policy=resolve_product_mcmc_policy(
            dimension=1,
            walkers=8,
            steps=1,
            root_seed=774,
        ),
        coordinate_units=(
            (workflow.problem.controlled_ids[0], ParameterUnit.FRACTION),
        ),
    )
    kernel = native_mcmc_module._LogDensityKernel(
        native_mcmc_module._McmcWorkerContext.from_plan(plan)
    )

    assert math.isfinite(kernel.evaluate(np.asarray((0.5,))).value)
    assert kernel.evaluate(np.asarray((0.5001,))).value == -math.inf

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    initial = operation.evidence.states[0]
    assert all(math.isfinite(value) for value in initial.log_densities)
    attempts = native_mcmc_module._initialization_attempts_for_positions(
        plan,
        initial.positions,
    )
    assert attempts is not None
    assert any(attempt > 0 for attempt in attempts)


@pytest.mark.parametrize(
    ("model_name", "populations", "root_seed"),
    (
        (
            "5st_fork",
            {"pb": 0.1, "pc": 0.2, "pd": 0.3, "pe": 0.4},
            775,
        ),
        (
            "6st_fork",
            {"pb": 0.1, "pc": 0.15, "pd": 0.2, "pe": 0.25, "pf": 0.3},
            776,
        ),
    ),
)
def test_higher_state_mcmc_initializer_retries_to_valid_simplex_walkers(
    model_name: str,
    populations: dict[str, float],
    root_seed: int,
) -> None:
    fitted_names = tuple(name.upper() for name in populations)
    workflow = _build_workflow(
        start=populations,
        target=populations,
        fitted_name=fitted_names,
        model_name=model_name,
    )
    template = _fit_workflow(workflow)
    accepted = _qualified_accepted_at(
        workflow,
        template,
        tuple(populations[name] for name in populations),
    )
    dimension = len(accepted.vector)
    plan = McmcPlan.for_accepted(
        accepted,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        policy=resolve_product_mcmc_policy(
            dimension=dimension,
            walkers=2 * dimension + 2,
            steps=1,
            root_seed=root_seed,
        ),
        coordinate_units=tuple(
            (param_id, ParameterUnit.FRACTION)
            for param_id in workflow.problem.controlled_ids
        ),
    )

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    initial = operation.evidence.states[0]
    attempts = native_mcmc_module._initialization_attempts_for_positions(
        plan,
        initial.positions,
    )
    assert attempts is not None
    assert any(attempt > 0 for attempt in attempts)
    for vector in initial.positions:
        resolved = workflow.parameterization.resolve(
            workflow.problem.lifecycle_frame(vector, workflow.parameterization)
        )
        values = tuple(
            resolved[workflow.local_ids[name]] for name in ("pa", *populations)
        )
        assert all(value >= 0.0 for value in values)
        assert math.fsum(values) == pytest.approx(1.0)


def test_mcmc_retries_public_population_proposals_outside_the_simplex() -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.4},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.5},
        fitted_name="PD",
    )
    accepted = _fit_workflow(workflow)
    plan = McmcPlan.for_accepted(
        accepted,
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
            (workflow.problem.controlled_ids[0], ParameterUnit.FRACTION),
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
        and "function_id='population_complement'" in evaluation.message
        for evaluation in initial_evaluations
    )

    operation = execute_mcmc_evidence(accepted, plan)

    assert operation.terminal is McmcOperationTerminal.COMPLETED
    assert operation.evidence is not None
    assert all(
        math.isfinite(value) for value in operation.evidence.states[0].log_densities
    )


def test_monte_carlo_reuses_direct_trf_simplex_coordinates() -> None:
    workflow = _build_workflow(
        start={"pb": 0.15, "pc": 0.25, "pd": 0.35},
        target={"pb": 0.1, "pc": 0.25, "pd": 0.35},
        fitted_name="PB",
    )
    accepted = _fit_workflow(workflow)
    size = workflow.engine.plan.observation_count
    dataset = ResamplingDatasetManifest(
        workflow.engine.plan,
        tuple(
            float(value) for value in accepted.evaluation_result.normalized_calculations
        ),
        tuple(False for _ in range(size)),
        tuple("G23" for _ in range(size)),
        tuple(f"nstate-observation-{index}" for index in range(size)),
    )
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=dataset,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        scheme=ResamplingScheme.MONTE_CARLO,
        replicate_count=2,
        replicate_structural_identities=("nstate-mc-0", "nstate-mc-1"),
        replicate_component_identities=(("nstate-trf-0",), ("nstate-trf-1",)),
        root_seed=772,
        output_scope=workflow.problem.commit_scope,
        output_units=("native",) * len(workflow.problem.commit_scope),
        minimum_successful_count=2,
        strategy=OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "100"),),
    )

    operation = execute_resampling_evidence(accepted, plan)

    assert operation.terminal is OperationTerminal.COMPLETED
    assert operation.evidence is not None
    assert operation.evidence.successful_count == 2
    for outcome in operation.evidence.outcomes:
        assert outcome.success is not None
        resolved = dict(outcome.success.resolved_items)
        populations = tuple(
            resolved[workflow.local_ids[name]] for name in ("pa", "pb", "pc", "pd")
        )
        assert all(value >= 0.0 for value in populations)
        assert sum(populations) == pytest.approx(1.0)


def test_bootstrap_near_simplex_face_never_commits_negative_pa() -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.45},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.499},
        fitted_name="PD",
    )
    accepted = _fit_workflow(workflow)
    size = workflow.engine.plan.observation_count
    dataset = ResamplingDatasetManifest(
        workflow.engine.plan,
        tuple(
            float(value) for value in accepted.evaluation_result.normalized_calculations
        ),
        tuple(False for _ in range(size)),
        tuple("G23" for _ in range(size)),
        tuple(f"nstate-bootstrap-{index}" for index in range(size)),
    )
    plan = ResamplingPlan.for_accepted(
        accepted,
        dataset=dataset,
        source_problem=workflow.problem,
        parameterization=workflow.parameterization,
        source_engine=workflow.engine,
        scheme=ResamplingScheme.BOOTSTRAP,
        replicate_count=2,
        replicate_structural_identities=("nstate-bs-0", "nstate-bs-1"),
        replicate_component_identities=(("nstate-bs-trf-0",), ("nstate-bs-trf-1",)),
        root_seed=777,
        output_scope=workflow.problem.commit_scope,
        output_units=("native",) * len(workflow.problem.commit_scope),
        minimum_successful_count=1,
        strategy=OptimizationStrategy.DIRECT_TRF,
        strategy_settings=(("objective_request_budget", "100"),),
    )
    before = workflow.session.analysis_values.snapshot()

    operation = execute_resampling_evidence(accepted, plan)

    assert operation.terminal is OperationTerminal.COMPLETED
    assert operation.evidence is not None
    assert workflow.session.analysis_values.snapshot() == before
    for outcome in operation.evidence.outcomes:
        if outcome.success is None:
            assert outcome.failure is not None
            continue
        resolved = dict(outcome.success.resolved_items)
        populations = tuple(
            resolved[workflow.local_ids[name]] for name in ("pa", "pb", "pc", "pd")
        )
        assert all(value >= 0.0 for value in populations)
        assert math.fsum(populations) == pytest.approx(1.0)


def test_post_fit_restart_round_trip_preserves_public_nstate_values_and_roles(
    tmp_path: Path,
) -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.25, "pd": 0.3},
        target={"pb": 0.13, "pc": 0.25, "pd": 0.3},
        fitted_name="PB",
        include_rates=True,
    )
    starting = workflow.session.analysis_values.snapshot()
    parameter_model = workflow.session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    args = Namespace(
        experiments=[],
        parameters=[],
        method=None,
        output=tmp_path / "Output",
        model="4st_fork",
        include=None,
        exclude=None,
        workers=1,
        native_threads="auto",
    )
    with patch("chemex.run_info._git_metadata", return_value=None):
        run_info = write_run_info(
            args,
            parameter_model=parameter_model,
            starting_values=starting,
            input_files=(),
            argv=("chemex", "fit"),
            working_directory=tmp_path,
        )
    accepted = _fit_workflow(workflow)
    committed = workflow.session.analysis_values.commit(
        dict(accepted.commit_items),
        expected=starting,
        scope=accepted.commit_scope,
    )

    run_info.publish_restart(committed)

    restart_path = tmp_path / "Output" / "run_info" / "restart.toml"
    restart = tomllib.loads(restart_path.read_text(encoding="utf-8"))
    restart_names = {name.split(",", 1)[0] for name in restart["GLOBAL"]}
    assert {"PB", "PC", "PD", "KEX_AB", "KEX_AC", "KEX_AD"} <= restart_names
    assert "PA" not in restart_names
    assert not {
        name
        for name in restart_names
        if name.startswith("K") and not name.startswith("KEX_")
    }

    continued = _build_workflow(
        start={"pb": 0.2, "pc": 0.25, "pd": 0.3},
        target={"pb": 0.13, "pc": 0.25, "pd": 0.3},
        fitted_name="PB",
        model_name="4st_fork",
        defaults=read_defaults([restart_path]),
    )
    continued_model = continued.session.parameter_factory.sealed_parameter_model
    assert continued_model is not None
    all_ids = {definition.param_id for definition in continued_model.definitions}
    ordinary = continued.session.compile_parameterization(Method(), all_ids)
    resolved = ordinary.resolve(
        ordinary.frame_from_snapshot(continued.session.analysis_values.snapshot())
    )
    for name in ("pb", "pc", "pd", "kex_ab", "kex_ac", "kex_ad"):
        param_id = continued.local_ids[name]
        assert resolved[param_id] == committed[workflow.local_ids[name]]
    for name in ("kex_ab", "kex_ac", "kex_ad"):
        assert ordinary.role(continued.local_ids[name]) is ParameterRole.FIT
    for name in ("pa", "kab", "kba", "kac", "kca", "kad", "kda"):
        assert ordinary.role(continued.local_ids[name]) is ParameterRole.DERIVED


def test_directional_rates_and_pa_receive_coupled_deterministic_uncertainty() -> None:
    workflow = _build_workflow(
        start={"pb": 0.15, "pc": 0.25, "pd": 0.35},
        target={"pb": 0.1, "pc": 0.25, "pd": 0.35},
        fitted_name="PB",
        include_rates=True,
    )
    accepted = _fit_workflow(workflow)
    uncertainty = _derive_uncertainty(workflow, accepted)
    standard_errors = {
        name: uncertainty.parameter(workflow.local_ids[name]).standard_error
        for name in ("pb", "pa", "kab", "kba")
    }

    assert all(value is not None for value in standard_errors.values())
    pb_error = standard_errors["pb"]
    assert pb_error is not None
    assert standard_errors["pa"] == pytest.approx(pb_error)
    assert standard_errors["kab"] == pytest.approx(500.0 * pb_error)
    assert standard_errors["kba"] == pytest.approx(500.0 * pb_error)


def test_pa_zero_qualifies_symmetric_deterministic_uncertainty_as_boundary() -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.4},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.5},
        fitted_name="PD",
        include_rates=True,
    )
    accepted = _fit_workflow(workflow)
    uncertainty = _derive_uncertainty(workflow, accepted)
    conclusion = uncertainty.parameter(workflow.local_ids["pd"])

    assert conclusion is not None
    assert conclusion.standard_error is not None
    assert conclusion.unavailable_kind is None
    assert conclusion.boundary_warning


@pytest.mark.parametrize("interrupt_second_block", (False, True))
def test_block_recovery_retains_pa_face_warning(
    interrupt_second_block: bool,
) -> None:
    workflow = _build_workflow(
        start={"pb": 0.2, "pc": 0.3, "pd": 0.4},
        target={"pb": 0.2, "pc": 0.3, "pd": 0.5},
        fitted_name="PD",
        temperatures=(25.0, 30.0),
    )
    accepted = _fit_workflow(workflow)
    original_svd = uncertainty_module.svd
    call_count = 0

    def rank_deficient_root_then_maybe_interrupt(*args: Any, **kwargs: Any):
        nonlocal call_count
        call_count += 1
        if interrupt_second_block and call_count == 3:
            raise KeyboardInterrupt
        left, singular, right = original_svd(*args, **kwargs)
        if call_count == 1:
            singular = np.array(singular, copy=True)
            singular[-1] = 0.0
        return left, singular, right

    with patch(
        "chemex.optimize.uncertainty.svd",
        side_effect=rank_deficient_root_then_maybe_interrupt,
    ):
        uncertainty = _derive_uncertainty(workflow, accepted)

    operation = uncertainty.block_operation
    assert operation is not None
    expected_completed = 1 if interrupt_second_block else 2
    assert len(operation.completed_blocks) == expected_completed
    warned = operation.simple_bound_warning_ids
    if operation.complete_evidence is not None:
        assert operation.complete_evidence.simple_bound_warning_ids == warned
    for local_ids in workflow.component_local_ids[:expected_completed]:
        assert local_ids["pd"] in warned
