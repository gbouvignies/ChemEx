"""Behavioral qualification for exact native grouped Direct TRF (#594)."""

from __future__ import annotations

import dataclasses
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    CancellationToken,
    DirectTrfConstructionError,
    DirectTrfInvocation,
    OptimizationProblem,
    canonical_chi_square,
    commit_accepted_fit,
    execute_direct_trf,
)
from chemex.optimize.grouped_direct_trf import (
    ComponentDisposition,
    FitDecomposition,
    GroupedDirectTrfInvocation,
    GroupedDirectTrfTerminal,
    execute_direct_trf_components,
    execute_grouped_direct_trf,
    materialize_grouped_direct_trf,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


def _grouped_problem(
    method: Method | None = None,
) -> tuple[
    AnalysisSession,
    Experiments,
    ActiveParameterization,
    EvaluationEngine,
    OptimizationProblem,
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    selected_method = read_methods([METHOD])["DEFAULT"] if method is None else method
    parameterization = session.compile_parameterization(
        selected_method,
        experiments.param_ids,
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
    return session, experiments, parameterization, engine, problem


def test_exact_decomposition_is_deterministic_and_distinguishes_shared_coordinates() -> (
    None
):
    _session, _experiments, parameterization, engine, problem = _grouped_problem()

    first = FitDecomposition.from_root(problem, parameterization, engine)
    second = FitDecomposition.from_root(problem, parameterization, engine)

    expected = (
        "__R1A_A_G2N_H_800_0MHZ",
        "__R1A_A_H3N_H_800_0MHZ",
        "__R1A_A_K4N_H_800_0MHZ",
        "__R1A_A_L6N_H_800_0MHZ",
        "__R1A_A_S5N_H_800_0MHZ",
    )
    assert tuple(component.controlled_ids[0] for component in first.components) == (
        expected
    )
    assert tuple(component.root_profile_indices for component in first.components) == (
        (0,),
        (1,),
        (2,),
        (4,),
        (3,),
    )
    assert first.constant_profile_indices == ()
    assert first.non_objective_profile_indices == ()
    assert first.partition_proof.root_plan_identity == engine.plan.identity
    assert tuple(
        record[2] for record in first.partition_proof.profile_records
    ) == tuple((param_id,) for param_id in problem.controlled_ids)
    assert tuple(
        record[2] for record in first.partition_proof.component_records
    ) == tuple(component.root_profile_indices for component in first.components)
    assert all(
        component.problem.source_snapshot is problem.source_snapshot
        for component in first.components
    )
    assert all(
        component.problem.derivation is not None
        and component.problem.derivation.root_problem_identity == problem.identity
        and component.problem.derivation.component_identity == component.identity
        and component.problem.derivation.projected_plan_identity
        == component.problem.evaluation_plan_identity
        and not component.problem.acceptance_authority
        for component in first.components
    )
    assert all(
        tuple(param_id for param_id, _value in component.problem.held_items)
        == tuple(
            param_id
            for param_id, _value in problem.independent_items
            if param_id not in component.controlled_ids
        )
        for component in first.components
    )
    assert first.identity == second.identity
    assert tuple(component.identity for component in first.components) == tuple(
        component.identity for component in second.components
    )
    child = first.components[0]
    child_invocation = DirectTrfInvocation.for_problem(
        child.problem,
        objective_request_budget=10,
    )
    with pytest.raises(DirectTrfConstructionError, match="no acceptance authority"):
        execute_direct_trf(
            child.problem,
            child_invocation,
            parameterization,
            engine.project_profiles(child.root_profile_indices),
        )

    _session, _experiments, shared_parameterization, shared_engine, shared_problem = (
        _grouped_problem(Method(fit=["PB"], fix=["KEX_AB"]))
    )
    shared = FitDecomposition.from_root(
        shared_problem,
        shared_parameterization,
        shared_engine,
    )
    assert len(shared.components) == 1
    assert set(shared.components[0].controlled_ids) == {"__PB", *expected}
    assert shared.components[0].root_profile_indices == (0, 1, 2, 3, 4)

    profile_record = first.partition_proof.profile_records[0]
    first_path = profile_record[3][0]
    invalid_paths = (
        (first_path[0], first_path[1], (first_path[0], "__UNDECLARED", first_path[1])),
        *profile_record[3][1:],
    )
    invalid_proof = dataclasses.replace(
        first.partition_proof,
        profile_records=(
            (profile_record[0], profile_record[1], profile_record[2], invalid_paths),
            *first.partition_proof.profile_records[1:],
        ),
    )
    invalid_decomposition = dataclasses.replace(first, partition_proof=invalid_proof)
    invalid_invocation = GroupedDirectTrfInvocation.for_decomposition(
        invalid_decomposition,
        objective_request_budgets=(10,) * len(invalid_decomposition.components),
    )
    invalid_outcome = execute_grouped_direct_trf(
        problem,
        invalid_decomposition,
        invalid_invocation,
        parameterization,
        engine,
    )
    assert invalid_outcome.terminal is GroupedDirectTrfTerminal.EXECUTION_FAILURE
    assert all(
        component.disposition is ComponentDisposition.NOT_STARTED
        for component in invalid_outcome.components
    )


def test_all_masked_unscaled_profile_cannot_supply_a_control_dependency() -> None:
    session, experiments, parameterization, _engine, _problem = _grouped_problem()
    profile = next(iter(experiments)).profiles[0]
    profile.is_scaled = False
    profile.data.mask[:] = False
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )

    assert engine.plan.profiles[0].retained_observation_indices == ()
    with pytest.raises(
        DirectTrfConstructionError,
        match="lack a retained-objective dependency",
    ):
        FitDecomposition.from_root(problem, parameterization, engine)


def test_successful_components_compose_one_fresh_root_accepted_result_and_commit_once() -> (
    None
):
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()

    outcome = execute_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    assert all(
        component.disposition is ComponentDisposition.SUCCEEDED
        for component in outcome.components
    )
    assert all(
        not hasattr(component, "accepted_result") for component in outcome.components
    )
    assert outcome.accepted_result is not None
    assert outcome.commit_authority is not None
    accepted = outcome.accepted_result
    selected = {
        param_id: component.candidate.vector[index]
        for component in outcome.components
        if component.candidate is not None
        for index, param_id in enumerate(component.controlled_ids)
    }
    assert accepted.vector == tuple(
        selected[param_id] for param_id in problem.controlled_ids
    )
    assert accepted.chi_square == canonical_chi_square(
        accepted.evaluation_result.residuals
    )
    assert session.analysis_values.snapshot() == before

    receipt = commit_accepted_fit(
        accepted,
        outcome.commit_authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=session.analysis_values,
    )
    assert receipt.old_revision == before.revision
    assert receipt.new_revision == before.revision + 1
    committed = session.analysis_values.snapshot()
    with pytest.raises(
        DirectTrfConstructionError,
        match="exact live Direct TRF commit authority",
    ):
        commit_accepted_fit(
            accepted,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == committed


def test_fresh_root_validation_rejects_an_inconsistent_component_result() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
    )
    first = components[0]
    assert first.candidate is not None
    inconsistent_evaluation = dataclasses.replace(
        first.candidate.evaluation_result,
        normalized_calculations=np.asarray(
            first.candidate.evaluation_result.normalized_calculations
        )
        + 1.0,
    )
    inconsistent_candidate = dataclasses.replace(
        first.candidate,
        evaluation_result=inconsistent_evaluation,
    )
    inconsistent_components = (
        dataclasses.replace(first, candidate=inconsistent_candidate),
        *components[1:],
    )

    outcome = materialize_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        inconsistent_components,
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before

    inconsistent_objective = dataclasses.replace(
        first.candidate,
        chi_square=first.candidate.chi_square + 1.0,
    )
    objective_outcome = materialize_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        (dataclasses.replace(first, candidate=inconsistent_objective), *components[1:]),
    )
    assert objective_outcome.failure is not None
    assert objective_outcome.failure.category == "decomposition_projection_mismatch"
    assert objective_outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before

    out_of_bounds_candidate = dataclasses.replace(first.candidate, vector=(-1.0,))
    out_of_bounds = materialize_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        (
            dataclasses.replace(first, candidate=out_of_bounds_candidate),
            *components[1:],
        ),
    )
    assert out_of_bounds.failure is not None
    assert out_of_bounds.failure.category == "aggregate_feasibility_failure"
    assert out_of_bounds.accepted_result is None
    assert session.analysis_values.snapshot() == before


@pytest.mark.parametrize(
    "provenance_field",
    (
        "problem_identity",
        "invocation_identity",
        "execution_identity",
        "evaluator_compatibility_identity",
        "evaluation_identity",
    ),
)
def test_component_materialization_provenance_is_non_retargetable(
    provenance_field: str,
) -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
    )
    first = components[0]
    assert first.candidate is not None
    foreign_materialization = dataclasses.replace(
        first.candidate.materialization,
        **{provenance_field: f"foreign-{provenance_field}"},
    )
    foreign_candidate = dataclasses.replace(
        first.candidate,
        materialization=foreign_materialization,
    )

    outcome = materialize_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        (dataclasses.replace(first, candidate=foreign_candidate), *components[1:]),
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE
    assert outcome.failure is not None
    assert outcome.failure.category == "decomposition_projection_mismatch"
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_component_vector_cannot_be_retargeted_to_transposed_controlled_ids() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem(
        Method(fit=["PB"], fix=["KEX_AB"]),
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,),
    )
    before = session.analysis_values.snapshot()

    def successful_backend(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64)
        assert np.isfinite(candidate).all()
        return _backend_result(candidate, fun(candidate), success=True)

    with patch("chemex.optimize.direct_trf.least_squares", successful_backend):
        components = execute_direct_trf_components(
            decomposition,
            invocation,
            parameterization,
            engine,
        )
    assert len(components) == 1
    component = components[0]
    assert component.candidate is not None
    first, second, *remaining = component.controlled_ids
    transposed_ownership = (second, first, *remaining)

    outcome = materialize_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        (dataclasses.replace(component, controlled_ids=transposed_ownership),),
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE
    assert outcome.failure is not None
    assert outcome.failure.category == "aggregate_assignment_failure"
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def _backend_result(
    candidate: Array,
    residuals: Array,
    *,
    success: bool,
) -> object:
    return SimpleNamespace(
        x=candidate,
        fun=residuals,
        status=1 if success else 0,
        success=success,
        message="converged" if success else "required component did not converge",
        nfev=1,
        njev=0,
        cost=0.5 * float(np.dot(residuals, residuals)),
        optimality=0.0 if success else 1.0,
        active_mask=np.zeros_like(candidate, dtype=np.int64),
    )


def test_required_component_failure_collects_other_outcomes_and_prevents_aggregate() -> (
    None
):
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    invocation_count = 0

    def fail_first(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        nonlocal invocation_count
        invocation_count += 1
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        return _backend_result(
            candidate,
            residuals,
            success=invocation_count != 1,
        )

    with patch("chemex.optimize.direct_trf.least_squares", fail_first):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert invocation_count == len(decomposition.components)
    assert outcome.terminal is GroupedDirectTrfTerminal.COMPONENT_FAILURE
    assert outcome.components[0].disposition is ComponentDisposition.FAILED
    assert all(
        component.disposition is ComponentDisposition.SUCCEEDED
        for component in outcome.components[1:]
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_component_contract_failure_stops_later_components() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()

    with patch(
        "chemex.optimize.direct_trf.least_squares",
        return_value=SimpleNamespace(),
    ):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.EXECUTION_FAILURE
    assert tuple(component.disposition for component in outcome.components) == (
        ComponentDisposition.EXECUTION_FAILURE,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_component_materialization_failure_retains_typed_detail() -> None:
    session, experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = cast("Callable[..., Array]", pulse_type.calculate)
    fail_materialization = False

    def controlled_calculate(*args: object, **kwargs: object) -> Array:
        if fail_materialization:
            raise RuntimeError("grouped component materialization failed")
        return original(*args, **kwargs)

    def converge_then_arm(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        nonlocal fail_materialization
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        fail_materialization = True
        return _backend_result(candidate, residuals, success=True)

    with (
        patch.object(pulse_type, "calculate", controlled_calculate),
        patch("chemex.optimize.direct_trf.least_squares", converge_then_arm),
    ):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.EXECUTION_FAILURE
    assert outcome.components[0].failure is not None
    assert outcome.components[0].failure.category == "kernel_exception"
    assert "grouped component materialization failed" in (
        outcome.components[0].failure.message
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_interruption_stops_later_components_without_exposing_partial_authority() -> (
    None
):
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    invocation_count = 0

    def interrupt_second(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        nonlocal invocation_count
        invocation_count += 1
        if invocation_count == 2:
            raise KeyboardInterrupt
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        return _backend_result(candidate, residuals, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", interrupt_second):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.INTERRUPTED
    assert tuple(component.disposition for component in outcome.components) == (
        ComponentDisposition.SUCCEEDED,
        ComponentDisposition.INTERRUPTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_projection_interruption_gives_every_component_a_disposition() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()

    with patch.object(engine, "project_profiles", side_effect=KeyboardInterrupt):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.INTERRUPTED
    assert all(
        component.disposition is ComponentDisposition.NOT_STARTED
        for component in outcome.components
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_cancellation_stops_after_in_flight_component_and_commits_nothing() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(10,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    token = CancellationToken()

    def cancel_first(
        fun: Callable[[Array], Array],
        x0: Array,
        **_settings: object,
    ) -> object:
        candidate = np.asarray(x0, dtype=np.float64) * 0.9
        residuals = fun(candidate)
        token.cancel()
        return _backend_result(candidate, residuals, success=True)

    with patch("chemex.optimize.direct_trf.least_squares", cancel_first):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            cancellation=token,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.CANCELLED
    assert tuple(component.disposition for component in outcome.components) == (
        ComponentDisposition.CANCELLED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
    )
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_cancellation_during_first_context_projection_stops_validation() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
    )
    token = CancellationToken()
    projection_count = 0
    project_profiles = engine.project_profiles

    def cancel_during_first_context_projection(
        profile_indices: tuple[int, ...],
    ) -> EvaluationEngine:
        nonlocal projection_count
        projection_count += 1
        projected = project_profiles(profile_indices)
        token.cancel()
        return projected

    with (
        patch.object(
            engine,
            "project_profiles",
            side_effect=cancel_during_first_context_projection,
        ),
        patch(
            "chemex.optimize.grouped_direct_trf._aggregate_vector",
            side_effect=AssertionError("aggregate reconstruction started"),
        ) as aggregate_vector,
        patch.object(
            engine,
            "new_evaluator",
            wraps=engine.new_evaluator,
        ) as new_root_evaluator,
    ):
        outcome = materialize_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            components,
            cancellation=token,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.CANCELLED
    assert projection_count == 1
    aggregate_vector.assert_not_called()
    new_root_evaluator.assert_not_called()
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before


def test_cancellation_during_first_evidence_comparison_stops_reconstruction() -> None:
    session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    before = session.analysis_values.snapshot()
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
    )
    token = CancellationToken()
    projection_count = 0
    comparison_count = 0
    context_projection_count = len(decomposition.components)
    project_profiles = engine.project_profiles

    def record_projection(
        profile_indices: tuple[int, ...],
    ) -> EvaluationEngine:
        nonlocal projection_count
        projection_count += 1
        return project_profiles(profile_indices)

    def cancel_during_first_evidence_comparison(residuals: Array) -> float:
        nonlocal comparison_count
        comparison_count += 1
        result = canonical_chi_square(residuals)
        token.cancel()
        return result

    with (
        patch.object(
            engine,
            "project_profiles",
            side_effect=record_projection,
        ),
        patch(
            "chemex.optimize.grouped_direct_trf.canonical_chi_square",
            side_effect=cancel_during_first_evidence_comparison,
        ),
        patch(
            "chemex.optimize.grouped_direct_trf._aggregate_vector",
            side_effect=AssertionError("aggregate reconstruction started"),
        ) as aggregate_vector,
        patch.object(
            engine,
            "new_evaluator",
            wraps=engine.new_evaluator,
        ) as new_root_evaluator,
    ):
        outcome = materialize_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            components,
            cancellation=token,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.CANCELLED
    assert projection_count == context_projection_count + 1
    assert comparison_count == 1
    aggregate_vector.assert_not_called()
    new_root_evaluator.assert_not_called()
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
