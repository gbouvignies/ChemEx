"""Behavioral qualification for exact native grouped Direct TRF (#594)."""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest

import chemex.optimize.grouped_direct_trf as grouped_direct_trf_owner
import chemex.parameters.feasible_coordinates as feasible_coordinates_owner
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.optimize import direct_trf as direct_trf_owner
from chemex.optimize.direct_trf import (
    AffineHalfSpace,
    CancellationToken,
    ComponentProblemDerivation,
    DirectTrfConstructionError,
    DirectTrfInvocation,
    OptimizationProblem,
    ProblemDerivation,
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
from chemex.optimize.progress import (
    FitProgressContext,
    ProgressEvent,
    ProgressPhase,
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
    if method is None:
        plan = read_method_plan([METHOD])
        parameterization = session.compile_parameterization_from_actions(
            plan.effective_role_actions()["DEFAULT"], experiments.param_ids
        )
    else:
        parameterization = session.compile_parameterization(
            method,
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

    derivations: list[ProblemDerivation] = []
    original_derive_child = OptimizationProblem.derive_grouped_component

    def track_derive_child(
        self: OptimizationProblem,
        *,
        controlled_ids: tuple[str, ...],
        start: tuple[float, ...],
        derivation: ComponentProblemDerivation,
    ) -> OptimizationProblem:
        derivations.append(derivation)
        return original_derive_child(
            self,
            controlled_ids=controlled_ids,
            start=start,
            derivation=derivation,
        )

    with patch.object(
        OptimizationProblem,
        "derive_grouped_component",
        track_derive_child,
    ):
        first = FitDecomposition.from_root(problem, parameterization, engine)
        second = FitDecomposition.from_root(problem, parameterization, engine)

    assert len(derivations) == 2 * len(first.components)
    assert all(
        isinstance(derivation, ComponentProblemDerivation) for derivation in derivations
    )

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
    root_chart = problem.feasible_coordinates
    assert root_chart is not None
    for component in first.components:
        frame = root_chart.frame_with_updates(
            dict(zip(component.controlled_ids, component.problem.start, strict=True))
        )
        freshly_compiled = feasible_coordinates_owner.compile_feasible_coordinates(
            parameterization,
            frame,
            component.controlled_ids,
            component.problem.lower_bounds,
            component.problem.upper_bounds,
        )
        assert component.problem.feasible_coordinates == freshly_compiled
        assert component.problem.feasible_coordinates is not None
        assert (
            component.problem.feasible_coordinates.identity == freshly_compiled.identity
        )
        assert (
            component.problem.feasible_coordinates.solver_domain
            == freshly_compiled.solver_domain
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
    assert (
        shared.components[0].problem.feasible_coordinates
        is shared_problem.feasible_coordinates
    )

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


def test_grouped_execution_does_not_rederive_a_supplied_decomposition() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )

    with patch.object(
        FitDecomposition,
        "from_root",
        side_effect=AssertionError("grouped execution rederived its decomposition"),
    ):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED


def test_grouped_validation_rejects_a_self_consistent_forged_component_identity() -> (
    None
):
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    original = decomposition.components[0]
    assert isinstance(original.problem.derivation, ComponentProblemDerivation)
    forged_identity = "forged-component-identity"
    forged_derivation = dataclasses.replace(
        original.problem.derivation,
        component_identity=forged_identity,
    )
    forged_problem = dataclasses.replace(
        original.problem,
        derivation=forged_derivation,
    )
    forged_component = dataclasses.replace(
        original,
        problem=forged_problem,
        identity=forged_identity,
    )
    forged_components = (forged_component, *decomposition.components[1:])
    first_record = decomposition.partition_proof.component_records[0]
    forged_proof = dataclasses.replace(
        decomposition.partition_proof,
        component_records=(
            (forged_identity, *first_record[1:]),
            *decomposition.partition_proof.component_records[1:],
        ),
    )
    forged_decomposition = dataclasses.replace(
        decomposition,
        components=forged_components,
        partition_proof=forged_proof,
    )
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        forged_decomposition,
        objective_request_budgets=(10,) * len(forged_components),
    )

    outcome = execute_grouped_direct_trf(
        problem,
        forged_decomposition,
        invocation,
        parameterization,
        engine,
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.EXECUTION_FAILURE
    assert all(
        component.disposition is ComponentDisposition.NOT_STARTED
        for component in outcome.components
    )


@pytest.mark.parametrize(
    "mutation",
    (
        "missing_chart",
        "altered_held_frame",
        "altered_block",
        "changed_start",
        "changed_policy",
    ),
)
def test_grouped_validation_rejects_forged_component_feasibility(
    mutation: str,
) -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    original = decomposition.components[0]
    chart = original.problem.feasible_coordinates
    assert chart is not None
    if mutation == "missing_chart":
        forged_problem = dataclasses.replace(
            original.problem,
            feasible_coordinates=None,
        )
    elif mutation == "altered_held_frame":
        held_id, held_value = original.problem.held_items[0]
        forged_chart = dataclasses.replace(
            chart,
            base_frame=chart.base_frame.with_updates({held_id: held_value + 0.25}),
        )
        forged_problem = dataclasses.replace(
            original.problem,
            feasible_coordinates=forged_chart,
        )
    elif mutation == "altered_block":
        block = chart.blocks[0]
        altered_block = dataclasses.replace(
            block,
            off_diagonal_ids=(
                (*block.off_diagonal_ids[0][:2], block.diagonal_ids[0]),
                *block.off_diagonal_ids[1:],
            ),
        )
        forged_problem = dataclasses.replace(
            original.problem,
            feasible_coordinates=dataclasses.replace(
                chart,
                blocks=(altered_block, *chart.blocks[1:]),
            ),
        )
    elif mutation == "changed_policy":
        derivation = original.problem.derivation
        assert isinstance(derivation, ComponentProblemDerivation)
        forged_problem = dataclasses.replace(
            original.problem,
            derivation=dataclasses.replace(
                derivation,
                projection_policy="forged-projection-policy",
            ),
        )
    else:
        forged_problem = dataclasses.replace(
            original.problem,
            start=(original.problem.start[0] + 0.01, *original.problem.start[1:]),
        )
    forged_component = dataclasses.replace(original, problem=forged_problem)
    forged_decomposition = dataclasses.replace(
        decomposition,
        components=(forged_component, *decomposition.components[1:]),
    )
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        forged_decomposition,
        objective_request_budgets=(10,) * len(forged_decomposition.components),
    )

    outcome = execute_grouped_direct_trf(
        problem,
        forged_decomposition,
        invocation,
        parameterization,
        engine,
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.EXECUTION_FAILURE
    assert all(
        component.disposition is ComponentDisposition.NOT_STARTED
        for component in outcome.components
    )


def test_standalone_materialization_rejects_a_foreign_component_plan() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    original = decomposition.components[0]
    original_profile = original.evaluation_plan.profiles[0]
    altered_profile = dataclasses.replace(
        original_profile,
        experimental_values=(
            original_profile.experimental_values[0] + 1.0,
            *original_profile.experimental_values[1:],
        ),
    )
    altered_plan = dataclasses.replace(
        original.evaluation_plan,
        profiles=(altered_profile, *original.evaluation_plan.profiles[1:]),
    )
    derivation = original.problem.derivation
    assert isinstance(derivation, ComponentProblemDerivation)
    forged_derivation = dataclasses.replace(
        derivation,
        projected_plan_identity=altered_plan.identity,
    )
    forged_problem = dataclasses.replace(
        original.problem,
        evaluation_plan_identity=altered_plan.identity,
        derivation=forged_derivation,
    )
    forged_identity = grouped_direct_trf_owner._component_identity(
        problem,
        engine,
        original.controlled_ids,
        original.root_profile_indices,
        altered_plan.identity,
        original.problem.lower_bounds,
        original.problem.upper_bounds,
    )
    forged_derivation = dataclasses.replace(
        forged_derivation,
        component_identity=forged_identity,
    )
    forged_problem = dataclasses.replace(
        forged_problem,
        derivation=forged_derivation,
    )
    forged_component = dataclasses.replace(
        original,
        problem=forged_problem,
        evaluation_plan=altered_plan,
        identity=forged_identity,
    )
    first_record = decomposition.partition_proof.component_records[0]
    forged_proof = dataclasses.replace(
        decomposition.partition_proof,
        component_records=(
            (forged_identity, *first_record[1:]),
            *decomposition.partition_proof.component_records[1:],
        ),
    )
    forged_decomposition = dataclasses.replace(
        decomposition,
        components=(forged_component, *decomposition.components[1:]),
        partition_proof=forged_proof,
    )
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        forged_decomposition,
        objective_request_budgets=(80,) * len(forged_decomposition.components),
    )

    outcome = materialize_grouped_direct_trf(
        problem,
        forged_decomposition,
        invocation,
        parameterization,
        engine,
        (),
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE


def test_grouped_validation_does_not_rebuild_projected_plans() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    project_plan = engine.project_plan
    project_profiles = engine.project_profiles

    with (
        patch.object(engine, "project_plan", wraps=project_plan) as plan_calls,
        patch.object(
            engine,
            "project_profiles",
            wraps=project_profiles,
        ) as engine_calls,
    ):
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    assert plan_calls.call_count == len(decomposition.components)
    assert engine_calls.call_count == len(decomposition.components)


def test_combined_grouped_execution_validates_context_once() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )

    with patch.object(
        grouped_direct_trf_owner,
        "_validate_grouped_context",
        wraps=grouped_direct_trf_owner._validate_grouped_context,
    ) as validation_calls:
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    validation_calls.assert_called_once()


def test_grouped_decomposition_rejects_a_derived_component_root() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    component = decomposition.components[0]
    child_engine = engine.project_profiles(component.root_profile_indices)

    with (
        patch.object(
            grouped_direct_trf_owner,
            "_profile_dependencies",
            side_effect=AssertionError("non-root decomposition inspected profiles"),
        ),
        pytest.raises(DirectTrfConstructionError, match="complete root"),
    ):
        FitDecomposition.from_root(
            component.problem,
            parameterization,
            child_engine,
        )


def test_grouped_component_projection_does_not_recompile_feasibility() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()

    with (
        patch.object(
            feasible_coordinates_owner,
            "compile_feasible_coordinates",
            side_effect=AssertionError("component feasibility was recompiled"),
        ),
        patch.object(
            problem.feasible_coordinates.__class__,
            "derive",
            side_effect=AssertionError("component feasibility used general derivation"),
        ),
        patch.object(
            engine,
            "project_profiles",
            side_effect=AssertionError("decomposition compiled child engines"),
        ),
    ):
        decomposition = FitDecomposition.from_root(problem, parameterization, engine)

    assert len(decomposition.components) == len(problem.controlled_ids)


def test_general_component_derivation_recompiles_feasibility_at_root_start() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    root_feasible = problem.feasible_coordinates
    assert root_feasible is not None
    derivation = ComponentProblemDerivation(
        problem.identity,
        problem.affine_feasibility_identity,
        "resampling-complete-scope",
        "native-resampling-complete-scope-v1",
        engine.plan.identity,
        problem.controlled_ids,
        problem.held_items,
    )
    derive_calls = 0
    original_derive = type(root_feasible).derive

    def track_derive(*args, **kwargs):
        nonlocal derive_calls
        derive_calls += 1
        return original_derive(*args, **kwargs)

    with patch.object(type(root_feasible), "derive", track_derive):
        child = problem.derive_child(
            controlled_ids=problem.controlled_ids,
            start=problem.start,
            derivation=derivation,
        )

    assert derive_calls == 1
    assert child.feasible_coordinates is not root_feasible


def test_grouped_exact_domain_rejects_projected_chart_as_forged_root() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    component = decomposition.components[0]
    assert component.problem.feasible_coordinates is not None
    forged_root = dataclasses.replace(component.problem, derivation=None)
    derivation = component.problem.derivation
    assert isinstance(derivation, ComponentProblemDerivation)

    with pytest.raises(
        feasible_coordinates_owner.FeasibleCoordinateConstructionError,
        match="compiled root chart",
    ):
        forged_root.derive_grouped_component(
            controlled_ids=forged_root.controlled_ids,
            start=forged_root.start,
            derivation=derivation,
        )


def test_grouped_components_preserve_root_affine_feasibility() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    controlled_id = problem.controlled_ids[0]
    coefficients = tuple(
        1.0 if param_id == controlled_id else 0.0
        for param_id, _value in problem.independent_items
    )
    restriction = AffineHalfSpace(
        "grouped-root-upper",
        coefficients,
        problem.start[0] + 1.0,
    )
    constrained = dataclasses.replace(
        problem,
        affine_half_spaces=(restriction,),
    )
    decomposition = FitDecomposition.from_root(
        constrained,
        parameterization,
        engine,
    )
    assert decomposition.components
    for component in decomposition.components:
        assert component.problem.affine_half_spaces == (restriction,)
        with pytest.raises(DirectTrfConstructionError, match="feasibility"):
            dataclasses.replace(component.problem, affine_half_spaces=())


def test_relaxation_feasibility_hyperedge_merges_separate_profile_controls() -> None:
    components = grouped_direct_trf_owner._ordered_component_controls(
        ("r2", "eta", "pb"),
        (frozenset(("r2",)), frozenset(("eta",)), frozenset(("pb",))),
        (frozenset(("r2", "eta")),),
    )

    assert components == (("pb",), ("r2", "eta"))


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

    with patch.object(
        direct_trf_owner,
        "materialize_root_candidate",
        wraps=direct_trf_owner.materialize_root_candidate,
    ) as shared_materialization:
        outcome = execute_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    assert shared_materialization.call_count == 1
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
    retained = accepted.final_residual_jacobian
    assert retained is not None
    assert retained.source.value == "fit-partition-composition"
    assert retained.controlled_ids == problem.controlled_ids
    base = np.asarray(accepted.evaluation_result.residuals, dtype=np.float64)
    independent_columns: list[Array] = []
    evaluator = engine.new_evaluator()
    for index, value in enumerate(accepted.vector):
        step = math.sqrt(np.finfo(np.float64).eps) * max(1.0, abs(value))
        candidate = list(accepted.vector)
        candidate[index] += step
        frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            problem.lifecycle_frame(tuple(candidate), parameterization),
        )
        evaluated = evaluator.evaluate(frame)
        assert isinstance(evaluated, EvaluationResult)
        independent_columns.append(
            (np.asarray(evaluated.residuals, dtype=np.float64) - base) / step
        )
    np.testing.assert_allclose(
        retained.matrix,
        np.column_stack(independent_columns),
        rtol=2.0e-6,
        atol=2.0e-8,
    )
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
        match="exact live fit commit authority",
    ):
        commit_accepted_fit(
            accepted,
            outcome.commit_authority,
            problem=problem,
            parameterization=parameterization,
            analysis_values=session.analysis_values,
        )
    assert session.analysis_values.snapshot() == committed


def test_root_jacobian_composition_converts_each_component_matrix_once() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem(
        Method(fit=["PB"], fix=["KEX_AB"])
    )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    assert len(decomposition.components) == 1
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(100,),
    )
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
    )
    retained_matrices = tuple(
        component.final_residual_jacobian.matrix
        for component in components
        if component.final_residual_jacobian is not None
    )
    assert len(retained_matrices) == 1
    conversions: list[object] = []
    numpy_owner = grouped_direct_trf_owner.np

    class NumpyConversionProbe:
        def __getattr__(self, name: str) -> object:
            return getattr(numpy_owner, name)

        def asarray(
            self,
            value: object,
            dtype: type[np.float64] | np.dtype[np.float64] | None = None,
        ) -> Array:
            if any(value is matrix for matrix in retained_matrices):
                conversions.append(value)
            return numpy_owner.asarray(value, dtype=dtype)

    with patch.object(grouped_direct_trf_owner, "np", NumpyConversionProbe()):
        outcome = materialize_grouped_direct_trf(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            components,
        )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    assert conversions == list(retained_matrices)


def test_grouped_direct_progress_identifies_each_component_in_order() -> None:
    _session, _experiments, parameterization, engine, problem = _grouped_problem()
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    invocation = GroupedDirectTrfInvocation.for_decomposition(
        decomposition,
        objective_request_budgets=(80,) * len(decomposition.components),
    )
    observed: list[tuple[FitProgressContext, ProgressEvent]] = []

    outcome = execute_grouped_direct_trf(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        progress_observer=lambda context, event: observed.append((context, event)),
    )

    assert outcome.terminal is GroupedDirectTrfTerminal.ACCEPTED
    started = [item for item in observed if item[1].phase is ProgressPhase.STARTED]
    terminated = [
        item for item in observed if item[1].phase is ProgressPhase.TERMINATED
    ]
    assert [context.component_ordinal for context, _event in started] == list(
        range(1, len(decomposition.components) + 1)
    )
    assert all(
        context.component_total == len(decomposition.components)
        for context, _event in started
    )
    assert [context.controlled_ids for context, _event in started] == [
        component.controlled_ids for component in decomposition.components
    ]
    assert [context for context, _event in terminated] == [
        context for context, _event in started
    ]


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
            dataclasses.replace(
                first,
                candidate=out_of_bounds_candidate,
                final_residual_jacobian=None,
            ),
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
        (
            dataclasses.replace(
                component,
                controlled_ids=transposed_ownership,
                final_residual_jacobian=None,
            ),
        ),
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
        jac=np.zeros((residuals.size, candidate.size), dtype=np.float64),
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


def test_projection_interruption_marks_the_in_flight_component() -> None:
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
    assert tuple(component.disposition for component in outcome.components) == (
        ComponentDisposition.INTERRUPTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
        ComponentDisposition.NOT_STARTED,
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


def test_cancellation_during_first_carried_projection_stops_validation() -> None:
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
    materialized_projection = (
        grouped_direct_trf_owner._materialized_component_projection
    )

    def cancel_during_first_carried_projection(*args, **kwargs):
        nonlocal projection_count
        projection_count += 1
        projected = materialized_projection(*args, **kwargs)
        token.cancel()
        return projected

    with (
        patch.object(
            grouped_direct_trf_owner,
            "_materialized_component_projection",
            side_effect=cancel_during_first_carried_projection,
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
    comparison_count = 0

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
            side_effect=AssertionError("materialization rebuilt a child engine"),
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
    assert comparison_count == 1
    aggregate_vector.assert_not_called()
    new_root_evaluator.assert_not_called()
    assert outcome.accepted_result is None
    assert session.analysis_values.snapshot() == before
