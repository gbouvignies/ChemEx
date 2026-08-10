"""Behavioral qualification tests for native uncertainty evidence (#598).

The public seam is one immutable accepted fit plus its exact native problem,
parameterization, and evaluator lineage.  Evidence derivation must not mutate
the fit, its commit authority, or session-owned Analysis Values.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    DirectTrfInvocation,
    OptimizationProblem,
    execute_direct_trf,
)
from chemex.optimize.uncertainty import (
    ClaimState,
    FunctionAnalyticPartialCapability,
    FunctionFiniteDifferenceCapability,
    OperationTerminal,
    ParameterUnit,
    ResidualVarianceScaling,
    UncertaintyPolicy,
    compile_constraint_linearization_capabilities,
    derive_uncertainty_evidence,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
DCEST_EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
DCEST_PARAMETERS = (
    ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"
)


def _accepted_relaxation_fit() -> tuple[
    AnalysisSession,
    ActiveParameterization,
    EvaluationEngine,
    OptimizationProblem,
    AcceptedFitResult,
]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(
            include=[SpinSystem.from_name("G2N-HN")],
            exclude=None,
        ),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    method = read_methods([METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    invocation = DirectTrfInvocation.for_problem(
        problem,
        objective_request_budget=80,
    )
    outcome = execute_direct_trf(
        problem,
        invocation,
        parameterization,
        engine,
    )
    assert outcome.accepted_result is not None
    return session, parameterization, engine, problem, outcome.accepted_result


def _qualification_policy(controlled_id: str) -> UncertaintyPolicy:
    # These are explicit qualification inputs, not inferred production defaults.
    return UncertaintyPolicy(
        calibration_identity="issue-598-relaxation-qualification-v1",
        numerical_compatibility_requirement=(
            "binary64-scipy-economical-svd-fixed-pairwise-v1"
        ),
        coordinate_scales=((controlled_id, 1.0),),
        coordinate_units=((controlled_id, ParameterUnit.RATE_PER_SECOND),),
        residual_variance_scaling=(
            ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE
        ),
        relative_step_tolerance=1.0e-4,
        roundoff_multiplier=64.0,
        smaller_step_extent=8,
        larger_step_extent=8,
        svd_driver="gesdd",
        rank_absolute_tolerance=0.0,
        rank_relative_tolerance=1.0e-12,
        weak_relative_tolerance=1.0e-6,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="box-domain-no-affine-restrictions-v1",
    )


def _hd_problem(
    method: Method,
) -> tuple[ActiveParameterization, EvaluationEngine, OptimizationProblem]:
    session = AnalysisSession.create()
    session.set_model("2st_hd")
    experiments = build_experiments(
        [DCEST_EXPERIMENT],
        Selection(include=[SpinSystem.from_name("1N")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([DCEST_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    return parameterization, engine, problem


def _accepted_hd_fit() -> tuple[
    ActiveParameterization,
    EvaluationEngine,
    OptimizationProblem,
    AcceptedFitResult,
]:
    method = Method(
        fit=("D2O",),
        fix=("CS_A", "DW_AB", "KDH", "PHI", "R1_A", "R2_A", "R2_B"),
    )
    parameterization, engine, problem = _hd_problem(method)
    outcome = execute_direct_trf(
        problem,
        DirectTrfInvocation.for_problem(problem, objective_request_budget=80),
        parameterization,
        engine,
    )
    assert outcome.accepted_result is not None
    return parameterization, engine, problem, outcome.accepted_result


def _accepted_reference_anchor(
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    problem: OptimizationProblem,
    vector: tuple[float, ...] | None = None,
) -> AcceptedFitResult:
    anchor = problem.start if vector is None else vector
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame(anchor, parameterization),
    )
    evaluation = engine.new_evaluator().evaluate(frame)
    assert isinstance(evaluation, EvaluationResult)
    chi_square = float(np.dot(evaluation.residuals, evaluation.residuals))
    return AcceptedFitResult(
        "qualification-anchor-occurrence",
        problem.identity,
        "qualification-anchor-invocation",
        "qualification-anchor-execution",
        "qualification-anchor-materialization",
        problem.parameterization_identity,
        problem.evaluator_parameterization_identity,
        problem.source_snapshot.occurrence_identity,
        problem.source_snapshot.revision,
        problem.controlled_ids,
        anchor,
        chi_square,
        evaluation,
        problem.controlled_ids,
        tuple(zip(problem.controlled_ids, anchor, strict=True)),
        "qualification-anchor-origin",
    )


def test_accepted_fit_yields_typed_covariance_and_joint_constraint_propagation() -> (
    None
):
    session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    before = session.analysis_values.snapshot()
    controlled_id = problem.controlled_ids[0]
    derived_id = "__R1A_B_G2N_H_800_0MHZ"
    compiled = compile_constraint_linearization_capabilities(
        parameterization,
        (controlled_id, derived_id),
        (),
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        constrained_scope=(controlled_id, derived_id),
        constrained_units=(
            (controlled_id, ParameterUnit.RATE_PER_SECOND),
            (derived_id, ParameterUnit.RATE_PER_SECOND),
        ),
        constrained_scales=((controlled_id, 1.0), (derived_id, 1.0)),
        compiled_constraint_linearization=compiled,
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.accepted_result_identity == accepted.identity
    assert evidence.residual_jacobian is not None
    assert evidence.residual_jacobian.controlled_ids == (controlled_id,)
    assert evidence.residual_jacobian.residual_count == 7
    assert evidence.residual_jacobian.coordinate_count == 1
    assert evidence.residual_jacobian.complete_reliable
    # Frozen from an independently evaluated five-point stencil (h=1e-4).  The
    # tolerance covers its O(h^4) truncation plus the production two-scale
    # stencil's binary64 roundoff without masking a scientifically material shift.
    expected_jacobian = np.asarray(
        (
            6.877571026489022,
            3.791023675708857,
            1.0915054108627373,
            -1.2596720503206598,
            -3.297643939993577,
            -5.054307903075824,
            -6.558611676446162,
        )
    )
    np.testing.assert_allclose(
        np.asarray(evidence.residual_jacobian.matrix)[:, 0],
        expected_jacobian,
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    assert evidence.covariance is not None
    factor = np.asarray(evidence.covariance.factor)
    for column in range(factor.shape[1]):
        pivot = int(np.argmax(np.abs(factor[:, column])))
        assert factor[pivot, column] >= 0.0
    assert evidence.covariance.controlled_ids == (controlled_id,)
    assert evidence.covariance.retained_residual_count == 7
    assert evidence.covariance.controlled_coordinate_count == 1
    assert evidence.covariance.profiled_normalization_count == 1
    assert evidence.covariance.nominal_residual_degrees_of_freedom == 5
    assert evidence.covariance.rank == 1
    assert evidence.covariance.factorization == "scaled-svd-gram-v1"
    np.testing.assert_allclose(
        evidence.covariance.covariance[0][0],
        0.018371558813164147,
        rtol=4.0e-7,
        atol=2.0e-12,
    )
    assert evidence.marginal_errors is not None
    assert evidence.marginal_errors.entries[0].value is not None
    assert evidence.correlations is not None
    assert evidence.correlations.entries[0][0].value == 1.0
    assert evidence.correlations.units == (ParameterUnit.RATE_PER_SECOND,)
    assert any(
        item.name == "CORRELATION_SCOPE_REPORTABILITY"
        for item in evidence.correlations.claims
    )

    assert evidence.constraint_jacobian is not None
    assert evidence.constraint_jacobian.output_ids == (controlled_id, derived_id)
    assert evidence.constraint_jacobian.matrix == ((1.0,), (1.0,))
    assert (
        evidence.constraint_jacobian.claims[0].name
        == "CONSTRAINT_OUTPUT_LINEARIZATION_COMPLETE"
    )
    assert evidence.constrained_propagation is not None
    propagated = evidence.constrained_propagation.covariance
    source_variance = evidence.covariance.covariance[0][0]
    np.testing.assert_allclose(
        propagated,
        ((source_variance, source_variance),) * 2,
        rtol=0.0,
        atol=0.0,
    )
    assert (
        evidence.constrained_propagation.claim("OUTPUT_RANK_DEFICIENCY_EXPECTED")
        is ClaimState.SATISFIED
    )
    assert any(
        item.name.startswith("SOURCE_COVARIANCE::")
        for item in evidence.constrained_propagation.claims
    )
    assert evidence.residual_jacobian.evaluation_plan_identity == engine.plan.identity
    assert evidence.covariance.source_jacobian_identity == (
        evidence.residual_jacobian.identity
    )
    assert evidence.constrained_propagation.source_covariance_identity == (
        evidence.covariance.identity
    )
    assert evidence.resolved_environment_identity == ("local-qualification-environment")
    assert all(
        item.resolved_environment_identity == "local-qualification-environment"
        for item in evidence.operations
    )
    assert evidence.rank_diagnostic is not None
    assert evidence.rank_diagnostic.identifiable_projector == ((1.0,),)
    assert evidence.rank_diagnostic.null_projector == ((0.0,),)
    assert evidence.rank_diagnostic.weak_projector == ((0.0,),)
    record = evidence.to_record()
    json.dumps(record, allow_nan=False, sort_keys=True)
    payload = record["payload"]
    assert isinstance(payload, dict)
    typed_payload = cast("dict[str, object]", payload)
    assert typed_payload["resolved_environment_identity"] == (
        "local-qualification-environment"
    )
    other_environment = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        constrained_scope=(controlled_id, derived_id),
        constrained_units=(
            (controlled_id, ParameterUnit.RATE_PER_SECOND),
            (derived_id, ParameterUnit.RATE_PER_SECOND),
        ),
        constrained_scales=((controlled_id, 1.0), (derived_id, 1.0)),
        compiled_constraint_linearization=compiled,
        resolved_environment_identity="different-local-environment",
    )
    assert other_environment.identity == evidence.identity
    assert tuple(item.identity for item in other_environment.operations) == tuple(
        item.identity for item in evidence.operations
    )
    assert tuple(
        item.occurrence_identity for item in other_environment.operations
    ) != tuple(item.occurrence_identity for item in evidence.operations)
    assert session.analysis_values.snapshot() == before


def test_scientific_constraint_requires_and_uses_declared_numerical_capability() -> (
    None
):
    parameterization, engine, problem, accepted = _accepted_hd_fit()
    controlled_id = problem.controlled_ids[0]
    output_id = "__PA_1_0_1000"
    base_policy = dataclasses.replace(
        _qualification_policy(controlled_id),
        coordinate_units=((controlled_id, ParameterUnit.FRACTION),),
    )

    with pytest.raises(ValueError, match="lacks capabilities"):
        compile_constraint_linearization_capabilities(
            parameterization,
            (output_id,),
            (),
        )

    capability = FunctionFiniteDifferenceCapability(
        function_id="pop_2st",
        component="pa",
        argument_scales=(10.0, 10.0),
        output_scale=1.0,
    )
    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=base_policy,
        constrained_scope=(output_id,),
        constrained_units=((output_id, ParameterUnit.FRACTION),),
        constrained_scales=((output_id, 1.0),),
        compiled_constraint_linearization=(
            compile_constraint_linearization_capabilities(
                parameterization,
                (output_id,),
                (capability,),
            )
        ),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.constraint_jacobian is not None
    d2o = accepted.evaluation_result.resolved_values[controlled_id]
    phi = accepted.evaluation_result.resolved_values["__PHI_1"]
    expected = -phi / (1.0 + d2o * (phi - 1.0)) ** 2
    np.testing.assert_allclose(
        evidence.constraint_jacobian.matrix[0][0],
        expected,
        rtol=2.0e-6,
        atol=2.0e-8,
    )
    assert len(evidence.constraint_jacobian.function_partial_diagnostics) == 2
    assert all(
        item.method == "centered_two_scale_numerical"
        for item in evidence.constraint_jacobian.function_partial_diagnostics
    )

    unreliable = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=dataclasses.replace(
            base_policy,
            relative_step_tolerance=0.0,
            roundoff_multiplier=0.0,
            smaller_step_extent=0,
            larger_step_extent=0,
        ),
        constrained_scope=(output_id,),
        constrained_units=((output_id, ParameterUnit.FRACTION),),
        constrained_scales=((output_id, 1.0),),
        compiled_constraint_linearization=(
            compile_constraint_linearization_capabilities(
                parameterization,
                (output_id,),
                (capability,),
            )
        ),
        resolved_environment_identity="local-qualification-environment",
    )
    assert any(
        item.category == "unreliable_or_nondifferentiable_function_partial"
        and item.trajectory_fingerprint is not None
        for item in unreliable.failures
    )

    analytic = FunctionAnalyticPartialCapability(
        function_id="pop_2st",
        component="pa",
        implementation_identity="pop-2st-pa-analytic-partials-v1",
        partials=(
            lambda kab, kba: -kba / (kab + kba) ** 2,
            lambda kab, kba: kab / (kab + kba) ** 2,
        ),
    )
    analytic_evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=base_policy,
        constrained_scope=(output_id,),
        constrained_units=((output_id, ParameterUnit.FRACTION),),
        constrained_scales=((output_id, 1.0),),
        compiled_constraint_linearization=(
            compile_constraint_linearization_capabilities(
                parameterization,
                (output_id,),
                (analytic,),
            )
        ),
        resolved_environment_identity="local-qualification-environment",
    )
    assert analytic_evidence.constraint_jacobian is not None
    np.testing.assert_allclose(
        analytic_evidence.constraint_jacobian.matrix[0][0],
        expected,
        rtol=5.0e-15,
        atol=5.0e-15,
    )
    assert all(
        item.method == "versioned_analytic_partial"
        for item in analytic_evidence.constraint_jacobian.function_partial_diagnostics
    )


def test_multivariate_scaled_svd_correlation_and_joint_propagation_reference() -> None:
    method = Method(
        fit=("D2O", "KDH"),
        fix=("CS_A", "DW_AB", "PHI", "R1_A", "R2_A", "R2_B"),
    )
    parameterization, engine, problem = _hd_problem(method)
    accepted = _accepted_reference_anchor(parameterization, engine, problem)
    d2o_id, kdh_id = problem.controlled_ids
    kab_id = "__KAB_1_0_1000"
    policy = dataclasses.replace(
        _qualification_policy(d2o_id),
        coordinate_scales=((d2o_id, 0.1), (kdh_id, 10.0)),
        coordinate_units=(
            (d2o_id, ParameterUnit.FRACTION),
            (kdh_id, ParameterUnit.RATE_PER_SECOND),
        ),
        residual_variance_scaling=(
            ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
        ),
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        constrained_scope=(d2o_id, kdh_id, kab_id),
        constrained_units=(
            (d2o_id, ParameterUnit.FRACTION),
            (kdh_id, ParameterUnit.RATE_PER_SECOND),
            (kab_id, ParameterUnit.RATE_PER_SECOND),
        ),
        constrained_scales=((d2o_id, 0.1), (kdh_id, 10.0), (kab_id, 10.0)),
        compiled_constraint_linearization=(
            compile_constraint_linearization_capabilities(
                parameterization,
                (d2o_id, kdh_id, kab_id),
                (),
            )
        ),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.failures == ()
    assert evidence.residual_jacobian is not None
    np.testing.assert_allclose(
        np.asarray(evidence.residual_jacobian.matrix)[:3],
        (
            (666.6712446920574, 1.910664913855726),
            (0.27370067627634853, -0.01357184843345749),
            (-5.879327716305852, -0.3296351577155292),
        ),
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    assert evidence.rank_diagnostic is not None
    np.testing.assert_allclose(
        evidence.rank_diagnostic.singular_values,
        (154.19754083595208, 24.21131174143699),
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    np.testing.assert_allclose(
        evidence.rank_diagnostic.identifiable_projector,
        np.eye(2),
        rtol=0.0,
        atol=3.0e-16,
    )
    assert evidence.covariance is not None
    expected_covariance = np.asarray(
        (
            (4.078404361040173e-06, -0.0006890732177197917),
            (-0.0006890732177197917, 0.13401556849617918),
        )
    )
    np.testing.assert_allclose(
        evidence.covariance.covariance,
        expected_covariance,
        rtol=4.0e-7,
        atol=2.0e-12,
    )
    assert evidence.correlations is not None
    np.testing.assert_allclose(
        np.asarray(
            [
                [cast("float", entry.value) for entry in row]
                for row in evidence.correlations.entries
            ]
        ),
        ((1.0, -0.9320572795396814), (-0.9320572795396814, 1.0)),
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    assert evidence.constraint_jacobian is not None
    d2o, kdh = accepted.vector
    phi = accepted.evaluation_result.resolved_values["__PHI_1"]
    expected_gradient = np.asarray(((1.0, 0.0), (0.0, 1.0), (kdh * phi, d2o * phi)))
    np.testing.assert_allclose(
        evidence.constraint_jacobian.matrix,
        expected_gradient,
        rtol=0.0,
        atol=0.0,
    )
    assert evidence.constrained_propagation is not None
    np.testing.assert_allclose(
        evidence.constrained_propagation.covariance,
        expected_gradient @ expected_covariance @ expected_gradient.T,
        rtol=4.0e-15,
        atol=2.0e-17,
    )


def test_rank_deficiency_is_diagnostic_and_never_manufactures_uncertainty() -> None:
    session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    policy = _qualification_policy(controlled_id)
    rank_rejecting_policy = dataclasses.replace(
        policy,
        rank_absolute_tolerance=1.0e100,
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=rank_rejecting_policy,
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is not None
    assert evidence.rank_diagnostic is not None
    assert evidence.rank_diagnostic.rank == 0
    assert evidence.covariance is None
    assert evidence.marginal_errors is None
    assert evidence.correlations is None
    assert [(item.stage, item.category) for item in evidence.failures] == [
        ("covariance", "rank_deficient")
    ]
    assert all(
        operation.artifact_identity is None
        for operation in evidence.operations
        if operation.stage == "covariance"
    )
    assert session.analysis_values.snapshot().revision == 0


def test_stencil_exhaustion_retains_trajectory_and_each_search_extent() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    policy = dataclasses.replace(
        _qualification_policy(problem.controlled_ids[0]),
        relative_step_tolerance=0.0,
        roundoff_multiplier=0.0,
        smaller_step_extent=0,
        larger_step_extent=0,
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="local-qualification-environment",
    )

    failure = evidence.failures[0]
    assert failure.category == "exhausted_declared_step_search_extent"
    assert failure.trajectory_fingerprint is not None
    assert failure.termination_details == (
        "exhausted_declared_smaller_step_extent",
        "exhausted_declared_larger_step_extent",
    )


def test_boundary_covariance_remains_diagnostic_but_reporting_fails_closed() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    boundary_problem = dataclasses.replace(
        problem,
        lower_bounds=(accepted.vector[0],),
    )
    boundary_accepted = dataclasses.replace(
        accepted,
        problem_identity=boundary_problem.identity,
    )

    evidence = derive_uncertainty_evidence(
        boundary_accepted,
        problem=boundary_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is not None
    assert evidence.residual_jacobian.columns[0].orientation == "one_sided_positive"
    assert evidence.covariance is not None
    assert evidence.covariance.claim("INTERIOR_POINT") is ClaimState.VIOLATED
    assert evidence.covariance.claim("BOUNDARY_SEPARATION") is ClaimState.VIOLATED
    assert evidence.covariance.claim("USABLE_LOCAL_COVARIANCE") is ClaimState.VIOLATED
    assert evidence.marginal_errors is not None
    assert evidence.marginal_errors.entries[0].value is not None
    assert not evidence.marginal_errors.scope_reportable
    assert evidence.correlations is not None
    assert evidence.correlations.entries[0][0].value == 1.0
    assert not evidence.correlations.scope_reportable


def test_unrepresentable_centered_interval_falls_back_to_inward_stencil() -> None:
    _session, parameterization, engine, problem, _accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    interior = (float(np.nextafter(0.0, 1.0)),)
    anchor = _accepted_reference_anchor(
        parameterization,
        engine,
        problem,
        vector=interior,
    )

    evidence = derive_uncertainty_evidence(
        anchor,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is not None
    assert evidence.residual_jacobian.columns[0].orientation == "one_sided_positive"


def test_nonfinite_boundary_separation_is_indeterminate() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    finite_extreme_problem = dataclasses.replace(
        problem,
        upper_bounds=(float(np.finfo(np.float64).max),),
    )
    compatible_accepted = dataclasses.replace(
        accepted,
        problem_identity=finite_extreme_problem.identity,
    )

    evidence = derive_uncertainty_evidence(
        compatible_accepted,
        problem=finite_extreme_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.covariance is not None
    assert evidence.covariance.claim("BOUNDARY_SEPARATION") is ClaimState.INDETERMINATE


def test_non_finite_accepted_objective_yields_only_typed_failed_covariance() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    non_finite = dataclasses.replace(accepted, chi_square=float("inf"))

    evidence = derive_uncertainty_evidence(
        non_finite,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is not None
    assert evidence.covariance is None
    assert evidence.rank_diagnostic is not None
    assert [(item.stage, item.category) for item in evidence.failures] == [
        ("covariance", "non_finite_chi_square")
    ]
    assert evidence.marginal_errors is None
    assert evidence.correlations is None


def test_foreign_coordinate_scope_fails_before_any_numerical_evidence() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    policy = dataclasses.replace(
        _qualification_policy(problem.controlled_ids[0]),
        coordinate_scales=(("__FOREIGN", 1.0),),
        coordinate_units=(("__FOREIGN", ParameterUnit.RATE_PER_SECOND),),
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is None
    assert evidence.rank_diagnostic is None
    assert evidence.covariance is None
    assert evidence.constraint_jacobian is None
    assert evidence.failures[0].category == "coordinate_scale_scope_mismatch"
    assert len(evidence.operations) == 1


def test_held_constraint_output_is_explicitly_outside_uncertainty_basis() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        constrained_scope=("__PB",),
        constrained_units=(("__PB", ParameterUnit.FRACTION),),
        constrained_scales=(("__PB", 1.0),),
        compiled_constraint_linearization=(
            compile_constraint_linearization_capabilities(
                parameterization,
                ("__PB",),
                (),
            )
        ),
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.covariance is not None
    assert evidence.constraint_jacobian is None
    assert evidence.constrained_propagation is None
    assert [item.category for item in evidence.failures[-2:]] == [
        "held_independent_outside_uncertainty_basis",
        "source_artifact_unavailable",
    ]


def test_absolute_observation_uncertainties_do_not_apply_residual_scaling() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    policy = dataclasses.replace(
        _qualification_policy(problem.controlled_ids[0]),
        residual_variance_scaling=(
            ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
        ),
    )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.covariance is not None
    assert evidence.covariance.residual_variance_scale == 1.0
    assert evidence.covariance.covariance == (evidence.covariance.unscaled_covariance)
    assert (
        evidence.covariance.claim("RESIDUAL_VARIANCE_NONDEGENERACY")
        is ClaimState.NOT_APPLICABLE
    )


def test_malformed_non_finite_svd_output_is_a_typed_covariance_failure() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    malformed_svd = (
        np.ones((7, 1), dtype=np.float64),
        np.array([np.nan], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
    )

    with patch("chemex.optimize.uncertainty.svd", return_value=malformed_svd):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            resolved_environment_identity="local-qualification-environment",
        )

    assert evidence.residual_jacobian is not None
    assert evidence.covariance is None
    assert evidence.rank_diagnostic is None
    assert evidence.failures[-1].category == "invalid_covariance_arithmetic"


def test_negative_svd_spectrum_is_typed_invalid_evidence() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    malformed_svd = (
        np.ones((7, 1), dtype=np.float64),
        np.array([-1.0], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
    )

    with patch("chemex.optimize.uncertainty.svd", return_value=malformed_svd):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            resolved_environment_identity="local-qualification-environment",
        )

    assert evidence.covariance is None
    assert evidence.failures[-1].category == "invalid_singular_spectrum"


def test_positive_scale_covariance_factor_underflow_is_typed_failure() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    underflowing_svd = (
        np.ones((7, 1), dtype=np.float64),
        np.array([1.0e200], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
    )

    with patch("chemex.optimize.uncertainty.svd", return_value=underflowing_svd):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            resolved_environment_identity="local-qualification-environment",
        )

    assert evidence.covariance is None
    assert evidence.failures[-1].category == "covariance_factor_underflow"


def test_cancellation_freezes_terminal_operation_without_artifact() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    probe_count = 0

    def cancel_inside_residual_operation() -> OperationTerminal | None:
        nonlocal probe_count
        probe_count += 1
        return None if probe_count == 1 else OperationTerminal.CANCELLED

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        cancellation_probe=cancel_inside_residual_operation,
        resolved_environment_identity="local-qualification-environment",
    )

    assert evidence.residual_jacobian is None
    assert evidence.operations[0].terminal is OperationTerminal.CANCELLED
    assert evidence.operations[0].artifact_identity is None
    assert evidence.failures[0].category == "cancelled"
    assert len(evidence.operations) == 1
