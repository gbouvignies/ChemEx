"""Behavioral qualification tests for native uncertainty evidence (#598).

The public seam is one immutable accepted fit plus its exact native problem,
parameterization, and evaluator lineage.  Evidence derivation must not mutate
the fit, its commit authority, or session-owned Analysis Values.
"""

from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path
from typing import cast
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFrame,
    EvaluationResult,
)
from chemex.experiments.builder import build_experiments
from chemex.optimize import native_deterministic as native_deterministic_module
from chemex.optimize import uncertainty as uncertainty_module
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    AffineEquality,
    AffineHalfSpace,
    DirectTrfConstructionError,
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
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PARAMETERS = ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"
DCEST_EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
DCEST_PARAMETERS = (
    ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"
)


def _analytic_pa_kab(kab: float, kba: float) -> float:
    return -kba / (kab + kba) ** 2


def _analytic_pa_kab_alternative(kab: float, kba: float) -> float:
    return -(kba + 1.0) / (kab + kba + 1.0) ** 2


def _analytic_pa_kba(kab: float, kba: float) -> float:
    return kab / (kab + kba) ** 2


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
        singular_value_cluster_relative_tolerance=1.0e-10,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1",
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
    return AcceptedFitResult.for_qualification(
        occurrence_identity="qualification-anchor-occurrence",
        problem_identity=problem.identity,
        invocation_identity="qualification-anchor-invocation",
        execution_identity="qualification-anchor-execution",
        materialization_identity="qualification-anchor-materialization",
        parameterization_identity=problem.parameterization_identity,
        evaluator_parameterization_identity=problem.evaluator_parameterization_identity,
        source_occurrence_identity=problem.source_snapshot.occurrence_identity,
        source_revision=problem.source_snapshot.revision,
        controlled_ids=problem.controlled_ids,
        vector=anchor,
        chi_square=chi_square,
        evaluation_result=evaluation,
        commit_scope=problem.controlled_ids,
        commit_items=tuple(zip(problem.controlled_ids, anchor, strict=True)),
        origin_context_identity="qualification-anchor-origin",
    )


def _qualified_accepted_copy(
    accepted: AcceptedFitResult,
    *,
    occurrence_identity: str | None = None,
    problem_identity: str | None = None,
    chi_square: float | None = None,
) -> AcceptedFitResult:
    """Mint a distinct authoritative occurrence for isolated qualification."""
    return AcceptedFitResult.for_qualification(
        occurrence_identity=(
            accepted.occurrence_identity
            if occurrence_identity is None
            else occurrence_identity
        ),
        problem_identity=(
            accepted.problem_identity if problem_identity is None else problem_identity
        ),
        invocation_identity=accepted.invocation_identity,
        execution_identity=accepted.execution_identity,
        materialization_identity=accepted.materialization_identity,
        parameterization_identity=accepted.parameterization_identity,
        evaluator_parameterization_identity=accepted.evaluator_parameterization_identity,
        source_occurrence_identity=accepted.source_occurrence_identity,
        source_revision=accepted.source_revision,
        controlled_ids=accepted.controlled_ids,
        vector=accepted.vector,
        chi_square=accepted.chi_square if chi_square is None else chi_square,
        evaluation_result=accepted.evaluation_result,
        commit_scope=accepted.commit_scope,
        commit_items=accepted.commit_items,
        origin_context_identity=accepted.origin_context_identity,
    )


def _independent_five_point_jacobian(
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    problem: OptimizationProblem,
    vector: tuple[float, ...],
    steps: tuple[float, ...],
    *,
    forward_columns: frozenset[int] = frozenset(),
) -> Array:
    """Auditable five-point reference independent of uncertainty.py."""

    def residuals(candidate: tuple[float, ...]) -> Array:
        frame = EvaluationFrame.from_lifecycle_frame(
            parameterization,
            problem.lifecycle_frame(candidate, parameterization),
        )
        evaluated = engine.new_evaluator().evaluate(frame)
        assert isinstance(evaluated, EvaluationResult)
        return np.asarray(evaluated.residuals, dtype=np.float64)

    base = residuals(vector)
    columns: list[Array] = []
    for column, step in enumerate(steps):
        values: list[Array] = []
        if column in forward_columns:
            for multiplier in (1.0, 2.0, 3.0, 4.0):
                candidate = list(vector)
                candidate[column] += multiplier * step
                values.append(residuals(tuple(candidate)))
            derivative = (
                -25.0 * base
                + 48.0 * values[0]
                - 36.0 * values[1]
                + 16.0 * values[2]
                - 3.0 * values[3]
            ) / (12.0 * step)
        else:
            for multiplier in (-2.0, -1.0, 1.0, 2.0):
                candidate = list(vector)
                candidate[column] += multiplier * step
                values.append(residuals(tuple(candidate)))
            derivative = (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (
                12.0 * step
            )
        columns.append(derivative)
    return np.column_stack(columns)


def _independent_svd_covariance(
    jacobian: Array,
    coordinate_scales: tuple[float, ...],
    residual_variance_scale: float,
) -> tuple[Array, Array, Array]:
    """Direct closed-form SVD reference independent of uncertainty.py."""
    scales = np.asarray(coordinate_scales, dtype=np.float64)
    scaled = np.asarray(jacobian, dtype=np.float64) * scales[np.newaxis, :]
    _left, singular_values, right_transpose = np.linalg.svd(
        scaled,
        full_matrices=False,
    )
    right = right_transpose.T
    inverse_singular = 1.0 / singular_values
    unscaled_factor = scales[:, np.newaxis] * right * inverse_singular[np.newaxis, :]
    factor = math.sqrt(residual_variance_scale) * unscaled_factor
    return unscaled_factor @ unscaled_factor.T, factor, factor @ factor.T


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
    # Independent five-point O(h^4) derivative, followed below by the scalar
    # closed form C = (chi²/nu) / (J^T J).  Neither calls uncertainty.py.
    # The reference J is [6.8775710242, 3.7910236726, 1.0915054084,
    # -1.2596720492, -3.2976439424, -5.0543079057, -6.5586116745]; the
    # observed production/reference maxima are 3.16e-9 absolute and 2.23e-9
    # relative.  Its independently derived variance is 0.01837155881904635;
    # production differs by 5.89e-12, supporting the declared tolerances.
    expected_jacobian = _independent_five_point_jacobian(
        parameterization,
        engine,
        problem,
        accepted.vector,
        (1.0e-4,),
    )
    np.testing.assert_allclose(
        np.asarray(evidence.residual_jacobian.matrix)[:, 0],
        expected_jacobian[:, 0],
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    assert evidence.covariance is not None
    factor = np.asarray(evidence.covariance.factor)
    np.testing.assert_allclose(
        factor @ factor.T,
        evidence.covariance.covariance,
        rtol=0.0,
        atol=0.0,
    )
    assert evidence.covariance.controlled_ids == (controlled_id,)
    assert evidence.covariance.retained_residual_count == 7
    assert evidence.covariance.controlled_coordinate_count == 1
    assert evidence.covariance.profiled_normalization_count == 1
    assert evidence.covariance.nominal_residual_degrees_of_freedom == 5
    assert evidence.covariance.rank == 1
    assert evidence.covariance.factorization == ("direct-scaled-svd-factor-gram-v3")
    expected_variance = (accepted.chi_square / 5.0) / float(
        expected_jacobian[:, 0] @ expected_jacobian[:, 0]
    )
    np.testing.assert_allclose(
        evidence.covariance.covariance[0][0],
        expected_variance,
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
    reference_jacobian = _independent_five_point_jacobian(
        parameterization,
        engine,
        problem,
        accepted.vector,
        (1.0e-4, 1.0e-2),
    )
    np.testing.assert_allclose(
        np.asarray(evidence.residual_jacobian.matrix),
        reference_jacobian,
        # The independent five-point and production adaptive finite-difference
        # calculations differ by up to 4.674e-6 across the qualified local,
        # GitHub Actions Python 3.13, and Python 3.14 platform envelope.
        rtol=6.0e-6,
        atol=2.0e-9,
    )
    assert evidence.rank_diagnostic is not None
    scaled_reference = reference_jacobian * np.asarray((0.1, 10.0))[np.newaxis, :]
    expected_singular_values = np.linalg.svd(scaled_reference, compute_uv=False)
    # Independent expected singular values are [154.1975407925,
    # 24.2113113641].  The SVD comparison therefore measures the production
    # derivative error, not a copied production decomposition.
    np.testing.assert_allclose(
        evidence.rank_diagnostic.singular_values,
        expected_singular_values,
        rtol=2.0e-7,
        atol=2.0e-9,
    )
    np.testing.assert_allclose(
        evidence.rank_diagnostic.identifiable_projector,
        np.eye(2),
        rtol=0.0,
        atol=3.0e-16,
    )
    with pytest.raises(ValueError, match="classified subspaces"):
        dataclasses.replace(
            evidence.rank_diagnostic,
            identifiable_projector=((1.0, 0.0), (0.0, 0.0)),
            null_projector=((0.0, 0.0), (0.0, 1.0)),
        )
    assert evidence.covariance is not None
    expected_unscaled, _expected_factor, expected_covariance = (
        _independent_svd_covariance(
            reference_jacobian,
            (0.1, 10.0),
            1.0,
        )
    )
    # Direct SVD reference C = [[4.0784044539e-6, -6.8907323809e-4],
    # [-6.8907323809e-4, 1.3401557289e-1]].  Observed production/reference
    # covariance differences are 4.40e-9 absolute and 3.28e-8 relative, well
    # inside the 4e-7 tolerance inherited from the finite-difference error.
    np.testing.assert_allclose(
        evidence.covariance.covariance,
        expected_covariance,
        rtol=4.0e-7,
        atol=2.0e-12,
    )
    np.testing.assert_allclose(
        evidence.covariance.unscaled_covariance,
        expected_unscaled,
        rtol=4.0e-7,
        atol=2.0e-12,
    )
    assert evidence.correlations is not None
    expected_correlation = expected_covariance[0, 1] / math.sqrt(
        expected_covariance[0, 0] * expected_covariance[1, 1]
    )
    np.testing.assert_allclose(
        np.asarray(
            [
                [cast("float", entry.value) for entry in row]
                for row in evidence.correlations.entries
            ]
        ),
        ((1.0, expected_correlation), (expected_correlation, 1.0)),
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
        rtol=4.0e-7,
        atol=2.0e-12,
    )


def test_covariance_uses_direct_svd_not_normal_equations() -> None:
    jacobian = np.asarray(
        (
            (1.0, 1.0),
            (1.0, 1.0 + 1.0e-6),
            (1.0, 1.0 - 1.0e-6),
            (2.0, 2.0 + 0.5e-6),
        ),
        dtype=np.float64,
    )
    scales = (1.0e-3, 1.0e3)
    expected_unscaled, expected_factor, expected_covariance = (
        _independent_svd_covariance(jacobian, scales, 0.75)
    )
    scaled = jacobian * np.asarray(scales)[np.newaxis, :]
    _left, singular_values, right_transpose = np.linalg.svd(
        scaled,
        full_matrices=False,
    )
    actual_unscaled, actual_factor, actual_covariance = (
        uncertainty_module._canonical_covariance_reduction(
            singular_values,
            right_transpose,
            scales,
            0.75,
        )
    )
    np.testing.assert_allclose(
        actual_unscaled,
        expected_unscaled,
        rtol=2.0e-15,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.abs(actual_factor),
        np.abs(expected_factor),
        rtol=2.0e-15,
        atol=0.0,
    )
    np.testing.assert_allclose(
        actual_covariance,
        expected_covariance,
        rtol=2.0e-15,
        atol=0.0,
    )

    scale_matrix = np.diag(scales)
    normal_equation_covariance = 0.75 * (
        scale_matrix @ np.linalg.solve(scaled.T @ scaled, np.eye(2)) @ scale_matrix
    )
    relative_difference = np.max(
        np.abs((normal_equation_covariance - expected_covariance) / expected_covariance)
    )
    assert relative_difference > 1.0e-5


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
    boundary_accepted = _qualified_accepted_copy(
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
    independent_forward = _independent_five_point_jacobian(
        parameterization,
        engine,
        boundary_problem,
        boundary_accepted.vector,
        (1.0e-4,),
        forward_columns=frozenset({0}),
    )
    # The independent one-sided O(h^4) reference J is [6.8775710282,
    # 3.7910236723, 1.0915054092, -1.2596720491, -3.2976439423,
    # -5.0543079056, -6.5586116752].  Observed maxima are 2.02e-8 absolute
    # and 9.41e-9 relative, supporting the same 2e-7 / 2e-9 policy.
    np.testing.assert_allclose(
        np.asarray(evidence.residual_jacobian.matrix),
        independent_forward,
        rtol=2.0e-7,
        atol=2.0e-9,
    )
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
    forged_reportability_claims = tuple(
        dataclasses.replace(item, state=ClaimState.SATISFIED)
        if item.name == "MARGINAL_SCOPE_REPORTABILITY"
        else item
        for item in evidence.marginal_errors.claims
    )
    with pytest.raises(ValueError, match="reportability"):
        dataclasses.replace(
            evidence.marginal_errors,
            scope_reportable=True,
            claims=forged_reportability_claims,
        )


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
    compatible_accepted = _qualified_accepted_copy(
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
    non_finite = _qualified_accepted_copy(
        accepted,
        chi_square=float("inf"),
    )

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
    np.testing.assert_allclose(
        evidence.covariance.covariance,
        evidence.covariance.unscaled_covariance,
        rtol=3.0e-16,
        atol=0.0,
    )
    assert (
        evidence.covariance.claim("RESIDUAL_VARIANCE_NONDEGENERACY")
        is ClaimState.NOT_APPLICABLE
    )


def test_product_request_excludes_unsupported_scientific_function_propagation() -> None:
    parameterization, _engine, problem = _hd_problem(Method())

    request, unsupported = native_deterministic_module._product_uncertainty_request(
        problem,
        parameterization,
    )

    assert unsupported
    assert set(unsupported).isdisjoint(request.constrained_scope)
    assert request.policy.residual_variance_scaling is (
        ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
    )
    assert request.policy.coordinate_scales


def test_product_request_does_not_hide_structural_capability_errors() -> None:
    parameterization, _engine, problem = _hd_problem(Method())

    with (
        patch.object(
            native_deterministic_module,
            "compile_constraint_linearization_capabilities",
            side_effect=uncertainty_module.UncertaintyConstructionError("malformed"),
        ),
        pytest.raises(
            uncertainty_module.UncertaintyConstructionError, match="malformed"
        ),
    ):
        native_deterministic_module._product_uncertainty_request(
            problem,
            parameterization,
        )


def test_scaled_svd_preserves_full_rank_when_information_condition_overflows() -> None:
    policy = dataclasses.replace(
        _qualification_policy("a"),
        coordinate_scales=(("a", 1.0), ("b", 1.0)),
        coordinate_units=(
            ("a", ParameterUnit.UNSPECIFIED),
            ("b", ParameterUnit.UNSPECIFIED),
        ),
        rank_relative_tolerance=0.0,
    )

    with patch.object(
        uncertainty_module,
        "svd",
        return_value=(np.eye(2), np.array((1.0e200, 1.0)), np.eye(2)),
    ):
        analysis, failure = uncertainty_module._analyze_scaled_jacobian(
            np.eye(2),
            (1.0, 1.0),
            policy=policy,
            source_identity="overflowing-information-condition",
            cancellation_probe=None,
        )

    assert analysis is not None
    assert analysis.rank == 2
    assert analysis.jacobian_condition == pytest.approx(1.0e200)
    assert analysis.information_condition is not None
    assert math.isinf(analysis.information_condition)
    assert failure is not None
    assert failure.category == "non_finite_information_condition"


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


def test_exact_accepted_occurrence_is_not_reconstructible_or_mixable() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    reconstructed = dataclasses.replace(
        accepted, occurrence_identity="replacement-occurrence"
    )
    assert reconstructed.identity == accepted.identity
    with pytest.raises(ValueError, match="authoritative accepted occurrence"):
        derive_uncertainty_evidence(
            reconstructed,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            resolved_environment_identity="reconstructed-occurrence",
        )
    no_witness = dataclasses.replace(accepted, occurrence_witness=None)
    with pytest.raises(ValueError, match="authoritative accepted occurrence"):
        derive_uncertainty_evidence(
            no_witness,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            resolved_environment_identity="missing-occurrence-witness",
        )

    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="occurrence-integrity",
    )
    assert evidence.residual_jacobian is not None
    assert evidence.rank_diagnostic is not None
    assert evidence.covariance is not None
    reconstructed_anchor = dataclasses.replace(
        accepted,
        occurrence_identity="direct-artifact-reconstruction",
        occurrence_witness=None,
    )
    with pytest.raises(ValueError, match="exact accepted occurrence"):
        dataclasses.replace(
            evidence.residual_jacobian,
            accepted_anchor=reconstructed_anchor,
            accepted_occurrence_identity=reconstructed_anchor.occurrence_identity,
        )
    with pytest.raises(ValueError, match="exact accepted occurrence"):
        dataclasses.replace(
            evidence.residual_jacobian,
            accepted_occurrence_identity="foreign-occurrence",
        )
    with pytest.raises(ValueError, match="exact accepted occurrence"):
        dataclasses.replace(
            evidence.covariance,
            accepted_occurrence_identity="foreign-occurrence",
        )
    with pytest.raises(ValueError, match="exact accepted occurrence"):
        dataclasses.replace(
            evidence,
            accepted_occurrence_identity="foreign-occurrence",
            accepted_anchor=accepted,
        )
    with pytest.raises(ValueError, match="request"):
        dataclasses.replace(evidence, request_identity="foreign-request")
    with pytest.raises(ValueError, match="request"):
        dataclasses.replace(evidence, policy_identity="foreign-policy")

    equivalent_occurrence = _qualified_accepted_copy(
        accepted,
        occurrence_identity="equivalent-authoritative-occurrence",
    )
    assert equivalent_occurrence.identity == accepted.identity
    other_evidence = derive_uncertainty_evidence(
        equivalent_occurrence,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="other-equivalent-occurrence",
    )
    assert other_evidence.residual_jacobian is not None
    assert other_evidence.rank_diagnostic is not None
    assert other_evidence.covariance is not None
    with pytest.raises(ValueError, match="accepted fit occurrences"):
        dataclasses.replace(
            evidence,
            residual_jacobian=other_evidence.residual_jacobian,
        )
    with pytest.raises(ValueError, match="accepted fit occurrences"):
        dataclasses.replace(
            evidence,
            rank_diagnostic=other_evidence.rank_diagnostic,
        )
    with pytest.raises(ValueError, match="accepted fit occurrences"):
        dataclasses.replace(evidence, covariance=other_evidence.covariance)


def test_authoritative_artifacts_reject_inconsistent_reconstruction() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    derived_id = "__R1A_B_G2N_H_800_0MHZ"
    evidence = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        constrained_scope=(derived_id,),
        constrained_units=((derived_id, ParameterUnit.RATE_PER_SECOND),),
        constrained_scales=((derived_id, 1.0),),
        compiled_constraint_linearization=compile_constraint_linearization_capabilities(
            parameterization, (derived_id,), ()
        ),
        resolved_environment_identity="artifact-integrity",
    )
    jacobian = evidence.residual_jacobian
    rank = evidence.rank_diagnostic
    covariance = evidence.covariance
    constraint = evidence.constraint_jacobian
    propagation = evidence.constrained_propagation
    assert jacobian is not None
    assert rank is not None
    assert covariance is not None
    assert constraint is not None
    assert propagation is not None
    assert evidence.marginal_errors is not None

    replacements = (
        (jacobian, {"matrix": ((42.0,),) * 7}),
        (jacobian, {"constraint_program_identity": "foreign-program"}),
        (jacobian, {"evaluation_plan_identity": "foreign-plan"}),
        (jacobian, {"policy_identity": "foreign-policy"}),
        (jacobian, {"calibration_identity": "foreign-calibration"}),
        (
            jacobian,
            {"numerical_compatibility_requirement": "foreign-requirement"},
        ),
        (rank, {"source_jacobian_identity": "foreign-jacobian"}),
        (rank, {"singular_values": (42.0,)}),
        (rank, {"weak_threshold": 42.0}),
        (rank, {"rank": 0}),
        (rank, {"identifiable_projector": ((0.0,),)}),
        (rank, {"weak_projector": ((1.0,),)}),
        (covariance, {"source_jacobian_identity": "foreign-jacobian"}),
        (covariance, {"rank_diagnostic_identity": "foreign-svd"}),
        (covariance, {"coordinate_scales": (42.0,)}),
        (covariance, {"residual_variance_scale": 42.0}),
        (covariance, {"factor": ((0.0,),)}),
        (covariance, {"unscaled_covariance": ((42.0,),)}),
        (covariance, {"covariance": ((42.0,),)}),
        (constraint, {"accepted_evaluation_identity": "foreign-evaluation"}),
        (constraint, {"matrix": ((42.0,),)}),
        (constraint, {"output_ids": ("__FOREIGN",)}),
        (propagation, {"source_constraint_jacobian_identity": "foreign-G"}),
        (propagation, {"output_ids": ("__FOREIGN",)}),
        (propagation, {"covariance": ((42.0,),)}),
    )
    for artifact, changes in replacements:
        with pytest.raises(ValueError):
            dataclasses.replace(artifact, **changes)

    forged_column = dataclasses.replace(
        jacobian.columns[0],
        fine_estimate=tuple(value + 1.0 for value in jacobian.columns[0].fine_estimate),
    )
    forged_matrix = tuple((value,) for value in forged_column.fine_estimate)
    with pytest.raises(ValueError, match="replayed"):
        dataclasses.replace(
            jacobian,
            columns=(forged_column,),
            matrix=forged_matrix,
        )
    forged_subspace = dataclasses.replace(
        rank.subspaces[0],
        classification="isolated_null",
    )
    with pytest.raises(ValueError, match="classification"):
        dataclasses.replace(rank, subspaces=(forged_subspace,))

    altered_claims = tuple(
        dataclasses.replace(item, state=ClaimState.SATISFIED)
        if item.name == "USABLE_LOCAL_COVARIANCE"
        else item
        for item in covariance.claims
    )
    if altered_claims == covariance.claims:
        altered_claims = tuple(
            dataclasses.replace(item, state=ClaimState.VIOLATED)
            if item.name == "USABLE_LOCAL_COVARIANCE"
            else item
            for item in covariance.claims
        )
    with pytest.raises(ValueError, match="Covariance claims"):
        dataclasses.replace(covariance, claims=altered_claims)
    forged_variance_claims = tuple(
        dataclasses.replace(item, detail="forged residual variance claim")
        if item.name == "RESIDUAL_VARIANCE_NONDEGENERACY"
        else item
        for item in covariance.claims
    )
    with pytest.raises(ValueError, match="Covariance claims"):
        dataclasses.replace(covariance, claims=forged_variance_claims)
    forged_rank_claims = tuple(
        dataclasses.replace(item, detail="forged propagated rank claim")
        if item.name == "OUTPUT_RANK_DEFICIENCY_EXPECTED"
        else item
        for item in propagation.claims
    )
    with pytest.raises(ValueError, match="propagation claims"):
        dataclasses.replace(propagation, claims=forged_rank_claims)
    with pytest.raises(ValueError):
        dataclasses.replace(
            evidence.marginal_errors,
            entries=(
                dataclasses.replace(evidence.marginal_errors.entries[0], value=42.0),
            ),
        )


def _full_frame_coefficients(
    problem: OptimizationProblem,
    coefficients: dict[str, float],
) -> tuple[float, ...]:
    return tuple(
        coefficients.get(param_id, 0.0)
        for param_id, _value in problem.independent_items
    )


@pytest.mark.parametrize(
    ("zeta", "expected"),
    (
        (4.0, ClaimState.SATISFIED),
        (2.0, ClaimState.VIOLATED),
        (0.0, ClaimState.VIOLATED),
    ),
)
def test_affine_directional_separation_uses_exact_full_frame_slack(
    zeta: float,
    expected: ClaimState,
) -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    policy = _qualification_policy(problem.controlled_ids[0])
    baseline = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="affine-baseline",
    )
    assert baseline.covariance is not None
    standard_error = math.sqrt(baseline.covariance.covariance[0][0])
    restriction = AffineHalfSpace(
        "controlled-upper",
        _full_frame_coefficients(problem, {problem.controlled_ids[0]: 1.0}),
        accepted.vector[0] + zeta * standard_error,
    )
    affine_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(restriction,),
    )
    affine_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=affine_problem.identity,
        occurrence_identity=f"affine-{zeta.hex()}",
    )
    evidence = derive_uncertainty_evidence(
        affine_accepted,
        problem=affine_problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="affine-case",
    )
    assert evidence.covariance is not None
    assert evidence.covariance.claim("AFFINE_FEASIBILITY") is expected
    assert evidence.covariance.usable is (expected is ClaimState.SATISFIED)


def test_active_affine_upper_uses_inward_stencil_without_off_feasible_calls() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    restriction = AffineHalfSpace(
        "active-stencil-upper",
        _full_frame_coefficients(problem, {controlled_id: 1.0}),
        accepted.vector[0],
    )
    constrained_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(restriction,),
    )
    constrained_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=constrained_problem.identity,
        occurrence_identity="active-stencil-upper-occurrence",
    )
    original_evaluate = BoundEvaluator.evaluate
    evaluated_values: list[float] = []

    def guarded_evaluate(
        bound: BoundEvaluator,
        frame: EvaluationFrame,
    ) -> object:
        items = cast("list[list[object]]", frame.to_record()["items"])
        encoded = cast(
            "str",
            next(item[1] for item in items if item[0] == controlled_id),
        )
        value = float.fromhex(encoded)
        assert value <= accepted.vector[0]
        evaluated_values.append(value)
        return original_evaluate(bound, frame)

    with patch.object(BoundEvaluator, "evaluate", guarded_evaluate):
        evidence = derive_uncertainty_evidence(
            constrained_accepted,
            problem=constrained_problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(controlled_id),
            resolved_environment_identity="active-affine-stencil",
        )

    assert evidence.residual_jacobian is not None
    column = evidence.residual_jacobian.columns[0]
    assert column.orientation == "one_sided_negative"
    assert column.feasible_displacement_interval[1] == 0.0
    assert column.upper_feasibility_limiters == ("affine:active-stencil-upper",)
    assert evaluated_values


def test_affine_lower_inactive_box_intersection_and_coupled_held_stencils() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    lower_restriction = AffineHalfSpace(
        "active-stencil-lower",
        _full_frame_coefficients(problem, {controlled_id: -1.0}),
        -accepted.vector[0],
    )
    lower_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(lower_restriction,),
    )
    lower_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=lower_problem.identity,
        occurrence_identity="active-stencil-lower-occurrence",
    )
    lower_evidence = derive_uncertainty_evidence(
        lower_accepted,
        problem=lower_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="active-affine-lower-stencil",
    )
    assert lower_evidence.residual_jacobian is not None
    lower_column = lower_evidence.residual_jacobian.columns[0]
    assert lower_column.orientation == "one_sided_positive"
    assert lower_column.feasible_displacement_interval[0] == 0.0
    assert lower_column.lower_feasibility_limiters == ("affine:active-stencil-lower",)

    inactive = AffineHalfSpace(
        "inactive-stencil-upper",
        _full_frame_coefficients(problem, {controlled_id: 1.0}),
        accepted.vector[0] + 0.5,
    )
    inactive_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(inactive,),
    )
    interval = inactive_problem.coordinate_line_feasibility(accepted.vector, 0)
    assert interval.minimum_displacement == problem.lower_bounds[0] - accepted.vector[0]
    assert interval.maximum_displacement == 0.5
    assert interval.upper_limiters == ("affine:inactive-stencil-upper",)
    inactive_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=inactive_problem.identity,
        occurrence_identity="inactive-stencil-occurrence",
    )
    inactive_evidence = derive_uncertainty_evidence(
        inactive_accepted,
        problem=inactive_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="inactive-affine-stencil",
    )
    assert inactive_evidence.residual_jacobian is not None
    assert inactive_evidence.residual_jacobian.columns[0].orientation == "centered"

    held_id, held_value = problem.held_items[0]
    coupled = AffineHalfSpace(
        "coupled-held-stencil-upper",
        _full_frame_coefficients(
            problem,
            {controlled_id: 1.0, held_id: 1.0},
        ),
        accepted.vector[0] + held_value,
    )
    coupled_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(coupled,),
    )
    coupled_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=coupled_problem.identity,
        occurrence_identity="coupled-held-stencil-occurrence",
    )
    coupled_evidence = derive_uncertainty_evidence(
        coupled_accepted,
        problem=coupled_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="coupled-held-affine-stencil",
    )
    assert coupled_evidence.residual_jacobian is not None
    coupled_column = coupled_evidence.residual_jacobian.columns[0]
    assert coupled_column.orientation == "one_sided_negative"
    assert coupled_column.upper_feasibility_limiters == (
        "affine:coupled-held-stencil-upper",
    )


def test_affine_equality_has_no_single_coordinate_stencil() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    equality = AffineEquality(
        "fixed-stencil-coordinate",
        _full_frame_coefficients(problem, {controlled_id: 1.0}),
        accepted.vector[0],
    )
    equality_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_equalities=(equality,),
    )
    equality_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=equality_problem.identity,
        occurrence_identity="fixed-stencil-coordinate-occurrence",
    )
    evidence = derive_uncertainty_evidence(
        equality_accepted,
        problem=equality_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(controlled_id),
        resolved_environment_identity="fixed-affine-stencil",
    )
    assert evidence.residual_jacobian is None
    assert evidence.covariance is None
    assert evidence.failures[0].category == "exhausted_exact_feasible_distance"


def test_affine_held_only_and_equality_semantics_fail_closed() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    held_id, held_value = problem.held_items[0]
    held_only = AffineHalfSpace(
        "held-only",
        _full_frame_coefficients(problem, {held_id: 1.0}),
        held_value,
    )
    equality = AffineEquality(
        "controlled-equality",
        _full_frame_coefficients(problem, {problem.controlled_ids[0]: 1.0}),
        accepted.vector[0],
    )
    held_problem = dataclasses.replace(problem, affine_half_spaces=(held_only,))
    held_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=held_problem.identity,
        occurrence_identity="held-only-occurrence",
    )
    held_evidence = derive_uncertainty_evidence(
        held_accepted,
        problem=held_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="held-affine",
    )
    assert held_evidence.covariance is not None
    assert held_evidence.covariance.claim("AFFINE_FEASIBILITY") is ClaimState.SATISFIED
    assert held_evidence.covariance.usable

    held_equality = AffineEquality(
        "held-equality",
        _full_frame_coefficients(problem, {held_id: 1.0}),
        held_value,
    )
    held_equality_problem = dataclasses.replace(
        problem, affine_equalities=(held_equality,)
    )
    held_equality_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=held_equality_problem.identity,
        occurrence_identity="held-equality-occurrence",
    )
    held_equality_evidence = derive_uncertainty_evidence(
        held_equality_accepted,
        problem=held_equality_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="held-equality-affine",
    )
    assert held_equality_evidence.covariance is not None
    assert held_equality_evidence.covariance.usable
    assert (
        held_equality_evidence.covariance.claim("CONTROLLED_AFFINE_SEPARATION")
        is ClaimState.NOT_APPLICABLE
    )

    equality_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_equalities=(equality,),
    )
    equality_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=equality_problem.identity,
        occurrence_identity="equality-occurrence",
    )
    equality_evidence = derive_uncertainty_evidence(
        equality_accepted,
        problem=equality_problem,
        parameterization=parameterization,
        engine=engine,
        policy=_qualification_policy(problem.controlled_ids[0]),
        resolved_environment_identity="equality-affine",
    )
    assert equality_evidence.residual_jacobian is None
    assert equality_evidence.covariance is None
    assert equality_evidence.failures[0].category == "exhausted_exact_feasible_distance"


def test_affine_coupled_negative_degenerate_and_nonfinite_cases() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    held_id, held_value = problem.held_items[0]
    policy = _qualification_policy(controlled_id)
    baseline = derive_uncertainty_evidence(
        accepted,
        problem=problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="affine-exception-baseline",
    )
    assert baseline.covariance is not None
    standard_error = math.sqrt(baseline.covariance.covariance[0][0])

    coupled = AffineHalfSpace(
        "coupled",
        _full_frame_coefficients(problem, {controlled_id: 1.0, held_id: 1.0}),
        accepted.vector[0] + held_value + 4.0 * standard_error,
    )
    coupled_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(coupled,),
    )
    coupled_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=coupled_problem.identity,
        occurrence_identity="coupled-affine-occurrence",
    )
    coupled_evidence = derive_uncertainty_evidence(
        coupled_accepted,
        problem=coupled_problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="affine-coupled",
    )
    assert coupled_evidence.covariance is not None, coupled_evidence.failures
    assert (
        coupled_evidence.covariance.claim("AFFINE_FEASIBILITY") is ClaimState.SATISFIED
    )

    invalid_restrictions = (
        AffineHalfSpace(
            "negative",
            _full_frame_coefficients(problem, {controlled_id: 1.0}),
            accepted.vector[0] - 1.0,
        ),
        AffineHalfSpace(
            "nonfinite",
            _full_frame_coefficients(
                problem,
                {
                    controlled_id: float(np.finfo(np.float64).max),
                    held_id: float(np.finfo(np.float64).max),
                },
            ),
            float(np.finfo(np.float64).max),
        ),
    )
    for restriction in invalid_restrictions:
        with pytest.raises(DirectTrfConstructionError, match="affine half-space"):
            dataclasses.replace(
                problem,
                start=accepted.vector,
                affine_half_spaces=(restriction,),
            )

    positive = AffineHalfSpace(
        "positive-with-zero-variance",
        _full_frame_coefficients(problem, {controlled_id: 1.0}),
        accepted.vector[0] + 1.0,
    )
    degenerate_problem = dataclasses.replace(
        problem,
        start=accepted.vector,
        affine_half_spaces=(positive,),
    )
    degenerate_accepted = _qualified_accepted_copy(
        accepted,
        problem_identity=degenerate_problem.identity,
        occurrence_identity="degenerate-affine-occurrence",
        chi_square=0.0,
    )
    degenerate = derive_uncertainty_evidence(
        degenerate_accepted,
        problem=degenerate_problem,
        parameterization=parameterization,
        engine=engine,
        policy=policy,
        resolved_environment_identity="affine-degenerate",
    )
    assert degenerate.covariance is not None
    assert degenerate.covariance.claim("AFFINE_FEASIBILITY") is ClaimState.INDETERMINATE


def test_analytic_capability_identity_binds_actual_callable_semantics() -> None:
    first = FunctionAnalyticPartialCapability(
        "pop_2st",
        "pa",
        "same-user-facing-label",
        (_analytic_pa_kab, _analytic_pa_kba),
    )
    different = FunctionAnalyticPartialCapability(
        "pop_2st",
        "pa",
        "same-user-facing-label",
        (_analytic_pa_kab_alternative, _analytic_pa_kba),
    )
    rebound = dataclasses.replace(
        first,
        partials=(_analytic_pa_kab_alternative, _analytic_pa_kba),
    )
    assert first.identity != different.identity
    assert first.identity != rebound.identity
    assert first.implementation_fingerprints != different.implementation_fingerprints

    parameterization, _engine, _problem, _accepted = _accepted_hd_fit()
    output_id = "__PA_1_0_1000"

    def wrong_domain(_single: float) -> float:
        return 1.0

    incompatible = FunctionAnalyticPartialCapability(
        "pop_2st",
        "pa",
        "incompatible-domain",
        (wrong_domain, _analytic_pa_kba),
    )
    with pytest.raises(ValueError, match="domain arity"):
        compile_constraint_linearization_capabilities(
            parameterization,
            (output_id,),
            (incompatible,),
        )


def test_clustered_singular_evidence_retains_only_invariant_projectors() -> None:
    identity = np.eye(3)
    angle = 0.731
    rotated = np.asarray(
        (
            (math.cos(angle), -math.sin(angle), 0.0),
            (math.sin(angle), math.cos(angle), 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    exact = uncertainty_module._invariant_singular_subspaces(
        identity,
        (2.0, 2.0, 1.0),
        rank_threshold=1.0e-12,
        weak_threshold=1.0e-6,
        cluster_relative_tolerance=1.0e-10,
    )
    rotated_exact = uncertainty_module._invariant_singular_subspaces(
        rotated,
        (2.0, 2.0, 1.0),
        rank_threshold=1.0e-12,
        weak_threshold=1.0e-6,
        cluster_relative_tolerance=1.0e-10,
    )
    assert exact[0].classification == "clustered_identifiable"
    assert exact[1].classification == "isolated_identifiable"
    np.testing.assert_allclose(
        exact[0].projector, rotated_exact[0].projector, atol=2e-16
    )

    near = uncertainty_module._invariant_singular_subspaces(
        identity,
        (2.0, 2.0 * (1.0 - 5.0e-11), 1.0),
        rank_threshold=1.0e-12,
        weak_threshold=1.0e-6,
        cluster_relative_tolerance=1.0e-10,
    )
    assert near[0].classification == "clustered_identifiable"
    classified = uncertainty_module._invariant_singular_subspaces(
        identity,
        (2.0, 1.0e-7, 0.0),
        rank_threshold=1.0e-9,
        weak_threshold=1.0e-6,
        cluster_relative_tolerance=1.0e-10,
    )
    assert tuple(item.classification for item in classified) == (
        "isolated_identifiable",
        "isolated_weak",
        "isolated_null",
    )


def test_cancellation_after_svd_prevents_covariance_assembly() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    observed_svd = False
    original_svd = uncertainty_module.svd

    def traced_svd(*args: object, **kwargs: object) -> object:
        nonlocal observed_svd
        result = original_svd(*args, **kwargs)
        observed_svd = True
        return result

    with (
        patch("chemex.optimize.uncertainty.svd", side_effect=traced_svd),
        patch(
            "chemex.optimize.uncertainty._canonical_covariance_reduction",
            wraps=uncertainty_module._canonical_covariance_reduction,
        ) as covariance_reduction,
    ):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if observed_svd else None
            ),
            resolved_environment_identity="cancel-after-svd",
        )
    covariance_reduction.assert_not_called()
    assert evidence.rank_diagnostic is None
    assert evidence.covariance is None
    assert evidence.operations[-1].terminal is OperationTerminal.CANCELLED


def test_cancellation_after_rank_retains_diagnostic_but_not_covariance() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    rank_subspaces_complete = False
    original_subspaces = uncertainty_module._invariant_singular_subspaces

    def traced_subspaces(*args: object, **kwargs: object) -> object:
        nonlocal rank_subspaces_complete
        result = original_subspaces(*args, **kwargs)
        rank_subspaces_complete = True
        return result

    with (
        patch(
            "chemex.optimize.uncertainty._invariant_singular_subspaces",
            side_effect=traced_subspaces,
        ),
        patch(
            "chemex.optimize.uncertainty._canonical_covariance_reduction"
        ) as covariance_reduction,
    ):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(problem.controlled_ids[0]),
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if rank_subspaces_complete else None
            ),
            resolved_environment_identity="cancel-after-rank",
        )
    covariance_reduction.assert_not_called()
    assert evidence.rank_diagnostic is not None
    assert evidence.covariance is None
    assert evidence.operations[-1].terminal is OperationTerminal.CANCELLED


def test_cancellation_before_marginals_and_constraint_phases_is_ordered() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    derived_id = "__R1A_B_G2N_H_800_0MHZ"
    original_reduction = uncertainty_module._canonical_covariance_reduction
    covariance_complete = False

    def traced_reduction(*args: object, **kwargs: object) -> object:
        nonlocal covariance_complete
        result = original_reduction(*args, **kwargs)
        covariance_complete = True
        return result

    with (
        patch(
            "chemex.optimize.uncertainty._canonical_covariance_reduction",
            side_effect=traced_reduction,
        ),
        patch("chemex.optimize.uncertainty._marginal_errors") as marginal,
    ):
        before_marginal = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(controlled_id),
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if covariance_complete else None
            ),
            resolved_environment_identity="cancel-before-marginal",
        )
    marginal.assert_not_called()
    assert before_marginal.covariance is not None
    assert before_marginal.marginal_errors is None

    original_correlations = uncertainty_module._correlations
    correlations_complete = False

    def traced_correlations(*args: object, **kwargs: object) -> object:
        nonlocal correlations_complete
        result = original_correlations(*args, **kwargs)
        correlations_complete = True
        return result

    with (
        patch(
            "chemex.optimize.uncertainty._correlations",
            side_effect=traced_correlations,
        ),
        patch("chemex.optimize.uncertainty._linearize_constraints") as linearize_g,
    ):
        before_g = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(controlled_id),
            constrained_scope=(derived_id,),
            constrained_units=((derived_id, ParameterUnit.RATE_PER_SECOND),),
            constrained_scales=((derived_id, 1.0),),
            compiled_constraint_linearization=compile_constraint_linearization_capabilities(
                parameterization, (derived_id,), ()
            ),
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if correlations_complete else None
            ),
            resolved_environment_identity="cancel-before-g",
        )
    linearize_g.assert_not_called()
    assert before_g.constraint_jacobian is None


def test_cancellation_during_g_and_before_propagation_emits_no_downstream() -> None:
    _session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    derived_id = "__R1A_B_G2N_H_800_0MHZ"
    compiled = compile_constraint_linearization_capabilities(
        parameterization, (controlled_id, derived_id), ()
    )
    original_target = uncertainty_module._differentiate_target
    target_started = False

    def traced_target(*args: object, **kwargs: object) -> object:
        nonlocal target_started
        result = original_target(*args, **kwargs)
        target_started = True
        return result

    common = {
        "problem": problem,
        "parameterization": parameterization,
        "engine": engine,
        "policy": _qualification_policy(controlled_id),
        "constrained_scope": (controlled_id, derived_id),
        "constrained_units": (
            (controlled_id, ParameterUnit.RATE_PER_SECOND),
            (derived_id, ParameterUnit.RATE_PER_SECOND),
        ),
        "constrained_scales": ((controlled_id, 1.0), (derived_id, 1.0)),
        "compiled_constraint_linearization": compiled,
    }
    with (
        patch(
            "chemex.optimize.uncertainty._differentiate_target",
            side_effect=traced_target,
        ),
        patch("chemex.optimize.uncertainty._propagate_constraints") as propagate,
    ):
        during_g = derive_uncertainty_evidence(
            accepted,
            **common,
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if target_started else None
            ),
            resolved_environment_identity="cancel-during-g",
        )
    propagate.assert_not_called()
    assert during_g.constraint_jacobian is None
    assert during_g.constrained_propagation is None

    original_linearize = uncertainty_module._linearize_constraints
    g_complete = False

    def traced_linearize(*args: object, **kwargs: object) -> object:
        nonlocal g_complete
        result = original_linearize(*args, **kwargs)
        g_complete = True
        return result

    with (
        patch(
            "chemex.optimize.uncertainty._linearize_constraints",
            side_effect=traced_linearize,
        ),
        patch("chemex.optimize.uncertainty._propagate_constraints") as propagate,
    ):
        before_propagation = derive_uncertainty_evidence(
            accepted,
            **common,
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if g_complete else None
            ),
            resolved_environment_identity="cancel-before-propagation",
        )
    propagate.assert_not_called()
    assert before_propagation.constraint_jacobian is not None
    assert before_propagation.constrained_propagation is None


def test_cancellation_between_constrained_marginals_and_correlations() -> None:
    session, parameterization, engine, problem, accepted = _accepted_relaxation_fit()
    controlled_id = problem.controlled_ids[0]
    derived_id = "__R1A_B_G2N_H_800_0MHZ"
    before = session.analysis_values.snapshot()
    original_marginals = uncertainty_module._marginal_errors
    original_correlations = uncertainty_module._correlations
    constrained_marginals_complete = False
    constrained_correlations_started = False

    def traced_marginals(*args: object, **kwargs: object) -> object:
        nonlocal constrained_marginals_complete
        result = original_marginals(*args, **kwargs)
        if kwargs.get("source_family") == "constrained_propagation":
            constrained_marginals_complete = True
        return result

    def traced_correlations(*args: object, **kwargs: object) -> object:
        nonlocal constrained_correlations_started
        if kwargs.get("source_family") == "constrained_propagation":
            constrained_correlations_started = True
        return original_correlations(*args, **kwargs)

    with (
        patch(
            "chemex.optimize.uncertainty._marginal_errors",
            side_effect=traced_marginals,
        ),
        patch(
            "chemex.optimize.uncertainty._correlations",
            side_effect=traced_correlations,
        ),
    ):
        evidence = derive_uncertainty_evidence(
            accepted,
            problem=problem,
            parameterization=parameterization,
            engine=engine,
            policy=_qualification_policy(controlled_id),
            constrained_scope=(derived_id,),
            constrained_units=((derived_id, ParameterUnit.RATE_PER_SECOND),),
            constrained_scales=((derived_id, 1.0),),
            compiled_constraint_linearization=compile_constraint_linearization_capabilities(
                parameterization, (derived_id,), ()
            ),
            cancellation_probe=lambda: (
                OperationTerminal.CANCELLED if constrained_marginals_complete else None
            ),
            resolved_environment_identity="cancel-before-constrained-correlations",
        )

    assert not constrained_correlations_started
    assert evidence.constrained_propagation is not None
    assert evidence.constrained_marginal_errors is not None
    assert evidence.constrained_correlations is None
    assert evidence.operations[-1].stage == "constrained_correlations"
    assert evidence.operations[-1].terminal is OperationTerminal.CANCELLED
    assert all(
        operation.stage != "final_bundle_assembly" for operation in evidence.operations
    )
    assert session.analysis_values.snapshot() == before
