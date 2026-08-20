"""Qualification tests for the #585 native evaluation boundary.

The public seam is an immutable evaluation plan bound to an isolated evaluator.
Legacy Profile residuals remain the independent scientific oracle at this
checkpoint.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from threading import Thread
from unittest.mock import patch

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
    _identity,
    _parameterization_failure_validity,
)
from chemex.experiments.builder import build_experiments
from chemex.parameters.parameterization import (
    AmbiguousParameterReferenceError,
    ConstraintCycleError,
    ConstraintDomainError,
    ConstraintEvaluationError,
    ConstraintProgramMismatchError,
    ConstraintSelfReferenceError,
    IncompatibleParameterizationInputError,
    IncompleteParameterDependenciesError,
    ModelDerivationOverrideError,
    NonFiniteParameterValueError,
    NoParameterMatchError,
    UnsupportedConstraintExpressionError,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
PARAMETERS = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"
CEST_EXPERIMENT = ROOT / "examples/Experiments/CEST_13C_LABEL_CN/Experiments/23hz.toml"
CEST_PARAMETERS = (
    ROOT / "examples/Experiments/CEST_13C_LABEL_CN/Parameters/parameters.toml"
)


def _shipped_dcest(*names: str) -> tuple[AnalysisSession, Experiments]:
    session = AnalysisSession.create()
    session.set_model("2st_hd")
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(
            include=[SpinSystem.from_name(name) for name in (names or ("1N",))],
            exclude=None,
        ),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    return session, experiments


def _evaluation_frame(
    session: AnalysisSession, parameterization: object
) -> EvaluationFrame:
    return EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )


def test_cest_infinite_sweep_width_sentinel_has_a_stable_native_plan() -> None:
    identities: list[str] = []
    for _ in range(2):
        session = AnalysisSession.create()
        session.set_model("2st")
        experiments = build_experiments(
            [CEST_EXPERIMENT],
            Selection(
                include=[SpinSystem.from_name("L18CB")],
                exclude=None,
            ),
            session=session,
        )
        session.parameters.set_defaults(read_defaults([CEST_PARAMETERS]))
        assert session.try_build_analysis_values()
        parameterization = session.compile_current_parameterization(
            experiments.param_ids
        )
        engine = EvaluationEngine.from_experiments(
            experiments,
            parameterization,
        )
        identities.append(engine.plan.identity)
        outcome = engine.new_evaluator().evaluate(
            _evaluation_frame(session, parameterization)
        )
        assert isinstance(outcome, EvaluationResult)

    assert identities[0] == identities[1]


def test_shipped_two_profile_dcest_plan_matches_legacy_completely() -> None:
    native_session, native_experiments = _shipped_dcest("1N", "2N")
    legacy_session, legacy_experiments = _shipped_dcest("1N", "2N")
    parameterization = native_session.compile_parameterization(
        Method(), native_experiments.param_ids
    )
    native = EvaluationEngine.from_experiments(native_experiments, parameterization)
    evaluator = native.new_evaluator()
    outcome = evaluator.evaluate(_evaluation_frame(native_session, parameterization))
    assert isinstance(outcome, EvaluationResult)
    legacy_parameters = legacy_session.parameters.build_lmfit_params(
        legacy_experiments.param_ids
    )
    legacy_residuals = legacy_experiments.residuals(legacy_parameters)
    legacy_profiles = next(iter(legacy_experiments)).profiles
    assert [item.profile_ordinal for item in native.plan.profiles] == [0, 1]
    assert outcome.resolved_values == pytest.approx(legacy_parameters.valuesdict())
    # The isolated native workspace has bitwise-identical pulse-entry inputs,
    # but the SciPy/BLAS propagator path is allowed a few binary64 ulps across
    # independent workspace executions on Linux.  Keep this separate from the
    # deliberate ordered-normalization tolerance below.
    np.testing.assert_allclose(
        outcome.unscaled_calculations,
        np.concatenate([profile.data.calc_unscaled for profile in legacy_profiles]),
        rtol=6.0e-16,
        atol=0.0,
    )
    # #572 deliberately replaces legacy BLAS dot accumulation with the frozen
    # ordered binary64 reduction.  The DCEST difference is confined to scale.
    np.testing.assert_allclose(
        outcome.normalized_calculations,
        np.concatenate([profile.data.calc for profile in legacy_profiles]),
        rtol=6.0e-16,
        atol=1.2e-8,
    )
    np.testing.assert_allclose(outcome.residuals, legacy_residuals, rtol=6.0e-13)
    np.testing.assert_allclose(
        [item.normalization_factor for item in outcome.profiles],
        [profile.data.scale for profile in legacy_profiles],
        rtol=6.0e-16,
    )
    assert [item.retained_observation_indices for item in outcome.profiles] == [
        tuple(np.flatnonzero(profile.data.mask)) for profile in legacy_profiles
    ]
    repeated = evaluator.evaluate(_evaluation_frame(native_session, parameterization))
    assert isinstance(repeated, EvaluationResult)
    np.testing.assert_array_equal(repeated.residuals, outcome.residuals)
    assert evaluator.cache_statistics.hits == 2


@pytest.mark.parametrize(
    "field,value",
    (
        ("kernel_identity", "other:kernel:v2"),
        ("kernel_configuration", "{}"),
        ("spectrometer_configuration", "{}"),
        ("local_inputs", (("wrong", "__wrong"),)),
        ("observation_metadata", "[]"),
        ("output_shape", (999,)),
        ("observation_offset", 1),
        ("experiment_ordinal", 1),
        ("profile_ordinal", 1),
        ("identity", "foreign-profile-identity"),
    ),
)
def test_trusted_rebinding_rejects_each_mutated_kernel_descriptor(
    field: str, value: object
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    payload = json.loads(json.dumps(engine.plan.to_record()))
    serialized_value = (
        [list(item) for item in value]
        if field == "local_inputs"
        else list(value)
        if field == "output_shape"
        else value
    )
    payload["profiles"][0][field] = serialized_value
    descriptor = dataclasses.replace(engine.plan.profiles[0], **{field: value})
    signed_plan = EvaluationPlan(
        engine.plan.parameterization_identity,
        engine.plan.constraint_program_identity,
        (descriptor,),
        resolved_ids=engine.plan.resolved_ids,
    )
    payload["identity"] = signed_plan.identity
    if field in {
        "output_shape",
        "observation_offset",
        "experiment_ordinal",
        "profile_ordinal",
        "identity",
    }:
        with pytest.raises(ValueError):
            EvaluationPlan.from_record(payload)
    else:
        plan = EvaluationPlan.from_record(payload)
        with pytest.raises(ValueError, match="Trusted"):
            EvaluationEngine.bind(plan, parameterization, experiments)


def test_trusted_rebinding_rejects_nested_serialized_basis_descriptor_mutation() -> (
    None
):
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    payload = json.loads(json.dumps(engine.plan.to_record()))
    spectrometer = json.loads(payload["profiles"][0]["spectrometer_configuration"])
    spectrometer["basis_type"] = "other-basis"
    configuration = json.dumps(spectrometer, separators=(",", ":"), sort_keys=True)
    payload["profiles"][0]["spectrometer_configuration"] = configuration
    descriptor = dataclasses.replace(
        engine.plan.profiles[0], spectrometer_configuration=configuration
    )
    payload["identity"] = EvaluationPlan(
        engine.plan.parameterization_identity,
        engine.plan.constraint_program_identity,
        (descriptor,),
        resolved_ids=engine.plan.resolved_ids,
    ).identity
    plan = EvaluationPlan.from_record(payload)
    with pytest.raises(ValueError, match="Trusted"):
        EvaluationEngine.bind(plan, parameterization, experiments)


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("normalization_version", "other-normalization"),
        ("residual_version", "other-residual"),
    ),
)
def test_serialized_plan_execution_semantics_must_match_local_contract(
    field: str, value: str
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    payload = json.loads(json.dumps(engine.plan.to_record()))
    payload[field] = value
    signed_plan = dataclasses.replace(engine.plan, **{field: value})
    payload["identity"] = signed_plan.identity
    plan = EvaluationPlan.from_record(payload)
    with pytest.raises(ValueError, match="execution semantics"):
        EvaluationEngine.bind(plan, parameterization, experiments)


def test_spectrometer_native_kernel_descriptor_is_immutable_by_value() -> None:
    _session, experiments = _shipped_dcest()
    descriptor = (
        next(iter(experiments)).profiles[0].spectrometer.native_kernel_descriptor()
    )
    with pytest.raises(TypeError):
        descriptor["basis_type"] = "other-basis"
    with pytest.raises(TypeError):
        descriptor["jeff_i"]["dephasing"] = True


def test_native_snapshot_failure_does_not_disable_legacy_dcest() -> None:
    normal_session, normal_experiments = _shipped_dcest()
    normal_parameters = normal_session.parameters.build_lmfit_params(
        normal_experiments.param_ids
    )
    normal_residuals = normal_experiments.residuals(normal_parameters)
    normal_profile = next(iter(normal_experiments)).profiles[0]
    normal_calculation = normal_profile.data.calc.copy()

    with patch(
        "chemex.nmr.spectrometer.deepcopy",
        side_effect=RuntimeError("native snapshot failed"),
    ):
        failed_session, failed_experiments = _shipped_dcest()
    failed_profile = next(iter(failed_experiments)).profiles[0]
    failed_spectrometer = failed_profile.spectrometer
    assert failed_spectrometer.native_construction_diagnostic == (
        "RuntimeError: native snapshot failed"
    )
    with pytest.raises(ValueError, match="Native construction is unavailable"):
        failed_spectrometer.native_kernel_descriptor()
    with pytest.raises(ValueError, match="Native construction is unavailable"):
        failed_spectrometer.new_native_workspace()

    failed_parameters = failed_session.parameters.build_lmfit_params(
        failed_experiments.param_ids
    )
    failed_residuals = failed_experiments.residuals(failed_parameters)
    np.testing.assert_array_equal(failed_profile.data.calc, normal_calculation)
    np.testing.assert_array_equal(failed_residuals, normal_residuals)
    parameterization = failed_session.compile_parameterization(
        Method(), failed_experiments.param_ids
    )
    with pytest.raises(ValueError, match="Native construction is unavailable"):
        EvaluationEngine.from_experiments(failed_experiments, parameterization)

    failed_profile.calculate(failed_parameters)
    assert not failed_spectrometer.try_finalize_native_construction()
    with pytest.raises(ValueError, match="native snapshot failed"):
        failed_spectrometer.native_kernel_descriptor()


def test_lifecycle_frames_reject_stale_and_foreign_snapshots_before_projection() -> (
    None
):
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    snapshot = session.analysis_values.snapshot()
    with pytest.raises(ValueError):
        parameterization.frame_from_snapshot(dataclasses.replace(snapshot, revision=1))
    foreign, _ = _shipped_dcest()
    with pytest.raises(ValueError):
        parameterization.frame_from_snapshot(foreign.analysis_values.snapshot())


def test_evaluation_frame_validates_canonical_binary64_input_and_projection_order() -> (
    None
):
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    lifecycle = parameterization.frame_from_snapshot(session.analysis_values.snapshot())

    with pytest.raises(ValueError, match="unique"):
        EvaluationFrame(parameterization.identity, (("a", 1.0), ("a", 2.0)))
    with pytest.raises(TypeError, match="real binary64"):
        EvaluationFrame(parameterization.identity, (("a", True),))
    with pytest.raises(ValueError, match="finite"):
        EvaluationFrame(parameterization.identity, (("a", np.inf),))
    with pytest.raises(TypeError, match="real binary64"):
        EvaluationFrame(parameterization.identity, (("a", 10**10_000),))
    with pytest.raises(ValueError, match="canonical independent ID order"):
        EvaluationFrame.from_lifecycle_frame(
            parameterization,
            dataclasses.replace(lifecycle, _items=tuple(reversed(lifecycle._items))),
        )
    with pytest.raises(ValueError, match="incompatible"):
        EvaluationFrame.from_lifecycle_frame(
            parameterization,
            dataclasses.replace(lifecycle, program_fingerprint="foreign"),
        )

    frame = EvaluationFrame.from_lifecycle_frame(parameterization, lifecycle)
    assert tuple(field.name for field in dataclasses.fields(frame)) == (
        "parameterization_identity",
        "_items",
    )


def test_later_profile_invalid_trial_is_atomic_resets_workspace_and_isolates_legacy() -> (
    None
):
    session, experiments = _shipped_dcest("1N", "2N")
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    evaluator = engine.new_evaluator()
    frame = _evaluation_frame(session, parameterization)
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate
    calls: list[object] = []

    def fail_second(self: object, spectrometer: object, data: object) -> Array:
        calls.append(self)
        if len(calls) == 2:
            return np.full_like(data.metadata, np.nan)
        return original(self, spectrometer, data)

    with patch.object(pulse_type, "calculate", fail_second):
        failed = evaluator.evaluate(frame)
    assert isinstance(failed, EvaluationFailure)
    assert failed.validity == "INVALID_TRIAL"
    assert len(calls) == 2
    assert evaluator.cache_statistics.entries == 0
    recovered = evaluator.evaluate(frame)
    fresh = engine.new_evaluator().evaluate(frame)
    assert isinstance(recovered, EvaluationResult)
    assert isinstance(fresh, EvaluationResult)
    np.testing.assert_array_equal(recovered.residuals, fresh.residuals)
    legacy = experiments.residuals(
        session.parameters.build_lmfit_params(experiments.param_ids)
    )
    np.testing.assert_allclose(legacy, recovered.residuals, rtol=6.0e-13)


def test_later_profile_implementation_failure_poisoning_requires_fresh_evaluator() -> (
    None
):
    session, experiments = _shipped_dcest("1N", "2N")
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = _evaluation_frame(session, parameterization)
    evaluator = engine.new_evaluator()
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate
    calls: list[object] = []

    def fail_second(self: object, spectrometer: object, data: object) -> Array:
        calls.append(self)
        if len(calls) == 2:
            return np.array([0.0])
        return original(self, spectrometer, data)

    with patch.object(pulse_type, "calculate", fail_second):
        failed = evaluator.evaluate(frame)
    assert isinstance(failed, EvaluationFailure)
    assert failed.validity == "IMPLEMENTATION_FAILURE"
    assert failed.category == "unexpected_shape"
    assert len(calls) == 2
    assert evaluator.cache_statistics.entries == 0
    poisoned = evaluator.evaluate(frame)
    assert isinstance(poisoned, EvaluationFailure)
    assert poisoned.category == "workspace_poisoned"
    recovered = engine.new_evaluator().evaluate(frame)
    assert isinstance(recovered, EvaluationResult)
    reference = engine.new_evaluator().evaluate(frame)
    assert isinstance(reference, EvaluationResult)
    assert recovered.identity == reference.identity


def test_one_evaluation_frame_has_identical_complete_results_across_evaluators() -> (
    None
):
    session, experiments = _shipped_dcest("1N", "2N")
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = _evaluation_frame(session, parameterization)
    first_engine = EvaluationEngine.from_experiments(experiments, parameterization)
    second_engine = EvaluationEngine.bind(
        EvaluationPlan.from_record(first_engine.plan.to_record()),
        parameterization,
        experiments,
    )
    first = first_engine.new_evaluator().evaluate(frame)
    second = second_engine.new_evaluator().evaluate(frame)
    assert isinstance(first, EvaluationResult)
    assert isinstance(second, EvaluationResult)
    assert first.identity == second.identity
    assert first.resolved_values == second.resolved_values
    assert first.profiles == second.profiles
    np.testing.assert_array_equal(
        first.unscaled_calculations, second.unscaled_calculations
    )
    np.testing.assert_array_equal(
        first.normalized_calculations, second.normalized_calculations
    )
    np.testing.assert_array_equal(first.residuals, second.residuals)


def test_native_plan_matches_legacy_complete_profile_evaluation() -> None:
    native_session, native_experiments = _shipped_dcest()
    legacy_session, legacy_experiments = _shipped_dcest()
    parameterization = native_session.compile_parameterization(
        Method(), native_experiments.param_ids
    )
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(native_session.analysis_values.snapshot()),
    )
    engine = EvaluationEngine.from_experiments(native_experiments, parameterization)

    assert not hasattr(frame, "occurrence_identity")
    assert not hasattr(frame, "revision")

    outcome = engine.new_evaluator().evaluate(frame)
    assert isinstance(outcome, EvaluationResult)

    legacy_parameters = legacy_session.parameters.build_lmfit_params(
        legacy_experiments.param_ids
    )
    legacy_residuals = legacy_experiments.residuals(legacy_parameters)
    legacy_profile = next(iter(legacy_experiments)).profiles[0]

    assert outcome.plan_identity == engine.plan.identity
    assert outcome.resolved_values == pytest.approx(legacy_parameters.valuesdict())
    # Native and legacy calculate through separate isolated Spectrometers.
    # Their pulse-entry scalar values, matrices, weights, and detection vector
    # are bitwise-identical; Linux SciPy/BLAS propagator reductions may differ
    # by a few ulps between those independent executions.
    np.testing.assert_allclose(
        outcome.unscaled_calculations,
        legacy_profile.data.calc_unscaled,
        rtol=6.0e-16,
        atol=0.0,
    )
    np.testing.assert_allclose(
        outcome.normalized_calculations,
        legacy_profile.data.calc,
        rtol=6.0e-16,
        atol=1.2e-8,
    )
    np.testing.assert_allclose(outcome.residuals, legacy_residuals, rtol=6.0e-13)
    assert outcome.profiles[0].normalization_factor == pytest.approx(
        legacy_profile.data.scale, rel=6.0e-16
    )
    assert outcome.profiles[0].retained_observation_indices == tuple(
        np.flatnonzero(legacy_profile.data.mask)
    )


def test_native_private_workspace_has_bitwise_identical_dcest_pulse_inputs() -> None:
    native_session, native_experiments = _shipped_dcest()
    legacy_session, legacy_experiments = _shipped_dcest()
    parameterization = native_session.compile_parameterization(
        Method(), native_experiments.param_ids
    )
    frame = _evaluation_frame(native_session, parameterization)
    captured: list[dict[str, object]] = []
    pulse_type = type(next(iter(native_experiments)).profiles[0].pulse_sequence)
    original = pulse_type.calculate

    def capture(self: object, spectrometer: object, data: object) -> Array:
        engine = spectrometer._engine
        captured.append(
            {
                "settings": self.settings.model_dump(),
                "metadata_dtype": data.metadata.dtype,
                "metadata": data.metadata.copy(),
                "scalars": {
                    name: getattr(spectrometer, name)
                    for name in (
                        "ppm_i",
                        "ppm_s",
                        "carrier_i",
                        "carrier_s",
                        "offset_i",
                        "offset_s",
                        "b1_i",
                        "b1_s",
                        "detection",
                    )
                }
                | {"h_larmor_frq": engine.h_frq},
                "b1_i_values": spectrometer.b1_i_distribution.values.copy(),
                "b1_i_weights": spectrometer.b1_i_distribution.weights.copy(),
                "jeff_i_values": spectrometer.jeff_i.values.copy(),
                "jeff_i_weights": spectrometer.jeff_i.weights.copy(),
                "parameter_values": dict(engine.par_values),
                "matrix_keys": frozenset(engine._matrices),
                "matrices": {
                    name: matrix.copy() for name, matrix in engine._matrices.items()
                },
                "l_free": engine.l_free.copy(),
                "weights": engine.weights.copy(),
                "detection_vector": engine._readout._detect_vector.copy(),
            }
        )
        return original(self, spectrometer, data)

    with patch.object(pulse_type, "calculate", capture):
        outcome = (
            EvaluationEngine.from_experiments(native_experiments, parameterization)
            .new_evaluator()
            .evaluate(frame)
        )
        legacy_experiments.residuals(
            legacy_session.parameters.build_lmfit_params(legacy_experiments.param_ids)
        )
    assert isinstance(outcome, EvaluationResult)
    native, legacy = captured
    assert native["settings"] == legacy["settings"]
    assert native["metadata_dtype"] == legacy["metadata_dtype"]
    assert native["scalars"] == legacy["scalars"]
    assert native["parameter_values"] == legacy["parameter_values"]
    assert native["matrix_keys"] == legacy["matrix_keys"]
    for name in (
        "metadata",
        "b1_i_values",
        "b1_i_weights",
        "jeff_i_values",
        "jeff_i_weights",
        "l_free",
        "weights",
        "detection_vector",
    ):
        np.testing.assert_array_equal(native[name], legacy[name])
    for matrix_name in native["matrix_keys"]:
        np.testing.assert_array_equal(
            native["matrices"][matrix_name], legacy["matrices"][matrix_name]
        )


def test_plan_identity_is_deterministic_and_rebinds_after_serialization() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)

    record = engine.plan.to_record()
    restored = EvaluationPlan.from_record(record)
    rebound = EvaluationEngine.bind(restored, parameterization, experiments)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )

    assert restored.identity == engine.plan.identity
    assert restored.profiles == engine.plan.profiles
    original = engine.new_evaluator().evaluate(frame)
    rebound_result = rebound.new_evaluator().evaluate(frame)
    assert isinstance(original, EvaluationResult)
    assert isinstance(rebound_result, EvaluationResult)
    np.testing.assert_array_equal(
        original.unscaled_calculations, rebound_result.unscaled_calculations
    )
    np.testing.assert_array_equal(
        original.normalized_calculations, rebound_result.normalized_calculations
    )
    np.testing.assert_array_equal(original.residuals, rebound_result.residuals)
    restored_result = EvaluationResult.from_record(
        json.loads(json.dumps(original.to_record())), restored
    )
    np.testing.assert_array_equal(
        original.unscaled_calculations, restored_result.unscaled_calculations
    )
    np.testing.assert_array_equal(
        original.normalized_calculations, restored_result.normalized_calculations
    )
    np.testing.assert_array_equal(original.residuals, restored_result.residuals)


def test_workspace_reuse_cache_and_fresh_workspace_are_scientifically_identical() -> (
    None
):
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    reused = engine.new_evaluator()

    first = reused.evaluate(frame)
    second = reused.evaluate(frame)
    fresh = engine.new_evaluator().evaluate(frame)
    changed_id = parameterization.independent_ids[0]
    changed_frame = EvaluationFrame(
        frame.parameterization_identity,
        tuple(
            (
                param_id,
                session.analysis_values.snapshot()[param_id] * 1.000001
                if param_id == changed_id
                else value,
            )
            for param_id, value in frame._items
        ),
    )
    changed = reused.evaluate(changed_frame)

    assert isinstance(first, EvaluationResult)
    assert isinstance(second, EvaluationResult)
    assert isinstance(fresh, EvaluationResult)
    assert isinstance(changed, EvaluationResult)
    assert reused.cache_statistics.hits == 1
    assert reused.cache_statistics.misses == 2
    np.testing.assert_array_equal(
        first.unscaled_calculations, second.unscaled_calculations
    )
    np.testing.assert_array_equal(
        first.normalized_calculations, second.normalized_calculations
    )
    np.testing.assert_array_equal(first.residuals, fresh.residuals)
    assert not first.unscaled_calculations.flags.writeable
    assert not first.normalized_calculations.flags.writeable
    assert not first.residuals.flags.writeable
    with pytest.raises(ValueError):
        first.residuals[0] = 0.0


def test_native_normalization_uses_the_no_epsilon_572_formula() -> None:
    session, experiments = _shipped_dcest()
    profile = next(iter(experiments)).profiles[0]
    profile.data.err = np.linspace(1.0e-18, 2.0e-15, profile.data.exp.size)
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    evaluator = EvaluationEngine.from_experiments(
        experiments, parameterization
    ).new_evaluator()
    outcome = evaluator.evaluate(_evaluation_frame(session, parameterization))
    assert isinstance(outcome, EvaluationResult)

    unscaled = outcome.unscaled_calculations
    expected = float(
        np.dot(profile.data.exp / profile.data.err, unscaled / profile.data.err)
        / np.dot(unscaled / profile.data.err, unscaled / profile.data.err)
    )
    assert outcome.profiles[0].normalization_factor == pytest.approx(expected)
    profile.calculate(session.parameters.build_lmfit_params(experiments.param_ids))
    assert profile.data.scale != expected


def test_distinct_evaluators_are_reentrant_while_workspace_owner_is_enforced() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    first = engine.new_evaluator()
    thread_outcomes: list[EvaluationResult | EvaluationFailure] = []

    def evaluate_second() -> None:
        thread_outcomes.append(engine.new_evaluator().evaluate(frame))

    worker = Thread(target=evaluate_second)
    worker.start()
    worker.join()
    owner_violation: list[EvaluationResult | EvaluationFailure] = []

    def misuse_first() -> None:
        owner_violation.append(first.evaluate(frame))

    misuse = Thread(target=misuse_first)
    misuse.start()
    misuse.join()

    assert isinstance(thread_outcomes[0], EvaluationResult)
    assert isinstance(owner_violation[0], EvaluationFailure)
    assert owner_violation[0].validity == "INVALID_REQUEST"


@pytest.mark.parametrize(
    ("values", "category", "validity"),
    (
        (np.array([np.nan]), "unexpected_shape", "IMPLEMENTATION_FAILURE"),
        (np.array([]), "invalid_normalization", "INVALID_TRIAL"),
    ),
)
def test_native_kernel_and_normalization_failures_are_atomic(
    values: Array,
    category: str,
    validity: str,
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    pulse_type = type(next(iter(experiments)).profiles[0].pulse_sequence)
    if not values.size:
        values = np.zeros_like(next(iter(experiments)).profiles[0].data.exp)

    def calculate(_self: object, _spectrometer: object, _data: object) -> Array:
        return values

    with patch.object(pulse_type, "calculate", calculate):
        outcome = (
            EvaluationEngine.from_experiments(experiments, parameterization)
            .new_evaluator()
            .evaluate(frame)
        )

    assert isinstance(outcome, EvaluationFailure)
    assert outcome.category == category
    assert outcome.validity == validity
    assert not hasattr(outcome, "residuals")


def test_non_finite_kernel_output_is_an_invalid_trial_without_partial_result() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )
    profile = next(iter(experiments)).profiles[0]
    pulse_type = type(profile.pulse_sequence)

    def calculate(_self: object, _spectrometer: object, _data: object) -> Array:
        return np.full_like(profile.data.exp, np.nan)

    with patch.object(pulse_type, "calculate", calculate):
        outcome = (
            EvaluationEngine.from_experiments(experiments, parameterization)
            .new_evaluator()
            .evaluate(frame)
        )

    assert isinstance(outcome, EvaluationFailure)
    assert outcome.category == "non_finite_calculation"
    assert outcome.validity == "INVALID_TRIAL"
    restored_plan = EvaluationPlan.from_record(
        EvaluationEngine.from_experiments(
            experiments, parameterization
        ).plan.to_record()
    )
    assert EvaluationFailure.from_record(outcome.to_record(), restored_plan) == outcome


def test_native_failure_is_qualification_only_and_cannot_veto_legacy() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    with (
        patch(
            "chemex.evaluation.native.EvaluationEngine.from_experiments",
            side_effect=ValueError("native qualification failure"),
        ),
        pytest.raises(ValueError, match="native qualification failure"),
    ):
        EvaluationEngine.from_experiments(experiments, parameterization)
    legacy = experiments.residuals(
        session.parameters.build_lmfit_params(experiments.param_ids)
    )
    assert legacy.size and np.all(np.isfinite(legacy))


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("experimental_values", 1.0),
        ("uncertainties", 1.0),
        ("mask", False),
        ("is_scaled", False),
        ("source_identity", "foreign-source"),
    ),
)
def test_trusted_rebinding_rejects_complete_profile_population_tampering(
    field: str, value: object
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    descriptor = engine.plan.profiles[0]
    if field in {"experimental_values", "uncertainties", "mask"}:
        changed = list(getattr(descriptor, field))
        changed[0] = not changed[0] if field == "mask" else value
        value = tuple(changed)
    tampered = dataclasses.replace(descriptor, **{field: value})
    # A hostile serializer can recompute every fingerprint; trusted rebinding
    # must still compare the complete local scientific descriptor.
    if field in {"source_identity"}:
        tampered = dataclasses.replace(
            tampered,
            identity=_identity(
                "profile-plan",
                (
                    tampered.source_identity,
                    tampered.experiment_ordinal,
                    tampered.profile_ordinal,
                    tampered.observation_offset,
                ),
            ),
        )
    plan = EvaluationPlan.from_record(
        json.loads(
            json.dumps(
                dataclasses.replace(engine.plan, profiles=(tampered,)).to_record()
            )
        )
    )
    with pytest.raises(ValueError):
        EvaluationEngine.bind(plan, parameterization, experiments)


def test_frame_result_and_failure_records_fail_closed_after_tampering() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = _evaluation_frame(session, parameterization)
    result = engine.new_evaluator().evaluate(frame)
    assert isinstance(result, EvaluationResult)
    for mutate in (
        lambda record: record["items"].__setitem__(0, ["bad", "0x1.0000000000000p+0"]),
        lambda record: record.__setitem__("identity", "bad"),
    ):
        record = json.loads(json.dumps(frame.to_record()))
        mutate(record)
        with pytest.raises((TypeError, ValueError)):
            EvaluationFrame.from_record(record)
    for mutate in (
        lambda record: record["profiles"][0].__setitem__("observation_offset", 999),
        lambda record: record["profiles"][0].__setitem__(
            "normalization_factor", "0x1.0000000000000p+0"
        ),
        lambda record: record["residuals"].__setitem__(0, "0x1.0000000000000p+0"),
        lambda record: record.__setitem__("parameterization_identity", "foreign"),
    ):
        record = json.loads(json.dumps(result.to_record()))
        mutate(record)
        with pytest.raises((TypeError, ValueError)):
            EvaluationResult.from_record(record, engine.plan)
    failure = EvaluationFailure(
        engine.plan.identity,
        engine.plan.parameterization_identity,
        "kernel",
        "non_finite_calculation",
        "INVALID_TRIAL",
        engine.plan.profiles[0].identity,
    )
    record = failure.to_record()
    record["identity"] = "foreign"
    with pytest.raises(ValueError):
        EvaluationFailure.from_record(record, engine.plan)


def test_result_buffers_are_owned_immutable_bytes() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    result = (
        EvaluationEngine.from_experiments(experiments, parameterization)
        .new_evaluator()
        .evaluate(_evaluation_frame(session, parameterization))
    )
    assert isinstance(result, EvaluationResult)
    for values in (
        result.unscaled_calculations,
        result.normalized_calculations,
        result.residuals,
    ):
        before = values.copy()
        with pytest.raises(ValueError):
            values[0] = 0.0
        with pytest.raises(ValueError):
            values.flags.writeable = True
        assert isinstance(values.base, bytes)
        np.testing.assert_array_equal(values, before)


def test_native_evaluation_is_stable_under_hostile_numpy_error_state() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    previous = np.seterr(all="raise")
    try:
        outcome = (
            EvaluationEngine.from_experiments(experiments, parameterization)
            .new_evaluator()
            .evaluate(_evaluation_frame(session, parameterization))
        )
    finally:
        np.seterr(**previous)
    assert isinstance(outcome, EvaluationResult)


def test_semantic_parameterization_identity_excludes_lifecycle_identity() -> None:
    first_session, first_experiments = _shipped_dcest()
    second_session, second_experiments = _shipped_dcest()
    first = first_session.compile_parameterization(
        Method(), first_experiments.param_ids
    )
    second = second_session.compile_parameterization(
        Method(), second_experiments.param_ids
    )
    assert first.identity != second.identity
    assert first.evaluator_identity == second.evaluator_identity
    assert (
        EvaluationEngine.from_experiments(
            first_experiments, first
        ).plan.parameterization_identity
        == EvaluationEngine.from_experiments(
            second_experiments, second
        ).plan.parameterization_identity
    )


def test_masked_and_empty_unscaled_profiles_preserve_complete_calculations() -> None:
    session, experiments = _shipped_dcest()
    profile = next(iter(experiments)).profiles[0]
    profile.data.mask[0] = False
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    masked_engine = EvaluationEngine.from_experiments(experiments, parameterization)
    restored_plan = EvaluationPlan.from_record(
        json.loads(json.dumps(masked_engine.plan.to_record()))
    )
    masked = (
        EvaluationEngine.bind(restored_plan, parameterization, experiments)
        .new_evaluator()
        .evaluate(_evaluation_frame(session, parameterization))
    )
    assert isinstance(masked, EvaluationResult)
    assert masked.unscaled_calculations.size == profile.data.exp.size
    assert masked.normalized_calculations.size == profile.data.exp.size
    assert masked.residuals.size == profile.data.exp.size - 1
    assert masked.profiles[0].retained_observation_indices == tuple(
        np.flatnonzero(profile.data.mask)
    )
    assert masked.profiles[0].observation_offset == 0
    assert masked.profiles[0].residual_offset == 0
    assert masked.profiles[0].residual_count == profile.data.exp.size - 1
    retained = profile.data.mask
    unscaled = masked.unscaled_calculations
    numerator = 0.0
    denominator = 0.0
    for exp, calc, err in zip(
        profile.data.exp[retained],
        unscaled[retained],
        profile.data.err[retained],
        strict=True,
    ):
        numerator += (calc / err) * (exp / err)
        denominator += (calc / err) * (calc / err)
    assert masked.profiles[0].normalization_factor == numerator / denominator
    profile.is_scaled = False
    profile.data.mask[:] = False
    empty = (
        EvaluationEngine.from_experiments(experiments, parameterization)
        .new_evaluator()
        .evaluate(_evaluation_frame(session, parameterization))
    )
    assert isinstance(empty, EvaluationResult)
    assert empty.profiles[0].normalization_factor == 1.0
    assert empty.residuals.size == 0
    assert empty.unscaled_calculations.size == profile.data.exp.size


@pytest.mark.parametrize("sign", (0.0, -1.0))
def test_finite_zero_and_negative_normalization_are_valid(sign: float) -> None:
    session, experiments = _shipped_dcest()
    profile = next(iter(experiments)).profiles[0]
    profile.data.exp = sign * np.abs(profile.data.exp)
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    result = (
        EvaluationEngine.from_experiments(experiments, parameterization)
        .new_evaluator()
        .evaluate(_evaluation_frame(session, parameterization))
    )
    assert isinstance(result, EvaluationResult)
    factor = result.profiles[0].normalization_factor
    assert factor == 0.0 if sign == 0.0 else factor < 0.0


def test_descriptor_is_history_independent_and_evaluator_rejects_reentry_and_pid() -> (
    None
):
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    first = EvaluationEngine.from_experiments(experiments, parameterization).plan
    profile = next(iter(experiments)).profiles[0]
    profile.calculate(session.parameters.build_lmfit_params(experiments.param_ids))
    second = EvaluationEngine.from_experiments(experiments, parameterization).plan
    assert first.identity == second.identity
    evaluator = EvaluationEngine.from_experiments(
        experiments, parameterization
    ).new_evaluator()
    frame = _evaluation_frame(session, parameterization)
    pulse_type = type(profile.pulse_sequence)
    original = pulse_type.calculate
    nested: list[EvaluationFailure] = []

    def calculate(self: object, spectrometer: object, data: object) -> Array:
        nested_result = evaluator.evaluate(frame)
        assert isinstance(nested_result, EvaluationFailure)
        nested.append(nested_result)
        return original(self, spectrometer, data)

    with patch.object(pulse_type, "calculate", calculate):
        assert isinstance(evaluator.evaluate(frame), EvaluationResult)
    assert nested[0].category == "workspace_reentrant"
    evaluator._owner_pid = -1
    failed = evaluator.evaluate(frame)
    assert isinstance(failed, EvaluationFailure)
    assert failed.category == "workspace_process_violation"


def test_dcest_construction_descriptor_precedes_and_survives_pulse_history() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    profile = next(iter(experiments)).profiles[0]
    initial = dict(profile.spectrometer.native_kernel_descriptor())
    profile.calculate(session.parameters.build_lmfit_params(experiments.param_ids))
    after_legacy = dict(profile.spectrometer.native_kernel_descriptor())
    before_plan = EvaluationEngine.from_experiments(experiments, parameterization).plan
    profile.calculate(session.parameters.build_lmfit_params(experiments.param_ids))
    after_plan = EvaluationEngine.from_experiments(experiments, parameterization).plan
    assert initial == after_legacy
    assert before_plan.identity == after_plan.identity


def test_construction_setting_changes_native_descriptor_and_plan_identity() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    original = EvaluationEngine.from_experiments(experiments, parameterization).plan
    changed = next(iter(experiments)).profiles[0].spectrometer
    changed._native_kernel_descriptor = None
    changed._native_workspace_template = None
    changed._native_construction_attempted = False
    changed._native_construction_diagnostic = None
    changed.detection = "[iz_b]"
    changed.finalize_native_construction()
    changed_plan = EvaluationEngine.from_experiments(experiments, parameterization).plan
    assert changed.native_kernel_descriptor()["detection"] == "[iz_b]"
    assert original.identity != changed_plan.identity


def test_native_workspace_cannot_mutate_authoritative_profile_or_metadata() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    profile = next(iter(experiments)).profiles[0]
    metadata = profile.data.metadata.copy()
    profile.data.mask[0] = False
    legacy_params = session.parameters.build_lmfit_params(experiments.param_ids)
    profile.calculate(legacy_params)
    legacy_unscaled = profile.data.calc_unscaled.copy()
    offset = profile.spectrometer.offset_i
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    baseline = engine.new_evaluator().evaluate(
        _evaluation_frame(session, parameterization)
    )
    assert isinstance(baseline, EvaluationResult)
    observed_parameter_ids: list[set[str]] = []
    observed_data: list[tuple[Array, Array, Array]] = []
    pulse_type = type(profile.pulse_sequence)
    original = pulse_type.calculate
    spectrometer_type = type(profile.spectrometer)
    original_update = spectrometer_type.update

    def update(self: object, values: dict[str, float]) -> None:
        observed_parameter_ids.append(set(values))
        original_update(self, values)

    def mutate(self: object, spectrometer: object, data: object) -> Array:
        observed_data.append((data.exp.copy(), data.err.copy(), data.mask.copy()))
        data.metadata[:] = 999.0
        spectrometer.offset_i = 999.0
        return original(self, spectrometer, data)

    with (
        patch.object(pulse_type, "calculate", mutate),
        patch.object(spectrometer_type, "update", update),
    ):
        result = engine.new_evaluator().evaluate(
            _evaluation_frame(session, parameterization)
        )
    assert isinstance(result, EvaluationResult)
    np.testing.assert_array_equal(profile.data.metadata, metadata)
    assert profile.spectrometer.offset_i == offset
    assert observed_parameter_ids
    assert observed_parameter_ids[0] == set(profile.name_map)
    assert all(
        np.array_equal(exp, np.zeros_like(exp)) for exp, _err, _mask in observed_data
    )
    assert all(
        np.array_equal(err, np.ones_like(err)) for _exp, err, _mask in observed_data
    )
    assert all(np.all(mask) for _exp, _err, mask in observed_data)
    repeated = engine.new_evaluator().evaluate(
        _evaluation_frame(session, parameterization)
    )
    assert isinstance(repeated, EvaluationResult)
    np.testing.assert_array_equal(
        baseline.unscaled_calculations, repeated.unscaled_calculations
    )
    profile.calculate(legacy_params)
    np.testing.assert_array_equal(profile.data.calc_unscaled, legacy_unscaled)


def test_normalization_overflow_is_typed_under_hostile_numpy_settings() -> None:
    session, experiments = _shipped_dcest()
    profile = next(iter(experiments)).profiles[0]
    profile.data.err[:] = np.finfo(np.float64).tiny
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = _evaluation_frame(session, parameterization)
    default = (
        EvaluationEngine.from_experiments(experiments, parameterization)
        .new_evaluator()
        .evaluate(frame)
    )
    previous = np.seterr(all="raise")
    try:
        outcome = (
            EvaluationEngine.from_experiments(experiments, parameterization)
            .new_evaluator()
            .evaluate(frame)
        )
        assert np.geterr() == {
            "divide": "raise",
            "over": "raise",
            "under": "raise",
            "invalid": "raise",
        }
    finally:
        np.seterr(**previous)
    assert isinstance(default, EvaluationFailure)
    assert (default.stage, default.category, default.validity) == (
        outcome.stage,
        outcome.category,
        outcome.validity,
    )
    assert isinstance(outcome, EvaluationFailure)
    assert (outcome.stage, outcome.category, outcome.validity) == (
        "normalization",
        "invalid_normalization",
        "INVALID_TRIAL",
    )


def test_normalization_is_independent_of_blas_thread_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    frame = _evaluation_frame(session, parameterization)
    outcomes: list[EvaluationResult] = []
    for threads in ("1", "4"):
        monkeypatch.setenv("OPENBLAS_NUM_THREADS", threads)
        outcome = (
            EvaluationEngine.from_experiments(experiments, parameterization)
            .new_evaluator()
            .evaluate(frame)
        )
        assert isinstance(outcome, EvaluationResult)
        outcomes.append(outcome)
    np.testing.assert_array_equal(
        outcomes[0].normalized_calculations, outcomes[1].normalized_calculations
    )
    assert [item.normalization_factor for item in outcomes[0].profiles] == [
        item.normalization_factor for item in outcomes[1].profiles
    ]


@pytest.mark.parametrize(
    ("error_type", "expected"),
    (
        (NoParameterMatchError, "INVALID_PLAN_OR_BINDING"),
        (AmbiguousParameterReferenceError, "INVALID_PLAN_OR_BINDING"),
        (ConstraintSelfReferenceError, "INVALID_PLAN_OR_BINDING"),
        (ConstraintCycleError, "INVALID_PLAN_OR_BINDING"),
        (ConstraintDomainError, "INVALID_TRIAL"),
        (ConstraintEvaluationError, "INVALID_TRIAL"),
        (NonFiniteParameterValueError, "INVALID_TRIAL"),
        (IncompleteParameterDependenciesError, "INVALID_PLAN_OR_BINDING"),
        (IncompatibleParameterizationInputError, "INVALID_PLAN_OR_BINDING"),
        (ConstraintProgramMismatchError, "INVALID_PLAN_OR_BINDING"),
        (UnsupportedConstraintExpressionError, "INVALID_PLAN_OR_BINDING"),
        (ModelDerivationOverrideError, "INVALID_PLAN_OR_BINDING"),
    ),
)
def test_parameterization_failure_taxonomy_is_closed(
    error_type: type[Exception], expected: str
) -> None:
    error = error_type("qualification")
    assert _parameterization_failure_validity(error) == expected


@pytest.mark.parametrize(
    ("error_type", "expected"),
    (
        (ConstraintDomainError, "INVALID_TRIAL"),
        (NonFiniteParameterValueError, "INVALID_TRIAL"),
        (ConstraintProgramMismatchError, "INVALID_PLAN_OR_BINDING"),
    ),
)
def test_evaluator_maps_declared_parameterization_failures(
    error_type: type[Exception], expected: str
) -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    evaluator = EvaluationEngine.from_experiments(
        experiments, parameterization
    ).new_evaluator()
    with patch.object(
        type(parameterization), "resolve", side_effect=error_type("test")
    ):
        outcome = evaluator.evaluate(_evaluation_frame(session, parameterization))
    assert isinstance(outcome, EvaluationFailure)
    assert outcome.validity == expected


def test_projection_and_result_assembly_failures_are_typed() -> None:
    session, experiments = _shipped_dcest()
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    profile = next(iter(experiments)).profiles[0]
    profile.name_map["missing"] = "__missing"
    try:
        with pytest.raises(ValueError, match="outside its closure"):
            EvaluationEngine.from_experiments(experiments, parameterization)
    finally:
        del profile.name_map["missing"]
    evaluator = engine.new_evaluator()
    with patch("chemex.evaluation.native.EvaluationResult", side_effect=RuntimeError):
        outcome = evaluator.evaluate(_evaluation_frame(session, parameterization))
    assert isinstance(outcome, EvaluationFailure)
    assert (outcome.stage, outcome.validity) == ("result", "IMPLEMENTATION_FAILURE")
