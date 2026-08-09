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
)
from chemex.experiments.builder import build_experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).parent.parent
EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
PARAMETERS = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"


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
    np.testing.assert_array_equal(
        outcome.unscaled_calculations,
        np.concatenate([profile.data.calc_unscaled for profile in legacy_profiles]),
    )
    np.testing.assert_array_equal(
        outcome.normalized_calculations,
        np.concatenate([profile.data.calc for profile in legacy_profiles]),
    )
    np.testing.assert_array_equal(outcome.residuals, legacy_residuals)
    assert [item.normalization_factor for item in outcome.profiles] == [
        profile.data.scale for profile in legacy_profiles
    ]
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
    )
    payload["identity"] = signed_plan.identity
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
    with pytest.raises(ValueError, match="finite"):
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
    np.testing.assert_array_equal(legacy, recovered.residuals)


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
    np.testing.assert_allclose(
        outcome.unscaled_calculations,
        legacy_profile.data.calc_unscaled,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_array_equal(
        outcome.normalized_calculations, legacy_profile.data.calc
    )
    np.testing.assert_array_equal(outcome.residuals, legacy_residuals)
    assert outcome.profiles[0].normalization_factor == legacy_profile.data.scale
    assert outcome.profiles[0].retained_observation_indices == tuple(
        np.flatnonzero(legacy_profile.data.mask)
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
    assert outcome.profiles[0].normalization_factor == expected
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
