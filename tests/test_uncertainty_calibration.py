"""Verification for the executable prospective issue #601 calibration."""

# ruff: noqa: E701, E702, I001 -- compact verification for the frozen artifact.
from __future__ import annotations
import hashlib
import importlib
import json
import math
from dataclasses import replace
from types import SimpleNamespace
from typing import Any, cast
import pytest
from tests.qualification import capture_uncertainty_calibration as calibration


# fmt: off
def selected_policy(rank: float = 2.0**-40) -> Any:
    calibration._import_machinery()
    return calibration.compose_uncertainty_policy("A1", (1.0, 100.0), (2.0**-24, 256.0, 8, 8), "gesvd", (0.0, rank), max(2.0**-22, rank), 2.0**-33, 2.0**19, 64.0)
def empty_counts(role: str = "canonical") -> dict[str, int]:
    return dict.fromkeys(calibration.CEILINGS[role], 0)
def canonical_record(policy: Any) -> dict[str, object]:
    serialized = calibration.policy_record(policy)
    phases = {"scales": {"state_population": 1.0, "exchange_rate": 100.0}, "finite_difference": (2.0**-24, 256.0, 8, 8), "svd_driver": "gesvd", "rank": (0.0, 2.0**-40), "weak": 2.0**-22, "cluster": 2.0**-33, "conditioning": 2.0**19, "correlation": 64.0}
    metrics = {"scale": {"status": "SELECTED"}, "finite_difference": {"status": "SELECTED", "policy": phases["finite_difference"]}, "driver": {"status": "SELECTED", "driver": "gesvd"}, "rank": {"status": "SELECTED", "policy": phases["rank"]}, **{name: {"status": "SELECTED", "value": phases[name]} for name in ("weak", "cluster", "conditioning", "correlation")}}
    return {"source_digest": hashlib.sha256(calibration.Path(calibration.__file__).read_bytes()).hexdigest(), "specification_id": calibration.SPECIFICATION_ID, "specification_digest": calibration.frozen_digest(), "status": "QUALIFIED", "canonical_provenance": {}, "selected_phase_policies": phases, "decisive_metrics": metrics, "policy": serialized, "policy_digest": hashlib.sha256(json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode()).hexdigest(), "composed": ({"passed": True},), "holdouts": ({"passed": True},) * 6}
def test_frozen_wayfinder_populations_and_digest() -> None:
    assert calibration.COUNTS == (17, 10_296, 2, 400, 27, 27, 27, 14)
    assert len(tuple(calibration.finite_difference_policies())) == 10_296
    assert calibration.EXTENTS == (0, 2, 4, 8, 12, 16) and (0.0, 1.0, 4, 12) in calibration.finite_difference_policies() and (0.0, 1.0, 12, 4) in calibration.finite_difference_policies()
    assert len(calibration.RANK_GRID) ** 2 == 400
    assert calibration.frozen_digest() == "9c34b947a024a81cd53785192fc129f90b616ca99d9f5302744c7f2b63f59a3e"
def test_catalogue_partition_is_executable_complete_and_disjoint() -> None:
    actual = calibration.actual_catalogue_names(); supported, unsupported = calibration.catalogue_partition()
    assert actual == supported | unsupported and not supported & unsupported and "l2_free" in unsupported
    assert len(actual) == 277 and {item[1] for item in calibration.SUPPORTED} == {"dimensionless", "fraction", "rate", "frequency"}
    with pytest.raises(KeyError, match="not catalogued"):
        calibration.classify_scale_name("invented_parameter")
def test_all_frozen_truths_are_executable() -> None:
    calibration._import_machinery()
    names = ("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "C4", "C5", "F2", "H1", "H2", "H3", "H4", "H5", "H6")
    assert all(calibration._calculation(name, calibration._setup(name)[0], 1.0) is not None and calibration.truth_jacobian(name) for name in names)
    assert calibration.scientific_partials(0.5, 2.0) == pytest.approx((math.exp(0.25) / 2.0, -0.5 * math.exp(0.25) / 4.0 + 0.8))
def test_a0_order_threshold_ambiguity_and_rank_neighbors(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(calibration, "SCALE_EXPONENTS", (-1, 0, 1)); monkeypatch.setattr(calibration, "finite_difference_policies", lambda: iter(((0.0, 1.0, 0, 0),)))
    monkeypatch.setattr(calibration, "_score", lambda _t, q, _o, _p, _c: (0, 9.0, 99.0, 2 if q[0] == 0.5 else 1))
    value, record = calibration._select_one_scale("A0", 1.0, ((None,),), (1.0,), empty_counts())
    assert value == 1.0 and record["exponent"] == 0
    assert calibration._threshold_outcome((1.0, 2.0, 3.0), [1.0, 3.0], 3.0) == "AMBIGUOUS_SELECTION" and calibration._threshold_outcome((1.0, 2.0), [2.0], 2.0) == "GRID_SATURATED"
    observations = tuple(SimpleNamespace(singular_values=spectrum) for spectrum in ((1.0e-9, 2.0e-21), (1.0, 5.0e-13), (math.sqrt(28.0), 0.0)))
    rank, rank_record = calibration.select_rank_policy(observations, empty_counts())
    neighbors = cast("tuple[dict[str, object], ...]", rank_record["rejected_neighbors"])
    assert rank is not None and len(neighbors) == 4 and all({"status", "decisive_reasons", "error", "cost"} <= set(item) for item in neighbors)
def test_b_consumes_selected_a_and_c_consumes_selected_ab() -> None:
    phase_a = selected_policy(0.0); counts = empty_counts()
    b, b_record = calibration._spectral_observation("B1", phase_a, "gesvd", (0.0, 0.0), "environment", counts)
    phase_ab = replace(phase_a, rank_relative_tolerance=2.0**-40, weak_relative_tolerance=2.0**-40)
    assert b is not None and b.source_policy.relative_step_tolerance == phase_a.relative_step_tolerance and b_record["policy_identity"] == b.source_policy.identity
    for name in ("C1", "C2", "C3", "C4"):
        diagnostic, record = calibration._spectral_observation(name, phase_ab, "gesvd", (0.0, 2.0**-40), "environment", counts)
        assert diagnostic is not None and diagnostic.source_policy.rank_relative_tolerance == phase_ab.rank_relative_tolerance and record["policy_identity"] == diagnostic.source_policy.identity
    c5, c5_record = calibration._derive_case("C5", phase_ab, "environment", counts, "canonical")
    assert c5_record["passed"] and c5.constrained_propagation is not None and c5.constrained_correlations is not None
def test_both_scalings_and_h3_h5_use_real_typed_evidence() -> None:
    policy = selected_policy(2.0**-20)
    for scaling in (calibration.ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES, calibration.ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE):
        evidence, record = calibration._derive_case("A3", policy, "environment", empty_counts(), "canonical", scaling=scaling)
        assert record["passed"] and record["m"] == 4 and record["n"] == 2 and record["g"] == 1 and record["nu"] == 1 and record["covariance_truth"] is None
        assert evidence.rank_diagnostic is not None and evidence.source_policy.residual_variance_scaling is scaling
    h3, h3_record = calibration._derive_case("H3", policy, "environment", empty_counts(), "canonical", scaling=calibration.ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE)
    h5, _record = calibration._derive_case("H5", selected_policy(), "environment", empty_counts(), "canonical", numerical_partial=True)
    assert h3_record["passed"] and h3.residual_jacobian is not None and h3.rank_diagnostic is not None
    assert h5.constraint_jacobian is not None and h5.constrained_propagation is not None and h5.constrained_marginal_errors is not None
    zero = h5.constraint_jacobian.output_ids.index("ZERO")
    assert h5.constraint_jacobian.matrix[zero] == (0.0, 0.0) and all(value == 0.0 for value in h5.constrained_propagation.factor[zero])
    assert all(h5.constrained_propagation.covariance[zero][i] == h5.constrained_propagation.covariance[i][zero] == 0.0 for i in range(len(h5.constrained_propagation.output_ids)))
    assert h5.constrained_propagation.claim("LOCAL_FIRST_ORDER_DEGENERACY") is calibration.ClaimState.VIOLATED and h5.constrained_propagation.claim("OUTPUT_FIRST_ORDER_NONDEGENERACY::ZERO") is calibration.ClaimState.VIOLATED and not h5.constrained_marginal_errors.scope_reportable
def test_every_budget_refuses_before_the_first_over_budget_operation(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    class Evaluator:
        def evaluate(self, _frame: object) -> None: calls.append("evaluation")
    counts = dict(calibration.CEILINGS["canonical"]); owner = SimpleNamespace(charging=True, counts=counts, role="canonical"); wrapped = calibration._CountingEvaluator(Evaluator(), owner)
    with pytest.raises(RuntimeError): wrapped.evaluate(object())
    monkeypatch.setattr(calibration.uncertainty_module, "svd", lambda *_a, **_k: calls.append("svd"))
    with pytest.raises(RuntimeError), calibration._charged_svd(counts, "canonical"): calibration.uncertainty_module.svd(object())
    monkeypatch.setattr(calibration, "_correlations", lambda **_k: calls.append("correlation"))
    with pytest.raises(RuntimeError): calibration._correlation_evidence(object(), object(), counts, "canonical")
    calibration._ACTIVE_COUNTS = counts; calibration._ACTIVE_ROLE = "canonical"
    with pytest.raises(RuntimeError): calibration.counted_scientific_value(0.5, 2.0)
    assert calls == []
def test_holdout_failure_is_irreversible(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(calibration, "_derive_case", lambda name, *_a, **_k: (calls.append(name), {"case": name, "passed": False}))
    passed, records = calibration.run_holdouts(object(), "environment", empty_counts())
    assert not passed and calls == ["H1"] and tuple(item["case"] for item in records) == ("H1",)
def test_attestation_precedes_import(monkeypatch: pytest.MonkeyPatch) -> None:
    reloaded = importlib.reload(calibration); calls: list[str] = []
    monkeypatch.setattr(reloaded, "attest_lane", lambda *_: calls.append("attest") or {}); monkeypatch.setattr(reloaded, "_import_machinery", lambda: calls.append("import"))
    monkeypatch.setattr(reloaded, "validate_lane_records", lambda *_: (_ for _ in ()).throw(RuntimeError("stop")))
    with pytest.raises(RuntimeError, match="stop"): reloaded.acquire("canonical", "digest")
    assert calls == ["attest", "import"]
def test_compatibility_is_bound_to_canonical_selected_policy(monkeypatch: pytest.MonkeyPatch) -> None:
    policy = selected_policy(); record = canonical_record(policy); monkeypatch.setattr(calibration, "validate_lane_records", lambda *_: "environment")
    rebuilt = calibration.policy_from_canonical_record(record)
    evidence, result = calibration._derive_case("A2", rebuilt, "compatibility-environment", empty_counts("compatibility"), "compatibility")
    assert isinstance(rebuilt, calibration.UncertaintyPolicy) and result["passed"] and evidence.residual_jacobian is not None
    runner_up = replace(policy, correlation_roundoff_multiplier=128.0); record["policy"] = calibration.policy_record(runner_up); record["policy_digest"] = hashlib.sha256(json.dumps(record["policy"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    with pytest.raises(RuntimeError, match="not the canonical selected policy"): calibration.policy_from_canonical_record(record)
