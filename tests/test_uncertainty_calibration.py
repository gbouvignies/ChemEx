"""Verification for the executable prospective issue #601 calibration."""

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


def selected_policy(rank: float = 2.0**-40) -> Any:
    calibration._import_machinery()
    return calibration.compose_uncertainty_policy(
        "A1",
        (1.0, 100.0),
        (2.0**-24, 256.0, 8, 8),
        "gesvd",
        (0.0, rank),
        max(2.0**-22, rank),
        2.0**-32,
        2.0**20,
        64.0,
    )


def empty_counts(role: str = "canonical") -> dict[str, int]:
    return dict.fromkeys(calibration.CEILINGS[role], 0)


def canonical_record(policy: Any) -> dict[str, object]:
    scales = {
        family: q0 for family, _stratum, q0, _unit, _patterns in calibration.SUPPORTED
    } | {"A4-argument-0": 1.0, "A4-argument-1": 2.0, "A4-output": 4.0}
    policy = replace(
        policy,
        calibration_identity=f"{calibration.SPECIFICATION_ID}:{calibration._scale_catalogue_digest(scales)}",
    )
    serialized = calibration.policy_record(policy)
    phases = {
        "scales": scales,
        "finite_difference": (2.0**-24, 256.0, 8, 8),
        "svd_driver": "gesvd",
        "rank": (0.0, 2.0**-40),
        "weak": 2.0**-22,
        "cluster": 2.0**-32,
        "conditioning": 2.0**20,
        "correlation": 64.0,
    }
    families = {
        family: {
            "status": "SELECTED",
            "family": family,
            "q0": q0,
            "exponent": 0,
            "scale": q0,
            "candidate_identity": calibration._scale_candidate_identity(
                family, q0, 0, q0
            ),
        }
        for family, _stratum, q0, _unit, _patterns in calibration.SUPPORTED
    }
    metrics = {
        "scale": {"status": "SELECTED", "families": families},
        "finite_difference": {
            "status": "SELECTED",
            "policy": phases["finite_difference"],
        },
        "driver": {
            "status": "SELECTED",
            "driver": "gesvd",
            "cases": tuple(
                {"case": name, "status": "QUALIFIED"} for name in ("B1", "B2", "B3")
            ),
        },
        "rank": {
            "status": "SELECTED",
            "policy": phases["rank"],
            "typed_truth_cases": tuple(
                {"case": name, "passed": True} for name in ("B2", "B3")
            ),
        },
        "phase_c_cases": tuple(
            (
                {"case": name, "status": "ACQUIRED"}
                if name != "C5"
                else {"case": name, "passed": True}
            )
            for name in ("C1", "C2", "C3", "C4", "C5")
        ),
        **{
            name: {"status": "SELECTED", "value": phases[name]}
            for name in ("weak", "cluster", "conditioning", "correlation")
        },
    }
    return {
        "source_digest": hashlib.sha256(
            calibration.Path(calibration.__file__).read_bytes()
        ).hexdigest(),
        "specification_id": calibration.SPECIFICATION_ID,
        "specification_digest": calibration.frozen_digest(),
        "status": "QUALIFIED",
        "canonical_provenance": {},
        "selected_phase_policies": phases,
        "decisive_metrics": metrics,
        "policy": serialized,
        "policy_digest": hashlib.sha256(
            json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "composed": tuple(
            {"case": name, "passed": True} for name in calibration.COMPOSED_CASES
        ),
        "holdouts": tuple(
            {"case": name, "passed": True} for name in calibration.HOLDOUT_CASES
        ),
    }


def test_frozen_wayfinder_populations_and_digest() -> None:
    assert (
        calibration.COUNTS == (17, 10_296, 2, 400, 27, 27, 27, 14)
        and len(tuple(calibration.finite_difference_policies())) == 10_296
        and calibration.PLANNED_BY_ROLE["canonical"]["svd"] == 30
        and calibration.PLANNED_BY_ROLE["canonical"]["correlation"] == 22
        and calibration.PLANNED_BY_ROLE["compatibility"]["svd"] == 4
    )
    assert (
        calibration.EXTENTS == (0, 2, 4, 8, 12, 16)
        and (0.0, 1.0, 4, 12) in calibration.finite_difference_policies()
        and (0.0, 1.0, 12, 4) in calibration.finite_difference_policies()
    )
    assert len(calibration.RANK_GRID) ** 2 == 400
    source = calibration.Path(calibration.__file__).read_bytes()
    digest = calibration.frozen_digest()
    changed_source = source.replace(
        b"1.0 + 128.0 * 2.0**-52", b"1.0 + 256.0 * 2.0**-52", 1
    )
    assert (
        digest == calibration._specification_digest(source)
        and changed_source != source
        and calibration._specification_digest(changed_source) != digest
    )


def test_canonical_v2_fixture_records_corrected_fail_closed_result() -> None:
    fixture_path = (
        calibration.Path(__file__).parent
        / "fixtures"
        / "canonical_uncertainty_calibration_v2.json"
    )
    record = json.loads(fixture_path.read_text())
    correlation = record["decisive_metrics"]["correlation"]
    assert (
        record["source_digest"]
        == hashlib.sha256(
            calibration.Path(calibration.__file__).read_bytes()
        ).hexdigest()
        and record["specification_digest"] == calibration.frozen_digest()
        and record["status"] == "COMPOSED_VALIDATION_FAILED"
        and record["selected_phase_policies"]["correlation"] == 64.0
        and correlation["status"] == "SELECTED"
        and correlation["passed"] is True
        and correlation["truth_raw"]
        == [
            0.75,
            1.0 + 64.0 * 2.0**-52,
            -(1.0 + 64.0 * 2.0**-52),
            1.0 + 128.0 * 2.0**-52,
        ]
        and correlation["outcomes"]
        == [
            "AVAILABLE",
            "ENDPOINT_CANONICALIZED",
            "ENDPOINT_CANONICALIZED",
            "CORRELATION_RANGE_VIOLATION",
        ]
        and record["composed"][-1]["case"] == "F1-absolute"
        and record["composed"][-1]["passed"] is False
        and all(item["status"] == "NOT_RUN" for item in record["holdouts"])
        and record["compatibility"]
        == {"reason": "composed_validation_failed", "status": "NOT_RUN"}
        and "raw_truth_mismatch" not in fixture_path.read_text()
    )


def test_catalogue_partition_is_executable_complete_and_disjoint() -> None:
    actual = calibration.actual_catalogue_names()
    supported, unsupported = calibration.catalogue_partition()
    assert (
        actual == supported | unsupported
        and not supported & unsupported
        and "l2_free" in unsupported
        and len(actual) == 277
        and {item[1] for item in calibration.SUPPORTED}
        == {"dimensionless", "fraction", "rate", "frequency"}
    )
    with pytest.raises(KeyError, match="not catalogued"):
        calibration.classify_scale_name("invented_parameter")


def test_all_frozen_truths_are_executable() -> None:
    calibration._import_machinery()
    names = (
        "A1",
        "A2",
        "A3",
        "B1",
        "B2",
        "B3",
        "C1",
        "C2",
        "C3",
        "C4",
        "F2",
        "H1",
        "H2",
        "H3",
        "H4",
        "H5",
        "H6",
    )
    assert all(
        calibration._calculation(name, calibration._setup(name)[0], 1.0) is not None
        and calibration.truth_jacobian(name)
        for name in names
    ) and calibration.scientific_partials(0.5, 2.0) == pytest.approx(
        (math.exp(0.25) / 2.0, -0.5 * math.exp(0.25) / 4.0 + 0.8)
    )


def test_a0_order_threshold_ambiguity_and_rank_neighbors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(calibration, "SCALE_EXPONENTS", (-1, 0, 1))
    monkeypatch.setattr(
        calibration, "finite_difference_policies", lambda: iter(((0.0, 1.0, 0, 0),))
    )
    monkeypatch.setattr(
        calibration,
        "_score",
        lambda _t, q, _o, _p, _c: (0, 9.0, 99.0, 2 if q[0] == 0.5 else 1),
    )
    value, record = calibration._select_one_scale(
        "A0", 1.0, ((None,),), (1.0,), empty_counts()
    )
    assert value == 1.0 and record["exponent"] == 0
    assert (
        calibration._qualifier_failure(
            {
                (0,): (True, (0,), "policy-0"),
                (1,): (False, None, None),
                (2,): (True, (2,), "policy-2"),
            },
            ((0, 1, 2),),
        )
        == "NON_MONOTONE_ADEQUACY"
    )
    assert (
        calibration._qualifier_failure(
            {(0,): (True, (1,), "same-policy"), (1,): (True, (1,), "same-policy")},
            ((0, 1),),
        )
        == "AMBIGUOUS_SELECTION"
        and calibration._qualifier_failure(
            {(0,): (True, (1,), "policy-0"), (1,): (True, (1,), "policy-1")}, ((0, 1),)
        )
        is None
    )
    assert (
        calibration._threshold_outcome((1.0, 2.0, 3.0), [1.0, 3.0], 3.0)
        == "AMBIGUOUS_SELECTION"
        and calibration._threshold_outcome((1.0, 2.0), [2.0], 2.0) == "GRID_SATURATED"
    )
    observations = tuple(
        SimpleNamespace(singular_values=spectrum)
        for spectrum in ((1.0e-9, 2.0e-21), (1.0, 5.0e-13), (math.sqrt(28.0), 0.0))
    )
    rank, rank_record = calibration.select_rank_policy(observations, empty_counts())
    neighbors = cast("tuple[dict[str, object], ...]", rank_record["rejected_neighbors"])
    assert (
        rank == (0.0, 2.0**-40)
        and rank_record["policy"] == rank
        and len(neighbors) == 4
        and all(
            {"status", "decisive_reasons", "error", "cost"} <= set(item)
            for item in neighbors
        )
    )


def test_b_consumes_selected_a_and_c_consumes_selected_ab(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phase_a = selected_policy(0.0)
    counts = empty_counts()
    scales = {calibration.BC_SCALE_FAMILY: 1.0}
    b, b_record = calibration._spectral_observation(
        "B1", scales, phase_a, "gesvd", (0.0, 0.0), "environment", counts
    )
    scaled_b, scaled_record = calibration._spectral_observation(
        "B1",
        {calibration.BC_SCALE_FAMILY: 2.0},
        phase_a,
        "gesvd",
        (0.0, 0.0),
        "environment",
        counts,
    )
    phase_ab = replace(
        phase_a, rank_relative_tolerance=2.0**-40, weak_relative_tolerance=2.0**-40
    )
    final_weak_b2, _record = calibration._spectral_observation(
        "B2",
        scales,
        replace(phase_ab, weak_relative_tolerance=2.0**-22),
        "gesvd",
        (0.0, 2.0**-40),
        "environment",
        counts,
    )
    provisional_b2, _record = calibration._spectral_observation(
        "B2", scales, phase_ab, "gesvd", (0.0, 2.0**-40), "environment", counts
    )
    assert (
        b is not None
        and scaled_b is not None
        and b.source_policy.relative_step_tolerance == phase_a.relative_step_tolerance
        and b_record["policy_identity"] == b.source_policy.identity
    )
    assert scaled_record["policy_identity"] != b_record[
        "policy_identity"
    ] and scaled_b.singular_values[0] == pytest.approx(2.0 * b.singular_values[0])
    assert (
        provisional_b2 is not None
        and final_weak_b2 is not None
        and (
            provisional_b2.singular_values,
            provisional_b2.threshold,
            provisional_b2.rank,
            provisional_b2.identifiable_projector,
            provisional_b2.null_projector,
        )
        == (
            final_weak_b2.singular_values,
            final_weak_b2.threshold,
            final_weak_b2.rank,
            final_weak_b2.identifiable_projector,
            final_weak_b2.null_projector,
        )
    )
    for name in ("C1", "C2", "C3", "C4"):
        diagnostic, record = calibration._spectral_observation(
            name, scales, phase_ab, "gesvd", (0.0, 2.0**-40), "environment", counts
        )
        assert (
            diagnostic is not None
            and diagnostic.source_policy.rank_relative_tolerance
            == phase_ab.rank_relative_tolerance
            and record["policy_identity"] == diagnostic.source_policy.identity
        )
    passed, rank_cases = calibration.validate_rank_truth_cases(
        scales, phase_ab, "environment", counts
    )
    assert (
        passed
        and tuple(item["case"] for item in rank_cases) == ("B2", "B3")
        and all(not item["covariance_available"] for item in rank_cases)
    )
    perturbation = cast("dict[str, Any]", rank_cases[0]["perturbation"])
    assert (
        perturbation["relative_envelope"] == 1.0e-12
        and perturbation["unperturbed_error"] <= calibration.TOL["derivative"]
        and perturbation["weak_entries"] == pytest.approx((-5.0e-13, 1.5e-12))
        and perturbation["ranks"] == (1, 2)
        and perturbation["unstable_modes"] == (1,)
    )
    assert calibration._h4_spectral_case(
        selected_policy(), "environment", counts, "canonical"
    )["passed"]
    driver, qualified_record, _observations = calibration.select_svd_driver(
        scales, phase_a, "environment", counts
    )
    assert driver in calibration.SVD_DRIVERS and all(
        item["status"] == "QUALIFIED"
        for item in cast("tuple[dict[str, object], ...]", qualified_record["cases"])
    )
    monkeypatch.setattr(
        calibration,
        "_spectral_observation",
        lambda name, *_a: (
            SimpleNamespace(singular_values=(1.0, 0.0)),
            {"case": name, "status": "DISQUALIFIED"},
        ),
    )
    driver, record, _observations = calibration.select_svd_driver(
        scales, phase_a, "environment", empty_counts()
    )
    assert driver is None and record["status"] == "UNSUPPORTED"


def test_phase_c_replaces_provisional_weak_before_canonical_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provisional, final = 2.0**-40, 2.0**-22
    seen: dict[str, Any] = {}
    scales = {
        family: q0 for family, _stratum, q0, _unit, _patterns in calibration.SUPPORTED
    } | {"A4-argument-0": 1.0, "A4-argument-1": 2.0, "A4-output": 4.0}
    fd = (2.0**-24, 256.0, 8, 8)
    scale_record = {"status": "SELECTED", "families": tuple(scales)}
    driver_record = {
        "status": "SELECTED",
        "driver": "gesvd",
        "cases": tuple(
            {"case": name, "status": "QUALIFIED"} for name in ("B1", "B2", "B3")
        ),
    }
    monkeypatch.setattr(
        calibration, "select_scales", lambda _counts: (scales, scale_record, (), ())
    )
    monkeypatch.setattr(
        calibration,
        "calibrate_finite_differences",
        lambda *_a: (fd, {"status": "SELECTED", "policy": fd}),
    )
    monkeypatch.setattr(
        calibration, "select_svd_driver", lambda *_a: ("gesvd", driver_record, ())
    )
    monkeypatch.setattr(
        calibration,
        "select_rank_policy",
        lambda *_a: (
            (0.0, provisional),
            {"status": "SELECTED", "policy": (0.0, provisional)},
        ),
    )
    b2 = {
        "case": "B2",
        "passed": True,
        "perturbation": {"ranks": (1, 2), "unstable_modes": (1,)},
    }

    def rank_truth(
        _scales: object, policy: Any, *_a: object
    ) -> tuple[bool, tuple[dict[str, object], ...]]:
        seen["phase_ab"] = policy
        return True, (b2, {"case": "B3", "passed": True})

    monkeypatch.setattr(calibration, "validate_rank_truth_cases", rank_truth)
    monkeypatch.setattr(
        calibration,
        "_spectral_observation",
        lambda name, *_a: (object(), {"case": name, "status": "ACQUIRED"}),
    )
    monkeypatch.setattr(
        calibration,
        "_derive_case",
        lambda name, *_a, **_k: (
            SimpleNamespace(constrained_propagation=object()),
            {"case": name, "passed": True},
        ),
    )
    monkeypatch.setattr(
        calibration,
        "select_weak_policy",
        lambda *_a: (final, {"status": "SELECTED", "value": final}),
    )
    monkeypatch.setattr(
        calibration,
        "select_cluster_policy",
        lambda *_a: (2.0**-32, {"status": "SELECTED", "value": 2.0**-32}),
    )
    monkeypatch.setattr(
        calibration,
        "select_conditioning_policy",
        lambda *_a: (2.0**20, {"status": "SELECTED", "value": 2.0**20}),
    )
    monkeypatch.setattr(
        calibration,
        "select_correlation_policy",
        lambda *_a: (64.0, {"status": "SELECTED", "value": 64.0}),
    )

    def composed(
        policy: Any, *_a: object
    ) -> tuple[bool, tuple[dict[str, object], ...]]:
        seen["composed"] = policy
        return True, tuple(
            {"case": name, "passed": True} for name in calibration.COMPOSED_CASES
        )

    monkeypatch.setattr(calibration, "validate_composed_cases", composed)
    monkeypatch.setattr(
        calibration,
        "run_holdouts",
        lambda policy, *_a: (
            seen.__setitem__("holdout", policy) is None,
            tuple({"case": name, "passed": True} for name in calibration.HOLDOUT_CASES),
        ),
    )
    result = calibration.acquire_canonical("environment")
    policy = cast("dict[str, object]", result["policy"])
    phases = cast("dict[str, object]", result["selected_phase_policies"])
    metrics = cast("dict[str, dict[str, object]]", result["decisive_metrics"])
    assert (
        seen["phase_ab"].weak_relative_tolerance == provisional
        and provisional in calibration.WEAK_GRID
        and final in calibration.WEAK_GRID
        and final != provisional
    )
    assert (
        seen["composed"].weak_relative_tolerance
        == seen["holdout"].weak_relative_tolerance
        == policy["weak_relative_tolerance"]
        == phases["weak"]
        == metrics["weak"]["value"]
        == final
    )
    correlation_record = {
        "status": "UNSUPPORTED",
        "qualified_count": 0,
        "rejected_neighbors": (
            {
                "candidate": 0.0,
                "status": "DISQUALIFIED",
                "reasons": ("raw_truth_mismatch",),
            },
        ),
    }
    monkeypatch.setattr(
        calibration, "select_correlation_policy", lambda *_a: (None, correlation_record)
    )
    monkeypatch.setattr(
        calibration,
        "validate_composed_cases",
        lambda *_a: (_ for _ in ()).throw(
            AssertionError("composed validation must not run")
        ),
    )
    monkeypatch.setattr(
        calibration,
        "run_holdouts",
        lambda *_a: (_ for _ in ()).throw(AssertionError("holdouts must not run")),
    )
    terminal = calibration.acquire_canonical("environment")
    terminal_phases = cast("dict[str, object]", terminal["selected_phase_policies"])
    terminal_metrics = cast("dict[str, Any]", terminal["decisive_metrics"])
    assert terminal["status"] == "UNSUPPORTED_OR_SATURATED" and terminal_phases == {
        "scales": scales,
        "finite_difference": fd,
        "svd_driver": "gesvd",
        "rank": (0.0, provisional),
        "weak": final,
        "cluster": 2.0**-32,
        "conditioning": 2.0**20,
        "correlation": None,
    }
    assert (
        terminal_metrics["scale"] is scale_record
        and terminal_metrics["driver"] is driver_record
        and terminal_metrics["rank"]["typed_truth_cases"][0]["perturbation"]
        == b2["perturbation"]
        and terminal_metrics["correlation"] is correlation_record
    )
    assert (
        terminal["policy"]
        == {"status": "UNAVAILABLE", "reason": "incomplete_phase_c_selection"}
        and terminal["composed"]["status"] == "NOT_RUN"
        and tuple(item["case"] for item in terminal["holdouts"])
        == calibration.HOLDOUT_CASES
        and all(item["status"] == "NOT_RUN" for item in terminal["holdouts"])
        and terminal["compatibility"]
        == {"status": "NOT_RUN", "reason": "no_qualified_complete_policy"}
    )
    monkeypatch.setattr(
        calibration,
        "select_correlation_policy",
        lambda *_a: (
            64.0,
            {"case": "C5", "status": "SELECTED", "passed": True, "value": 64.0},
        ),
    )
    composed_failure = ({"case": "A1", "passed": False},)
    monkeypatch.setattr(
        calibration,
        "validate_composed_cases",
        lambda *_a: (False, composed_failure),
    )
    failed = calibration.acquire_canonical("environment")
    assert (
        failed["status"] == "COMPOSED_VALIDATION_FAILED"
        and failed["selected_phase_policies"]["correlation"] == 64.0
        and failed["decisive_metrics"]["correlation"]["passed"] is True
        and failed["composed"] == composed_failure
        and all(item["status"] == "NOT_RUN" for item in failed["holdouts"])
        and failed["compatibility"]
        == {"status": "NOT_RUN", "reason": "composed_validation_failed"}
        and isinstance(failed["policy_digest"], str)
    )


def test_both_scalings_and_h3_h5_use_real_typed_evidence() -> None:
    policy = selected_policy(2.0**-20)
    for scaling in (
        calibration.ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES,
        calibration.ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE,
    ):
        evidence, record = calibration._derive_case(
            "A3", policy, "environment", empty_counts(), "canonical", scaling=scaling
        )
        assert (
            record["passed"]
            and record["m"] == 4
            and record["n"] == 2
            and record["g"] == 1
            and record["nu"] == 1
            and record["covariance_truth"] is None
        )
        assert (
            evidence.rank_diagnostic is not None
            and evidence.source_policy.residual_variance_scaling is scaling
        )
    h3, h3_record = calibration._derive_case(
        "H3",
        policy,
        "environment",
        empty_counts(),
        "canonical",
        scaling=calibration.ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE,
    )
    h5_counts = empty_counts()
    h5, _record = calibration._derive_case(
        "H5",
        selected_policy(),
        "environment",
        h5_counts,
        "canonical",
        numerical_partial=True,
    )
    assert (
        h3_record["passed"]
        and h3.residual_jacobian is not None
        and h3.rank_diagnostic is not None
    )
    assert (
        h5.constraint_jacobian is not None
        and h5.constrained_propagation is not None
        and h5.constrained_marginal_errors is not None
    )
    zero = h5.constraint_jacobian.output_ids.index("ZERO")
    assert h5.constraint_jacobian.matrix[zero] == (0.0, 0.0) and all(
        value == 0.0 for value in h5.constrained_propagation.factor[zero]
    )
    assert all(
        h5.constrained_propagation.covariance[zero][i]
        == h5.constrained_propagation.covariance[i][zero]
        == 0.0
        for i in range(len(h5.constrained_propagation.output_ids))
    )
    assert (
        h5_counts["svd"] == 2
        and h5.constrained_correlations is not None
        and h5.correlations is not None
        and h5.constrained_propagation.claim("LOCAL_FIRST_ORDER_DEGENERACY")
        is calibration.ClaimState.VIOLATED
        and h5.constrained_propagation.claim("OUTPUT_FIRST_ORDER_NONDEGENERACY::ZERO")
        is calibration.ClaimState.VIOLATED
        and not h5.constrained_marginal_errors.scope_reportable
    )


def test_every_budget_refuses_before_the_first_over_budget_operation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    class Evaluator:
        def evaluate(self, _frame: object) -> None:
            calls.append("evaluation")

    counts = empty_counts()
    counts["evaluation"] = calibration.CEILINGS["canonical"]["evaluation"]
    owner = SimpleNamespace(charging=True, counts=counts, role="canonical")
    wrapped = calibration._CountingEvaluator(Evaluator(), owner)
    with pytest.raises(RuntimeError):
        wrapped.evaluate(object())
    monkeypatch.setattr(
        calibration.uncertainty_module, "svd", lambda *_a, **_k: calls.append("svd")
    )
    counts = empty_counts()
    counts["svd"] = calibration.CEILINGS["canonical"]["svd"] - 1
    with pytest.raises(RuntimeError), calibration._charged_svd(counts, "canonical"):
        calibration.uncertainty_module.svd(object())
        calibration.uncertainty_module.svd(object())
    monkeypatch.setattr(
        calibration, "_correlations", lambda **_k: calls.append("correlation")
    )
    counts = empty_counts()
    counts["correlation"] = calibration.CEILINGS["canonical"]["correlation"]
    with pytest.raises(RuntimeError):
        calibration._correlation_evidence(object(), object(), counts, "canonical")
    counts = empty_counts()
    counts["scalar"] = calibration.CEILINGS["canonical"]["scalar"]
    calibration._ACTIVE_COUNTS = counts
    calibration._ACTIVE_ROLE = "canonical"
    with pytest.raises(RuntimeError):
        calibration.counted_scientific_value(0.5, 2.0)
    assert calls == ["svd"]


def test_correlation_selection_uses_production_endpoint_calculation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expected_raw = (
        0.75,
        1.0 + 64.0 * 2.0**-52,
        -(1.0 + 64.0 * 2.0**-52),
        1.0 + 128.0 * 2.0**-52,
    )
    expected_outcomes = (
        "AVAILABLE",
        "ENDPOINT_CANONICALIZED",
        "ENDPOINT_CANONICALIZED",
        "CORRELATION_RANGE_VIOLATION",
    )
    counts = empty_counts()
    chosen, record = calibration.select_correlation_policy(selected_policy(), counts)
    assert (
        chosen == 64.0
        and record["case"] == "C5"
        and record["passed"] is True
        and record["truth_raw"] == expected_raw
        and record["outcomes"] == expected_outcomes
        and counts["correlation"] == len(calibration.CORRELATION_GRID)
    )

    original = calibration._expected_correlation_entries

    def mismatched(*args: object, **kwargs: object) -> tuple[tuple[object, ...], ...]:
        entries = [list(row) for row in original(*args, **kwargs)]
        entries[0][1] = replace(entries[0][1], raw_value=0.5)
        return tuple(tuple(row) for row in entries)

    monkeypatch.setattr(calibration, "_expected_correlation_entries", mismatched)
    failed_counts = empty_counts()
    chosen, failed = calibration.select_correlation_policy(
        selected_policy(), failed_counts
    )
    assert (
        chosen is None
        and failed["case"] == "C5"
        and failed["status"] == "TRUTH_PROBE_FAILED"
        and failed["passed"] is False
        and failed["reason"] == "raw_truth_mismatch"
        and "candidate_summaries" not in failed
        and failed_counts["correlation"] == 1
    )


def test_holdout_failure_is_irreversible(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        calibration,
        "_derive_case",
        lambda name, *_a, **_k: (calls.append(name), {"case": name, "passed": False}),
    )
    passed, records = calibration.run_holdouts(object(), "environment", empty_counts())
    assert (
        not passed
        and calls == ["H1"]
        and tuple(item["case"] for item in records) == ("H1",)
    )


def test_attestation_precedes_import(monkeypatch: pytest.MonkeyPatch) -> None:
    reloaded = importlib.reload(calibration)
    calls: list[str] = []
    monkeypatch.setattr(
        reloaded, "attest_lane", lambda *_: calls.append("attest") or {}
    )
    monkeypatch.setattr(reloaded, "_import_machinery", lambda: calls.append("import"))
    monkeypatch.setattr(
        reloaded,
        "validate_lane_records",
        lambda *_: (_ for _ in ()).throw(RuntimeError("stop")),
    )
    with pytest.raises(RuntimeError, match="stop"):
        reloaded.acquire("canonical", "digest")
    assert calls == ["attest", "import"]


def test_compatibility_is_bound_to_canonical_selected_policy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    policy = selected_policy()
    record = canonical_record(policy)
    monkeypatch.setattr(calibration, "validate_lane_records", lambda *_: "environment")
    rebuilt = calibration.policy_from_canonical_record(record)
    evidence, result = calibration._derive_case(
        "A2",
        rebuilt,
        "compatibility-environment",
        empty_counts("compatibility"),
        "compatibility",
    )
    assert (
        isinstance(rebuilt, calibration.UncertaintyPolicy)
        and result["passed"]
        and evidence.residual_jacobian is not None
    )
    runner_up = replace(policy, correlation_roundoff_multiplier=128.0)
    record["policy"] = calibration.policy_record(runner_up)
    record["policy_digest"] = hashlib.sha256(
        json.dumps(record["policy"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    with pytest.raises(RuntimeError, match="not the canonical selected policy"):
        calibration.policy_from_canonical_record(record)
    anonymous = canonical_record(policy)
    anonymous["composed"] = ({"passed": True},) * len(calibration.COMPOSED_CASES)
    with pytest.raises(RuntimeError, match="composed validation"):
        calibration.policy_from_canonical_record(anonymous)
    runner_up_scale = canonical_record(policy)
    runner_scales = cast(
        "dict[str, float]",
        cast("dict[str, Any]", runner_up_scale["selected_phase_policies"])["scales"],
    )
    runner_scales[calibration.BC_SCALE_FAMILY] = 2.0
    runner_policy = replace(
        calibration.policy_from_record(
            cast("dict[str, object]", runner_up_scale["policy"])
        ),
        calibration_identity=f"{calibration.SPECIFICATION_ID}:{calibration._scale_catalogue_digest(runner_scales)}",
    )
    runner_up_scale["policy"] = calibration.policy_record(runner_policy)
    runner_up_scale["policy_digest"] = hashlib.sha256(
        json.dumps(
            runner_up_scale["policy"], sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    with pytest.raises(RuntimeError, match="not the canonical selected policy"):
        calibration.policy_from_canonical_record(runner_up_scale)
    failed_driver = canonical_record(policy)
    cast(
        "dict[str, Any]",
        cast("dict[str, Any]", failed_driver["decisive_metrics"])["driver"],
    )["cases"] = (
        {"case": "B1", "status": "QUALIFIED"},
        {"case": "B2", "status": "DISQUALIFIED"},
        {"case": "B3", "status": "QUALIFIED"},
    )
    with pytest.raises(RuntimeError, match="typed phase evidence"):
        calibration.policy_from_canonical_record(failed_driver)
    monkeypatch.setattr(
        calibration,
        "_derive_case",
        lambda name, *_a, **_k: (
            SimpleNamespace(constrained_propagation=None),
            {"case": name, "passed": True},
        ),
    )
    monkeypatch.setattr(
        calibration, "_h4_spectral_case", lambda *_a: {"case": "H4", "passed": True}
    )
    monkeypatch.setattr(
        calibration,
        "_correlation_candidate_record",
        lambda *_a: ({"status": "QUALIFIED", "candidate": 64.0}, True),
    )
    replay = calibration.run_compatibility_replay(
        rebuilt,
        {calibration.BC_SCALE_FAMILY: 1.0},
        "environment",
        empty_counts("compatibility"),
    )
    assert (
        replay["status"] == "COMPATIBLE"
        and tuple(
            item["case"]
            for item in cast("tuple[dict[str, object], ...]", replay["results"])
        )
        == calibration.COMPATIBILITY_CASES
        and replay["retuning"] is False
    )
