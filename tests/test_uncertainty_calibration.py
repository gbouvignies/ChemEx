"""Prospective checks for the frozen issue #601 calibration specification."""

from __future__ import annotations

import importlib
import inspect
from collections.abc import Mapping
from typing import Any, cast

import pytest

from chemex.numerical_lanes import (
    LaneAttestation,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)
from tests.qualification import capture_uncertainty_calibration as calibration


def test_import_defers_scientific_calibration_machinery() -> None:
    reloaded = importlib.reload(calibration)
    assert not reloaded._SCIENTIFIC_MACHINERY_IMPORTED
    assert reloaded.UncertaintyPolicy is None
    assert reloaded.migration_core_authority_selection is None


def test_wayfinder_candidate_grids_are_literal_and_complete() -> None:
    assert tuple(range(-8, 9)) == calibration.SCALE_EXPONENTS
    assert (0.0,) + tuple(
        2.0**-power for power in range(8, 49, 2)
    ) == calibration.RELATIVE_STEP_TOLERANCES
    assert tuple(2.0**power for power in range(13)) == calibration.ROUNDOFF_MULTIPLIERS
    assert calibration.STEP_EXTENTS == (0, 2, 4, 8, 12, 16)
    assert calibration.SVD_DRIVERS == ("gesvd", "gesdd")
    assert (0.0,) + tuple(
        2.0**-power for power in range(16, 53, 2)
    ) == calibration.RANK_TOLERANCES
    assert (0.0,) + tuple(
        2.0**-power for power in range(2, 53, 2)
    ) == calibration.WEAK_MODE_THRESHOLDS
    assert calibration.CLUSTER_THRESHOLDS == calibration.WEAK_MODE_THRESHOLDS
    assert (1.0,) + tuple(
        2.0**power for power in range(2, 53, 2)
    ) == calibration.CONDITIONING_LIMITS
    assert (0.0,) + tuple(
        2.0**power for power in range(13)
    ) == calibration.CORRELATION_MULTIPLIERS

    assert calibration.CANDIDATE_COUNTS == {
        "characteristic_scale_per_family": 17,
        "finite_difference": 10_296,
        "svd_driver": 2,
        "rank": 400,
        "weak_mode": 27,
        "cluster": 27,
        "conditioning": 27,
        "correlation": 14,
    }
    policies = calibration.finite_difference_policies()
    assert sum(1 for _ in policies) == 10_296
    assert calibration.finite_difference_policies().__next__() == (0.0, 1.0, 0, 0)
    assert calibration.FINITE_DIFFERENCE_SEMANTICS["independent_extents"] is True
    assert calibration.FINITE_DIFFERENCE_SEMANTICS["exponents"] == (
        "range(-E_minus,E_plus+1)"
    )
    assert calibration.SCALE_SELECTION_SCOPE["fallback"] is None


def test_smaller_and_larger_extents_are_independent_grid_axes() -> None:
    policies = calibration.finite_difference_policies()
    assert (0.0, 1.0, 4, 12) in policies
    assert (0.0, 1.0, 12, 4) in policies


def test_characteristic_scale_catalogue_is_closed_and_has_no_fallback() -> None:
    supported = calibration.SUPPORTED_SCALE_FAMILIES
    unsupported = calibration.UNSUPPORTED_SCALE_FAMILIES
    assert len(supported) == 16
    assert {entry[1] for entry in supported} == {
        "dimensionless",
        "fraction",
        "rate",
        "frequency",
    }
    assert all(entry[3] in {"1", "s^-1", "ppm", "Hz"} for entry in supported)
    assert unsupported
    assert {entry[0] for entry in supported}.isdisjoint(unsupported)
    assert calibration.characteristic_scale_family("kab")[0] == (
        "directional_exchange_rate"
    )
    assert calibration.characteristic_scale_family("kex_ef")[0] == "exchange_rate"
    assert calibration.characteristic_scale_family("etaxy_a")[0] == "cross_relaxation"
    assert calibration.scale_candidates(100.0) == tuple(
        100.0 * 2.0**exponent for exponent in range(-8, 9)
    )
    with pytest.raises(KeyError, match="unsupported characteristic-scale family"):
        calibration.characteristic_scale_family("kd")
    with pytest.raises(KeyError, match="not catalogued"):
        calibration.characteristic_scale_family("invented_parameter")


def test_phase_ownership_is_irreversible_and_holdout_cannot_reopen_selection() -> None:
    assert calibration.PHASE_OWNERSHIP == (
        ("A0", "characteristic_scale"),
        ("A1", "finite_difference"),
        ("B1", "svd_driver"),
        ("B2", "rank_threshold"),
        ("C1", "weak_mode"),
        ("C2", "cluster"),
        ("C3", "conditioning"),
        ("C4", "correlation_endpoint"),
        ("F", "composed_policy"),
        ("H", "holdout_validation_only"),
    )
    selected = {"A0": "scale", "A1": "fd", "B1": "gesvd"}
    assert calibration.advance_phase(selected, "B2", "rank") == {
        **selected,
        "B2": "rank",
    }
    with pytest.raises(RuntimeError, match="cannot retune"):
        calibration.advance_phase(selected, "A1", "other-fd")
    assert calibration.holdout_decision((True,) * 6) == "QUALIFIED"
    assert calibration.holdout_decision((True, True, False, True, True, True)) == (
        "HOLDOUT_FAILED_POLICY_UNAVAILABLE"
    )
    with pytest.raises(ValueError, match="exactly six"):
        calibration.holdout_decision((True,) * 5)


def test_saturation_and_non_monotone_outcomes_fail_closed() -> None:
    assert calibration.phase_decision([], edge_only=False, ambiguous=False) == (
        "UNSUPPORTED"
    )
    assert calibration.phase_decision(
        [{"candidate": "edge"}], edge_only=True, ambiguous=False
    ) == ("GRID_SATURATED")
    assert calibration.phase_decision(
        [{"candidate": "x"}], edge_only=False, ambiguous=True
    ) == ("AMBIGUOUS_NON_MONOTONE")
    assert calibration.phase_decision(
        [{"candidate": "x"}], edge_only=False, ambiguous=False
    ) == ("SELECTABLE")


def test_resource_ledger_separates_canonical_selection_from_compatibility_replay() -> (
    None
):
    ledger = calibration.RESOURCE_LEDGER
    breakdown = calibration.RESOURCE_BREAKDOWN
    assert ledger["planned_maximum"] == {
        "evaluation_engine_requests": 26_360,
        "scalar_function_requests": 3_586,
        "svd_kernels": 26,
        "correlation_kernels": 18,
        "offline_candidate_comparisons": 12_000_000,
    }
    assert ledger["ceiling"] == {
        "evaluation_engine_requests": 28_000,
        "scalar_function_requests": 4_000,
        "svd_kernels": 32,
        "correlation_kernels": 24,
        "offline_candidate_comparisons": 12_000_000,
    }
    assert ledger["canonical_python_313"]["evaluation_engine_requests"] == 25_609
    assert ledger["compatibility_python_314"]["evaluation_engine_requests"] == 751
    assert ledger["compatibility_python_314"]["candidate_selection"] == 0
    assert ledger["compatibility_python_314"]["holdouts"] == 0
    assert (
        sum(
            value
            for name, value in breakdown["canonical_python_313"].items()
            if name.endswith("evaluation_engine")
        )
        == 25_609
    )
    assert (
        sum(
            value
            for name, value in breakdown["canonical_python_313"].items()
            if name.endswith("scalar_function")
        )
        == 3_320
    )
    assert (
        breakdown["compatibility_python_314"]["selected_policy_evaluation_engine"]
        == 751
    )


def test_truth_and_holdout_catalogues_are_frozen() -> None:
    assert tuple(calibration.TRUTH_CASES) == (
        "A1",
        "A2",
        "A3",
        "A4",
        "B1",
        "B2",
        "B3",
        "C1",
        "C2",
        "C3",
        "C4",
        "C5",
        "F1",
        "F2",
    )
    assert tuple(calibration.HOLDOUT_CASES) == ("H1", "H2", "H3", "H4", "H5", "H6")
    assert calibration.COMPATIBILITY_REPLAY_CASES == ("A1", "A2", "H4", "C5", "F2")
    assert calibration.METRIC_TOLERANCES["derivative_relative_envelope"] == 2.0**-24
    assert calibration.METRIC_TOLERANCES["truth_uncertainty_fraction"] == 1.0 / 16.0
    assert calibration.METRIC_TOLERANCES["boundary_zeta"] == 3.0
    assert calibration.METRIC_TOLERANCES["correlation_interior_epsilon"] == 8.0
    assert set(calibration.METRIC_FORMULAS) == {
        "derivative",
        "reliability",
        "covariance",
        "rank",
        "weak",
        "cluster",
        "conditioning",
        "correlation",
        "boundary",
        "normalization_dof",
        "constraint_partials",
        "propagation",
        "zero_gradient",
    }
    assert calibration.EVIDENCE_FIELDS == (
        "selected_phase_policies",
        "composed_uncertainty_policy_identity",
        "qualification_scope",
        "decisive_candidate_metrics",
        "rejected_neighbors",
        "holdout_results",
        "canonical_provenance",
    )


def test_lane_specific_acquisition_plan_forbids_compatibility_selection() -> None:
    canonical = calibration.prospective_record("canonical")
    compatibility = calibration.prospective_record(
        "compatibility", {"identity": "selected-on-canonical"}
    )
    canonical_plan = cast("tuple[tuple[str, str], ...]", canonical["acquisition_plan"])
    compatibility_plan = cast(
        "tuple[tuple[str, str], ...]", compatibility["acquisition_plan"]
    )
    compatibility_resources = cast(
        "Mapping[str, Mapping[str, int]]", compatibility["resource_ledger"]
    )
    assert canonical["status"] == "FROZEN_AWAITING_INDEPENDENT_REVIEW"
    assert tuple(step[0] for step in canonical_plan) == (
        "A0",
        "A1",
        "B1",
        "B2",
        "C",
        "F",
        "H",
    )
    assert tuple(step[0] for step in compatibility_plan) == ("P", "D")
    assert (
        compatibility_resources["compatibility_python_314"]["candidate_selection"] == 0
    )
    assert compatibility_resources["compatibility_python_314"]["holdouts"] == 0


def test_live_attestation_precedes_scientific_imports_and_uses_typed_records(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lane = canonical_lanes()[0]
    environment = RuntimeEnvironment(lane.semantics)

    def observe_current_process(
        cls: type[RuntimeEnvironment], image_digest: str, provenance_path: object = None
    ) -> RuntimeEnvironment:
        _ = cls, provenance_path
        assert image_digest == lane.semantics.image_digest
        return environment

    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(observe_current_process),
    )
    records = calibration.attest_lane("canonical", lane.semantics.image_digest)
    reconstructed_lane = NumericalLane.from_record(
        cast("Any", records["numerical_lane"])
    )
    attestation = LaneAttestation.from_record(cast("Any", records["lane_attestation"]))
    reconstructed_environment = RuntimeEnvironment.from_record(
        cast("Any", records["runtime_environment"])
    )
    assert reconstructed_lane == lane
    assert attestation.lane_identity == lane.identity
    assert attestation.environment_identity == reconstructed_environment.identity

    calls: list[str] = []

    def accept_attestation(role: str, image_digest: str) -> dict[str, object]:
        calls.append(f"attest:{role}:{image_digest}")
        return {}

    def import_machinery() -> None:
        calls.append("import")

    def reject_records(records: Mapping[str, object], role: str) -> None:
        _ = records
        calls.append(f"validate:{role}")
        raise RuntimeError("typed lane reconstruction rejected")

    monkeypatch.setattr(calibration, "attest_lane", accept_attestation)
    monkeypatch.setattr(calibration, "_import_scientific_machinery", import_machinery)
    monkeypatch.setattr(calibration, "validate_lane_records", reject_records)
    with pytest.raises(RuntimeError, match="typed lane reconstruction rejected"):
        calibration.acquire("canonical", "externally-observed-image-digest")
    assert calls == [
        "attest:canonical:externally-observed-image-digest",
        "import",
        "validate:canonical",
    ]
    assert inspect.signature(calibration.acquire).parameters[
        "image_digest"
    ].default is (inspect.Parameter.empty)


def test_compatibility_lane_cannot_select_or_run_holdouts() -> None:
    with pytest.raises(RuntimeError, match="selected canonical policy"):
        calibration.acquire("compatibility", "digest")
    with pytest.raises(ValueError, match="canonical or compatibility"):
        calibration.acquire("other", "digest")


def test_streaming_reducer_keeps_only_decisive_summary_and_neighbors() -> None:
    observations = (
        {
            "candidate": (0,),
            "qualified": False,
            "truth_error": 3.0,
            "cost": 1,
            "reasons": ("truth",),
        },
        {
            "candidate": (1,),
            "qualified": True,
            "truth_error": 1.0,
            "cost": 2,
            "reasons": (),
        },
        {
            "candidate": (2,),
            "qualified": True,
            "truth_error": 2.0,
            "cost": 3,
            "reasons": (),
        },
    )
    summary = calibration.reduce_candidates(iter(observations))
    assert summary["selected"] == (1,)
    assert summary["frontier"] == (((1,), 1.0, 2),)
    assert summary["rejected_neighbors"] == (
        ((0,), ("truth",)),
        ((2,), ("pareto_dominated",)),
    )
    assert "records" not in summary
