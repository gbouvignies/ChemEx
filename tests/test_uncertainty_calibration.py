"""Verification for the executable prospective issue #601 calibration."""

from __future__ import annotations

import importlib
import math
from typing import Any, cast

import pytest

from tests.qualification import capture_uncertainty_calibration as calibration

# The direct verification mirrors the compact qualification artifact.
# fmt: off

def selected_policy() -> Any:
    calibration._import_machinery()
    return calibration.compose_uncertainty_policy(
        "A1", (1.0, 100.0), (2.0**-24, 256.0, 8, 8), "gesvd",
        (0.0, 2.0**-40), 2.0**-22, 2.0**-33, 2.0**19, 64.0,
    )


def test_frozen_wayfinder_populations_and_digest() -> None:
    assert calibration.COUNTS == (17, 10_296, 2, 400, 27, 27, 27, 14)
    assert len(tuple(calibration.finite_difference_policies())) == 10_296
    assert calibration.EXTENTS == (0, 2, 4, 8, 12, 16)
    assert (0.0, 1.0, 4, 12) in calibration.finite_difference_policies()
    assert (0.0, 1.0, 12, 4) in calibration.finite_difference_policies()
    assert len(calibration.RANK_GRID) ** 2 == 400
    assert calibration.frozen_digest() == "7a68b507192d22a725f290c003394ccd5291cd9dd91e4592bd1f65695d1b913c"


def test_catalogue_is_closed_and_generated_names_are_classified() -> None:
    assert len(calibration.SUPPORTED) == 16
    assert {item[1] for item in calibration.SUPPORTED} == {"dimensionless", "fraction", "rate", "frequency"}
    expected = {"kab": "directional_exchange_rate", "kex_ef": "exchange_rate", "r2_i_f": "transverse_relaxation", "dwm_i_af": None, "tauc_f": None}
    for name, family in expected.items():
        if family is None:
            with pytest.raises(KeyError, match="catalogued but unsupported"):
                calibration.classify_scale_name(name)
        else:
            assert calibration.classify_scale_name(name)[0] == family
    for name in calibration.UNSUPPORTED:
        with pytest.raises(KeyError, match="catalogued but unsupported"):
            calibration.classify_scale_name(name)
    with pytest.raises(KeyError, match="not catalogued"):
        calibration.classify_scale_name("invented_parameter")


def test_truth_functions_are_independent_and_executable() -> None:
    calibration._import_machinery()
    anchor = calibration._setup("A1")[0]
    step = 1.0e-6
    for column in range(2):
        plus, minus = list(anchor), list(anchor)
        plus[column] += step
        minus[column] -= step
        observed = tuple((a - b) / (2.0 * step) for a, b in zip(calibration._calculation("A1", tuple(plus), 1.0), calibration._calculation("A1", tuple(minus), 1.0), strict=True))
        expected = tuple(row[column] for row in calibration.truth_jacobian("A1"))
        assert observed == pytest.approx(expected, rel=2.0e-9, abs=2.0e-9)
    assert calibration.scientific_partials(0.5, 2.0) == pytest.approx((math.exp(0.25) / 2.0, -0.5 * math.exp(0.25) / 4.0 + 0.8))
    assert all(calibration.truth_jacobian(name) for name in ("A2", "A3", "F2", "H1", "H2", "H3", "H4", "H5"))


def test_explicit_selectors_saturate_and_resources_refuse() -> None:
    assert calibration._smallest_threshold((1.0, 2.0), lambda value: value == 2.0, 2.0)[1]["status"] == "GRID_SATURATED"
    counts = dict(calibration.CEILINGS["canonical"])
    with pytest.raises(RuntimeError, match="ceiling exceeded"):
        calibration._spend(counts, "canonical", "svd")


def test_holdout_failure_is_irreversible(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    def fail_first(name: str, *_args: object, **_kwargs: object) -> tuple[None, dict[str, object]]:
        calls.append(name)
        return None, {"case": name, "passed": False}
    monkeypatch.setattr(calibration, "_derive_case", fail_first)
    passed, records = calibration.run_holdouts(object(), "environment", dict.fromkeys(calibration.CEILINGS["canonical"], 0))
    assert not passed and calls == ["H1"] and tuple(item["case"] for item in records) == ("H1",)


def test_attestation_precedes_calibration_import(monkeypatch: pytest.MonkeyPatch) -> None:
    reloaded = importlib.reload(calibration)
    assert not reloaded._MACHINERY_IMPORTED
    calls: list[str] = []
    monkeypatch.setattr(reloaded, "attest_lane", lambda *_: calls.append("attest") or {})
    monkeypatch.setattr(reloaded, "_import_machinery", lambda: calls.append("import"))
    def stop(*_args: object) -> str:
        calls.append("validate")
        raise RuntimeError("stop")
    monkeypatch.setattr(reloaded, "validate_lane_records", stop)
    with pytest.raises(RuntimeError, match="stop"):
        reloaded.acquire("canonical", "digest")
    assert calls == ["attest", "import", "validate"]


def test_real_typed_policy_compatibility_replay_and_compact_reconstruction() -> None:
    policy = selected_policy()
    rebuilt = calibration.policy_from_record(calibration.policy_record(policy))
    assert rebuilt.identity == policy.identity
    counts = dict.fromkeys(calibration.CEILINGS["compatibility"], 0)
    result = calibration.run_compatibility_replay(rebuilt, "qualification-environment", counts)
    assert result["status"] == "COMPATIBLE" and result["retuning"] is False
    results = cast("tuple[dict[str, object], ...]", result["results"])
    assert tuple(item["case"] for item in results) == ("A1", "A2", "F2", "H4", "C5")
    assert counts == {"evaluation": 751, "scalar": 266, "svd": 4, "correlation": 4, "offline": 0}
