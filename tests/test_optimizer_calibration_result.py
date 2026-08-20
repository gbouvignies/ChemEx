"""Qualification disposition for the single prospective #602 v2 occurrence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import cast

import tests.qualification.capture_optimizer_calibration_v2 as calibration

RESULT = Path(__file__).parent / "fixtures/canonical_optimizer_calibration_v2.json"


def test_authoritative_v2_result_is_exactly_content_identified() -> None:
    content = RESULT.read_bytes()
    assert hashlib.sha256(content).hexdigest() == (
        "1eb81770df3ccac0e0cf247c409439a51942836f4354aac62a233c2e31489e0d"
    )
    record = json.loads(content)
    identity = record.pop("identity")

    assert (
        identity == "44ac3bc386e080f8b6d854bd2c9a50f7717b6df2634cfedcd9ce7b82199b341b"
    )
    assert identity == calibration._identity("canonical-optimizer-calibration", record)


def test_authoritative_v2_result_uses_the_frozen_source_lane_and_specification() -> (
    None
):
    record = json.loads(RESULT.read_bytes())

    assert record["specification"]["identity"] == (
        "c8db7972f1d43262b3935e33c2b711d6f392b5f3fcb9e78339079de31ca0c9c3"
    )
    assert record["source"]["manifest_sha256"] == (
        "8235f78c16e587302f692a2f6f407cd239a5893e912474c2492267e35ba18f23"
    )
    assert record["source"]["implementation_wheel_sha256"] == (
        "fc210667a013886c43b05cce3d84d4a83d2e59df6bdfd1cc8acc901fc5e842c2"
    )
    assert record["canonical_lane"]["numerical_lane"]["identity"] == (
        calibration.EXPECTED_CANONICAL_LANE_IDENTITY
    )
    assert record["authoritative_acquisition"] == {
        "adaptive_extension": False,
        "retuned": False,
        "retry_count": 0,
        "round": 1,
        "specification_commit": "6b13bad304cdce4a0e7944c3a3ce8befd16cc719",
        "v1_observations_used": False,
        "version": 2,
    }


def test_authoritative_v2_result_is_bounded_unsupported_and_adds_no_claims() -> None:
    record = json.loads(RESULT.read_bytes())
    truth = cast("dict[str, object]", record["truth_probes"])
    strata = cast("dict[str, dict[str, object]]", record["strata"])

    assert truth["status"] == "ARCHITECTURAL_FAILURE"
    assert truth["stage"] == "truth_probes"
    assert cast("dict[str, object]", truth["qualification"])["status"] == (
        "UNSUPPORTED"
    )
    assert strata["coupled_grid_trf"]["status"] == "ARCHITECTURAL_FAILURE"

    grouped = strata["decomposed_grouped_grid_trf"]
    assert cast("dict[str, object]", grouped["selection"]) == {
        "budget_per_component": 16,
        "status": "SELECTED",
    }
    candidates = cast("list[dict[str, object]]", grouped["candidates"])
    assert [candidate["budget_per_component"] for candidate in candidates] == [
        16,
        32,
        64,
        128,
    ]
    assert all(candidate["status"] == "PASS" for candidate in candidates)
    assert cast("dict[str, object]", grouped["qualification"]) == {
        "holdout_passed": True,
        "independent_holdout": True,
        "replay_passed": True,
        "resource_ceiling_passed": True,
        "status": "UNSUPPORTED",
        "truth_authority_matched": False,
    }
    assert cast("dict[str, object]", grouped["holdout"])["case"] == "RELAXATION_NZ"

    assert strata["routine_direct_trf"]["qualification"] == {
        "reason": (
            "routine Direct TRF representative case lacks pre-existing case-specific "
            "truth authority; v1 observations are ineligible"
        ),
        "status": "UNSUPPORTED",
    }
    assert strata["difficult_direct_trf"]["qualification"] == {
        "reason": (
            "difficult-start Direct TRF representative case lacks pre-existing "
            "case-specific truth authority; v1 observations are ineligible"
        ),
        "status": "UNSUPPORTED",
    }
    assert cast("dict[str, object]", strata["de_trf"]["selection"])["status"] == (
        "DIRECT_REFERENCE_UNQUALIFIED"
    )
    assert not any(
        cast("dict[str, object]", stratum.get("qualification", {})).get("status")
        == "QUALIFIED"
        for stratum in strata.values()
    )
