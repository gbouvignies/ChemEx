from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import pytest

import chemex.numerical_lanes as numerical_lanes
from chemex.numerical_lanes import (
    CrossLaneNumericalPolicy,
    LaneAuthorityError,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
    compare_values,
    comparison_scope,
)


def _environment(lane: NumericalLane) -> RuntimeEnvironment:
    return RuntimeEnvironment(
        python_implementation="CPython",
        python_version=lane.python_version,
        python_abi=lane.python_abi,
        python_source_hash=lane.python_source_hash,
        platform=lane.platform,
        platform_manifest=lane.platform_manifest,
        dependency_lock_hash=lane.dependency_lock_hash,
        build_recipe_hash=lane.build_recipe_hash,
        numerical_library=lane.numerical_library,
        isa=lane.isa,
        workers=lane.workers,
        native_threads=lane.native_threads,
        floating_point_mode=lane.floating_point_mode,
        imported_packages=lane.required_packages,
    )


def test_canonical_lanes_are_content_addressed_and_limit_compatibility_delta() -> None:
    python_313, python_314 = canonical_lanes()

    assert python_313.identity == canonical_lanes()[0].identity
    assert python_313.identity != python_314.identity
    assert python_313.platform == python_314.platform
    assert python_313.dependency_lock_hash == python_314.dependency_lock_hash
    assert python_313.numerical_library == python_314.numerical_library
    assert python_313.isa == python_314.isa
    assert python_313.workers == python_314.workers == 1
    assert python_313.native_threads == python_314.native_threads == 1
    assert python_313.floating_point_mode == python_314.floating_point_mode
    assert python_313.compatibility_delta(python_314) == {
        "python_abi": ("cp313", "cp314"),
        "python_source_hash": (
            "e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385",
            "9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc",
        ),
        "python_version": ("3.13.5", "3.14.5"),
    }


def test_post_import_attestation_requires_every_declared_claim() -> None:
    lane = canonical_lanes()[0]
    environment = _environment(lane)

    attestation = lane.attest(environment)

    assert attestation.lane_identity == lane.identity
    assert attestation.environment_identity == environment.identity
    assert attestation.identity != lane.identity

    unavailable = replace(environment, imported_packages=("numpy==2.5.1",))
    with pytest.raises(LaneAuthorityError, match="post-import package"):
        lane.attest(unavailable)

    wrong_threads = replace(environment, native_threads=2)
    with pytest.raises(LaneAuthorityError, match="native numerical threads"):
        lane.attest(wrong_threads)

    wrong_manifest = replace(environment, platform_manifest=f"debian@sha256:{'0' * 64}")
    with pytest.raises(LaneAuthorityError, match="platform manifest"):
        lane.attest(wrong_manifest)

    with pytest.raises(LaneAuthorityError):
        lane.attest_current_process()


def test_post_import_probe_normalizes_the_linux_amd64_platform_name(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numerical_lanes.sys, "platform", "linux")
    monkeypatch.setattr(numerical_lanes.platform_module, "machine", lambda: "x86_64")

    environment = RuntimeEnvironment.from_current_process()

    assert environment.platform == "linux/amd64"


def test_post_import_probe_reads_image_provenance_not_launcher_environment(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    lane = canonical_lanes()[0]
    provenance = tmp_path / "provenance.json"
    provenance.write_text(
        json.dumps(
            {
                "build_recipe_hash": lane.build_recipe_hash,
                "dependency_lock_hash": lane.dependency_lock_hash,
                "platform_manifest": lane.platform_manifest,
                "python_source_hash": lane.python_source_hash,
            }
        ),
        encoding="ascii",
    )
    monkeypatch.setenv("CHEMEX_NUMERICAL_LANE_RECIPE_HASH", "0" * 64)

    environment = RuntimeEnvironment.from_current_process(provenance)

    assert environment.build_recipe_hash == lane.build_recipe_hash


def test_replay_scope_never_uses_cross_lane_comparison_for_within_lane_replay() -> None:
    python_313, python_314 = canonical_lanes()

    assert comparison_scope(python_313, python_313).kind == "WITHIN_LANE_BITWISE"
    assert comparison_scope(python_313, python_313).left_lane_identity == (
        comparison_scope(python_313, python_313).right_lane_identity
    )

    cross_lane = comparison_scope(python_313, python_314)
    assert cross_lane.kind == "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON"
    assert cross_lane.left_lane_identity != cross_lane.right_lane_identity

    changed = math.nextafter(1.0, math.inf)
    assert not compare_values(
        comparison_scope(python_313, python_313), [1.0], [changed]
    ).equivalent
    with pytest.raises(ValueError, match="tolerance"):
        compare_values(
            comparison_scope(python_313, python_313),
            [1.0],
            [1.0],
            policy=CrossLaneNumericalPolicy(1.0e-6, 0.0),
        )
    with pytest.raises(ValueError, match="explicit numerical policy"):
        compare_values(cross_lane, [1.0], [changed])
    assert compare_values(
        cross_lane,
        [1.0],
        [changed],
        policy=CrossLaneNumericalPolicy(1.0e-12, 0.0),
    ).equivalent
