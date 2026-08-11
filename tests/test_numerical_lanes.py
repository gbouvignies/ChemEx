from __future__ import annotations

import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

import pytest

import chemex.numerical_lanes as numerical_lanes
from chemex.numerical_lanes import (
    ComparisonOutcome,
    ComparisonScope,
    CrossLaneNumericalPolicy,
    LaneAttestation,
    LaneAuthorityError,
    LaneSemantics,
    LiveLaneAuthority,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
    compare_values,
    comparison_scope,
)


def _semantics(*, compatibility: bool = False) -> LaneSemantics:
    python_version = "3.14.5" if compatibility else "3.13.5"
    python_digits = "314" if compatibility else "313"
    python_source_hash = (
        "9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc"
        if compatibility
        else "e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385"
    )
    differing = "2" if compatibility else "1"
    return LaneSemantics(
        python_implementation="CPython",
        python_version=python_version,
        python_abi=f"cpython-{python_digits}-x86_64-linux-gnu",
        python_source_hash=python_source_hash,
        python_executable_hash=differing * 64,
        platform="linux/amd64",
        platform_manifest=numerical_lanes._PLATFORM_MANIFEST,
        dependency_lock_hash=numerical_lanes._LOCKFILE_HASH,
        build_recipe_hash=numerical_lanes._recipe_hash(),
        build_context_hash="3" * 64,
        uv_version=numerical_lanes._UV_VERSION,
        uv_wheel_hash=numerical_lanes._UV_WHEEL_HASH,
        wheel_manifest_hash=differing * 64,
        os_package_manifest_hash="4" * 64,
        image_digest=f"sha256:{differing * 64}",
        numpy_version="2.5.1",
        numpy_installation_hash=differing * 64,
        scipy_version="1.18.0",
        scipy_installation_hash=differing * 64,
        numpy_blas=numerical_lanes._NUMPY_BLAS,
        numpy_lapack=numerical_lanes._NUMPY_BLAS,
        scipy_blas=(
            numerical_lanes._SCIPY_BLAS_314
            if compatibility
            else numerical_lanes._SCIPY_BLAS_313
        ),
        scipy_lapack=(
            numerical_lanes._SCIPY_LAPACK_314
            if compatibility
            else numerical_lanes._SCIPY_LAPACK_313
        ),
        openblas_configuration=numerical_lanes._OPENBLAS_RUNTIME,
        openblas_core="Haswell",
        loaded_library_manifest_hash=differing * 64,
        numpy_dispatch_restrictions=",".join(
            numerical_lanes._NUMPY_DISPATCH_RESTRICTIONS
        ),
        isa="x86-64-v3",
        workers=1,
        native_threads=1,
        floating_point_mode="binary64-round-nearest-gradual-underflow",
    )


def _lane(*, compatibility: bool = False) -> NumericalLane:
    return NumericalLane(
        (
            "compatibility-linux-amd64-python-3.14-v1"
            if compatibility
            else "canonical-linux-amd64-python-3.13-v1"
        ),
        "PYTHON_COMPATIBILITY" if compatibility else "CANONICAL_NUMERICAL",
        _semantics(compatibility=compatibility),
    )


def _attestation(
    monkeypatch: pytest.MonkeyPatch, lane: NumericalLane
) -> LiveLaneAuthority:
    environment = RuntimeEnvironment(lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: environment),
    )
    return lane.attest_current_process(lane.semantics.image_digest)


def _fabricated_evidence(lane_identity: str) -> LaneAttestation:
    environment_identity = "0" * 64
    method = "POST_IMPORT_CURRENT_PROCESS"
    identity_payload = {
        "kind": "numerical-lane-attestation",
        "record": [lane_identity, environment_identity, 1, 1, method],
        "schema_version": 2,
        "semantic_version": "chemex-numerical-lane-v2",
    }
    identity = hashlib.sha256(
        json.dumps(
            identity_payload,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
    ).hexdigest()
    return LaneAttestation.from_record(
        {
            "environment_identity": environment_identity,
            "identity": identity,
            "lane_identity": lane_identity,
            "method": method,
            "native_threads": 1,
            "schema_version": 2,
            "semantic_version": "chemex-numerical-lane-v2",
            "workers": 1,
        }
    )


def test_lane_identity_covers_complete_semantics_and_round_trips_strictly() -> None:
    lane = _lane()

    assert NumericalLane.from_record(lane.to_record()) == lane
    assert LaneSemantics.from_record(lane.semantics.to_record()) == lane.semantics
    assert (
        replace(
            lane,
            semantics=replace(lane.semantics, image_digest=f"sha256:{'f' * 64}"),
        ).identity
        != lane.identity
    )
    assert (
        replace(
            lane,
            semantics=replace(lane.semantics, loaded_library_manifest_hash="f" * 64),
        ).identity
        != lane.identity
    )

    malformed = lane.to_record()
    malformed["unexpected"] = True
    with pytest.raises(ValueError, match="not canonical"):
        NumericalLane.from_record(malformed)


def test_compatibility_lane_records_only_python_abi_structural_differences() -> None:
    python_313 = _lane()
    python_314 = _lane(compatibility=True)

    assert python_313.identity != python_314.identity
    python_313.validate_compatibility_lane(python_314)
    assert set(python_313.compatibility_delta(python_314)) == {
        "image_digest",
        "loaded_library_manifest_hash",
        "numpy_installation_hash",
        "python_abi",
        "python_executable_hash",
        "python_source_hash",
        "python_version",
        "scipy_blas",
        "scipy_installation_hash",
        "scipy_lapack",
        "wheel_manifest_hash",
    }

    drifted = replace(
        python_314,
        semantics=replace(python_314.semantics, os_package_manifest_hash="f" * 64),
    )
    with pytest.raises(LaneAuthorityError, match="os_package_manifest_hash"):
        python_313.validate_compatibility_lane(drifted)


def test_canonical_manifest_loader_is_strict_and_role_ordered(tmp_path: Path) -> None:
    python_313 = _lane()
    python_314 = _lane(compatibility=True)
    (tmp_path / "canonical-linux-amd64-python-3.13-v1.json").write_text(
        json.dumps(python_313.to_record()), encoding="ascii"
    )
    (tmp_path / "compatibility-linux-amd64-python-3.14-v1.json").write_text(
        json.dumps(python_314.to_record()), encoding="ascii"
    )

    assert canonical_lanes(tmp_path) == (python_313, python_314)


def test_shipped_lane_manifests_are_the_declared_compatibility_pair() -> None:
    python_313, python_314 = canonical_lanes()

    python_313.validate_compatibility_lane(python_314)
    assert python_313.role == "CANONICAL_NUMERICAL"
    assert python_314.role == "PYTHON_COMPATIBILITY"


def test_canonical_manifest_loader_rejects_duplicate_json_members(
    tmp_path: Path,
) -> None:
    duplicate = f'{{"identity":"{"0" * 64}","identity":"{"1" * 64}"}}'
    (tmp_path / "canonical-linux-amd64-python-3.13-v1.json").write_text(
        duplicate, encoding="ascii"
    )

    with pytest.raises(LaneAuthorityError, match="unavailable"):
        canonical_lanes(tmp_path)


def test_only_current_process_probe_can_mint_live_authority_and_evidence_round_trips(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lane = _lane()

    with pytest.raises(TypeError, match="minted only"):
        LiveLaneAuthority()

    authority = _attestation(monkeypatch, lane)
    evidence = LaneAttestation.from_record(authority.to_record())

    assert isinstance(authority, LiveLaneAuthority)
    assert not isinstance(evidence, LiveLaneAuthority)
    assert evidence.to_record() == authority.to_record()
    assert evidence.lane_identity == lane.identity
    assert evidence.environment_identity == RuntimeEnvironment(lane.semantics).identity


def test_process_registry_binding_cannot_be_mutated_in_place(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    authority = _attestation(monkeypatch, _lane())
    binding = numerical_lanes._LIVE_LANE_BINDINGS[authority]

    assert isinstance(binding, tuple)
    for field, value in (
        ("lane_identity", "a" * 64),
        ("lane_role", "PYTHON_COMPATIBILITY"),
        ("attestation_identity", "b" * 64),
        ("environment_identity", "c" * 64),
        ("workers", 8),
        ("native_threads", 8),
    ):
        with pytest.raises(AttributeError):
            object.__setattr__(binding, field, value)

    assert comparison_scope(authority, authority).kind == "WITHIN_LANE_BITWISE"


def test_post_import_attestation_rejects_any_actual_process_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lane = _lane()
    wrong_environment = RuntimeEnvironment(replace(lane.semantics, native_threads=2))
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: wrong_environment),
    )

    with pytest.raises(LaneAuthorityError, match="native_threads"):
        lane.attest_current_process(lane.semantics.image_digest)


def test_current_process_probe_fails_closed_outside_the_lane() -> None:
    with pytest.raises(LaneAuthorityError, match="provenance"):
        RuntimeEnvironment.from_current_process(f"sha256:{'0' * 64}")


def test_comparison_scopes_require_attested_lanes_and_round_trip(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    python_313 = _lane()
    python_314 = _lane(compatibility=True)
    attestation_313 = _attestation(monkeypatch, python_313)
    attestation_314 = _attestation(monkeypatch, python_314)

    within = comparison_scope(attestation_313, attestation_313)
    cross = comparison_scope(attestation_313, attestation_314)

    assert within.kind == "WITHIN_LANE_BITWISE"
    assert cross.kind == "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON"
    assert ComparisonScope.from_record(within.to_record()) == within
    with pytest.raises(TypeError, match="live current-process lane authority"):
        comparison_scope(python_313, python_313)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Unknown comparison scope"):
        ComparisonScope(
            "UNKNOWN",  # type: ignore[arg-type]
            python_313.identity,
            python_314.identity,
            LaneAttestation.from_record(attestation_313.to_record()).identity,
            LaneAttestation.from_record(attestation_314.to_record()).identity,
        )


def test_deserialized_evidence_cannot_create_an_authoritative_comparison_scope() -> (
    None
):
    evidence = _fabricated_evidence(_lane().identity)

    with pytest.raises(TypeError, match="live current-process lane authority"):
        comparison_scope(evidence, evidence)


def test_within_lane_replay_is_bitwise_and_never_accepts_tolerances(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lane = _lane()
    attestation = _attestation(monkeypatch, lane)
    scope = comparison_scope(attestation, attestation)
    changed = math.nextafter(1.0, math.inf)

    outcome = compare_values(scope, [1.0], [changed])

    assert not outcome.equivalent
    assert ComparisonOutcome.from_record(outcome.to_record()) == outcome
    assert not compare_values(scope, [0.0], [-0.0]).equivalent
    with pytest.raises(ValueError, match="tolerance"):
        compare_values(
            scope,
            [1.0],
            [1.0],
            policy=CrossLaneNumericalPolicy("profile-intensity-v1", 1.0e-6, 0.0),
        )


def test_cross_lane_comparison_is_structural_finite_and_content_identified(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    left_lane = _lane()
    right_lane = _lane(compatibility=True)
    scope = comparison_scope(
        _attestation(monkeypatch, left_lane),
        _attestation(monkeypatch, right_lane),
    )
    policy = CrossLaneNumericalPolicy("profile-intensity-v1", 1.0e-12, 0.0)
    changed = math.nextafter(1.0, math.inf)

    outcome = compare_values(scope, [1.0], [changed], policy=policy)

    assert outcome.equivalent
    assert outcome.policy_identity == policy.identity
    assert CrossLaneNumericalPolicy.from_record(policy.to_record()) == policy
    with pytest.raises(ValueError, match="explicit numerical policy"):
        compare_values(scope, [1.0], [changed])
    with pytest.raises(ValueError, match="matching structure"):
        compare_values(scope, [1.0], [[1.0]], policy=policy)  # type: ignore[list-item]
    with pytest.raises(ValueError, match="finite"):
        compare_values(scope, [float("nan")], [float("nan")], policy=policy)
