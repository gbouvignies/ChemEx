from __future__ import annotations

import hashlib
import io
import json
import tarfile
from argparse import Namespace
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

import pytest

import chemex.baselines as baselines
import chemex.migration_core as migration_core
import chemex.migration_core_lifecycle as migration_core_lifecycle
from chemex.baselines import (
    ApprovedAnchorName,
    ArtifactContent,
    BaselinePublisher,
    Occurrence,
    ResultBundle,
    ScientificAnchorPublisher,
    capture_approved_scientific_anchor_legacy_observation,
)
from chemex.migration_core import (
    MigrationCoreCoverageError,
    MigrationCoreCoverageManifest,
    compile_migration_core_coverage,
    migration_core_coverage_manifest,
)
from chemex.numerical_lanes import (
    LaneAttestation,
    LiveLaneAuthority,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)
from chemex.optimize.numerical_probes import NumericalProbeBaseline

_REPOSITORY_ROOT = Path(__file__).parents[1]
_ANCHOR_RELEASE = (
    _REPOSITORY_ROOT
    / "tests/fixtures/migration_core_canonical_anchor_release_v1.tar.xz"
)
_ANCHOR_RELEASE_HASH = (
    "95df63b1c22dcee32ee57b08e6921deb492920f3b4e4ec02604946cdbb77ea15"
)
_ANCHOR_RELEASE_IDENTITY = (
    "7a82791ecf7d44252cfb502ce0738e614ba95333dc2ffe925bb9d938699e18c4"
)


def _canonical_bytes(record: Mapping[str, object]) -> bytes:
    return json.dumps(
        record, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")


def _canonical_probe_baseline() -> NumericalProbeBaseline:
    record = json.loads(
        (
            Path(__file__).parent
            / "fixtures"
            / "canonical_numerical_probe_baseline.json"
        ).read_text(encoding="utf-8")
    )
    assert isinstance(record, dict)
    return NumericalProbeBaseline.from_record(record)


def _extract_anchor_release(destination: Path) -> None:
    migration_core._safe_extract_release(_ANCHOR_RELEASE, destination)


def _repack_anchor_release(root: Path, archive: Path) -> tuple[str, str]:
    release_path = root / "release.json"
    release = json.loads(release_path.read_bytes())
    for member in release["payload_members"]:
        path = root / member["path"]
        if path.is_file():
            content = path.read_bytes()
            member["content_hash"] = hashlib.sha256(content).hexdigest()
            member["size"] = len(content)
    release.pop("identity")
    release_identity = hashlib.sha256(_canonical_bytes(release)).hexdigest()
    release["identity"] = release_identity
    release_path.write_bytes(_canonical_bytes(release))
    with tarfile.open(archive, mode="w:xz", preset=0) as opened:
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            content = path.read_bytes()
            info = tarfile.TarInfo(path.relative_to(root).as_posix())
            info.size = len(content)
            info.mtime = 0
            info.mode = 0o444
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            opened.addfile(info, io.BytesIO(content))
    return hashlib.sha256(archive.read_bytes()).hexdigest(), release_identity


def _cpmg_release_records(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    release = json.loads((root / "release.json").read_bytes())
    anchor = next(
        item for item in release["anchors"] if item["evidence"] == "anchor:cpmg-15n-ip"
    )
    manifest_path = root / anchor["publication_root"] / "manifest.json"
    return release, json.loads(manifest_path.read_bytes())


def _write_cpmg_release_records(
    root: Path, release: Mapping[str, object], manifest: Mapping[str, object]
) -> None:
    anchors = release["anchors"]
    assert isinstance(anchors, list)
    anchor = next(item for item in anchors if item["evidence"] == "anchor:cpmg-15n-ip")
    manifest_path = root / anchor["publication_root"] / "manifest.json"
    manifest_path.write_bytes(_canonical_bytes(manifest))
    (root / "release.json").write_bytes(_canonical_bytes(release))


def _live_authority(
    monkeypatch: pytest.MonkeyPatch, lane: NumericalLane
) -> LiveLaneAuthority:
    environment = RuntimeEnvironment(lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: environment),
    )
    return lane.attest_current_process(lane.semantics.image_digest)


def _write_complete_anchor_output(
    output: Path, snapshot: Path, anchor_name: ApprovedAnchorName
) -> None:
    anchor = baselines._approved_anchor(anchor_name)
    inputs = baselines._capture_anchor_inputs(anchor, snapshot)
    inventory = baselines._anchor_artifact_inventory(anchor, inputs)
    roles = inventory["required_roles"]
    assert isinstance(roles, list)
    for role in roles:
        assert isinstance(role, str)
        if role.startswith("legacy-output:"):
            path = output / role.removeprefix("legacy-output:")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(role.encode("ascii"))


def _published_migration_core(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> tuple[ScientificAnchorPublisher, dict[str, str]]:
    lane = canonical_lanes()[0]
    authority = _live_authority(monkeypatch, lane)
    publisher = ScientificAnchorPublisher(tmp_path / "evidence")
    names_by_first_experiment: dict[str, ApprovedAnchorName] = {
        "500mhz": "cpmg-15n-ip",
        "23hz": "cest-13c-label-cn",
        "cest_10hz_10p_1": "2st-binding",
        "1.25hz": "dcest-fifu-drd",
    }

    def complete_fake_run(args: Namespace, **_kwargs: object) -> None:
        assert isinstance(args.experiments, list)
        assert isinstance(args.output, Path)
        first = args.experiments[0]
        _write_complete_anchor_output(
            args.output, first.parents[1], names_by_first_experiment[first.stem]
        )

    monkeypatch.setattr(baselines, "run", complete_fake_run)
    anchors: tuple[tuple[ApprovedAnchorName, str], ...] = (
        ("cpmg-15n-ip", "examples/Experiments/CPMG_15N_IP"),
        ("cest-13c-label-cn", "examples/Experiments/CEST_13C_LABEL_CN"),
        ("2st-binding", "examples/Combinations/2stBinding"),
        ("dcest-fifu-drd", "examples/Experiments/DCEST_15N_3States"),
    )
    selections: dict[str, str] = {}
    for anchor_name, directory in anchors:
        published = capture_approved_scientific_anchor_legacy_observation(
            anchor_name,
            publisher=publisher,
            anchor_directory=Path(directory),
            lane_authority=authority,
            attempt_token=f"coverage-{anchor_name}",
        )
        selections[f"anchor:{anchor_name}"] = published.bundle.identity

    return publisher, selections


def test_versioned_manifest_covers_the_complete_migration_core_audit() -> None:
    manifest = migration_core_coverage_manifest()
    identifiers = {requirement.identifier for requirement in manifest.requirements}

    assert manifest.manifest_version == "migration-core-coverage-v2"
    assert len(manifest.anchors) == 4
    assert len(manifest.supporting_evidence) == 11
    assert {
        "migration-core.direct-trf.routine-behavior",
        "migration-core.grid-trf.decomposed-execution",
        "migration-core.de-trf.reduced-poor-start-dcest",
        "migration-core.grouped.aggregate-validation",
        "migration-core.covariance.constrained-propagation",
        "migration-core.resampling.serial-two-worker-replay",
        "migration-core.mcmc.truth-backed-statistical-comparison",
        "migration-core.failure.publication",
        "migration-core.operational.cache",
        "migration-core.performance.fixed-work-mcmc",
    } <= identifiers


def test_compiler_rejects_self_asserted_supporting_evidence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    publisher, selections = _published_migration_core(monkeypatch, tmp_path)
    manifest = migration_core_coverage_manifest()
    labels = {
        expected.evidence: {"evidence": expected.evidence}
        for expected in manifest.supporting_evidence
    }

    with pytest.raises(MigrationCoreCoverageError, match="not eligible"):
        compile_migration_core_coverage(manifest, publisher, selections, labels)


def test_numerical_probe_requires_the_frozen_reference_manifest(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    baseline = _canonical_probe_baseline()

    assert migration_core._validated_numerical_probe_identity(baseline) == (
        "NumericalProbeBaseline",
        baseline.identity,
    )
    monkeypatch.setattr(
        migration_core,
        "CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY",
        "0" * 64,
    )
    assert migration_core._validated_numerical_probe_identity(baseline) is None


def test_numerical_probe_rejects_self_consistent_alternate_authority() -> None:
    baseline = _canonical_probe_baseline()
    altered = NumericalProbeBaseline(
        baseline.definitions,
        baseline.artifacts,
        baseline.observed_lane_identity,
        baseline.observed_lane_role,
        "a" * 64,
        "b" * 64,
        baseline.reference_manifest_identity,
        "REFERENCE_MATCHED",
    )

    assert altered.manifest_identity == baseline.manifest_identity
    assert altered.identity != baseline.identity
    assert migration_core._validated_numerical_probe_identity(altered) is None


def test_numerical_probe_rejects_self_consistent_alternate_artifact() -> None:
    baseline = _canonical_probe_baseline()
    evidence = replace(
        baseline.artifacts[0].evidence,
        trajectory_fingerprint="a" * 64,
    )
    artifact = replace(baseline.artifacts[0], evidence=evidence)
    artifacts = (artifact, *baseline.artifacts[1:])
    capture_only = NumericalProbeBaseline(baseline.definitions, artifacts)
    altered = NumericalProbeBaseline(
        baseline.definitions,
        artifacts,
        baseline.observed_lane_identity,
        baseline.observed_lane_role,
        baseline.observed_attestation_identity,
        baseline.observed_environment_identity,
        capture_only.manifest_identity,
        "REFERENCE_MATCHED",
    )

    assert artifact.identity != baseline.artifacts[0].identity
    assert altered.identity != baseline.identity
    assert migration_core._validated_numerical_probe_identity(altered) is None


def test_selected_numerical_probe_fixture_hash_fails_closed() -> None:
    current = replace(
        migration_core.migration_core_current_release_selection(),
        probe_file_hash="a" * 64,
    )
    status = migration_core.MigrationCorePhasedStatus.from_bytes(
        (_REPOSITORY_ROOT / current.status_locator).read_bytes()
    )

    with pytest.raises(MigrationCoreCoverageError, match="probe hash does not match"):
        migration_core._load_selected_probe(_REPOSITORY_ROOT, current, status)


def test_repository_frozen_authority_selection_is_exact() -> None:
    selection = migration_core.migration_core_authority_selection()

    assert selection.identity == (
        "001b6b9e791e66677b265d8b1dd3c2d8151f4a9e0b33bdcfe5196782cecdd066"
    )
    assert selection.lane_identity == canonical_lanes()[0].identity
    assert selection.attestation_identity == (
        "3b3e0bc184826d61ec6652194486c907a5faa4c64b68ec03bbe60b63c660d687"
    )
    assert selection.environment_identity == (
        "cc5359f90df35ec9b60fd56e483911745209519ecce37d69e17c4edd6ea3604f"
    )


def test_durable_release_recursively_reconstructs_all_four_publications() -> None:
    with migration_core.resolve_migration_core_anchor_release(
        _ANCHOR_RELEASE,
        expected_archive_hash=_ANCHOR_RELEASE_HASH,
        expected_release_identity=_ANCHOR_RELEASE_IDENTITY,
    ) as resolved:
        assert resolved.release.repository_commit == (
            "2263bd9d162323c1b26949b3c3d7428f52b1c697"
        )
        assert resolved.release.source_archive_hash == (
            "20df597c93ecf965a582de907e62d95436c1dda581fb0ddd42dd5a3c48a1f50f"
        )
        assert {
            evidence: (
                published.case.identity,
                published.specification.identity,
                published.occurrence.identity,
                published.occurrence.lifecycle_identity,
                published.bundle.identity,
                published.bundle.manifest_identity,
                published.manifest_identity,
                len(published.bundle.members),
            )
            for evidence, published in resolved.published.items()
        } == {
            "anchor:2st-binding": (
                "35af1e0705653511fc3d20f7e9800deb9383af0236d965e34df85c5e16d670c6",
                "8acab6a7ef2010ff9187cf9203e68b092e1d03e74433214849961ad8ef52c839",
                "2c87cd4769195e59860646f380adb0d72fb05b6844b41c51e40969f922664b1a",
                "b3f38a19395e19b1eb70398d8a9657d60ac7dc24756fa9008fbb890f5b881a91",
                "bd7c8424f1b24b090d4e87946b45b017eb850450cd1ff6b4d582ed66757bbd44",
                "c37e4e427530646486c761d1a695e818eca1b25b8ea1fdea0ec8912c69be5dcb",
                "9bf207521c5c42c0317163fdd9841dda6ed3736b77ca9fe7f69a6c81d47757c9",
                896,
            ),
            "anchor:cest-13c-label-cn": (
                "e221428fe2bed8694f6400f5c10a0ae096763f008949c70eecdb685f5afe7b01",
                "16999ec4bfe695df11dd9752de9b275daaee1c59a4574c10e88539f85a673e8d",
                "b2aa92f23ed505402bd9419b717c8a245f06ebd525d8ff8be11fce3465f2f614",
                "e88400280b8389de7a7def4ca736273bf707b5908a8155575b98a5653246f620",
                "4dd0452ec492aaf920f24eea5b773723f240520bff3d3dbdf1ad5efad0b660db",
                "9ca652014856faa7b738997517653de8b5e69a4f4576eefaf46995796956e331",
                "74e325ad35dcc237abcdbf2434a1d11e78b6ce1bddace5d0b84aa074d5b95ac4",
                92,
            ),
            "anchor:cpmg-15n-ip": (
                "18d9a4fb78e300a474132313c898c78240be7aa6ae7eb2f3385d2fed69e6ebd3",
                "97cde21a4efa9a2bc6ab1bb458f517b3d7ca8db70641f8b6c4ce10072b0d90f9",
                "96fdb5d9f0425b92a892e7ea784d95481ebeebbfa62655785b79fc1ef3a3ba4b",
                "fcd8ed5417f2106bea646e680e64e70b85f0ade60e94886545ef4c196e1c6920",
                "cc8e5072746712f6fc1fd2f976739d1fe1d9d11c436ae9f9fdd6f34404a12eec",
                "76c8d4507f862db64d14e9f08244d7dc28d21c955f807031ea1728a2e18c49cd",
                "6f6f7afd96ec66c4852736928302229568ea77e6db2668d540effeb4e66a57c2",
                342,
            ),
            "anchor:dcest-fifu-drd": (
                "3c379d070b9aea9fdf9675511df717b4208f4d0d09979e62f304216a6147c241",
                "a8aad410db1c559faf0fc220d32e6aaa247ef45a9d03c6c4de3de9f9d2d8489f",
                "714f523c0fedb5a4bc422ddfa51e0a1fbf45bcb2ecc81444e1ce0a70ba87fbb5",
                "a0b389ee27338fec2e4c437d22fdb3c07a773625a308b9d85a2f33f6e9004638",
                "80d52a017bda14233ce5009f5dc22a8e76b113b8edcf28c6bb67845cc90368e3",
                "98620feb28660dddaaa974dcb411d7876c453c7a78064e4a40edbf8aa07579a2",
                "e01ed45d43be92bd59d53fca4299126aa2a8b686ff9c7773c079d954cc939cf8",
                24,
            ),
        }
        for anchor in resolved.release.anchors:
            assert (
                resolved.publisher.read(anchor.result_bundle_identity)
                == (resolved.published[anchor.evidence])
            )


@pytest.mark.parametrize(
    ("layer", "manifest_path", "release_field"),
    (
        ("case", ("case", "identity"), "case_identity"),
        (
            "execution specification",
            ("specification", "identity"),
            "specification_identity",
        ),
        ("occurrence", ("occurrence", "identity"), "occurrence_identity"),
        (
            "occurrence lifecycle",
            ("occurrence", "lifecycle_identity"),
            "occurrence_lifecycle_identity",
        ),
        ("result bundle", ("bundle", "identity"), "result_bundle_identity"),
        (
            "member manifest",
            ("bundle", "manifest_identity"),
            "member_manifest_identity",
        ),
        (
            "publication",
            ("manifest_identity",),
            "publication_manifest_identity",
        ),
    ),
)
def test_release_rejects_changed_recursive_identity_after_outer_rehash(
    tmp_path: Path,
    layer: str,
    manifest_path: tuple[str, ...],
    release_field: str,
) -> None:
    root = tmp_path / "release"
    root.mkdir()
    _extract_anchor_release(root)
    release, manifest = _cpmg_release_records(root)
    anchors = release["anchors"]
    assert isinstance(anchors, list)
    anchor = next(item for item in anchors if item["evidence"] == "anchor:cpmg-15n-ip")
    target = manifest
    for component in manifest_path[:-1]:
        nested = target[component]
        assert isinstance(nested, dict)
        target = nested
    target[manifest_path[-1]] = "a" * 64
    if layer == "member manifest":
        manifest["member_manifest_identity"] = "a" * 64
    anchor[release_field] = "a" * 64
    _write_cpmg_release_records(root, release, manifest)
    archive_hash, release_identity = _repack_anchor_release(
        root, tmp_path / "altered.tar.xz"
    )

    with (
        pytest.raises(MigrationCoreCoverageError, match="cannot be reconstructed"),
        migration_core.resolve_migration_core_anchor_release(
            tmp_path / "altered.tar.xz",
            expected_archive_hash=archive_hash,
            expected_release_identity=release_identity,
        ),
    ):
        pass


@pytest.mark.parametrize("member_kind", ("captured input", "output artifact"))
def test_release_rejects_changed_member_byte_after_outer_rehash(
    tmp_path: Path, member_kind: str
) -> None:
    root = tmp_path / "release"
    root.mkdir()
    _extract_anchor_release(root)
    release, manifest = _cpmg_release_records(root)
    anchor = next(
        item for item in release["anchors"] if item["evidence"] == "anchor:cpmg-15n-ip"
    )
    if member_kind == "captured input":
        target = next((root / anchor["input_root"]).iterdir())
    else:
        member = next(
            item
            for item in manifest["bundle"]["members"]
            if item["role"].startswith("legacy-output:")
        )
        target = root / anchor["publication_root"] / "members" / member["content_hash"]
    content = bytearray(target.read_bytes())
    content[0] ^= 1
    target.write_bytes(content)
    archive_hash, release_identity = _repack_anchor_release(
        root, tmp_path / "altered.tar.xz"
    )

    with (
        pytest.raises(MigrationCoreCoverageError),
        migration_core.resolve_migration_core_anchor_release(
            tmp_path / "altered.tar.xz",
            expected_archive_hash=archive_hash,
            expected_release_identity=release_identity,
        ),
    ):
        pass


def test_release_rejects_wrong_hash_missing_corrupt_and_missing_member(
    tmp_path: Path,
) -> None:
    with (
        pytest.raises(MigrationCoreCoverageError, match="hash does not match"),
        migration_core.resolve_migration_core_anchor_release(
            _ANCHOR_RELEASE,
            expected_archive_hash="a" * 64,
            expected_release_identity=_ANCHOR_RELEASE_IDENTITY,
        ),
    ):
        pass
    with (
        pytest.raises(MigrationCoreCoverageError, match="missing or unavailable"),
        migration_core.resolve_migration_core_anchor_release(
            tmp_path / "missing.tar.xz",
            expected_archive_hash=_ANCHOR_RELEASE_HASH,
            expected_release_identity=_ANCHOR_RELEASE_IDENTITY,
        ),
    ):
        pass

    corrupt = tmp_path / "corrupt.tar.xz"
    corrupt.write_bytes(_ANCHOR_RELEASE.read_bytes()[:1024])
    with (
        pytest.raises(MigrationCoreCoverageError, match="corrupt or truncated"),
        migration_core.resolve_migration_core_anchor_release(
            corrupt,
            expected_archive_hash=hashlib.sha256(corrupt.read_bytes()).hexdigest(),
            expected_release_identity=_ANCHOR_RELEASE_IDENTITY,
        ),
    ):
        pass

    root = tmp_path / "release"
    root.mkdir()
    _extract_anchor_release(root)
    release, manifest = _cpmg_release_records(root)
    anchor = next(
        item for item in release["anchors"] if item["evidence"] == "anchor:cpmg-15n-ip"
    )
    member = manifest["bundle"]["members"][0]
    (root / anchor["publication_root"] / "members" / member["content_hash"]).unlink()
    archive_hash, release_identity = _repack_anchor_release(
        root, tmp_path / "missing-member.tar.xz"
    )
    with (
        pytest.raises(MigrationCoreCoverageError, match="payload is not closed"),
        migration_core.resolve_migration_core_anchor_release(
            tmp_path / "missing-member.tar.xz",
            expected_archive_hash=archive_hash,
            expected_release_identity=release_identity,
        ),
    ):
        pass


def test_release_locator_cannot_escape_trusted_repository_root(tmp_path: Path) -> None:
    with pytest.raises(MigrationCoreCoverageError, match="repository-relative"):
        migration_core._repository_path(tmp_path, "../outside", "anchor release")


def test_current_status_is_derived_from_resolved_evidence() -> None:
    result = migration_core.compile_current_phased_migration_core_status(
        _REPOSITORY_ROOT
    )

    assert result.authority_selection_identity == (
        "001b6b9e791e66677b265d8b1dd3c2d8151f4a9e0b33bdcfe5196782cecdd066"
    )
    assert result.anchor_release_identity == _ANCHOR_RELEASE_IDENTITY
    assert len(result.anchor_evidence) == 4
    assert result.numerical_probe.identity == (
        "d11f9caa404b8ce3fa96e041659939446e205aae9376b08d6a137ed40cf0bdb0"
    )
    assert (
        len(result.eligible_requirement_ids),
        len(result.uncovered_requirement_ids),
    ) == (30, 21)
    failure_coverage = {
        identifier
        for identifier in result.eligible_requirement_ids
        if identifier.startswith("migration-core.failure.")
    }
    assert failure_coverage == set(migration_core_lifecycle.FAILURE_REQUIREMENTS)
    assert set(result.unqualified_requirement_ids) == {
        "migration-core.covariance.boundary",
        "migration-core.covariance.constrained-propagation",
        "migration-core.covariance.finite-difference-reliability",
        "migration-core.resampling.nucleus-bootstrap-truth-probe",
    }
    assert result.compiler_status == "FAILED_CLOSED"
    assert result.compiled_coverage_identity is None


def test_bounded_calibrations_qualify_only_their_demonstrated_claims() -> None:
    uncertainty = (
        _REPOSITORY_ROOT / "tests/fixtures/canonical_uncertainty_calibration_v3.json"
    ).read_bytes()
    resampling = (
        _REPOSITORY_ROOT / "tests/fixtures/canonical_resampling_calibration_v2.json"
    ).read_bytes()

    assert migration_core.qualified_migration_core_requirements(
        "execution:covariance-constrained-uncertainty", uncertainty
    ) == {
        "migration-core.covariance.analytic-normalization",
        "migration-core.covariance.both-scaling-policies",
        "migration-core.covariance.conditioning",
        "migration-core.covariance.correlation",
        "migration-core.covariance.propagation-degeneracy",
        "migration-core.covariance.rank",
    }
    assert migration_core.qualified_migration_core_requirements(
        "execution:resampling", resampling
    ) == {
        "migration-core.resampling.bootstrap-truth-probe",
        "migration-core.resampling.mc-truth-probe",
        "migration-core.resampling.serial-two-worker-replay",
    }
    assert not migration_core.qualified_migration_core_requirements(
        "execution:covariance-constrained-uncertainty", uncertainty + b" "
    )


def test_selected_supporting_evidence_hash_fails_closed() -> None:
    current = migration_core.migration_core_current_release_selection()
    altered = replace(current.supporting_evidence[0], file_hash="a" * 64)
    current = replace(
        current,
        supporting_evidence=(altered, *current.supporting_evidence[1:]),
    )

    with pytest.raises(
        MigrationCoreCoverageError, match="selection hash does not match"
    ):
        migration_core._load_selected_supporting_evidence(_REPOSITORY_ROOT, current)


def test_changed_requirement_selection_is_not_trusted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    status = migration_core.MigrationCorePhasedStatus.from_bytes(
        (
            _REPOSITORY_ROOT
            / migration_core.migration_core_current_release_selection().status_locator
        ).read_bytes()
    )
    selected = dict(status.selected_evidence)
    selected["anchor:cpmg-15n-ip"] = "a" * 64
    altered = replace(status, selected_evidence=tuple(sorted(selected.items())))
    authority = migration_core.migration_core_authority_selection()
    monkeypatch.setattr(
        migration_core,
        "_load_selected_authority_and_status",
        lambda repository_root, current: (authority, altered),
    )

    with pytest.raises(
        MigrationCoreCoverageError,
        match="requirement-to-evidence selection is incompatible",
    ):
        migration_core.compile_current_phased_migration_core_status(_REPOSITORY_ROOT)


def test_anchor_publications_retain_exact_attestation_and_environment_records(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    publisher, selections = _published_migration_core(monkeypatch, tmp_path)
    authority = migration_core.migration_core_authority_selection()

    for bundle_identity in selections.values():
        published = publisher.read(bundle_identity)
        members = {member.role: member for member in published.bundle.members}
        attestation = LaneAttestation.from_record(
            json.loads(
                (
                    published.location
                    / "members"
                    / members["environment:lane-attestation.json"].content_hash
                ).read_bytes()
            )
        )
        environment = RuntimeEnvironment.from_record(
            json.loads(
                (
                    published.location
                    / "members"
                    / members["environment:runtime-environment.json"].content_hash
                ).read_bytes()
            )
        )

        assert attestation.identity == authority.attestation_identity
        assert attestation.environment_identity == authority.environment_identity
        assert environment.identity == authority.environment_identity


def test_structurally_valid_alternate_anchor_authority_is_ineligible(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    publisher, selections = _published_migration_core(monkeypatch, tmp_path)
    original = publisher.read(selections["anchor:cpmg-15n-ip"])
    alternative_environment = RuntimeEnvironment(canonical_lanes()[1].semantics)
    alternative_attestation = LaneAttestation(
        canonical_lanes()[0].identity,
        alternative_environment.identity,
        1,
        1,
        "POST_IMPORT_CURRENT_PROCESS",
    )
    artifacts = []
    for member in original.bundle.members:
        content = (original.location / "members" / member.content_hash).read_bytes()
        if member.role == "environment:lane-attestation.json":
            content = _canonical_bytes(alternative_attestation.to_record())
        elif member.role == "environment:runtime-environment.json":
            content = _canonical_bytes(alternative_environment.to_record())
        artifacts.append(ArtifactContent(member.role, content))
    requested = Occurrence(
        original.specification.identity,
        original.case.identity,
        original.specification.implementation.identity,
        original.specification.lane_reference,
        alternative_attestation.identity,
        tuple(sorted(member.identity for member in original.case.inputs)),
        "alternate-valid-authority",
    )
    bundle = ResultBundle.create(
        requested.identity,
        original.specification.identity,
        original.specification.implementation,
        tuple(artifact.member for artifact in artifacts),
    )
    alternate_publisher = BaselinePublisher(tmp_path / "alternate")
    alternate_publisher.reserve(original.case, original.specification, requested)
    alternate = alternate_publisher.publish(
        original.case,
        original.specification,
        requested.succeeded(bundle),
        bundle,
        artifacts,
    )
    reconstructed = alternate_publisher.read(alternate.bundle.identity)
    expected = migration_core_coverage_manifest().anchors[2]

    assert reconstructed == alternate
    with pytest.raises(MigrationCoreCoverageError, match="exact authority mismatch"):
        migration_core._validate_anchor_evidence(expected, reconstructed)


def test_compiler_fails_closed_on_missing_unresolved_or_unavailable_evidence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    publisher, selections = _published_migration_core(monkeypatch, tmp_path)
    manifest = migration_core_coverage_manifest()
    supporting = {item.evidence: object() for item in manifest.supporting_evidence}
    missing = dict(supporting)
    missing.pop("probe:lifecycle-publication-safety")
    with pytest.raises(MigrationCoreCoverageError, match="missing required evidence"):
        compile_migration_core_coverage(manifest, publisher, selections, missing)

    unresolved = dict(selections)
    unresolved["anchor:dcest-fifu-drd"] = "a" * 64
    with pytest.raises(MigrationCoreCoverageError, match="unresolved or unavailable"):
        compile_migration_core_coverage(manifest, publisher, unresolved, supporting)

    with pytest.raises(MigrationCoreCoverageError, match="unavailable or not eligible"):
        compile_migration_core_coverage(manifest, publisher, selections, supporting)


def test_manifest_validator_rejects_incompatible_or_empty_coverage() -> None:
    manifest = migration_core_coverage_manifest()
    incompatible = json.loads(manifest.to_bytes())
    incompatible["manifest_version"] = "migration-core-coverage-v3"
    with pytest.raises(MigrationCoreCoverageError, match="incompatible version"):
        MigrationCoreCoverageManifest.from_bytes(_canonical_bytes(incompatible))

    uncovered = json.loads(manifest.to_bytes())
    uncovered["requirements"][0]["evidence"] = []
    with pytest.raises(MigrationCoreCoverageError, match="eligible evidence"):
        MigrationCoreCoverageManifest.from_bytes(_canonical_bytes(uncovered))

    truncated = json.loads(manifest.to_bytes())
    truncated["requirements"].pop()
    with pytest.raises(MigrationCoreCoverageError, match="fixed audited requirements"):
        MigrationCoreCoverageManifest.from_bytes(_canonical_bytes(truncated))


def test_canonical_capture_status_records_the_fail_closed_empirical_result() -> None:
    record = json.loads(
        (
            Path(__file__).parent
            / "fixtures"
            / "migration_core_canonical_capture_status_v1.json"
        ).read_text(encoding="ascii")
    )
    identity = record.pop("identity")
    anchors = {anchor["anchor"]: anchor for anchor in record["anchors"]}
    coverage = record["coverage"]

    assert hashlib.sha256(_canonical_bytes(record)).hexdigest() == identity
    assert record["status"] == "NOT_READY"
    assert record["source"]["repository_commit"] == (
        "09269a1d48e79d1684a69fa51ef3f603f7d23c95"
    )
    assert record["lane"] == {
        "attestation_identity": (
            "3b3e0bc184826d61ec6652194486c907a5faa4c64b68ec03bbe60b63c660d687"
        ),
        "attestation_method": "POST_IMPORT_CURRENT_PROCESS",
        "environment_identity": (
            "cc5359f90df35ec9b60fd56e483911745209519ecce37d69e17c4edd6ea3604f"
        ),
        "identity": (
            "4f1e7dab3384f88149fd33befb09ba3e96730dc336e9427b224a39a8a7167e4f"
        ),
        "image_digest": (
            "sha256:ed4f97e00bd6fe494f46772cc31a338d756b5bb5fe6e4e480d987af51af85550"
        ),
        "name": "canonical-linux-amd64-python-3.13-v1",
        "native_threads": 1,
        "workers": 1,
    }
    assert set(anchors) == {
        "cpmg-15n-ip",
        "cest-13c-label-cn",
        "2st-binding",
        "dcest-fifu-drd",
    }
    assert {
        name for name, anchor in anchors.items() if anchor["status"] == "ELIGIBLE"
    } == {"cpmg-15n-ip", "2st-binding"}
    assert anchors["cpmg-15n-ip"]["artifact_count"] == 341
    assert anchors["2st-binding"]["artifact_count"] == 895
    assert anchors["cest-13c-label-cn"]["result_bundle_identity"] is None
    assert anchors["dcest-fifu-drd"]["result_bundle_identity"] is None
    assert coverage["manifest_identity"] == (
        "c79d318b753eed1e863c6e5e99bbee7f7c1f7cc19ad8c00dc15acf983a6e69c8"
    )
    assert coverage["requirement_count"] == 51
    assert coverage["eligible_requirement_count"] == 4
    assert coverage["uncovered_requirement_count"] == 47
    assert len(coverage["uncovered_requirement_ids"]) == 47
    assert coverage["compiled_coverage_identity"] is None
    assert coverage["compiler"]["status"] == "FAILED_CLOSED"


def test_v2_status_is_historical_asserted_unresolved_evidence_only() -> None:
    content = (
        Path(__file__).parent
        / "fixtures"
        / "migration_core_canonical_coverage_status_v2.json"
    ).read_bytes()
    record = json.loads(content)
    identity = record.pop("identity")
    anchors = {anchor["anchor"]: anchor for anchor in record["anchors"]}
    coverage = record["coverage"]

    assert hashlib.sha256(_canonical_bytes(record)).hexdigest() == identity
    assert record["status"] == "READY_TO_REVIEW_AS_PHASED_590_INFRASTRUCTURE"
    assert record["completion_claim"] == "PHASED_INFRASTRUCTURE_ONLY"
    assert record["source"] == {
        "archive_sha256": (
            "7e5677cf87dee315f05c5ea9b56ad904578e5cc904f39964c2f6d16ab14a2a1e"
        ),
        "case_source_authority": {
            "authority_version": "case-source-authority-v1",
            "identity": (
                "2151f8713fe10a6f2bf6f9b6918c8aa4b1c1b503eb2ba3a3e925a7f35f0e4d7e"
            ),
            "lockfile_hash": (
                "f0fb2ffc7b1a5ecd1bf7ac43956fc4861b96c058d158948b68b4e97027a6086a"
            ),
            "schema_version": 1,
            "source_commit": "d5ed0c87e8ce7a7f17745feea346af4dfbae7ecf",
        },
        "execution_source_label": (
            "cc82cba558b885cdb8f26f173a08a09294195031+working-tree@"
            "7e5677cf87dee315f05c5ea9b56ad904578e5cc904f39964c2f6d16ab14a2a1e"
        ),
        "implementation_authority_identity": (
            "f284deeb88aba9041764d7a471ad962746f92c857f3f1dc0d96d6bafc1d10eac"
        ),
        "implementation_source_manifest_identity": (
            "d9829b9bb00d7403630d2416dafa0fcb2a33a8cfd7558f213d267b2081ac00a8"
        ),
        "repository_base_commit": "cc82cba558b885cdb8f26f173a08a09294195031",
        "tree_sha256": (
            "819447cbe245b015acd063992450fce3e7bdb811fb13a513233d0bb2075951bc"
        ),
    }
    assert record["lane"]["identity"] == canonical_lanes()[0].identity
    assert record["lane"]["attestation_identity"] == (
        "3b3e0bc184826d61ec6652194486c907a5faa4c64b68ec03bbe60b63c660d687"
    )
    assert record["lane"]["environment_identity"] == (
        "cc5359f90df35ec9b60fd56e483911745209519ecce37d69e17c4edd6ea3604f"
    )
    assert set(anchors) == {
        "cpmg-15n-ip",
        "cest-13c-label-cn",
        "2st-binding",
        "dcest-fifu-drd",
    }
    assert all(anchor["status"] == "ELIGIBLE" for anchor in anchors.values())
    assert all(anchor["fresh_occurrence"] for anchor in anchors.values())
    assert all(anchor["shipped_input_bytes_revalidated"] for anchor in anchors.values())
    assert {
        name: (anchor["output_artifact_count"], anchor["artifact_count"])
        for name, anchor in anchors.items()
    } == {
        "cpmg-15n-ip": (340, 341),
        "cest-13c-label-cn": (90, 91),
        "2st-binding": (894, 895),
        "dcest-fifu-drd": (22, 23),
    }

    assert coverage["manifest_identity"] == (
        "c79d318b753eed1e863c6e5e99bbee7f7c1f7cc19ad8c00dc15acf983a6e69c8"
    )
    assert coverage["requirement_count"] == 51
    assert coverage["eligible_requirement_count"] == 10
    assert len(coverage["eligible_requirement_ids"]) == 10
    assert coverage["uncovered_requirement_count"] == 41
    uncovered = set(coverage["uncovered_requirement_ids"])
    assert len(uncovered) == 41
    missing_contracts = coverage["required_evidence_not_yet_acquired"]
    assert {contract["evidence"] for contract in missing_contracts} == set(
        coverage["compiler"]["reason"]
        .removeprefix("Migration-core is missing required evidence: ")
        .split(", ")
    )
    assert (
        set().union(
            *(set(contract["requirement_ids"]) for contract in missing_contracts)
        )
        == uncovered
    )
    assert all(
        contract["supplying_or_qualification_tickets"] for contract in missing_contracts
    )
    assert coverage["compiled_coverage_identity"] is None
    assert coverage["compiler"]["status"] == "FAILED_CLOSED"
    with pytest.raises(MigrationCoreCoverageError):
        migration_core.MigrationCorePhasedStatus.from_bytes(content)
    assert (
        migration_core.migration_core_current_release_selection().status_locator
        != "tests/fixtures/migration_core_canonical_coverage_status_v2.json"
    )
