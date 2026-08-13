from __future__ import annotations

import hashlib
import json
from argparse import Namespace
from collections.abc import Mapping
from pathlib import Path

import pytest

import chemex.baselines as baselines
import chemex.migration_core as migration_core
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

    assert manifest.manifest_version == "migration-core-coverage-v1"
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
    incompatible["manifest_version"] = "migration-core-coverage-v2"
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
    assert coverage["manifest_identity"] == migration_core_coverage_manifest().identity
    assert coverage["requirement_count"] == 51
    assert coverage["eligible_requirement_count"] == 4
    assert coverage["uncovered_requirement_count"] == 47
    assert len(coverage["uncovered_requirement_ids"]) == 47
    assert coverage["compiled_coverage_identity"] is None
    assert coverage["compiler"]["status"] == "FAILED_CLOSED"


def test_current_canonical_coverage_status_records_phased_anchor_eligibility() -> None:
    record = json.loads(
        (
            Path(__file__).parent
            / "fixtures"
            / "migration_core_canonical_coverage_status_v2.json"
        ).read_text(encoding="ascii")
    )
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

    assert coverage["manifest_identity"] == migration_core_coverage_manifest().identity
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
