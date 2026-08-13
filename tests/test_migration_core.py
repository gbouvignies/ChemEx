from __future__ import annotations

import json
from argparse import Namespace
from collections.abc import Mapping
from pathlib import Path

import pytest

import chemex.baselines as baselines
import chemex.migration_core as migration_core
from chemex.baselines import (
    ApprovedAnchorName,
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

    assert migration_core._numerical_probe_identity(baseline) == (
        "NumericalProbeBaseline",
        baseline.identity,
    )
    monkeypatch.setattr(
        migration_core,
        "CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY",
        "0" * 64,
    )
    assert migration_core._numerical_probe_identity(baseline) is None


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
