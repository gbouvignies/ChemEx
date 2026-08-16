"""Fixed #590 migration-core requirement-to-evidence composition."""

from __future__ import annotations

import hashlib
import json
import tarfile
import tempfile
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path, PurePosixPath
from typing import cast

from chemex.baselines import (
    BaselinePublisher,
    CanonicalBaselineValue,
    PublishedEvidence,
)
from chemex.migration_core_lifecycle import (
    FAILURE_REQUIREMENTS,
    LifecycleProbeCapture,
    eligible_failure_requirements,
)
from chemex.migration_core_operational import (
    OPERATIONAL_REQUIREMENTS,
    OperationalReplayCapture,
    eligible_operational_requirements,
)
from chemex.numerical_lanes import (
    LaneAttestation,
    RuntimeEnvironment,
    canonical_lanes,
)
from chemex.optimize.numerical_probes import (
    CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    NumericalProbeBaseline,
)

_SCHEMA_VERSION = 1
_MANIFEST_VERSION = "migration-core-coverage-v2"
_BASELINE_SEMANTIC_VERSION = "chemex-baseline-v1"
_MANIFEST_RESOURCE = "migration_core_coverage_v2.json"
_AUTHORITY_SELECTION_RESOURCE = "migration_core_authority_selection_v1.json"
_CURRENT_RELEASE_SELECTION_RESOURCES = (
    "migration_core_current_release_v3.json",
    "migration_core_current_release_v2.json",
    "migration_core_current_release_v1.json",
)
_CURRENT_STATUS_VERSIONS = {
    "migration-core-canonical-coverage-status-v3": (10, 41),
    "migration-core-canonical-coverage-status-v4": (17, 34),
    "migration-core-canonical-coverage-status-v5": (30, 21),
}
_ANCHOR_RELEASE_VERSION = "migration-core-canonical-anchor-release-v1"
_CURRENT_RELEASE_VERSIONS = {
    "migration-core-current-release-v1",
    "migration-core-current-release-v2",
    "migration-core-current-release-v3",
}
_AUTHORITY_SELECTION_VERSION = "migration-core-authority-selection-v1"
_CANONICAL_ATTESTATION_IDENTITY = (
    "3b3e0bc184826d61ec6652194486c907a5faa4c64b68ec03bbe60b63c660d687"
)
_CANONICAL_ENVIRONMENT_IDENTITY = (
    "cc5359f90df35ec9b60fd56e483911745209519ecce37d69e17c4edd6ea3604f"
)
_CANONICAL_NUMERICAL_PROBE_BASELINE_IDENTITY = (
    "d11f9caa404b8ce3fa96e041659939446e205aae9376b08d6a137ed40cf0bdb0"
)
_UNCERTAINTY_CALIBRATION_FILE_HASH = (
    "cb149ef426d4c32bfaa2b84b5c26170ada1a3fb1b66b5221a447c9fd4a178c5c"
)
_RESAMPLING_CALIBRATION_FILE_HASH = (
    "f1ccbe45222351a2f9fd138d4250dfba80c8d84f3c658de2d33a2acd0bd5622b"
)
_OPERATIONAL_CAPTURE_IDENTITY = (
    "8a0394e50bef697dce880ee35bf7d77d2f8ffa3c8fc8e81ab0a6afbb4d15a59a"
)
_OPERATIONAL_SOURCE_COMMIT = "fb156f86f431a90c65d1a7285bdb6532ab2c51ec"
_OPERATIONAL_LOCKFILE_HASH = (
    "cc7a8e08d8fb8f1ea4255b63452598f6dbe041a8b4024de0f3af065020088004"
)
_APPROVED_EVIDENCE = (
    "anchor:2st-binding",
    "anchor:cest-13c-label-cn",
    "anchor:cpmg-15n-ip",
    "anchor:dcest-fifu-drd",
)
_APPROVED_SUPPORTING_EVIDENCE = (
    "execution:covariance-constrained-uncertainty",
    "execution:de-trf-search",
    "execution:fit-component-aggregation",
    "execution:grid-decomposed",
    "execution:grid-single-component",
    "execution:mcmc",
    "execution:resampling",
    "probe:lifecycle-publication-safety",
    "probe:numerical-optimizer-evaluator",
    "probe:serialization-multiprocessing-cache-replay",
    "workload:qualified-performance",
)
_AUDITED_REQUIREMENTS = (
    "migration-core.anchor.2st-binding.canonical-completion",
    "migration-core.anchor.cest-13c-label-cn.canonical-completion",
    "migration-core.anchor.cpmg-15n-ip.canonical-completion",
    "migration-core.anchor.dcest-fifu-drd.canonical-completion",
    "migration-core.anchor-set.canonical-environment-authority",
    "migration-core.anchor-set.closed-artifact-integrity",
    "migration-core.anchor-set.exact-lineage-identities",
    "migration-core.anchor-set.shipped-workflow-identities",
    "migration-core.direct-trf.routine-behavior",
    "migration-core.direct-trf.legacy-outlier-behavior",
    "migration-core.grid-trf.coupled-execution",
    "migration-core.grid-trf.immutable-27-seed-physical-grid",
    "migration-core.grid-trf.decomposed-execution",
    "migration-core.de-trf.reduced-poor-start-dcest",
    "migration-core.de-trf.deterministic-compact-roots",
    "migration-core.grouped.exact-component-execution",
    "migration-core.grouped.fresh-aggregate-materialization",
    "migration-core.grouped.aggregate-validation",
    "migration-core.covariance.both-scaling-policies",
    "migration-core.covariance.finite-difference-reliability",
    "migration-core.covariance.rank",
    "migration-core.covariance.conditioning",
    "migration-core.covariance.boundary",
    "migration-core.covariance.correlation",
    "migration-core.covariance.analytic-normalization",
    "migration-core.covariance.constrained-propagation",
    "migration-core.covariance.propagation-degeneracy",
    "migration-core.resampling.mc-truth-probe",
    "migration-core.resampling.bootstrap-truth-probe",
    "migration-core.resampling.nucleus-bootstrap-truth-probe",
    "migration-core.resampling.serial-two-worker-replay",
    "migration-core.mcmc.bounded-analytic-probe",
    "migration-core.mcmc.realistic-probe",
    "migration-core.mcmc.complete-chain-replay",
    "migration-core.mcmc.truth-backed-statistical-comparison",
    "migration-core.failure.construction",
    "migration-core.failure.execution",
    "migration-core.failure.materialization",
    "migration-core.failure.commit",
    "migration-core.failure.partial-stochastic-evidence",
    "migration-core.failure.workflow-stop",
    "migration-core.failure.publication",
    "migration-core.operational.serialization",
    "migration-core.operational.multiprocessing",
    "migration-core.operational.cache",
    "migration-core.operational.deterministic-replay",
    "migration-core.performance.fixed-work-profile",
    "migration-core.performance.fixed-work-grid",
    "migration-core.performance.fixed-work-serialization",
    "migration-core.performance.fixed-work-resampling",
    "migration-core.performance.fixed-work-mcmc",
)

_QUALIFIED_UNCERTAINTY_REQUIREMENTS = frozenset(
    {
        "migration-core.covariance.analytic-normalization",
        "migration-core.covariance.both-scaling-policies",
        "migration-core.covariance.conditioning",
        "migration-core.covariance.correlation",
        "migration-core.covariance.propagation-degeneracy",
        "migration-core.covariance.rank",
    }
)
_QUALIFIED_RESAMPLING_REQUIREMENTS = frozenset(
    {
        "migration-core.resampling.bootstrap-truth-probe",
        "migration-core.resampling.mc-truth-probe",
        "migration-core.resampling.serial-two-worker-replay",
    }
)
_QUALIFIED_NUMERICAL_PROBE_REQUIREMENTS = frozenset(
    {
        "migration-core.covariance.finite-difference-reliability",
        "migration-core.direct-trf.legacy-outlier-behavior",
        "migration-core.direct-trf.routine-behavior",
    }
)
_EXPECTED_UNQUALIFIED_REQUIREMENTS = frozenset(
    {
        "migration-core.covariance.boundary",
        "migration-core.covariance.constrained-propagation",
        "migration-core.covariance.finite-difference-reliability",
        "migration-core.resampling.nucleus-bootstrap-truth-probe",
    }
)


class MigrationCoreCoverageError(ValueError):
    """The fixed migration-core coverage cannot be compiled authoritatively."""


def _duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise MigrationCoreCoverageError(
                "Migration-core manifest contains duplicate keys"
            )
        result[key] = value
    return result


def _constant(value: str) -> object:
    raise MigrationCoreCoverageError(f"Migration-core manifest contains {value}")


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("ascii")


def _identity(kind: str, value: object) -> str:
    return hashlib.sha256(
        _canonical_bytes(
            {
                "kind": kind,
                "schema_version": _SCHEMA_VERSION,
                "value": value,
            }
        )
    ).hexdigest()


def _content_identity(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _json_record_from_bytes(content: bytes, name: str) -> Mapping[str, object]:
    try:
        decoded = json.loads(
            content.decode("utf-8"),
            object_pairs_hook=_duplicates,
            parse_constant=_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise MigrationCoreCoverageError(
            f"Migration-core {name} is malformed"
        ) from error
    return _record(decoded, name)


def _canonical_provenance_matches(record: Mapping[str, object]) -> bool:
    try:
        lane = _record(record.get("numerical_lane"), "numerical lane")
        attestation = _record(record.get("lane_attestation"), "lane attestation")
        environment = _record(record.get("runtime_environment"), "runtime environment")
    except MigrationCoreCoverageError:
        return False
    authority = migration_core_authority_selection()
    return (
        lane.get("identity") == authority.lane_identity
        and lane.get("role") == authority.lane_role
        and attestation.get("identity") == authority.attestation_identity
        and attestation.get("environment_identity") == authority.environment_identity
        and attestation.get("workers") == authority.workers
        and attestation.get("native_threads") == authority.native_threads
        and environment.get("identity") == authority.environment_identity
    )


def _qualified_uncertainty_requirements(content: bytes) -> frozenset[str]:
    if hashlib.sha256(content).hexdigest() != _UNCERTAINTY_CALIBRATION_FILE_HASH:
        return frozenset()
    try:
        record = _json_record_from_bytes(content, "uncertainty calibration")
        provenance = _record(
            record.get("canonical_provenance"), "uncertainty provenance"
        )
        phases = _record(
            record.get("selected_phase_policies"), "uncertainty selections"
        )
        metrics = _record(record.get("decisive_metrics"), "uncertainty metrics")
        composed = record.get("composed")
        holdouts = record.get("holdouts")
    except MigrationCoreCoverageError:
        return frozenset()
    if not isinstance(composed, list) or not isinstance(holdouts, list):
        return frozenset()
    composed_by_case = {
        item.get("case"): item
        for value in composed
        if isinstance(value, Mapping)
        for item in (cast("Mapping[str, object]", value),)
    }
    selected_metrics = (
        "scale",
        "finite_difference",
        "driver",
        "rank",
        "weak",
        "cluster",
        "conditioning",
        "correlation",
    )
    try:
        all_selected = all(
            _record(metrics.get(name), f"uncertainty {name}").get("status")
            == "SELECTED"
            for name in selected_metrics
        )
    except MigrationCoreCoverageError:
        return frozenset()
    if (
        record.get("specification_id") != "chemex-uncertainty-calibration-v3"
        or record.get("source_digest")
        != "8015729a614ad1128a8d30a9dbba5557d86c6a3864feeaa5514324397182d752"
        or record.get("specification_digest")
        != "4cae4629318008d1f8fff41657bd4cb01784d03247b700025c57f89dee76fbf9"
        or record.get("policy_digest")
        != "d5a42efd796d4cb5afc31af069157b9fa0d700f61534b94226bdcd4ccd44840d"
        or record.get("status") != "COMPOSED_VALIDATION_FAILED"
        or not _canonical_provenance_matches(provenance)
        or not all_selected
        or phases.get("finite_difference") != [0.0, 2.0, 2, 0]
        or phases.get("rank") != [0.0, 2.0**-36]
        or phases.get("conditioning") != 2.0**20
        or phases.get("correlation") != 64.0
        or tuple(composed_by_case)
        != (
            "A1",
            "A2",
            "F1-absolute",
            "F1-estimated",
            "F2",
        )
        or not all(
            composed_by_case[name].get("passed") is True
            for name in ("A1", "A2", "F1-absolute", "F1-estimated")
        )
        or composed_by_case["F2"].get("passed") is not False
        or composed_by_case["F2"].get("partial_error") != 4.733702319015265e-11
        or not all(
            isinstance(item, Mapping) and item.get("status") == "NOT_RUN"
            for item in holdouts
        )
        or record.get("compatibility")
        != {"reason": "composed_validation_failed", "status": "NOT_RUN"}
    ):
        return frozenset()
    return _QUALIFIED_UNCERTAINTY_REQUIREMENTS


def _qualified_resampling_requirements(content: bytes) -> frozenset[str]:
    if hashlib.sha256(content).hexdigest() != _RESAMPLING_CALIBRATION_FILE_HASH:
        return frozenset()
    try:
        record = _json_record_from_bytes(content, "resampling calibration")
        lane = _record(record.get("canonical_lane"), "resampling lane")
        source = _record(record.get("source"), "resampling source")
        families = _record(record.get("families"), "resampling families")
        mc = _record(families.get("mc"), "MC calibration")
        bs = _record(families.get("bs"), "BS calibration")
        bsn = _record(families.get("bsn"), "BSN calibration")
    except MigrationCoreCoverageError:
        return frozenset()

    def qualified_family(family: Mapping[str, object]) -> bool:
        selection = family.get("selection")
        holdout = family.get("holdout")
        replay = family.get("replay")
        return (
            isinstance(selection, Mapping)
            and selection.get("status") == "SELECTED"
            and selection.get("selected_count") == 128
            and isinstance(holdout, Mapping)
            and holdout.get("status") == "PASS"
            and isinstance(replay, Mapping)
            and replay.get("status") == "PASS"
            and replay.get("serial_workers") == 1
            and replay.get("parallel_workers") == 2
            and replay.get("exact_identity_match") is True
        )

    bsn_selection = bsn.get("selection")
    if (
        record.get("schema_version") != 1
        or record.get("record_version") != "canonical-resampling-calibration-v2"
        or record.get("identity")
        != "3adb78c1c9910d16bd5201f9ef0984bf8974417632fbae2c3af5518f78fa70c4"
        or source.get("dependency_lock_sha256")
        != "cc7a8e08d8fb8f1ea4255b63452598f6dbe041a8b4024de0f3af065020088004"
        or source.get("qualification_script_sha256")
        != "6b6c6709dba53316da120992f199d9701e4a7f6e3e2e3c674e88e548b2e1dcec"
        or source.get("specification_commit")
        != "65517c1242a387635bdaeaa96b05f547e0a72fcd"
        or not _canonical_provenance_matches(lane)
        or not qualified_family(mc)
        or not qualified_family(bs)
        or not isinstance(bsn_selection, Mapping)
        or bsn_selection.get("status") != "NO_ADEQUATE_CANDIDATE"
        or bsn_selection.get("selected_count") is not None
        or bsn.get("holdout") is not None
        or bsn.get("replay") is not None
    ):
        return frozenset()
    return _QUALIFIED_RESAMPLING_REQUIREMENTS


def qualified_migration_core_requirements(
    evidence: str, content: bytes
) -> frozenset[str]:
    """Return only claims qualified by an exact immutable evidence payload."""

    validators = {
        "execution:covariance-constrained-uncertainty": (
            _qualified_uncertainty_requirements
        ),
        "execution:resampling": _qualified_resampling_requirements,
    }
    validator = validators.get(evidence)
    return frozenset() if validator is None else validator(content)


def _record(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or not all(isinstance(key, str) for key in value):
        raise MigrationCoreCoverageError(f"Migration-core {name} must be a record")
    return cast("Mapping[str, object]", value)


def _exact_keys(record: Mapping[str, object], expected: set[str], name: str) -> None:
    if set(record) != expected:
        raise MigrationCoreCoverageError(
            f"Migration-core {name} has unknown or missing fields"
        )


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must be a non-empty string"
        )
    return value


def _digest(value: object, name: str) -> str:
    result = _text(value, name)
    if len(result) != 64 or any(
        character not in "0123456789abcdef" for character in result
    ):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must be a SHA-256 digest"
        )
    return result


def _git_commit(value: object, name: str) -> str:
    result = _text(value, name)
    if len(result) != 40 or any(
        character not in "0123456789abcdef" for character in result
    ):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must be a full Git commit identity"
        )
    return result


def _ordered_texts(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise MigrationCoreCoverageError(f"Migration-core {name} must be a list")
    result = tuple(_text(item, name) for item in value)
    if not result or result != tuple(sorted(result)) or len(set(result)) != len(result):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must contain ordered unique values"
        )
    return result


def _relative_locator(value: object, name: str) -> str:
    locator = _text(value, name)
    path = PurePosixPath(locator)
    if (
        path.is_absolute()
        or path == PurePosixPath(".")
        or any(part in {"", ".", ".."} for part in path.parts)
        or "\\" in locator
    ):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must be a safe repository-relative locator"
        )
    return locator


def _repository_path(root: Path, locator: str, name: str) -> Path:
    trusted_root = root.resolve()
    candidate = (trusted_root / _relative_locator(locator, name)).resolve()
    if not candidate.is_relative_to(trusted_root):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} escapes the trusted repository root"
        )
    return candidate


@dataclass(frozen=True, slots=True)
class MigrationCoreAuthoritySelection:
    """The one repository-frozen #588 authority selected by this release."""

    schema_version: int
    selection_version: str
    lane_identity: str
    lane_role: str
    attestation_identity: str
    environment_identity: str
    image_digest: str
    workers: int
    native_threads: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        lane = canonical_lanes()[0]
        if (
            self.schema_version != _SCHEMA_VERSION
            or self.selection_version != _AUTHORITY_SELECTION_VERSION
            or self.lane_identity != lane.identity
            or self.lane_role != lane.role
            or self.attestation_identity != _CANONICAL_ATTESTATION_IDENTITY
            or self.environment_identity != _CANONICAL_ENVIRONMENT_IDENTITY
            or self.image_digest != lane.semantics.image_digest
            or self.workers != lane.semantics.workers
            or self.native_threads != lane.semantics.native_threads
        ):
            raise MigrationCoreCoverageError(
                "Migration-core authority selection is not the exact frozen lane"
            )
        object.__setattr__(self, "identity", _content_identity(self.to_record()))

    def to_record(self) -> dict[str, object]:
        return {
            "attestation_identity": self.attestation_identity,
            "environment_identity": self.environment_identity,
            "image_digest": self.image_digest,
            "lane_identity": self.lane_identity,
            "lane_role": self.lane_role,
            "native_threads": self.native_threads,
            "schema_version": self.schema_version,
            "selection_version": self.selection_version,
            "workers": self.workers,
        }

    @classmethod
    def from_bytes(cls, content: bytes) -> MigrationCoreAuthoritySelection:
        record = _json_record_from_bytes(content, "authority selection")
        _exact_keys(
            record,
            {
                "attestation_identity",
                "environment_identity",
                "identity",
                "image_digest",
                "lane_identity",
                "lane_role",
                "native_threads",
                "schema_version",
                "selection_version",
                "workers",
            },
            "authority selection",
        )
        schema_version = record.get("schema_version")
        workers = record.get("workers")
        native_threads = record.get("native_threads")
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in (schema_version, workers, native_threads)
        ):
            raise MigrationCoreCoverageError(
                "Migration-core authority selection has malformed integer fields"
            )
        selection = cls(
            cast("int", schema_version),
            _text(record.get("selection_version"), "authority selection version"),
            _digest(record.get("lane_identity"), "selected lane identity"),
            _text(record.get("lane_role"), "selected lane role"),
            _digest(
                record.get("attestation_identity"),
                "selected attestation identity",
            ),
            _digest(
                record.get("environment_identity"),
                "selected environment identity",
            ),
            _text(record.get("image_digest"), "selected image digest"),
            cast("int", workers),
            cast("int", native_threads),
        )
        if record.get("identity") != selection.identity:
            raise MigrationCoreCoverageError(
                "Migration-core authority-selection identity does not match payload"
            )
        return selection


def migration_core_authority_selection() -> MigrationCoreAuthoritySelection:
    """Load the exact repository-frozen authority; callers cannot replace it."""

    content = files("chemex").joinpath(_AUTHORITY_SELECTION_RESOURCE).read_bytes()
    return MigrationCoreAuthoritySelection.from_bytes(content)


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreReleasePayloadMember:
    """One byte-bearing member of the closed canonical anchor release."""

    path: str
    content_hash: str
    size: int

    def __post_init__(self) -> None:
        _relative_locator(self.path, "release payload path")
        _digest(self.content_hash, "release payload hash")
        if isinstance(self.size, bool) or self.size < 0:
            raise MigrationCoreCoverageError(
                "Migration-core release payload size must be non-negative"
            )

    def to_record(self) -> dict[str, object]:
        return {"content_hash": self.content_hash, "path": self.path, "size": self.size}

    @classmethod
    def from_record(
        cls, record: Mapping[str, object]
    ) -> MigrationCoreReleasePayloadMember:
        _exact_keys(record, {"content_hash", "path", "size"}, "release member")
        size = record.get("size")
        if isinstance(size, bool) or not isinstance(size, int):
            raise MigrationCoreCoverageError(
                "Migration-core release payload size must be an integer"
            )
        return cls(
            _relative_locator(record.get("path"), "release payload path"),
            _digest(record.get("content_hash"), "release payload hash"),
            size,
        )


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreReleaseAnchor:
    """Exact publication and input locators for one approved anchor."""

    evidence: str
    publication_root: str
    input_root: str
    case_identity: str
    specification_identity: str
    occurrence_identity: str
    occurrence_lifecycle_identity: str
    result_bundle_identity: str
    member_manifest_identity: str
    publication_manifest_identity: str

    def __post_init__(self) -> None:
        _text(self.evidence, "release anchor evidence")
        _relative_locator(self.publication_root, "release publication root")
        _relative_locator(self.input_root, "release input root")
        for value, name in (
            (self.case_identity, "release case identity"),
            (self.specification_identity, "release specification identity"),
            (self.occurrence_identity, "release occurrence identity"),
            (self.occurrence_lifecycle_identity, "release lifecycle identity"),
            (self.result_bundle_identity, "release bundle identity"),
            (self.member_manifest_identity, "release member-manifest identity"),
            (self.publication_manifest_identity, "release publication identity"),
        ):
            _digest(value, name)

    def to_record(self) -> dict[str, object]:
        return {
            "case_identity": self.case_identity,
            "evidence": self.evidence,
            "input_root": self.input_root,
            "member_manifest_identity": self.member_manifest_identity,
            "occurrence_identity": self.occurrence_identity,
            "occurrence_lifecycle_identity": self.occurrence_lifecycle_identity,
            "publication_manifest_identity": self.publication_manifest_identity,
            "publication_root": self.publication_root,
            "result_bundle_identity": self.result_bundle_identity,
            "specification_identity": self.specification_identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> MigrationCoreReleaseAnchor:
        expected = {
            "case_identity",
            "evidence",
            "input_root",
            "member_manifest_identity",
            "occurrence_identity",
            "occurrence_lifecycle_identity",
            "publication_manifest_identity",
            "publication_root",
            "result_bundle_identity",
            "specification_identity",
        }
        _exact_keys(record, expected, "release anchor")
        return cls(
            _text(record.get("evidence"), "release anchor evidence"),
            _relative_locator(
                record.get("publication_root"), "release publication root"
            ),
            _relative_locator(record.get("input_root"), "release input root"),
            _digest(record.get("case_identity"), "release case identity"),
            _digest(
                record.get("specification_identity"),
                "release specification identity",
            ),
            _digest(record.get("occurrence_identity"), "release occurrence identity"),
            _digest(
                record.get("occurrence_lifecycle_identity"),
                "release lifecycle identity",
            ),
            _digest(record.get("result_bundle_identity"), "release bundle identity"),
            _digest(
                record.get("member_manifest_identity"),
                "release member-manifest identity",
            ),
            _digest(
                record.get("publication_manifest_identity"),
                "release publication identity",
            ),
        )


@dataclass(frozen=True, slots=True)
class MigrationCoreAnchorRelease:
    """Closed content manifest for the durable four-anchor publication release."""

    schema_version: int
    release_version: str
    repository_commit: str
    publisher_root: str
    source_archive_locator: str
    source_archive_hash: str
    anchors: tuple[MigrationCoreReleaseAnchor, ...]
    payload_members: tuple[MigrationCoreReleasePayloadMember, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            self.schema_version != _SCHEMA_VERSION
            or self.release_version != _ANCHOR_RELEASE_VERSION
        ):
            raise MigrationCoreCoverageError(
                "Migration-core anchor release has incompatible version"
            )
        _git_commit(self.repository_commit, "release repository commit")
        _relative_locator(self.publisher_root, "release publisher root")
        _relative_locator(self.source_archive_locator, "release source locator")
        _digest(self.source_archive_hash, "release source archive hash")
        if (
            tuple(item.evidence for item in self.anchors) != _APPROVED_EVIDENCE
            or tuple(sorted(self.anchors)) != self.anchors
        ):
            raise MigrationCoreCoverageError(
                "Migration-core release must contain the ordered approved anchors"
            )
        if (
            not self.payload_members
            or tuple(sorted(self.payload_members)) != self.payload_members
            or len({item.path for item in self.payload_members})
            != len(self.payload_members)
        ):
            raise MigrationCoreCoverageError(
                "Migration-core release payload manifest is not closed and ordered"
            )
        source = next(
            (
                item
                for item in self.payload_members
                if item.path == self.source_archive_locator
            ),
            None,
        )
        if source is None or source.content_hash != self.source_archive_hash:
            raise MigrationCoreCoverageError(
                "Migration-core release source archive is unresolved"
            )
        object.__setattr__(self, "identity", _content_identity(self.to_record()))

    def to_record(self) -> dict[str, object]:
        return {
            "anchors": [item.to_record() for item in self.anchors],
            "payload_members": [item.to_record() for item in self.payload_members],
            "publisher_root": self.publisher_root,
            "release_version": self.release_version,
            "repository_commit": self.repository_commit,
            "schema_version": self.schema_version,
            "source_archive_hash": self.source_archive_hash,
            "source_archive_locator": self.source_archive_locator,
        }

    @classmethod
    def from_bytes(cls, content: bytes) -> MigrationCoreAnchorRelease:
        record = _json_record_from_bytes(content, "anchor release")
        _exact_keys(
            record,
            {
                "anchors",
                "identity",
                "payload_members",
                "publisher_root",
                "release_version",
                "repository_commit",
                "schema_version",
                "source_archive_hash",
                "source_archive_locator",
            },
            "anchor release",
        )
        raw_anchors = record.get("anchors")
        raw_members = record.get("payload_members")
        schema_version = record.get("schema_version")
        if (
            not isinstance(raw_anchors, list)
            or not isinstance(raw_members, list)
            or isinstance(schema_version, bool)
            or not isinstance(schema_version, int)
        ):
            raise MigrationCoreCoverageError(
                "Migration-core anchor release has malformed fields"
            )
        release = cls(
            schema_version,
            _text(record.get("release_version"), "anchor release version"),
            _git_commit(record.get("repository_commit"), "release repository commit"),
            _relative_locator(record.get("publisher_root"), "release publisher root"),
            _relative_locator(
                record.get("source_archive_locator"), "release source locator"
            ),
            _digest(record.get("source_archive_hash"), "release source hash"),
            tuple(
                MigrationCoreReleaseAnchor.from_record(_record(item, "release anchor"))
                for item in raw_anchors
            ),
            tuple(
                MigrationCoreReleasePayloadMember.from_record(
                    _record(item, "release member")
                )
                for item in raw_members
            ),
        )
        if record.get("identity") != release.identity:
            raise MigrationCoreCoverageError(
                "Migration-core anchor-release identity does not match payload"
            )
        return release


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreAnchorEvidence:
    """One exact approved anchor execution eligible for manifest coverage."""

    evidence: str
    case_name: str
    case_identity: str
    workflow: CanonicalBaselineValue
    artifact_inventory_version: str

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> MigrationCoreAnchorEvidence:
        _exact_keys(
            record,
            {
                "artifact_inventory_version",
                "case_identity",
                "case_name",
                "evidence",
                "workflow",
            },
            "anchor evidence",
        )
        workflow = _record(record.get("workflow"), "anchor workflow")
        return cls(
            _text(record.get("evidence"), "anchor evidence identifier"),
            _text(record.get("case_name"), "anchor case name"),
            _digest(record.get("case_identity"), "anchor case identity"),
            CanonicalBaselineValue.from_record(workflow, "anchor workflow"),
            _text(
                record.get("artifact_inventory_version"),
                "anchor artifact inventory version",
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "artifact_inventory_version": self.artifact_inventory_version,
            "case_identity": self.case_identity,
            "case_name": self.case_name,
            "evidence": self.evidence,
            "workflow": self.workflow.to_record_value(),
        }


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreEvidenceContract:
    """One explicitly eligible supporting execution selected by the audit."""

    evidence: str
    evidence_type: str
    source_issue: int

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> MigrationCoreEvidenceContract:
        _exact_keys(
            record,
            {"evidence", "evidence_type", "source_issue"},
            "supporting evidence",
        )
        source_issue = record.get("source_issue")
        if isinstance(source_issue, bool) or not isinstance(source_issue, int):
            raise MigrationCoreCoverageError(
                "Migration-core supporting evidence issue must be an integer"
            )
        return cls(
            _text(record.get("evidence"), "supporting evidence identifier"),
            _text(record.get("evidence_type"), "supporting evidence type"),
            source_issue,
        )

    def to_record(self) -> dict[str, object]:
        return {
            "evidence": self.evidence,
            "evidence_type": self.evidence_type,
            "source_issue": self.source_issue,
        }


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreRequirement:
    """One required claim and its exact eligible evidence composition."""

    identifier: str
    evidence: tuple[str, ...]

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> MigrationCoreRequirement:
        _exact_keys(record, {"evidence", "identifier"}, "coverage requirement")
        return cls(
            _text(record.get("identifier"), "requirement identifier"),
            _ordered_texts(record.get("evidence"), "eligible evidence"),
        )

    def to_record(self) -> dict[str, object]:
        return {"evidence": list(self.evidence), "identifier": self.identifier}


@dataclass(frozen=True, slots=True)
class MigrationCoreCoverageManifest:
    """The single reviewed requirement-to-evidence manifest for #590."""

    schema_version: int
    manifest_version: str
    baseline_semantic_version: str
    anchors: tuple[MigrationCoreAnchorEvidence, ...]
    supporting_evidence: tuple[MigrationCoreEvidenceContract, ...]
    requirements: tuple[MigrationCoreRequirement, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA_VERSION:
            raise MigrationCoreCoverageError(
                "Migration-core schema has incompatible version"
            )
        if self.manifest_version != _MANIFEST_VERSION:
            raise MigrationCoreCoverageError(
                "Migration-core manifest has incompatible version"
            )
        if self.baseline_semantic_version != _BASELINE_SEMANTIC_VERSION:
            raise MigrationCoreCoverageError(
                "Migration-core baseline semantics have incompatible version"
            )
        if tuple(item.evidence for item in self.anchors) != _APPROVED_EVIDENCE:
            raise MigrationCoreCoverageError(
                "Migration-core manifest must contain exactly the four approved anchors"
            )
        if (
            tuple(item.evidence for item in self.supporting_evidence)
            != _APPROVED_SUPPORTING_EVIDENCE
        ):
            raise MigrationCoreCoverageError(
                "Migration-core manifest has incomplete supporting evidence"
            )
        requirement_ids = tuple(item.identifier for item in self.requirements)
        if requirement_ids != _AUDITED_REQUIREMENTS:
            raise MigrationCoreCoverageError(
                "Migration-core manifest does not contain the fixed audited requirements"
            )
        approved = set(_APPROVED_EVIDENCE) | set(_APPROVED_SUPPORTING_EVIDENCE)
        if any(
            not set(requirement.evidence) <= approved
            for requirement in self.requirements
        ):
            raise MigrationCoreCoverageError(
                "Migration-core requirement names ineligible evidence"
            )
        object.__setattr__(
            self, "identity", _identity("migration-core-manifest", self.to_record())
        )

    @classmethod
    def from_bytes(cls, content: bytes) -> MigrationCoreCoverageManifest:
        try:
            decoded = json.loads(
                content.decode("utf-8"),
                object_pairs_hook=_duplicates,
                parse_constant=_constant,
            )
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise MigrationCoreCoverageError(
                "Migration-core manifest is malformed"
            ) from error
        record = _record(decoded, "manifest")
        _exact_keys(
            record,
            {
                "anchors",
                "baseline_semantic_version",
                "manifest_version",
                "requirements",
                "schema_version",
                "supporting_evidence",
            },
            "manifest",
        )
        raw_anchors = record.get("anchors")
        raw_supporting = record.get("supporting_evidence")
        raw_requirements = record.get("requirements")
        if (
            not isinstance(raw_anchors, list)
            or not isinstance(raw_supporting, list)
            or not isinstance(raw_requirements, list)
        ):
            raise MigrationCoreCoverageError(
                "Migration-core manifest evidence and requirements must be lists"
            )
        schema_version = record.get("schema_version")
        if isinstance(schema_version, bool) or not isinstance(schema_version, int):
            raise MigrationCoreCoverageError(
                "Migration-core schema version must be an integer"
            )
        return cls(
            schema_version,
            _text(record.get("manifest_version"), "manifest version"),
            _text(
                record.get("baseline_semantic_version"),
                "baseline semantic version",
            ),
            tuple(
                MigrationCoreAnchorEvidence.from_record(_record(item, "anchor"))
                for item in raw_anchors
            ),
            tuple(
                MigrationCoreEvidenceContract.from_record(
                    _record(item, "supporting evidence")
                )
                for item in raw_supporting
            ),
            tuple(
                MigrationCoreRequirement.from_record(_record(item, "requirement"))
                for item in raw_requirements
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "anchors": [item.to_record() for item in self.anchors],
            "baseline_semantic_version": self.baseline_semantic_version,
            "manifest_version": self.manifest_version,
            "requirements": [item.to_record() for item in self.requirements],
            "schema_version": self.schema_version,
            "supporting_evidence": [
                item.to_record() for item in self.supporting_evidence
            ],
        }

    def to_bytes(self) -> bytes:
        return _canonical_bytes(self.to_record())


@dataclass(frozen=True, slots=True)
class CompiledMigrationCoreEvidence:
    """Exact publisher-validated lineage selected for one evidence contract."""

    evidence: str
    case_identity: str
    specification_identity: str
    occurrence_identity: str
    lane_attestation_identity: str
    workflow: CanonicalBaselineValue
    result_bundle_identity: str
    artifact_member_identities: tuple[str, ...]
    implementation_authority_identity: str
    publication_manifest_identity: str

    def __post_init__(self) -> None:
        for value, name in (
            (self.case_identity, "compiled case identity"),
            (self.specification_identity, "compiled specification identity"),
            (self.occurrence_identity, "compiled occurrence identity"),
            (self.lane_attestation_identity, "compiled lane attestation identity"),
            (self.result_bundle_identity, "compiled result bundle identity"),
            (
                self.implementation_authority_identity,
                "compiled implementation authority identity",
            ),
            (
                self.publication_manifest_identity,
                "compiled publication manifest identity",
            ),
        ):
            _digest(value, name)
        if not self.artifact_member_identities:
            raise MigrationCoreCoverageError(
                "Compiled evidence must retain artifact identities"
            )
        for identity in self.artifact_member_identities:
            _digest(identity, "compiled artifact member identity")

    def to_record(self) -> dict[str, object]:
        return {
            "artifact_member_identities": list(self.artifact_member_identities),
            "case_identity": self.case_identity,
            "evidence": self.evidence,
            "implementation_authority_identity": (
                self.implementation_authority_identity
            ),
            "lane_attestation_identity": self.lane_attestation_identity,
            "occurrence_identity": self.occurrence_identity,
            "publication_manifest_identity": self.publication_manifest_identity,
            "result_bundle_identity": self.result_bundle_identity,
            "specification_identity": self.specification_identity,
            "workflow": self.workflow.to_record_value(),
        }


@dataclass(frozen=True, slots=True, order=True)
class CompiledMigrationCoreSupportingEvidence:
    """Concrete native evidence accepted by one fixed audit contract."""

    evidence: str
    evidence_type: str
    identity: str

    def __post_init__(self) -> None:
        _text(self.evidence, "compiled supporting evidence")
        _text(self.evidence_type, "compiled supporting evidence type")
        _digest(self.identity, "compiled supporting evidence identity")


@dataclass(frozen=True, slots=True)
class CompiledMigrationCoreCoverage:
    """Publisher-verified evidence identities composed per manifest requirement."""

    manifest_identity: str
    coverage: tuple[tuple[str, tuple[str, ...]], ...]
    evidence: tuple[CompiledMigrationCoreEvidence, ...]
    supporting_evidence: tuple[CompiledMigrationCoreSupportingEvidence, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _digest(self.manifest_identity, "compiled manifest identity")
        if not self.coverage or not self.evidence or not self.supporting_evidence:
            raise MigrationCoreCoverageError(
                "Compiled migration-core coverage is empty"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "compiled-migration-core-coverage",
                (
                    self.manifest_identity,
                    self.coverage,
                    tuple(item.to_record() for item in self.evidence),
                    tuple(
                        (item.evidence, item.evidence_type, item.identity)
                        for item in self.supporting_evidence
                    ),
                ),
            ),
        )

    @property
    def bundle_identities(self) -> tuple[str, ...]:
        return tuple(sorted(item.result_bundle_identity for item in self.evidence))


def migration_core_coverage_manifest() -> MigrationCoreCoverageManifest:
    """Load the one shipped, versioned #590 coverage manifest."""
    content = files("chemex").joinpath(_MANIFEST_RESOURCE).read_bytes()
    return MigrationCoreCoverageManifest.from_bytes(content)


def _file_hash(path: Path, name: str) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise MigrationCoreCoverageError(
            f"Migration-core {name} is missing or unavailable"
        ) from error


def _safe_extract_release(archive: Path, destination: Path) -> None:
    seen: set[str] = set()
    try:
        with tarfile.open(archive, mode="r:xz") as opened:
            for member in opened:
                locator = _relative_locator(member.name, "archive member")
                if locator in seen:
                    raise MigrationCoreCoverageError(
                        "Migration-core archive contains duplicate members"
                    )
                seen.add(locator)
                target = _repository_path(destination, locator, "archive member")
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                if not member.isfile():
                    raise MigrationCoreCoverageError(
                        "Migration-core archive contains an unsafe member type"
                    )
                source = opened.extractfile(member)
                if source is None:
                    raise MigrationCoreCoverageError(
                        "Migration-core archive member is unavailable"
                    )
                content = source.read()
                if len(content) != member.size:
                    raise MigrationCoreCoverageError(
                        "Migration-core archive member is truncated"
                    )
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
    except (EOFError, OSError, tarfile.TarError) as error:
        raise MigrationCoreCoverageError(
            "Migration-core anchor archive is corrupt or truncated"
        ) from error


@dataclass(frozen=True, slots=True)
class ResolvedMigrationCoreAnchorRelease:
    """One safely extracted release with recursively validated publications."""

    release: MigrationCoreAnchorRelease
    publisher: BaselinePublisher
    published: Mapping[str, PublishedEvidence]
    extraction_root: Path


def _validate_release_anchor(
    anchor: MigrationCoreReleaseAnchor,
    published: PublishedEvidence,
    extraction_root: Path,
) -> None:
    actual = {
        "case_identity": published.case.identity,
        "specification_identity": published.specification.identity,
        "occurrence_identity": published.occurrence.identity,
        "occurrence_lifecycle_identity": published.occurrence.lifecycle_identity,
        "result_bundle_identity": published.bundle.identity,
        "member_manifest_identity": published.bundle.manifest_identity,
        "publication_manifest_identity": published.manifest_identity,
    }
    expected = {
        "case_identity": anchor.case_identity,
        "specification_identity": anchor.specification_identity,
        "occurrence_identity": anchor.occurrence_identity,
        "occurrence_lifecycle_identity": anchor.occurrence_lifecycle_identity,
        "result_bundle_identity": anchor.result_bundle_identity,
        "member_manifest_identity": anchor.member_manifest_identity,
        "publication_manifest_identity": anchor.publication_manifest_identity,
    }
    if actual != expected:
        raise MigrationCoreCoverageError(
            f"Evidence {anchor.evidence} differs from its frozen release lineage"
        )
    expected_location = _repository_path(
        extraction_root, anchor.publication_root, "release publication root"
    )
    if published.location.resolve() != expected_location:
        raise MigrationCoreCoverageError(
            f"Evidence {anchor.evidence} has an incompatible publication locator"
        )
    input_root = _repository_path(
        extraction_root, anchor.input_root, "release input root"
    )
    expected_inputs = {member.content_hash: member for member in published.case.inputs}
    try:
        actual_inputs = {path.name for path in input_root.iterdir() if path.is_file()}
    except OSError as error:
        raise MigrationCoreCoverageError(
            f"Evidence {anchor.evidence} has unavailable captured inputs"
        ) from error
    if actual_inputs != set(expected_inputs):
        raise MigrationCoreCoverageError(
            f"Evidence {anchor.evidence} has missing or unexpected captured inputs"
        )
    for content_hash, member in expected_inputs.items():
        content_path = input_root / content_hash
        try:
            content = content_path.read_bytes()
        except OSError as error:
            raise MigrationCoreCoverageError(
                f"Evidence {anchor.evidence} has an unavailable captured input"
            ) from error
        if (
            hashlib.sha256(content).hexdigest() != content_hash
            or len(content) != member.size
        ):
            raise MigrationCoreCoverageError(
                f"Evidence {anchor.evidence} captured input hash mismatch"
            )


@contextmanager
def resolve_migration_core_anchor_release(
    archive: Path,
    *,
    expected_archive_hash: str,
    expected_release_identity: str,
) -> Iterator[ResolvedMigrationCoreAnchorRelease]:
    """Resolve the one closed release without granting it scientific authority."""

    _digest(expected_archive_hash, "expected anchor archive hash")
    _digest(expected_release_identity, "expected anchor release identity")
    if _file_hash(archive, "anchor archive") != expected_archive_hash:
        raise MigrationCoreCoverageError(
            "Migration-core anchor archive hash does not match selection"
        )
    with tempfile.TemporaryDirectory(prefix="chemex-migration-core-release-") as raw:
        extraction_root = Path(raw)
        _safe_extract_release(archive, extraction_root)
        try:
            release_content = (extraction_root / "release.json").read_bytes()
        except OSError as error:
            raise MigrationCoreCoverageError(
                "Migration-core anchor release record is missing"
            ) from error
        release = MigrationCoreAnchorRelease.from_bytes(release_content)
        if release.identity != expected_release_identity:
            raise MigrationCoreCoverageError(
                "Migration-core anchor release identity does not match selection"
            )
        extracted_files = {
            path.relative_to(extraction_root).as_posix()
            for path in extraction_root.rglob("*")
            if path.is_file()
        }
        expected_files = {"release.json"} | {
            item.path for item in release.payload_members
        }
        if extracted_files != expected_files:
            raise MigrationCoreCoverageError(
                "Migration-core anchor release payload is not closed"
            )
        for member in release.payload_members:
            path = _repository_path(
                extraction_root, member.path, "release payload member"
            )
            try:
                content = path.read_bytes()
            except OSError as error:
                raise MigrationCoreCoverageError(
                    "Migration-core anchor release member is missing"
                ) from error
            if (
                hashlib.sha256(content).hexdigest() != member.content_hash
                or len(content) != member.size
            ):
                raise MigrationCoreCoverageError(
                    "Migration-core anchor release payload hash mismatch"
                )
        publisher_root = _repository_path(
            extraction_root, release.publisher_root, "release publisher root"
        )
        publisher = BaselinePublisher(publisher_root)
        published: dict[str, PublishedEvidence] = {}
        for anchor in release.anchors:
            try:
                evidence = publisher.read(anchor.result_bundle_identity)
            except (TypeError, ValueError) as error:
                raise MigrationCoreCoverageError(
                    f"Evidence {anchor.evidence} cannot be reconstructed"
                ) from error
            _validate_release_anchor(anchor, evidence, extraction_root)
            published[anchor.evidence] = evidence
        yield ResolvedMigrationCoreAnchorRelease(
            release, publisher, published, extraction_root
        )


def _validate_anchor_evidence(
    expected: MigrationCoreAnchorEvidence, published: PublishedEvidence
) -> None:
    authority = migration_core_authority_selection()
    canonical_lane = canonical_lanes()[0]
    workflow = published.specification.workflow
    inventory = published.specification.artifact_inventory.to_record_value()
    settings = published.specification.execution_settings.to_record_value()
    if not isinstance(inventory, Mapping) or not isinstance(settings, Mapping):
        raise MigrationCoreCoverageError("Anchor specification records are malformed")
    if (
        published.case.name != expected.case_name
        or published.case.identity != expected.case_identity
        or workflow != expected.workflow
        or inventory.get("version") != expected.artifact_inventory_version
        or published.specification.roles != ("Scientific anchor",)
        or published.specification.claims != ("legacy-observation-continuity",)
        or published.specification.lane_reference != authority.lane_identity
        or published.occurrence.lane_reference != authority.lane_identity
        or settings.get("workers") != authority.workers
        or settings.get("native_threads") != authority.native_threads
        or published.occurrence.lifecycle != "SUCCEEDED"
        or published.bundle.implementation.authority_role
        != "LegacyObservationImplementation"
    ):
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} is incompatible with the coverage manifest"
        )
    if published.occurrence.lane_attestation_identity != authority.attestation_identity:
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} has an exact authority mismatch"
        )
    attestation_members = tuple(
        member
        for member in published.bundle.members
        if member.role == "environment:lane-attestation.json"
    )
    environment_members = tuple(
        member
        for member in published.bundle.members
        if member.role == "environment:runtime-environment.json"
    )
    if len(attestation_members) != 1 or len(environment_members) != 1:
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} lacks exact environment evidence"
        )
    try:
        attestation_content = (
            published.location / "members" / attestation_members[0].content_hash
        ).read_bytes()
        attestation = LaneAttestation.from_record(
            _json_record_from_bytes(attestation_content, "lane attestation")
        )
        environment_content = (
            published.location / "members" / environment_members[0].content_hash
        ).read_bytes()
        environment = RuntimeEnvironment.from_record(
            _json_record_from_bytes(environment_content, "runtime environment")
        )
    except (MigrationCoreCoverageError, OSError, TypeError, ValueError) as error:
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} has unavailable environment evidence"
        ) from error
    if (
        attestation.identity != authority.attestation_identity
        or attestation.identity != published.occurrence.lane_attestation_identity
        or attestation.lane_identity != authority.lane_identity
        or attestation.environment_identity != authority.environment_identity
        or attestation.workers != authority.workers
        or attestation.native_threads != authority.native_threads
        or environment.identity != authority.environment_identity
        or environment.semantics != canonical_lane.semantics
        or environment.semantics.image_digest != authority.image_digest
        or environment.semantics.workers != authority.workers
        or environment.semantics.native_threads != authority.native_threads
    ):
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} has an exact authority mismatch"
        )


def _validated_numerical_probe_identity(
    evidence: object,
) -> tuple[str, str] | None:
    authority = migration_core_authority_selection()
    if (
        isinstance(evidence, NumericalProbeBaseline)
        and evidence.historical_qualification == "REFERENCE_MATCHED"
        and evidence.manifest_identity == CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
        and evidence.reference_manifest_identity
        == CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
        and evidence.observed_lane_role == authority.lane_role
        and evidence.observed_lane_identity == authority.lane_identity
        and evidence.observed_attestation_identity == authority.attestation_identity
        and evidence.observed_environment_identity == authority.environment_identity
        and evidence.identity == _CANONICAL_NUMERICAL_PROBE_BASELINE_IDENTITY
    ):
        return type(evidence).__name__, evidence.identity
    return None


def _validated_lifecycle_publication_safety_identity(
    evidence: object,
) -> tuple[str, str] | None:
    current = migration_core_current_release_selection()
    authority = migration_core_authority_selection()
    if (
        isinstance(evidence, LifecycleProbeCapture)
        and current.lifecycle_probe_identity is not None
        and evidence.identity == current.lifecycle_probe_identity
        and current.lifecycle_source_commit is not None
        and current.lifecycle_lockfile_hash is not None
        and eligible_failure_requirements(
            evidence,
            source_commit=current.lifecycle_source_commit,
            lockfile_hash=current.lifecycle_lockfile_hash,
            lane_identity=authority.lane_identity,
            attestation_identity=authority.attestation_identity,
            environment_identity=authority.environment_identity,
        )
        == FAILURE_REQUIREMENTS
    ):
        return "LifecyclePublicationSafetyEvidence", evidence.identity
    return None


def _validated_uncertainty_calibration_identity(
    evidence: object,
) -> tuple[str, str] | None:
    if isinstance(evidence, bytes) and qualified_migration_core_requirements(
        "execution:covariance-constrained-uncertainty", evidence
    ):
        return "UncertaintyCalibrationEvidence", hashlib.sha256(evidence).hexdigest()
    return None


def _validated_resampling_calibration_identity(
    evidence: object,
) -> tuple[str, str] | None:
    if isinstance(evidence, bytes) and qualified_migration_core_requirements(
        "execution:resampling", evidence
    ):
        record = _json_record_from_bytes(evidence, "resampling calibration")
        return "ResamplingCalibrationEvidence", _digest(
            record.get("identity"), "resampling evidence identity"
        )
    return None


def _validated_operational_replay_identity(
    evidence: object,
) -> tuple[str, str] | None:
    authority = migration_core_authority_selection()
    if (
        isinstance(evidence, OperationalReplayCapture)
        and evidence.identity == _OPERATIONAL_CAPTURE_IDENTITY
        and eligible_operational_requirements(
            evidence,
            source_commit=_OPERATIONAL_SOURCE_COMMIT,
            lockfile_hash=_OPERATIONAL_LOCKFILE_HASH,
            lane_identity=authority.lane_identity,
            attestation_identity=authority.attestation_identity,
            environment_identity=authority.environment_identity,
        )
        == OPERATIONAL_REQUIREMENTS
    ):
        return "SerializationMultiprocessingCacheReplayEvidence", evidence.identity
    return None


def _unavailable_identity(_evidence: object) -> tuple[str, str] | None:
    return None


_SUPPORTING_VALIDATORS: dict[str, Callable[[object], tuple[str, str] | None]] = {
    # Exact reduced-case and policy fingerprints are empirical release inputs.
    # Broad successful native objects deliberately cannot mint these selections.
    "execution:covariance-constrained-uncertainty": (
        _validated_uncertainty_calibration_identity
    ),
    "execution:de-trf-search": _unavailable_identity,
    "execution:fit-component-aggregation": _unavailable_identity,
    "execution:grid-decomposed": _unavailable_identity,
    "execution:grid-single-component": _unavailable_identity,
    "execution:mcmc": _unavailable_identity,
    "execution:resampling": _validated_resampling_calibration_identity,
    "probe:lifecycle-publication-safety": (
        _validated_lifecycle_publication_safety_identity
    ),
    "probe:numerical-optimizer-evaluator": _validated_numerical_probe_identity,
    "probe:serialization-multiprocessing-cache-replay": (
        _validated_operational_replay_identity
    ),
    "workload:qualified-performance": _unavailable_identity,
}


def _supporting_qualified_requirements(
    evidence_name: str, evidence: object
) -> frozenset[str]:
    if evidence_name == "probe:numerical-optimizer-evaluator":
        return (
            _QUALIFIED_NUMERICAL_PROBE_REQUIREMENTS
            if _validated_numerical_probe_identity(evidence) is not None
            else frozenset()
        )
    if evidence_name == "probe:lifecycle-publication-safety":
        return (
            frozenset(FAILURE_REQUIREMENTS)
            if _validated_lifecycle_publication_safety_identity(evidence) is not None
            else frozenset()
        )
    if evidence_name in {
        "execution:covariance-constrained-uncertainty",
        "execution:resampling",
    } and isinstance(evidence, bytes):
        return qualified_migration_core_requirements(evidence_name, evidence)
    if (
        evidence_name == "probe:serialization-multiprocessing-cache-replay"
        and isinstance(evidence, OperationalReplayCapture)
    ):
        authority = migration_core_authority_selection()
        return eligible_operational_requirements(
            evidence,
            source_commit=_OPERATIONAL_SOURCE_COMMIT,
            lockfile_hash=_OPERATIONAL_LOCKFILE_HASH,
            lane_identity=authority.lane_identity,
            attestation_identity=authority.attestation_identity,
            environment_identity=authority.environment_identity,
        )
    return frozenset()


def _supporting_native_identity(
    evidence_name: str, evidence: object
) -> tuple[str, str]:
    validated = _SUPPORTING_VALIDATORS[evidence_name](evidence)
    if validated is None:
        raise MigrationCoreCoverageError(
            f"Evidence {evidence_name} is unavailable or not eligible"
        )
    evidence_type, identity = validated
    return evidence_type, _digest(identity, "supporting native evidence identity")


def _compiled_evidence(
    evidence: str, published: PublishedEvidence
) -> CompiledMigrationCoreEvidence:
    lane_attestation_identity = published.occurrence.lane_attestation_identity
    if lane_attestation_identity is None:
        raise MigrationCoreCoverageError(
            f"Evidence {evidence} lacks environment identity"
        )
    return CompiledMigrationCoreEvidence(
        evidence,
        published.case.identity,
        published.specification.identity,
        published.occurrence.identity,
        lane_attestation_identity,
        published.specification.workflow,
        published.bundle.identity,
        tuple(member.identity for member in published.bundle.members),
        published.bundle.implementation.identity,
        published.manifest_identity,
    )


def _qualified_coverage_items(
    manifest: MigrationCoreCoverageManifest,
    evidence_identities: Mapping[str, str],
    qualified_requirements: Mapping[str, frozenset[str]],
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    items: list[tuple[str, tuple[str, ...]]] = []
    for requirement in manifest.requirements:
        unqualified = tuple(
            evidence
            for evidence in requirement.evidence
            if evidence in qualified_requirements
            and requirement.identifier not in qualified_requirements[evidence]
        )
        if unqualified:
            raise MigrationCoreCoverageError(
                f"Evidence {', '.join(unqualified)} does not qualify requirement "
                f"{requirement.identifier}"
            )
        items.append(
            (
                requirement.identifier,
                tuple(
                    evidence_identities[evidence] for evidence in requirement.evidence
                ),
            )
        )
    return tuple(items)


def compile_migration_core_coverage(
    manifest: MigrationCoreCoverageManifest,
    publisher: BaselinePublisher,
    selections: Mapping[str, str],
    supporting_evidence: Mapping[str, object],
) -> CompiledMigrationCoreCoverage:
    """Compile #590 coverage, failing closed before any claim is marked covered."""
    if not isinstance(manifest, MigrationCoreCoverageManifest):
        raise TypeError("Migration-core compilation requires its versioned manifest")
    selected_keys = set(selections)
    supporting_keys = set(supporting_evidence)
    missing = sorted(
        (set(_APPROVED_EVIDENCE) - selected_keys)
        | (set(_APPROVED_SUPPORTING_EVIDENCE) - supporting_keys)
    )
    unknown = sorted(
        (selected_keys - set(_APPROVED_EVIDENCE))
        | (supporting_keys - set(_APPROVED_SUPPORTING_EVIDENCE))
    )
    if missing:
        raise MigrationCoreCoverageError(
            "Migration-core is missing required evidence: " + ", ".join(missing)
        )
    if unknown:
        raise MigrationCoreCoverageError(
            "Migration-core contains unapproved evidence: " + ", ".join(unknown)
        )

    resolved: dict[str, PublishedEvidence] = {}
    for expected in manifest.anchors:
        identity = selections.get(expected.evidence)
        try:
            bundle_identity = _digest(identity, "selected bundle identity")
            published = publisher.read(bundle_identity)
        except (MigrationCoreCoverageError, TypeError, ValueError) as error:
            raise MigrationCoreCoverageError(
                f"Evidence {expected.evidence} is unresolved or unavailable"
            ) from error
        _validate_anchor_evidence(expected, published)
        resolved[expected.evidence] = published

    native_items: list[CompiledMigrationCoreSupportingEvidence] = []
    for expected in manifest.supporting_evidence:
        evidence_type, identity = _supporting_native_identity(
            expected.evidence, supporting_evidence[expected.evidence]
        )
        if evidence_type != expected.evidence_type:
            raise MigrationCoreCoverageError(
                f"Evidence {expected.evidence} has an incompatible native type"
            )
        native_items.append(
            CompiledMigrationCoreSupportingEvidence(
                expected.evidence, evidence_type, identity
            )
        )
    native = tuple(native_items)
    qualified_requirements = {
        name: _supporting_qualified_requirements(name, evidence)
        for name, evidence in supporting_evidence.items()
    }
    evidence_identities = {
        **{name: published.bundle.identity for name, published in resolved.items()},
        **{item.evidence: item.identity for item in native},
    }

    coverage = _qualified_coverage_items(
        manifest, evidence_identities, qualified_requirements
    )
    if len(coverage) != len(manifest.requirements):
        raise MigrationCoreCoverageError(
            "Migration-core has missing requirement coverage"
        )
    return CompiledMigrationCoreCoverage(
        manifest.identity,
        coverage,
        tuple(
            _compiled_evidence(evidence, resolved[evidence])
            for evidence in sorted(resolved)
        ),
        native,
    )


@dataclass(frozen=True, slots=True, order=True)
class MigrationCoreSupportingEvidenceSelection:
    """One immutable repository file selected as supporting evidence."""

    evidence: str
    locator: str
    file_hash: str
    identity: str

    def __post_init__(self) -> None:
        if self.evidence not in _APPROVED_SUPPORTING_EVIDENCE:
            raise MigrationCoreCoverageError(
                "Migration-core selects unapproved supporting evidence"
            )
        _relative_locator(self.locator, "supporting evidence locator")
        _digest(self.file_hash, "supporting evidence file hash")
        _digest(self.identity, "supporting evidence identity")

    def to_record(self) -> dict[str, object]:
        return {
            "evidence": self.evidence,
            "file_hash": self.file_hash,
            "identity": self.identity,
            "locator": self.locator,
        }

    @classmethod
    def from_record(
        cls, record: Mapping[str, object]
    ) -> MigrationCoreSupportingEvidenceSelection:
        _exact_keys(
            record,
            {"evidence", "file_hash", "identity", "locator"},
            "supporting evidence selection",
        )
        return cls(
            _text(record.get("evidence"), "supporting evidence name"),
            _relative_locator(record.get("locator"), "supporting evidence locator"),
            _digest(record.get("file_hash"), "supporting evidence file hash"),
            _digest(record.get("identity"), "supporting evidence identity"),
        )


@dataclass(frozen=True, slots=True)
class MigrationCoreCurrentReleaseSelection:
    """Repository-frozen locators and hashes for the one current phased release."""

    schema_version: int
    selection_version: str
    authority_locator: str
    authority_file_hash: str
    authority_identity: str
    status_locator: str
    status_file_hash: str
    status_identity: str
    anchor_release_locator: str
    anchor_release_hash: str
    anchor_release_identity: str
    probe_locator: str
    probe_file_hash: str
    probe_identity: str
    lifecycle_probe_locator: str | None = None
    lifecycle_probe_file_hash: str | None = None
    lifecycle_probe_identity: str | None = None
    lifecycle_source_commit: str | None = None
    lifecycle_lockfile_hash: str | None = None
    supporting_evidence: tuple[MigrationCoreSupportingEvidenceSelection, ...] = ()
    operational_source_commit: str | None = None
    operational_lockfile_hash: str | None = None
    identity: str = field(init=False)

    def _validate_extensions(self, lifecycle: tuple[str | None, ...]) -> None:
        expected_supporting = (
            "execution:covariance-constrained-uncertainty",
            "execution:resampling",
            "probe:serialization-multiprocessing-cache-replay",
        )
        operational_source = (
            self.operational_source_commit,
            self.operational_lockfile_hash,
        )
        if self.selection_version == "migration-core-current-release-v3":
            if (
                tuple(item.evidence for item in self.supporting_evidence)
                != expected_supporting
                or tuple(sorted(self.supporting_evidence)) != self.supporting_evidence
                or not all(operational_source)
                or not all(lifecycle)
                or self.operational_source_commit != _OPERATIONAL_SOURCE_COMMIT
                or self.operational_lockfile_hash != _OPERATIONAL_LOCKFILE_HASH
            ):
                raise MigrationCoreCoverageError(
                    "Migration-core current release has incomplete qualified evidence"
                )
            _git_commit(self.operational_source_commit, "operational source commit")
            _digest(self.operational_lockfile_hash, "operational lockfile hash")
        elif self.supporting_evidence or any(operational_source):
            raise MigrationCoreCoverageError(
                "Historical migration-core release selects later evidence"
            )

    def __post_init__(self) -> None:
        if (
            self.schema_version != _SCHEMA_VERSION
            or self.selection_version not in _CURRENT_RELEASE_VERSIONS
        ):
            raise MigrationCoreCoverageError(
                "Migration-core current release has incompatible version"
            )
        for value, name in (
            (self.authority_locator, "authority locator"),
            (self.status_locator, "status locator"),
            (self.anchor_release_locator, "anchor release locator"),
            (self.probe_locator, "probe locator"),
        ):
            _relative_locator(value, name)
        for value, name in (
            (self.authority_file_hash, "authority file hash"),
            (self.authority_identity, "authority identity"),
            (self.status_file_hash, "status file hash"),
            (self.status_identity, "status identity"),
            (self.anchor_release_hash, "anchor release hash"),
            (self.anchor_release_identity, "anchor release identity"),
            (self.probe_file_hash, "probe file hash"),
            (self.probe_identity, "probe identity"),
        ):
            _digest(value, name)
        lifecycle = (
            self.lifecycle_probe_locator,
            self.lifecycle_probe_file_hash,
            self.lifecycle_probe_identity,
            self.lifecycle_source_commit,
            self.lifecycle_lockfile_hash,
        )
        if any(value is not None for value in lifecycle):
            if not all(value is not None for value in lifecycle):
                raise MigrationCoreCoverageError(
                    "Migration-core lifecycle selection is incomplete"
                )
            _relative_locator(self.lifecycle_probe_locator, "lifecycle probe locator")
            _digest(self.lifecycle_probe_file_hash, "lifecycle probe file hash")
            _digest(self.lifecycle_probe_identity, "lifecycle probe identity")
            _git_commit(self.lifecycle_source_commit, "lifecycle source commit")
            _digest(self.lifecycle_lockfile_hash, "lifecycle lockfile hash")
        if (self.selection_version.endswith("v2")) != bool(
            lifecycle[0]
        ) and self.selection_version != "migration-core-current-release-v3":
            raise MigrationCoreCoverageError(
                "Migration-core release version differs from lifecycle selection"
            )
        self._validate_extensions(lifecycle)
        object.__setattr__(self, "identity", _content_identity(self.to_record()))

    def to_record(self) -> dict[str, object]:
        record: dict[str, object] = {
            "anchor_release": {
                "archive_hash": self.anchor_release_hash,
                "identity": self.anchor_release_identity,
                "locator": self.anchor_release_locator,
            },
            "authority_selection": {
                "file_hash": self.authority_file_hash,
                "identity": self.authority_identity,
                "locator": self.authority_locator,
            },
            "numerical_probe": {
                "file_hash": self.probe_file_hash,
                "identity": self.probe_identity,
                "locator": self.probe_locator,
            },
            "phased_status": {
                "file_hash": self.status_file_hash,
                "identity": self.status_identity,
                "locator": self.status_locator,
            },
            "schema_version": self.schema_version,
            "selection_version": self.selection_version,
        }
        if self.lifecycle_probe_locator is not None:
            record["lifecycle_probe"] = {
                "file_hash": self.lifecycle_probe_file_hash,
                "identity": self.lifecycle_probe_identity,
                "locator": self.lifecycle_probe_locator,
                "lockfile_hash": self.lifecycle_lockfile_hash,
                "source_commit": self.lifecycle_source_commit,
            }
        if self.supporting_evidence:
            record["supporting_evidence"] = [
                item.to_record() for item in self.supporting_evidence
            ]
            record["operational_source"] = {
                "lockfile_hash": self.operational_lockfile_hash,
                "source_commit": self.operational_source_commit,
            }
        return record

    @classmethod
    def from_bytes(cls, content: bytes) -> MigrationCoreCurrentReleaseSelection:
        record = _json_record_from_bytes(content, "current release selection")
        expected = {
            "anchor_release",
            "authority_selection",
            "identity",
            "numerical_probe",
            "phased_status",
            "schema_version",
            "selection_version",
        }
        if "lifecycle_probe" in record:
            expected.add("lifecycle_probe")
        if "supporting_evidence" in record or "operational_source" in record:
            expected.update({"supporting_evidence", "operational_source"})
        _exact_keys(record, expected, "current release selection")
        authority = _record(record.get("authority_selection"), "authority selection")
        status = _record(record.get("phased_status"), "phased status selection")
        release = _record(record.get("anchor_release"), "anchor release selection")
        probe = _record(record.get("numerical_probe"), "probe selection")
        lifecycle = (
            _record(record.get("lifecycle_probe"), "lifecycle probe selection")
            if "lifecycle_probe" in record
            else None
        )
        raw_supporting = record.get("supporting_evidence", [])
        operational_source = (
            _record(record.get("operational_source"), "operational source")
            if "operational_source" in record
            else None
        )
        if not isinstance(raw_supporting, list):
            raise MigrationCoreCoverageError(
                "Migration-core supporting evidence selection must be a list"
            )
        for nested, name in (
            (authority, "authority selection"),
            (status, "phased status selection"),
            (release, "anchor release selection"),
            (probe, "probe selection"),
        ):
            _exact_keys(
                nested,
                {"file_hash", "identity", "locator"}
                if name != "anchor release selection"
                else {"archive_hash", "identity", "locator"},
                name,
            )
        if lifecycle is not None:
            _exact_keys(
                lifecycle,
                {"file_hash", "identity", "locator", "lockfile_hash", "source_commit"},
                "lifecycle probe selection",
            )
        if operational_source is not None:
            _exact_keys(
                operational_source,
                {"lockfile_hash", "source_commit"},
                "operational source",
            )
        schema_version = record.get("schema_version")
        if isinstance(schema_version, bool) or not isinstance(schema_version, int):
            raise MigrationCoreCoverageError(
                "Migration-core current release schema must be an integer"
            )
        selection = cls(
            schema_version,
            _text(record.get("selection_version"), "current release version"),
            _relative_locator(authority.get("locator"), "authority locator"),
            _digest(authority.get("file_hash"), "authority file hash"),
            _digest(authority.get("identity"), "authority identity"),
            _relative_locator(status.get("locator"), "status locator"),
            _digest(status.get("file_hash"), "status file hash"),
            _digest(status.get("identity"), "status identity"),
            _relative_locator(release.get("locator"), "anchor release locator"),
            _digest(release.get("archive_hash"), "anchor release hash"),
            _digest(release.get("identity"), "anchor release identity"),
            _relative_locator(probe.get("locator"), "probe locator"),
            _digest(probe.get("file_hash"), "probe file hash"),
            _digest(probe.get("identity"), "probe identity"),
            None
            if lifecycle is None
            else _relative_locator(lifecycle.get("locator"), "lifecycle probe locator"),
            None
            if lifecycle is None
            else _digest(lifecycle.get("file_hash"), "lifecycle probe file hash"),
            None
            if lifecycle is None
            else _digest(lifecycle.get("identity"), "lifecycle probe identity"),
            None
            if lifecycle is None
            else _git_commit(lifecycle.get("source_commit"), "lifecycle source commit"),
            None
            if lifecycle is None
            else _digest(lifecycle.get("lockfile_hash"), "lifecycle lockfile hash"),
            tuple(
                MigrationCoreSupportingEvidenceSelection.from_record(
                    _record(item, "supporting evidence selection")
                )
                for item in raw_supporting
            ),
            None
            if operational_source is None
            else _git_commit(
                operational_source.get("source_commit"), "operational source commit"
            ),
            None
            if operational_source is None
            else _digest(
                operational_source.get("lockfile_hash"), "operational lockfile hash"
            ),
        )
        if record.get("identity") != selection.identity:
            raise MigrationCoreCoverageError(
                "Migration-core current-release identity does not match payload"
            )
        return selection


@dataclass(frozen=True, slots=True)
class MigrationCorePhasedStatus:
    """Index and expected terminal reproduced from durable scientific evidence."""

    schema_version: int
    record_version: str
    authority_locator: str
    authority_identity: str
    anchor_release_locator: str
    anchor_release_hash: str
    anchor_release_identity: str
    source_commit: str
    source_archive_locator: str
    source_archive_hash: str
    anchors: tuple[MigrationCoreReleaseAnchor, ...]
    probe_locator: str
    probe_file_hash: str
    probe_identity: str
    probe_manifest_identity: str
    selected_evidence: tuple[tuple[str, str], ...]
    eligible_requirement_ids: tuple[str, ...]
    uncovered_requirement_ids: tuple[str, ...]
    unqualified_requirement_ids: tuple[str, ...]
    compiler_status: str
    compiler_reason: str
    compiled_coverage_identity: None
    completion_claim: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        expected_claim = (
            "MIGRATION_CORE_MATRIX_COMPLETE_WITH_EXPLICIT_UNCOVERED_SCOPE"
            if self.record_version == "migration-core-canonical-coverage-status-v5"
            else "PHASED_INFRASTRUCTURE_ONLY"
        )
        if (
            self.schema_version != _SCHEMA_VERSION
            or self.record_version not in _CURRENT_STATUS_VERSIONS
            or self.completion_claim != expected_claim
            or self.compiler_status != "FAILED_CLOSED"
            or self.compiled_coverage_identity is not None
        ):
            raise MigrationCoreCoverageError(
                "Migration-core phased status has incompatible semantics"
            )
        for value, name in (
            (self.authority_locator, "status authority locator"),
            (self.anchor_release_locator, "status release locator"),
            (self.source_archive_locator, "status source locator"),
            (self.probe_locator, "status probe locator"),
        ):
            _relative_locator(value, name)
        for value, name in (
            (self.authority_identity, "status authority identity"),
            (self.anchor_release_hash, "status release hash"),
            (self.anchor_release_identity, "status release identity"),
            (self.source_archive_hash, "status source hash"),
            (self.probe_file_hash, "status probe hash"),
            (self.probe_identity, "status probe identity"),
            (self.probe_manifest_identity, "status probe manifest identity"),
        ):
            _digest(value, name)
        _git_commit(self.source_commit, "status source commit")
        if (
            tuple(item.evidence for item in self.anchors) != _APPROVED_EVIDENCE
            or tuple(sorted(self.selected_evidence)) != self.selected_evidence
            or len({name for name, _ in self.selected_evidence})
            != len(self.selected_evidence)
            or (
                len(self.eligible_requirement_ids),
                len(self.uncovered_requirement_ids),
            )
            != _CURRENT_STATUS_VERSIONS[self.record_version]
            or set(self.eligible_requirement_ids) & set(self.uncovered_requirement_ids)
            or not set(self.unqualified_requirement_ids)
            <= set(self.uncovered_requirement_ids)
            or (
                self.record_version == "migration-core-canonical-coverage-status-v5"
                and set(self.unqualified_requirement_ids)
                != _EXPECTED_UNQUALIFIED_REQUIREMENTS
            )
            or (
                self.record_version != "migration-core-canonical-coverage-status-v5"
                and self.unqualified_requirement_ids
            )
        ):
            raise MigrationCoreCoverageError(
                "Migration-core phased status has malformed coverage selections"
            )
        object.__setattr__(self, "identity", _content_identity(self.to_record()))

    def to_record(self) -> dict[str, object]:
        coverage: dict[str, object] = {
            "compiled_coverage_identity": self.compiled_coverage_identity,
            "compiler": {
                "reason": self.compiler_reason,
                "status": self.compiler_status,
            },
            "eligible_requirement_ids": list(self.eligible_requirement_ids),
            "requirement_count": len(_AUDITED_REQUIREMENTS),
            "selected_evidence": [
                {"evidence": name, "identity": identity}
                for name, identity in self.selected_evidence
            ],
            "uncovered_requirement_ids": list(self.uncovered_requirement_ids),
        }
        if self.unqualified_requirement_ids:
            coverage["unqualified_requirement_ids"] = list(
                self.unqualified_requirement_ids
            )
        return {
            "anchor_release": {
                "anchors": [item.to_record() for item in self.anchors],
                "archive_hash": self.anchor_release_hash,
                "identity": self.anchor_release_identity,
                "locator": self.anchor_release_locator,
                "source_archive_hash": self.source_archive_hash,
                "source_archive_locator": self.source_archive_locator,
                "source_commit": self.source_commit,
            },
            "authority_selection": {
                "identity": self.authority_identity,
                "locator": self.authority_locator,
            },
            "completion_claim": self.completion_claim,
            "coverage": coverage,
            "numerical_probe": {
                "file_hash": self.probe_file_hash,
                "identity": self.probe_identity,
                "locator": self.probe_locator,
                "manifest_identity": self.probe_manifest_identity,
            },
            "record_version": self.record_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_bytes(cls, content: bytes) -> MigrationCorePhasedStatus:
        record = _json_record_from_bytes(content, "phased status")
        _exact_keys(
            record,
            {
                "anchor_release",
                "authority_selection",
                "completion_claim",
                "coverage",
                "identity",
                "numerical_probe",
                "record_version",
                "schema_version",
            },
            "phased status",
        )
        authority = _record(record.get("authority_selection"), "status authority")
        release = _record(record.get("anchor_release"), "status anchor release")
        probe = _record(record.get("numerical_probe"), "status probe")
        coverage = _record(record.get("coverage"), "status coverage")
        compiler = _record(coverage.get("compiler"), "status compiler")
        record_version = _text(record.get("record_version"), "status version")
        _exact_keys(authority, {"identity", "locator"}, "status authority")
        _exact_keys(
            release,
            {
                "anchors",
                "archive_hash",
                "identity",
                "locator",
                "source_archive_hash",
                "source_archive_locator",
                "source_commit",
            },
            "status anchor release",
        )
        _exact_keys(
            probe,
            {"file_hash", "identity", "locator", "manifest_identity"},
            "status probe",
        )
        coverage_keys = {
            "compiled_coverage_identity",
            "compiler",
            "eligible_requirement_ids",
            "requirement_count",
            "selected_evidence",
            "uncovered_requirement_ids",
        }
        if record_version == "migration-core-canonical-coverage-status-v5":
            coverage_keys.add("unqualified_requirement_ids")
        _exact_keys(coverage, coverage_keys, "status coverage")
        _exact_keys(compiler, {"reason", "status"}, "status compiler")
        raw_anchors = release.get("anchors")
        raw_selected = coverage.get("selected_evidence")
        raw_eligible = coverage.get("eligible_requirement_ids")
        raw_uncovered = coverage.get("uncovered_requirement_ids")
        raw_unqualified = coverage.get("unqualified_requirement_ids", [])
        schema_version = record.get("schema_version")
        if (
            not isinstance(raw_anchors, list)
            or not isinstance(raw_selected, list)
            or not isinstance(raw_eligible, list)
            or not isinstance(raw_uncovered, list)
            or not isinstance(raw_unqualified, list)
            or isinstance(schema_version, bool)
            or not isinstance(schema_version, int)
            or coverage.get("requirement_count") != len(_AUDITED_REQUIREMENTS)
        ):
            raise MigrationCoreCoverageError(
                "Migration-core phased status has malformed fields"
            )
        selections: list[tuple[str, str]] = []
        for item in raw_selected:
            selected = _record(item, "status evidence selection")
            _exact_keys(selected, {"evidence", "identity"}, "status evidence selection")
            selections.append(
                (
                    _text(selected.get("evidence"), "selected evidence"),
                    _digest(selected.get("identity"), "selected evidence identity"),
                )
            )
        status = cls(
            schema_version,
            record_version,
            _relative_locator(authority.get("locator"), "status authority locator"),
            _digest(authority.get("identity"), "status authority identity"),
            _relative_locator(release.get("locator"), "status release locator"),
            _digest(release.get("archive_hash"), "status release hash"),
            _digest(release.get("identity"), "status release identity"),
            _git_commit(release.get("source_commit"), "status source commit"),
            _relative_locator(
                release.get("source_archive_locator"), "status source locator"
            ),
            _digest(release.get("source_archive_hash"), "status source hash"),
            tuple(
                MigrationCoreReleaseAnchor.from_record(
                    _record(item, "status release anchor")
                )
                for item in raw_anchors
            ),
            _relative_locator(probe.get("locator"), "status probe locator"),
            _digest(probe.get("file_hash"), "status probe hash"),
            _digest(probe.get("identity"), "status probe identity"),
            _digest(probe.get("manifest_identity"), "status probe manifest identity"),
            tuple(selections),
            tuple(_text(item, "eligible requirement") for item in raw_eligible),
            tuple(_text(item, "uncovered requirement") for item in raw_uncovered),
            tuple(_text(item, "unqualified requirement") for item in raw_unqualified),
            _text(compiler.get("status"), "status compiler terminal"),
            _text(compiler.get("reason"), "status compiler reason"),
            cast("None", coverage.get("compiled_coverage_identity")),
            _text(record.get("completion_claim"), "status completion claim"),
        )
        if record.get("identity") != status.identity:
            raise MigrationCoreCoverageError(
                "Migration-core phased-status identity does not match payload"
            )
        return status


@dataclass(frozen=True, slots=True)
class DerivedMigrationCorePhasedCoverage:
    """Coverage facts derived only after durable evidence reconstruction."""

    current_release_identity: str
    authority_selection_identity: str
    anchor_release_identity: str
    anchor_evidence: tuple[CompiledMigrationCoreEvidence, ...]
    numerical_probe: CompiledMigrationCoreSupportingEvidence
    eligible_requirement_ids: tuple[str, ...]
    uncovered_requirement_ids: tuple[str, ...]
    unqualified_requirement_ids: tuple[str, ...]
    compiler_status: str
    compiler_reason: str
    compiled_coverage_identity: None = None


def migration_core_current_release_selection() -> MigrationCoreCurrentReleaseSelection:
    """Load the sole repository-frozen current phased release selection."""

    try:
        package = files("chemex")
        resource = next(
            item
            for name in _CURRENT_RELEASE_SELECTION_RESOURCES
            if (item := package.joinpath(name)).is_file()
        )
        content = resource.read_bytes()
    except (OSError, StopIteration) as error:
        raise MigrationCoreCoverageError(
            "Migration-core current release selection is unavailable"
        ) from error
    return MigrationCoreCurrentReleaseSelection.from_bytes(content)


def _read_selected_file(
    repository_root: Path, locator: str, expected_hash: str, name: str
) -> tuple[Path, bytes]:
    path = _repository_path(repository_root, locator, name)
    try:
        content = path.read_bytes()
    except OSError as error:
        raise MigrationCoreCoverageError(
            f"Migration-core selected {name} is missing"
        ) from error
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise MigrationCoreCoverageError(
            f"Migration-core selected {name} hash does not match"
        )
    return path, content


def _load_selected_authority_and_status(
    repository_root: Path, current: MigrationCoreCurrentReleaseSelection
) -> tuple[MigrationCoreAuthoritySelection, MigrationCorePhasedStatus]:
    _, authority_content = _read_selected_file(
        repository_root,
        current.authority_locator,
        current.authority_file_hash,
        "authority selection",
    )
    authority = MigrationCoreAuthoritySelection.from_bytes(authority_content)
    frozen_authority = migration_core_authority_selection()
    if (
        authority != frozen_authority
        or authority.identity != current.authority_identity
    ):
        raise MigrationCoreCoverageError(
            "Migration-core current release selects incompatible authority"
        )
    _, status_content = _read_selected_file(
        repository_root,
        current.status_locator,
        current.status_file_hash,
        "phased status",
    )
    status = MigrationCorePhasedStatus.from_bytes(status_content)
    if status.identity != current.status_identity:
        raise MigrationCoreCoverageError(
            "Migration-core phased status differs from current selection"
        )
    if (
        status.authority_locator != current.authority_locator
        or status.authority_identity != authority.identity
        or status.anchor_release_locator != current.anchor_release_locator
        or status.anchor_release_hash != current.anchor_release_hash
        or status.anchor_release_identity != current.anchor_release_identity
        or status.probe_locator != current.probe_locator
        or status.probe_file_hash != current.probe_file_hash
        or status.probe_identity != current.probe_identity
    ):
        raise MigrationCoreCoverageError(
            "Migration-core phased status has incompatible frozen selections"
        )
    return authority, status


def _load_selected_probe(
    repository_root: Path,
    current: MigrationCoreCurrentReleaseSelection,
    status: MigrationCorePhasedStatus,
) -> tuple[NumericalProbeBaseline, tuple[str, str]]:
    _, probe_content = _read_selected_file(
        repository_root,
        current.probe_locator,
        current.probe_file_hash,
        "numerical probe",
    )
    try:
        probe_record = _json_record_from_bytes(probe_content, "numerical probe")
        probe = NumericalProbeBaseline.from_record(probe_record)
    except (TypeError, ValueError) as error:
        raise MigrationCoreCoverageError(
            "Migration-core numerical probe cannot be reconstructed"
        ) from error
    validated_probe = _validated_numerical_probe_identity(probe)
    if (
        validated_probe is None
        or probe.identity != current.probe_identity
        or probe.identity != status.probe_identity
        or probe.manifest_identity != status.probe_manifest_identity
    ):
        raise MigrationCoreCoverageError(
            "Migration-core numerical probe has an exact authority mismatch"
        )
    return probe, validated_probe


def _load_selected_lifecycle_probe(
    repository_root: Path,
    current: MigrationCoreCurrentReleaseSelection,
) -> LifecycleProbeCapture | None:
    if current.lifecycle_probe_locator is None:
        return None
    lifecycle_hash = current.lifecycle_probe_file_hash
    if lifecycle_hash is None:
        raise MigrationCoreCoverageError("Lifecycle probe selection is incomplete")
    _, content = _read_selected_file(
        repository_root,
        current.lifecycle_probe_locator,
        lifecycle_hash,
        "lifecycle probe capture",
    )
    try:
        capture = LifecycleProbeCapture.from_bytes(content)
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise MigrationCoreCoverageError(
            "Migration-core lifecycle probe capture cannot be reconstructed"
        ) from error
    validated = _validated_lifecycle_publication_safety_identity(capture)
    if validated is None or capture.identity != current.lifecycle_probe_identity:
        raise MigrationCoreCoverageError(
            "Migration-core lifecycle probe capture is not eligible"
        )
    return capture


def _load_selected_supporting_evidence(
    repository_root: Path,
    current: MigrationCoreCurrentReleaseSelection,
) -> dict[str, object]:
    supporting: dict[str, object] = {}
    for selection in current.supporting_evidence:
        _, content = _read_selected_file(
            repository_root,
            selection.locator,
            selection.file_hash,
            f"{selection.evidence} selection",
        )
        if selection.evidence == "probe:serialization-multiprocessing-cache-replay":
            try:
                evidence: object = OperationalReplayCapture.from_bytes(content)
            except (TypeError, ValueError, json.JSONDecodeError) as error:
                raise MigrationCoreCoverageError(
                    "Migration-core operational replay cannot be reconstructed"
                ) from error
        else:
            evidence = content
        _evidence_type, identity = _supporting_native_identity(
            selection.evidence, evidence
        )
        if identity != selection.identity:
            raise MigrationCoreCoverageError(
                f"Evidence {selection.evidence} differs from its frozen selection"
            )
        supporting[selection.evidence] = evidence
    return supporting


def _validate_resolved_release_status(
    resolved: ResolvedMigrationCoreAnchorRelease,
    status: MigrationCorePhasedStatus,
) -> None:
    if (
        resolved.release.repository_commit != status.source_commit
        or resolved.release.source_archive_locator != status.source_archive_locator
        or resolved.release.source_archive_hash != status.source_archive_hash
        or resolved.release.anchors != status.anchors
    ):
        raise MigrationCoreCoverageError(
            "Migration-core anchor release differs from phased status"
        )


def _derive_release_anchor_evidence(
    resolved: ResolvedMigrationCoreAnchorRelease,
    manifest: MigrationCoreCoverageManifest,
) -> tuple[dict[str, str], tuple[CompiledMigrationCoreEvidence, ...]]:
    expected_anchors = {item.evidence: item for item in manifest.anchors}
    selections: dict[str, str] = {}
    compiled_anchors: list[CompiledMigrationCoreEvidence] = []
    for anchor in resolved.release.anchors:
        published = resolved.published[anchor.evidence]
        _validate_anchor_evidence(expected_anchors[anchor.evidence], published)
        selections[anchor.evidence] = published.bundle.identity
        compiled_anchors.append(_compiled_evidence(anchor.evidence, published))
    return selections, tuple(compiled_anchors)


def _derive_requirement_partitions(
    manifest: MigrationCoreCoverageManifest,
    selections: Mapping[str, str],
    supporting_evidence: Mapping[str, object],
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    selected_names = set(selections) | set(supporting_evidence)
    qualified = {
        name: _supporting_qualified_requirements(name, evidence)
        for name, evidence in supporting_evidence.items()
    }
    eligible_items: list[str] = []
    uncovered_items: list[str] = []
    unqualified_items: list[str] = []
    for requirement in manifest.requirements:
        required = set(requirement.evidence)
        if not required <= selected_names:
            uncovered_items.append(requirement.identifier)
            continue
        unqualified = tuple(
            evidence
            for evidence in requirement.evidence
            if evidence in qualified
            and requirement.identifier not in qualified[evidence]
        )
        if unqualified:
            uncovered_items.append(requirement.identifier)
            unqualified_items.append(requirement.identifier)
        else:
            eligible_items.append(requirement.identifier)
    return tuple(eligible_items), tuple(uncovered_items), tuple(unqualified_items)


def _fail_closed_compiler_reason(
    manifest: MigrationCoreCoverageManifest,
    resolved: ResolvedMigrationCoreAnchorRelease,
    selections: Mapping[str, str],
    supporting_evidence: Mapping[str, object],
) -> str:
    try:
        compile_migration_core_coverage(
            manifest,
            resolved.publisher,
            selections,
            supporting_evidence,
        )
    except MigrationCoreCoverageError as error:
        return str(error)
    raise MigrationCoreCoverageError(
        "Migration-core incomplete phased evidence unexpectedly compiled"
    )


def compile_current_phased_migration_core_status(
    repository_root: Path,
) -> DerivedMigrationCorePhasedCoverage:
    """Reconstruct durable evidence and derive the current phased #590 status."""

    current = migration_core_current_release_selection()
    authority, status = _load_selected_authority_and_status(repository_root, current)
    probe, validated_probe = _load_selected_probe(repository_root, current, status)
    lifecycle_probe = _load_selected_lifecycle_probe(repository_root, current)
    supporting_evidence = _load_selected_supporting_evidence(repository_root, current)
    supporting_evidence["probe:numerical-optimizer-evaluator"] = probe
    if lifecycle_probe is not None:
        supporting_evidence["probe:lifecycle-publication-safety"] = lifecycle_probe
    archive_path = _repository_path(
        repository_root, current.anchor_release_locator, "anchor release"
    )
    manifest = migration_core_coverage_manifest()
    with resolve_migration_core_anchor_release(
        archive_path,
        expected_archive_hash=current.anchor_release_hash,
        expected_release_identity=current.anchor_release_identity,
    ) as resolved:
        _validate_resolved_release_status(resolved, status)
        selections, compiled_anchors = _derive_release_anchor_evidence(
            resolved, manifest
        )
        supporting_identities = [
            (name, _supporting_native_identity(name, evidence)[1])
            for name, evidence in supporting_evidence.items()
        ]
        selected_identities = tuple(
            sorted((*selections.items(), *supporting_identities))
        )
        if selected_identities != status.selected_evidence:
            raise MigrationCoreCoverageError(
                "Migration-core requirement-to-evidence selection is incompatible"
            )
        eligible, uncovered, unqualified = _derive_requirement_partitions(
            manifest, selections, supporting_evidence
        )
        if (
            eligible != status.eligible_requirement_ids
            or uncovered != status.uncovered_requirement_ids
            or unqualified != status.unqualified_requirement_ids
        ):
            raise MigrationCoreCoverageError(
                "Migration-core phased coverage counts are not evidence-derived"
            )
        compiler_reason = _fail_closed_compiler_reason(
            manifest, resolved, selections, supporting_evidence
        )
        if (
            status.compiler_status != "FAILED_CLOSED"
            or status.compiler_reason != compiler_reason
        ):
            raise MigrationCoreCoverageError(
                "Migration-core compiler terminal differs from phased status"
            )
        probe_reference = CompiledMigrationCoreSupportingEvidence(
            "probe:numerical-optimizer-evaluator",
            validated_probe[0],
            validated_probe[1],
        )
        return DerivedMigrationCorePhasedCoverage(
            current.identity,
            authority.identity,
            resolved.release.identity,
            compiled_anchors,
            probe_reference,
            eligible,
            uncovered,
            unqualified,
            "FAILED_CLOSED",
            compiler_reason,
        )
