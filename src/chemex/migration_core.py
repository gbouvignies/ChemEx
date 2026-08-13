"""Fixed #590 migration-core requirement-to-evidence composition."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from importlib.resources import files
from typing import cast

from chemex.baselines import (
    BaselinePublisher,
    CanonicalBaselineValue,
    PublishedEvidence,
)
from chemex.numerical_lanes import LaneAttestation, canonical_lanes
from chemex.optimize.numerical_probes import (
    CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    NumericalProbeBaseline,
)

_SCHEMA_VERSION = 1
_MANIFEST_VERSION = "migration-core-coverage-v1"
_BASELINE_SEMANTIC_VERSION = "chemex-baseline-v1"
_MANIFEST_RESOURCE = "migration_core_coverage_v1.json"
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


def _ordered_texts(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise MigrationCoreCoverageError(f"Migration-core {name} must be a list")
    result = tuple(_text(item, name) for item in value)
    if not result or result != tuple(sorted(result)) or len(set(result)) != len(result):
        raise MigrationCoreCoverageError(
            f"Migration-core {name} must contain ordered unique values"
        )
    return result


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


def _validate_anchor_evidence(
    expected: MigrationCoreAnchorEvidence, published: PublishedEvidence
) -> None:
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
        or published.specification.lane_reference != canonical_lane.identity
        or published.occurrence.lane_reference != canonical_lane.identity
        or published.occurrence.lane_attestation_identity is None
        or settings.get("workers") != 1
        or settings.get("native_threads") != 1
        or published.occurrence.lifecycle != "SUCCEEDED"
        or published.bundle.implementation.authority_role
        != "LegacyObservationImplementation"
    ):
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} is incompatible with the coverage manifest"
        )
    attestation_members = tuple(
        member
        for member in published.bundle.members
        if member.role == "environment:lane-attestation.json"
    )
    if len(attestation_members) != 1:
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} lacks exact environment evidence"
        )
    try:
        content = (
            published.location / "members" / attestation_members[0].content_hash
        ).read_bytes()
        attestation = LaneAttestation.from_record(
            _record(json.loads(content), "lane attestation")
        )
    except (OSError, TypeError, ValueError) as error:
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} has unavailable environment evidence"
        ) from error
    if (
        attestation.identity != published.occurrence.lane_attestation_identity
        or attestation.lane_identity != canonical_lane.identity
    ):
        raise MigrationCoreCoverageError(
            f"Evidence {expected.evidence} has incompatible environment identity"
        )


def _numerical_probe_identity(evidence: object) -> tuple[str, str] | None:
    if (
        isinstance(evidence, NumericalProbeBaseline)
        and evidence.historical_qualification == "REFERENCE_MATCHED"
        and evidence.manifest_identity == CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
        and evidence.reference_manifest_identity
        == CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY
        and evidence.observed_lane_role == "CANONICAL_NUMERICAL"
        and evidence.observed_lane_identity == canonical_lanes()[0].identity
    ):
        return type(evidence).__name__, evidence.identity
    return None


def _unavailable_identity(_evidence: object) -> tuple[str, str] | None:
    return None


_SUPPORTING_VALIDATORS: dict[str, Callable[[object], tuple[str, str] | None]] = {
    # Exact reduced-case and policy fingerprints are empirical release inputs.
    # Broad successful native objects deliberately cannot mint these selections.
    "execution:covariance-constrained-uncertainty": _unavailable_identity,
    "execution:de-trf-search": _unavailable_identity,
    "execution:fit-component-aggregation": _unavailable_identity,
    "execution:grid-decomposed": _unavailable_identity,
    "execution:grid-single-component": _unavailable_identity,
    "execution:mcmc": _unavailable_identity,
    "execution:resampling": _unavailable_identity,
    "probe:lifecycle-publication-safety": _unavailable_identity,
    "probe:numerical-optimizer-evaluator": _numerical_probe_identity,
    "probe:serialization-multiprocessing-cache-replay": _unavailable_identity,
    "workload:qualified-performance": _unavailable_identity,
}


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
    evidence_identities = {
        **{name: published.bundle.identity for name, published in resolved.items()},
        **{item.evidence: item.identity for item in native},
    }

    coverage = tuple(
        (
            requirement.identifier,
            tuple(evidence_identities[evidence] for evidence in requirement.evidence),
        )
        for requirement in manifest.requirements
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
