"""Truth-backed reduced numerical probes for migration baseline qualification.

The probes in this module are calibration evidence, not product fit entry points.
They freeze independently derived truth, closed eligibility claims, numerical
policy, work budgets, seeds, and failure expectations before any probe runs.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import product
from typing import Literal, Never, cast

import numpy as np
from scipy.linalg import svd
from scipy.optimize import Bounds, OptimizeResult, differential_evolution, least_squares

from chemex.baselines import CanonicalBaselineValue
from chemex.numerical_lanes import LaneRole, LiveLaneAuthority, canonical_lanes
from chemex.typing import Array

_SCHEMA_VERSION = 2
_SEMANTIC_VERSION = "chemex-numerical-probes-v2"
_SOURCE_REVISION = "700cb71df950fc54c268c0bca199403e5621269d"
CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY = (
    "c8052d59d8b648bd0ba55effdd91d86cc18bf68b948faa420450ede3d47b93bf"
)

type ProbeCategory = Literal[
    "TRF_ROUTINE",
    "TRF_DIFFICULT",
    "GRID",
    "GRID_ORDERING",
    "DE_SEARCH",
    "FINITE_DIFFERENCE",
    "BOUNDS",
    "RANK",
    "CONDITIONING",
]
type TruthAuthorityKind = Literal["ANALYTIC_DERIVATION", "EXACT_LINEAR_ALGEBRA"]
type RequestFailureKind = Literal[
    "INVALID_INPUT",
    "EVALUATION_FAILURE",
    "NON_FINITE_RESULT",
]


def _identity(kind: str, record: object) -> str:
    content = json.dumps(
        {
            "kind": kind,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "record": record,
        },
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(content).hexdigest()


def _float_token(value: float) -> str:
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError("Probe values must be finite binary64 values")
    return (0.0 if scalar == 0.0 else scalar).hex()


def _float_from_record(value: object, field_name: str) -> float:
    if not isinstance(value, str):
        raise TypeError(f"Probe {field_name} must be a binary64 token")
    try:
        scalar = float.fromhex(value)
    except ValueError as error:
        raise ValueError(f"Probe {field_name} is not binary64") from error
    if _float_token(scalar) != value:
        raise ValueError(f"Probe {field_name} is not canonical binary64")
    return scalar


def _vector_tokens(values: Sequence[float]) -> tuple[str, ...]:
    return tuple(_float_token(value) for value in values)


def _vector_from_record(value: object, field_name: str) -> tuple[float, ...]:
    if not isinstance(value, list) or not value:
        raise TypeError(f"Probe {field_name} must be a non-empty list")
    return tuple(
        _float_from_record(item, f"{field_name}[{index}]")
        for index, item in enumerate(value)
    )


def _record_nonnegative_int(value: object, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Probe {field_name} must be a non-negative integer")
    return value


def _nonempty_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"Probe {field_name} must be a non-empty string")
    return value


def _canonical_strings(values: Sequence[str], field_name: str) -> tuple[str, ...]:
    result = tuple(values)
    if (
        not result
        or any(not isinstance(value, str) or not value for value in result)
        or tuple(sorted(result)) != result
        or len(set(result)) != len(result)
    ):
        raise ValueError(
            f"Probe {field_name} must be non-empty, unique, and canonically ordered"
        )
    return result


def _exact_keys(record: Mapping[str, object], expected: set[str], name: str) -> None:
    if set(record) != expected:
        raise ValueError(f"{name} record has unknown or missing fields")


@dataclass(frozen=True, slots=True)
class ProbeDerivation:
    """Exact source and reduction provenance for one synthetic probe."""

    source_issue: int
    source_revision: str
    parent_case: str
    reduction: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if isinstance(self.source_issue, bool) or self.source_issue < 1:
            raise ValueError("Probe source issue must be a positive integer")
        for value, name in (
            (self.source_revision, "source revision"),
            (self.parent_case, "parent case"),
            (self.reduction, "reduction"),
        ):
            _nonempty_string(value, name)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "probe-derivation",
                (
                    self.source_issue,
                    self.source_revision,
                    self.parent_case,
                    self.reduction,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "source_issue": self.source_issue,
            "source_revision": self.source_revision,
            "parent_case": self.parent_case,
            "reduction": self.reduction,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ProbeDerivation:
        _exact_keys(
            record,
            {
                "source_issue",
                "source_revision",
                "parent_case",
                "reduction",
                "identity",
            },
            "Probe derivation",
        )
        source_issue = record.get("source_issue")
        if not isinstance(source_issue, int) or isinstance(source_issue, bool):
            raise TypeError("Probe source issue must be an integer")
        derivation = cls(
            source_issue,
            _nonempty_string(record.get("source_revision"), "source revision"),
            _nonempty_string(record.get("parent_case"), "parent case"),
            _nonempty_string(record.get("reduction"), "reduction"),
        )
        if record.get("identity") != derivation.identity:
            raise ValueError("Probe derivation identity does not match its payload")
        return derivation


@dataclass(frozen=True, slots=True)
class ProbeTruthAuthority:
    """Non-observational authority that makes a reduced-case claim eligible."""

    kind: TruthAuthorityKind
    reference: str
    derivation_digest: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.kind not in {"ANALYTIC_DERIVATION", "EXACT_LINEAR_ALGEBRA"}:
            raise ValueError("Unknown probe truth authority")
        _nonempty_string(self.reference, "truth reference")
        digest = hashlib.sha256(self.reference.encode("utf-8")).hexdigest()
        object.__setattr__(self, "derivation_digest", digest)
        object.__setattr__(
            self,
            "identity",
            _identity("probe-truth-authority", (self.kind, digest)),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "reference": self.reference,
            "derivation_digest": self.derivation_digest,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ProbeTruthAuthority:
        _exact_keys(
            record,
            {"kind", "reference", "derivation_digest", "identity"},
            "Probe truth authority",
        )
        authority = cls(
            cast("TruthAuthorityKind", record.get("kind")),
            _nonempty_string(record.get("reference"), "truth reference"),
        )
        if (
            record.get("derivation_digest") != authority.derivation_digest
            or record.get("identity") != authority.identity
        ):
            raise ValueError(
                "Probe truth authority identity does not match its payload"
            )
        return authority


@dataclass(frozen=True, slots=True, order=True)
class ProbeFailureExpectation:
    """Required terminal and risk truth for one probe run."""

    terminal: str
    required_risks: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _nonempty_string(self.terminal, "expected terminal")
        risks = tuple(self.required_risks)
        if risks and (
            tuple(sorted(risks)) != risks
            or len(set(risks)) != len(risks)
            or any(not risk for risk in risks)
        ):
            raise ValueError("Probe expected risks must be unique and ordered")

    def to_record(self) -> dict[str, object]:
        return {"terminal": self.terminal, "required_risks": list(self.required_risks)}

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ProbeFailureExpectation:
        _exact_keys(record, {"terminal", "required_risks"}, "Probe expectation")
        risks = record.get("required_risks")
        if not isinstance(risks, list) or any(
            not isinstance(value, str) for value in risks
        ):
            raise TypeError("Probe expected risks must be strings")
        return cls(
            _nonempty_string(record.get("terminal"), "expected terminal"),
            tuple(cast("list[str]", risks)),
        )


@dataclass(frozen=True, slots=True)
class NumericalProbeDefinition:
    """Immutable predeclared contract for one truth-backed reduced case."""

    probe_id: str
    category: ProbeCategory
    derivation: ProbeDerivation
    truth: ProbeTruthAuthority
    eligible_claims: tuple[str, ...]
    policy: CanonicalBaselineValue
    budget: CanonicalBaselineValue
    seed: CanonicalBaselineValue
    failure_expectations: tuple[ProbeFailureExpectation, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _nonempty_string(self.probe_id, "identifier")
        if self.category not in {
            "TRF_ROUTINE",
            "TRF_DIFFICULT",
            "GRID",
            "GRID_ORDERING",
            "DE_SEARCH",
            "FINITE_DIFFERENCE",
            "BOUNDS",
            "RANK",
            "CONDITIONING",
        }:
            raise ValueError("Unknown numerical probe category")
        claims = _canonical_strings(self.eligible_claims, "eligible claims")
        if (
            not isinstance(self.policy, CanonicalBaselineValue)
            or not isinstance(self.budget, CanonicalBaselineValue)
            or not isinstance(self.seed, CanonicalBaselineValue)
        ):
            raise TypeError("Probe policy, budget, and seed must be canonical values")
        expectations = tuple(self.failure_expectations)
        if not expectations or tuple(sorted(expectations)) != expectations:
            raise ValueError("Probe failure expectations must be non-empty and ordered")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-probe-definition",
                (
                    self.probe_id,
                    self.category,
                    self.derivation.identity,
                    self.truth.identity,
                    claims,
                    self.policy.encoded,
                    self.budget.encoded,
                    self.seed.encoded,
                    tuple(
                        (expectation.terminal, expectation.required_risks)
                        for expectation in expectations
                    ),
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "probe_id": self.probe_id,
            "category": self.category,
            "derivation": self.derivation.to_record(),
            "truth": self.truth.to_record(),
            "eligible_claims": list(self.eligible_claims),
            "policy": self.policy.to_record_value(),
            "budget": self.budget.to_record_value(),
            "seed": self.seed.to_record_value(),
            "failure_expectations": [
                expectation.to_record() for expectation in self.failure_expectations
            ],
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> NumericalProbeDefinition:
        _exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "probe_id",
                "category",
                "derivation",
                "truth",
                "eligible_claims",
                "policy",
                "budget",
                "seed",
                "failure_expectations",
                "identity",
            },
            "Numerical probe definition",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported numerical probe semantics")
        raw_derivation = record.get("derivation")
        raw_truth = record.get("truth")
        raw_claims = record.get("eligible_claims")
        raw_expectations = record.get("failure_expectations")
        if (
            not isinstance(raw_derivation, Mapping)
            or not isinstance(raw_truth, Mapping)
            or not isinstance(raw_claims, list)
            or not isinstance(raw_expectations, list)
            or any(not isinstance(value, str) for value in raw_claims)
            or any(not isinstance(value, Mapping) for value in raw_expectations)
        ):
            raise TypeError("Malformed numerical probe definition")
        definition = cls(
            _nonempty_string(record.get("probe_id"), "identifier"),
            cast("ProbeCategory", record.get("category")),
            ProbeDerivation.from_record(cast("Mapping[str, object]", raw_derivation)),
            ProbeTruthAuthority.from_record(cast("Mapping[str, object]", raw_truth)),
            _canonical_strings(cast("Sequence[str]", raw_claims), "eligible claims"),
            CanonicalBaselineValue.from_record(record.get("policy"), "probe policy"),
            CanonicalBaselineValue.from_record(record.get("budget"), "probe budget"),
            CanonicalBaselineValue.from_record(record.get("seed"), "probe seed"),
            tuple(
                ProbeFailureExpectation.from_record(
                    cast("Mapping[str, object]", expectation)
                )
                for expectation in raw_expectations
            ),
        )
        if record.get("identity") != definition.identity:
            raise ValueError("Numerical probe identity does not match its payload")
        return definition


@dataclass(frozen=True, slots=True)
class AuthoritativeObjectiveAccounting:
    """Authoritative evaluator requests, distinct from backend diagnostics."""

    workflow_requests: int
    materialization_requests: int
    fresh_evaluations: int
    cache_served_requests: int
    requests_refused: int
    failure_kinds: tuple[RequestFailureKind, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        values = (
            self.workflow_requests,
            self.materialization_requests,
            self.fresh_evaluations,
            self.cache_served_requests,
            self.requests_refused,
        )
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in values
        ):
            raise ValueError("Objective-request counters must be non-negative integers")
        if any(
            kind not in {"INVALID_INPUT", "EVALUATION_FAILURE", "NON_FINITE_RESULT"}
            for kind in self.failure_kinds
        ):
            raise ValueError("Unknown authoritative request failure kind")
        if (
            self.fresh_evaluations
            + self.cache_served_requests
            + self.requests_refused
            + self.requests_failed
            != self.requests_received
        ):
            raise ValueError("Objective requests do not reconcile with terminal counts")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "probe-authoritative-objective-accounting",
                (values, self.failure_kinds),
            ),
        )

    @property
    def requests_received(self) -> int:
        return self.workflow_requests + self.materialization_requests

    @property
    def requests_completed(self) -> int:
        return self.fresh_evaluations + self.cache_served_requests

    @property
    def requests_failed(self) -> int:
        return len(self.failure_kinds)

    @property
    def cache_hits(self) -> int:
        return self.cache_served_requests

    @property
    def cache_misses(self) -> int:
        return self.fresh_evaluations

    def to_record(self) -> dict[str, object]:
        return {
            "workflow_requests": self.workflow_requests,
            "materialization_requests": self.materialization_requests,
            "requests_received": self.requests_received,
            "requests_completed": self.requests_completed,
            "fresh_evaluations": self.fresh_evaluations,
            "cache_served_requests": self.cache_served_requests,
            "requests_refused": self.requests_refused,
            "requests_failed": self.requests_failed,
            "failure_kinds": list(self.failure_kinds),
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls, record: Mapping[str, object]
    ) -> AuthoritativeObjectiveAccounting:
        _exact_keys(
            record,
            {
                "workflow_requests",
                "materialization_requests",
                "requests_received",
                "requests_completed",
                "fresh_evaluations",
                "cache_served_requests",
                "requests_refused",
                "requests_failed",
                "failure_kinds",
                "identity",
            },
            "Objective request accounting",
        )
        raw_failures = record.get("failure_kinds")
        if not isinstance(raw_failures, list) or any(
            not isinstance(value, str) for value in raw_failures
        ):
            raise TypeError("Objective request failure kinds must be strings")
        accounting = cls(
            _record_nonnegative_int(
                record.get("workflow_requests"), "workflow requests"
            ),
            _record_nonnegative_int(
                record.get("materialization_requests"),
                "materialization requests",
            ),
            _record_nonnegative_int(record.get("fresh_evaluations"), "evaluations"),
            _record_nonnegative_int(
                record.get("cache_served_requests"), "cache-served requests"
            ),
            _record_nonnegative_int(record.get("requests_refused"), "refused requests"),
            tuple(cast("list[RequestFailureKind]", raw_failures)),
        )
        if (
            record.get("requests_received") != accounting.requests_received
            or record.get("requests_completed") != accounting.requests_completed
            or record.get("requests_failed") != accounting.requests_failed
            or record.get("identity") != accounting.identity
        ):
            raise ValueError("Objective request accounting does not match its payload")
        return accounting


# Compatibility name for the initial #591 draft; new code should use the explicit name.
ObjectiveRequestAccounting = AuthoritativeObjectiveAccounting


@dataclass(frozen=True, slots=True)
class BackendDiagnosticCounters:
    """Non-authoritative backend/callback counters retained for diagnosis."""

    nfev: int
    njev: int | None
    iterations: int | None
    callback_calls: int
    diagnostic_evaluations: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        values = (self.nfev, self.callback_calls, self.diagnostic_evaluations)
        if (
            any(
                isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in values
            )
            or (
                self.njev is not None
                and (
                    isinstance(self.njev, bool)
                    or not isinstance(self.njev, int)
                    or self.njev < 0
                )
            )
            or (
                self.iterations is not None
                and (
                    isinstance(self.iterations, bool)
                    or not isinstance(self.iterations, int)
                    or self.iterations < 0
                )
            )
        ):
            raise ValueError("Backend diagnostic counters must be non-negative")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "probe-backend-diagnostics",
                (
                    self.nfev,
                    self.njev,
                    self.iterations,
                    self.callback_calls,
                    self.diagnostic_evaluations,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "nfev": self.nfev,
            "njev": self.njev,
            "iterations": self.iterations,
            "callback_calls": self.callback_calls,
            "diagnostic_evaluations": self.diagnostic_evaluations,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> BackendDiagnosticCounters:
        _exact_keys(
            record,
            {
                "nfev",
                "njev",
                "iterations",
                "callback_calls",
                "diagnostic_evaluations",
                "identity",
            },
            "Backend diagnostics",
        )
        raw_njev = record.get("njev")
        if raw_njev is not None:
            raw_njev = _record_nonnegative_int(raw_njev, "backend njev")
        raw_iterations = record.get("iterations")
        if raw_iterations is not None:
            raw_iterations = _record_nonnegative_int(
                raw_iterations,
                "backend iterations",
            )
        diagnostics = cls(
            _record_nonnegative_int(record.get("nfev"), "backend nfev"),
            raw_njev,
            raw_iterations,
            _record_nonnegative_int(record.get("callback_calls"), "callback calls"),
            _record_nonnegative_int(
                record.get("diagnostic_evaluations"), "diagnostic evaluations"
            ),
        )
        if record.get("identity") != diagnostics.identity:
            raise ValueError("Backend diagnostics identity does not match its payload")
        return diagnostics


@dataclass(frozen=True, slots=True)
class SolverProbeEvidence:
    """Typed bounded least-squares result with detached accounting evidence."""

    start: tuple[float, ...]
    accepted: tuple[float, ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    residuals: tuple[float, ...]
    chi_square: float
    active_mask: tuple[int, ...]
    objective_accounting: ObjectiveRequestAccounting
    backend_diagnostics: BackendDiagnosticCounters
    trajectory_fingerprint: str
    candidate_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        dimension = len(self.start)
        if (
            dimension < 1
            or any(
                len(values) != dimension
                for values in (self.accepted, self.lower_bounds, self.upper_bounds)
            )
            or len(self.active_mask) != dimension
        ):
            raise ValueError("Solver probe vectors have inconsistent dimensions")
        for values in (
            self.start,
            self.accepted,
            self.lower_bounds,
            self.upper_bounds,
            self.residuals,
        ):
            _vector_tokens(values)
        chi_square = float(self.chi_square)
        if not math.isfinite(chi_square) or chi_square < 0.0:
            raise ValueError("Solver probe chi square must be finite and non-negative")
        if chi_square != math.fsum(value * value for value in self.residuals):
            raise ValueError("Solver probe chi square does not match residuals")
        if any(mask not in {-1, 0, 1} for mask in self.active_mask):
            raise ValueError("Solver active-mask entries must be -1, 0, or 1")
        for fingerprint in (
            self.trajectory_fingerprint,
            self.candidate_fingerprint,
        ):
            if len(fingerprint) != 64:
                raise ValueError("Solver probe fingerprints must be SHA-256 digests")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "solver-probe-evidence",
                (
                    _vector_tokens(self.start),
                    _vector_tokens(self.accepted),
                    _vector_tokens(self.lower_bounds),
                    _vector_tokens(self.upper_bounds),
                    _vector_tokens(self.residuals),
                    _float_token(chi_square),
                    self.active_mask,
                    self.objective_accounting.identity,
                    self.backend_diagnostics.identity,
                    self.trajectory_fingerprint,
                    self.candidate_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_accounting(self) -> AuthoritativeObjectiveAccounting:
        return self.objective_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "SOLVER",
            "start": list(_vector_tokens(self.start)),
            "accepted": list(_vector_tokens(self.accepted)),
            "lower_bounds": list(_vector_tokens(self.lower_bounds)),
            "upper_bounds": list(_vector_tokens(self.upper_bounds)),
            "residuals": list(_vector_tokens(self.residuals)),
            "chi_square": _float_token(self.chi_square),
            "active_mask": list(self.active_mask),
            "objective_accounting": self.objective_accounting.to_record(),
            "backend_diagnostics": self.backend_diagnostics.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "candidate_fingerprint": self.candidate_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> SolverProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "start",
                "accepted",
                "lower_bounds",
                "upper_bounds",
                "residuals",
                "chi_square",
                "active_mask",
                "objective_accounting",
                "backend_diagnostics",
                "trajectory_fingerprint",
                "candidate_fingerprint",
                "identity",
            },
            "Solver probe evidence",
        )
        if record.get("kind") != "SOLVER":
            raise ValueError("Unknown solver probe evidence kind")
        raw_mask = record.get("active_mask")
        raw_accounting = record.get("objective_accounting")
        raw_diagnostics = record.get("backend_diagnostics")
        if (
            not isinstance(raw_mask, list)
            or any(not isinstance(value, int) for value in raw_mask)
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
        ):
            raise TypeError("Malformed solver probe evidence")
        evidence = cls(
            _vector_from_record(record.get("start"), "start"),
            _vector_from_record(record.get("accepted"), "accepted"),
            _vector_from_record(record.get("lower_bounds"), "lower bounds"),
            _vector_from_record(record.get("upper_bounds"), "upper bounds"),
            _vector_from_record(record.get("residuals"), "residuals"),
            _float_from_record(record.get("chi_square"), "chi square"),
            tuple(cast("list[int]", raw_mask)),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            _nonempty_string(
                record.get("trajectory_fingerprint"), "trajectory fingerprint"
            ),
            _nonempty_string(
                record.get("candidate_fingerprint"), "candidate fingerprint"
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("Solver probe evidence identity does not match payload")
        return evidence


@dataclass(frozen=True, slots=True)
class GridSeedEvidence:
    """One physical-coordinate seed and its bounded local-solver evidence."""

    ordinal: int
    seed: tuple[float, ...]
    terminal: str
    solver: SolverProbeEvidence
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _record_nonnegative_int(self.ordinal, "GRID seed ordinal")
        _nonempty_string(self.terminal, "GRID seed terminal")
        if len(self.seed) != len(self.solver.start) or self.seed != self.solver.start:
            raise ValueError("GRID seed does not match its solver start")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "grid-seed-probe-evidence",
                (
                    self.ordinal,
                    _vector_tokens(self.seed),
                    self.terminal,
                    self.solver.identity,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "ordinal": self.ordinal,
            "seed": list(_vector_tokens(self.seed)),
            "terminal": self.terminal,
            "solver": self.solver.to_record(),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> GridSeedEvidence:
        _exact_keys(
            record,
            {"ordinal", "seed", "terminal", "solver", "identity"},
            "GRID seed",
        )
        raw_solver = record.get("solver")
        if not isinstance(raw_solver, Mapping):
            raise TypeError("GRID seed solver evidence must be a record")
        evidence = cls(
            _record_nonnegative_int(record.get("ordinal"), "GRID seed ordinal"),
            _vector_from_record(record.get("seed"), "GRID seed"),
            _nonempty_string(record.get("terminal"), "GRID seed terminal"),
            SolverProbeEvidence.from_record(cast("Mapping[str, object]", raw_solver)),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("GRID seed identity does not match its payload")
        return evidence


@dataclass(frozen=True, slots=True)
class GridProbeEvidence:
    """Typed immutable 27-seed execution and candidate-ordering evidence."""

    seeds: tuple[GridSeedEvidence, ...]
    ordered_seed_ordinals: tuple[int, ...]
    selected_seed_ordinal: int
    objective_accounting: ObjectiveRequestAccounting
    backend_diagnostics: BackendDiagnosticCounters
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if tuple(seed.ordinal for seed in self.seeds) != tuple(range(len(self.seeds))):
            raise ValueError("GRID probe seeds must be in canonical ordinal order")
        expected_order = tuple(
            seed.ordinal
            for seed in sorted(
                self.seeds,
                key=lambda seed: (
                    seed.solver.chi_square,
                    seed.solver.accepted,
                    seed.ordinal,
                ),
            )
        )
        if self.ordered_seed_ordinals != expected_order:
            raise ValueError("GRID candidate ordering does not match evidence")
        if not expected_order or self.selected_seed_ordinal != expected_order[0]:
            raise ValueError("GRID selected seed is not the first ordered candidate")
        if len(self.trajectory_fingerprint) != 64:
            raise ValueError("GRID trajectory fingerprint must be a SHA-256 digest")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "grid-probe-evidence",
                (
                    tuple(seed.identity for seed in self.seeds),
                    self.ordered_seed_ordinals,
                    self.selected_seed_ordinal,
                    self.objective_accounting.identity,
                    self.backend_diagnostics.identity,
                    self.trajectory_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_accounting(self) -> AuthoritativeObjectiveAccounting:
        return self.objective_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "GRID",
            "seeds": [seed.to_record() for seed in self.seeds],
            "ordered_seed_ordinals": list(self.ordered_seed_ordinals),
            "selected_seed_ordinal": self.selected_seed_ordinal,
            "objective_accounting": self.objective_accounting.to_record(),
            "backend_diagnostics": self.backend_diagnostics.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> GridProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "seeds",
                "ordered_seed_ordinals",
                "selected_seed_ordinal",
                "objective_accounting",
                "backend_diagnostics",
                "trajectory_fingerprint",
                "identity",
            },
            "GRID probe evidence",
        )
        raw_seeds = record.get("seeds")
        raw_order = record.get("ordered_seed_ordinals")
        raw_accounting = record.get("objective_accounting")
        raw_diagnostics = record.get("backend_diagnostics")
        if (
            record.get("kind") != "GRID"
            or not isinstance(raw_seeds, list)
            or not isinstance(raw_order, list)
            or any(not isinstance(value, Mapping) for value in raw_seeds)
            or any(not isinstance(value, int) for value in raw_order)
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
        ):
            raise TypeError("Malformed GRID probe evidence")
        evidence = cls(
            tuple(
                GridSeedEvidence.from_record(cast("Mapping[str, object]", seed))
                for seed in raw_seeds
            ),
            tuple(cast("list[int]", raw_order)),
            _record_nonnegative_int(
                record.get("selected_seed_ordinal"), "selected GRID seed"
            ),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            _nonempty_string(
                record.get("trajectory_fingerprint"), "GRID trajectory fingerprint"
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("GRID evidence identity does not match its payload")
        return evidence


@dataclass(frozen=True, slots=True, order=True)
class GridOrderingCandidate:
    """One analytically declared candidate for the GRID ordering contract."""

    ordinal: int
    final_vector: tuple[float, ...]
    chi_square: float
    identity: str = field(init=False, compare=False)

    def __post_init__(self) -> None:
        _record_nonnegative_int(self.ordinal, "GRID ordering ordinal")
        _vector_tokens(self.final_vector)
        if not math.isfinite(self.chi_square) or self.chi_square < 0.0:
            raise ValueError("GRID ordering chi square must be finite and non-negative")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "grid-ordering-candidate",
                (
                    self.ordinal,
                    _vector_tokens(self.final_vector),
                    _float_token(self.chi_square),
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "ordinal": self.ordinal,
            "final_vector": list(_vector_tokens(self.final_vector)),
            "chi_square": _float_token(self.chi_square),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> GridOrderingCandidate:
        _exact_keys(
            record,
            {"ordinal", "final_vector", "chi_square", "identity"},
            "GRID ordering candidate",
        )
        candidate = cls(
            _record_nonnegative_int(record.get("ordinal"), "GRID ordering ordinal"),
            _vector_from_record(record.get("final_vector"), "GRID final vector"),
            _float_from_record(record.get("chi_square"), "GRID chi square"),
        )
        if record.get("identity") != candidate.identity:
            raise ValueError("GRID ordering candidate identity does not match payload")
        return candidate


@dataclass(frozen=True, slots=True)
class GridOrderingProbeEvidence:
    """Observed ordering retained separately from its predeclared expected truth."""

    candidates: tuple[GridOrderingCandidate, ...]
    observed_order: tuple[int, ...]
    expected_order: tuple[int, ...]
    selected_ordinal: int
    objective_accounting: ObjectiveRequestAccounting
    backend_diagnostics: BackendDiagnosticCounters
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        ordinals = tuple(candidate.ordinal for candidate in self.candidates)
        if not ordinals or ordinals != tuple(range(len(ordinals))):
            raise ValueError("GRID ordering candidates must retain ordinal order")
        if sorted(self.observed_order) != list(ordinals):
            raise ValueError("Observed GRID ordering must be a candidate permutation")
        if sorted(self.expected_order) != list(ordinals):
            raise ValueError("Expected GRID ordering must be a candidate permutation")
        if self.selected_ordinal != self.observed_order[0]:
            raise ValueError(
                "GRID ordering selection must match observed first candidate"
            )
        if len(self.trajectory_fingerprint) != 64:
            raise ValueError("GRID ordering trajectory must be a SHA-256 digest")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "grid-ordering-probe-evidence",
                (
                    tuple(candidate.identity for candidate in self.candidates),
                    self.observed_order,
                    self.expected_order,
                    self.selected_ordinal,
                    self.objective_accounting.identity,
                    self.backend_diagnostics.identity,
                    self.trajectory_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_accounting(self) -> AuthoritativeObjectiveAccounting:
        return self.objective_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "GRID_ORDERING",
            "candidates": [candidate.to_record() for candidate in self.candidates],
            "observed_order": list(self.observed_order),
            "expected_order": list(self.expected_order),
            "selected_ordinal": self.selected_ordinal,
            "objective_accounting": self.objective_accounting.to_record(),
            "backend_diagnostics": self.backend_diagnostics.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> GridOrderingProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "candidates",
                "observed_order",
                "expected_order",
                "selected_ordinal",
                "objective_accounting",
                "backend_diagnostics",
                "trajectory_fingerprint",
                "identity",
            },
            "GRID ordering probe evidence",
        )
        raw_candidates = record.get("candidates")
        raw_observed = record.get("observed_order")
        raw_expected = record.get("expected_order")
        raw_accounting = record.get("objective_accounting")
        raw_diagnostics = record.get("backend_diagnostics")
        if (
            record.get("kind") != "GRID_ORDERING"
            or not isinstance(raw_candidates, list)
            or any(not isinstance(value, Mapping) for value in raw_candidates)
            or not isinstance(raw_observed, list)
            or any(not isinstance(value, int) for value in raw_observed)
            or not isinstance(raw_expected, list)
            or any(not isinstance(value, int) for value in raw_expected)
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
        ):
            raise TypeError("Malformed GRID ordering probe evidence")
        evidence = cls(
            tuple(
                GridOrderingCandidate.from_record(
                    cast("Mapping[str, object]", candidate)
                )
                for candidate in raw_candidates
            ),
            tuple(cast("list[int]", raw_observed)),
            tuple(cast("list[int]", raw_expected)),
            _record_nonnegative_int(
                record.get("selected_ordinal"), "selected GRID ordering candidate"
            ),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            _nonempty_string(
                record.get("trajectory_fingerprint"),
                "GRID ordering trajectory fingerprint",
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("GRID ordering evidence identity does not match payload")
        return evidence


@dataclass(frozen=True, slots=True, order=True)
class DePopulationCandidate:
    """One terminal DE population member in stable backend index order."""

    population_index: int
    vector: tuple[float, ...]
    objective: float
    backend_objective: float
    fingerprint: str = field(init=False, compare=False)

    def __post_init__(self) -> None:
        _record_nonnegative_int(self.population_index, "DE population index")
        _vector_tokens(self.vector)
        objective = float(self.objective)
        if not math.isfinite(objective) or objective < 0.0:
            raise ValueError("DE population objective must be finite and non-negative")
        backend_objective = float(self.backend_objective)
        if not math.isfinite(backend_objective) or backend_objective < 0.0:
            raise ValueError(
                "DE backend population objective must be finite and non-negative"
            )
        object.__setattr__(
            self,
            "fingerprint",
            _identity(
                "de-population-candidate",
                (
                    self.population_index,
                    _vector_tokens(self.vector),
                    _float_token(objective),
                    _float_token(backend_objective),
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "population_index": self.population_index,
            "vector": list(_vector_tokens(self.vector)),
            "objective": _float_token(self.objective),
            "backend_objective": _float_token(self.backend_objective),
            "fingerprint": self.fingerprint,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> DePopulationCandidate:
        _exact_keys(
            record,
            {
                "population_index",
                "vector",
                "objective",
                "backend_objective",
                "fingerprint",
            },
            "DE population candidate",
        )
        candidate = cls(
            _record_nonnegative_int(
                record.get("population_index"), "DE population index"
            ),
            _vector_from_record(record.get("vector"), "DE candidate vector"),
            _float_from_record(record.get("objective"), "DE candidate objective"),
            _float_from_record(
                record.get("backend_objective"),
                "DE backend candidate objective",
            ),
        )
        if record.get("fingerprint") != candidate.fingerprint:
            raise ValueError("DE candidate fingerprint does not match its payload")
        return candidate


@dataclass(frozen=True, slots=True)
class DeProbeEvidence:
    """Seeded DE population ordering plus independently accounted TRF polish."""

    root_seed: int
    final_population: tuple[DePopulationCandidate, ...]
    ordered_population_indices: tuple[int, ...]
    selected_population_index: int
    search_accounting: ObjectiveRequestAccounting
    search_diagnostics: BackendDiagnosticCounters
    polish: SolverProbeEvidence
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.root_seed, bool)
            or not isinstance(self.root_seed, int)
            or not 0 <= self.root_seed < 2**64
        ):
            raise ValueError("DE root seed must be an unsigned 64-bit integer")
        if tuple(
            candidate.population_index for candidate in self.final_population
        ) != tuple(range(len(self.final_population))):
            raise ValueError("DE population must retain backend index order")
        expected_order = tuple(
            candidate.population_index
            for candidate in sorted(
                self.final_population,
                key=lambda candidate: (candidate.objective, candidate.population_index),
            )
        )
        if self.ordered_population_indices != expected_order:
            raise ValueError("DE population ordering does not match candidate evidence")
        if not expected_order or self.selected_population_index != expected_order[0]:
            raise ValueError("DE selection is not the first ordered population member")
        if len(self.trajectory_fingerprint) != 64:
            raise ValueError("DE trajectory fingerprint must be a SHA-256 digest")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "de-probe-evidence",
                (
                    self.root_seed,
                    tuple(candidate.fingerprint for candidate in self.final_population),
                    self.ordered_population_indices,
                    self.selected_population_index,
                    self.search_accounting.identity,
                    self.search_diagnostics.identity,
                    self.polish.identity,
                    self.trajectory_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_search_accounting(
        self,
    ) -> AuthoritativeObjectiveAccounting:
        return self.search_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "DE_SEARCH",
            "root_seed": self.root_seed,
            "final_population": [
                candidate.to_record() for candidate in self.final_population
            ],
            "ordered_population_indices": list(self.ordered_population_indices),
            "selected_population_index": self.selected_population_index,
            "search_accounting": self.search_accounting.to_record(),
            "search_diagnostics": self.search_diagnostics.to_record(),
            "polish": self.polish.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> DeProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "root_seed",
                "final_population",
                "ordered_population_indices",
                "selected_population_index",
                "search_accounting",
                "search_diagnostics",
                "polish",
                "trajectory_fingerprint",
                "identity",
            },
            "DE probe evidence",
        )
        raw_population = record.get("final_population")
        raw_order = record.get("ordered_population_indices")
        raw_accounting = record.get("search_accounting")
        raw_diagnostics = record.get("search_diagnostics")
        raw_polish = record.get("polish")
        if (
            record.get("kind") != "DE_SEARCH"
            or not isinstance(raw_population, list)
            or not isinstance(raw_order, list)
            or any(not isinstance(value, Mapping) for value in raw_population)
            or any(not isinstance(value, int) for value in raw_order)
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
            or not isinstance(raw_polish, Mapping)
        ):
            raise TypeError("Malformed DE probe evidence")
        evidence = cls(
            _record_nonnegative_int(record.get("root_seed"), "DE root seed"),
            tuple(
                DePopulationCandidate.from_record(
                    cast("Mapping[str, object]", candidate)
                )
                for candidate in raw_population
            ),
            tuple(cast("list[int]", raw_order)),
            _record_nonnegative_int(
                record.get("selected_population_index"), "selected DE candidate"
            ),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            SolverProbeEvidence.from_record(cast("Mapping[str, object]", raw_polish)),
            _nonempty_string(
                record.get("trajectory_fingerprint"), "DE trajectory fingerprint"
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("DE evidence identity does not match its payload")
        return evidence


@dataclass(frozen=True, slots=True)
class FiniteDifferenceProbeEvidence:
    """Chosen finite-difference estimate with exact steps and request evidence."""

    point: float
    actual_steps: tuple[float, ...]
    weights: tuple[float, ...]
    sampled_values: tuple[float, ...]
    truth: float
    estimate: float
    absolute_error: float
    absolute_tolerance: float
    objective_accounting: ObjectiveRequestAccounting
    backend_diagnostics: BackendDiagnosticCounters
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _float_token(self.point)
        if (
            len(self.actual_steps) != 5
            or len(self.weights) != 5
            or len(self.sampled_values) != 5
        ):
            raise ValueError("Finite-difference evidence requires five stencil points")
        _vector_tokens(self.actual_steps)
        _vector_tokens(self.weights)
        _vector_tokens(self.sampled_values)
        for value in (
            self.truth,
            self.estimate,
            self.absolute_error,
            self.absolute_tolerance,
        ):
            _float_token(value)
        if self.absolute_error != abs(self.estimate - self.truth):
            raise ValueError(
                "Finite-difference error does not match estimate and truth"
            )
        reconstructed = math.fsum(
            weight * value
            for weight, value in zip(self.weights, self.sampled_values, strict=True)
        )
        if self.estimate != reconstructed:
            raise ValueError("Finite-difference estimate does not match its weights")
        if self.absolute_tolerance <= 0.0:
            raise ValueError("Finite-difference tolerance must be positive")
        if len(self.trajectory_fingerprint) != 64:
            raise ValueError("Stencil trajectory fingerprint must be a SHA-256 digest")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "finite-difference-probe-evidence",
                (
                    _float_token(self.point),
                    _vector_tokens(self.actual_steps),
                    _vector_tokens(self.weights),
                    _vector_tokens(self.sampled_values),
                    _float_token(self.truth),
                    _float_token(self.estimate),
                    _float_token(self.absolute_error),
                    _float_token(self.absolute_tolerance),
                    self.objective_accounting.identity,
                    self.backend_diagnostics.identity,
                    self.trajectory_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_accounting(self) -> AuthoritativeObjectiveAccounting:
        return self.objective_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "FINITE_DIFFERENCE",
            "point": _float_token(self.point),
            "actual_steps": list(_vector_tokens(self.actual_steps)),
            "weights": list(_vector_tokens(self.weights)),
            "sampled_values": list(_vector_tokens(self.sampled_values)),
            "truth": _float_token(self.truth),
            "estimate": _float_token(self.estimate),
            "absolute_error": _float_token(self.absolute_error),
            "absolute_tolerance": _float_token(self.absolute_tolerance),
            "objective_accounting": self.objective_accounting.to_record(),
            "backend_diagnostics": self.backend_diagnostics.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> FiniteDifferenceProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "point",
                "actual_steps",
                "weights",
                "sampled_values",
                "truth",
                "estimate",
                "absolute_error",
                "absolute_tolerance",
                "objective_accounting",
                "backend_diagnostics",
                "trajectory_fingerprint",
                "identity",
            },
            "Finite-difference probe evidence",
        )
        raw_accounting = record.get("objective_accounting")
        raw_diagnostics = record.get("backend_diagnostics")
        if (
            record.get("kind") != "FINITE_DIFFERENCE"
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
        ):
            raise TypeError("Malformed finite-difference probe evidence")
        evidence = cls(
            _float_from_record(record.get("point"), "finite-difference point"),
            _vector_from_record(record.get("actual_steps"), "actual steps"),
            _vector_from_record(record.get("weights"), "stencil weights"),
            _vector_from_record(record.get("sampled_values"), "sampled values"),
            _float_from_record(record.get("truth"), "derivative truth"),
            _float_from_record(record.get("estimate"), "derivative estimate"),
            _float_from_record(record.get("absolute_error"), "absolute error"),
            _float_from_record(record.get("absolute_tolerance"), "absolute tolerance"),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            _nonempty_string(
                record.get("trajectory_fingerprint"), "stencil trajectory fingerprint"
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError(
                "Finite-difference evidence identity does not match payload"
            )
        return evidence


@dataclass(frozen=True, slots=True)
class SpectralRiskProbeEvidence:
    """Exact small-matrix rank and conditioning risk evidence."""

    jacobian: tuple[tuple[float, ...], ...]
    singular_values: tuple[float, ...]
    rank: int
    dimension: int
    condition: float | None
    condition_limit: float
    risks: tuple[str, ...]
    objective_accounting: ObjectiveRequestAccounting
    backend_diagnostics: BackendDiagnosticCounters
    trajectory_fingerprint: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.jacobian or any(
            len(row) != self.dimension for row in self.jacobian
        ):
            raise ValueError("Spectral probe Jacobian has the wrong dimension")
        for row in self.jacobian:
            _vector_tokens(row)
        if len(self.singular_values) != self.dimension:
            raise ValueError(
                "Spectral probe must retain the complete singular spectrum"
            )
        _vector_tokens(self.singular_values)
        if not 0 <= self.rank <= self.dimension:
            raise ValueError("Spectral rank is outside its matrix dimension")
        if self.condition is not None:
            _float_token(self.condition)
            if self.condition < 1.0:
                raise ValueError("Spectral condition must be at least one")
        _float_token(self.condition_limit)
        if self.condition_limit <= 0.0:
            raise ValueError("Conditioning limit must be positive")
        _canonical_strings(self.risks, "spectral risks")
        if len(self.trajectory_fingerprint) != 64:
            raise ValueError("Spectral trajectory fingerprint must be a SHA-256 digest")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "spectral-risk-probe-evidence",
                (
                    tuple(_vector_tokens(row) for row in self.jacobian),
                    _vector_tokens(self.singular_values),
                    self.rank,
                    self.dimension,
                    None if self.condition is None else _float_token(self.condition),
                    _float_token(self.condition_limit),
                    self.risks,
                    self.objective_accounting.identity,
                    self.backend_diagnostics.identity,
                    self.trajectory_fingerprint,
                ),
            ),
        )

    @property
    def authoritative_accounting(self) -> AuthoritativeObjectiveAccounting:
        return self.objective_accounting

    def to_record(self) -> dict[str, object]:
        return {
            "kind": "SPECTRAL_RISK",
            "jacobian": [list(_vector_tokens(row)) for row in self.jacobian],
            "singular_values": list(_vector_tokens(self.singular_values)),
            "rank": self.rank,
            "dimension": self.dimension,
            "condition": (
                None if self.condition is None else _float_token(self.condition)
            ),
            "condition_limit": _float_token(self.condition_limit),
            "risks": list(self.risks),
            "objective_accounting": self.objective_accounting.to_record(),
            "backend_diagnostics": self.backend_diagnostics.to_record(),
            "trajectory_fingerprint": self.trajectory_fingerprint,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> SpectralRiskProbeEvidence:
        _exact_keys(
            record,
            {
                "kind",
                "jacobian",
                "singular_values",
                "rank",
                "dimension",
                "condition",
                "condition_limit",
                "risks",
                "objective_accounting",
                "backend_diagnostics",
                "trajectory_fingerprint",
                "identity",
            },
            "Spectral risk probe evidence",
        )
        raw_jacobian = record.get("jacobian")
        raw_risks = record.get("risks")
        raw_accounting = record.get("objective_accounting")
        raw_diagnostics = record.get("backend_diagnostics")
        if (
            record.get("kind") != "SPECTRAL_RISK"
            or not isinstance(raw_jacobian, list)
            or any(not isinstance(row, list) for row in raw_jacobian)
            or not isinstance(raw_risks, list)
            or any(not isinstance(risk, str) for risk in raw_risks)
            or not isinstance(raw_accounting, Mapping)
            or not isinstance(raw_diagnostics, Mapping)
        ):
            raise TypeError("Malformed spectral risk probe evidence")
        raw_condition = record.get("condition")
        evidence = cls(
            tuple(
                _vector_from_record(row, f"Jacobian row {index}")
                for index, row in enumerate(raw_jacobian)
            ),
            _vector_from_record(record.get("singular_values"), "singular values"),
            _record_nonnegative_int(record.get("rank"), "rank"),
            _record_nonnegative_int(record.get("dimension"), "dimension"),
            (
                None
                if raw_condition is None
                else _float_from_record(raw_condition, "condition")
            ),
            _float_from_record(record.get("condition_limit"), "condition limit"),
            _canonical_strings(cast("Sequence[str]", raw_risks), "spectral risks"),
            ObjectiveRequestAccounting.from_record(
                cast("Mapping[str, object]", raw_accounting)
            ),
            BackendDiagnosticCounters.from_record(
                cast("Mapping[str, object]", raw_diagnostics)
            ),
            _nonempty_string(
                record.get("trajectory_fingerprint"), "spectral fingerprint"
            ),
        )
        if record.get("identity") != evidence.identity:
            raise ValueError("Spectral evidence identity does not match its payload")
        return evidence


type NumericalProbeEvidence = (
    SolverProbeEvidence
    | GridProbeEvidence
    | GridOrderingProbeEvidence
    | DeProbeEvidence
    | FiniteDifferenceProbeEvidence
    | SpectralRiskProbeEvidence
)


@dataclass(frozen=True, slots=True)
class NumericalProbeArtifact:
    """Content-addressed result of one predeclared probe definition."""

    definition_identity: str
    probe_id: str
    terminal: str
    risks: tuple[str, ...]
    satisfied_claims: tuple[str, ...]
    evidence: NumericalProbeEvidence
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for value, name in (
            (self.definition_identity, "definition identity"),
            (self.probe_id, "identifier"),
            (self.terminal, "terminal"),
        ):
            _nonempty_string(value, name)
        if len(self.definition_identity) != 64:
            raise ValueError("Probe definition identity must be a SHA-256 digest")
        if self.risks and (
            tuple(sorted(self.risks)) != self.risks
            or len(set(self.risks)) != len(self.risks)
        ):
            raise ValueError("Probe risks must be unique and ordered")
        if self.satisfied_claims and (
            tuple(sorted(self.satisfied_claims)) != self.satisfied_claims
            or len(set(self.satisfied_claims)) != len(self.satisfied_claims)
            or any(not claim for claim in self.satisfied_claims)
        ):
            raise ValueError("Probe satisfied claims must be unique and ordered")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-probe-artifact",
                (
                    self.definition_identity,
                    self.probe_id,
                    self.terminal,
                    self.risks,
                    self.satisfied_claims,
                    self.evidence.identity,
                ),
            ),
        )

    def validate_definition(self, definition: NumericalProbeDefinition) -> None:
        if (
            self.definition_identity != definition.identity
            or self.probe_id != definition.probe_id
            or not set(self.satisfied_claims).issubset(definition.eligible_claims)
        ):
            raise ValueError("Numerical probe artifact violates its definition")
        expected_evidence: type[NumericalProbeEvidence]
        if definition.category in {"TRF_ROUTINE", "TRF_DIFFICULT", "BOUNDS"}:
            expected_evidence = SolverProbeEvidence
        elif definition.category == "GRID":
            expected_evidence = GridProbeEvidence
        elif definition.category == "GRID_ORDERING":
            expected_evidence = GridOrderingProbeEvidence
        elif definition.category == "DE_SEARCH":
            expected_evidence = DeProbeEvidence
        elif definition.category == "FINITE_DIFFERENCE":
            expected_evidence = FiniteDifferenceProbeEvidence
        else:
            expected_evidence = SpectralRiskProbeEvidence
        if not isinstance(self.evidence, expected_evidence):
            raise TypeError("Numerical probe evidence does not match its category")

    def require_qualification(self, definition: NumericalProbeDefinition) -> None:
        self.validate_definition(definition)
        expectation = definition.failure_expectations[0]
        if (
            self.terminal != expectation.terminal
            or not set(expectation.required_risks).issubset(self.risks)
            or self.satisfied_claims != definition.eligible_claims
        ):
            raise ValueError("Numerical probe artifact does not qualify its definition")
        if (
            definition.category == "GRID_ORDERING"
            and isinstance(self.evidence, GridOrderingProbeEvidence)
            and self.evidence.observed_order != self.evidence.expected_order
        ):
            raise ValueError("GRID ordering evidence differs from independent truth")

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "definition_identity": self.definition_identity,
            "probe_id": self.probe_id,
            "terminal": self.terminal,
            "risks": list(self.risks),
            "satisfied_claims": list(self.satisfied_claims),
            "evidence": self.evidence.to_record(),
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        definition: NumericalProbeDefinition,
    ) -> NumericalProbeArtifact:
        _exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "definition_identity",
                "probe_id",
                "terminal",
                "risks",
                "satisfied_claims",
                "evidence",
                "identity",
            },
            "Numerical probe artifact",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported numerical probe artifact semantics")
        raw_risks = record.get("risks")
        raw_claims = record.get("satisfied_claims")
        raw_evidence = record.get("evidence")
        if (
            not isinstance(raw_risks, list)
            or not isinstance(raw_claims, list)
            or any(not isinstance(value, str) for value in (*raw_risks, *raw_claims))
            or not isinstance(raw_evidence, Mapping)
        ):
            raise TypeError("Malformed numerical probe artifact")
        evidence_kind = raw_evidence.get("kind")
        if evidence_kind == "SOLVER":
            evidence: NumericalProbeEvidence = SolverProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        elif evidence_kind == "GRID":
            evidence = GridProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        elif evidence_kind == "GRID_ORDERING":
            evidence = GridOrderingProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        elif evidence_kind == "DE_SEARCH":
            evidence = DeProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        elif evidence_kind == "FINITE_DIFFERENCE":
            evidence = FiniteDifferenceProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        elif evidence_kind == "SPECTRAL_RISK":
            evidence = SpectralRiskProbeEvidence.from_record(
                cast("Mapping[str, object]", raw_evidence)
            )
        else:
            raise ValueError("Unknown numerical probe evidence kind")
        artifact = cls(
            _nonempty_string(record.get("definition_identity"), "definition identity"),
            _nonempty_string(record.get("probe_id"), "identifier"),
            _nonempty_string(record.get("terminal"), "terminal"),
            tuple(cast("list[str]", raw_risks)),
            tuple(cast("list[str]", raw_claims)),
            evidence,
        )
        artifact.validate_definition(definition)
        if record.get("identity") != artifact.identity:
            raise ValueError("Numerical probe artifact identity does not match payload")
        return artifact


class NumericalProbeQualificationError(RuntimeError):
    """Typed rejection of evidence at the live lane-qualification seam."""

    def __init__(self, terminal: str) -> None:
        self.terminal = terminal
        super().__init__(f"Numerical probe qualification rejected: {terminal}")


@dataclass(frozen=True, slots=True)
class NumericalProbeBaseline:
    """Serializable observations; live lane authority is always supplied separately."""

    definitions: tuple[NumericalProbeDefinition, ...]
    artifacts: tuple[NumericalProbeArtifact, ...]
    observed_lane_identity: str | None = None
    observed_lane_role: LaneRole | None = None
    observed_attestation_identity: str | None = None
    observed_environment_identity: str | None = None
    reference_manifest_identity: str | None = None
    historical_qualification: Literal["CAPTURE_ONLY", "REFERENCE_MATCHED"] = (
        "CAPTURE_ONLY"
    )
    manifest_identity: str = field(init=False)
    identity: str = field(init=False)

    def _validate_historical_observation(self, manifest: str) -> None:
        observation = (
            self.observed_lane_identity,
            self.observed_lane_role,
            self.observed_attestation_identity,
            self.observed_environment_identity,
            self.reference_manifest_identity,
        )
        if self.historical_qualification == "CAPTURE_ONLY":
            if any(value is not None for value in observation):
                raise ValueError(
                    "Capture-only numerical evidence cannot claim lane qualification"
                )
            return
        if self.historical_qualification != "REFERENCE_MATCHED":
            raise ValueError("Unknown numerical probe historical qualification")
        for value, name in (
            (self.observed_lane_identity, "observed lane identity"),
            (self.observed_attestation_identity, "observed attestation identity"),
            (self.observed_environment_identity, "observed environment identity"),
            (self.reference_manifest_identity, "reference manifest identity"),
        ):
            if not isinstance(value, str) or len(value) != 64:
                raise ValueError(f"Matched numerical evidence requires {name}")
        if self.observed_lane_role not in {
            "CANONICAL_NUMERICAL",
            "PYTHON_COMPATIBILITY",
        }:
            raise ValueError("Matched numerical evidence requires a lane role")
        if self.reference_manifest_identity != manifest:
            raise ValueError("Historical reference does not match probe manifest")

    def __post_init__(self) -> None:
        definitions = tuple(self.definitions)
        artifacts = tuple(self.artifacts)
        if (
            not definitions
            or len(definitions) != len(artifacts)
            or len({definition.probe_id for definition in definitions})
            != len(definitions)
        ):
            raise ValueError(
                "Probe baseline requires one artifact per unique definition"
            )
        if definitions != numerical_probe_definitions():
            raise ValueError(
                "Probe baseline does not qualify the canonical probe catalogue"
            )
        for definition, artifact in zip(definitions, artifacts, strict=True):
            artifact.require_qualification(definition)
        manifest = _identity(
            "numerical-probe-manifest",
            tuple(
                (definition.probe_id, definition.identity, artifact.identity)
                for definition, artifact in zip(
                    definitions,
                    artifacts,
                    strict=True,
                )
            ),
        )
        observation = (
            self.observed_lane_identity,
            self.observed_lane_role,
            self.observed_attestation_identity,
            self.observed_environment_identity,
            self.reference_manifest_identity,
        )
        self._validate_historical_observation(manifest)
        object.__setattr__(self, "manifest_identity", manifest)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-probe-baseline",
                (
                    manifest,
                    observation,
                    self.historical_qualification,
                    tuple(definition.identity for definition in definitions),
                    tuple(artifact.identity for artifact in artifacts),
                ),
            ),
        )

    @property
    def historically_satisfied_claims(self) -> tuple[str, ...]:
        """Claims recorded at capture time, never a replacement for live authority."""

        if self.historical_qualification == "REFERENCE_MATCHED":
            return ("trajectory-manifest-compatible",)
        return ()

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "definitions": [definition.to_record() for definition in self.definitions],
            "artifacts": [artifact.to_record() for artifact in self.artifacts],
            "observed_lane_identity": self.observed_lane_identity,
            "observed_lane_role": self.observed_lane_role,
            "observed_attestation_identity": self.observed_attestation_identity,
            "observed_environment_identity": self.observed_environment_identity,
            "reference_manifest_identity": self.reference_manifest_identity,
            "historical_qualification": self.historical_qualification,
            "manifest_identity": self.manifest_identity,
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> NumericalProbeBaseline:
        _exact_keys(
            record,
            {
                "schema_version",
                "semantic_version",
                "definitions",
                "artifacts",
                "observed_lane_identity",
                "observed_lane_role",
                "observed_attestation_identity",
                "observed_environment_identity",
                "reference_manifest_identity",
                "historical_qualification",
                "manifest_identity",
                "identity",
            },
            "Numerical probe baseline",
        )
        if (
            record.get("schema_version") != _SCHEMA_VERSION
            or record.get("semantic_version") != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported numerical probe baseline semantics")
        raw_definitions = record.get("definitions")
        raw_artifacts = record.get("artifacts")
        if (
            not isinstance(raw_definitions, list)
            or not isinstance(raw_artifacts, list)
            or any(not isinstance(value, Mapping) for value in raw_definitions)
            or any(not isinstance(value, Mapping) for value in raw_artifacts)
        ):
            raise TypeError("Malformed numerical probe baseline")
        definitions = tuple(
            NumericalProbeDefinition.from_record(
                cast("Mapping[str, object]", definition)
            )
            for definition in raw_definitions
        )
        if len(definitions) != len(raw_artifacts):
            raise ValueError("Probe baseline definition and artifact counts differ")
        artifacts = tuple(
            NumericalProbeArtifact.from_record(
                cast("Mapping[str, object]", artifact),
                definition,
            )
            for definition, artifact in zip(
                definitions,
                raw_artifacts,
                strict=True,
            )
        )
        historical = record.get("historical_qualification")
        if historical not in {"CAPTURE_ONLY", "REFERENCE_MATCHED"}:
            raise ValueError("Unknown numerical probe historical qualification")
        lane_role = record.get("observed_lane_role")
        if lane_role not in {None, "CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"}:
            raise ValueError("Unknown observed numerical lane role")
        baseline = cls(
            definitions,
            artifacts,
            cast("str | None", record.get("observed_lane_identity")),
            cast("LaneRole | None", lane_role),
            cast("str | None", record.get("observed_attestation_identity")),
            cast("str | None", record.get("observed_environment_identity")),
            cast("str | None", record.get("reference_manifest_identity")),
            cast("Literal['CAPTURE_ONLY', 'REFERENCE_MATCHED']", historical),
        )
        if (
            record.get("manifest_identity") != baseline.manifest_identity
            or record.get("identity") != baseline.identity
        ):
            raise ValueError("Numerical probe baseline identity does not match payload")
        return baseline

    def require_live_qualification(
        self,
        authority: LiveLaneAuthority,
        *,
        expected_manifest_identity: str,
        required_lane_role: LaneRole,
    ) -> None:
        """Validate this evidence using #588's exact process-owned capability."""

        if not isinstance(authority, LiveLaneAuthority):
            raise TypeError("Probe qualification requires live lane authority")
        evidence = authority.evidence
        if len(expected_manifest_identity) != 64:
            raise ValueError("Expected probe manifest must be a SHA-256 digest")
        if required_lane_role != "CANONICAL_NUMERICAL":
            raise NumericalProbeQualificationError("EXPECTED_REFERENCE_UNAVAILABLE")
        if expected_manifest_identity != CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY:
            raise NumericalProbeQualificationError("EXPECTED_REFERENCE_MISMATCH")
        if self.manifest_identity != expected_manifest_identity:
            raise NumericalProbeQualificationError("MANIFEST_MISMATCH")
        if authority.lane_role != required_lane_role:
            raise NumericalProbeQualificationError("LANE_ROLE_MISMATCH")
        required_lane = next(
            lane for lane in canonical_lanes() if lane.role == required_lane_role
        )
        if authority.lane_identity != required_lane.identity:
            raise NumericalProbeQualificationError("LANE_MISMATCH")
        if self.historical_qualification == "REFERENCE_MATCHED" and (
            self.observed_lane_identity != evidence.lane_identity
            or self.observed_lane_role != authority.lane_role
            or self.observed_attestation_identity != evidence.identity
            or self.observed_environment_identity != evidence.environment_identity
            or self.reference_manifest_identity != expected_manifest_identity
        ):
            raise NumericalProbeQualificationError("LANE_MISMATCH")


class NumericalProbeExecutionError(RuntimeError):
    """Typed fail-closed terminal before a complete probe artifact exists."""

    def __init__(
        self,
        definition_identity: str,
        terminal: str,
        accounting: ObjectiveRequestAccounting,
        trajectory_fingerprint: str,
        failure_kind: RequestFailureKind | None = None,
    ) -> None:
        self.definition_identity = definition_identity
        self.terminal = terminal
        self.accounting = accounting
        self.trajectory_fingerprint = trajectory_fingerprint
        self.failure_kind = failure_kind
        super().__init__(f"Numerical probe terminated: {terminal}")


class _ProbeBudgetExhausted(RuntimeError):
    def __init__(
        self,
        accounting: ObjectiveRequestAccounting,
        trajectory_fingerprint: str,
    ) -> None:
        self.accounting = accounting
        self.trajectory_fingerprint = trajectory_fingerprint
        super().__init__("Probe objective-request budget exhausted")


class _ProbeRequestFailed(RuntimeError):
    def __init__(
        self,
        accounting: ObjectiveRequestAccounting,
        trajectory_fingerprint: str,
        failure_kind: RequestFailureKind,
    ) -> None:
        self.accounting = accounting
        self.trajectory_fingerprint = trajectory_fingerprint
        self.failure_kind = failure_kind
        super().__init__("Authoritative probe objective request failed")


class _ProbeBackendFailed(RuntimeError):
    pass


class _ProbeBackendExecutionFailed(RuntimeError):
    def __init__(
        self,
        accounting: ObjectiveRequestAccounting,
        trajectory_fingerprint: str,
    ) -> None:
        self.accounting = accounting
        self.trajectory_fingerprint = trajectory_fingerprint
        super().__init__("Probe backend execution failed")


class _ObjectiveRecorder:
    def __init__(
        self,
        residual: Callable[[tuple[float, ...]], tuple[float, ...]],
        budget: int,
    ) -> None:
        self._residual = residual
        self._budget = budget
        self._cache: dict[tuple[str, ...], tuple[float, ...]] = {}
        self._trace: list[object] = []
        self.workflow_requests = 0
        self.materialization_requests = 0
        self.fresh_evaluations = 0
        self.cache_served_requests = 0
        self.refused = 0
        self.failure_kinds: list[RequestFailureKind] = []

    def _fail(self, source: str, failure_kind: RequestFailureKind) -> Never:
        self.failure_kinds.append(failure_kind)
        self._trace.append((source, "failed", failure_kind))
        raise _ProbeRequestFailed(
            self.accounting,
            self.trajectory_fingerprint,
            failure_kind,
        )

    def evaluate(self, vector: Array | Sequence[float], source: str) -> Array:
        if source == "solver":
            self.workflow_requests += 1
        elif source == "materialization":
            self.materialization_requests += 1
        else:
            raise ValueError("Unknown probe objective request source")
        try:
            candidate = tuple(float(value) for value in vector)
            key = _vector_tokens(candidate)
        except (TypeError, ValueError, OverflowError):
            self._fail(source, "INVALID_INPUT")
        if self.requests_completed >= self._budget:
            self.refused += 1
            self._trace.append((source, key, "budget_refused"))
            raise _ProbeBudgetExhausted(
                self.accounting,
                self.trajectory_fingerprint,
            )
        cached = self._cache.get(key)
        disposition = "cache_hit"
        if cached is None:
            try:
                raw_result = self._residual(candidate)
            except Exception as error:  # noqa: BLE001 - evaluator boundary
                try:
                    self._fail(source, "EVALUATION_FAILURE")
                except _ProbeRequestFailed as failure:
                    raise failure from error
            try:
                cached = tuple(float(value) for value in raw_result)
                _vector_tokens(cached)
            except (TypeError, ValueError, OverflowError):
                self._fail(source, "NON_FINITE_RESULT")
            self._cache[key] = cached
            self.fresh_evaluations += 1
            disposition = "cache_miss"
        else:
            self.cache_served_requests += 1
        self._trace.append((source, key, disposition, _vector_tokens(cached)))
        return np.asarray(cached, dtype=np.float64)

    @property
    def requests_completed(self) -> int:
        return self.fresh_evaluations + self.cache_served_requests

    @property
    def accounting(self) -> ObjectiveRequestAccounting:
        return ObjectiveRequestAccounting(
            self.workflow_requests,
            self.materialization_requests,
            self.fresh_evaluations,
            self.cache_served_requests,
            self.refused,
            tuple(self.failure_kinds),
        )

    @property
    def trajectory_fingerprint(self) -> str:
        return _identity("probe-objective-trajectory", tuple(self._trace))


class _BackendObjectiveAdapter:
    """Count SciPy calls without granting them ChemEx request authority."""

    def __init__(
        self,
        residual: Callable[[tuple[float, ...]], tuple[float, ...]],
        budget: int,
    ) -> None:
        self._residual = residual
        self._budget = budget
        self.calls = 0

    def evaluate(self, vector: Array | Sequence[float]) -> Array:
        if self.calls >= self._budget:
            raise _ProbeBackendFailed
        self.calls += 1
        try:
            candidate = tuple(float(value) for value in vector)
            values = tuple(float(value) for value in self._residual(candidate))
            _vector_tokens(candidate)
            _vector_tokens(values)
        except Exception as error:
            raise _ProbeBackendFailed from error
        return np.asarray(values, dtype=np.float64)


def _budget_value(definition: NumericalProbeDefinition, field_name: str) -> int:
    budget = definition.budget.to_record_value()
    if not isinstance(budget, Mapping):
        raise TypeError("Probe budget must be a record")
    value = budget.get(field_name)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"Probe budget {field_name!r} must be positive")
    return value


@dataclass(frozen=True, slots=True)
class _TrfSettings:
    ftol: float
    xtol: float
    gtol: float
    x_scale: float

    @classmethod
    def from_definition(
        cls,
        definition: NumericalProbeDefinition,
        prefix: str = "",
    ) -> _TrfSettings:
        if prefix == "" and (
            _policy_string(definition, "backend") != "scipy.optimize.least_squares"
            or _policy_string(definition, "method") != "trf"
        ):
            raise ValueError("TRF probe policy selects an unsupported backend")
        return cls(
            _policy_float(definition, f"{prefix}ftol"),
            _policy_float(definition, f"{prefix}xtol"),
            _policy_float(definition, f"{prefix}gtol"),
            _policy_float(definition, f"{prefix}x_scale"),
        )


def _least_squares_evidence(
    *,
    start: tuple[float, ...],
    lower: tuple[float, ...],
    upper: tuple[float, ...],
    residual: Callable[[tuple[float, ...]], tuple[float, ...]],
    authoritative_budget: int,
    backend_max_nfev: int,
    settings: _TrfSettings,
) -> tuple[SolverProbeEvidence, bool]:
    recorder = _ObjectiveRecorder(residual, authoritative_budget)
    recorder.evaluate(start, "solver")
    backend = _BackendObjectiveAdapter(residual, backend_max_nfev * (len(start) + 1))

    try:
        result = least_squares(
            backend.evaluate,
            np.asarray(start, dtype=np.float64),
            bounds=(
                np.asarray(lower, dtype=np.float64),
                np.asarray(upper, dtype=np.float64),
            ),
            method="trf",
            ftol=settings.ftol,
            xtol=settings.xtol,
            gtol=settings.gtol,
            x_scale=settings.x_scale,
            max_nfev=backend_max_nfev,
        )
    except Exception as error:
        raise _ProbeBackendExecutionFailed(
            recorder.accounting,
            recorder.trajectory_fingerprint,
        ) from error
    accepted = tuple(float(value) for value in result.x)
    residuals = tuple(
        float(value) for value in recorder.evaluate(accepted, "materialization")
    )
    chi_square = math.fsum(value * value for value in residuals)
    candidate_fingerprint = _identity(
        "probe-candidate",
        (_vector_tokens(accepted), _vector_tokens(residuals), _float_token(chi_square)),
    )
    return (
        SolverProbeEvidence(
            start,
            accepted,
            lower,
            upper,
            residuals,
            chi_square,
            tuple(int(value) for value in result.active_mask),
            recorder.accounting,
            BackendDiagnosticCounters(
                int(result.nfev),
                None if result.njev is None else int(result.njev),
                None,
                0,
                backend.calls,
            ),
            recorder.trajectory_fingerprint,
            candidate_fingerprint,
        ),
        bool(result.success),
    )


def _sum_accounting(
    accounting: Sequence[ObjectiveRequestAccounting],
) -> ObjectiveRequestAccounting:
    return ObjectiveRequestAccounting(
        sum(item.workflow_requests for item in accounting),
        sum(item.materialization_requests for item in accounting),
        sum(item.fresh_evaluations for item in accounting),
        sum(item.cache_served_requests for item in accounting),
        sum(item.requests_refused for item in accounting),
        tuple(kind for item in accounting for kind in item.failure_kinds),
    )


def _sum_diagnostics(
    diagnostics: Sequence[BackendDiagnosticCounters],
) -> BackendDiagnosticCounters:
    njev = (
        None
        if any(item.njev is None for item in diagnostics)
        else sum(cast("int", item.njev) for item in diagnostics)
    )
    iterations = (
        None
        if any(item.iterations is None for item in diagnostics)
        else sum(cast("int", item.iterations) for item in diagnostics)
    )
    return BackendDiagnosticCounters(
        sum(item.nfev for item in diagnostics),
        njev,
        iterations,
        sum(item.callback_calls for item in diagnostics),
        sum(item.diagnostic_evaluations for item in diagnostics),
    )


def _solver_matches_truth(
    definition: NumericalProbeDefinition,
    evidence: SolverProbeEvidence,
) -> bool:
    expected = _policy_vector(definition, "truth")
    tolerance = _policy_float(definition, "truth_absolute_tolerance")
    if len(expected) != len(evidence.accepted) or any(
        not math.isclose(actual, truth, rel_tol=0.0, abs_tol=tolerance)
        for actual, truth in zip(evidence.accepted, expected, strict=True)
    ):
        return False
    if definition.category == "BOUNDS":
        return evidence.active_mask == _policy_int_vector(
            definition,
            "expected_active_mask",
        )
    return evidence.chi_square <= _policy_float(
        definition,
        "maximum_chi_square",
    )


def _run_solver_probe(
    definition: NumericalProbeDefinition,
) -> NumericalProbeArtifact:
    start = _policy_vector(definition, "start")
    lower, upper = _policy_bounds(definition)
    if definition.category == "TRF_ROUTINE":

        def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
            x_value, y_value = vector
            return (x_value - 1.5, y_value + 0.5)

    elif definition.category == "TRF_DIFFICULT":

        def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
            x_value, y_value = vector
            return (10.0 * (y_value - x_value * x_value), 1.0 - x_value)

    elif definition.category == "BOUNDS":

        def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
            return (vector[0] - 2.0,)

    else:
        raise ValueError(f"Probe {definition.probe_id!r} is not a solver probe")
    evidence, succeeded = _least_squares_evidence(
        start=start,
        lower=lower,
        upper=upper,
        residual=residual,
        authoritative_budget=_budget_value(
            definition, "authoritative_objective_requests"
        ),
        backend_max_nfev=_budget_value(definition, "backend_max_nfev"),
        settings=_TrfSettings.from_definition(definition),
    )
    qualified = succeeded and _solver_matches_truth(definition, evidence)
    terminal = "CONVERGED" if qualified else "TRUTH_MISMATCH"
    risks = ("ACTIVE_BOUND",) if qualified and definition.category == "BOUNDS" else ()
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        terminal,
        risks,
        definition.eligible_claims if qualified else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def _run_grid_probe(definition: NumericalProbeDefinition) -> NumericalProbeArtifact:
    lower, upper = _policy_bounds(definition)

    def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
        x_value, y_value, z_value = vector
        return (x_value - 1.0, y_value + 1.0, z_value - 0.5)

    axes = _policy_matrix(definition, "physical_axes")
    starts = tuple(product(*axes))
    expected_seed_count = _budget_value(definition, "seed_count")
    if len(starts) != expected_seed_count:
        raise ValueError("GRID policy axes do not match its declared seed budget")
    seeds: list[GridSeedEvidence] = []
    for ordinal, start in enumerate(starts):
        solver, succeeded = _least_squares_evidence(
            start=start,
            lower=lower,
            upper=upper,
            residual=residual,
            authoritative_budget=_budget_value(
                definition, "authoritative_objective_requests_per_seed"
            ),
            backend_max_nfev=_budget_value(definition, "backend_max_nfev_per_seed"),
            settings=_TrfSettings.from_definition(definition),
        )
        seeds.append(
            GridSeedEvidence(
                ordinal,
                start,
                "CONVERGED" if succeeded else "NON_CONVERGED",
                solver,
            )
        )
    seed_evidence = tuple(seeds)
    order = tuple(
        seed.ordinal
        for seed in sorted(
            seed_evidence,
            key=lambda seed: (
                seed.solver.chi_square,
                seed.solver.accepted,
                seed.ordinal,
            ),
        )
    )
    evidence = GridProbeEvidence(
        seed_evidence,
        order,
        order[0],
        _sum_accounting(
            tuple(seed.solver.objective_accounting for seed in seed_evidence)
        ),
        _sum_diagnostics(
            tuple(seed.solver.backend_diagnostics for seed in seed_evidence)
        ),
        _identity(
            "grid-probe-trajectory",
            tuple(seed.solver.trajectory_fingerprint for seed in seed_evidence),
        ),
    )
    selected = seed_evidence[order[0]].solver
    qualified = (
        all(seed.terminal == "CONVERGED" for seed in seed_evidence)
        and _solver_matches_truth(definition, selected)
        and _policy_string(definition, "candidate_tie_break")
        == "chi-square-then-seed-ordinal"
    )
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        "SELECTED" if qualified else "TRUTH_MISMATCH",
        (),
        definition.eligible_claims if qualified else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def _run_grid_ordering_probe(
    definition: NumericalProbeDefinition,
) -> NumericalProbeArtifact:
    raw_candidates = _policy_record(definition).get("candidates")
    if not isinstance(raw_candidates, list) or any(
        not isinstance(value, Mapping) for value in raw_candidates
    ):
        raise TypeError("GRID ordering policy candidates must be records")
    candidates: list[GridOrderingCandidate] = []
    for index, raw_candidate in enumerate(raw_candidates):
        candidate_record = cast("Mapping[str, object]", raw_candidate)
        _exact_keys(
            candidate_record,
            {"ordinal", "final_vector", "chi_square"},
            f"GRID ordering policy candidate {index}",
        )
        candidates.append(
            GridOrderingCandidate(
                _record_nonnegative_int(
                    candidate_record.get("ordinal"), "GRID ordering ordinal"
                ),
                _semantic_vector(
                    candidate_record.get("final_vector"),
                    f"GRID ordering final_vector[{index}]",
                ),
                _semantic_float(
                    candidate_record.get("chi_square"),
                    f"GRID ordering chi_square[{index}]",
                ),
            )
        )
    candidate_evidence = tuple(candidates)
    if len(candidate_evidence) != _budget_value(definition, "candidate_count"):
        raise ValueError("GRID ordering policy differs from its declared budget")
    if _policy_string(definition, "candidate_tie_break") != (
        "chi-square-then-final-vector-then-seed-ordinal"
    ):
        raise ValueError("GRID ordering policy selects an unsupported tie break")
    observed_order = tuple(
        candidate.ordinal
        for candidate in sorted(
            candidate_evidence,
            key=lambda candidate: (
                candidate.chi_square,
                candidate.final_vector,
                candidate.ordinal,
            ),
        )
    )
    expected_order = _policy_int_vector(definition, "expected_order")
    evidence = GridOrderingProbeEvidence(
        candidate_evidence,
        observed_order,
        expected_order,
        observed_order[0],
        ObjectiveRequestAccounting(0, 0, 0, 0, 0, ()),
        BackendDiagnosticCounters(0, None, None, 0, 0),
        _identity(
            "grid-ordering-trajectory",
            tuple(candidate.identity for candidate in candidate_evidence),
        ),
    )
    qualified = observed_order == expected_order
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        "ORDERED" if qualified else "ORDERING_MISMATCH",
        (),
        definition.eligible_claims if qualified else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def _run_de_probe(definition: NumericalProbeDefinition) -> NumericalProbeArtifact:
    root_seed = definition.seed.to_record_value()
    if isinstance(root_seed, bool) or not isinstance(root_seed, int):
        raise TypeError("DE probe root seed must be an integer")
    oscillation_amplitude = _policy_float(definition, "oscillation_amplitude")
    oscillation_frequency = _policy_float(definition, "oscillation_frequency")

    def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
        displacement = vector[0] - 0.25
        return (
            displacement,
            oscillation_amplitude * math.sin(oscillation_frequency * displacement),
        )

    search_recorder = _ObjectiveRecorder(
        residual,
        _budget_value(definition, "de_authoritative_objective_requests"),
    )
    backend = _BackendObjectiveAdapter(
        residual,
        _budget_value(definition, "de_backend_objective_calls"),
    )

    def objective(vector: Array) -> float:
        values = backend.evaluate(vector)
        return float(np.dot(values, values))

    callback_calls = 0

    def count_callback(intermediate_result: OptimizeResult) -> None:
        nonlocal callback_calls
        _ = intermediate_result
        callback_calls += 1

    lower, upper = _policy_bounds(definition)
    mutation = _policy_vector(definition, "mutation")
    if len(mutation) != 2:
        raise ValueError("DE mutation policy requires one declared range")

    strategy = _policy_string(definition, "de_strategy")
    if strategy != "best1bin":
        raise ValueError("DE probe policy selects an unsupported strategy")
    try:
        result = differential_evolution(
            objective,
            Bounds(np.asarray(lower), np.asarray(upper)),
            strategy="best1bin",
            maxiter=_policy_int(definition, "maximum_generations"),
            popsize=_policy_int(definition, "population_multiplier"),
            mutation=mutation,
            recombination=_policy_float(definition, "recombination"),
            tol=_policy_float(definition, "relative_tolerance"),
            atol=_policy_float(definition, "absolute_tolerance"),
            polish=_policy_bool(definition, "backend_polish"),
            rng=np.random.default_rng(root_seed),
            updating=cast(
                "Literal['immediate', 'deferred']",
                _policy_string(definition, "updating"),
            ),
            workers=_policy_int(definition, "workers"),
            callback=count_callback,
        )
    except Exception as error:
        raise _ProbeBackendExecutionFailed(
            search_recorder.accounting,
            search_recorder.trajectory_fingerprint,
        ) from error
    population = tuple(
        DePopulationCandidate(
            index,
            tuple(float(value) for value in vector),
            float(
                np.dot(
                    values := search_recorder.evaluate(vector, "materialization"),
                    values,
                )
            ),
            float(result.population_energies[index]),
        )
        for index, vector in enumerate(result.population)
    )
    order = tuple(
        candidate.population_index
        for candidate in sorted(
            population,
            key=lambda candidate: (candidate.objective, candidate.population_index),
        )
    )
    selected_vector = population[order[0]].vector
    polish, polish_succeeded = _least_squares_evidence(
        start=selected_vector,
        lower=lower,
        upper=upper,
        residual=residual,
        authoritative_budget=_budget_value(
            definition, "trf_authoritative_objective_requests"
        ),
        backend_max_nfev=_budget_value(definition, "trf_backend_max_nfev"),
        settings=_TrfSettings.from_definition(definition, "polish_"),
    )
    evidence = DeProbeEvidence(
        root_seed,
        population,
        order,
        order[0],
        search_recorder.accounting,
        BackendDiagnosticCounters(
            int(result.nfev),
            None,
            int(result.nit),
            callback_calls,
            backend.calls,
        ),
        polish,
        _identity(
            "de-probe-trajectory",
            (
                search_recorder.trajectory_fingerprint,
                tuple(candidate.fingerprint for candidate in population),
                polish.trajectory_fingerprint,
            ),
        ),
    )
    backend_selected_vector = tuple(float(value) for value in result.x)
    qualified = (
        polish_succeeded
        and _solver_matches_truth(definition, polish)
        and all(
            math.isclose(actual, selected, rel_tol=0.0, abs_tol=0.0)
            for actual, selected in zip(
                backend_selected_vector,
                selected_vector,
                strict=True,
            )
        )
        and _policy_string(definition, "workflow_polish") == "trf"
    )
    terminal = "POLISHED" if qualified else "TRUTH_MISMATCH"
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        terminal,
        (),
        definition.eligible_claims if qualified else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def _policy_record(definition: NumericalProbeDefinition) -> Mapping[str, object]:
    policy = definition.policy.to_record_value()
    if not isinstance(policy, Mapping):
        raise TypeError("Probe policy must be a record")
    return cast("Mapping[str, object]", policy)


def _semantic_float(value: object, field_name: str) -> float:
    raw_value = value
    if not isinstance(raw_value, Mapping) or set(raw_value) != {"binary64"}:
        raise TypeError(f"Probe policy {field_name!r} must be binary64")
    return _float_from_record(raw_value.get("binary64"), field_name)


def _policy_float(definition: NumericalProbeDefinition, field_name: str) -> float:
    return _semantic_float(_policy_record(definition).get(field_name), field_name)


def _policy_int(definition: NumericalProbeDefinition, field_name: str) -> int:
    value = _policy_record(definition).get(field_name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"Probe policy {field_name!r} must be an integer")
    return value


def _policy_string(definition: NumericalProbeDefinition, field_name: str) -> str:
    return _nonempty_string(
        _policy_record(definition).get(field_name),
        f"policy {field_name}",
    )


def _policy_bool(definition: NumericalProbeDefinition, field_name: str) -> bool:
    value = _policy_record(definition).get(field_name)
    if not isinstance(value, bool):
        raise TypeError(f"Probe policy {field_name!r} must be boolean")
    return value


def _semantic_vector(value: object, field_name: str) -> tuple[float, ...]:
    if not isinstance(value, list) or not value:
        raise TypeError(f"Probe policy {field_name!r} must be a non-empty list")
    return tuple(
        _semantic_float(item, f"{field_name}[{index}]")
        for index, item in enumerate(value)
    )


def _policy_vector(
    definition: NumericalProbeDefinition,
    field_name: str,
) -> tuple[float, ...]:
    return _semantic_vector(_policy_record(definition).get(field_name), field_name)


def _policy_matrix(
    definition: NumericalProbeDefinition,
    field_name: str,
) -> tuple[tuple[float, ...], ...]:
    value = _policy_record(definition).get(field_name)
    if not isinstance(value, list) or not value:
        raise TypeError(f"Probe policy {field_name!r} must be a non-empty matrix")
    return tuple(
        _semantic_vector(row, f"{field_name}[{index}]")
        for index, row in enumerate(value)
    )


def _policy_int_vector(
    definition: NumericalProbeDefinition,
    field_name: str,
) -> tuple[int, ...]:
    value = _policy_record(definition).get(field_name)
    if not isinstance(value, list) or any(
        isinstance(item, bool) or not isinstance(item, int) for item in value
    ):
        raise TypeError(f"Probe policy {field_name!r} must be an integer list")
    return tuple(cast("list[int]", value))


def _policy_bounds(
    definition: NumericalProbeDefinition,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    rows = _policy_matrix(definition, "bounds")
    if any(len(row) != 2 or row[0] >= row[1] for row in rows):
        raise ValueError("Probe bounds must be ordered lower/upper pairs")
    return tuple(row[0] for row in rows), tuple(row[1] for row in rows)


def _run_finite_difference_probe(
    definition: NumericalProbeDefinition,
) -> NumericalProbeArtifact:
    if _policy_string(definition, "stencil") != "centered-five-point":
        raise ValueError("Finite-difference probe requires its declared stencil")
    if (
        _policy_string(definition, "weight_derivation")
        != "lagrange-first-derivative-at-zero"
    ):
        raise ValueError("Finite-difference probe requires Lagrange weights")
    point = _policy_float(definition, "point")
    step = _policy_float(definition, "step")
    nominal_steps = (-2.0 * step, -step, 0.0, step, 2.0 * step)
    sample_points = tuple(float(point + displacement) for displacement in nominal_steps)
    actual_steps = tuple(sample - point for sample in sample_points)

    def residual(vector: tuple[float, ...]) -> tuple[float, ...]:
        return (math.exp(vector[0]),)

    recorder = _ObjectiveRecorder(
        residual,
        _budget_value(definition, "authoritative_objective_requests"),
    )
    sampled = tuple(
        float(recorder.evaluate((sample,), "solver")[0]) for sample in sample_points
    )
    weights = tuple(
        math.fsum(
            math.prod(
                -actual_steps[other]
                for other in range(len(actual_steps))
                if other not in (index, omitted)
            )
            / math.prod(
                actual_steps[index] - actual_steps[other]
                for other in range(len(actual_steps))
                if other != index
            )
            for omitted in range(len(actual_steps))
            if omitted != index
        )
        for index in range(len(actual_steps))
    )
    estimate = math.fsum(
        weight * value for weight, value in zip(weights, sampled, strict=True)
    )
    truth = math.exp(point)
    absolute_error = abs(estimate - truth)
    absolute_tolerance = _policy_float(definition, "absolute_tolerance")
    evidence = FiniteDifferenceProbeEvidence(
        point,
        actual_steps,
        weights,
        sampled,
        truth,
        estimate,
        absolute_error,
        absolute_tolerance,
        recorder.accounting,
        BackendDiagnosticCounters(0, None, None, 0, 0),
        recorder.trajectory_fingerprint,
    )
    reliable = absolute_error <= absolute_tolerance
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        "RELIABLE" if reliable else "UNRELIABLE",
        (),
        definition.eligible_claims if reliable else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def _run_spectral_risk_probe(
    definition: NumericalProbeDefinition,
) -> NumericalProbeArtifact:
    if definition.category == "RANK":
        jacobian = ((1.0, 1.0), (2.0, 2.0))
        relative_tolerance = _policy_float(definition, "rank_relative_tolerance")
        condition_limit = _policy_float(definition, "conditioning_limit")
    elif definition.category == "CONDITIONING":
        jacobian = ((1.0, 0.0), (0.0, 1.0e-8))
        relative_tolerance = 0.0
        condition_limit = _policy_float(definition, "conditioning_limit")
    else:
        raise ValueError(f"Probe {definition.probe_id!r} is not spectral")
    singular_values = tuple(
        float(value)
        for value in svd(
            np.asarray(jacobian, dtype=np.float64),
            compute_uv=False,
            lapack_driver=cast(
                "Literal['gesdd', 'gesvd']",
                _policy_string(definition, "svd_driver"),
            ),
        )
    )
    rank_threshold = singular_values[0] * relative_tolerance
    rank = sum(value > rank_threshold for value in singular_values)
    dimension = len(jacobian[0])
    condition = None if rank < dimension else singular_values[0] / singular_values[-1]
    risks = (
        ("RANK_DEFICIENT",)
        if rank < dimension
        else (
            ("ILL_CONDITIONED",) if cast("float", condition) > condition_limit else ()
        )
    )
    evidence = SpectralRiskProbeEvidence(
        jacobian,
        singular_values,
        rank,
        dimension,
        condition,
        condition_limit,
        risks,
        ObjectiveRequestAccounting(0, 0, 0, 0, 0, ()),
        BackendDiagnosticCounters(0, None, None, 0, 1),
        _identity(
            "spectral-probe-trajectory",
            (
                tuple(_vector_tokens(row) for row in jacobian),
                _vector_tokens(singular_values),
                _float_token(rank_threshold),
                rank,
                risks,
            ),
        ),
    )
    expected_rank = _policy_int(definition, "expected_rank")
    qualified = rank == expected_rank and bool(risks)
    if definition.category == "CONDITIONING":
        expected_condition = _policy_float(definition, "expected_condition")
        qualified = (
            qualified
            and condition is not None
            and math.isclose(
                condition,
                expected_condition,
                rel_tol=_policy_float(definition, "condition_relative_tolerance"),
                abs_tol=0.0,
            )
        )
    artifact = NumericalProbeArtifact(
        definition.identity,
        definition.probe_id,
        "RISK_DETECTED" if qualified else "TRUTH_MISMATCH",
        risks if qualified else (),
        definition.eligible_claims if qualified else (),
        evidence,
    )
    artifact.validate_definition(definition)
    return artifact


def run_numerical_probe(
    definition: NumericalProbeDefinition,
) -> NumericalProbeArtifact:
    """Execute one predeclared reduced case and return typed baseline evidence."""
    canonical = next(
        (
            candidate
            for candidate in numerical_probe_definitions()
            if candidate.probe_id == definition.probe_id
        ),
        None,
    )
    if canonical != definition:
        raise ValueError("Only an exact predeclared numerical probe may run")
    runners: Mapping[
        ProbeCategory,
        Callable[[NumericalProbeDefinition], NumericalProbeArtifact],
    ] = {
        "TRF_ROUTINE": _run_solver_probe,
        "TRF_DIFFICULT": _run_solver_probe,
        "GRID": _run_grid_probe,
        "GRID_ORDERING": _run_grid_ordering_probe,
        "DE_SEARCH": _run_de_probe,
        "FINITE_DIFFERENCE": _run_finite_difference_probe,
        "BOUNDS": _run_solver_probe,
        "RANK": _run_spectral_risk_probe,
        "CONDITIONING": _run_spectral_risk_probe,
    }
    try:
        return runners[definition.category](definition)
    except _ProbeBudgetExhausted as error:
        raise NumericalProbeExecutionError(
            definition.identity,
            "BUDGET_EXHAUSTED",
            error.accounting,
            error.trajectory_fingerprint,
        ) from error
    except _ProbeRequestFailed as error:
        raise NumericalProbeExecutionError(
            definition.identity,
            "REQUEST_FAILED",
            error.accounting,
            error.trajectory_fingerprint,
            error.failure_kind,
        ) from error
    except _ProbeBackendExecutionFailed as error:
        raise NumericalProbeExecutionError(
            definition.identity,
            "BACKEND_FAILURE",
            error.accounting,
            error.trajectory_fingerprint,
        ) from error


def run_numerical_probe_baseline(
    authority: LiveLaneAuthority | None = None,
    *,
    expected_manifest_identity: str | None = None,
    required_lane_role: LaneRole = "CANONICAL_NUMERICAL",
) -> NumericalProbeBaseline:
    """Capture the closed catalogue and optionally record a live-qualified replay."""
    definitions = numerical_probe_definitions()
    artifacts = tuple(run_numerical_probe(definition) for definition in definitions)
    capture = NumericalProbeBaseline(definitions, artifacts)
    if expected_manifest_identity is None:
        return capture
    if authority is None:
        raise TypeError(
            "Qualified numerical replay requires live authority and expected manifest"
        )
    capture.require_live_qualification(
        authority,
        expected_manifest_identity=expected_manifest_identity,
        required_lane_role=required_lane_role,
    )
    evidence = authority.evidence
    return NumericalProbeBaseline(
        definitions,
        artifacts,
        evidence.lane_identity,
        authority.lane_role,
        evidence.identity,
        evidence.environment_identity,
        expected_manifest_identity,
        "REFERENCE_MATCHED",
    )


def _definition(
    probe_id: str,
    category: ProbeCategory,
    *,
    parent_case: str,
    reduction: str,
    truth_kind: TruthAuthorityKind,
    truth_reference: str,
    eligible_claims: Sequence[str],
    policy: Mapping[str, object],
    budget: Mapping[str, object],
    seed: object,
    terminal: str,
    required_risks: Sequence[str] = (),
) -> NumericalProbeDefinition:
    return NumericalProbeDefinition(
        probe_id,
        category,
        ProbeDerivation(591, _SOURCE_REVISION, parent_case, reduction),
        ProbeTruthAuthority(truth_kind, truth_reference),
        tuple(sorted(eligible_claims)),
        CanonicalBaselineValue.from_value(policy),
        CanonicalBaselineValue.from_value(budget),
        CanonicalBaselineValue.from_value(seed),
        (
            ProbeFailureExpectation(
                terminal,
                tuple(sorted(required_risks)),
            ),
        ),
    )


def numerical_probe_definitions() -> tuple[NumericalProbeDefinition, ...]:
    """Return the closed v1 reduced-case catalogue in canonical run order."""
    trf_policy = {
        "backend": "scipy.optimize.least_squares",
        "method": "trf",
        "ftol": 1.0e-10,
        "xtol": 1.0e-10,
        "gtol": 1.0e-10,
        "x_scale": 1.0,
    }
    return (
        _definition(
            "trf-routine-quadratic-v1",
            "TRF_ROUTINE",
            parent_case="analytic bounded linear residual",
            reduction="Two independent coordinates retain routine TRF convergence.",
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference="r(x,y)=(x-1.5,y+0.5); unique minimizer=(1.5,-0.5).",
            eligible_claims=(
                "deterministic-trajectory-fingerprint",
                "routine-trf-convergence",
                "solver-request-accounting",
            ),
            policy={
                **trf_policy,
                "start": [0.0, 0.0],
                "bounds": [[-4.0, 4.0], [-4.0, 4.0]],
                "truth": [1.5, -0.5],
                "truth_absolute_tolerance": 1.0e-8,
                "maximum_chi_square": 1.0e-16,
            },
            budget={
                "authoritative_objective_requests": 2,
                "backend_max_nfev": 64,
            },
            seed="not-applicable",
            terminal="CONVERGED",
        ),
        _definition(
            "trf-difficult-rosenbrock-v1",
            "TRF_DIFFICULT",
            parent_case="Rosenbrock least-squares valley",
            reduction="Canonical poor start retains curved-valley TRF behavior.",
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference="r(x,y)=(10*(y-x^2),1-x); unique zero=(1,1).",
            eligible_claims=(
                "deterministic-trajectory-fingerprint",
                "difficult-trf-convergence",
                "solver-request-accounting",
            ),
            policy={
                **trf_policy,
                "start": [-1.2, 1.0],
                "bounds": [[-3.0, 3.0], [-3.0, 3.0]],
                "truth": [1.0, 1.0],
                "truth_absolute_tolerance": 1.0e-8,
                "maximum_chi_square": 1.0e-16,
            },
            budget={
                "authoritative_objective_requests": 2,
                "backend_max_nfev": 256,
            },
            seed="not-applicable",
            terminal="CONVERGED",
        ),
        _definition(
            "grid-27-seed-coverage-v1",
            "GRID",
            parent_case="immutable 27-seed physical-coordinate grid",
            reduction="Three three-point axes retain complete Cartesian seed coverage.",
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference="r(x,y,z)=(x-1,y+1,z-0.5); unique zero=(1,-1,0.5).",
            eligible_claims=(
                "grid-seed-coverage",
                "solver-request-accounting",
            ),
            policy={
                **trf_policy,
                "physical_axes": [
                    [-2.0, 0.0, 2.0],
                    [-2.0, 0.0, 2.0],
                    [-2.0, 0.0, 2.0],
                ],
                "bounds": [[-3.0, 3.0], [-3.0, 3.0], [-3.0, 3.0]],
                "candidate_tie_break": "chi-square-then-seed-ordinal",
                "truth": [1.0, -1.0, 0.5],
                "truth_absolute_tolerance": 1.0e-8,
                "maximum_chi_square": 1.0e-16,
            },
            budget={
                "authoritative_objective_requests_per_seed": 2,
                "backend_max_nfev_per_seed": 48,
                "seed_count": 27,
            },
            seed="not-applicable",
            terminal="SELECTED",
        ),
        _definition(
            "grid-candidate-ordering-v1",
            "GRID_ORDERING",
            parent_case="four-candidate deterministic GRID selection",
            reduction=(
                "Four analytically scored candidates isolate chi-square, final-vector, "
                "and seed-ordinal ordering."
            ),
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference=(
                "Candidates (ordinal,chi2,vector)=(0,1,(1,0)),(1,0,(0,0)),"
                "(2,1,(0,1)),(3,1,(-1,0)); lexicographic ordering by "
                "(chi2,vector,ordinal) is exactly (1,3,2,0)."
            ),
            eligible_claims=("canonical-grid-candidate-ordering",),
            policy={
                "candidate_tie_break": (
                    "chi-square-then-final-vector-then-seed-ordinal"
                ),
                "candidates": [
                    {"ordinal": 0, "chi_square": 1.0, "final_vector": [1.0, 0.0]},
                    {"ordinal": 1, "chi_square": 0.0, "final_vector": [0.0, 0.0]},
                    {"ordinal": 2, "chi_square": 1.0, "final_vector": [0.0, 1.0]},
                    {"ordinal": 3, "chi_square": 1.0, "final_vector": [-1.0, 0.0]},
                ],
                "expected_order": [1, 3, 2, 0],
            },
            budget={"candidate_count": 4},
            seed="not-applicable",
            terminal="ORDERED",
        ),
        _definition(
            "de-bounded-search-v1",
            "DE_SEARCH",
            parent_case="bounded multimodal basin search",
            reduction=(
                "One selected coordinate retains oscillatory local basins, seeded DE, "
                "and separate TRF polish."
            ),
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference=(
                "f(x)=(x-0.25)^2+0.25 sin^2(5(x-0.25)); its non-negative "
                "first term makes x=0.25 the unique bounded global minimizer."
            ),
            eligible_claims=(
                "de-seeded-replay",
                "de-to-trf-candidate-ordering",
                "solver-request-accounting",
            ),
            policy={
                "de_strategy": "best1bin",
                "bounds": [[-2.0, 2.0]],
                "oscillation_amplitude": 0.5,
                "oscillation_frequency": 5.0,
                "population_multiplier": 5,
                "maximum_generations": 6,
                "mutation": [0.5, 1.0],
                "recombination": 0.7,
                "relative_tolerance": 0.0,
                "absolute_tolerance": 0.0,
                "updating": "immediate",
                "workers": 1,
                "backend_polish": False,
                "workflow_polish": "trf",
                "polish_ftol": 1.0e-10,
                "polish_xtol": 1.0e-10,
                "polish_gtol": 1.0e-10,
                "polish_x_scale": 1.0,
                "truth": [0.25],
                "truth_absolute_tolerance": 1.0e-8,
                "maximum_chi_square": 1.0e-16,
            },
            budget={
                "de_authoritative_objective_requests": 5,
                "de_backend_objective_calls": 128,
                "trf_authoritative_objective_requests": 2,
                "trf_backend_max_nfev": 48,
            },
            seed=591,
            terminal="POLISHED",
        ),
        _definition(
            "finite-difference-reliability-v1",
            "FINITE_DIFFERENCE",
            parent_case="analytic exponential derivative",
            reduction="Five-point centered stencil retains actual-step reliability.",
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference="d exp(x)/dx at x=0.5 equals exp(0.5).",
            eligible_claims=(
                "finite-difference-reliability",
                "stencil-request-accounting",
                "stencil-trajectory-fingerprint",
            ),
            policy={
                "stencil": "centered-five-point",
                "weight_derivation": "lagrange-first-derivative-at-zero",
                "point": 0.5,
                "step": 1.0e-4,
                "absolute_tolerance": 1.0e-10,
            },
            budget={"authoritative_objective_requests": 5},
            seed="not-applicable",
            terminal="RELIABLE",
        ),
        _definition(
            "trf-active-bound-v1",
            "BOUNDS",
            parent_case="bounded scalar residual",
            reduction="An infeasible unconstrained optimum retains active-bound truth.",
            truth_kind="ANALYTIC_DERIVATION",
            truth_reference="r(x)=x-2 on [0,1]; constrained minimizer x=1.",
            eligible_claims=(
                "active-bound-detection",
                "bound-safe-trf-convergence",
                "solver-request-accounting",
            ),
            policy={
                **trf_policy,
                "start": [0.5],
                "bounds": [[0.0, 1.0]],
                "truth": [1.0],
                "truth_absolute_tolerance": 1.0e-8,
                "expected_active_mask": [1],
            },
            budget={
                "authoritative_objective_requests": 2,
                "backend_max_nfev": 64,
            },
            seed="not-applicable",
            terminal="CONVERGED",
            required_risks=("ACTIVE_BOUND",),
        ),
        _definition(
            "rank-deficient-linearization-v1",
            "RANK",
            parent_case="exact rank-one Jacobian",
            reduction="A duplicated column retains fail-closed rank diagnosis.",
            truth_kind="EXACT_LINEAR_ALGEBRA",
            truth_reference="J=((1,1),(2,2)); exact rank=1 < 2.",
            eligible_claims=("rank-risk-detection",),
            policy={
                "rank_relative_tolerance": 1.0e-12,
                "conditioning_limit": 1.0e12,
                "svd_driver": "gesdd",
                "expected_rank": 1,
            },
            budget={"svd_decompositions": 1},
            seed="not-applicable",
            terminal="RISK_DETECTED",
            required_risks=("RANK_DEFICIENT",),
        ),
        _definition(
            "ill-conditioned-linearization-v1",
            "CONDITIONING",
            parent_case="exact diagonal Jacobian spectrum",
            reduction="Two separated singular values retain conditioning diagnosis.",
            truth_kind="EXACT_LINEAR_ALGEBRA",
            truth_reference="J=diag(1,1e-8); exact kappa_2(J)=1e8.",
            eligible_claims=("conditioning-risk-detection",),
            policy={
                "conditioning_limit": 1.0e6,
                "svd_driver": "gesdd",
                "expected_rank": 2,
                "expected_condition": 1.0e8,
                "condition_relative_tolerance": 1.0e-12,
            },
            budget={"svd_decompositions": 1},
            seed="not-applicable",
            terminal="RISK_DETECTED",
            required_risks=("ILL_CONDITIONED",),
        ),
    )
