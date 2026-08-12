"""Deterministic native resampling evidence qualification (#599).

This module is an isolated qualification seam.  It plans identity-derived
Monte Carlo, bootstrap, and nucleus-bootstrap replicates; delegates each draw
to an evidence-only projected optimization callback; freezes every terminal
replicate outcome; and derives summaries from one common complete-scope sample
set.  Production fitting and reporting do not consume these artifacts yet.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from collections.abc import Callable, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real

import numpy as np

from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    accepted_occurrence_is_authoritative,
)
from chemex.runtime import ExecutionSettings
from chemex.runtime.execution import native_thread_environment
from chemex.typing import Array

_SCHEMA_VERSION = 1
_SEED_POLICY_VERSION = "sha256-structural-u64-v1"
_RNG_ALGORITHM = "numpy-pcg64-v1"
_GENERATION_POLICY_VERSION = "chemex-mc-bs-bsn-v1"
_SUMMARY_POLICY_VERSION = "complete-scope-percentile-v1"
_MAX_U64 = (1 << 64) - 1


class ResamplingConstructionError(ValueError):
    """Raised when a native resampling policy or artifact is malformed."""


class ResamplingScheme(StrEnum):
    """Closed scientific replicate-generation schemes retained by ChemEx."""

    MONTE_CARLO = "mc"
    BOOTSTRAP = "bs"
    NUCLEUS_BOOTSTRAP = "bsn"


class OptimizationStrategy(StrEnum):
    """Closed projected optimization strategies available to a replicate."""

    DIRECT_TRF = "direct_trf"
    GRID_DIRECT_TRF = "grid_direct_trf"
    DE_DIRECT_TRF = "de_direct_trf"


class ReplicateDisposition(StrEnum):
    """Terminal disposition of one actually executed replicate."""

    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


class ResamplingLifecycle(StrEnum):
    """Lifecycle of an evidence artifact containing genuine outcomes."""

    COMPLETED = "completed"
    PARTIAL = "partial"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


class OperationTerminal(StrEnum):
    """Terminal disposition of the overall resampling operation."""

    COMPLETED = "completed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


class SummaryTerminal(StrEnum):
    """Terminal disposition of a scope-atomic resampling summary request."""

    COMPLETED = "completed"
    INSUFFICIENT_COVERAGE = "insufficient_coverage"
    SOURCE_INVALID = "source_invalid"


class CorrelationAvailability(StrEnum):
    """Availability of one empirical correlation entry."""

    AVAILABLE = "available"
    ZERO_VARIANCE = "zero_variance"


class ClaimState(StrEnum):
    """Closed assessment state for a resampling validity claim."""

    SATISFIED = "satisfied"
    VIOLATED = "violated"
    INDETERMINATE = "indeterminate"
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, slots=True)
class ClaimAssessment:
    """One versioned lifecycle, integrity, or scientific validity claim."""

    name: str
    state: ClaimState
    policy_identity: str
    details: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name or not self.policy_identity:
            raise ResamplingConstructionError(
                "Claim assessment requires a name and policy identity"
            )


def _claim_state(claims: Sequence[ClaimAssessment], name: str) -> ClaimState:
    for claim in claims:
        if claim.name == name:
            return claim.state
    return ClaimState.INDETERMINATE


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ResamplingConstructionError(f"{name} must be a real binary64 scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ResamplingConstructionError(f"{name} must be finite")
    return 0.0 if result == 0.0 else result


def _float_token(value: float) -> str:
    return float(value).hex()


def _items_tokens(items: Sequence[tuple[str, float]]) -> tuple[tuple[str, str], ...]:
    return tuple((param_id, _float_token(value)) for param_id, value in items)


def _non_empty_identity(value: str, *, name: str) -> str:
    if not value:
        raise ResamplingConstructionError(f"{name} must not be empty")
    return value


def _canonical_scope(scope: Sequence[str]) -> tuple[str, ...]:
    result = tuple(scope)
    if (
        not result
        or any(not item for item in result)
        or len(set(result)) != len(result)
    ):
        raise ResamplingConstructionError(
            "Resampling output scope must be non-empty, ordered, and unique"
        )
    return result


def _positive_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ResamplingConstructionError(f"{name} must be a positive integer")
    return value


def _unsigned_seed(value: object, *, name: str) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < 0
        or value > _MAX_U64
    ):
        raise ResamplingConstructionError(f"{name} must be an unsigned 64-bit seed")
    return value


def _canonical_field(value: bytes) -> bytes:
    return len(value).to_bytes(8, "big") + value


def _derive_seed(
    *,
    root_seed: int,
    stage_identity: str,
    work_unit_kind: str,
    structural_identity: str,
) -> int:
    fields = (
        b"chemex-native-seed",
        _SEED_POLICY_VERSION.encode("ascii"),
        root_seed.to_bytes(8, "big"),
        stage_identity.encode("ascii"),
        work_unit_kind.encode("ascii"),
        structural_identity.encode("ascii"),
    )
    digest = hashlib.sha256(
        b"".join(_canonical_field(item) for item in fields)
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


@dataclass(frozen=True, slots=True)
class ReplicateRequest:
    """One canonical structural replicate identity and its derived seed."""

    plan_identity: str
    ordinal: int
    structural_identity: str
    optimization_projection_identity: str
    component_identities: tuple[str, ...]
    seed: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.ordinal < 0:
            raise ResamplingConstructionError("Replicate ordinal must be non-negative")
        _non_empty_identity(
            self.optimization_projection_identity,
            name="optimization projection identity",
        )
        if (
            not self.component_identities
            or any(not item for item in self.component_identities)
            or len(set(self.component_identities)) != len(self.component_identities)
        ):
            raise ResamplingConstructionError(
                "Replicate request requires unique projected component identities"
            )
        _unsigned_seed(self.seed, name="replicate seed")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-replicate-request",
                (
                    self.plan_identity,
                    self.ordinal,
                    self.structural_identity,
                    self.optimization_projection_identity,
                    self.component_identities,
                    self.seed,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ResamplingPlan:
    """Immutable accepted-anchor resampling population and validity contract."""

    accepted_result_identity: str
    accepted_occurrence_identity: str
    scheme: ResamplingScheme
    replicate_count: int
    replicate_structural_identities: tuple[str, ...]
    replicate_component_identities: tuple[tuple[str, ...], ...]
    root_seed: int
    source_dataset_identity: str
    output_scope: tuple[str, ...]
    optimization_projection_identity: str
    minimum_successful_count: int
    seed_policy_version: str = _SEED_POLICY_VERSION
    rng_algorithm: str = _RNG_ALGORITHM
    generation_policy_version: str = _GENERATION_POLICY_VERSION
    identity: str = field(init=False)
    replicates: tuple[ReplicateRequest, ...] = field(init=False)

    def __post_init__(self) -> None:
        count = _positive_integer(self.replicate_count, name="replicate count")
        minimum = _positive_integer(
            self.minimum_successful_count,
            name="minimum successful replicate count",
        )
        if minimum > count:
            raise ResamplingConstructionError(
                "Minimum successful replicate count cannot exceed the population"
            )
        supplied_structural_identities = tuple(self.replicate_structural_identities)
        supplied_component_identities = tuple(
            tuple(items) for items in self.replicate_component_identities
        )
        if (
            len(supplied_structural_identities) != count
            or any(not item for item in supplied_structural_identities)
            or len(set(supplied_structural_identities)) != count
        ):
            raise ResamplingConstructionError(
                "Resampling requires one unique structural identity per replicate"
            )
        if len(supplied_component_identities) != count or any(
            not items
            or any(not item for item in items)
            or len(set(items)) != len(items)
            for items in supplied_component_identities
        ):
            raise ResamplingConstructionError(
                "Resampling requires complete unique component identities per replicate"
            )
        canonical_work_units = tuple(
            sorted(
                zip(
                    supplied_structural_identities,
                    supplied_component_identities,
                    strict=True,
                )
            )
        )
        structural_identities = tuple(item[0] for item in canonical_work_units)
        component_identities = tuple(item[1] for item in canonical_work_units)
        root_seed = _unsigned_seed(self.root_seed, name="root seed")
        scope = _canonical_scope(self.output_scope)
        for name, value in (
            ("accepted result identity", self.accepted_result_identity),
            ("accepted occurrence identity", self.accepted_occurrence_identity),
            ("source dataset identity", self.source_dataset_identity),
            ("optimization projection identity", self.optimization_projection_identity),
            ("seed policy version", self.seed_policy_version),
            ("RNG algorithm", self.rng_algorithm),
            ("generation policy version", self.generation_policy_version),
        ):
            _non_empty_identity(value, name=name)
        identity = _identity(
            "native-resampling-plan",
            (
                self.accepted_result_identity,
                self.accepted_occurrence_identity,
                self.scheme.value,
                count,
                structural_identities,
                component_identities,
                root_seed,
                self.source_dataset_identity,
                scope,
                self.optimization_projection_identity,
                minimum,
                self.seed_policy_version,
                self.rng_algorithm,
                self.generation_policy_version,
            ),
        )
        requests: list[ReplicateRequest] = []
        for ordinal, (structural_identity, components) in enumerate(
            zip(structural_identities, component_identities, strict=True)
        ):
            seed = _derive_seed(
                root_seed=root_seed,
                stage_identity=identity,
                work_unit_kind="resampling-replicate",
                structural_identity=structural_identity,
            )
            requests.append(
                ReplicateRequest(
                    identity,
                    ordinal,
                    structural_identity,
                    self.optimization_projection_identity,
                    components,
                    seed,
                )
            )
        seeds = tuple(item.seed for item in requests)
        if len(set(seeds)) != len(seeds):
            raise ResamplingConstructionError(
                "Distinct resampling work units produced a derived-seed collision"
            )
        object.__setattr__(self, "replicate_count", count)
        object.__setattr__(
            self,
            "replicate_structural_identities",
            structural_identities,
        )
        object.__setattr__(
            self,
            "replicate_component_identities",
            component_identities,
        )
        object.__setattr__(self, "minimum_successful_count", minimum)
        object.__setattr__(self, "root_seed", root_seed)
        object.__setattr__(self, "output_scope", scope)
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "replicates", tuple(requests))

    @classmethod
    def for_accepted(
        cls,
        accepted: AcceptedFitResult,
        *,
        scheme: ResamplingScheme,
        replicate_count: int,
        replicate_structural_identities: Sequence[str],
        replicate_component_identities: Sequence[Sequence[str]],
        root_seed: int,
        source_dataset_identity: str,
        output_scope: Sequence[str],
        optimization_projection_identity: str,
        minimum_successful_count: int,
    ) -> ResamplingPlan:
        """Bind a deterministic replicate population to one accepted occurrence."""
        if not accepted_occurrence_is_authoritative(accepted):
            raise ResamplingConstructionError(
                "Resampling planning requires an exact authoritative accepted occurrence"
            )
        return cls(
            accepted.identity,
            accepted.occurrence_identity,
            scheme,
            replicate_count,
            tuple(replicate_structural_identities),
            tuple(tuple(items) for items in replicate_component_identities),
            root_seed,
            source_dataset_identity,
            tuple(output_scope),
            optimization_projection_identity,
            minimum_successful_count,
        )


@dataclass(frozen=True, slots=True)
class ResamplingPopulation:
    """Canonical source values needed by all three retained schemes."""

    source_dataset_identity: str
    observed: tuple[float, ...]
    calculated: tuple[float, ...]
    standard_errors: tuple[float, ...]
    mask: tuple[bool, ...]
    references: tuple[bool, ...]
    nucleus_groups: tuple[str, ...]
    profile_blocks: tuple[str, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _non_empty_identity(
            self.source_dataset_identity,
            name="source dataset identity",
        )
        size = len(self.observed)
        if size < 1 or any(
            len(values) != size
            for values in (
                self.calculated,
                self.standard_errors,
                self.mask,
                self.references,
                self.nucleus_groups,
                self.profile_blocks,
            )
        ):
            raise ResamplingConstructionError(
                "Resampling population fields must have one common non-zero size"
            )
        observed = tuple(
            _finite(value, name=f"observed[{index}]")
            for index, value in enumerate(self.observed)
        )
        calculated = tuple(
            _finite(value, name=f"calculated[{index}]")
            for index, value in enumerate(self.calculated)
        )
        errors = tuple(
            _finite(value, name=f"standard_errors[{index}]")
            for index, value in enumerate(self.standard_errors)
        )
        if any(value < 0.0 for value in errors):
            raise ResamplingConstructionError(
                "Resampling standard errors must be non-negative"
            )
        if any(type(value) is not bool for value in (*self.mask, *self.references)):
            raise ResamplingConstructionError(
                "Mask and reference flags must be Boolean"
            )
        if any(not group for group in self.nucleus_groups):
            raise ResamplingConstructionError(
                "Nucleus group identities must not be empty"
            )
        if any(not block for block in self.profile_blocks):
            raise ResamplingConstructionError(
                "Profile block identities must not be empty"
            )
        object.__setattr__(self, "observed", observed)
        object.__setattr__(self, "calculated", calculated)
        object.__setattr__(self, "standard_errors", errors)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-population",
                (
                    self.source_dataset_identity,
                    tuple(_float_token(value) for value in observed),
                    tuple(_float_token(value) for value in calculated),
                    tuple(_float_token(value) for value in errors),
                    self.mask,
                    self.references,
                    self.nucleus_groups,
                    self.profile_blocks,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ResamplingDraw:
    """One immutable seeded transformation with exact source selection."""

    population_identity: str
    scheme: ResamplingScheme
    seed: int
    observations: tuple[float, ...]
    standard_errors: tuple[float, ...]
    mask: tuple[bool, ...]
    references: tuple[bool, ...]
    nucleus_groups: tuple[str, ...]
    profile_blocks: tuple[str, ...]
    source_indices: tuple[int, ...]
    sampled_nucleus_groups: tuple[str, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _unsigned_seed(self.seed, name="draw seed")
        size = len(self.observations)
        if any(
            len(values) != size
            for values in (
                self.standard_errors,
                self.mask,
                self.references,
                self.nucleus_groups,
                self.profile_blocks,
                self.source_indices,
            )
        ):
            raise ResamplingConstructionError(
                "Resampling draw fields must preserve one common row extent"
            )
        observations = tuple(
            _finite(value, name=f"draw observation[{index}]")
            for index, value in enumerate(self.observations)
        )
        errors = tuple(
            _finite(value, name=f"draw error[{index}]")
            for index, value in enumerate(self.standard_errors)
        )
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "standard_errors", errors)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-draw",
                (
                    self.population_identity,
                    self.scheme.value,
                    self.seed,
                    tuple(_float_token(value) for value in observations),
                    tuple(_float_token(value) for value in errors),
                    self.mask,
                    self.references,
                    self.nucleus_groups,
                    self.profile_blocks,
                    self.source_indices,
                    self.sampled_nucleus_groups,
                ),
            ),
        )


def generate_resampling_draw(
    population: ResamplingPopulation,
    scheme: ResamplingScheme,
    *,
    seed: int,
) -> ResamplingDraw:
    """Apply one retained scientific generation scheme with explicit PCG64 state."""
    normalized_seed = _unsigned_seed(seed, name="draw seed")
    rng = np.random.Generator(np.random.PCG64(normalized_seed))
    size = len(population.observed)
    if scheme is ResamplingScheme.MONTE_CARLO:
        observations = tuple(
            float(value)
            for value in rng.normal(population.calculated, population.standard_errors)
        )
        return ResamplingDraw(
            population.identity,
            scheme,
            normalized_seed,
            observations,
            population.standard_errors,
            population.mask,
            population.references,
            population.nucleus_groups,
            population.profile_blocks,
            tuple(range(size)),
        )
    if scheme is ResamplingScheme.BOOTSTRAP:
        source_indices = list(range(size))
        for block in dict.fromkeys(population.profile_blocks):
            masked = tuple(
                index
                for index, active in enumerate(population.mask)
                if active and population.profile_blocks[index] == block
            )
            reference_pool = tuple(
                index for index in masked if population.references[index]
            )
            non_reference_pool = tuple(
                index for index in masked if not population.references[index]
            )
            selected = sorted(
                (
                    *(
                        int(value)
                        for value in rng.choice(
                            reference_pool,
                            size=len(reference_pool),
                            replace=True,
                        )
                    ),
                    *(
                        int(value)
                        for value in rng.choice(
                            non_reference_pool,
                            size=len(non_reference_pool),
                            replace=True,
                        )
                    ),
                )
            )
            for target, source in zip(masked, selected, strict=True):
                source_indices[target] = source
        indexes = tuple(source_indices)
        return ResamplingDraw(
            population.identity,
            scheme,
            normalized_seed,
            tuple(population.observed[index] for index in indexes),
            tuple(population.standard_errors[index] for index in indexes),
            population.mask,
            tuple(population.references[index] for index in indexes),
            population.nucleus_groups,
            population.profile_blocks,
            indexes,
        )
    canonical_groups = tuple(sorted(set(population.nucleus_groups)))
    sampled_groups = tuple(
        str(value)
        for value in rng.choice(
            canonical_groups,
            size=len(canonical_groups),
            replace=True,
        )
    )
    indexes = tuple(
        index
        for group in sampled_groups
        for index, source_group in enumerate(population.nucleus_groups)
        if source_group == group
    )
    return ResamplingDraw(
        population.identity,
        scheme,
        normalized_seed,
        tuple(population.observed[index] for index in indexes),
        tuple(population.standard_errors[index] for index in indexes),
        tuple(population.mask[index] for index in indexes),
        tuple(population.references[index] for index in indexes),
        tuple(population.nucleus_groups[index] for index in indexes),
        tuple(population.profile_blocks[index] for index in indexes),
        indexes,
        sampled_groups,
    )


@dataclass(frozen=True, slots=True)
class ProjectedOptimizationSuccess:
    """Complete projected candidate evidence with no acceptance or commit authority."""

    transformed_data_identity: str
    optimization_projection_identity: str
    evaluation_plan_identity: str
    problem_identity: str
    invocation_identity: str
    execution_identity: str
    strategy: OptimizationStrategy
    component_identities: tuple[str, ...]
    component_outcome_identities: tuple[str, ...]
    resolved_items: tuple[tuple[str, float], ...]
    chi_square: float
    expected_output_scope: tuple[str, ...] | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for name, value in (
            ("transformed data identity", self.transformed_data_identity),
            ("optimization projection identity", self.optimization_projection_identity),
            ("evaluation plan identity", self.evaluation_plan_identity),
            ("problem identity", self.problem_identity),
            ("invocation identity", self.invocation_identity),
            ("execution identity", self.execution_identity),
        ):
            _non_empty_identity(value, name=name)
        if (
            not self.component_outcome_identities
            or any(not item for item in self.component_outcome_identities)
            or len(self.component_outcome_identities) != len(self.component_identities)
        ):
            raise ResamplingConstructionError(
                "Projected optimization requires complete component outcomes"
            )
        if (
            not self.component_identities
            or any(not item for item in self.component_identities)
            or len(set(self.component_identities)) != len(self.component_identities)
        ):
            raise ResamplingConstructionError(
                "Projected optimization requires unique component identities"
            )
        items = tuple(
            (param_id, _finite(value, name=f"resolved value {param_id!r}"))
            for param_id, value in self.resolved_items
        )
        ids = tuple(param_id for param_id, _value in items)
        if (
            not ids
            or any(not param_id for param_id in ids)
            or len(set(ids)) != len(ids)
        ):
            raise ResamplingConstructionError(
                "Projected optimization requires unique resolved output IDs"
            )
        if self.expected_output_scope is not None and ids != tuple(
            self.expected_output_scope
        ):
            raise ResamplingConstructionError(
                "Projected optimization must resolve the complete output scope"
            )
        chi_square = _finite(self.chi_square, name="projected chi-square")
        if chi_square < 0.0:
            raise ResamplingConstructionError(
                "Projected chi-square must be non-negative"
            )
        object.__setattr__(self, "resolved_items", items)
        object.__setattr__(self, "chi_square", chi_square)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-projected-optimization-success",
                (
                    self.transformed_data_identity,
                    self.optimization_projection_identity,
                    self.evaluation_plan_identity,
                    self.problem_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.strategy.value,
                    self.component_identities,
                    self.component_outcome_identities,
                    _items_tokens(items),
                    _float_token(chi_square),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ProjectedOptimizationFailure:
    """Typed unsuccessful projected path retaining all available lineage."""

    transformed_data_identity: str
    disposition: ReplicateDisposition
    category: str
    message: str
    optimization_projection_identity: str = field(kw_only=True)
    evaluation_plan_identity: str | None = None
    problem_identity: str | None = None
    execution_identity: str | None = None
    component_outcome_identities: tuple[str, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.disposition is ReplicateDisposition.SUCCEEDED:
            raise ResamplingConstructionError(
                "Projected failure cannot use the successful disposition"
            )
        if not self.transformed_data_identity or not self.category or not self.message:
            raise ResamplingConstructionError(
                "Projected failure requires data identity, category, and message"
            )
        _non_empty_identity(
            self.optimization_projection_identity,
            name="optimization projection identity",
        )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-projected-optimization-failure",
                (
                    self.transformed_data_identity,
                    self.disposition.value,
                    self.category,
                    self.message,
                    self.optimization_projection_identity,
                    self.evaluation_plan_identity,
                    self.problem_identity,
                    self.execution_identity,
                    self.component_outcome_identities,
                ),
            ),
        )


type ProjectedOptimizationResult = (
    ProjectedOptimizationSuccess | ProjectedOptimizationFailure
)


type ReplicateExecutor = Callable[
    [ReplicateRequest, ResamplingDraw],
    ProjectedOptimizationResult,
]


@dataclass(frozen=True, slots=True)
class ReplicateOutcome:
    """One genuine terminal outcome in canonical replicate order."""

    plan_identity: str
    request_identity: str
    ordinal: int
    seed: int
    draw_identity: str
    disposition: ReplicateDisposition
    success: ProjectedOptimizationSuccess | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    failure: ProjectedOptimizationFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        succeeded = self.disposition is ReplicateDisposition.SUCCEEDED
        if succeeded != (self.success is not None) or succeeded == (
            self.failure is not None
        ):
            raise ResamplingConstructionError(
                "Replicate outcome must contain exactly its typed terminal payload"
            )
        payload = self.success if self.success is not None else self.failure
        if payload is None:
            raise ResamplingConstructionError(
                "Replicate outcome lacks a terminal payload"
            )
        if payload.transformed_data_identity != self.draw_identity:
            raise ResamplingConstructionError(
                "Replicate outcome and projected path use different transformed data"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-replicate-outcome",
                (
                    self.plan_identity,
                    self.request_identity,
                    self.ordinal,
                    self.seed,
                    self.draw_identity,
                    self.disposition.value,
                    payload.identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class ResamplingEvidence:
    """Exact frozen outcome set with separate coverage and lifecycle claims."""

    plan: ResamplingPlan = field(repr=False, compare=False)
    population_identity: str
    outcomes: tuple[ReplicateOutcome, ...]
    operation_terminal: OperationTerminal = OperationTerminal.COMPLETED
    lifecycle: ResamplingLifecycle = field(init=False)
    completed_count: int = field(init=False)
    successful_count: int = field(init=False)
    failed_count: int = field(init=False)
    coverage_satisfied: bool = field(init=False)
    claims: tuple[ClaimAssessment, ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.outcomes:
            raise ResamplingConstructionError(
                "Resampling evidence requires at least one genuine replicate outcome"
            )
        ordinals = tuple(outcome.ordinal for outcome in self.outcomes)
        if ordinals != tuple(sorted(ordinals)) or len(set(ordinals)) != len(ordinals):
            raise ResamplingConstructionError(
                "Resampling outcomes must use unique canonical ordinal order"
            )
        if any(
            ordinal >= self.plan.replicate_count
            or outcome.plan_identity != self.plan.identity
            or outcome.request_identity != self.plan.replicates[ordinal].identity
            or outcome.seed != self.plan.replicates[ordinal].seed
            for ordinal, outcome in zip(ordinals, self.outcomes, strict=True)
        ):
            raise ResamplingConstructionError(
                "Resampling outcomes differ from their planned replicate identities"
            )
        for outcome in self.outcomes:
            if (
                outcome.success is not None
                and tuple(
                    param_id for param_id, _value in outcome.success.resolved_items
                )
                != self.plan.output_scope
            ):
                raise ResamplingConstructionError(
                    "Successful replicate must resolve the complete output scope"
                )
        completed = len(self.outcomes)
        successful = sum(
            outcome.disposition is ReplicateDisposition.SUCCEEDED
            for outcome in self.outcomes
        )
        failed = sum(
            outcome.disposition is ReplicateDisposition.FAILED
            for outcome in self.outcomes
        )
        if self.operation_terminal is OperationTerminal.INTERRUPTED:
            lifecycle = ResamplingLifecycle.INTERRUPTED
        elif self.operation_terminal is OperationTerminal.CANCELLED:
            lifecycle = ResamplingLifecycle.CANCELLED
        elif completed == self.plan.replicate_count:
            lifecycle = ResamplingLifecycle.COMPLETED
        else:
            lifecycle = ResamplingLifecycle.PARTIAL
        object.__setattr__(self, "completed_count", completed)
        object.__setattr__(self, "successful_count", successful)
        object.__setattr__(self, "failed_count", failed)
        object.__setattr__(
            self,
            "coverage_satisfied",
            successful >= self.plan.minimum_successful_count,
        )
        claims = (
            ClaimAssessment(
                "STRUCTURAL_INTEGRITY",
                ClaimState.SATISFIED,
                self.plan.identity,
                (f"terminal_outcomes={completed}",),
            ),
            ClaimAssessment(
                "INTENDED_POPULATION_TERMINAL",
                (
                    ClaimState.SATISFIED
                    if completed == self.plan.replicate_count
                    else ClaimState.VIOLATED
                ),
                self.plan.identity,
                (
                    f"intended={self.plan.replicate_count}",
                    f"terminal={completed}",
                ),
            ),
            ClaimAssessment(
                "MINIMUM_SUCCESSFUL_COVERAGE",
                (
                    ClaimState.SATISFIED
                    if successful >= self.plan.minimum_successful_count
                    else ClaimState.VIOLATED
                ),
                self.plan.identity,
                (
                    f"required={self.plan.minimum_successful_count}",
                    f"successful={successful}",
                ),
            ),
            ClaimAssessment(
                "COMPLETE_SCOPE_SUCCESS_ROWS",
                (ClaimState.SATISFIED if successful > 0 else ClaimState.NOT_APPLICABLE),
                self.plan.identity,
                (f"scope={','.join(self.plan.output_scope)}",),
            ),
        )
        object.__setattr__(self, "claims", claims)
        object.__setattr__(self, "lifecycle", lifecycle)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-evidence",
                (
                    self.plan.identity,
                    self.population_identity,
                    tuple(outcome.identity for outcome in self.outcomes),
                    lifecycle.value,
                    completed,
                    successful,
                    failed,
                    successful >= self.plan.minimum_successful_count,
                    tuple(
                        (
                            claim.name,
                            claim.state.value,
                            claim.policy_identity,
                            claim.details,
                        )
                        for claim in claims
                    ),
                ),
            ),
        )

    def claim(self, name: str) -> ClaimState:
        """Return one exact named claim; unknown claims fail closed."""
        return _claim_state(self.claims, name)


@dataclass(frozen=True, slots=True)
class ResamplingOperation:
    """Operation record plus optional truthfully frozen evidence artifact."""

    plan_identity: str
    terminal: OperationTerminal
    evidence: ResamplingEvidence | None
    unstarted_ordinals: tuple[int, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.evidence is None and not self.unstarted_ordinals:
            raise ResamplingConstructionError(
                "An empty operation must retain its planned unstarted population"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-operation",
                (
                    self.plan_identity,
                    self.terminal.value,
                    None if self.evidence is None else self.evidence.identity,
                    self.unstarted_ordinals,
                ),
            ),
        )


def _execute_replicate(
    plan: ResamplingPlan,
    population: ResamplingPopulation,
    executor: ReplicateExecutor,
    request: ReplicateRequest,
) -> ReplicateOutcome:
    try:
        draw = generate_resampling_draw(population, plan.scheme, seed=request.seed)
        projected = executor(request, draw)
    except KeyboardInterrupt:
        failure = ProjectedOptimizationFailure(
            draw.identity,
            ReplicateDisposition.INTERRUPTED,
            "asynchronous_interruption",
            "Replicate projected optimization was interrupted",
            optimization_projection_identity=request.optimization_projection_identity,
        )
        return ReplicateOutcome(
            plan.identity,
            request.identity,
            request.ordinal,
            request.seed,
            draw.identity,
            failure.disposition,
            failure=failure,
        )
    except Exception as error:  # noqa: BLE001 - freeze typed replicate failure
        draw_identity = (
            draw.identity
            if "draw" in locals()
            else _identity(
                "native-resampling-draw-failure",
                (population.identity, plan.scheme.value, request.seed),
            )
        )
        failure = ProjectedOptimizationFailure(
            draw_identity,
            ReplicateDisposition.FAILED,
            "projected_execution_failure",
            f"{type(error).__name__}: {error}",
            optimization_projection_identity=request.optimization_projection_identity,
        )
        return ReplicateOutcome(
            plan.identity,
            request.identity,
            request.ordinal,
            request.seed,
            draw_identity,
            failure.disposition,
            failure=failure,
        )
    if isinstance(projected, ProjectedOptimizationSuccess):
        if tuple(param_id for param_id, _value in projected.resolved_items) != (
            plan.output_scope
        ):
            failure = ProjectedOptimizationFailure(
                draw.identity,
                ReplicateDisposition.FAILED,
                "incomplete_output_scope",
                "Projected optimization did not resolve the complete output scope",
                projected.evaluation_plan_identity,
                projected.problem_identity,
                projected.execution_identity,
                projected.component_outcome_identities,
                optimization_projection_identity=(
                    request.optimization_projection_identity
                ),
            )
            return ReplicateOutcome(
                plan.identity,
                request.identity,
                request.ordinal,
                request.seed,
                draw.identity,
                failure.disposition,
                failure=failure,
            )
        payload: ProjectedOptimizationResult = projected
    elif isinstance(projected, ProjectedOptimizationFailure):
        payload = projected
    else:
        failure = ProjectedOptimizationFailure(
            draw.identity,
            ReplicateDisposition.FAILED,
            "invalid_projected_outcome",
            "Projected optimization returned an unsupported outcome type",
            optimization_projection_identity=request.optimization_projection_identity,
        )
        payload = failure
    if payload.transformed_data_identity != draw.identity:
        mismatch = ProjectedOptimizationFailure(
            draw.identity,
            ReplicateDisposition.FAILED,
            "transformed_data_lineage_mismatch",
            "Projected optimization used a different transformed-data identity",
            optimization_projection_identity=request.optimization_projection_identity,
        )
        payload = mismatch
    if (
        payload.optimization_projection_identity
        != request.optimization_projection_identity
    ):
        payload = ProjectedOptimizationFailure(
            draw.identity,
            ReplicateDisposition.FAILED,
            "optimization_projection_lineage_mismatch",
            "Projected optimization used a different declared projection",
            optimization_projection_identity=request.optimization_projection_identity,
        )
    if isinstance(payload, ProjectedOptimizationSuccess) and (
        payload.component_identities != request.component_identities
    ):
        payload = ProjectedOptimizationFailure(
            draw.identity,
            ReplicateDisposition.FAILED,
            "incomplete_component_projection",
            "Projected optimization did not return every planned component",
            optimization_projection_identity=request.optimization_projection_identity,
            evaluation_plan_identity=payload.evaluation_plan_identity,
            problem_identity=payload.problem_identity,
            execution_identity=payload.execution_identity,
            component_outcome_identities=payload.component_outcome_identities,
        )
    if isinstance(payload, ProjectedOptimizationSuccess):
        return ReplicateOutcome(
            plan.identity,
            request.identity,
            request.ordinal,
            request.seed,
            draw.identity,
            ReplicateDisposition.SUCCEEDED,
            success=payload,
        )
    return ReplicateOutcome(
        plan.identity,
        request.identity,
        request.ordinal,
        request.seed,
        draw.identity,
        payload.disposition,
        failure=payload,
    )


def _execute_serial_replicates(
    plan: ResamplingPlan,
    population: ResamplingPopulation,
    executor: ReplicateExecutor,
    cancellation_probe: Callable[[], bool] | None,
) -> tuple[list[ReplicateOutcome], OperationTerminal]:
    outcomes: list[ReplicateOutcome] = []
    terminal = OperationTerminal.COMPLETED
    for request in plan.replicates:
        if cancellation_probe is not None and cancellation_probe():
            return outcomes, OperationTerminal.CANCELLED
        outcome = _execute_replicate(plan, population, executor, request)
        outcomes.append(outcome)
        if outcome.disposition is ReplicateDisposition.CANCELLED:
            return outcomes, OperationTerminal.CANCELLED
        if outcome.disposition is ReplicateDisposition.INTERRUPTED:
            return outcomes, OperationTerminal.INTERRUPTED
    return outcomes, terminal


def _replicate_terminal(
    outcomes: Sequence[ReplicateOutcome],
) -> OperationTerminal:
    if any(
        outcome.disposition is ReplicateDisposition.INTERRUPTED for outcome in outcomes
    ):
        return OperationTerminal.INTERRUPTED
    if any(
        outcome.disposition is ReplicateDisposition.CANCELLED for outcome in outcomes
    ):
        return OperationTerminal.CANCELLED
    return OperationTerminal.COMPLETED


def _execute_parallel_replicates(
    plan: ResamplingPlan,
    population: ResamplingPopulation,
    executor: ReplicateExecutor,
    execution: ExecutionSettings,
    cancellation_probe: Callable[[], bool] | None,
) -> tuple[list[ReplicateOutcome], OperationTerminal]:
    worker_count = min(execution.workers, plan.replicate_count)
    outcomes: list[ReplicateOutcome] = []
    terminal = OperationTerminal.COMPLETED
    with native_thread_environment(execution.native_thread_env(parallel=True)):
        pool = ThreadPoolExecutor(max_workers=worker_count)
        pending: set[Future[ReplicateOutcome]] = {
            pool.submit(_execute_replicate, plan, population, executor, request)
            for request in plan.replicates
        }
        try:
            while pending:
                done, pending = wait(
                    pending,
                    timeout=0.01 if cancellation_probe is not None else None,
                    return_when=FIRST_COMPLETED,
                )
                outcomes.extend(future.result() for future in done)
                terminal = _replicate_terminal(outcomes)
                if cancellation_probe is not None and cancellation_probe():
                    terminal = OperationTerminal.CANCELLED
                if terminal is not OperationTerminal.COMPLETED:
                    break
        except KeyboardInterrupt:
            terminal = OperationTerminal.INTERRUPTED
        finally:
            cancelled_before_start: set[Future[ReplicateOutcome]] = set()
            for future in pending:
                if future.cancel():
                    cancelled_before_start.add(future)
            pool.shutdown(wait=True, cancel_futures=True)
            outcomes.extend(
                future.result()
                for future in pending - cancelled_before_start
                if not future.cancelled()
            )
    outcomes.sort(key=lambda outcome: outcome.ordinal)
    return outcomes, terminal


def execute_resampling_evidence(
    accepted: AcceptedFitResult,
    plan: ResamplingPlan,
    population: ResamplingPopulation,
    executor: ReplicateExecutor,
    *,
    execution: ExecutionSettings | None = None,
    cancellation_probe: Callable[[], bool] | None = None,
) -> ResamplingOperation:
    """Execute canonical evidence-only replicates without analysis-state authority."""
    if (
        not accepted_occurrence_is_authoritative(accepted)
        or accepted.identity != plan.accepted_result_identity
        or accepted.occurrence_identity != plan.accepted_occurrence_identity
    ):
        raise ResamplingConstructionError(
            "Resampling execution requires the plan's exact authoritative anchor"
        )
    if population.source_dataset_identity != plan.source_dataset_identity:
        raise ResamplingConstructionError(
            "Resampling population belongs to a different source dataset"
        )
    settings = ExecutionSettings() if execution is None else execution
    if cancellation_probe is not None and cancellation_probe():
        return ResamplingOperation(
            plan.identity,
            OperationTerminal.CANCELLED,
            None,
            tuple(range(plan.replicate_count)),
        )
    if settings.workers == 1:
        outcomes, terminal = _execute_serial_replicates(
            plan,
            population,
            executor,
            cancellation_probe,
        )
    else:
        outcomes, terminal = _execute_parallel_replicates(
            plan,
            population,
            executor,
            settings,
            cancellation_probe,
        )
    evidence = (
        ResamplingEvidence(plan, population.identity, tuple(outcomes), terminal)
        if outcomes
        else None
    )
    completed_ordinals = {outcome.ordinal for outcome in outcomes}
    unstarted = tuple(
        request.ordinal
        for request in plan.replicates
        if request.ordinal not in completed_ordinals
    )
    return ResamplingOperation(plan.identity, terminal, evidence, unstarted)


@dataclass(frozen=True, slots=True)
class SummaryFailure:
    """Typed reason why no resampling summary artifact was produced."""

    category: str
    message: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.category or not self.message:
            raise ResamplingConstructionError(
                "Summary failure requires category and message"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-summary-failure", (self.category, self.message)
            ),
        )


@dataclass(frozen=True, slots=True)
class SummarySample:
    """One complete-scope row included in a joint empirical distribution."""

    ordinal: int
    outcome_identity: str
    items: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class ResamplingExclusion:
    """One explicit unsuccessful source outcome excluded by policy."""

    ordinal: int
    outcome_identity: str
    disposition: ReplicateDisposition
    category: str


@dataclass(frozen=True, slots=True)
class ParameterDistribution:
    """One parameter's summary over the exact common included sample set."""

    parameter_id: str
    mean: float
    standard_deviation: float
    median: float
    percentile_95_lower: float
    percentile_95_upper: float

    def __post_init__(self) -> None:
        if not self.parameter_id:
            raise ResamplingConstructionError(
                "Parameter distribution requires a stable parameter ID"
            )
        for name in (
            "mean",
            "standard_deviation",
            "median",
            "percentile_95_lower",
            "percentile_95_upper",
        ):
            object.__setattr__(
                self,
                name,
                _finite(getattr(self, name), name=f"distribution {name}"),
            )


@dataclass(frozen=True, slots=True)
class CovarianceEntry:
    """One finite empirical covariance entry over the common sample set."""

    parameter_a: str
    parameter_b: str
    value: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "value",
            _finite(self.value, name="empirical covariance"),
        )


@dataclass(frozen=True, slots=True)
class CorrelationEntry:
    """One empirical correlation or an explicit typed unavailability."""

    parameter_a: str
    parameter_b: str
    availability: CorrelationAvailability
    value: float | None

    def __post_init__(self) -> None:
        if (self.availability is CorrelationAvailability.AVAILABLE) != (
            self.value is not None
        ):
            raise ResamplingConstructionError(
                "Correlation availability and value are inconsistent"
            )
        if self.value is not None:
            object.__setattr__(
                self,
                "value",
                _finite(self.value, name="empirical correlation"),
            )


@dataclass(frozen=True, slots=True)
class ResamplingSummary:
    """Qualified joint empirical summary with one inclusion denominator."""

    evidence_identity: str
    output_scope: tuple[str, ...]
    included_ordinals: tuple[int, ...]
    unstarted_ordinals: tuple[int, ...]
    exclusions: tuple[ResamplingExclusion, ...]
    samples: tuple[SummarySample, ...]
    distributions: tuple[ParameterDistribution, ...]
    covariance: tuple[CovarianceEntry, ...]
    correlations: tuple[CorrelationEntry, ...]
    policy_version: str = _SUMMARY_POLICY_VERSION
    claims: tuple[ClaimAssessment, ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        policy_identity = _identity(
            "native-resampling-summary-policy",
            (
                self.policy_version,
                self.evidence_identity,
                self.output_scope,
                self.included_ordinals,
                self.unstarted_ordinals,
            ),
        )
        claims = (
            ClaimAssessment(
                "COMMON_SCOPE_INCLUSION",
                ClaimState.SATISFIED,
                policy_identity,
                (
                    f"scope={','.join(self.output_scope)}",
                    f"included={len(self.included_ordinals)}",
                ),
            ),
            ClaimAssessment(
                "MINIMUM_SUCCESSFUL_COVERAGE",
                ClaimState.SATISFIED,
                policy_identity,
                (f"included={len(self.included_ordinals)}",),
            ),
            ClaimAssessment(
                "FINITE_EMPIRICAL_SUMMARY",
                ClaimState.SATISFIED,
                policy_identity,
            ),
        )
        object.__setattr__(self, "claims", claims)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-summary",
                (
                    self.evidence_identity,
                    self.output_scope,
                    self.included_ordinals,
                    self.unstarted_ordinals,
                    tuple(
                        (
                            item.ordinal,
                            item.outcome_identity,
                            _items_tokens(item.items),
                        )
                        for item in self.samples
                    ),
                    tuple(
                        (
                            item.ordinal,
                            item.outcome_identity,
                            item.disposition.value,
                            item.category,
                        )
                        for item in self.exclusions
                    ),
                    tuple(
                        (
                            item.parameter_id,
                            _float_token(item.mean),
                            _float_token(item.standard_deviation),
                            _float_token(item.median),
                            _float_token(item.percentile_95_lower),
                            _float_token(item.percentile_95_upper),
                        )
                        for item in self.distributions
                    ),
                    tuple(
                        (
                            item.parameter_a,
                            item.parameter_b,
                            _float_token(item.value),
                        )
                        for item in self.covariance
                    ),
                    tuple(
                        (
                            item.parameter_a,
                            item.parameter_b,
                            item.availability.value,
                            None if item.value is None else _float_token(item.value),
                        )
                        for item in self.correlations
                    ),
                    self.policy_version,
                    tuple(
                        (
                            claim.name,
                            claim.state.value,
                            claim.policy_identity,
                            claim.details,
                        )
                        for claim in claims
                    ),
                ),
            ),
        )

    def claim(self, name: str) -> ClaimState:
        """Return one exact named claim; unknown claims fail closed."""
        return _claim_state(self.claims, name)


@dataclass(frozen=True, slots=True)
class ResamplingSummaryOutcome:
    """Closed result union for one summary derivation request."""

    terminal: SummaryTerminal
    evidence_identity: str
    summary: ResamplingSummary | None = None
    failure: SummaryFailure | None = None

    def __post_init__(self) -> None:
        completed = self.terminal is SummaryTerminal.COMPLETED
        if completed != (self.summary is not None) or completed == (
            self.failure is not None
        ):
            raise ResamplingConstructionError(
                "Summary outcome must contain exactly its typed terminal payload"
            )


def _build_summary_statistics(
    matrix: Array,
    output_scope: tuple[str, ...],
) -> tuple[
    tuple[ParameterDistribution, ...],
    tuple[CovarianceEntry, ...],
    tuple[CorrelationEntry, ...],
]:
    distributions: list[ParameterDistribution] = []
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for index, param_id in enumerate(output_scope):
                values = matrix[:, index]
                lower, median, upper = np.percentile(values, (2.5, 50.0, 97.5))
                distributions.append(
                    ParameterDistribution(
                        param_id,
                        float(np.mean(values)),
                        float(np.std(values, ddof=1)),
                        float(median),
                        float(lower),
                        float(upper),
                    )
                )
            centered = matrix - np.mean(matrix, axis=0)
            covariance_matrix = (centered.T @ centered) / float(len(matrix) - 1)
    covariance = tuple(
        CovarianceEntry(param_a, param_b, float(covariance_matrix[row, column]))
        for row, param_a in enumerate(output_scope)
        for column, param_b in enumerate(output_scope)
    )
    variances = np.diag(covariance_matrix)
    correlations: list[CorrelationEntry] = []
    for row, param_a in enumerate(output_scope):
        for column, param_b in enumerate(output_scope):
            denominator = math.sqrt(float(variances[row])) * math.sqrt(
                float(variances[column])
            )
            if not math.isfinite(denominator):
                raise ResamplingConstructionError(
                    "Empirical correlation denominator must be finite"
                )
            if denominator == 0.0:
                correlations.append(
                    CorrelationEntry(
                        param_a,
                        param_b,
                        CorrelationAvailability.ZERO_VARIANCE,
                        None,
                    )
                )
            else:
                correlations.append(
                    CorrelationEntry(
                        param_a,
                        param_b,
                        CorrelationAvailability.AVAILABLE,
                        float(covariance_matrix[row, column] / denominator),
                    )
                )
    return tuple(distributions), covariance, tuple(correlations)


def summarize_resampling_evidence(
    evidence: ResamplingEvidence,
) -> ResamplingSummaryOutcome:
    """Derive one summary from all and only complete successful scope rows."""
    successful = tuple(
        outcome for outcome in evidence.outcomes if outcome.success is not None
    )
    if len(successful) < evidence.plan.minimum_successful_count:
        failure = SummaryFailure(
            "insufficient_successful_coverage",
            "Successful complete-scope replicate coverage is below the declared minimum",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.INSUFFICIENT_COVERAGE,
            evidence.identity,
            failure=failure,
        )
    samples = tuple(
        SummarySample(
            outcome.ordinal,
            outcome.identity,
            outcome.success.resolved_items,
        )
        for outcome in successful
        if outcome.success is not None
    )
    if any(
        tuple(param_id for param_id, _value in sample.items)
        != evidence.plan.output_scope
        for sample in samples
    ):
        failure = SummaryFailure(
            "source_scope_mismatch",
            "Successful source outcomes do not share the declared complete scope",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.SOURCE_INVALID,
            evidence.identity,
            failure=failure,
        )
    matrix = np.asarray(
        [[value for _param_id, value in sample.items] for sample in samples],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(matrix)):
        failure = SummaryFailure(
            "non_finite_source_sample",
            "A source sample contains a non-finite complete-scope value",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.SOURCE_INVALID,
            evidence.identity,
            failure=failure,
        )
    try:
        distributions, covariance, correlations = _build_summary_statistics(
            matrix,
            evidence.plan.output_scope,
        )
    except (FloatingPointError, RuntimeWarning, OverflowError, ValueError) as error:
        failure = SummaryFailure(
            "non_finite_summary_arithmetic",
            f"Summary arithmetic did not produce finite binary64 evidence: {error}",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.SOURCE_INVALID,
            evidence.identity,
            failure=failure,
        )
    exclusions = tuple(
        ResamplingExclusion(
            outcome.ordinal,
            outcome.identity,
            outcome.disposition,
            outcome.failure.category if outcome.failure is not None else "unknown",
        )
        for outcome in evidence.outcomes
        if outcome.failure is not None
    )
    summary = ResamplingSummary(
        evidence.identity,
        evidence.plan.output_scope,
        tuple(outcome.ordinal for outcome in successful),
        tuple(
            request.ordinal
            for request in evidence.plan.replicates
            if request.ordinal not in {outcome.ordinal for outcome in evidence.outcomes}
        ),
        exclusions,
        samples,
        distributions,
        covariance,
        correlations,
    )
    return ResamplingSummaryOutcome(
        SummaryTerminal.COMPLETED,
        evidence.identity,
        summary=summary,
    )
