"""Deterministic native resampling evidence qualification (#599).

This module plans identity-derived Monte Carlo, bootstrap, and
nucleus-bootstrap replicates; delegates each draw to an internally bound
evidence-only native optimizer; freshly evaluates every eligible candidate in a
separate replicate-local workspace; freezes every terminal outcome; and derives
summaries from one common complete-scope sample set. The production adapter
consumes the execution artifacts without exposing qualification policy objects
through the method-file surface.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field, replace
from enum import StrEnum
from numbers import Real
from typing import Literal, cast

import numpy as np

from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
    EvaluationPlan,
    EvaluationResult,
    ProfilePlan,
    ResampledProfileBinding,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    ComponentProblemDerivation,
    DirectTrfCandidateTerminal,
    DirectTrfInvocation,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
    canonical_chi_square,
    execute_direct_trf_candidate,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.runtime import ExecutionSettings
from chemex.runtime.execution import native_thread_environment
from chemex.typing import Array

_SCHEMA_VERSION = 1
_SEED_POLICY_VERSION = "sha256-product-stable-structural-u64-v2"
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


class ReplicateDisposition(StrEnum):
    """Terminal disposition of one actually executed replicate."""

    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    NOT_STARTED = "not_started"


class ReplicateStage(StrEnum):
    """Last truthful execution stage reached by a planned replicate."""

    PLANNED = "planned"
    GENERATING = "generating"
    DRAW_READY = "draw_ready"
    EXECUTOR_CREATING = "executor_creating"
    EXECUTING = "executing"
    VALIDATING = "validating"
    TERMINAL = "terminal"


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


def _evaluation_identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"schema": 1, "kind": kind, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
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
    source_dataset_identity: str
    scheme: ResamplingScheme
    generation_policy_version: str
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
                    self.source_dataset_identity,
                    self.scheme.value,
                    self.generation_policy_version,
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
    accepted_vector: tuple[float, ...]
    dataset: ResamplingDatasetManifest = field(repr=False, compare=False)
    source_problem: OptimizationProblem = field(repr=False, compare=False)
    parameterization: ActiveParameterization = field(repr=False, compare=False)
    source_engine: EvaluationEngine = field(repr=False, compare=False)
    scheme: ResamplingScheme
    replicate_count: int
    replicate_structural_identities: tuple[str, ...]
    replicate_component_identities: tuple[tuple[str, ...], ...]
    root_seed: int
    output_scope: tuple[str, ...]
    output_units: tuple[str, ...]
    minimum_successful_count: int
    strategy: OptimizationStrategy = OptimizationStrategy.DIRECT_TRF
    strategy_settings: tuple[tuple[str, str], ...] = ()
    seed_policy_version: str = _SEED_POLICY_VERSION
    rng_algorithm: str = _RNG_ALGORITHM
    generation_policy_version: str = _GENERATION_POLICY_VERSION
    identity: str = field(init=False)
    replicates: tuple[ReplicateRequest, ...] = field(init=False)
    source_dataset_identity: str = field(init=False)
    optimization_projection_identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - complete plan invariant
        count = _positive_integer(self.replicate_count, name="replicate count")
        accepted_vector = tuple(
            _finite(value, name=f"accepted vector[{index}]")
            for index, value in enumerate(self.accepted_vector)
        )
        if len(accepted_vector) != len(self.source_problem.controlled_ids):
            raise ResamplingConstructionError(
                "Accepted resampling vector must cover the controlled coordinates"
            )
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
        units = tuple(self.output_units)
        if len(units) != len(scope) or any(not unit for unit in units):
            raise ResamplingConstructionError(
                "Resampling output units must cover the exact ordered scope"
            )
        settings = tuple(sorted(self.strategy_settings))
        if any(not key or not value for key, value in settings) or len(
            {key for key, _value in settings}
        ) != len(settings):
            raise ResamplingConstructionError(
                "Optimization strategy settings require unique non-empty keys"
            )
        if (
            self.dataset.evaluation_plan.identity
            != self.source_problem.evaluation_plan_identity
            or self.source_problem.parameterization_identity
            != self.parameterization.identity
            or self.source_problem.evaluator_parameterization_identity
            != self.parameterization.evaluator_identity
            or self.source_problem.constraint_program_identity
            != self.parameterization.program.fingerprint
            or self.source_engine.plan.identity != self.dataset.evaluation_plan.identity
        ):
            raise ResamplingConstructionError(
                "Source dataset, problem, and parameterization lineage differ"
            )
        if tuple(scope) != self.source_problem.commit_scope:
            raise ResamplingConstructionError(
                "Resampling output scope must equal the exact source problem scope"
            )
        projection_identity = _identity(
            "native-resampling-optimization-projection",
            (
                self.source_problem.identity,
                self.parameterization.identity,
                self.strategy.value,
                settings,
                component_identities,
                scope,
                units,
            ),
        )
        for name, value in (
            ("accepted result identity", self.accepted_result_identity),
            ("accepted occurrence identity", self.accepted_occurrence_identity),
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
                tuple(_float_token(value) for value in accepted_vector),
                self.scheme.value,
                count,
                structural_identities,
                component_identities,
                root_seed,
                self.dataset.identity,
                self.source_problem.identity,
                self.parameterization.identity,
                scope,
                units,
                projection_identity,
                minimum,
                self.strategy.value,
                settings,
                self.seed_policy_version,
                self.rng_algorithm,
                self.generation_policy_version,
            ),
        )
        seed_stage_identity = _identity(
            "native-resampling-seed-stage",
            (
                self.scheme.value,
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
                stage_identity=seed_stage_identity,
                work_unit_kind="resampling-replicate",
                structural_identity=structural_identity,
            )
            requests.append(
                ReplicateRequest(
                    identity,
                    self.dataset.identity,
                    self.scheme,
                    self.generation_policy_version,
                    ordinal,
                    structural_identity,
                    projection_identity,
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
        object.__setattr__(self, "accepted_vector", accepted_vector)
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
        object.__setattr__(self, "output_units", units)
        object.__setattr__(self, "strategy_settings", settings)
        object.__setattr__(self, "source_dataset_identity", self.dataset.identity)
        object.__setattr__(
            self, "optimization_projection_identity", projection_identity
        )
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "replicates", tuple(requests))

    @classmethod
    def for_accepted(
        cls,
        accepted: AcceptedFitResult,
        *,
        dataset: ResamplingDatasetManifest,
        source_problem: OptimizationProblem,
        parameterization: ActiveParameterization,
        source_engine: EvaluationEngine,
        scheme: ResamplingScheme,
        replicate_count: int,
        replicate_structural_identities: Sequence[str],
        replicate_component_identities: Sequence[Sequence[str]],
        root_seed: int,
        output_scope: Sequence[str],
        output_units: Sequence[str],
        minimum_successful_count: int,
        strategy: OptimizationStrategy = OptimizationStrategy.DIRECT_TRF,
        strategy_settings: Sequence[tuple[str, str]] = (),
    ) -> ResamplingPlan:
        """Bind a deterministic replicate population to one accepted occurrence."""
        if not accepted_occurrence_is_authoritative(accepted):
            raise ResamplingConstructionError(
                "Resampling planning requires an exact authoritative accepted occurrence"
            )
        _validate_accepted_source(
            accepted,
            dataset,
            source_problem,
            parameterization,
        )
        if source_engine.plan.identity != dataset.evaluation_plan.identity:
            raise ResamplingConstructionError(
                "Source evaluator belongs to another EvaluationPlan"
            )
        return cls(
            accepted.identity,
            accepted.occurrence_identity,
            accepted.vector,
            dataset,
            source_problem,
            parameterization,
            source_engine,
            scheme,
            replicate_count,
            tuple(replicate_structural_identities),
            tuple(tuple(items) for items in replicate_component_identities),
            root_seed,
            tuple(output_scope),
            tuple(output_units),
            minimum_successful_count,
            strategy,
            tuple(strategy_settings),
        )


def _validate_accepted_source(
    accepted: AcceptedFitResult,
    dataset: ResamplingDatasetManifest,
    source_problem: OptimizationProblem,
    parameterization: ActiveParameterization,
) -> None:
    snapshot = source_problem.source_snapshot
    if (
        accepted.problem_identity != source_problem.identity
        or accepted.parameterization_identity != parameterization.identity
        or accepted.evaluator_parameterization_identity
        != parameterization.evaluator_identity
        or accepted.source_occurrence_identity != snapshot.occurrence_identity
        or accepted.source_revision != snapshot.revision
        or accepted.evaluation_result.plan_identity != dataset.evaluation_plan.identity
        or accepted.evaluation_result.parameterization_identity
        != parameterization.evaluator_identity
        or tuple(
            float(value) for value in accepted.evaluation_result.normalized_calculations
        )
        != dataset.calculated
        or accepted.controlled_ids != source_problem.controlled_ids
        or accepted.commit_scope != source_problem.commit_scope
        or accepted.commit_items
        != accepted.evaluation_result.resolved_values.ordered_items()
    ):
        raise ResamplingConstructionError(
            "Resampling source does not match the exact accepted native lineage"
        )


@dataclass(frozen=True, slots=True)
class ResamplingDatasetManifest:
    """Canonical by-value source population owned by one resampling request."""

    evaluation_plan: EvaluationPlan
    calculated: tuple[float, ...]
    references: tuple[bool, ...]
    nucleus_groups: tuple[str, ...]
    observation_descriptors: tuple[str, ...]
    identity: str = field(init=False)

    @property
    def observed(self) -> tuple[float, ...]:
        return tuple(
            value
            for profile in self.evaluation_plan.profiles
            for value in profile.experimental_values
        )

    @property
    def standard_errors(self) -> tuple[float, ...]:
        return tuple(
            value
            for profile in self.evaluation_plan.profiles
            for value in profile.uncertainties
        )

    @property
    def mask(self) -> tuple[bool, ...]:
        return tuple(
            value for profile in self.evaluation_plan.profiles for value in profile.mask
        )

    @property
    def profile_blocks(self) -> tuple[str, ...]:
        return tuple(
            profile.identity
            for profile in self.evaluation_plan.profiles
            for _ in range(profile.observation_count)
        )

    def __post_init__(self) -> None:
        try:
            evaluation_plan = EvaluationPlan.from_record(
                self.evaluation_plan.to_record()
            )
        except (TypeError, ValueError) as error:
            raise ResamplingConstructionError(
                "Source dataset requires a canonical EvaluationPlan descriptor"
            ) from error
        size = evaluation_plan.observation_count
        if size < 1 or any(
            len(values) != size
            for values in (
                self.calculated,
                self.references,
                self.nucleus_groups,
                self.observation_descriptors,
            )
        ):
            raise ResamplingConstructionError(
                "Source dataset fields must match the EvaluationPlan population"
            )
        calculated = tuple(
            _finite(value, name=f"calculated[{index}]")
            for index, value in enumerate(self.calculated)
        )
        if any(type(value) is not bool for value in self.references):
            raise ResamplingConstructionError("Reference flags must be Boolean")
        if any(value < 0.0 for value in self.standard_errors):
            raise ResamplingConstructionError(
                "Source dataset uncertainties must be non-negative"
            )
        if any(not group for group in self.nucleus_groups):
            raise ResamplingConstructionError(
                "Nucleus group identities must not be empty"
            )
        if any(not descriptor for descriptor in self.observation_descriptors):
            raise ResamplingConstructionError(
                "Observation descriptors must not be empty"
            )
        object.__setattr__(self, "evaluation_plan", evaluation_plan)
        object.__setattr__(self, "calculated", calculated)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-source-dataset",
                (
                    self.evaluation_plan.identity,
                    tuple(_float_token(value) for value in calculated),
                    self.references,
                    self.nucleus_groups,
                    self.observation_descriptors,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        """Return a canonical record whose digest can be independently verified."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "evaluation_plan": self.evaluation_plan.to_record(),
            "calculated": [_float_token(value) for value in self.calculated],
            "references": list(self.references),
            "nucleus_groups": list(self.nucleus_groups),
            "observation_descriptors": list(self.observation_descriptors),
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ResamplingDatasetManifest:
        """Restore a source manifest only after recomputing its exact identity."""
        expected_keys = {
            "schema_version",
            "evaluation_plan",
            "calculated",
            "references",
            "nucleus_groups",
            "observation_descriptors",
            "identity",
        }
        if (
            set(record) != expected_keys
            or record.get("schema_version") != _SCHEMA_VERSION
        ):
            raise ValueError("Malformed native resampling dataset record")
        plan_record = record.get("evaluation_plan")
        if not isinstance(plan_record, Mapping):
            raise TypeError("Resampling dataset EvaluationPlan must be a mapping")
        try:
            calculated_record = cast("list[object]", record["calculated"])
            references_record = cast("list[bool]", record["references"])
            groups_record = cast("list[object]", record["nucleus_groups"])
            descriptors_record = cast("list[object]", record["observation_descriptors"])
            if not all(
                isinstance(value, list)
                for value in (
                    calculated_record,
                    references_record,
                    groups_record,
                    descriptors_record,
                )
            ):
                raise TypeError("Resampling dataset vectors must be lists")
            manifest = cls(
                EvaluationPlan.from_record(cast("Mapping[str, object]", plan_record)),
                tuple(float.fromhex(str(value)) for value in calculated_record),
                tuple(references_record),
                tuple(str(value) for value in groups_record),
                tuple(str(value) for value in descriptors_record),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Malformed native resampling dataset payload") from error
        if record.get("identity") != manifest.identity:
            raise ValueError("Resampling dataset identity does not match its content")
        return manifest


_DRAW_CONSTRUCTION_KEY = object()


@dataclass(frozen=True, slots=True)
class ResamplingDraw:
    """One immutable seeded transformation with exact source selection."""

    _construction_key: object = field(repr=False, compare=False)
    request_identity: str
    ordinal: int
    source_dataset_identity: str
    scheme: ResamplingScheme
    generation_policy_version: str
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
        if self._construction_key is not _DRAW_CONSTRUCTION_KEY:
            raise ResamplingConstructionError(
                "Resampling draws must be produced by canonical generation"
            )
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
                    self.request_identity,
                    self.ordinal,
                    self.source_dataset_identity,
                    self.scheme.value,
                    self.generation_policy_version,
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
    dataset: ResamplingDatasetManifest,
    request: ReplicateRequest,
) -> ResamplingDraw:
    """Apply one retained scientific generation scheme with explicit PCG64 state."""
    if request.source_dataset_identity != dataset.identity:
        raise ResamplingConstructionError(
            "Replicate request belongs to another source dataset"
        )
    normalized_seed = _unsigned_seed(request.seed, name="draw seed")
    rng = np.random.Generator(np.random.PCG64(normalized_seed))
    size = len(dataset.observed)
    if request.scheme is ResamplingScheme.MONTE_CARLO:
        observations = tuple(
            float(value)
            for value in rng.normal(dataset.calculated, dataset.standard_errors)
        )
        return ResamplingDraw(
            _DRAW_CONSTRUCTION_KEY,
            request.identity,
            request.ordinal,
            dataset.identity,
            request.scheme,
            request.generation_policy_version,
            normalized_seed,
            observations,
            dataset.standard_errors,
            dataset.mask,
            dataset.references,
            dataset.nucleus_groups,
            dataset.profile_blocks,
            tuple(range(size)),
        )
    if request.scheme is ResamplingScheme.BOOTSTRAP:
        source_indices = list(range(size))
        for block in dict.fromkeys(dataset.profile_blocks):
            masked = tuple(
                index
                for index, active in enumerate(dataset.mask)
                if active and dataset.profile_blocks[index] == block
            )
            reference_pool = tuple(
                index for index in masked if dataset.references[index]
            )
            non_reference_pool = tuple(
                index for index in masked if not dataset.references[index]
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
            _DRAW_CONSTRUCTION_KEY,
            request.identity,
            request.ordinal,
            dataset.identity,
            request.scheme,
            request.generation_policy_version,
            normalized_seed,
            tuple(dataset.observed[index] for index in indexes),
            tuple(dataset.standard_errors[index] for index in indexes),
            dataset.mask,
            tuple(dataset.references[index] for index in indexes),
            dataset.nucleus_groups,
            dataset.profile_blocks,
            indexes,
        )
    canonical_groups = tuple(sorted(set(dataset.nucleus_groups)))
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
        for index, source_group in enumerate(dataset.nucleus_groups)
        if source_group == group
    )
    return ResamplingDraw(
        _DRAW_CONSTRUCTION_KEY,
        request.identity,
        request.ordinal,
        dataset.identity,
        request.scheme,
        request.generation_policy_version,
        normalized_seed,
        tuple(dataset.observed[index] for index in indexes),
        tuple(dataset.standard_errors[index] for index in indexes),
        tuple(dataset.mask[index] for index in indexes),
        tuple(dataset.references[index] for index in indexes),
        tuple(dataset.nucleus_groups[index] for index in indexes),
        tuple(dataset.profile_blocks[index] for index in indexes),
        indexes,
        sampled_groups,
    )


def _transformed_evaluation_plan(
    dataset: ResamplingDatasetManifest,
    draw: ResamplingDraw,
    source_engine: EvaluationEngine,
) -> tuple[EvaluationPlan, tuple[ResampledProfileBinding, ...]]:
    source_profiles = {
        profile.identity: (profile_index, profile)
        for profile_index, profile in enumerate(dataset.evaluation_plan.profiles)
    }
    profiles: list[ProfilePlan] = []
    bindings: list[ResampledProfileBinding] = []
    offset = 0
    index = 0
    while index < len(draw.observations):
        block_identity = draw.profile_blocks[index]
        stop = index + 1
        while stop < len(draw.observations) and (
            draw.profile_blocks[stop] == block_identity
        ):
            stop += 1
        source_profile_index, source = source_profiles[block_identity]
        root_indices = tuple(
            source_index - source.observation_offset
            for source_index in draw.source_indices[index:stop]
        )
        binding = ResampledProfileBinding(source_profile_index, root_indices)
        source_identity = _identity(
            "native-resampling-transformed-profile-source",
            (
                dataset.identity,
                draw.identity,
                source.identity,
                len(profiles),
                tuple(draw.source_indices[index:stop]),
                tuple(_float_token(value) for value in draw.observations[index:stop]),
                tuple(
                    _float_token(value) for value in draw.standard_errors[index:stop]
                ),
                tuple(draw.mask[index:stop]),
                tuple(draw.references[index:stop]),
                tuple(draw.nucleus_groups[index:stop]),
            ),
        )
        profile_ordinal = len(profiles)
        profile_identity = _evaluation_identity(
            "profile-plan",
            (source_identity, source.experiment_ordinal, profile_ordinal, offset),
        )
        profiles.append(
            ProfilePlan(
                profile_identity,
                source_identity,
                source.experiment_ordinal,
                profile_ordinal,
                offset,
                source.local_inputs,
                source.is_scaled,
                draw.observations[index:stop],
                draw.standard_errors[index:stop],
                draw.mask[index:stop],
                source.kernel_identity,
                source.kernel_configuration,
                source.spectrometer_configuration,
                source_engine.resampled_observation_metadata(binding),
                (stop - index,),
            )
        )
        bindings.append(binding)
        offset += stop - index
        index = stop
    source_plan = dataset.evaluation_plan
    return EvaluationPlan(
        source_plan.parameterization_identity,
        source_plan.constraint_program_identity,
        tuple(profiles),
        source_plan.resolved_ids,
        source_plan.scalar_version,
        source_plan.normalization_version,
        source_plan.residual_version,
        source_plan.masking_version,
        source_plan.ordering_version,
        source_plan.failure_version,
        source_plan.diagnostics_version,
    ), tuple(bindings)


@dataclass(frozen=True, slots=True)
class ReplicateExecutionPlan:
    """Exact non-authoritative native problem prepared for one generated draw."""

    plan_identity: str
    request: ReplicateRequest
    draw: ResamplingDraw
    accepted_result_identity: str
    accepted_occurrence_identity: str
    source_snapshot_occurrence_identity: str
    source_snapshot_revision: int
    evaluation_plan: EvaluationPlan = field(repr=False, compare=False)
    profile_bindings: tuple[ResampledProfileBinding, ...] = field(
        repr=False, compare=False
    )
    problem: OptimizationProblem = field(repr=False, compare=False)
    parameterization: ActiveParameterization = field(repr=False, compare=False)
    invocation: DirectTrfInvocation = field(repr=False, compare=False)
    strategy: OptimizationStrategy
    workflow_identity: str
    invocation_identity: str
    component_identities: tuple[str, ...]
    output_scope: tuple[str, ...]
    output_units: tuple[str, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            self.draw.request_identity != self.request.identity
            or self.draw.source_dataset_identity != self.request.source_dataset_identity
            or self.evaluation_plan.identity != self.problem.evaluation_plan_identity
            or self.problem.parameterization_identity != self.parameterization.identity
            or self.problem.evaluator_parameterization_identity
            != self.parameterization.evaluator_identity
            or self.problem.constraint_program_identity
            != self.parameterization.program.fingerprint
            or self.invocation.problem_identity != self.problem.identity
            or self.invocation.identity != self.invocation_identity
            or self.problem.acceptance_authority
            or self.problem.commit_scope != self.output_scope
            or len(self.output_units) != len(self.output_scope)
            or len(self.profile_bindings) != len(self.evaluation_plan.profiles)
        ):
            raise ResamplingConstructionError(
                "Prepared replicate native lineage is internally inconsistent"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-prepared-replicate",
                (
                    self.plan_identity,
                    self.request.identity,
                    self.draw.identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_snapshot_occurrence_identity,
                    self.source_snapshot_revision,
                    self.evaluation_plan.identity,
                    tuple(
                        (
                            binding.root_profile_index,
                            binding.root_observation_indices,
                        )
                        for binding in self.profile_bindings
                    ),
                    self.problem.identity,
                    self.parameterization.identity,
                    self.parameterization.evaluator_identity,
                    self.parameterization.program.fingerprint,
                    self.strategy.value,
                    self.workflow_identity,
                    self.invocation_identity,
                    self.component_identities,
                    self.output_scope,
                    self.output_units,
                ),
            ),
        )

    @classmethod
    def prepare(
        cls,
        plan: ResamplingPlan,
        request: ReplicateRequest,
        draw: ResamplingDraw,
    ) -> ReplicateExecutionPlan:
        """Derive every projected identity before executor construction."""
        if (
            request.plan_identity != plan.identity
            or draw.request_identity != request.identity
        ):
            raise ResamplingConstructionError(
                "Cannot prepare a replicate from foreign request or draw lineage"
            )
        if draw != generate_resampling_draw(plan.dataset, request):
            raise ResamplingConstructionError(
                "Prepared replicate draw differs from canonical seeded generation"
            )
        evaluation_plan, profile_bindings = _transformed_evaluation_plan(
            plan.dataset, draw, plan.source_engine
        )
        derivation = ComponentProblemDerivation(
            plan.source_problem.identity,
            plan.source_problem.affine_feasibility_identity,
            request.structural_identity,
            "native-resampling-complete-scope-v1",
            evaluation_plan.identity,
            plan.source_problem.controlled_ids,
            plan.source_problem.held_items,
        )
        problem = plan.source_problem.derive_child(
            controlled_ids=plan.source_problem.controlled_ids,
            start=plan.accepted_vector,
            derivation=derivation,
        )
        if plan.strategy is not OptimizationStrategy.DIRECT_TRF:
            raise ResamplingConstructionError(
                "Qualified resampling currently supports the closed Direct-TRF binder"
            )
        settings = dict(plan.strategy_settings)
        if set(settings) - {"objective_request_budget"}:
            raise ResamplingConstructionError(
                "Direct-TRF resampling received unsupported strategy settings"
            )
        try:
            objective_request_budget = int(
                settings.get("objective_request_budget", "100")
            )
        except ValueError as error:
            raise ResamplingConstructionError(
                "Direct-TRF objective request budget must be an integer"
            ) from error
        invocation = DirectTrfInvocation.for_problem(
            problem,
            objective_request_budget=objective_request_budget,
        )
        workflow_identity = _identity(
            "native-resampling-projected-workflow",
            (
                plan.optimization_projection_identity,
                plan.strategy.value,
                plan.strategy_settings,
                problem.identity,
                request.component_identities,
            ),
        )
        invocation_identity = invocation.identity
        snapshot = problem.source_snapshot
        return cls(
            plan.identity,
            request,
            draw,
            plan.accepted_result_identity,
            plan.accepted_occurrence_identity,
            snapshot.occurrence_identity,
            snapshot.revision,
            evaluation_plan,
            profile_bindings,
            problem,
            plan.parameterization,
            invocation,
            plan.strategy,
            workflow_identity,
            invocation_identity,
            request.component_identities,
            plan.output_scope,
            plan.output_units,
        )

    def resolve_candidate(
        self, candidate_vector: Sequence[float]
    ) -> tuple[tuple[str, float], ...]:
        """Freshly apply the exact problem feasibility and constraint program."""
        frame = self.problem.lifecycle_frame(candidate_vector, self.parameterization)
        resolved = self.parameterization.resolve(frame)
        return tuple((param_id, resolved[param_id]) for param_id in self.output_scope)


@dataclass(frozen=True, slots=True)
class ProjectedOptimizationSuccess:
    """Non-authoritative optimizer candidate with exact prepared lineage."""

    prepared_replicate_identity: str
    accepted_result_identity: str
    accepted_occurrence_identity: str
    source_snapshot_occurrence_identity: str
    source_snapshot_revision: int
    transformed_data_identity: str
    optimization_projection_identity: str
    evaluation_plan_identity: str
    problem_identity: str
    parameterization_identity: str
    evaluator_parameterization_identity: str
    constraint_program_identity: str
    workflow_identity: str
    invocation_identity: str
    execution_identity: str
    strategy: OptimizationStrategy
    component_identities: tuple[str, ...]
    component_outcome_identities: tuple[str, ...]
    controlled_ids: tuple[str, ...]
    candidate_identity: str
    candidate_vector: tuple[float, ...]
    backend_chi_square: float | None
    backend_evaluation_identity: str | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for name, value in (
            ("prepared replicate identity", self.prepared_replicate_identity),
            ("accepted result identity", self.accepted_result_identity),
            ("accepted occurrence identity", self.accepted_occurrence_identity),
            (
                "source snapshot occurrence identity",
                self.source_snapshot_occurrence_identity,
            ),
            ("transformed data identity", self.transformed_data_identity),
            ("optimization projection identity", self.optimization_projection_identity),
            ("evaluation plan identity", self.evaluation_plan_identity),
            ("problem identity", self.problem_identity),
            ("parameterization identity", self.parameterization_identity),
            (
                "evaluator parameterization identity",
                self.evaluator_parameterization_identity,
            ),
            ("constraint program identity", self.constraint_program_identity),
            ("workflow identity", self.workflow_identity),
            ("invocation identity", self.invocation_identity),
            ("execution identity", self.execution_identity),
            ("candidate identity", self.candidate_identity),
        ):
            _non_empty_identity(value, name=name)
        if self.source_snapshot_revision < 0:
            raise ResamplingConstructionError(
                "Source snapshot revision must be non-negative"
            )
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
        vector = tuple(
            _finite(value, name=f"candidate vector[{index}]")
            for index, value in enumerate(self.candidate_vector)
        )
        controlled_ids = _canonical_scope(self.controlled_ids)
        if len(vector) != len(controlled_ids):
            raise ResamplingConstructionError(
                "Projected candidate must cover the exact controlled-coordinate order"
            )
        backend_chi_square = (
            None
            if self.backend_chi_square is None
            else _finite(self.backend_chi_square, name="backend diagnostic chi-square")
        )
        if backend_chi_square is not None and backend_chi_square < 0.0:
            raise ResamplingConstructionError(
                "Backend diagnostic chi-square must be non-negative"
            )
        object.__setattr__(self, "candidate_vector", vector)
        object.__setattr__(self, "controlled_ids", controlled_ids)
        object.__setattr__(self, "backend_chi_square", backend_chi_square)
        if self.backend_evaluation_identity is not None:
            _non_empty_identity(
                self.backend_evaluation_identity,
                name="backend evaluation diagnostic identity",
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-projected-optimization-success",
                (
                    self.prepared_replicate_identity,
                    self.accepted_result_identity,
                    self.accepted_occurrence_identity,
                    self.source_snapshot_occurrence_identity,
                    self.source_snapshot_revision,
                    self.transformed_data_identity,
                    self.optimization_projection_identity,
                    self.evaluation_plan_identity,
                    self.problem_identity,
                    self.parameterization_identity,
                    self.evaluator_parameterization_identity,
                    self.constraint_program_identity,
                    self.workflow_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.strategy.value,
                    self.component_identities,
                    self.component_outcome_identities,
                    controlled_ids,
                    self.candidate_identity,
                    tuple(_float_token(value) for value in vector),
                    (
                        None
                        if backend_chi_square is None
                        else _float_token(backend_chi_square)
                    ),
                    self.backend_evaluation_identity,
                ),
            ),
        )

    @classmethod
    def for_prepared(
        cls,
        prepared: ReplicateExecutionPlan,
        *,
        candidate_vector: Sequence[float],
        backend_chi_square: float | None,
        execution_identity: str,
        backend_evaluation_identity: str | None = None,
    ) -> ProjectedOptimizationSuccess:
        """Bind non-authoritative backend evidence to an exact prepared replicate."""
        vector = tuple(float(value) for value in candidate_vector)
        candidate_identity = _identity(
            "native-resampling-candidate",
            (
                prepared.identity,
                prepared.problem.controlled_ids,
                tuple(_float_token(value) for value in vector),
            ),
        )
        component_outcomes = tuple(
            _identity(
                "native-resampling-component-outcome",
                (prepared.identity, component_identity, candidate_identity),
            )
            for component_identity in prepared.component_identities
        )
        return cls(
            prepared.identity,
            prepared.accepted_result_identity,
            prepared.accepted_occurrence_identity,
            prepared.source_snapshot_occurrence_identity,
            prepared.source_snapshot_revision,
            prepared.draw.identity,
            prepared.request.optimization_projection_identity,
            prepared.evaluation_plan.identity,
            prepared.problem.identity,
            prepared.parameterization.identity,
            prepared.parameterization.evaluator_identity,
            prepared.parameterization.program.fingerprint,
            prepared.workflow_identity,
            prepared.invocation_identity,
            execution_identity,
            prepared.strategy,
            prepared.component_identities,
            component_outcomes,
            prepared.problem.controlled_ids,
            candidate_identity,
            vector,
            backend_chi_square,
            backend_evaluation_identity,
        )


@dataclass(frozen=True, slots=True)
class ProjectedOptimizationFailure:
    """Typed unsuccessful projected path retaining all available lineage."""

    transformed_data_identity: str | None
    disposition: ReplicateDisposition
    category: str
    message: str
    optimization_projection_identity: str = field(kw_only=True)
    stage: ReplicateStage = ReplicateStage.EXECUTING
    prepared_replicate_identity: str | None = None
    evaluation_plan_identity: str | None = None
    problem_identity: str | None = None
    invocation_identity: str | None = None
    execution_identity: str | None = None
    component_outcome_identities: tuple[str, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.disposition is ReplicateDisposition.SUCCEEDED:
            raise ResamplingConstructionError(
                "Projected failure cannot use the successful disposition"
            )
        if not self.category or not self.message:
            raise ResamplingConstructionError(
                "Projected failure requires category and message"
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
                    self.stage.value,
                    self.prepared_replicate_identity,
                    self.evaluation_plan_identity,
                    self.problem_identity,
                    self.invocation_identity,
                    self.execution_identity,
                    self.component_outcome_identities,
                ),
            ),
        )


type ProjectedOptimizationResult = (
    ProjectedOptimizationSuccess | ProjectedOptimizationFailure
)


type CandidateEvidenceTestHook = Callable[
    [ReplicateExecutionPlan, ProjectedOptimizationResult],
    ProjectedOptimizationResult,
]


@dataclass(frozen=True, slots=True, init=False)
class ValidatedReplicateSuccess:
    """Canonical complete sample derived from exact fresh materialization sources."""

    prepared: ReplicateExecutionPlan = field(init=False, repr=False, compare=False)
    projected: ProjectedOptimizationSuccess = field(
        init=False, repr=False, compare=False
    )
    evaluation_result: EvaluationResult = field(init=False, repr=False, compare=False)

    @classmethod
    def from_materialization(
        cls,
        prepared: ReplicateExecutionPlan,
        projected: ProjectedOptimizationSuccess,
        evaluation_result: EvaluationResult,
    ) -> ValidatedReplicateSuccess:
        """Construct one success solely from its exact subordinate evidence."""
        try:
            rebound_prepared = replace(prepared)
            canonical_projected = ProjectedOptimizationSuccess.for_prepared(
                rebound_prepared,
                candidate_vector=projected.candidate_vector,
                backend_chi_square=projected.backend_chi_square,
                execution_identity=projected.execution_identity,
                backend_evaluation_identity=(projected.backend_evaluation_identity),
            )
            canonical_evaluation = EvaluationResult.from_record(
                evaluation_result.to_record(), rebound_prepared.evaluation_plan
            )
            resolved_items = tuple(
                (param_id, canonical_evaluation.resolved_values[param_id])
                for param_id in rebound_prepared.output_scope
            )
            if (
                rebound_prepared.identity != prepared.identity
                or canonical_projected != projected
                or canonical_projected.identity != projected.identity
                or canonical_evaluation.identity != evaluation_result.identity
                or canonical_evaluation.plan_identity
                != rebound_prepared.evaluation_plan.identity
                or canonical_evaluation.parameterization_identity
                != rebound_prepared.parameterization.evaluator_identity
                or tuple(canonical_evaluation.resolved_values)
                != rebound_prepared.output_scope
                or resolved_items
                != rebound_prepared.resolve_candidate(projected.candidate_vector)
                or (
                    projected.backend_evaluation_identity is not None
                    and projected.backend_evaluation_identity
                    != canonical_evaluation.identity
                )
            ):
                raise ResamplingConstructionError(
                    "Fresh materialization differs from its exact prepared candidate"
                )
            canonical_chi_square(canonical_evaluation.residuals)
        except ResamplingConstructionError:
            raise
        except Exception as error:
            raise ResamplingConstructionError(
                "Successful replicate has an invalid fresh derivation chain"
            ) from error
        success = object.__new__(cls)
        object.__setattr__(success, "prepared", rebound_prepared)
        object.__setattr__(success, "projected", canonical_projected)
        object.__setattr__(success, "evaluation_result", canonical_evaluation)
        return success

    @property
    def prepared_replicate_identity(self) -> str:
        return self.prepared.identity

    @property
    def projected_outcome_identity(self) -> str:
        return self.projected.identity

    @property
    def evaluation_identity(self) -> str:
        return self.evaluation_result.identity

    @property
    def output_scope(self) -> tuple[str, ...]:
        return self.prepared.output_scope

    @property
    def output_units(self) -> tuple[str, ...]:
        return self.prepared.output_units

    @property
    def resolved_items(self) -> tuple[tuple[str, float], ...]:
        return tuple(
            (param_id, self.evaluation_result.resolved_values[param_id])
            for param_id in self.output_scope
        )

    @property
    def chi_square(self) -> float:
        return canonical_chi_square(self.evaluation_result.residuals)

    @property
    def backend_chi_square(self) -> float | None:
        return self.projected.backend_chi_square

    @property
    def backend_chi_square_agrees(self) -> bool | None:
        return (
            None
            if self.backend_chi_square is None
            else self.backend_chi_square == self.chi_square
        )

    @property
    def fresh_evaluation_count(self) -> int:
        return 1

    @property
    def materialization_identity(self) -> str:
        return _identity(
            "native-resampling-fresh-materialization",
            (
                self.prepared.identity,
                self.projected.candidate_identity,
                self.projected.execution_identity,
                self.evaluation_result.evaluator_compatibility_identity,
                self.evaluation_result.identity,
                _items_tokens(self.resolved_items),
            ),
        )

    @property
    def claims(self) -> tuple[ClaimAssessment, ...]:
        return (
            ClaimAssessment(
                "FRESH_MATERIALIZATION_INTEGRITY",
                ClaimState.SATISFIED,
                self.prepared.identity,
                (
                    f"evaluation={self.evaluation_identity}",
                    f"scope={','.join(self.output_scope)}",
                ),
            ),
        )

    @property
    def identity(self) -> str:
        return _identity(
            "native-resampling-validated-success",
            (
                self.prepared_replicate_identity,
                self.projected_outcome_identity,
                self.evaluation_identity,
                self.materialization_identity,
                _items_tokens(self.resolved_items),
                _float_token(self.chi_square),
                (
                    None
                    if self.backend_chi_square is None
                    else _float_token(self.backend_chi_square)
                ),
                self.backend_chi_square_agrees,
                self.fresh_evaluation_count,
                self.output_scope,
                self.output_units,
                tuple(
                    (
                        claim.name,
                        claim.state.value,
                        claim.policy_identity,
                        claim.details,
                    )
                    for claim in self.claims
                ),
            ),
        )

    @staticmethod
    def _projected_record(
        projected: ProjectedOptimizationSuccess,
    ) -> dict[str, object]:
        return {
            "prepared_replicate_identity": projected.prepared_replicate_identity,
            "accepted_result_identity": projected.accepted_result_identity,
            "accepted_occurrence_identity": projected.accepted_occurrence_identity,
            "source_snapshot_occurrence_identity": (
                projected.source_snapshot_occurrence_identity
            ),
            "source_snapshot_revision": projected.source_snapshot_revision,
            "transformed_data_identity": projected.transformed_data_identity,
            "optimization_projection_identity": (
                projected.optimization_projection_identity
            ),
            "evaluation_plan_identity": projected.evaluation_plan_identity,
            "problem_identity": projected.problem_identity,
            "parameterization_identity": projected.parameterization_identity,
            "evaluator_parameterization_identity": (
                projected.evaluator_parameterization_identity
            ),
            "constraint_program_identity": projected.constraint_program_identity,
            "workflow_identity": projected.workflow_identity,
            "invocation_identity": projected.invocation_identity,
            "execution_identity": projected.execution_identity,
            "strategy": projected.strategy.value,
            "component_identities": list(projected.component_identities),
            "component_outcome_identities": list(
                projected.component_outcome_identities
            ),
            "controlled_ids": list(projected.controlled_ids),
            "candidate_identity": projected.candidate_identity,
            "candidate_vector": [
                _float_token(value) for value in projected.candidate_vector
            ],
            "backend_chi_square": (
                None
                if projected.backend_chi_square is None
                else _float_token(projected.backend_chi_square)
            ),
            "backend_evaluation_identity": projected.backend_evaluation_identity,
            "identity": projected.identity,
        }

    def to_record(self) -> dict[str, object]:
        """Serialize sources plus checked readable projections."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "prepared_replicate_identity": self.prepared_replicate_identity,
            "request_identity": self.prepared.request.identity,
            "draw_identity": self.prepared.draw.identity,
            "evaluation_plan_identity": self.prepared.evaluation_plan.identity,
            "problem_identity": self.prepared.problem.identity,
            "projected": self._projected_record(self.projected),
            "evaluation_result": self.evaluation_result.to_record(),
            "output_scope": list(self.output_scope),
            "output_units": list(self.output_units),
            "resolved_items": [
                [param_id, _float_token(value)]
                for param_id, value in self.resolved_items
            ],
            "chi_square": _float_token(self.chi_square),
            "backend_chi_square": (
                None
                if self.backend_chi_square is None
                else _float_token(self.backend_chi_square)
            ),
            "backend_chi_square_agrees": self.backend_chi_square_agrees,
            "fresh_evaluation_count": self.fresh_evaluation_count,
            "materialization_identity": self.materialization_identity,
            "claims": [
                {
                    "name": claim.name,
                    "state": claim.state.value,
                    "policy_identity": claim.policy_identity,
                    "details": list(claim.details),
                }
                for claim in self.claims
            ],
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        plan: ResamplingPlan,
        request: ReplicateRequest,
    ) -> ValidatedReplicateSuccess:
        """Rebuild from exact source context and reject altered checked copies."""
        try:
            draw = generate_resampling_draw(plan.dataset, request)
            prepared = ReplicateExecutionPlan.prepare(plan, request, draw)
            projected_record = record.get("projected")
            evaluation_record = record.get("evaluation_result")
            if not isinstance(projected_record, Mapping) or not isinstance(
                evaluation_record, Mapping
            ):
                raise TypeError("Success source records must be mappings")
            vector_record = projected_record.get("candidate_vector")
            if not isinstance(vector_record, list):
                raise TypeError("Candidate vector must be a list")
            backend_record = projected_record.get("backend_chi_square")
            backend_chi_square = (
                None if backend_record is None else float.fromhex(str(backend_record))
            )
            execution_identity = projected_record.get("execution_identity")
            if not isinstance(execution_identity, str):
                raise TypeError("Execution identity must be a string")
            backend_evaluation_identity = projected_record.get(
                "backend_evaluation_identity"
            )
            if backend_evaluation_identity is not None and not isinstance(
                backend_evaluation_identity, str
            ):
                raise TypeError("Backend evaluation identity must be a string")
            projected = ProjectedOptimizationSuccess.for_prepared(
                prepared,
                candidate_vector=tuple(
                    float.fromhex(str(value)) for value in vector_record
                ),
                backend_chi_square=backend_chi_square,
                execution_identity=execution_identity,
                backend_evaluation_identity=backend_evaluation_identity,
            )
            if projected_record != cls._projected_record(projected):
                raise ResamplingConstructionError(
                    "Stored projected candidate differs from exact prepared lineage"
                )
            evaluation_result = EvaluationResult.from_record(
                cast("Mapping[str, object]", evaluation_record),
                prepared.evaluation_plan,
            )
            success = cls.from_materialization(prepared, projected, evaluation_result)
        except ResamplingConstructionError:
            raise
        except Exception as error:
            raise ResamplingConstructionError(
                "Malformed validated replicate success record"
            ) from error
        if record != success.to_record():
            raise ResamplingConstructionError(
                "Stored success differs from its canonical fresh derivation"
            )
        return success

    def validate_integrity(
        self,
        plan: ResamplingPlan,
        request: ReplicateRequest,
    ) -> ValidatedReplicateSuccess:
        """Reconstruct this success from the exact operation-owned derivation."""
        try:
            canonical_draw = generate_resampling_draw(plan.dataset, request)
            canonical_prepared = ReplicateExecutionPlan.prepare(
                plan, request, canonical_draw
            )
            stored_prepared = replace(self.prepared)
            if (
                stored_prepared.identity != self.prepared.identity
                or stored_prepared.identity != canonical_prepared.identity
                or self.prepared.request.identity != request.identity
                or self.prepared.draw.identity != canonical_draw.identity
            ):
                raise ResamplingConstructionError(
                    "Successful replicate belongs to another prepared request or draw"
                )
            validation_engine = plan.source_engine.resampled(
                canonical_prepared.evaluation_plan,
                canonical_prepared.profile_bindings,
            )
            if (
                self.evaluation_result.evaluator_compatibility_identity
                != validation_engine.compatibility_identity
            ):
                raise ResamplingConstructionError(
                    "Successful evaluation has foreign evaluator compatibility"
                )
            canonical = type(self).from_materialization(
                canonical_prepared,
                self.projected,
                self.evaluation_result,
            )
            if self.to_record() != canonical.to_record():
                raise ResamplingConstructionError(
                    "Successful replicate differs from its canonical derivation"
                )
        except ResamplingConstructionError:
            raise
        except Exception as error:
            raise ResamplingConstructionError(
                "Successful replicate integrity validation failed"
            ) from error
        return canonical


@dataclass(frozen=True, slots=True)
class ReplicateOutcome:
    """One genuine terminal outcome in canonical replicate order."""

    plan_identity: str
    request_identity: str
    ordinal: int
    seed: int
    draw_identity: str | None
    stage: ReplicateStage
    disposition: ReplicateDisposition
    success: ValidatedReplicateSuccess | None = field(
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
        if (
            self.failure is not None
            and self.failure.transformed_data_identity != self.draw_identity
        ):
            raise ResamplingConstructionError(
                "Replicate outcome and projected path use different transformed data"
            )
        if self.success is not None and self.draw_identity is None:
            raise ResamplingConstructionError(
                "Successful replicate outcome requires a generated draw"
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
                    self.stage.value,
                    self.disposition.value,
                    payload.identity,
                ),
            ),
        )

    def validate_integrity(self, plan: ResamplingPlan) -> None:  # noqa: C901
        """Recursively validate this terminal payload against its exact request."""
        try:
            if self.ordinal < 0 or self.ordinal >= plan.replicate_count:
                raise ResamplingConstructionError(
                    "Replicate outcome ordinal is outside its exact plan"
                )
            request = plan.replicates[self.ordinal]
            if (
                self.plan_identity != plan.identity
                or self.request_identity != request.identity
                or self.seed != request.seed
            ):
                raise ResamplingConstructionError(
                    "Replicate outcome differs from its exact planned request"
                )
            succeeded = self.disposition is ReplicateDisposition.SUCCEEDED
            if succeeded != (self.success is not None) or succeeded == (
                self.failure is not None
            ):
                raise ResamplingConstructionError(
                    "Replicate outcome has inconsistent terminal payload"
                )
            if self.success is not None:
                canonical_success = self.success.validate_integrity(plan, request)
                if (
                    self.draw_identity != canonical_success.prepared.draw.identity
                    or self.stage is not ReplicateStage.TERMINAL
                ):
                    raise ResamplingConstructionError(
                        "Successful outcome differs from its generated draw or stage"
                    )
                payload_identity = canonical_success.identity
            elif self.failure is not None:
                canonical_failure = replace(self.failure)
                if (
                    canonical_failure.identity != self.failure.identity
                    or canonical_failure.transformed_data_identity != self.draw_identity
                ):
                    raise ResamplingConstructionError(
                        "Failed outcome has altered terminal evidence"
                    )
                payload_identity = canonical_failure.identity
            else:
                raise ResamplingConstructionError(
                    "Replicate outcome lacks terminal evidence"
                )
            expected_identity = _identity(
                "native-resampling-replicate-outcome",
                (
                    self.plan_identity,
                    self.request_identity,
                    self.ordinal,
                    self.seed,
                    self.draw_identity,
                    self.stage.value,
                    self.disposition.value,
                    payload_identity,
                ),
            )
            if self.identity != expected_identity:
                raise ResamplingConstructionError(
                    "Replicate outcome identity differs from its terminal payload"
                )
        except ResamplingConstructionError:
            raise
        except Exception as error:
            raise ResamplingConstructionError(
                "Replicate outcome integrity validation failed"
            ) from error


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

    def __post_init__(self) -> None:  # noqa: C901
        if not self.outcomes:
            raise ResamplingConstructionError(
                "Resampling evidence requires at least one genuine replicate outcome"
            )
        if self.population_identity != self.plan.dataset.identity:
            raise ResamplingConstructionError(
                "Resampling evidence belongs to another source population"
            )
        ordinals = tuple(outcome.ordinal for outcome in self.outcomes)
        if ordinals != tuple(sorted(ordinals)) or len(set(ordinals)) != len(ordinals):
            raise ResamplingConstructionError(
                "Resampling outcomes must use unique canonical ordinal order"
            )
        if any(
            outcome.ordinal >= self.plan.replicate_count
            or outcome.plan_identity != self.plan.identity
            or outcome.request_identity
            != self.plan.replicates[outcome.ordinal].identity
            or outcome.seed != self.plan.replicates[outcome.ordinal].seed
            for outcome in self.outcomes
        ):
            raise ResamplingConstructionError(
                "Resampling outcomes differ from their planned replicate identities"
            )
        for outcome in self.outcomes:
            outcome.validate_integrity(self.plan)
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
        completed = sum(
            outcome.disposition is not ReplicateDisposition.NOT_STARTED
            for outcome in self.outcomes
        )
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

    def validate_integrity(self) -> None:
        """Recompute this evidence and recursively reject altered descendants."""
        try:
            canonical = type(self)(
                self.plan,
                self.population_identity,
                self.outcomes,
                self.operation_terminal,
            )
            if (
                self.lifecycle != canonical.lifecycle
                or self.completed_count != canonical.completed_count
                or self.successful_count != canonical.successful_count
                or self.failed_count != canonical.failed_count
                or self.coverage_satisfied != canonical.coverage_satisfied
                or self.claims != canonical.claims
                or self.identity != canonical.identity
            ):
                raise ResamplingConstructionError(
                    "Resampling evidence differs from its recursive derivation"
                )
        except ResamplingConstructionError:
            raise
        except Exception as error:
            raise ResamplingConstructionError(
                "Resampling evidence integrity validation failed"
            ) from error


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


def _execute_native_candidate(
    prepared: ReplicateExecutionPlan,
    engine: EvaluationEngine,
) -> ProjectedOptimizationResult:
    """Construct and execute the closed native strategy with replicate-local state."""
    outcome = execute_direct_trf_candidate(
        prepared.problem,
        prepared.invocation,
        prepared.parameterization,
        engine,
    )
    if (
        outcome.terminal is DirectTrfCandidateTerminal.SUCCESS
        and outcome.candidate is not None
    ):
        return ProjectedOptimizationSuccess.for_prepared(
            prepared,
            candidate_vector=outcome.candidate.vector,
            backend_chi_square=outcome.candidate.chi_square,
            execution_identity=outcome.execution.identity,
            backend_evaluation_identity=outcome.candidate.evaluation_result.identity,
        )
    return ProjectedOptimizationFailure(
        prepared.draw.identity,
        ReplicateDisposition.FAILED,
        f"direct_trf_{outcome.terminal.value}",
        "Replicate-local Direct-TRF did not produce an eligible candidate",
        optimization_projection_identity=(
            prepared.request.optimization_projection_identity
        ),
        stage=ReplicateStage.EXECUTING,
        prepared_replicate_identity=prepared.identity,
        evaluation_plan_identity=prepared.evaluation_plan.identity,
        problem_identity=prepared.problem.identity,
        invocation_identity=prepared.invocation.identity,
        execution_identity=outcome.execution.identity,
    )


def _execute_replicate(  # noqa: C901 - ordered staged ownership boundary
    plan: ResamplingPlan,
    request: ReplicateRequest,
    candidate_test_hook: CandidateEvidenceTestHook | None,
) -> ReplicateOutcome:
    draw: ResamplingDraw | None = None
    prepared: ReplicateExecutionPlan | None = None

    def failure_outcome(
        disposition: ReplicateDisposition,
        stage: ReplicateStage,
        category: str,
        message: str,
        *,
        projected: ProjectedOptimizationSuccess | None = None,
    ) -> ReplicateOutcome:
        failure = ProjectedOptimizationFailure(
            None if draw is None else draw.identity,
            disposition,
            category,
            message,
            optimization_projection_identity=(request.optimization_projection_identity),
            stage=stage,
            prepared_replicate_identity=(
                None if prepared is None else prepared.identity
            ),
            evaluation_plan_identity=(
                None if projected is None else projected.evaluation_plan_identity
            ),
            problem_identity=(
                None if projected is None else projected.problem_identity
            ),
            invocation_identity=(
                None if projected is None else projected.invocation_identity
            ),
            execution_identity=(
                None if projected is None else projected.execution_identity
            ),
            component_outcome_identities=(
                () if projected is None else projected.component_outcome_identities
            ),
        )
        return ReplicateOutcome(
            plan.identity,
            request.identity,
            request.ordinal,
            request.seed,
            None if draw is None else draw.identity,
            stage,
            disposition,
            failure=failure,
        )

    try:
        draw = generate_resampling_draw(plan.dataset, request)
    except KeyboardInterrupt:
        return failure_outcome(
            ReplicateDisposition.INTERRUPTED,
            ReplicateStage.GENERATING,
            "generation_interrupted",
            "Replicate draw generation was interrupted",
        )
    except Exception as error:  # noqa: BLE001 - typed generation failure
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.GENERATING,
            "generation_failure",
            f"{type(error).__name__}: {error}",
        )

    try:
        prepared = ReplicateExecutionPlan.prepare(plan, request, draw)
    except KeyboardInterrupt:
        return failure_outcome(
            ReplicateDisposition.INTERRUPTED,
            ReplicateStage.DRAW_READY,
            "preparation_interrupted",
            "Replicate native preparation was interrupted after draw creation",
        )
    except Exception as error:  # noqa: BLE001 - typed preparation failure
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.DRAW_READY,
            "preparation_failure",
            f"{type(error).__name__}: {error}",
        )

    try:
        optimizer_engine = plan.source_engine.resampled(
            prepared.evaluation_plan, prepared.profile_bindings
        )
    except KeyboardInterrupt:
        return failure_outcome(
            ReplicateDisposition.INTERRUPTED,
            ReplicateStage.EXECUTOR_CREATING,
            "executor_creation_interrupted",
            "Replicate-local optimizer binding was interrupted",
        )
    except Exception as error:  # noqa: BLE001 - typed local binding failure
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.EXECUTOR_CREATING,
            "executor_creation_failure",
            f"{type(error).__name__}: {error}",
        )

    try:
        native_projected = _execute_native_candidate(prepared, optimizer_engine)
        projected = native_projected
        if candidate_test_hook is not None:
            projected = candidate_test_hook(prepared, projected)
    except KeyboardInterrupt:
        return failure_outcome(
            ReplicateDisposition.INTERRUPTED,
            ReplicateStage.EXECUTING,
            "optimization_interrupted",
            "Replicate projected optimization was interrupted",
        )
    except Exception as error:  # noqa: BLE001 - freeze typed replicate failure
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.EXECUTING,
            "projected_execution_failure",
            f"{type(error).__name__}: {error}",
        )

    if isinstance(projected, ProjectedOptimizationFailure):
        if projected.transformed_data_identity != draw.identity:
            return failure_outcome(
                ReplicateDisposition.FAILED,
                ReplicateStage.VALIDATING,
                "transformed_data_lineage_mismatch",
                "Projected failure used a different transformed dataset",
            )
        return ReplicateOutcome(
            plan.identity,
            request.identity,
            request.ordinal,
            request.seed,
            draw.identity,
            projected.stage,
            projected.disposition,
            failure=projected,
        )
    if not isinstance(projected, ProjectedOptimizationSuccess):
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "invalid_projected_outcome",
            "Projected optimization returned an unsupported outcome type",
        )
    if not isinstance(native_projected, ProjectedOptimizationSuccess):
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "ineligible_native_candidate",
            "The closed native strategy did not produce an eligible candidate",
        )

    expected = native_projected
    lineage_fields = (
        (
            "prepared_replicate",
            projected.prepared_replicate_identity,
            expected.prepared_replicate_identity,
        ),
        (
            "accepted_result",
            projected.accepted_result_identity,
            expected.accepted_result_identity,
        ),
        (
            "accepted_occurrence",
            projected.accepted_occurrence_identity,
            expected.accepted_occurrence_identity,
        ),
        (
            "source_snapshot_occurrence",
            projected.source_snapshot_occurrence_identity,
            expected.source_snapshot_occurrence_identity,
        ),
        (
            "source_snapshot_revision",
            projected.source_snapshot_revision,
            expected.source_snapshot_revision,
        ),
        (
            "transformed_data",
            projected.transformed_data_identity,
            expected.transformed_data_identity,
        ),
        (
            "optimization_projection",
            projected.optimization_projection_identity,
            expected.optimization_projection_identity,
        ),
        (
            "evaluation_plan",
            projected.evaluation_plan_identity,
            expected.evaluation_plan_identity,
        ),
        ("optimization_problem", projected.problem_identity, expected.problem_identity),
        (
            "parameterization",
            projected.parameterization_identity,
            expected.parameterization_identity,
        ),
        (
            "evaluator_parameterization",
            projected.evaluator_parameterization_identity,
            expected.evaluator_parameterization_identity,
        ),
        (
            "constraint_program",
            projected.constraint_program_identity,
            expected.constraint_program_identity,
        ),
        ("workflow", projected.workflow_identity, expected.workflow_identity),
        ("invocation", projected.invocation_identity, expected.invocation_identity),
        ("strategy", projected.strategy, expected.strategy),
        ("components", projected.component_identities, expected.component_identities),
        (
            "component_outcomes",
            projected.component_outcome_identities,
            expected.component_outcome_identities,
        ),
        ("execution", projected.execution_identity, expected.execution_identity),
        ("controlled_coordinates", projected.controlled_ids, expected.controlled_ids),
        ("candidate_vector", projected.candidate_vector, expected.candidate_vector),
        ("candidate", projected.candidate_identity, expected.candidate_identity),
        (
            "backend_evaluation",
            projected.backend_evaluation_identity,
            expected.backend_evaluation_identity,
        ),
    )
    mismatch = next(
        (name for name, actual, exact in lineage_fields if actual != exact),
        None,
    )
    if mismatch is not None:
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            f"{mismatch}_lineage_mismatch",
            f"Projected optimization returned foreign {mismatch} lineage",
            projected=projected,
        )
    try:
        validation_engine = plan.source_engine.resampled(
            prepared.evaluation_plan, prepared.profile_bindings
        )
        validation_evaluator = validation_engine.new_evaluator()
        lifecycle_frame = prepared.problem.lifecycle_frame(
            projected.candidate_vector, prepared.parameterization
        )
        evaluation_frame = EvaluationFrame.from_lifecycle_frame(
            prepared.parameterization, lifecycle_frame
        )
        evaluated = validation_evaluator.evaluate(evaluation_frame)
    except KeyboardInterrupt:
        return failure_outcome(
            ReplicateDisposition.INTERRUPTED,
            ReplicateStage.VALIDATING,
            "validation_interrupted",
            "Fresh candidate evaluation was interrupted",
            projected=projected,
        )
    except Exception as error:  # noqa: BLE001 - typed evidence materialization failure
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "fresh_evaluation_failure",
            f"{type(error).__name__}: {error}",
            projected=projected,
        )
    if isinstance(evaluated, EvaluationFailure):
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "fresh_evaluation_failure",
            f"{evaluated.category}: {evaluated.message}",
            projected=projected,
        )
    if not isinstance(evaluated, EvaluationResult):
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "invalid_fresh_evaluation",
            "Fresh evaluator returned an unsupported result type",
            projected=projected,
        )
    try:
        canonical_evaluation = EvaluationResult.from_record(
            evaluated.to_record(), prepared.evaluation_plan
        )
        resolved_items = tuple(
            (param_id, canonical_evaluation.resolved_values[param_id])
            for param_id in prepared.output_scope
        )
        ordinary_resolved = prepared.resolve_candidate(projected.candidate_vector)
        canonical_chi_square(canonical_evaluation.residuals)
        if (
            canonical_evaluation.evaluator_compatibility_identity
            != validation_engine.compatibility_identity
            or canonical_evaluation.plan_identity != prepared.evaluation_plan.identity
            or resolved_items != ordinary_resolved
            or tuple(canonical_evaluation.resolved_values) != prepared.output_scope
            or (
                expected.backend_evaluation_identity is not None
                and canonical_evaluation.identity
                != expected.backend_evaluation_identity
            )
        ):
            raise ResamplingConstructionError(
                "Fresh evaluation differs from exact complete replicate scope"
            )
    except Exception as error:  # noqa: BLE001 - canonical result qualification
        return failure_outcome(
            ReplicateDisposition.FAILED,
            ReplicateStage.VALIDATING,
            "fresh_materialization_failure",
            f"{type(error).__name__}: {error}",
            projected=projected,
        )
    success = ValidatedReplicateSuccess.from_materialization(
        prepared,
        projected,
        canonical_evaluation,
    )
    return ReplicateOutcome(
        plan.identity,
        request.identity,
        request.ordinal,
        request.seed,
        draw.identity,
        ReplicateStage.TERMINAL,
        ReplicateDisposition.SUCCEEDED,
        success=success,
    )


def _execute_serial_replicates(
    plan: ResamplingPlan,
    cancellation_probe: Callable[[], bool] | None,
    candidate_test_hook: CandidateEvidenceTestHook | None,
) -> tuple[list[ReplicateOutcome], OperationTerminal]:
    outcomes: list[ReplicateOutcome] = []
    terminal = OperationTerminal.COMPLETED
    for request in plan.replicates:
        if cancellation_probe is not None and cancellation_probe():
            return outcomes, OperationTerminal.CANCELLED
        outcome = _execute_replicate(plan, request, candidate_test_hook)
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
    execution: ExecutionSettings,
    cancellation_probe: Callable[[], bool] | None,
    candidate_test_hook: CandidateEvidenceTestHook | None,
) -> tuple[list[ReplicateOutcome], OperationTerminal]:
    worker_count = min(execution.workers, plan.replicate_count)
    outcomes: list[ReplicateOutcome] = []
    terminal = OperationTerminal.COMPLETED
    with native_thread_environment(execution.native_thread_env(parallel=True)):
        pool = ThreadPoolExecutor(max_workers=worker_count)
        next_ordinal = 0
        pending: set[Future[ReplicateOutcome]] = set()
        while next_ordinal < worker_count:
            pending.add(
                pool.submit(
                    _execute_replicate,
                    plan,
                    plan.replicates[next_ordinal],
                    candidate_test_hook,
                )
            )
            next_ordinal += 1
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
                while (
                    len(pending) < worker_count and next_ordinal < plan.replicate_count
                ):
                    pending.add(
                        pool.submit(
                            _execute_replicate,
                            plan,
                            plan.replicates[next_ordinal],
                            candidate_test_hook,
                        )
                    )
                    next_ordinal += 1
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
    *,
    execution: ExecutionSettings | None = None,
    cancellation_probe: Callable[[], bool] | None = None,
    candidate_test_hook: CandidateEvidenceTestHook | None = None,
) -> ResamplingOperation:
    """Execute canonical evidence-only replicates without analysis-state authority."""
    if (
        not accepted_occurrence_is_authoritative(accepted)
        or accepted.identity != plan.accepted_result_identity
        or accepted.occurrence_identity != plan.accepted_occurrence_identity
        or accepted.vector != plan.accepted_vector
    ):
        raise ResamplingConstructionError(
            "Resampling execution requires the plan's exact authoritative anchor"
        )
    _validate_accepted_source(
        accepted,
        plan.dataset,
        plan.source_problem,
        plan.parameterization,
    )
    settings = ExecutionSettings() if execution is None else execution
    if cancellation_probe is not None and cancellation_probe():
        outcomes: list[ReplicateOutcome] = []
        terminal = OperationTerminal.CANCELLED
    elif settings.workers == 1:
        outcomes, terminal = _execute_serial_replicates(
            plan,
            cancellation_probe,
            candidate_test_hook,
        )
    else:
        outcomes, terminal = _execute_parallel_replicates(
            plan,
            settings,
            cancellation_probe,
            candidate_test_hook,
        )
    completed_ordinals = {outcome.ordinal for outcome in outcomes}
    unstarted = tuple(
        request.ordinal
        for request in plan.replicates
        if request.ordinal not in completed_ordinals
    )
    for ordinal in unstarted:
        request = plan.replicates[ordinal]
        failure = ProjectedOptimizationFailure(
            None,
            ReplicateDisposition.NOT_STARTED,
            "not_started_after_terminal",
            "Replicate was not started after operation termination",
            optimization_projection_identity=(request.optimization_projection_identity),
            stage=ReplicateStage.PLANNED,
        )
        outcomes.append(
            ReplicateOutcome(
                plan.identity,
                request.identity,
                request.ordinal,
                request.seed,
                None,
                ReplicateStage.PLANNED,
                ReplicateDisposition.NOT_STARTED,
                failure=failure,
            )
        )
    outcomes.sort(key=lambda outcome: outcome.ordinal)
    evidence = ResamplingEvidence(
        plan,
        plan.dataset.identity,
        tuple(outcomes),
        terminal,
    )
    return ResamplingOperation(plan.identity, terminal, evidence, unstarted)


@dataclass(frozen=True, slots=True)
class SummaryFailure:
    """Typed reason why no resampling summary artifact was produced."""

    source_evidence_identity: str
    category: str
    message: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.source_evidence_identity or not self.category or not self.message:
            raise ResamplingConstructionError(
                "Summary failure requires category and message"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-summary-failure",
                (self.source_evidence_identity, self.category, self.message),
            ),
        )


@dataclass(frozen=True, slots=True)
class ResamplingSummaryPolicy:
    """Explicit complete-case percentile/covariance convention."""

    percentile_lower: float = 2.5
    percentile_median: float = 50.0
    percentile_upper: float = 97.5
    percentile_method: Literal["median_unbiased"] = "median_unbiased"
    covariance_delta_degrees_of_freedom: int = 1
    version: str = _SUMMARY_POLICY_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        values = tuple(
            _finite(value, name="summary percentile")
            for value in (
                self.percentile_lower,
                self.percentile_median,
                self.percentile_upper,
            )
        )
        if not 0.0 <= values[0] < values[1] < values[2] <= 100.0:
            raise ResamplingConstructionError(
                "Summary percentiles must be strictly ordered in [0, 100]"
            )
        if self.covariance_delta_degrees_of_freedom != 1:
            raise ResamplingConstructionError(
                "Native resampling covariance uses the declared sample denominator"
            )
        if self.percentile_method != "median_unbiased":
            raise ResamplingConstructionError(
                "Native resampling percentiles use the declared median-unbiased convention"
            )
        _non_empty_identity(self.version, name="summary policy version")
        object.__setattr__(self, "percentile_lower", values[0])
        object.__setattr__(self, "percentile_median", values[1])
        object.__setattr__(self, "percentile_upper", values[2])
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-resampling-summary-policy",
                (
                    tuple(_float_token(value) for value in values),
                    self.percentile_method,
                    self.covariance_delta_degrees_of_freedom,
                    self.version,
                    "common-complete-sample-set",
                ),
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
    """Canonical joint summary reconstructed from one exact evidence artifact."""

    evidence: ResamplingEvidence = field(repr=False, compare=False)
    policy: ResamplingSummaryPolicy = field(repr=False, compare=False)
    evidence_identity: str = field(init=False)
    output_scope: tuple[str, ...] = field(init=False)
    output_units: tuple[str, ...] = field(init=False)
    included_ordinals: tuple[int, ...] = field(init=False)
    unstarted_ordinals: tuple[int, ...] = field(init=False)
    exclusions: tuple[ResamplingExclusion, ...] = field(init=False)
    samples: tuple[SummarySample, ...] = field(init=False)
    distributions: tuple[ParameterDistribution, ...] = field(init=False)
    covariance: tuple[CovarianceEntry, ...] = field(init=False)
    correlations: tuple[CorrelationEntry, ...] = field(init=False)
    policy_identity: str = field(init=False)
    policy_version: str = field(init=False)
    claims: tuple[ClaimAssessment, ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        self.evidence.validate_integrity()
        successful = tuple(
            outcome for outcome in self.evidence.outcomes if outcome.success is not None
        )
        if len(successful) < self.evidence.plan.minimum_successful_count:
            raise ResamplingConstructionError(
                "Successful complete-scope coverage is below the declared minimum"
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
            != self.evidence.plan.output_scope
            for sample in samples
        ):
            raise ResamplingConstructionError(
                "Successful source outcomes do not share the complete scope"
            )
        matrix = np.asarray(
            [[value for _param_id, value in sample.items] for sample in samples],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(matrix)):
            raise ResamplingConstructionError(
                "A source sample contains a non-finite complete-scope value"
            )
        distributions, covariance, correlations = _build_summary_statistics(
            matrix,
            self.evidence.plan.output_scope,
            self.policy,
        )
        unstarted = tuple(
            outcome.ordinal
            for outcome in self.evidence.outcomes
            if outcome.disposition is ReplicateDisposition.NOT_STARTED
        )
        exclusions = tuple(
            ResamplingExclusion(
                outcome.ordinal,
                outcome.identity,
                outcome.disposition,
                outcome.failure.category if outcome.failure is not None else "unknown",
            )
            for outcome in self.evidence.outcomes
            if outcome.failure is not None
            and outcome.disposition is not ReplicateDisposition.NOT_STARTED
        )
        object.__setattr__(self, "evidence_identity", self.evidence.identity)
        object.__setattr__(self, "output_scope", self.evidence.plan.output_scope)
        object.__setattr__(self, "output_units", self.evidence.plan.output_units)
        object.__setattr__(
            self, "included_ordinals", tuple(item.ordinal for item in successful)
        )
        object.__setattr__(self, "unstarted_ordinals", unstarted)
        object.__setattr__(self, "exclusions", exclusions)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "distributions", distributions)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(self, "correlations", correlations)
        object.__setattr__(self, "policy_identity", self.policy.identity)
        object.__setattr__(self, "policy_version", self.policy.version)
        claims = (
            ClaimAssessment(
                "COMMON_SCOPE_INCLUSION",
                ClaimState.SATISFIED,
                self.policy.identity,
                (
                    f"scope={','.join(self.output_scope)}",
                    f"units={','.join(self.output_units)}",
                    f"included={len(self.included_ordinals)}",
                ),
            ),
            ClaimAssessment(
                "MINIMUM_SUCCESSFUL_COVERAGE",
                ClaimState.SATISFIED,
                self.policy.identity,
                (f"included={len(self.included_ordinals)}",),
            ),
            ClaimAssessment(
                "FINITE_EMPIRICAL_SUMMARY",
                ClaimState.SATISFIED,
                self.policy.identity,
            ),
        )
        object.__setattr__(self, "claims", claims)
        object.__setattr__(
            self,
            "identity",
            _identity("native-resampling-summary", self._identity_record(claims)),
        )

    def _identity_record(
        self, claims: Sequence[ClaimAssessment] | None = None
    ) -> object:
        exact_claims = self.claims if claims is None else claims
        return (
            self.evidence_identity,
            self.output_scope,
            self.output_units,
            self.included_ordinals,
            self.unstarted_ordinals,
            tuple(
                (item.ordinal, item.outcome_identity, _items_tokens(item.items))
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
                (item.parameter_a, item.parameter_b, _float_token(item.value))
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
            self.policy_identity,
            self.policy_version,
            tuple(
                (claim.name, claim.state.value, claim.policy_identity, claim.details)
                for claim in exact_claims
            ),
        )

    def to_record(self) -> dict[str, object]:
        """Serialize every derived field for exact reconstruction checks."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "evidence_identity": self.evidence_identity,
            "output_scope": list(self.output_scope),
            "output_units": list(self.output_units),
            "included_ordinals": list(self.included_ordinals),
            "unstarted_ordinals": list(self.unstarted_ordinals),
            "samples": [
                {
                    "ordinal": item.ordinal,
                    "outcome_identity": item.outcome_identity,
                    "items": [[key, _float_token(value)] for key, value in item.items],
                }
                for item in self.samples
            ],
            "exclusions": [
                {
                    "ordinal": item.ordinal,
                    "outcome_identity": item.outcome_identity,
                    "disposition": item.disposition.value,
                    "category": item.category,
                }
                for item in self.exclusions
            ],
            "distributions": [
                {
                    "parameter_id": item.parameter_id,
                    "mean": _float_token(item.mean),
                    "standard_deviation": _float_token(item.standard_deviation),
                    "median": _float_token(item.median),
                    "percentile_95_lower": _float_token(item.percentile_95_lower),
                    "percentile_95_upper": _float_token(item.percentile_95_upper),
                }
                for item in self.distributions
            ],
            "covariance": [
                {
                    "parameter_a": item.parameter_a,
                    "parameter_b": item.parameter_b,
                    "value": _float_token(item.value),
                }
                for item in self.covariance
            ],
            "correlations": [
                {
                    "parameter_a": item.parameter_a,
                    "parameter_b": item.parameter_b,
                    "availability": item.availability.value,
                    "value": None if item.value is None else _float_token(item.value),
                }
                for item in self.correlations
            ],
            "policy_identity": self.policy_identity,
            "policy_version": self.policy_version,
            "claims": [
                {
                    "name": claim.name,
                    "state": claim.state.value,
                    "policy_identity": claim.policy_identity,
                    "details": list(claim.details),
                }
                for claim in self.claims
            ],
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        source_evidence: ResamplingEvidence,
        policy: ResamplingSummaryPolicy,
    ) -> ResamplingSummary:
        """Recompute and reject any independently altered derived summary field."""
        expected = cls.from_evidence(source_evidence, policy)
        if record != expected.to_record():
            raise ResamplingConstructionError(
                "Stored resampling summary differs from canonical source evidence"
            )
        return expected

    @classmethod
    def from_evidence(
        cls,
        evidence: ResamplingEvidence,
        policy: ResamplingSummaryPolicy,
    ) -> ResamplingSummary:
        """Compute all summary fields from one common complete successful set."""
        return cls(evidence, policy)

    def claim(self, name: str) -> ClaimState:
        """Return one exact named claim; unknown claims fail closed."""
        return _claim_state(self.claims, name)

    def validate_integrity(self) -> None:
        """Recompute the canonical summary and reject any altered derived field."""
        type(self).from_record(self.to_record(), self.evidence, self.policy)


@dataclass(frozen=True, slots=True)
class ResamplingSummaryOutcome:
    """Closed result whose evidence identity comes only from its exact payload."""

    terminal: SummaryTerminal
    summary: ResamplingSummary | None = None
    failure: SummaryFailure | None = None
    evidence_identity: str = field(init=False)

    def __post_init__(self) -> None:
        completed = self.terminal is SummaryTerminal.COMPLETED
        if completed != (self.summary is not None) or completed == (
            self.failure is not None
        ):
            raise ResamplingConstructionError(
                "Summary outcome must contain exactly its typed terminal payload"
            )
        identity = (
            self.summary.evidence_identity
            if self.summary is not None
            else self.failure.source_evidence_identity
            if self.failure is not None
            else ""
        )
        object.__setattr__(self, "evidence_identity", identity)


def _build_summary_statistics(
    matrix: Array,
    output_scope: tuple[str, ...],
    policy: ResamplingSummaryPolicy,
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
                lower, median, upper = np.percentile(
                    values,
                    (
                        policy.percentile_lower,
                        policy.percentile_median,
                        policy.percentile_upper,
                    ),
                    method=policy.percentile_method,
                )
                distributions.append(
                    ParameterDistribution(
                        param_id,
                        float(np.mean(values)),
                        float(
                            np.std(
                                values,
                                ddof=policy.covariance_delta_degrees_of_freedom,
                            )
                        ),
                        float(median),
                        float(lower),
                        float(upper),
                    )
                )
            centered = matrix - np.mean(matrix, axis=0)
            covariance_matrix = (centered.T @ centered) / float(
                len(matrix) - policy.covariance_delta_degrees_of_freedom
            )
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
    policy: ResamplingSummaryPolicy | None = None,
) -> ResamplingSummaryOutcome:
    """Derive one summary from all and only complete successful scope rows."""
    exact_policy = ResamplingSummaryPolicy() if policy is None else policy
    try:
        evidence.validate_integrity()
    except ResamplingConstructionError as error:
        failure = SummaryFailure(
            evidence.identity,
            "invalid_source_evidence",
            f"Source evidence failed recursive integrity validation: {error}",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.SOURCE_INVALID,
            failure=failure,
        )
    if evidence.successful_count < evidence.plan.minimum_successful_count:
        failure = SummaryFailure(
            evidence.identity,
            "insufficient_successful_coverage",
            "Successful complete-scope replicate coverage is below the declared minimum",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.INSUFFICIENT_COVERAGE,
            failure=failure,
        )
    try:
        summary = ResamplingSummary.from_evidence(evidence, exact_policy)
    except (
        FloatingPointError,
        RuntimeWarning,
        OverflowError,
        ValueError,
    ) as error:
        failure = SummaryFailure(
            evidence.identity,
            "non_finite_summary_arithmetic",
            f"Summary arithmetic did not produce finite binary64 evidence: {error}",
        )
        return ResamplingSummaryOutcome(
            SummaryTerminal.SOURCE_INVALID,
            failure=failure,
        )
    return ResamplingSummaryOutcome(
        SummaryTerminal.COMPLETED,
        summary=summary,
    )
