"""Prospective native fixed-topology MCMC evidence qualification (#600).

This module is an isolated qualification seam. Production method parsing,
workflow orchestration, reporting, and parameter-state mutation do not consume
these artifacts yet.
"""

from __future__ import annotations

import hashlib
import json
import math
import secrets
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from threading import RLock
from typing import SupportsIndex, cast
from weakref import WeakKeyDictionary

import emcee
import numpy as np

from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
    canonical_chi_square,
)
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import ActiveParameterization
from chemex.typing import Array

_SCHEMA_VERSION = 1
_SEED_POLICY_VERSION = "sha256-structural-u64-v1"
_RNG_ALGORITHM = "numpy-pcg64-lhs+seedsequence-mt19937-emcee-v1"
_SAMPLER_VERSION = "chemex-emcee-ensemble-v1"
_INITIALIZATION_VERSION = "bounded-latin-hypercube-v1"
_PROPOSAL_VERSION = "emcee-stretch-v1"
_MAX_U64 = (1 << 64) - 1


class McmcConstructionError(ValueError):
    """Raised when a native MCMC policy or artifact is malformed."""


class McmcPolicyKind(StrEnum):
    """Closed prospective native MCMC policy modes."""

    CALIBRATED = "calibrated"
    CALIBRATION_CANDIDATE = "calibration_candidate"
    EXPERT = "expert"


class InitializationKind(StrEnum):
    """Closed native ensemble initialization constructions."""

    BOUNDED_LATIN_HYPERCUBE = "bounded_latin_hypercube"


class ProposalKind(StrEnum):
    """Closed native ensemble proposal kernels."""

    STRETCH = "stretch"


class McmcEvidenceLifecycle(StrEnum):
    """Lifecycle of a primary chain artifact containing genuine states."""

    COMPLETED = "completed"
    PARTIAL = "partial"


class McmcOperationTerminal(StrEnum):
    """Terminal disposition of one native MCMC execution operation."""

    COMPLETED = "completed"
    CANCELLED = "cancelled"
    FAILED = "failed"
    INTERRUPTED = "interrupted"


class McmcExecutionStage(StrEnum):
    """Last immutable checkpoint reached by one operation."""

    BEFORE_INITIALIZATION = "before_initialization"
    INITIALIZING = "initializing"
    AFTER_INITIALIZATION = "after_initialization"
    BEFORE_TRANSITION = "before_transition"
    DURING_TRANSITION = "during_transition"
    AFTER_COMPLETE_STATE = "after_complete_state"
    FRESH_VALIDATION = "fresh_validation"
    BEFORE_FINAL_ASSEMBLY = "before_final_assembly"


class McmcInitializationOutcome(StrEnum):
    """Truthful bounded-initialization disposition."""

    NOT_STARTED = "not_started"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


class McmcDiagnosticStatus(StrEnum):
    """Typed availability of an observed or estimated diagnostic."""

    AVAILABLE = "available"
    UNAVAILABLE = "unavailable"


class McmcDiagnosticReason(StrEnum):
    """Closed reason for unavailable MCMC diagnostic evidence."""

    NO_TRANSITIONS = "no_transitions"
    INSUFFICIENT_STATES = "insufficient_states"
    INSUFFICIENT_RETAINED_WINDOW = "insufficient_retained_window"
    ZERO_VARIANCE = "zero_variance"
    INVALID_AUTOCORRELATION = "invalid_autocorrelation"
    PARTIAL_CHAIN = "partial_chain"
    BACKEND_TRANSITION_INCONSISTENT = "backend_transition_inconsistent"
    UNVALIDATED_BACKEND_TRANSITIONS = "unvalidated_backend_transitions"
    UNRELIABLE_AUTOCORRELATION = "unreliable_autocorrelation"


class BackendTransitionEvidenceKind(StrEnum):
    """Authority of one backend transition-mask artifact."""

    OBSERVED_EXECUTION = "observed_execution"
    HISTORICAL_OBSERVATION = "historical_observation"
    HYPOTHETICAL = "hypothetical"


class McmcTrajectoryClaim(StrEnum):
    """Authority level of deterministic trajectory observations."""

    ORDINARY_CAPTURE = "ordinary_capture"


class McmcExecutionError(RuntimeError):
    """Raised internally when native evaluation cannot produce log density."""


class _McmcCancelled(RuntimeError):
    pass


class _McmcEvidenceWitness:
    """Opaque process-local witness for one canonical evidence occurrence."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> _McmcEvidenceWitness:
        raise TypeError(
            "MCMC evidence witnesses are minted only by canonical factories"
        )

    def __copy__(self) -> _McmcEvidenceWitness:
        raise TypeError("MCMC evidence witnesses cannot be copied")

    def __deepcopy__(self, _memo: object) -> _McmcEvidenceWitness:
        raise TypeError("MCMC evidence witnesses cannot be copied")

    def __reduce__(self) -> tuple[object, ...]:
        raise TypeError("MCMC evidence witnesses cannot be serialized")

    def __reduce_ex__(self, _protocol: SupportsIndex, /) -> tuple[object, ...]:
        raise TypeError("MCMC evidence witnesses cannot be serialized")


@dataclass(frozen=True, slots=True)
class _McmcEvidenceBinding:
    kind: str
    identity: str | None = None
    object_identity: int | None = None


_MCMC_EVIDENCE_WITNESSES: WeakKeyDictionary[
    _McmcEvidenceWitness,
    _McmcEvidenceBinding,
] = WeakKeyDictionary()
_MCMC_EVIDENCE_WITNESSES_LOCK = RLock()


def _mint_mcmc_evidence_witness(kind: str) -> _McmcEvidenceWitness:
    witness = object.__new__(_McmcEvidenceWitness)
    with _MCMC_EVIDENCE_WITNESSES_LOCK:
        _MCMC_EVIDENCE_WITNESSES[witness] = _McmcEvidenceBinding(kind)
    return witness


def _bind_mcmc_evidence_witness(
    witness: _McmcEvidenceWitness,
    artifact: object,
    kind: str,
    identity: str,
) -> bool:
    with _MCMC_EVIDENCE_WITNESSES_LOCK:
        binding = _MCMC_EVIDENCE_WITNESSES.get(witness)
        if binding == _McmcEvidenceBinding(kind):
            _MCMC_EVIDENCE_WITNESSES[witness] = _McmcEvidenceBinding(
                kind,
                identity,
                id(artifact),
            )
            return True
        return binding == _McmcEvidenceBinding(kind, identity, id(artifact))


def _mcmc_evidence_occurrence_is_authoritative(
    witness: _McmcEvidenceWitness | None,
    artifact: object,
    kind: str,
    identity: str,
) -> bool:
    if witness is None:
        return False
    with _MCMC_EVIDENCE_WITNESSES_LOCK:
        binding = _MCMC_EVIDENCE_WITNESSES.get(witness)
    return binding == _McmcEvidenceBinding(kind, identity, id(artifact))


class _BackendExecutionObservation:
    """Opaque process-local owner of masks emitted by one sampler occurrence."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> _BackendExecutionObservation:
        raise TypeError(
            "Backend observations are minted only by native sampler execution"
        )

    def __copy__(self) -> _BackendExecutionObservation:
        raise TypeError("Backend execution observations cannot be copied")

    def __deepcopy__(self, _memo: object) -> _BackendExecutionObservation:
        raise TypeError("Backend execution observations cannot be copied")

    def __reduce__(self) -> tuple[object, ...]:
        raise TypeError("Backend execution observations cannot be serialized")

    def __reduce_ex__(self, _protocol: SupportsIndex, /) -> tuple[object, ...]:
        raise TypeError("Backend execution observations cannot be serialized")


@dataclass(frozen=True, slots=True)
class _BackendExecutionObservationBinding:
    plan_identity: str
    policy_identity: str
    walkers: int
    dimension: int
    backend_implementation_identity: str
    execution_occurrence_identity: str
    transitions: tuple[tuple[int, str, str, tuple[bool, ...]], ...] = ()
    raw_capture_identity: str | None = None
    raw_capture_object_identity: int | None = None
    objective_request_count: int | None = None
    evaluation_request_count: int | None = None
    mask_payload_identity: str | None = None


_BACKEND_EXECUTION_OBSERVATIONS: WeakKeyDictionary[
    _BackendExecutionObservation,
    _BackendExecutionObservationBinding,
] = WeakKeyDictionary()
_BACKEND_EXECUTION_OBSERVATIONS_LOCK = RLock()


def _mint_backend_execution_observation(
    plan: McmcPlan,
    *,
    backend_implementation_identity: str,
) -> _BackendExecutionObservation:
    observation = object.__new__(_BackendExecutionObservation)
    binding = _BackendExecutionObservationBinding(
        plan.identity,
        plan.policy.identity,
        plan.policy.walkers,
        plan.policy.dimension,
        backend_implementation_identity,
        _identity(
            "native-mcmc-backend-execution-occurrence",
            (plan.identity, secrets.token_hex(32)),
        ),
    )
    with _BACKEND_EXECUTION_OBSERVATIONS_LOCK:
        _BACKEND_EXECUTION_OBSERVATIONS[observation] = binding
    return observation


def _record_backend_execution_transition(
    observation: _BackendExecutionObservation,
    previous: EnsembleState,
    following: EnsembleState,
    mask: Sequence[bool],
) -> None:
    exact_mask = tuple(mask)
    with _BACKEND_EXECUTION_OBSERVATIONS_LOCK:
        binding = _BACKEND_EXECUTION_OBSERVATIONS.get(observation)
        if (
            binding is None
            or following.ordinal != len(binding.transitions) + 1
            or previous.ordinal + 1 != following.ordinal
        ):
            raise McmcExecutionError(
                "Backend mask does not belong to the active sampler occurrence"
            )
        _BACKEND_EXECUTION_OBSERVATIONS[observation] = (
            _BackendExecutionObservationBinding(
                binding.plan_identity,
                binding.policy_identity,
                binding.walkers,
                binding.dimension,
                binding.backend_implementation_identity,
                binding.execution_occurrence_identity,
                (
                    *binding.transitions,
                    (
                        following.ordinal,
                        previous.identity,
                        following.identity,
                        exact_mask,
                    ),
                ),
            )
        )


def _bind_backend_observation_to_raw_capture(
    observation: _BackendExecutionObservation,
    raw_capture: RawMcmcCapture,
) -> _BackendExecutionObservationBinding:
    with _BACKEND_EXECUTION_OBSERVATIONS_LOCK:
        binding = _BACKEND_EXECUTION_OBSERVATIONS.get(observation)
        expected_edges = tuple(
            (following.ordinal, previous.identity, following.identity)
            for previous, following in zip(
                raw_capture.states,
                raw_capture.states[1:],
                strict=False,
            )
        )
        observed_edges = (
            tuple(item[:3] for item in binding.transitions) if binding else ()
        )
        mask_payload_identity = (
            _identity(
                "native-mcmc-backend-mask-payload",
                binding.transitions,
            )
            if binding
            else ""
        )
        if (
            binding is None
            or binding.plan_identity != raw_capture.plan_identity
            or binding.policy_identity != raw_capture.policy_identity
            or binding.walkers != raw_capture.walkers
            or binding.dimension != raw_capture.dimension
            or binding.execution_occurrence_identity
            != raw_capture.backend_execution_occurrence_identity
            or observed_edges != expected_edges
        ):
            raise McmcConstructionError(
                "Backend masks cannot bind to a foreign execution occurrence"
            )
        if binding.raw_capture_identity is None:
            bound = _BackendExecutionObservationBinding(
                binding.plan_identity,
                binding.policy_identity,
                binding.walkers,
                binding.dimension,
                binding.backend_implementation_identity,
                binding.execution_occurrence_identity,
                binding.transitions,
                raw_capture.identity,
                id(raw_capture),
                raw_capture.objective_request_count,
                raw_capture.evaluation_request_count,
                mask_payload_identity,
            )
            _BACKEND_EXECUTION_OBSERVATIONS[observation] = bound
            return bound
        if (
            binding.raw_capture_identity != raw_capture.identity
            or binding.raw_capture_object_identity != id(raw_capture)
            or binding.objective_request_count != raw_capture.objective_request_count
            or binding.evaluation_request_count != raw_capture.evaluation_request_count
            or binding.mask_payload_identity != mask_payload_identity
        ):
            raise McmcConstructionError(
                "Backend masks cannot bind to a foreign execution occurrence"
            )
        return binding


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _nonempty(value: str, *, name: str) -> str:
    if not value:
        raise McmcConstructionError(f"{name} must not be empty")
    return value


def _exact_keys(record: Mapping[str, object], expected: set[str], name: str) -> None:
    if set(record) != expected:
        raise McmcConstructionError(f"{name} record has unexpected fields")


def _positive_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise McmcConstructionError(f"{name} must be a positive integer")
    return value


def _nonnegative_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise McmcConstructionError(f"{name} must be a non-negative integer")
    return value


def _unsigned_seed(value: object, *, name: str) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < 0
        or value > _MAX_U64
    ):
        raise McmcConstructionError(f"{name} must be an unsigned 64-bit seed")
    return value


def _finite_positive(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise McmcConstructionError(f"{name} must be a positive finite scalar")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise McmcConstructionError(f"{name} must be a positive finite scalar")
    return result


def _derive_seed(root_seed: int, namespace: str, work_unit: str) -> int:
    fields = (
        b"chemex-native-seed",
        _SEED_POLICY_VERSION.encode("ascii"),
        root_seed.to_bytes(8, "big"),
        namespace.encode("ascii"),
        work_unit.encode("ascii"),
    )
    encoded = b"".join(len(item).to_bytes(8, "big") + item for item in fields)
    return int.from_bytes(hashlib.sha256(encoded).digest()[:8], "big")


def build_bounded_latin_hypercube(
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    *,
    walkers: int,
    seed: int,
) -> tuple[tuple[float, ...], ...]:
    """Construct one reproducible open-bound Latin-hypercube ensemble."""
    walker_count = _positive_integer(walkers, name="MCMC walker count")
    rng_seed = _unsigned_seed(seed, name="MCMC initialization seed")
    lower = np.asarray(tuple(lower_bounds), dtype=np.float64)
    upper = np.asarray(tuple(upper_bounds), dtype=np.float64)
    if (
        lower.ndim != 1
        or upper.ndim != 1
        or lower.shape != upper.shape
        or lower.size < 1
        or not np.all(np.isfinite(lower))
        or not np.all(np.isfinite(upper))
        or not np.all(lower < upper)
    ):
        raise McmcConstructionError(
            "MCMC initialization requires finite strictly ordered bounds"
        )
    rng = np.random.Generator(np.random.PCG64(rng_seed))
    dimension = int(lower.size)
    unit = np.empty((walker_count, dimension), dtype=np.float64)
    epsilon = np.finfo(np.float64).eps
    for coordinate in range(dimension):
        jitter = np.clip(rng.random(walker_count), epsilon, 1.0 - epsilon)
        strata = rng.permutation(walker_count)
        unit[:, coordinate] = (strata + jitter) / walker_count
    positions = lower + unit * (upper - lower)
    if not np.all((positions > lower) & (positions < upper)):
        raise McmcConstructionError(
            "Bounded Latin-hypercube construction reached a closed bound"
        )
    return tuple(tuple(float(value) for value in row) for row in positions)


@dataclass(frozen=True, slots=True)
class ResolvedMcmcPolicy:
    """Complete immutable MCMC topology and work allocation fixed pre-run."""

    kind: McmcPolicyKind
    policy_version: str
    dimension: int
    walkers: int
    burn_steps: int
    retained_steps: int
    root_seed: int
    provenance_identity: str
    qualification_dimension_range: tuple[int, int] | None = None
    proposal_scale: float = 2.0
    thin: int = 1
    initialization: InitializationKind = InitializationKind.BOUNDED_LATIN_HYPERCUBE
    proposal: ProposalKind = ProposalKind.STRETCH
    sampler_version: str = _SAMPLER_VERSION
    initialization_version: str = _INITIALIZATION_VERSION
    proposal_version: str = _PROPOSAL_VERSION
    seed_policy_version: str = _SEED_POLICY_VERSION
    rng_algorithm: str = _RNG_ALGORITHM
    total_steps: int = field(init=False)
    objective_request_budget: int = field(init=False)
    initialization_seed: int = field(init=False)
    sampler_seed: int = field(init=False)
    has_calibrated_adequacy: bool = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        dimension = _positive_integer(self.dimension, name="MCMC dimension")
        walkers = _positive_integer(self.walkers, name="MCMC walker count")
        burn = _nonnegative_integer(self.burn_steps, name="MCMC burn extent")
        retained = _positive_integer(
            self.retained_steps,
            name="MCMC retained extent",
        )
        root_seed = _unsigned_seed(self.root_seed, name="MCMC root seed")
        qualification_range = self.qualification_dimension_range
        if qualification_range is not None and (
            len(qualification_range) != 2
            or isinstance(qualification_range[0], bool)
            or isinstance(qualification_range[1], bool)
            or not all(isinstance(value, int) for value in qualification_range)
            or qualification_range[0] < 1
            or qualification_range[0] > qualification_range[1]
            or not qualification_range[0] <= dimension <= qualification_range[1]
        ):
            raise McmcConstructionError("MCMC qualification dimension range is invalid")
        if self.kind is McmcPolicyKind.CALIBRATED and qualification_range is None:
            raise McmcConstructionError(
                "Calibrated MCMC policy requires an explicit qualification stratum"
            )
        if self.kind is McmcPolicyKind.CALIBRATED:
            raise McmcConstructionError(
                "No repository-frozen calibration reference is available for "
                "native MCMC"
            )
        if self.kind is McmcPolicyKind.EXPERT and qualification_range is not None:
            raise McmcConstructionError(
                "Expert MCMC policy cannot inherit calibrated qualification"
            )
        proposal_scale = _finite_positive(
            self.proposal_scale,
            name="MCMC proposal scale",
        )
        if proposal_scale <= 1.0:
            raise McmcConstructionError("MCMC stretch proposal scale must exceed one")
        if walkers < max(4, 2 * dimension):
            raise McmcConstructionError(
                "MCMC walker count must be at least max(4, 2 * dimension)"
            )
        if self.thin != 1:
            raise McmcConstructionError(
                "Native MCMC retained evidence has fixed THIN=1"
            )
        for value, name in (
            (self.policy_version, "MCMC policy version"),
            (self.provenance_identity, "MCMC policy provenance identity"),
            (self.sampler_version, "MCMC sampler version"),
            (self.initialization_version, "MCMC initialization version"),
            (self.proposal_version, "MCMC proposal version"),
            (self.seed_policy_version, "MCMC seed policy version"),
            (self.rng_algorithm, "MCMC RNG algorithm"),
        ):
            _nonempty(value, name=name)
        total_steps = burn + retained
        objective_request_budget = walkers * (1 + total_steps)
        namespace = _identity(
            "native-mcmc-seed-namespace",
            (
                self.kind.value,
                self.policy_version,
                self.provenance_identity,
                dimension,
                walkers,
                burn,
                retained,
                qualification_range,
                self.initialization.value,
                self.proposal.value,
                float(proposal_scale).hex(),
            ),
        )
        initialization_seed = _derive_seed(root_seed, namespace, "initialization")
        sampler_seed = _derive_seed(root_seed, namespace, "sampler")
        if initialization_seed == sampler_seed:
            raise McmcConstructionError("MCMC derived seeds collided")
        identity = _identity(
            "native-mcmc-resolved-policy",
            (
                self.kind.value,
                self.policy_version,
                self.provenance_identity,
                dimension,
                walkers,
                burn,
                retained,
                qualification_range,
                self.thin,
                self.initialization.value,
                self.initialization_version,
                self.proposal.value,
                self.proposal_version,
                float(proposal_scale).hex(),
                self.sampler_version,
                objective_request_budget,
                root_seed,
                initialization_seed,
                sampler_seed,
                self.seed_policy_version,
                self.rng_algorithm,
            ),
        )
        object.__setattr__(self, "dimension", dimension)
        object.__setattr__(self, "walkers", walkers)
        object.__setattr__(self, "burn_steps", burn)
        object.__setattr__(self, "retained_steps", retained)
        object.__setattr__(self, "root_seed", root_seed)
        object.__setattr__(self, "proposal_scale", proposal_scale)
        object.__setattr__(self, "total_steps", total_steps)
        object.__setattr__(
            self,
            "objective_request_budget",
            objective_request_budget,
        )
        object.__setattr__(self, "initialization_seed", initialization_seed)
        object.__setattr__(self, "sampler_seed", sampler_seed)
        object.__setattr__(
            self,
            "has_calibrated_adequacy",
            self.kind is McmcPolicyKind.CALIBRATED,
        )
        object.__setattr__(self, "identity", identity)


@dataclass(frozen=True, slots=True)
class McmcCalibrationReference:
    """Future #604 content reference; it carries no authority by itself."""

    calibration_identity: str
    baseline_release_identity: str
    numerical_lane_requirement: str
    policy_version: str
    minimum_dimension: int
    maximum_dimension: int
    walkers: int
    burn_steps: int
    retained_steps: int
    proposal_scale: float = 2.0
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        minimum = _positive_integer(
            self.minimum_dimension,
            name="calibration-reference minimum dimension",
        )
        maximum = _positive_integer(
            self.maximum_dimension,
            name="calibration-reference maximum dimension",
        )
        if minimum > maximum:
            raise McmcConstructionError(
                "MCMC calibration reference has an inverted dimension stratum"
            )
        for value, name in (
            (self.calibration_identity, "MCMC calibration identity"),
            (self.baseline_release_identity, "MCMC baseline release identity"),
            (self.numerical_lane_requirement, "MCMC numerical lane requirement"),
            (self.policy_version, "MCMC calibrated policy version"),
        ):
            _nonempty(value, name=name)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-calibration-reference",
                (
                    self.calibration_identity,
                    self.baseline_release_identity,
                    self.numerical_lane_requirement,
                    self.policy_version,
                    minimum,
                    maximum,
                    self.walkers,
                    self.burn_steps,
                    self.retained_steps,
                    float(self.proposal_scale).hex(),
                ),
            ),
        )


def _repository_frozen_calibration_matches(
    reference: McmcCalibrationReference,
    authority: object,
) -> bool:
    """#604 will bind exact frozen references; none exist at #600."""
    _ = reference, authority
    return False


@dataclass(frozen=True, slots=True)
class CalibratedMcmcPolicy:
    """Future exact-reference path, deliberately unavailable before #604."""

    reference: McmcCalibrationReference
    authority: object = field(repr=False, compare=False)

    def resolve(self, *, dimension: int, root_seed: int) -> ResolvedMcmcPolicy:
        if not _repository_frozen_calibration_matches(
            self.reference,
            self.authority,
        ):
            raise McmcConstructionError(
                "MCMC calibration reference has no exact repository-frozen authority"
            )
        raise AssertionError("#604 must implement calibrated MCMC policy resolution")


@dataclass(frozen=True, slots=True)
class CalibrationCandidateMcmcPolicy:
    """Executable prospective shape without calibrated adequacy authority."""

    candidate_identity: str
    minimum_dimension: int
    maximum_dimension: int
    walkers: int
    burn_steps: int
    retained_steps: int
    proposal_scale: float = 2.0

    def resolve(self, *, dimension: int, root_seed: int) -> ResolvedMcmcPolicy:
        minimum = _positive_integer(
            self.minimum_dimension,
            name="candidate minimum dimension",
        )
        maximum = _positive_integer(
            self.maximum_dimension,
            name="candidate maximum dimension",
        )
        if minimum > maximum or not minimum <= dimension <= maximum:
            raise McmcConstructionError(
                "MCMC calibration candidate does not cover the requested dimension"
            )
        return ResolvedMcmcPolicy(
            kind=McmcPolicyKind.CALIBRATION_CANDIDATE,
            policy_version=_nonempty(
                self.candidate_identity,
                name="MCMC calibration candidate identity",
            ),
            dimension=dimension,
            walkers=self.walkers,
            burn_steps=self.burn_steps,
            retained_steps=self.retained_steps,
            root_seed=root_seed,
            provenance_identity=self.candidate_identity,
            qualification_dimension_range=(minimum, maximum),
            proposal_scale=self.proposal_scale,
        )


@dataclass(frozen=True, slots=True)
class ExpertMcmcPolicy:
    """Prospective topology overrides without inherited adequacy claims."""

    burn_steps: int
    retained_steps: int
    walkers: int
    expert_provenance: str
    policy_version: str = field(init=False, default="expert-v1")

    def resolve(self, *, dimension: int, root_seed: int) -> ResolvedMcmcPolicy:
        return ResolvedMcmcPolicy(
            kind=McmcPolicyKind.EXPERT,
            policy_version=_nonempty(
                self.policy_version,
                name="MCMC expert policy version",
            ),
            dimension=dimension,
            walkers=self.walkers,
            burn_steps=self.burn_steps,
            retained_steps=self.retained_steps,
            root_seed=root_seed,
            provenance_identity=_nonempty(
                self.expert_provenance,
                name="MCMC expert provenance",
            ),
        )


@dataclass(frozen=True, slots=True)
class McmcAcceptedAnchor:
    """Exact occurrence-owned scientific context for one MCMC request."""

    accepted_semantic_identity: str
    accepted_occurrence_identity: str
    accepted_occurrence_witness: object = field(repr=False, compare=False)
    source_occurrence_identity: str
    source_revision: int
    problem_identity: str
    parameterization_identity: str
    evaluator_parameterization_identity: str
    constraint_program_identity: str
    evaluation_plan_identity: str
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    held_items: tuple[tuple[str, float], ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    affine_feasibility_identity: str
    accepted_vector: tuple[float, ...]
    accepted_evaluation_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-accepted-anchor",
                (
                    self.accepted_semantic_identity,
                    self.accepted_occurrence_identity,
                    self.source_occurrence_identity,
                    self.source_revision,
                    self.problem_identity,
                    self.parameterization_identity,
                    self.evaluator_parameterization_identity,
                    self.constraint_program_identity,
                    self.evaluation_plan_identity,
                    self.coordinate_units,
                    tuple((key, float(value).hex()) for key, value in self.held_items),
                    tuple(float(value).hex() for value in self.lower_bounds),
                    tuple(float(value).hex() for value in self.upper_bounds),
                    self.affine_feasibility_identity,
                    tuple(float(value).hex() for value in self.accepted_vector),
                    self.accepted_evaluation_identity,
                ),
            ),
        )

    def validate(
        self,
        accepted: AcceptedFitResult,
        problem: OptimizationProblem,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
    ) -> None:
        """Reject any live object outside the exact frozen accepted lineage."""
        expected = (
            self.accepted_semantic_identity,
            self.accepted_occurrence_identity,
            self.source_occurrence_identity,
            self.source_revision,
            self.problem_identity,
            self.parameterization_identity,
            self.evaluator_parameterization_identity,
            self.constraint_program_identity,
            self.evaluation_plan_identity,
            self.coordinate_units,
            self.held_items,
            self.lower_bounds,
            self.upper_bounds,
            self.affine_feasibility_identity,
            self.accepted_vector,
            self.accepted_evaluation_identity,
        )
        actual = (
            accepted.identity,
            accepted.occurrence_identity,
            problem.source_snapshot.occurrence_identity,
            problem.source_snapshot.revision,
            problem.identity,
            parameterization.identity,
            parameterization.evaluator_identity,
            parameterization.program.fingerprint,
            engine.plan.identity,
            tuple(self.coordinate_units),
            problem.held_items,
            problem.lower_bounds,
            problem.upper_bounds,
            problem.affine_feasibility_identity,
            accepted.vector,
            accepted.evaluation_result.identity,
        )
        if (
            not accepted_occurrence_is_authoritative(accepted)
            or accepted.occurrence_witness is not self.accepted_occurrence_witness
            or actual != expected
            or accepted.problem_identity != problem.identity
            or accepted.evaluation_result.plan_identity != engine.plan.identity
            or accepted.evaluation_result.parameterization_identity
            != parameterization.evaluator_identity
            or problem.evaluation_plan_identity != engine.plan.identity
            or problem.parameterization_identity != parameterization.identity
            or problem.evaluator_parameterization_identity
            != parameterization.evaluator_identity
            or problem.constraint_program_identity
            != parameterization.program.fingerprint
        ):
            raise McmcConstructionError(
                "MCMC accepted occurrence or execution context differs from its "
                "exact anchor"
            )


@dataclass(frozen=True, slots=True)
class McmcPlan:
    """One accepted-result MCMC request with all choices frozen pre-run."""

    accepted: AcceptedFitResult = field(repr=False, compare=False)
    source_problem: OptimizationProblem = field(repr=False, compare=False)
    parameterization: ActiveParameterization = field(repr=False, compare=False)
    source_engine: EvaluationEngine = field(repr=False, compare=False)
    policy: ResolvedMcmcPolicy
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    anchor: McmcAcceptedAnchor = field(init=False)
    accepted_result_identity: str = field(init=False)
    accepted_occurrence_identity: str = field(init=False)
    accepted_occurrence_witness: object = field(
        init=False,
        repr=False,
        compare=False,
    )
    problem_identity: str = field(init=False)
    coordinate_ids: tuple[str, ...] = field(init=False)
    lower_bounds: tuple[float, ...] = field(init=False)
    upper_bounds: tuple[float, ...] = field(init=False)
    initial_ensemble: tuple[tuple[float, ...], ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        accepted = self.accepted
        problem = self.source_problem
        parameterization = self.parameterization
        engine = self.source_engine
        if not accepted_occurrence_is_authoritative(accepted):
            raise McmcConstructionError(
                "MCMC planning requires an exact authoritative accepted occurrence"
            )
        snapshot = problem.source_snapshot
        if (
            accepted.problem_identity != problem.identity
            or accepted.parameterization_identity != parameterization.identity
            or accepted.evaluator_parameterization_identity
            != parameterization.evaluator_identity
            or accepted.source_occurrence_identity != snapshot.occurrence_identity
            or accepted.source_revision != snapshot.revision
            or accepted.evaluation_result.plan_identity != engine.plan.identity
            or accepted.evaluation_result.parameterization_identity
            != parameterization.evaluator_identity
            or accepted.controlled_ids != problem.controlled_ids
            or accepted.commit_scope != problem.commit_scope
            or accepted.commit_items
            != accepted.evaluation_result.resolved_values.ordered_items()
        ):
            raise McmcConstructionError(
                "MCMC source does not match the exact accepted native lineage"
            )
        if (
            engine.plan.identity != problem.evaluation_plan_identity
            or engine.plan.parameterization_identity
            != problem.evaluator_parameterization_identity
            or engine.plan.constraint_program_identity
            != problem.constraint_program_identity
        ):
            raise McmcConstructionError(
                "MCMC evaluator does not match the accepted optimization problem"
            )
        coordinate_ids = tuple(problem.controlled_ids)
        if len(coordinate_ids) != self.policy.dimension:
            raise McmcConstructionError(
                "Resolved MCMC policy dimension differs from the accepted problem"
            )
        resolved_items = dict(accepted.commit_items)
        if tuple(resolved_items[param_id] for param_id in coordinate_ids) != (
            accepted.vector
        ):
            raise McmcConstructionError(
                "MCMC accepted vector differs from its materialized central values"
            )
        units = tuple(self.coordinate_units)
        if tuple(param_id for param_id, _unit in units) != coordinate_ids or any(
            not isinstance(unit, ParameterUnit) for _param_id, unit in units
        ):
            raise McmcConstructionError(
                "MCMC coordinate units must bind the ordered sampled scope"
            )
        if problem.affine_half_spaces or problem.affine_equalities:
            raise McmcConstructionError(
                "Native MCMC qualification currently supports box-bounded problems"
            )
        lower = tuple(problem.lower_bounds)
        upper = tuple(problem.upper_bounds)
        initial = build_bounded_latin_hypercube(
            lower,
            upper,
            walkers=self.policy.walkers,
            seed=self.policy.initialization_seed,
        )
        anchor = McmcAcceptedAnchor(
            accepted.identity,
            accepted.occurrence_identity,
            accepted.occurrence_witness,
            snapshot.occurrence_identity,
            snapshot.revision,
            problem.identity,
            parameterization.identity,
            parameterization.evaluator_identity,
            parameterization.program.fingerprint,
            engine.plan.identity,
            units,
            problem.held_items,
            lower,
            upper,
            problem.affine_feasibility_identity,
            accepted.vector,
            accepted.evaluation_result.identity,
        )
        identity = _identity(
            "native-mcmc-plan",
            (
                anchor.identity,
                problem.identity,
                parameterization.identity,
                engine.plan.identity,
                self.policy.identity,
                coordinate_ids,
                units,
                tuple(float(value).hex() for value in lower),
                tuple(float(value).hex() for value in upper),
                tuple(tuple(float(value).hex() for value in row) for row in initial),
            ),
        )
        object.__setattr__(self, "coordinate_units", units)
        object.__setattr__(self, "anchor", anchor)
        object.__setattr__(self, "accepted_result_identity", accepted.identity)
        object.__setattr__(
            self,
            "accepted_occurrence_identity",
            accepted.occurrence_identity,
        )
        object.__setattr__(
            self,
            "accepted_occurrence_witness",
            accepted.occurrence_witness,
        )
        object.__setattr__(self, "problem_identity", problem.identity)
        object.__setattr__(self, "coordinate_ids", coordinate_ids)
        object.__setattr__(self, "lower_bounds", lower)
        object.__setattr__(self, "upper_bounds", upper)
        object.__setattr__(self, "initial_ensemble", initial)
        object.__setattr__(self, "identity", identity)

    def validate_integrity(self, accepted: AcceptedFitResult | None = None) -> None:
        """Recursively recheck all live plan objects and its exact anchor."""
        source = self.accepted if accepted is None else accepted
        self.anchor.validate(
            source,
            self.source_problem,
            self.parameterization,
            self.source_engine,
        )
        expected = McmcPlan.for_accepted(
            source,
            source_problem=self.source_problem,
            parameterization=self.parameterization,
            source_engine=self.source_engine,
            policy=self.policy,
            coordinate_units=self.coordinate_units,
        )
        if (
            expected.identity != self.identity
            or expected.anchor.identity != self.anchor.identity
        ):
            raise McmcConstructionError(
                "MCMC plan failed recursive integrity validation"
            )

    @classmethod
    def for_accepted(
        cls,
        accepted: AcceptedFitResult,
        *,
        source_problem: OptimizationProblem,
        parameterization: ActiveParameterization,
        source_engine: EvaluationEngine,
        policy: ResolvedMcmcPolicy,
        coordinate_units: Sequence[tuple[str, ParameterUnit]],
    ) -> McmcPlan:
        """Bind a prospective MCMC run to one accepted native occurrence."""
        return cls(
            accepted,
            source_problem,
            parameterization,
            source_engine,
            policy,
            tuple(coordinate_units),
        )


def _finite_state_scalar(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise McmcConstructionError(f"{name} must be a finite binary64 scalar")
    result = float(value)
    if not math.isfinite(result):
        raise McmcConstructionError(f"{name} must be a finite binary64 scalar")
    return 0.0 if result == 0.0 else result


@dataclass(frozen=True, slots=True)
class EnsembleState:
    """One atomically completed state with scientific position/log-density facts."""

    ordinal: int
    positions: tuple[tuple[float, ...], ...]
    log_densities: tuple[float, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.ordinal < 0:
            raise McmcConstructionError("MCMC state ordinal must be non-negative")
        positions = tuple(
            tuple(
                _finite_state_scalar(value, name="MCMC state coordinate")
                for value in row
            )
            for row in self.positions
        )
        if (
            not positions
            or not positions[0]
            or any(len(row) != len(positions[0]) for row in positions)
        ):
            raise McmcConstructionError(
                "MCMC state positions must be a non-empty rectangular ensemble"
            )
        log_densities = tuple(
            _finite_state_scalar(value, name="MCMC state log density")
            for value in self.log_densities
        )
        if len(log_densities) != len(positions):
            raise McmcConstructionError(
                "MCMC state payload must cover every walker exactly once"
            )
        identity = _identity(
            "native-mcmc-ensemble-state",
            (
                self.ordinal,
                tuple(tuple(float(value).hex() for value in row) for row in positions),
                tuple(float(value).hex() for value in log_densities),
            ),
        )
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "log_densities", log_densities)
        object.__setattr__(self, "identity", identity)

    def to_record(self) -> dict[str, object]:
        """Serialize one complete state with exact binary64 tokens."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "ordinal": self.ordinal,
            "positions": [
                [float(value).hex() for value in row] for row in self.positions
            ],
            "log_densities": [float(value).hex() for value in self.log_densities],
            "identity": self.identity,
        }

    def validate_integrity(self) -> None:
        expected = EnsembleState(
            self.ordinal,
            self.positions,
            self.log_densities,
        )
        if expected.identity != self.identity:
            raise McmcConstructionError("MCMC ensemble-state content identity mismatch")

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> EnsembleState:
        """Restore and content-validate one complete ensemble state."""
        _exact_keys(
            record,
            {
                "schema_version",
                "ordinal",
                "positions",
                "log_densities",
                "identity",
            },
            "MCMC ensemble-state",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise McmcConstructionError("Unsupported MCMC ensemble-state schema")
        ordinal = record.get("ordinal")
        raw_positions = record.get("positions")
        raw_log_densities = record.get("log_densities")
        identity = record.get("identity")
        if (
            isinstance(ordinal, bool)
            or not isinstance(ordinal, int)
            or not isinstance(raw_positions, list)
            or not all(isinstance(row, list) for row in raw_positions)
            or not isinstance(raw_log_densities, list)
            or not isinstance(identity, str)
        ):
            raise McmcConstructionError("Malformed MCMC ensemble-state record")
        position_rows = cast("list[list[object]]", raw_positions)
        if any(
            not isinstance(value, str) for row in position_rows for value in row
        ) or any(not isinstance(value, str) for value in raw_log_densities):
            raise McmcConstructionError("Malformed MCMC ensemble-state binary64 token")
        try:
            positions = tuple(
                tuple(float.fromhex(cast("str", value)) for value in row)
                for row in position_rows
            )
            log_densities = tuple(
                float.fromhex(cast("str", value)) for value in raw_log_densities
            )
        except ValueError as error:
            raise McmcConstructionError(
                "Malformed MCMC ensemble-state binary64 token"
            ) from error
        state = cls(
            ordinal,
            positions,
            log_densities,
        )
        if state.identity != identity:
            raise McmcConstructionError("MCMC ensemble-state identity mismatch")
        return state


@dataclass(frozen=True, slots=True)
class RawTransitionEvent:
    """One non-authoritative backend accept/reject observation."""

    ordinal: int
    previous_state_identity: str
    next_state_identity: str
    accepted: tuple[bool, ...]
    proposal_count: int
    accepted_count: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.ordinal < 1 or self.proposal_count < 0 or self.accepted_count < 0:
            raise McmcConstructionError("Raw MCMC transition counts are invalid")
        accepted = tuple(self.accepted)
        if (
            any(not isinstance(value, bool) for value in accepted)
            or self.proposal_count != len(accepted)
            or self.accepted_count != sum(accepted)
        ):
            raise McmcConstructionError(
                "Raw MCMC transition acceptance accounting is inconsistent"
            )
        object.__setattr__(self, "accepted", accepted)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-raw-transition",
                (
                    self.ordinal,
                    self.previous_state_identity,
                    self.next_state_identity,
                    accepted,
                    self.proposal_count,
                    self.accepted_count,
                ),
            ),
        )

    def validate_integrity(self) -> None:
        expected = RawTransitionEvent(
            self.ordinal,
            self.previous_state_identity,
            self.next_state_identity,
            self.accepted,
            self.proposal_count,
            self.accepted_count,
        )
        if expected.identity != self.identity:
            raise McmcConstructionError("Raw MCMC transition content identity mismatch")


@dataclass(frozen=True, slots=True)
class BackendTransitionConsistencyFailure:
    """Structural inconsistency in a backend-declared rejected move."""

    transition_ordinal: int
    walker_ordinal: int
    category: str


def _backend_transition_consistency_failures(
    source_capture: RawMcmcCapture,
    transitions: Sequence[RawTransitionEvent],
) -> tuple[BackendTransitionConsistencyFailure, ...]:
    states = source_capture.states
    if len(transitions) != max(0, len(states) - 1):
        raise McmcConstructionError(
            "Backend transition evidence must cover every complete state edge"
        )
    failures: list[BackendTransitionConsistencyFailure] = []
    for transition, previous, following in zip(
        transitions,
        states[:-1],
        states[1:],
        strict=True,
    ):
        transition.validate_integrity()
        if (
            transition.ordinal != following.ordinal
            or transition.previous_state_identity != previous.identity
            or transition.next_state_identity != following.identity
            or transition.proposal_count != source_capture.walkers
        ):
            raise McmcConstructionError(
                "Backend transition labels differ from their exact state edge"
            )
        for walker, accepted in enumerate(transition.accepted):
            if not accepted and (
                following.positions[walker] != previous.positions[walker]
                or following.log_densities[walker] != previous.log_densities[walker]
            ):
                failures.append(
                    BackendTransitionConsistencyFailure(
                        transition.ordinal,
                        walker,
                        "rejected_move_changed_state",
                    )
                )
    return tuple(failures)


@dataclass(frozen=True, slots=True, init=False)
class BackendTransitionEvidence:
    """Diagnostic-only backend masks bound to exact scientific raw states."""

    source_capture: RawMcmcCapture = field(repr=False, compare=False)
    kind: BackendTransitionEvidenceKind
    observation_provenance: str
    execution_occurrence_identity: str | None
    transitions: tuple[RawTransitionEvent, ...]
    _live_observation_authority: _BackendExecutionObservation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    source_capture_identity: str = field(init=False)
    observed_mask_payload_identity: str = field(init=False)
    consistency_failures: tuple[BackendTransitionConsistencyFailure, ...] = field(
        init=False
    )
    identity: str = field(init=False)

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        raise TypeError(
            "Backend transition evidence is created only by canonical factories"
        )

    def __post_init__(self) -> None:
        if self._occurrence_witness is None or not self.observation_provenance:
            raise McmcConstructionError(
                "Backend transition evidence requires canonical diagnostic provenance"
            )
        observed = self.kind in {
            BackendTransitionEvidenceKind.OBSERVED_EXECUTION,
            BackendTransitionEvidenceKind.HISTORICAL_OBSERVATION,
        }
        if observed != (self.execution_occurrence_identity is not None):
            raise McmcConstructionError(
                "Backend transition authority and execution occurrence disagree"
            )
        if (self.kind is BackendTransitionEvidenceKind.OBSERVED_EXECUTION) != (
            self._live_observation_authority is not None
        ):
            raise McmcConstructionError(
                "Live backend evidence requires its exact observation authority"
            )
        self.source_capture.validate_integrity()
        transitions = tuple(self.transitions)
        failures = _backend_transition_consistency_failures(
            self.source_capture,
            transitions,
        )
        object.__setattr__(self, "transitions", transitions)
        object.__setattr__(
            self,
            "source_capture_identity",
            self.source_capture.identity,
        )
        object.__setattr__(
            self,
            "observed_mask_payload_identity",
            _identity(
                "native-mcmc-backend-mask-payload",
                tuple(
                    (
                        item.ordinal,
                        item.previous_state_identity,
                        item.next_state_identity,
                        item.accepted,
                    )
                    for item in transitions
                ),
            ),
        )
        object.__setattr__(self, "consistency_failures", failures)
        if self.kind is BackendTransitionEvidenceKind.OBSERVED_EXECUTION:
            authority = cast(
                "_BackendExecutionObservation",
                self._live_observation_authority,
            )
            binding = _bind_backend_observation_to_raw_capture(
                authority,
                self.source_capture,
            )
            if (
                binding.execution_occurrence_identity
                != self.execution_occurrence_identity
                or binding.backend_implementation_identity
                != self.observation_provenance
                or binding.mask_payload_identity != self.observed_mask_payload_identity
            ):
                raise McmcConstructionError(
                    "Observed backend evidence differs from its authority registry"
                )
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "backend-transition-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Backend transition evidence occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-backend-transition-evidence",
            (
                self.source_capture.identity,
                self.kind.value,
                self.observation_provenance,
                self.execution_occurrence_identity,
                self.observed_mask_payload_identity,
                tuple(item.identity for item in self.transitions),
                tuple(
                    (item.transition_ordinal, item.walker_ordinal, item.category)
                    for item in self.consistency_failures
                ),
            ),
        )

    @property
    def observed_transition_count(self) -> int:
        return len(self.transitions)

    @property
    def is_structurally_consistent(self) -> bool:
        return not self.consistency_failures

    @property
    def supports_observed_diagnostics(self) -> bool:
        return self.kind in {
            BackendTransitionEvidenceKind.OBSERVED_EXECUTION,
            BackendTransitionEvidenceKind.HISTORICAL_OBSERVATION,
        }

    def validate_integrity(self) -> None:
        if not isinstance(self.source_capture, RawMcmcCapture):
            raise McmcConstructionError(
                "Backend transition evidence has a foreign source object"
            )
        self.source_capture.validate_integrity()
        if self.source_capture_identity != self.source_capture.identity:
            raise McmcConstructionError(
                "Backend transition source identity differs from its source object"
            )
        if self.supports_observed_diagnostics and (
            self.execution_occurrence_identity
            != self.source_capture.backend_execution_occurrence_identity
        ):
            raise McmcConstructionError(
                "Backend evidence differs from its source execution occurrence"
            )
        if self.kind is BackendTransitionEvidenceKind.OBSERVED_EXECUTION:
            authority = cast(
                "_BackendExecutionObservation",
                self._live_observation_authority,
            )
            binding = _bind_backend_observation_to_raw_capture(
                authority,
                self.source_capture,
            )
            if (
                binding.execution_occurrence_identity
                != self.execution_occurrence_identity
                or binding.backend_implementation_identity
                != self.observation_provenance
                or binding.mask_payload_identity != self.observed_mask_payload_identity
            ):
                raise McmcConstructionError(
                    "Observed backend evidence differs from its authority registry"
                )
        expected_failures = _backend_transition_consistency_failures(
            self.source_capture,
            self.transitions,
        )
        if (
            self.consistency_failures != expected_failures
            or self._content_identity() != self.identity
        ):
            raise McmcConstructionError(
                "Backend transition diagnostic content identity mismatch"
            )
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "backend-transition-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Backend transition evidence is not its canonical occurrence"
            )

    @classmethod
    def from_capture(
        cls,
        source_capture: RawMcmcCapture,
        masks: Sequence[Sequence[bool]],
        *,
        observation_provenance: str,
    ) -> BackendTransitionEvidence:
        """Create explicitly hypothetical masks without execution authority."""
        exact_masks = tuple(tuple(mask) for mask in masks)
        states = source_capture.states
        if len(exact_masks) != max(0, len(states) - 1):
            raise McmcConstructionError(
                "Backend masks must cover every observed transition"
            )
        transitions = tuple(
            RawTransitionEvent(
                following.ordinal,
                previous.identity,
                following.identity,
                mask,
                source_capture.walkers,
                sum(mask),
            )
            for previous, following, mask in zip(
                states[:-1],
                states[1:],
                exact_masks,
                strict=True,
            )
        )
        return _initialize_backend_transition_evidence(
            cls,
            source_capture,
            BackendTransitionEvidenceKind.HYPOTHETICAL,
            observation_provenance,
            None,
            transitions,
            None,
            _mint_mcmc_evidence_witness("backend-transition-evidence"),
        )

    def with_hypothetical_masks(
        self,
        masks: Sequence[Sequence[bool]],
        *,
        observation_provenance: str,
    ) -> BackendTransitionEvidence:
        """Create an unvalidated alternative without changing scientific states."""
        self.validate_integrity()
        return self.from_capture(
            self.source_capture,
            masks,
            observation_provenance=observation_provenance,
        )

    def to_record(self) -> dict[str, object]:
        """Serialize backend evidence without its process-local authority witness."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "kind": self.kind.value,
            "source_capture_identity": self.source_capture_identity,
            "observation_provenance": self.observation_provenance,
            "execution_occurrence_identity": self.execution_occurrence_identity,
            "observed_mask_payload_identity": self.observed_mask_payload_identity,
            "masks": [list(item.accepted) for item in self.transitions],
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        *,
        source: BackendTransitionEvidence,
    ) -> BackendTransitionEvidence:
        """Restore historical evidence only against its exact live observed source."""
        source.validate_integrity()
        if record.get("execution_occurrence_identity") != (
            source.execution_occurrence_identity
        ):
            raise McmcConstructionError(
                "Stored backend evidence belongs to a foreign execution occurrence"
            )
        if record != source.to_record():
            raise McmcConstructionError(
                "Stored backend evidence differs from its exact observed source"
            )
        if source.kind is not BackendTransitionEvidenceKind.OBSERVED_EXECUTION:
            raise McmcConstructionError(
                "Historical backend evidence requires a live observed source"
            )
        return _initialize_backend_transition_evidence(
            cls,
            source.source_capture,
            BackendTransitionEvidenceKind.HISTORICAL_OBSERVATION,
            source.observation_provenance,
            source.execution_occurrence_identity,
            source.transitions,
            None,
            _mint_mcmc_evidence_witness("backend-transition-evidence"),
        )


def _mint_observed_backend_transition_evidence(
    observation: _BackendExecutionObservation,
    raw_capture: RawMcmcCapture,
) -> BackendTransitionEvidence:
    """Mint observed masks solely from the execution-owned occurrence registry."""
    binding = _bind_backend_observation_to_raw_capture(observation, raw_capture)
    transitions = tuple(
        RawTransitionEvent(
            ordinal,
            previous_identity,
            following_identity,
            mask,
            raw_capture.walkers,
            sum(mask),
        )
        for ordinal, previous_identity, following_identity, mask in binding.transitions
    )
    return _initialize_backend_transition_evidence(
        BackendTransitionEvidence,
        raw_capture,
        BackendTransitionEvidenceKind.OBSERVED_EXECUTION,
        binding.backend_implementation_identity,
        binding.execution_occurrence_identity,
        transitions,
        observation,
        _mint_mcmc_evidence_witness("backend-transition-evidence"),
    )


def _initialize_backend_transition_evidence(
    cls: type[BackendTransitionEvidence],
    source_capture: RawMcmcCapture,
    kind: BackendTransitionEvidenceKind,
    observation_provenance: str,
    execution_occurrence_identity: str | None,
    transitions: tuple[RawTransitionEvent, ...],
    live_observation_authority: _BackendExecutionObservation | None,
    occurrence_witness: _McmcEvidenceWitness,
) -> BackendTransitionEvidence:
    evidence = object.__new__(cls)
    object.__setattr__(evidence, "source_capture", source_capture)
    object.__setattr__(evidence, "kind", kind)
    object.__setattr__(evidence, "observation_provenance", observation_provenance)
    object.__setattr__(
        evidence,
        "execution_occurrence_identity",
        execution_occurrence_identity,
    )
    object.__setattr__(evidence, "transitions", transitions)
    object.__setattr__(
        evidence,
        "_live_observation_authority",
        live_observation_authority,
    )
    object.__setattr__(evidence, "_occurrence_witness", occurrence_witness)
    evidence.__post_init__()
    return evidence


@dataclass(frozen=True, slots=True)
class RawMcmcCapture:
    """Execution-only backend observations with no scientific authority."""

    plan_identity: str
    policy_identity: str
    root_seed: int
    walkers: int
    dimension: int
    total_steps: int
    backend_execution_occurrence_identity: str = field(repr=False, compare=False)
    terminal: McmcOperationTerminal
    states: tuple[EnsembleState, ...]
    objective_request_count: int
    evaluation_request_count: int
    stage: McmcExecutionStage
    initialization_outcome: McmcInitializationOutcome
    failure_category: str | None = None
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        states = tuple(self.states)
        if tuple(state.ordinal for state in states) != tuple(range(len(states))):
            raise McmcConstructionError(
                "Raw MCMC capture must contain a contiguous complete-state prefix"
            )
        if not self.plan_identity or not self.backend_execution_occurrence_identity:
            raise McmcConstructionError("Raw MCMC capture requires a plan identity")
        if self.objective_request_count < 0 or self.evaluation_request_count < 0:
            raise McmcConstructionError(
                "Raw MCMC capture request counts must be non-negative"
            )
        if self._occurrence_witness is None:
            raise McmcConstructionError(
                "Raw MCMC capture must come from sampler execution"
            )
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "raw-capture",
            self.identity,
        ):
            raise McmcConstructionError(
                "Raw MCMC capture occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-raw-capture",
            (
                self.plan_identity,
                self.policy_identity,
                self.root_seed,
                self.walkers,
                self.dimension,
                self.total_steps,
                self.terminal.value,
                tuple(state.identity for state in self.states),
                self.objective_request_count,
                self.evaluation_request_count,
                self.stage.value,
                self.initialization_outcome.value,
                self.failure_category,
            ),
        )

    @property
    def complete_state_count(self) -> int:
        return len(self.states)

    def validate_integrity(self) -> None:
        for state in self.states:
            state.validate_integrity()
        if self._content_identity() != self.identity:
            raise McmcConstructionError("Raw MCMC capture content identity mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "raw-capture",
            self.identity,
        ):
            raise McmcConstructionError(
                "Raw MCMC capture is not the exact execution-owned occurrence"
            )


def _mint_raw_capture(
    plan: McmcPlan,
    terminal: McmcOperationTerminal,
    states: Sequence[EnsembleState],
    objective_request_count: int,
    evaluation_request_count: int,
    stage: McmcExecutionStage,
    initialization_outcome: McmcInitializationOutcome,
    backend_observation: _BackendExecutionObservation,
    failure_category: str | None = None,
) -> RawMcmcCapture:
    exact_states = tuple(states)
    with _BACKEND_EXECUTION_OBSERVATIONS_LOCK:
        observation_binding = _BACKEND_EXECUTION_OBSERVATIONS.get(backend_observation)
    if (
        observation_binding is None
        or observation_binding.plan_identity != plan.identity
        or observation_binding.policy_identity != plan.policy.identity
        or observation_binding.walkers != plan.policy.walkers
        or observation_binding.dimension != plan.policy.dimension
    ):
        raise McmcConstructionError(
            "Raw MCMC capture requires its exact backend execution occurrence"
        )
    capture = RawMcmcCapture(
        plan.identity,
        plan.policy.identity,
        plan.policy.root_seed,
        plan.policy.walkers,
        plan.policy.dimension,
        plan.policy.total_steps,
        observation_binding.execution_occurrence_identity,
        terminal,
        exact_states,
        objective_request_count,
        evaluation_request_count,
        stage,
        initialization_outcome,
        failure_category,
        _mint_mcmc_evidence_witness("raw-capture"),
    )
    return capture


@dataclass(frozen=True, slots=True)
class McmcEvidence:
    """Primary fixed-topology chain evidence preserving ensemble states."""

    plan: McmcPlan = field(repr=False, compare=False)
    lifecycle: McmcEvidenceLifecycle
    terminal: McmcOperationTerminal
    states: tuple[EnsembleState, ...]
    objective_request_count: int
    evaluation_request_count: int
    failure_category: str | None = None
    raw_capture_identity: str = ""
    backend_execution_occurrence_identity: str = ""
    backend_transition_evidence: BackendTransitionEvidence | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    plan_identity: str = field(init=False)
    accepted_result_identity: str = field(init=False)
    coordinate_ids: tuple[str, ...] = field(init=False)
    coordinate_units: tuple[tuple[str, ParameterUnit], ...] = field(init=False)
    completed_transition_count: int = field(init=False)
    trajectory_claim: McmcTrajectoryClaim = field(init=False)
    canonical_lane_qualified: bool = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:  # noqa: C901 - complete evidence invariant
        if self._occurrence_witness is None:
            raise McmcConstructionError(
                "Primary MCMC evidence must come from fresh native validation"
            )
        if not self.backend_execution_occurrence_identity:
            raise McmcConstructionError(
                "Primary MCMC evidence requires its backend execution occurrence"
            )
        states = tuple(self.states)
        expected_ordinals = tuple(range(len(states)))
        if not states or tuple(state.ordinal for state in states) != expected_ordinals:
            raise McmcConstructionError(
                "MCMC evidence must preserve a contiguous prefix of complete states"
            )
        if states[0].positions != self.plan.initial_ensemble:
            raise McmcConstructionError(
                "MCMC evidence initial state differs from its frozen plan"
            )
        for state in states:
            if len(state.positions) != self.plan.policy.walkers or any(
                len(position) != self.plan.policy.dimension
                for position in state.positions
            ):
                raise McmcConstructionError(
                    "MCMC evidence state topology differs from its frozen plan"
                )
        completed = len(states) - 1
        if completed > self.plan.policy.total_steps:
            raise McmcConstructionError("MCMC evidence exceeds its frozen topology")
        if self.lifecycle is McmcEvidenceLifecycle.COMPLETED:
            if (
                self.terminal is not McmcOperationTerminal.COMPLETED
                or completed != self.plan.policy.total_steps
                or self.failure_category is not None
            ):
                raise McmcConstructionError(
                    "Completed MCMC evidence must contain the complete frozen chain"
                )
        elif self.terminal is McmcOperationTerminal.COMPLETED:
            raise McmcConstructionError(
                "Partial MCMC evidence cannot have a completed terminal"
            )
        if (
            self.objective_request_count < len(states[0].positions)
            or self.objective_request_count > self.plan.policy.objective_request_budget
            or self.evaluation_request_count < 0
            or self.evaluation_request_count > self.objective_request_count
        ):
            raise McmcConstructionError("MCMC request accounting is inconsistent")
        if self.backend_transition_evidence is not None:
            self._validate_backend_occurrence()
        identity = _identity(
            "native-mcmc-primary-evidence",
            (
                self.plan.identity,
                self.plan.accepted_result_identity,
                self.lifecycle.value,
                self.terminal.value,
                tuple(state.identity for state in states),
                self.objective_request_count,
                self.evaluation_request_count,
                self.failure_category,
                self.raw_capture_identity,
            ),
        )
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "plan_identity", self.plan.identity)
        object.__setattr__(
            self,
            "accepted_result_identity",
            self.plan.accepted_result_identity,
        )
        object.__setattr__(self, "coordinate_ids", self.plan.coordinate_ids)
        object.__setattr__(self, "coordinate_units", self.plan.coordinate_units)
        object.__setattr__(self, "completed_transition_count", completed)
        object.__setattr__(
            self,
            "trajectory_claim",
            McmcTrajectoryClaim.ORDINARY_CAPTURE,
        )
        object.__setattr__(self, "canonical_lane_qualified", False)
        object.__setattr__(self, "identity", identity)
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "primary-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Primary MCMC evidence occurrence witness is already bound"
            )

    def validate_integrity(self) -> None:
        """Recursively reject reconstructed or altered primary evidence."""
        self.plan.validate_integrity()
        for state in self.states:
            state.validate_integrity()
        if self.backend_transition_evidence is not None:
            self._validate_backend_occurrence()
        expected_identity = _identity(
            "native-mcmc-primary-evidence",
            (
                self.plan.identity,
                self.plan.accepted_result_identity,
                self.lifecycle.value,
                self.terminal.value,
                tuple(state.identity for state in self.states),
                self.objective_request_count,
                self.evaluation_request_count,
                self.failure_category,
                self.raw_capture_identity,
            ),
        )
        if expected_identity != self.identity:
            raise McmcConstructionError(
                "Primary MCMC evidence content identity mismatch"
            )
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "primary-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Primary MCMC evidence is not its validation-owned occurrence"
            )

    def _validate_backend_occurrence(self) -> None:
        backend = cast(
            "BackendTransitionEvidence",
            self.backend_transition_evidence,
        )
        backend.validate_integrity()
        source_capture = backend.source_capture
        if (
            not backend.supports_observed_diagnostics
            or backend.execution_occurrence_identity
            != self.backend_execution_occurrence_identity
            or source_capture.backend_execution_occurrence_identity
            != self.backend_execution_occurrence_identity
            or source_capture.identity != self.raw_capture_identity
            or source_capture.plan_identity != self.plan.identity
            or source_capture.policy_identity != self.plan.policy.identity
            or source_capture.walkers != self.plan.policy.walkers
            or source_capture.dimension != self.plan.policy.dimension
        ):
            raise McmcConstructionError(
                "Primary backend diagnostics differ from its exact execution occurrence"
            )

    def to_record(self) -> dict[str, object]:
        """Serialize primary evidence without live evaluator or sampler objects."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "plan_identity": self.plan_identity,
            "accepted_result_identity": self.accepted_result_identity,
            "lifecycle": self.lifecycle.value,
            "terminal": self.terminal.value,
            "states": [state.to_record() for state in self.states],
            "objective_request_count": self.objective_request_count,
            "evaluation_request_count": self.evaluation_request_count,
            "failure_category": self.failure_category,
            "raw_capture_identity": self.raw_capture_identity,
            "backend_execution_occurrence_identity": (
                self.backend_execution_occurrence_identity
            ),
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        *,
        plan: McmcPlan,
        raw_capture: RawMcmcCapture | None = None,
    ) -> McmcEvidence:
        """Restore primary evidence against its exact trusted execution plan."""
        _exact_keys(
            record,
            {
                "schema_version",
                "plan_identity",
                "accepted_result_identity",
                "lifecycle",
                "terminal",
                "states",
                "objective_request_count",
                "evaluation_request_count",
                "failure_category",
                "raw_capture_identity",
                "backend_execution_occurrence_identity",
                "identity",
            },
            "MCMC evidence",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise McmcConstructionError("Unsupported MCMC evidence schema")
        if (
            record.get("plan_identity") != plan.identity
            or record.get("accepted_result_identity") != plan.accepted_result_identity
        ):
            raise McmcConstructionError("MCMC evidence lineage identity mismatch")
        raw_states = record.get("states")
        objective_count = record.get("objective_request_count")
        evaluation_count = record.get("evaluation_request_count")
        failure_category = record.get("failure_category")
        identity = record.get("identity")
        if (
            not isinstance(raw_states, list)
            or not all(isinstance(item, Mapping) for item in raw_states)
            or isinstance(objective_count, bool)
            or not isinstance(objective_count, int)
            or isinstance(evaluation_count, bool)
            or not isinstance(evaluation_count, int)
            or (failure_category is not None and not isinstance(failure_category, str))
            or not isinstance(identity, str)
        ):
            raise McmcConstructionError("Malformed MCMC evidence record")
        try:
            lifecycle = McmcEvidenceLifecycle(record.get("lifecycle"))
            terminal = McmcOperationTerminal(record.get("terminal"))
        except (TypeError, ValueError) as error:
            raise McmcConstructionError(
                "Malformed MCMC evidence lifecycle or terminal"
            ) from error
        _ = lifecycle, terminal, raw_states, objective_count, evaluation_count
        _ = failure_category, identity
        if raw_capture is None:
            raise McmcConstructionError(
                "Primary MCMC evidence reconstruction requires its raw capture and "
                "fresh validation"
            )
        validation = validate_raw_mcmc_capture(plan, raw_capture)
        expected = validation.primary_evidence
        if expected is None or record != expected.to_record():
            raise McmcConstructionError(
                "Stored MCMC evidence identity differs from fresh canonical validation"
            )
        return expected


def _mint_primary_evidence(
    plan: McmcPlan,
    raw_capture: RawMcmcCapture,
    lifecycle: McmcEvidenceLifecycle,
    states: Sequence[EnsembleState],
    backend_transition_evidence: BackendTransitionEvidence | None,
) -> McmcEvidence:
    evidence = McmcEvidence(
        plan,
        lifecycle,
        raw_capture.terminal,
        tuple(states),
        raw_capture.objective_request_count,
        raw_capture.evaluation_request_count,
        raw_capture.failure_category,
        raw_capture.identity,
        raw_capture.backend_execution_occurrence_identity,
        backend_transition_evidence,
        _mint_mcmc_evidence_witness("primary-evidence"),
    )
    return evidence


@dataclass(frozen=True, slots=True)
class McmcValidationFailure:
    """One typed fresh-validation failure for a state/walker observation."""

    state_ordinal: int
    walker_ordinal: int
    category: str
    message: str
    backend_log_density: float | None = None
    authoritative_log_density: float | None = None


@dataclass(frozen=True, slots=True)
class McmcChainValidation:
    """Fresh-workspace validation and any qualified primary-chain result."""

    raw_capture_identity: str
    validated_states: tuple[EnsembleState, ...]
    failures: tuple[McmcValidationFailure, ...]
    expected_state_count: int
    primary_evidence: McmcEvidence | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "identity", self._content_identity())

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-chain-validation",
            (
                self.raw_capture_identity,
                tuple(state.identity for state in self.validated_states),
                tuple(
                    (
                        item.state_ordinal,
                        item.walker_ordinal,
                        item.category,
                        item.message,
                        None
                        if item.backend_log_density is None
                        else float(item.backend_log_density).hex(),
                        None
                        if item.authoritative_log_density is None
                        else float(item.authoritative_log_density).hex(),
                    )
                    for item in self.failures
                ),
                self.expected_state_count,
                None
                if self.primary_evidence is None
                else self.primary_evidence.identity,
            ),
        )

    def validate_integrity(self) -> None:
        for state in self.validated_states:
            state.validate_integrity()
        if self.primary_evidence is not None:
            self.primary_evidence.validate_integrity()
        if self._content_identity() != self.identity:
            raise McmcConstructionError("MCMC chain validation content mismatch")

    @property
    def complete_state_count(self) -> int:
        return len(self.validated_states)

    @property
    def is_complete(self) -> bool:
        return (
            not self.failures
            and self.complete_state_count == self.expected_state_count
            and self.primary_evidence is not None
        )


def _fresh_authoritative_log_density(
    plan: McmcPlan,
    evaluator: BoundEvaluator,
    position: tuple[float, ...],
) -> float:
    """Re-derive one stored state through a validation-only evaluator."""
    candidate = np.asarray(position, dtype=np.float64)
    if (
        candidate.shape != (plan.policy.dimension,)
        or not np.all(np.isfinite(candidate))
        or np.any(candidate <= np.asarray(plan.lower_bounds))
        or np.any(candidate >= np.asarray(plan.upper_bounds))
    ):
        raise McmcExecutionError("stored MCMC state is outside its frozen bounds")
    lifecycle = plan.source_problem.lifecycle_frame(
        candidate,
        plan.parameterization,
    )
    frame = EvaluationFrame.from_lifecycle_frame(plan.parameterization, lifecycle)
    outcome = evaluator.evaluate(frame)
    if isinstance(outcome, EvaluationFailure):
        raise McmcExecutionError(
            f"fresh MCMC validation failed: {outcome.category}: {outcome.message}"
        )
    return -0.5 * canonical_chi_square(outcome.residuals)


def validate_raw_mcmc_capture(  # noqa: C901 - state/walker validation boundary
    plan: McmcPlan,
    raw_capture: RawMcmcCapture,
    *,
    backend_transition_evidence: BackendTransitionEvidence | None = None,
) -> McmcChainValidation:
    """Freshly validate raw backend observations before minting evidence."""
    expected_state_count = plan.policy.total_steps + 1
    try:
        plan.validate_integrity()
        raw_capture.validate_integrity()
    except McmcConstructionError as error:
        return McmcChainValidation(
            raw_capture.identity,
            (),
            (McmcValidationFailure(0, 0, "source_integrity_failure", str(error)),),
            expected_state_count,
            None,
        )
    if raw_capture.plan_identity != plan.identity:
        return McmcChainValidation(
            raw_capture.identity,
            (),
            (
                McmcValidationFailure(
                    0,
                    0,
                    "plan_identity_mismatch",
                    "raw MCMC capture is not bound to the supplied frozen plan",
                ),
            ),
            expected_state_count,
            None,
        )
    evaluator = plan.source_engine.new_evaluator()
    validated: list[EnsembleState] = []
    failures: list[McmcValidationFailure] = []
    for state in raw_capture.states:
        if (
            len(state.positions) != plan.policy.walkers
            or len(state.log_densities) != plan.policy.walkers
            or any(
                len(position) != plan.policy.dimension for position in state.positions
            )
        ):
            failures.append(
                McmcValidationFailure(
                    state.ordinal,
                    0,
                    "state_topology_mismatch",
                    "raw MCMC state differs from the frozen walker topology",
                )
            )
            break
        if state.ordinal == 0 and state.positions != plan.initial_ensemble:
            failures.append(
                McmcValidationFailure(
                    0,
                    0,
                    "initial_ensemble_mismatch",
                    "raw MCMC initial state differs from its frozen construction",
                )
            )
            break
        authoritative: list[float] = []
        for walker, (position, backend) in enumerate(
            zip(state.positions, state.log_densities, strict=True)
        ):
            try:
                value = _fresh_authoritative_log_density(plan, evaluator, position)
            except Exception as error:  # noqa: BLE001 - typed validation evidence
                failures.append(
                    McmcValidationFailure(
                        state.ordinal,
                        walker,
                        "fresh_native_evaluation_failed",
                        f"{type(error).__name__}: {error}",
                        backend,
                    )
                )
                break
            authoritative.append(value)
            if backend != value:
                failures.append(
                    McmcValidationFailure(
                        state.ordinal,
                        walker,
                        "backend_log_density_mismatch",
                        "backend log density differs from fresh native evaluation",
                        backend,
                        value,
                    )
                )
                break
        if failures:
            break
        validated.append(
            EnsembleState(
                state.ordinal,
                state.positions,
                tuple(authoritative),
            )
        )
    primary: McmcEvidence | None = None
    complete = (
        raw_capture.terminal is McmcOperationTerminal.COMPLETED
        and not failures
        and len(validated) == expected_state_count
    )
    if complete:
        primary = _mint_primary_evidence(
            plan,
            raw_capture,
            McmcEvidenceLifecycle.COMPLETED,
            tuple(validated),
            backend_transition_evidence,
        )
    elif not failures and validated:
        primary = _mint_primary_evidence(
            plan,
            raw_capture,
            McmcEvidenceLifecycle.PARTIAL,
            tuple(validated),
            backend_transition_evidence,
        )
    elif not failures and raw_capture.terminal is McmcOperationTerminal.COMPLETED:
        failures.append(
            McmcValidationFailure(
                len(validated),
                0,
                "incomplete_completed_capture",
                "completed raw MCMC capture does not contain its frozen topology",
            )
        )
    return McmcChainValidation(
        raw_capture.identity,
        tuple(validated),
        tuple(failures),
        expected_state_count,
        primary,
    )


@dataclass(frozen=True, slots=True)
class McmcOperation:
    """Frozen terminal record for one attempted MCMC evidence execution."""

    terminal: McmcOperationTerminal
    evidence: McmcEvidence | None
    failure_category: str | None = None
    failure_message: str = ""
    raw_capture: RawMcmcCapture | None = None
    validation: McmcChainValidation | None = None
    backend_transition_evidence: BackendTransitionEvidence | None = None
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.terminal is McmcOperationTerminal.COMPLETED:
            if (
                self.evidence is None
                or self.raw_capture is None
                or self.validation is None
                or not self.validation.is_complete
                or self.validation.primary_evidence is not self.evidence
                or self.failure_category is not None
            ):
                raise McmcConstructionError(
                    "Completed MCMC operation requires completed evidence"
                )
        elif self.failure_category is None:
            raise McmcConstructionError(
                "Non-completed MCMC operation requires typed failure evidence"
            )
        if (
            self.raw_capture is None
            or self.validation is None
            or self.backend_transition_evidence is None
        ):
            raise McmcConstructionError(
                "Every MCMC operation requires raw and backend diagnostic evidence"
            )
        if (
            self.backend_transition_evidence.kind
            is not BackendTransitionEvidenceKind.OBSERVED_EXECUTION
            or self.backend_transition_evidence.source_capture is not self.raw_capture
        ):
            raise McmcConstructionError(
                "MCMC operation requires execution-owned backend observations"
            )
        object.__setattr__(self, "identity", self._content_identity())
        if self._occurrence_witness is None or not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "operation",
            self.identity,
        ):
            raise McmcConstructionError(
                "MCMC operation must come from one canonical execution occurrence"
            )

    def _content_identity(self) -> str:
        capture = cast("RawMcmcCapture", self.raw_capture)
        validation = cast("McmcChainValidation", self.validation)
        backend = cast(
            "BackendTransitionEvidence",
            self.backend_transition_evidence,
        )
        return _identity(
            "native-mcmc-operation",
            (
                self.terminal.value,
                None if self.evidence is None else self.evidence.identity,
                self.failure_category,
                self.failure_message,
                capture.identity,
                validation.identity,
                backend.identity,
            ),
        )

    @property
    def complete_state_count(self) -> int:
        capture = cast("RawMcmcCapture", self.raw_capture)
        return capture.complete_state_count

    @property
    def stage(self) -> McmcExecutionStage:
        capture = cast("RawMcmcCapture", self.raw_capture)
        return capture.stage

    def validate_integrity(self) -> None:
        capture = cast("RawMcmcCapture", self.raw_capture)
        capture.validate_integrity()
        validation = cast("McmcChainValidation", self.validation)
        validation.validate_integrity()
        backend = cast(
            "BackendTransitionEvidence",
            self.backend_transition_evidence,
        )
        backend.validate_integrity()
        if backend.source_capture is not capture:
            raise McmcConstructionError(
                "Operation backend diagnostics differ from its raw capture"
            )
        if self.evidence is not None:
            self.evidence.validate_integrity()
            if self.evidence.backend_transition_evidence is not backend:
                raise McmcConstructionError(
                    "Operation backend diagnostics differ from primary attachment"
                )
        if self._content_identity() != self.identity:
            raise McmcConstructionError("MCMC operation content identity mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "operation",
            self.identity,
        ):
            raise McmcConstructionError(
                "MCMC operation occurrence is not authoritative"
            )


def _mint_mcmc_operation(
    terminal: McmcOperationTerminal,
    evidence: McmcEvidence | None,
    failure_category: str | None,
    failure_message: str,
    raw_capture: RawMcmcCapture,
    validation: McmcChainValidation,
    backend_transition_evidence: BackendTransitionEvidence,
) -> McmcOperation:
    return McmcOperation(
        terminal,
        evidence,
        failure_category,
        failure_message,
        raw_capture,
        validation,
        backend_transition_evidence,
        _mint_mcmc_evidence_witness("operation"),
    )


@dataclass(frozen=True, slots=True)
class RetainedSampleView:
    """Labeled step-major/walker-minor view derived from primary topology."""

    source_evidence_identity: str
    policy_identity: str
    coordinate_ids: tuple[str, ...]
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    selected_state_ordinals: tuple[int, ...]
    selected_walker_ordinals: tuple[int, ...]
    sample_indices: tuple[tuple[int, int], ...]
    samples: tuple[tuple[float, ...], ...]
    log_densities: tuple[float, ...]
    is_complete: bool
    source_evidence: McmcEvidence = field(repr=False, compare=False)
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self._occurrence_witness is None:
            raise McmcConstructionError(
                "Retained MCMC selection must be derived from primary evidence"
            )
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "retained-selection",
            self.identity,
        ):
            raise McmcConstructionError(
                "Retained MCMC selection occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-retained-sample-view",
            (
                self.source_evidence_identity,
                self.policy_identity,
                self.coordinate_ids,
                self.coordinate_units,
                self.selected_state_ordinals,
                self.selected_walker_ordinals,
                self.sample_indices,
                tuple(
                    tuple(float(value).hex() for value in sample)
                    for sample in self.samples
                ),
                tuple(float(value).hex() for value in self.log_densities),
                self.is_complete,
            ),
        )

    def validate_integrity(self) -> None:
        self.source_evidence.validate_integrity()
        if self._content_identity() != self.identity:
            raise McmcConstructionError("Retained MCMC selection content mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "retained-selection",
            self.identity,
        ):
            raise McmcConstructionError(
                "Retained MCMC selection is not its canonical derived occurrence"
            )


@dataclass(frozen=True, slots=True)
class AcceptanceDiagnostics:
    """Canonical diagnostics derived from one exact backend mask source."""

    source: BackendTransitionEvidence = field(repr=False, compare=False)
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    source_identity: str = field(init=False)
    state_ordinals: tuple[int, ...] = field(init=False)
    walker_ordinals: tuple[int, ...] = field(init=False)
    observed_transition_count: int = field(init=False)
    accepted_counts: tuple[int, ...] = field(init=False)
    status: McmcDiagnosticStatus = field(init=False)
    reason: McmcDiagnosticReason | None = field(init=False)
    acceptance_fractions: tuple[float, ...] | None = field(init=False)
    mean_acceptance_fraction: float | None = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            not isinstance(self.source, BackendTransitionEvidence)
            or self._occurrence_witness is None
        ):
            raise McmcConstructionError(
                "Acceptance diagnostics require canonical backend derivation"
            )
        self.source.validate_integrity()
        walker_ordinals = tuple(range(self.source.source_capture.walkers))
        transitions = self.source.transitions
        denominator = len(transitions)
        accepted_counts = tuple(
            sum(transition.accepted[walker] for transition in transitions)
            for walker in walker_ordinals
        )
        if denominator == 0:
            status = McmcDiagnosticStatus.UNAVAILABLE
            reason = McmcDiagnosticReason.NO_TRANSITIONS
            fractions = None
            mean = None
        elif not self.source.supports_observed_diagnostics:
            status = McmcDiagnosticStatus.UNAVAILABLE
            reason = McmcDiagnosticReason.UNVALIDATED_BACKEND_TRANSITIONS
            fractions = None
            mean = None
        elif not self.source.is_structurally_consistent:
            status = McmcDiagnosticStatus.UNAVAILABLE
            reason = McmcDiagnosticReason.BACKEND_TRANSITION_INCONSISTENT
            fractions = None
            mean = None
        else:
            status = McmcDiagnosticStatus.AVAILABLE
            reason = None
            fractions = tuple(count / denominator for count in accepted_counts)
            mean = sum(fractions) / len(fractions)
        object.__setattr__(self, "source_identity", self.source.identity)
        object.__setattr__(
            self,
            "state_ordinals",
            tuple(item.ordinal for item in transitions),
        )
        object.__setattr__(self, "walker_ordinals", walker_ordinals)
        object.__setattr__(self, "observed_transition_count", denominator)
        object.__setattr__(self, "accepted_counts", accepted_counts)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "reason", reason)
        object.__setattr__(self, "acceptance_fractions", fractions)
        object.__setattr__(self, "mean_acceptance_fraction", mean)
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "acceptance-diagnostics",
            self.identity,
        ):
            raise McmcConstructionError(
                "Acceptance diagnostic occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-acceptance-diagnostics",
            (
                self.source_identity,
                self.state_ordinals,
                self.walker_ordinals,
                self.observed_transition_count,
                self.accepted_counts,
                self.status.value,
                None if self.reason is None else self.reason.value,
                None
                if self.acceptance_fractions is None
                else tuple(float(value).hex() for value in self.acceptance_fractions),
                None
                if self.mean_acceptance_fraction is None
                else float(self.mean_acceptance_fraction).hex(),
            ),
        )

    @classmethod
    def from_backend_evidence(
        cls,
        source: BackendTransitionEvidence,
    ) -> AcceptanceDiagnostics:
        return cls(
            source,
            _mint_mcmc_evidence_witness("acceptance-diagnostics"),
        )

    def validate_integrity(self) -> None:
        if not isinstance(self.source, BackendTransitionEvidence):
            raise McmcConstructionError(
                "Acceptance diagnostics require an exact backend source object"
            )
        self.source.validate_integrity()
        if self.source_identity != self.source.identity:
            raise McmcConstructionError(
                "Acceptance diagnostic source identity differs from its source object"
            )
        expected = self.from_backend_evidence(self.source)
        if (
            self.state_ordinals != expected.state_ordinals
            or self.walker_ordinals != expected.walker_ordinals
            or self.observed_transition_count != expected.observed_transition_count
            or self.accepted_counts != expected.accepted_counts
            or self.status is not expected.status
            or self.reason is not expected.reason
            or self.acceptance_fractions != expected.acceptance_fractions
            or self.mean_acceptance_fraction != expected.mean_acceptance_fraction
            or self._content_identity() != self.identity
            or expected._content_identity() != self.identity
        ):
            raise McmcConstructionError("Acceptance diagnostic content mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "acceptance-diagnostics",
            self.identity,
        ):
            raise McmcConstructionError(
                "Acceptance diagnostics are not their canonical occurrence"
            )

    def to_record(self) -> dict[str, object]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "source_identity": self.source_identity,
            "state_ordinals": list(self.state_ordinals),
            "walker_ordinals": list(self.walker_ordinals),
            "observed_transition_count": self.observed_transition_count,
            "accepted_counts": list(self.accepted_counts),
            "status": self.status.value,
            "reason": None if self.reason is None else self.reason.value,
            "acceptance_fractions": None
            if self.acceptance_fractions is None
            else [float(value).hex() for value in self.acceptance_fractions],
            "mean_acceptance_fraction": None
            if self.mean_acceptance_fraction is None
            else float(self.mean_acceptance_fraction).hex(),
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        *,
        source: BackendTransitionEvidence,
    ) -> AcceptanceDiagnostics:
        expected = cls.from_backend_evidence(source)
        if record != expected.to_record():
            raise McmcConstructionError(
                "Stored acceptance diagnostics differ from canonical backend evidence"
            )
        return expected


McmcDiagnostics = AcceptanceDiagnostics


class PosteriorSampleDisposition(StrEnum):
    """Scope-atomic outcome for one retained state/walker label."""

    SUCCESS = "success"
    FAILED = "failed"


@dataclass(frozen=True, slots=True)
class PosteriorSampleFailure:
    """Typed exclusion of one entire requested posterior row."""

    category: str
    message: str


@dataclass(frozen=True, slots=True)
class PosteriorSampleOutcome:
    """One labeled complete-scope success or typed whole-row failure."""

    state_ordinal: int
    walker_ordinal: int
    disposition: PosteriorSampleDisposition
    independent_items: tuple[tuple[str, float], ...]
    resolved_items: tuple[tuple[str, float], ...] | None
    failure: PosteriorSampleFailure | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        success = self.disposition is PosteriorSampleDisposition.SUCCESS
        if success != (self.resolved_items is not None) or success == (
            self.failure is not None
        ):
            raise McmcConstructionError(
                "Posterior sample must contain exactly its scope-atomic payload"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-posterior-sample-outcome",
                (
                    self.state_ordinal,
                    self.walker_ordinal,
                    self.disposition.value,
                    tuple(
                        (key, float(value).hex())
                        for key, value in self.independent_items
                    ),
                    None
                    if self.resolved_items is None
                    else tuple(
                        (key, float(value).hex()) for key, value in self.resolved_items
                    ),
                    None
                    if self.failure is None
                    else (self.failure.category, self.failure.message),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class PosteriorSampleEvidence:
    """All retained labels resolved atomically to one requested output scope."""

    selection: RetainedSampleView = field(repr=False, compare=False)
    output_units: tuple[tuple[str, ParameterUnit], ...]
    outcomes: tuple[PosteriorSampleOutcome, ...]
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    selection_identity: str = field(init=False)
    output_scope: tuple[str, ...] = field(init=False)
    successful_labels: tuple[tuple[int, int], ...] = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self._occurrence_witness is None:
            raise McmcConstructionError(
                "Posterior sample evidence must come from canonical resolution"
            )
        units = tuple(self.output_units)
        scope = tuple(param_id for param_id, _unit in units)
        if (
            not scope
            or len(set(scope)) != len(scope)
            or any(not isinstance(unit, ParameterUnit) for _id, unit in units)
            or tuple(
                (outcome.state_ordinal, outcome.walker_ordinal)
                for outcome in self.outcomes
            )
            != self.selection.sample_indices
        ):
            raise McmcConstructionError(
                "Posterior evidence does not exactly cover its retained labels"
            )
        for outcome, sample in zip(
            self.outcomes,
            self.selection.samples,
            strict=True,
        ):
            if (
                tuple(key for key, _value in outcome.independent_items)
                != self.selection.coordinate_ids
                or tuple(value for _key, value in outcome.independent_items) != sample
                or (
                    outcome.resolved_items is not None
                    and tuple(key for key, _value in outcome.resolved_items) != scope
                )
            ):
                raise McmcConstructionError(
                    "Posterior outcome differs from its exact labeled complete scope"
                )
        successful = tuple(
            (item.state_ordinal, item.walker_ordinal)
            for item in self.outcomes
            if item.disposition is PosteriorSampleDisposition.SUCCESS
        )
        object.__setattr__(self, "output_units", units)
        object.__setattr__(self, "selection_identity", self.selection.identity)
        object.__setattr__(self, "output_scope", scope)
        object.__setattr__(self, "successful_labels", successful)
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "posterior-sample-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Posterior sample evidence occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-posterior-sample-evidence",
            (
                self.selection.identity,
                self.output_units,
                tuple(item.identity for item in self.outcomes),
            ),
        )

    def validate_integrity(self) -> None:
        self.selection.validate_integrity()
        for outcome in self.outcomes:
            expected = PosteriorSampleOutcome(
                outcome.state_ordinal,
                outcome.walker_ordinal,
                outcome.disposition,
                outcome.independent_items,
                outcome.resolved_items,
                outcome.failure,
            )
            if expected.identity != outcome.identity:
                raise McmcConstructionError(
                    "Posterior sample outcome content identity mismatch"
                )
        if self._content_identity() != self.identity:
            raise McmcConstructionError("Posterior sample evidence content mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "posterior-sample-evidence",
            self.identity,
        ):
            raise McmcConstructionError(
                "Posterior sample evidence is not its canonical resolved occurrence"
            )


def derive_posterior_sample_evidence(
    selection: RetainedSampleView,
    output_units: Sequence[tuple[str, ParameterUnit]],
    *,
    cancellation: CancellationToken | None = None,
) -> PosteriorSampleEvidence:
    """Resolve held and constrained values for every retained labeled sample."""
    selection.validate_integrity()
    if cancellation is not None and cancellation.is_cancelled:
        raise McmcConstructionError("Posterior sample derivation was cancelled")
    plan = selection.source_evidence.plan
    units = tuple(output_units)
    scope = tuple(param_id for param_id, _unit in units)
    if not scope or any(
        param_id not in plan.parameterization.scope_ids for param_id in scope
    ):
        raise McmcConstructionError(
            "Posterior output scope must be an ordered subset of the active scope"
        )
    outcomes: list[PosteriorSampleOutcome] = []
    for label, sample in zip(selection.sample_indices, selection.samples, strict=True):
        if cancellation is not None and cancellation.is_cancelled:
            raise McmcConstructionError("Posterior sample derivation was cancelled")
        independent_items = tuple(zip(plan.coordinate_ids, sample, strict=True))
        try:
            frame = plan.source_problem.lifecycle_frame(sample, plan.parameterization)
            resolved = plan.parameterization.resolve(frame)
            resolved_items = tuple(
                (param_id, float(resolved[param_id])) for param_id in scope
            )
            if any(not math.isfinite(value) for _key, value in resolved_items):
                raise ValueError("resolved posterior value is non-finite")
            outcome = PosteriorSampleOutcome(
                label[0],
                label[1],
                PosteriorSampleDisposition.SUCCESS,
                independent_items,
                resolved_items,
                None,
            )
        except Exception as error:  # noqa: BLE001 - typed whole-sample evidence
            outcome = PosteriorSampleOutcome(
                label[0],
                label[1],
                PosteriorSampleDisposition.FAILED,
                independent_items,
                None,
                PosteriorSampleFailure(type(error).__name__, str(error)),
            )
        outcomes.append(outcome)
    evidence = PosteriorSampleEvidence(
        selection,
        units,
        tuple(outcomes),
        _mint_mcmc_evidence_witness("posterior-sample-evidence"),
    )
    return evidence


@dataclass(frozen=True, slots=True)
class PosteriorSummaryPolicy:
    """Explicit equal-tailed complete-case summary convention."""

    quantiles: tuple[float, ...] = (0.025, 0.5, 0.975)
    covariance_delta_degrees_of_freedom: int = 1
    minimum_autocorrelation_steps: int = 16
    minimum_autocorrelation_walkers: int = 2
    autocorrelation_window_parameter: float = 5.0
    autocorrelation_adequacy_tolerance: float = 50.0
    version: str = "mcmc-complete-scope-equal-tailed-v2"
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        quantiles = tuple(float(value) for value in self.quantiles)
        minimum_steps = _positive_integer(
            self.minimum_autocorrelation_steps,
            name="Minimum autocorrelation retained steps",
        )
        minimum_walkers = _positive_integer(
            self.minimum_autocorrelation_walkers,
            name="Minimum autocorrelation retained walkers",
        )
        window_parameter = _finite_positive(
            self.autocorrelation_window_parameter,
            name="Autocorrelation window parameter",
        )
        adequacy_tolerance = _finite_positive(
            self.autocorrelation_adequacy_tolerance,
            name="Autocorrelation adequacy tolerance",
        )
        if (
            not quantiles
            or any(not math.isfinite(value) for value in quantiles)
            or tuple(sorted(set(quantiles))) != quantiles
            or quantiles[0] < 0.0
            or quantiles[-1] > 1.0
            or len(quantiles) != 3
            or quantiles[1] != 0.5
            or not math.isclose(
                quantiles[0],
                1.0 - quantiles[2],
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            or self.covariance_delta_degrees_of_freedom != 1
            or minimum_steps < 2
        ):
            raise McmcConstructionError("Posterior summary policy is invalid")
        object.__setattr__(self, "quantiles", quantiles)
        object.__setattr__(self, "minimum_autocorrelation_steps", minimum_steps)
        object.__setattr__(self, "minimum_autocorrelation_walkers", minimum_walkers)
        object.__setattr__(
            self,
            "autocorrelation_window_parameter",
            window_parameter,
        )
        object.__setattr__(
            self,
            "autocorrelation_adequacy_tolerance",
            adequacy_tolerance,
        )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-posterior-summary-policy",
                (
                    tuple(float(value).hex() for value in quantiles),
                    self.covariance_delta_degrees_of_freedom,
                    minimum_steps,
                    minimum_walkers,
                    float(window_parameter).hex(),
                    float(adequacy_tolerance).hex(),
                    self.version,
                ),
            ),
        )

    def validate_integrity(self) -> None:
        expected = PosteriorSummaryPolicy(
            self.quantiles,
            self.covariance_delta_degrees_of_freedom,
            self.minimum_autocorrelation_steps,
            self.minimum_autocorrelation_walkers,
            self.autocorrelation_window_parameter,
            self.autocorrelation_adequacy_tolerance,
            self.version,
        )
        if expected.identity != self.identity:
            raise McmcConstructionError(
                "Posterior summary policy content identity mismatch"
            )


@dataclass(frozen=True, slots=True)
class PosteriorScalarEstimate:
    """One available estimate or explicit typed unavailability."""

    status: McmcDiagnosticStatus
    reason: McmcDiagnosticReason | None
    value: float | None

    def __post_init__(self) -> None:
        available = self.status is McmcDiagnosticStatus.AVAILABLE
        if available != (self.value is not None) or available == (
            self.reason is not None
        ):
            raise McmcConstructionError(
                "Posterior estimate availability payload is inconsistent"
            )
        if self.value is not None and not math.isfinite(self.value):
            raise McmcConstructionError("Posterior estimate must be finite")


@dataclass(frozen=True, slots=True)
class PosteriorParameterSummary:
    parameter_id: str
    unit: ParameterUnit
    quantiles: tuple[tuple[float, float], ...]
    credible_interval_name: str
    credible_interval: tuple[float, float]
    posterior_standard_deviation: float
    autocorrelation_time: PosteriorScalarEstimate
    effective_sample_size: PosteriorScalarEstimate
    monte_carlo_standard_error: PosteriorScalarEstimate


@dataclass(frozen=True, slots=True)
class PosteriorMatrixEntry:
    parameter_a: str
    parameter_b: str
    value: float


@dataclass(frozen=True, slots=True)
class PosteriorCorrelationEntry:
    parameter_a: str
    parameter_b: str
    estimate: PosteriorScalarEstimate


def _posterior_estimate_record(item: PosteriorScalarEstimate) -> object:
    return (
        item.status.value,
        None if item.reason is None else item.reason.value,
        None if item.value is None else float(item.value).hex(),
    )


@dataclass(frozen=True, slots=True)
class PosteriorSummary:
    """Complete-case posterior summary over one exact successful sample set."""

    source: PosteriorSampleEvidence = field(repr=False, compare=False)
    policy: PosteriorSummaryPolicy = field(repr=False, compare=False)
    included_labels: tuple[tuple[int, int], ...]
    excluded_labels: tuple[tuple[int, int], ...]
    parameter_summaries: tuple[PosteriorParameterSummary, ...]
    covariance: tuple[PosteriorMatrixEntry, ...]
    correlations: tuple[PosteriorCorrelationEntry, ...]
    acceptance: McmcDiagnostics
    _occurrence_witness: _McmcEvidenceWitness | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    source_identity: str = field(init=False)
    policy_identity: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self._occurrence_witness is None:
            raise McmcConstructionError(
                "Posterior summary must come from canonical complete-case derivation"
            )
        object.__setattr__(self, "source_identity", self.source.identity)
        object.__setattr__(self, "policy_identity", self.policy.identity)
        object.__setattr__(self, "identity", self._content_identity())
        if not _bind_mcmc_evidence_witness(
            self._occurrence_witness,
            self,
            "posterior-summary",
            self.identity,
        ):
            raise McmcConstructionError(
                "Posterior summary occurrence witness is already bound"
            )

    def _content_identity(self) -> str:
        return _identity(
            "native-mcmc-posterior-summary",
            (
                self.source.identity,
                self.policy.identity,
                self.included_labels,
                self.excluded_labels,
                tuple(
                    (
                        item.parameter_id,
                        item.unit.value,
                        tuple(
                            (float(probability).hex(), float(value).hex())
                            for probability, value in item.quantiles
                        ),
                        item.credible_interval_name,
                        tuple(float(value).hex() for value in item.credible_interval),
                        float(item.posterior_standard_deviation).hex(),
                        _posterior_estimate_record(item.autocorrelation_time),
                        _posterior_estimate_record(item.effective_sample_size),
                        _posterior_estimate_record(item.monte_carlo_standard_error),
                    )
                    for item in self.parameter_summaries
                ),
                tuple(
                    (item.parameter_a, item.parameter_b, float(item.value).hex())
                    for item in self.covariance
                ),
                tuple(
                    (
                        item.parameter_a,
                        item.parameter_b,
                        _posterior_estimate_record(item.estimate),
                    )
                    for item in self.correlations
                ),
                self.acceptance.identity,
            ),
        )

    def validate_integrity(self) -> None:
        self.source.validate_integrity()
        self.policy.validate_integrity()
        self.acceptance.validate_integrity()
        expected_acceptance = derive_mcmc_diagnostics(
            self.source.selection.source_evidence
        )
        if (
            expected_acceptance != self.acceptance
            or self._content_identity() != self.identity
        ):
            raise McmcConstructionError("Posterior summary content mismatch")
        if not _mcmc_evidence_occurrence_is_authoritative(
            self._occurrence_witness,
            self,
            "posterior-summary",
            self.identity,
        ):
            raise McmcConstructionError(
                "Posterior summary is not its canonical derived occurrence"
            )


class PosteriorSummaryTerminal(StrEnum):
    """Closed result of one posterior summary request."""

    COMPLETED = "completed"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True, slots=True)
class PosteriorSummaryFailure:
    """Typed reason no complete posterior summary was available."""

    source_identity: str
    category: str
    message: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.source_identity or not self.category or not self.message:
            raise McmcConstructionError(
                "Posterior summary failure requires source, category, and message"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-posterior-summary-failure",
                (self.source_identity, self.category, self.message),
            ),
        )


@dataclass(frozen=True, slots=True)
class PosteriorSummaryOutcome:
    """Exactly one canonical summary or typed unavailability."""

    terminal: PosteriorSummaryTerminal
    summary: PosteriorSummary | None = None
    failure: PosteriorSummaryFailure | None = None

    def __post_init__(self) -> None:
        completed = self.terminal is PosteriorSummaryTerminal.COMPLETED
        if completed != (self.summary is not None) or completed == (
            self.failure is not None
        ):
            raise McmcConstructionError(
                "Posterior summary outcome has an inconsistent terminal payload"
            )


def _unavailable(reason: McmcDiagnosticReason) -> PosteriorScalarEstimate:
    return PosteriorScalarEstimate(McmcDiagnosticStatus.UNAVAILABLE, reason, None)


def derive_posterior_summary(  # noqa: C901 - joint typed summary assembly
    source: PosteriorSampleEvidence,
    policy: PosteriorSummaryPolicy | None = None,
    *,
    cancellation: CancellationToken | None = None,
) -> PosteriorSummary:
    """Summarize only the exact common complete successful posterior rows."""
    source.validate_integrity()
    if cancellation is not None and cancellation.is_cancelled:
        raise McmcConstructionError("Posterior summary derivation was cancelled")
    exact_policy = PosteriorSummaryPolicy() if policy is None else policy
    exact_policy.validate_integrity()
    successful = tuple(
        outcome
        for outcome in source.outcomes
        if outcome.disposition is PosteriorSampleDisposition.SUCCESS
    )
    if len(successful) < 2:
        raise McmcConstructionError(
            "Posterior summary requires at least two complete successful samples"
        )
    matrix = np.asarray(
        [
            [
                value
                for _key, value in cast(
                    "tuple[tuple[str, float], ...]", item.resolved_items
                )
            ]
            for item in successful
        ],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(matrix)):
        raise McmcConstructionError("Posterior complete-case matrix is non-finite")
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            covariance_matrix = np.cov(matrix, rowvar=False, ddof=1)
    covariance_matrix = np.atleast_2d(covariance_matrix)
    if not np.all(np.isfinite(covariance_matrix)):
        raise McmcConstructionError("Posterior covariance is unavailable")
    variances = np.diag(covariance_matrix)
    labels = tuple((item.state_ordinal, item.walker_ordinal) for item in successful)
    complete_grid = len(successful) == len(source.outcomes)
    step_count = len(source.selection.selected_state_ordinals)
    autocorrelation: tuple[PosteriorScalarEstimate, ...]
    ess: tuple[PosteriorScalarEstimate, ...]
    mcse: tuple[PosteriorScalarEstimate, ...]
    if (
        not source.selection.is_complete
        or not complete_grid
        or step_count < exact_policy.minimum_autocorrelation_steps
        or len(source.selection.selected_walker_ordinals)
        < exact_policy.minimum_autocorrelation_walkers
    ):
        reason = (
            McmcDiagnosticReason.PARTIAL_CHAIN
            if not source.selection.is_complete
            else McmcDiagnosticReason.INSUFFICIENT_RETAINED_WINDOW
        )
        autocorrelation = tuple(_unavailable(reason) for _ in source.output_scope)
        ess = tuple(_unavailable(reason) for _ in source.output_scope)
        mcse = tuple(_unavailable(reason) for _ in source.output_scope)
    else:
        chain = matrix.reshape(
            step_count,
            len(source.selection.selected_walker_ordinals),
            len(source.output_scope),
        )
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                tau_values = np.asarray(
                    emcee.autocorr.integrated_time(
                        chain,
                        c=exact_policy.autocorrelation_window_parameter,
                        tol=exact_policy.autocorrelation_adequacy_tolerance,
                        quiet=False,
                    ),
                    dtype=np.float64,
                )
        except emcee.autocorr.AutocorrError as error:
            tau_values = np.asarray(error.tau, dtype=np.float64)
        except Exception:  # noqa: BLE001 - typed unavailable estimate
            tau_values = np.full(len(source.output_scope), np.nan)
        if tau_values.shape != (len(source.output_scope),):
            tau_values = np.full(len(source.output_scope), np.nan)
        autocorrelation_values: list[PosteriorScalarEstimate] = []
        ess_values: list[PosteriorScalarEstimate] = []
        mcse_values: list[PosteriorScalarEstimate] = []
        for index, tau in enumerate(tau_values):
            if not math.isfinite(float(tau)) or tau <= 0.0:
                unavailable = _unavailable(McmcDiagnosticReason.INVALID_AUTOCORRELATION)
                autocorrelation_values.append(unavailable)
                ess_values.append(unavailable)
                mcse_values.append(unavailable)
                continue
            if (
                exact_policy.autocorrelation_adequacy_tolerance * float(tau)
                > step_count
            ):
                unavailable = _unavailable(
                    McmcDiagnosticReason.UNRELIABLE_AUTOCORRELATION
                )
                autocorrelation_values.append(unavailable)
                ess_values.append(unavailable)
                mcse_values.append(unavailable)
                continue
            effective = float(matrix.shape[0] / tau)
            standard_deviation = math.sqrt(float(variances[index]))
            autocorrelation_values.append(
                PosteriorScalarEstimate(
                    McmcDiagnosticStatus.AVAILABLE,
                    None,
                    float(tau),
                )
            )
            ess_values.append(
                PosteriorScalarEstimate(
                    McmcDiagnosticStatus.AVAILABLE,
                    None,
                    effective,
                )
            )
            mcse_values.append(
                PosteriorScalarEstimate(
                    McmcDiagnosticStatus.AVAILABLE,
                    None,
                    standard_deviation / math.sqrt(effective),
                )
            )
        autocorrelation = tuple(autocorrelation_values)
        ess = tuple(ess_values)
        mcse = tuple(mcse_values)
    summaries: list[PosteriorParameterSummary] = []
    for index, (param_id, unit) in enumerate(source.output_units):
        values = matrix[:, index]
        quantile_values = np.quantile(values, exact_policy.quantiles, method="linear")
        quantiles = tuple(
            (probability, float(value))
            for probability, value in zip(
                exact_policy.quantiles,
                quantile_values,
                strict=True,
            )
        )
        summaries.append(
            PosteriorParameterSummary(
                param_id,
                unit,
                quantiles,
                f"equal_tailed_{100.0 * (quantiles[-1][0] - quantiles[0][0]):g}_percent",
                (quantiles[0][1], quantiles[-1][1]),
                float(np.std(values, ddof=1)),
                autocorrelation[index],
                ess[index],
                mcse[index],
            )
        )
    covariance = tuple(
        PosteriorMatrixEntry(param_a, param_b, float(covariance_matrix[row, column]))
        for row, param_a in enumerate(source.output_scope)
        for column, param_b in enumerate(source.output_scope)
    )
    correlations: list[PosteriorCorrelationEntry] = []
    for row, param_a in enumerate(source.output_scope):
        for column, param_b in enumerate(source.output_scope):
            denominator = math.sqrt(float(variances[row])) * math.sqrt(
                float(variances[column])
            )
            estimate = (
                _unavailable(McmcDiagnosticReason.ZERO_VARIANCE)
                if denominator == 0.0
                else PosteriorScalarEstimate(
                    McmcDiagnosticStatus.AVAILABLE,
                    None,
                    float(covariance_matrix[row, column] / denominator),
                )
            )
            correlations.append(PosteriorCorrelationEntry(param_a, param_b, estimate))
    summary = PosteriorSummary(
        source,
        exact_policy,
        labels,
        tuple(
            (item.state_ordinal, item.walker_ordinal)
            for item in source.outcomes
            if item.disposition is PosteriorSampleDisposition.FAILED
        ),
        tuple(summaries),
        covariance,
        tuple(correlations),
        derive_mcmc_diagnostics(source.selection.source_evidence),
        _mint_mcmc_evidence_witness("posterior-summary"),
    )
    return summary


def summarize_posterior_evidence(
    source: PosteriorSampleEvidence,
    policy: PosteriorSummaryPolicy | None = None,
    *,
    cancellation: CancellationToken | None = None,
) -> PosteriorSummaryOutcome:
    """Return a summary or durable typed unavailability without fabrication."""
    try:
        summary = derive_posterior_summary(
            source,
            policy,
            cancellation=cancellation,
        )
    except McmcConstructionError as error:
        message = str(error)
        category = (
            "insufficient_complete_samples"
            if "at least two" in message
            else "cancelled"
            if "cancelled" in message
            else "invalid_source_or_summary"
        )
        return PosteriorSummaryOutcome(
            PosteriorSummaryTerminal.UNAVAILABLE,
            failure=PosteriorSummaryFailure(source.identity, category, message),
        )
    except (FloatingPointError, OverflowError, RuntimeWarning, ValueError) as error:
        return PosteriorSummaryOutcome(
            PosteriorSummaryTerminal.UNAVAILABLE,
            failure=PosteriorSummaryFailure(
                source.identity,
                "numerical_summary_unavailable",
                f"{type(error).__name__}: {error}",
            ),
        )
    return PosteriorSummaryOutcome(
        PosteriorSummaryTerminal.COMPLETED,
        summary=summary,
    )


class _LogDensityEvaluator:
    def __init__(self, plan: McmcPlan) -> None:
        self.plan = plan
        self.evaluator = plan.source_engine.new_evaluator()
        self.objective_request_count = 0
        self.evaluation_request_count = 0

    def __call__(self, vector: Array) -> float:
        if self.objective_request_count >= self.plan.policy.objective_request_budget:
            raise McmcExecutionError("MCMC objective-request budget exhausted")
        self.objective_request_count += 1
        candidate = np.asarray(vector, dtype=np.float64)
        if (
            candidate.shape != (self.plan.policy.dimension,)
            or not np.all(np.isfinite(candidate))
            or np.any(candidate <= np.asarray(self.plan.lower_bounds))
            or np.any(candidate >= np.asarray(self.plan.upper_bounds))
        ):
            return -np.inf
        try:
            lifecycle = self.plan.source_problem.lifecycle_frame(
                candidate,
                self.plan.parameterization,
            )
            frame = EvaluationFrame.from_lifecycle_frame(
                self.plan.parameterization,
                lifecycle,
            )
        except Exception as error:
            raise McmcExecutionError(
                f"MCMC candidate frame failed: {type(error).__name__}: {error}"
            ) from error
        outcome = self.evaluator.evaluate(frame)
        self.evaluation_request_count += 1
        if isinstance(outcome, EvaluationFailure):
            if outcome.validity == "INVALID_TRIAL":
                return -np.inf
            raise McmcExecutionError(
                f"MCMC native evaluation failed: {outcome.category}: {outcome.message}"
            )
        return -0.5 * canonical_chi_square(outcome.residuals)


class _RecordingStretchMove(emcee.moves.StretchMove):
    """Capture the backend's exact per-transition acceptance mask."""

    def __init__(self, *, scale: float) -> None:
        super().__init__(a=scale)
        self.last_accepted: Array | None = None

    def propose(self, model: object, state: object) -> tuple[object, Array]:
        proposed, accepted = super().propose(model, state)
        self.last_accepted = np.asarray(accepted, dtype=np.bool_).copy()
        return proposed, self.last_accepted


def _ensemble_state(
    ordinal: int,
    positions: Array,
    log_densities: Array,
) -> EnsembleState:
    current = np.asarray(positions, dtype=np.float64)
    return EnsembleState(
        ordinal,
        tuple(tuple(float(value) for value in row) for row in current),
        tuple(float(value) for value in np.asarray(log_densities, dtype=np.float64)),
    )


def _legacy_sampler_random_state(seed: int) -> object:
    keys = np.random.SeedSequence(seed).generate_state(624, dtype=np.uint32)
    return ("MT19937", keys, 624, 0, 0.0)


def execute_mcmc_evidence(  # noqa: C901 - checkpointed lifecycle boundary
    accepted: AcceptedFitResult,
    plan: McmcPlan,
    *,
    state_observer: Callable[[EnsembleState], None] | None = None,
    cancellation: CancellationToken | None = None,
    checkpoint_observer: Callable[[McmcExecutionStage, int], None] | None = None,
) -> McmcOperation:
    """Execute one fixed run while preserving only complete ensemble states."""
    plan.validate_integrity(accepted)
    signal = CancellationToken() if cancellation is None else cancellation
    log_density = _LogDensityEvaluator(plan)
    backend_observation = _mint_backend_execution_observation(
        plan,
        backend_implementation_identity="emcee-stretch-backend-v1",
    )
    states: list[EnsembleState] = []
    stage = McmcExecutionStage.BEFORE_INITIALIZATION
    initialization_outcome = McmcInitializationOutcome.NOT_STARTED

    def checkpoint(next_stage: McmcExecutionStage) -> None:
        nonlocal stage
        stage = next_stage
        if checkpoint_observer is not None:
            checkpoint_observer(stage, len(states))
        if signal.is_cancelled:
            raise _McmcCancelled(f"MCMC cancelled at {stage.value}")

    try:
        checkpoint(McmcExecutionStage.BEFORE_INITIALIZATION)
        initial = np.asarray(plan.initial_ensemble, dtype=np.float64)
        initial_values: list[float] = []
        for row in initial:
            checkpoint(McmcExecutionStage.INITIALIZING)
            initial_values.append(log_density(row))
        initial_log_density = np.asarray(initial_values, dtype=np.float64)
        if not np.all(np.isfinite(initial_log_density)):
            raise McmcExecutionError(
                "MCMC initial ensemble has unavailable native log density"
            )
        initialization_outcome = McmcInitializationOutcome.COMPLETED
        states.append(_ensemble_state(0, initial, initial_log_density))
        checkpoint(McmcExecutionStage.AFTER_INITIALIZATION)
        checkpoint(McmcExecutionStage.AFTER_COMPLETE_STATE)
        move = _RecordingStretchMove(scale=plan.policy.proposal_scale)
        sampler = emcee.EnsembleSampler(
            plan.policy.walkers,
            plan.policy.dimension,
            log_density,
            moves=move,
        )
        sampler.random_state = _legacy_sampler_random_state(plan.policy.sampler_seed)
        emcee_state = emcee.State(initial, log_prob=initial_log_density)
        iterator = sampler.sample(
            emcee_state,
            iterations=plan.policy.total_steps,
            tune=False,
            skip_initial_state_check=True,
            store=False,
        )
        for ordinal in range(1, plan.policy.total_steps + 1):
            checkpoint(McmcExecutionStage.BEFORE_TRANSITION)
            checkpoint(McmcExecutionStage.DURING_TRANSITION)
            sampled = next(iterator)
            if move.last_accepted is None:
                raise McmcExecutionError(
                    "MCMC backend omitted its transition diagnostic mask"
                )
            state = _ensemble_state(
                ordinal,
                sampled.coords,
                sampled.log_prob,
            )
            _record_backend_execution_transition(
                backend_observation,
                states[-1],
                state,
                tuple(bool(value) for value in move.last_accepted),
            )
            states.append(state)
            if state_observer is not None:
                state_observer(state)
            checkpoint(McmcExecutionStage.AFTER_COMPLETE_STATE)
        if log_density.objective_request_count != plan.policy.objective_request_budget:
            raise McmcExecutionError(
                "MCMC backend request count differs from the frozen budget"
            )
        checkpoint(McmcExecutionStage.FRESH_VALIDATION)
        raw_capture = _mint_raw_capture(
            plan,
            McmcOperationTerminal.COMPLETED,
            tuple(states),
            log_density.objective_request_count,
            log_density.evaluation_request_count,
            stage,
            initialization_outcome,
            backend_observation,
        )
        validation = validate_raw_mcmc_capture(plan, raw_capture)
        if not validation.is_complete or validation.primary_evidence is None:
            raise McmcExecutionError("fresh MCMC validation did not complete")
        checkpoint(McmcExecutionStage.BEFORE_FINAL_ASSEMBLY)
        backend_transition_evidence = _mint_observed_backend_transition_evidence(
            backend_observation,
            raw_capture,
        )
        evidence = validation.primary_evidence
        object.__setattr__(
            evidence,
            "backend_transition_evidence",
            backend_transition_evidence,
        )
        evidence.validate_integrity()
        return _mint_mcmc_operation(
            McmcOperationTerminal.COMPLETED,
            evidence,
            None,
            "",
            raw_capture,
            validation,
            backend_transition_evidence,
        )
    except _McmcCancelled as error:
        terminal = McmcOperationTerminal.CANCELLED
        category = "cancelled"
        message = str(error)
        if initialization_outcome is McmcInitializationOutcome.NOT_STARTED:
            initialization_outcome = McmcInitializationOutcome.CANCELLED
    except KeyboardInterrupt as error:
        terminal = McmcOperationTerminal.INTERRUPTED
        category = "interrupted"
        message = str(error)
        if initialization_outcome is McmcInitializationOutcome.NOT_STARTED:
            initialization_outcome = McmcInitializationOutcome.INTERRUPTED
    except Exception as error:  # noqa: BLE001 - freeze typed operation failure
        terminal = McmcOperationTerminal.FAILED
        category = type(error).__name__
        message = str(error)
        if initialization_outcome is McmcInitializationOutcome.NOT_STARTED:
            initialization_outcome = McmcInitializationOutcome.FAILED
    raw_capture = _mint_raw_capture(
        plan,
        terminal,
        tuple(states),
        log_density.objective_request_count,
        log_density.evaluation_request_count,
        stage,
        initialization_outcome,
        backend_observation,
        category,
    )
    backend_transition_evidence = _mint_observed_backend_transition_evidence(
        backend_observation,
        raw_capture,
    )
    validation = validate_raw_mcmc_capture(
        plan,
        raw_capture,
        backend_transition_evidence=backend_transition_evidence,
    )
    return _mint_mcmc_operation(
        terminal,
        validation.primary_evidence,
        category,
        message,
        raw_capture,
        validation,
        backend_transition_evidence,
    )


def derive_retained_sample_view(
    evidence: McmcEvidence,
    *,
    cancellation: CancellationToken | None = None,
) -> RetainedSampleView:
    """Select the prospective retained window without flattening primary evidence."""
    evidence.validate_integrity()
    if cancellation is not None and cancellation.is_cancelled:
        raise McmcConstructionError("Retained MCMC selection was cancelled")
    first = evidence.plan.policy.burn_steps + 1
    selected = tuple(state for state in evidence.states if state.ordinal >= first)
    ordinals = tuple(state.ordinal for state in selected)
    sample_indices = tuple(
        (state.ordinal, walker)
        for state in selected
        for walker in range(evidence.plan.policy.walkers)
    )
    samples = tuple(position for state in selected for position in state.positions)
    log_densities = tuple(value for state in selected for value in state.log_densities)
    retained = RetainedSampleView(
        evidence.identity,
        evidence.plan.policy.identity,
        evidence.coordinate_ids,
        evidence.coordinate_units,
        ordinals,
        tuple(range(evidence.plan.policy.walkers)),
        sample_indices,
        samples,
        log_densities,
        len(selected) == evidence.plan.policy.retained_steps,
        evidence,
        _mint_mcmc_evidence_witness("retained-selection"),
    )
    return retained


def derive_mcmc_diagnostics(
    evidence: McmcEvidence,
    *,
    cancellation: CancellationToken | None = None,
) -> McmcDiagnostics:
    """Derive diagnostic-only acceptance observations from their backend source."""
    evidence.validate_integrity()
    if cancellation is not None and cancellation.is_cancelled:
        raise McmcConstructionError("MCMC diagnostic derivation was cancelled")
    source = evidence.backend_transition_evidence
    if (
        source is None
        or source.source_capture_identity != evidence.raw_capture_identity
    ):
        raise McmcConstructionError(
            "Primary evidence has no matching backend transition diagnostic source"
        )
    return AcceptanceDiagnostics.from_backend_evidence(source)


def derive_mcmc_operation_diagnostics(
    operation: McmcOperation,
    *,
    cancellation: CancellationToken | None = None,
) -> McmcDiagnostics:
    """Derive truthful diagnostics even when no primary state exists."""
    operation.validate_integrity()
    if cancellation is not None and cancellation.is_cancelled:
        raise McmcConstructionError("MCMC diagnostic derivation was cancelled")
    source = cast("BackendTransitionEvidence", operation.backend_transition_evidence)
    if operation.evidence is not None:
        evidence_source = operation.evidence.backend_transition_evidence
        if evidence_source is not source:
            raise McmcConstructionError(
                "Operation and primary evidence backend diagnostics differ"
            )
    return AcceptanceDiagnostics.from_backend_evidence(source)
