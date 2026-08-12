"""Prospective native fixed-topology MCMC evidence qualification (#600).

This module is an isolated qualification seam. Production method parsing,
workflow orchestration, reporting, and parameter-state mutation do not consume
these artifacts yet.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import cast

import emcee
import numpy as np

from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
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
    FAILED = "failed"
    INTERRUPTED = "interrupted"


class McmcExecutionError(RuntimeError):
    """Raised internally when native evaluation cannot produce log density."""


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
        object.__setattr__(self, "identity", identity)


@dataclass(frozen=True, slots=True)
class CalibratedMcmcPolicy:
    """One versioned qualified topology for a bounded dimension stratum."""

    policy_version: str
    minimum_dimension: int
    maximum_dimension: int
    walkers: int
    burn_steps: int
    retained_steps: int
    proposal_scale: float = 2.0

    def resolve(self, *, dimension: int, root_seed: int) -> ResolvedMcmcPolicy:
        minimum = _positive_integer(
            self.minimum_dimension,
            name="calibrated minimum dimension",
        )
        maximum = _positive_integer(
            self.maximum_dimension,
            name="calibrated maximum dimension",
        )
        if minimum > maximum:
            raise McmcConstructionError(
                "Calibrated MCMC dimension stratum is strictly inverted"
            )
        if not minimum <= dimension <= maximum:
            raise McmcConstructionError(
                "No calibrated MCMC policy covers the requested dimension"
            )
        return ResolvedMcmcPolicy(
            kind=McmcPolicyKind.CALIBRATED,
            policy_version=_nonempty(
                self.policy_version,
                name="MCMC policy version",
            ),
            dimension=dimension,
            walkers=self.walkers,
            burn_steps=self.burn_steps,
            retained_steps=self.retained_steps,
            root_seed=root_seed,
            qualification_dimension_range=(minimum, maximum),
            proposal_scale=self.proposal_scale,
        )


@dataclass(frozen=True, slots=True)
class ExpertMcmcPolicy:
    """Prospective topology overrides without inherited adequacy claims."""

    burn_steps: int
    retained_steps: int
    walkers: int
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
    accepted_result_identity: str = field(init=False)
    accepted_occurrence_identity: str = field(init=False)
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
        identity = _identity(
            "native-mcmc-plan",
            (
                accepted.identity,
                accepted.occurrence_identity,
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
        object.__setattr__(self, "accepted_result_identity", accepted.identity)
        object.__setattr__(
            self,
            "accepted_occurrence_identity",
            accepted.occurrence_identity,
        )
        object.__setattr__(self, "problem_identity", problem.identity)
        object.__setattr__(self, "coordinate_ids", coordinate_ids)
        object.__setattr__(self, "lower_bounds", lower)
        object.__setattr__(self, "upper_bounds", upper)
        object.__setattr__(self, "initial_ensemble", initial)
        object.__setattr__(self, "identity", identity)

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
    """One atomically completed full-walker state and its log densities."""

    ordinal: int
    positions: tuple[tuple[float, ...], ...]
    log_densities: tuple[float, ...]
    accepted: tuple[bool, ...]
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
        accepted = tuple(self.accepted)
        if (
            len(log_densities) != len(positions)
            or len(accepted) != len(positions)
            or any(not isinstance(value, bool) for value in accepted)
        ):
            raise McmcConstructionError(
                "MCMC state payload must cover every walker exactly once"
            )
        identity = _identity(
            "native-mcmc-ensemble-state",
            (
                self.ordinal,
                tuple(tuple(float(value).hex() for value in row) for row in positions),
                tuple(float(value).hex() for value in log_densities),
                accepted,
            ),
        )
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "log_densities", log_densities)
        object.__setattr__(self, "accepted", accepted)
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
            "accepted": list(self.accepted),
            "identity": self.identity,
        }

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
                "accepted",
                "identity",
            },
            "MCMC ensemble-state",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise McmcConstructionError("Unsupported MCMC ensemble-state schema")
        ordinal = record.get("ordinal")
        raw_positions = record.get("positions")
        raw_log_densities = record.get("log_densities")
        raw_accepted = record.get("accepted")
        identity = record.get("identity")
        if (
            isinstance(ordinal, bool)
            or not isinstance(ordinal, int)
            or not isinstance(raw_positions, list)
            or not all(isinstance(row, list) for row in raw_positions)
            or not isinstance(raw_log_densities, list)
            or not isinstance(raw_accepted, list)
            or not all(isinstance(value, bool) for value in raw_accepted)
            or not isinstance(identity, str)
        ):
            raise McmcConstructionError("Malformed MCMC ensemble-state record")
        position_rows = cast("list[list[object]]", raw_positions)
        accepted_values = cast("list[bool]", raw_accepted)
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
            tuple(accepted_values),
        )
        if state.identity != identity:
            raise McmcConstructionError("MCMC ensemble-state identity mismatch")
        return state


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
    plan_identity: str = field(init=False)
    accepted_result_identity: str = field(init=False)
    coordinate_ids: tuple[str, ...] = field(init=False)
    coordinate_units: tuple[tuple[str, ParameterUnit], ...] = field(init=False)
    completed_transition_count: int = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
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
        object.__setattr__(self, "identity", identity)

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
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        *,
        plan: McmcPlan,
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
        evidence = cls(
            plan,
            lifecycle,
            terminal,
            tuple(
                EnsembleState.from_record(item)
                for item in cast("list[Mapping[str, object]]", raw_states)
            ),
            objective_count,
            evaluation_count,
            failure_category,
        )
        if evidence.identity != identity:
            raise McmcConstructionError("MCMC evidence identity mismatch")
        return evidence


@dataclass(frozen=True, slots=True)
class McmcOperation:
    """Frozen terminal record for one attempted MCMC evidence execution."""

    terminal: McmcOperationTerminal
    evidence: McmcEvidence | None
    failure_category: str | None = None
    failure_message: str = ""

    def __post_init__(self) -> None:
        if self.terminal is McmcOperationTerminal.COMPLETED:
            if self.evidence is None or self.failure_category is not None:
                raise McmcConstructionError(
                    "Completed MCMC operation requires completed evidence"
                )
        elif self.failure_category is None:
            raise McmcConstructionError(
                "Non-completed MCMC operation requires typed failure evidence"
            )


@dataclass(frozen=True, slots=True)
class RetainedSampleView:
    """Labeled step-major/walker-minor view derived from primary topology."""

    source_evidence_identity: str
    coordinate_ids: tuple[str, ...]
    coordinate_units: tuple[tuple[str, ParameterUnit], ...]
    selected_state_ordinals: tuple[int, ...]
    sample_indices: tuple[tuple[int, int], ...]
    samples: tuple[tuple[float, ...], ...]
    log_densities: tuple[float, ...]
    is_complete: bool
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-retained-sample-view",
                (
                    self.source_evidence_identity,
                    self.coordinate_ids,
                    self.coordinate_units,
                    self.selected_state_ordinals,
                    self.sample_indices,
                    tuple(
                        tuple(float(value).hex() for value in sample)
                        for sample in self.samples
                    ),
                    tuple(float(value).hex() for value in self.log_densities),
                    self.is_complete,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class McmcDiagnostics:
    """Acceptance diagnostics derived without altering primary chain evidence."""

    source_evidence_identity: str
    state_ordinals: tuple[int, ...]
    walker_ordinals: tuple[int, ...]
    acceptance_fractions: tuple[float, ...]
    mean_acceptance_fraction: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-mcmc-acceptance-diagnostics",
                (
                    self.source_evidence_identity,
                    self.state_ordinals,
                    self.walker_ordinals,
                    tuple(float(value).hex() for value in self.acceptance_fractions),
                    float(self.mean_acceptance_fraction).hex(),
                ),
            ),
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
            or np.any(candidate < np.asarray(self.plan.lower_bounds))
            or np.any(candidate > np.asarray(self.plan.upper_bounds))
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
    accepted_mask: Array | None,
) -> EnsembleState:
    current = np.asarray(positions, dtype=np.float64)
    accepted = (
        tuple(False for _ in current)
        if accepted_mask is None
        else tuple(bool(value) for value in accepted_mask)
    )
    return EnsembleState(
        ordinal,
        tuple(tuple(float(value) for value in row) for row in current),
        tuple(float(value) for value in np.asarray(log_densities, dtype=np.float64)),
        accepted,
    )


def _legacy_sampler_random_state(seed: int) -> object:
    keys = np.random.SeedSequence(seed).generate_state(624, dtype=np.uint32)
    return ("MT19937", keys, 624, 0, 0.0)


def execute_mcmc_evidence(
    accepted: AcceptedFitResult,
    plan: McmcPlan,
    *,
    state_observer: Callable[[EnsembleState], None] | None = None,
) -> McmcOperation:
    """Execute one fixed run while preserving only complete ensemble states."""
    if accepted.identity != plan.accepted_result_identity:
        raise McmcConstructionError(
            "MCMC execution accepted result differs from its frozen plan"
        )
    log_density = _LogDensityEvaluator(plan)
    states: list[EnsembleState] = []
    try:
        initial = np.asarray(plan.initial_ensemble, dtype=np.float64)
        initial_log_density = np.asarray(
            [log_density(row) for row in initial],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(initial_log_density)):
            raise McmcExecutionError(
                "MCMC initial ensemble has unavailable native log density"
            )
        states.append(_ensemble_state(0, initial, initial_log_density, None))
        move = _RecordingStretchMove(scale=plan.policy.proposal_scale)
        sampler = emcee.EnsembleSampler(
            plan.policy.walkers,
            plan.policy.dimension,
            log_density,
            moves=move,
        )
        sampler.random_state = _legacy_sampler_random_state(plan.policy.sampler_seed)
        emcee_state = emcee.State(initial, log_prob=initial_log_density)
        for ordinal, sampled in enumerate(
            sampler.sample(
                emcee_state,
                iterations=plan.policy.total_steps,
                tune=False,
                skip_initial_state_check=True,
                store=False,
            ),
            start=1,
        ):
            state = _ensemble_state(
                ordinal,
                sampled.coords,
                sampled.log_prob,
                move.last_accepted,
            )
            states.append(state)
            if state_observer is not None:
                state_observer(state)
        if log_density.objective_request_count != plan.policy.objective_request_budget:
            raise McmcExecutionError(
                "MCMC backend request count differs from the frozen budget"
            )
        evidence = McmcEvidence(
            plan,
            McmcEvidenceLifecycle.COMPLETED,
            McmcOperationTerminal.COMPLETED,
            tuple(states),
            log_density.objective_request_count,
            log_density.evaluation_request_count,
        )
        return McmcOperation(McmcOperationTerminal.COMPLETED, evidence)
    except KeyboardInterrupt as error:
        terminal = McmcOperationTerminal.INTERRUPTED
        category = "interrupted"
        message = str(error)
    except Exception as error:  # noqa: BLE001 - freeze typed operation failure
        terminal = McmcOperationTerminal.FAILED
        category = type(error).__name__
        message = str(error)
    evidence = (
        McmcEvidence(
            plan,
            McmcEvidenceLifecycle.PARTIAL,
            terminal,
            tuple(states),
            log_density.objective_request_count,
            log_density.evaluation_request_count,
            category,
        )
        if states
        else None
    )
    return McmcOperation(terminal, evidence, category, message)


def derive_retained_sample_view(evidence: McmcEvidence) -> RetainedSampleView:
    """Select the prospective retained window without flattening primary evidence."""
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
    return RetainedSampleView(
        evidence.identity,
        evidence.coordinate_ids,
        evidence.coordinate_units,
        ordinals,
        sample_indices,
        samples,
        log_densities,
        len(selected) == evidence.plan.policy.retained_steps,
    )


def derive_mcmc_diagnostics(evidence: McmcEvidence) -> McmcDiagnostics:
    """Derive acceptance diagnostics from complete transition states only."""
    transitions = evidence.states[1:]
    completed = len(transitions)
    accepted_counts = tuple(
        sum(state.accepted[walker] for state in transitions)
        for walker in range(evidence.plan.policy.walkers)
    )
    fractions = tuple(
        count / completed if completed else 0.0 for count in accepted_counts
    )
    mean = sum(fractions) / len(fractions)
    return McmcDiagnostics(
        evidence.identity,
        tuple(state.ordinal for state in transitions),
        tuple(range(evidence.plan.policy.walkers)),
        fractions,
        mean,
    )
