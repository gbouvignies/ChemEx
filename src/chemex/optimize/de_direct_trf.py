"""Qualified, non-authoritative selected-coordinate DE search machinery."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from numbers import Real
from typing import Protocol, cast
from uuid import uuid4

import numpy as np
from scipy.optimize import Bounds, differential_evolution

from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.optimize.direct_trf import (
    AttemptCounters,
    CancellationToken,
    DeSearchProblemDerivation,
    DirectTrfConstructionError,
    ObjectiveScalarizationError,
    OptimizationProblem,
    TerminalFailure,
    canonical_chi_square,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.typing import Array

_SCHEMA_VERSION = 1
_DE_BACKEND_POLICY_VERSION = "scipy-de-best1bin-lhs-deferred-pcg64-box-midpoint-v2"
_UINT64_MAX = 2**64 - 1
_PRODUCT_DE_POPULATION_MULTIPLIER = 15
_PRODUCT_DE_MAXIMUM_GENERATIONS = 1000


@dataclass(frozen=True, slots=True)
class _ExpectedSearchProjection:
    held_items: tuple[tuple[str, float], ...]
    start: tuple[float, ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    derivation: DeSearchProblemDerivation


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _finite_binary64(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise DirectTrfConstructionError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise DirectTrfConstructionError(f"{name} must be finite")
    return 0.0 if result == 0.0 else result


def _float_token(value: float) -> str:
    return float(value).hex()


class DeCoordinateSemantics(StrEnum):
    """Closed physical-to-search coordinate maps supported by DE v1."""

    LINEAR = "linear"
    LOG = "log"


@dataclass(frozen=True, slots=True)
class DeSearchCoordinate:
    """One selected independent physical coordinate and finite search range."""

    param_id: str
    physical_lower: float
    physical_upper: float
    semantics: DeCoordinateSemantics | str
    declaration_ordinal: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.param_id:
            raise DirectTrfConstructionError("DE coordinate ID cannot be empty")
        lower = _finite_binary64(self.physical_lower, name="DE search lower bound")
        upper = _finite_binary64(self.physical_upper, name="DE search upper bound")
        if lower >= upper:
            raise DirectTrfConstructionError(
                "DE search coordinates require strictly ordered finite bounds"
            )
        try:
            semantics = DeCoordinateSemantics(self.semantics)
        except ValueError as error:
            raise DirectTrfConstructionError(
                "DE coordinate semantics must be 'linear' or 'log'"
            ) from error
        if semantics is DeCoordinateSemantics.LOG and lower <= 0.0:
            raise DirectTrfConstructionError(
                "Logarithmic DE coordinates require positive physical bounds"
            )
        if (
            isinstance(self.declaration_ordinal, bool)
            or not isinstance(self.declaration_ordinal, int)
            or self.declaration_ordinal < 0
        ):
            raise DirectTrfConstructionError(
                "DE declaration ordinal must be a non-negative integer"
            )
        object.__setattr__(self, "physical_lower", lower)
        object.__setattr__(self, "physical_upper", upper)
        object.__setattr__(self, "semantics", semantics)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-search-coordinate",
                (
                    self.param_id,
                    _float_token(lower),
                    _float_token(upper),
                    semantics.value,
                    self.declaration_ordinal,
                ),
            ),
        )

    @property
    def physical_bounds(self) -> tuple[float, float]:
        return self.physical_lower, self.physical_upper

    @property
    def solver_bounds(self) -> tuple[float, float]:
        if self.semantics is DeCoordinateSemantics.LOG:
            return math.log(self.physical_lower), math.log(self.physical_upper)
        return self.physical_bounds

    @property
    def solver_initial(self) -> float:
        """Return a deterministic in-box initializer without changing captured state."""
        lower, upper = self.solver_bounds
        return math.fsum((0.5 * lower, 0.5 * upper))

    def to_solver(self, physical_value: float) -> float:
        value = _finite_binary64(
            physical_value,
            name=f"DE physical coordinate {self.param_id!r}",
        )
        if self.semantics is DeCoordinateSemantics.LOG:
            if value <= 0.0:
                raise DirectTrfConstructionError(
                    f"Logarithmic DE coordinate for {self.param_id!r} must be positive"
                )
            return math.log(value)
        return value

    def to_physical(self, solver_value: float) -> float:
        """Map one finite solver coordinate to its canonical physical value."""
        value = _finite_binary64(
            solver_value,
            name=f"DE solver coordinate {self.param_id!r}",
        )
        solver_lower, solver_upper = self.solver_bounds
        if not solver_lower <= value <= solver_upper:
            raise ValueError(
                f"DE solver coordinate for {self.param_id!r} is out of bounds"
            )
        if self.semantics is DeCoordinateSemantics.LINEAR:
            return value
        if value == solver_lower:
            return self.physical_lower
        if value == solver_upper:
            return self.physical_upper
        physical = math.exp(value)
        if not math.isfinite(physical):
            raise ValueError(f"DE logarithmic map for {self.param_id!r} is non-finite")
        return 0.0 if physical == 0.0 else physical


@dataclass(frozen=True, slots=True)
class DePopulation:
    """Serializable fixed population topology for one DE search attempt."""

    dimension: int
    multiplier: int
    size: int
    maximum_generations: int
    strategy: str = "best1bin"
    initialization: str = "latinhypercube"
    updating: str = "deferred"
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.dimension, bool)
            or not isinstance(self.dimension, int)
            or isinstance(self.multiplier, bool)
            or not isinstance(self.multiplier, int)
            or isinstance(self.size, bool)
            or not isinstance(self.size, int)
            or isinstance(self.maximum_generations, bool)
            or not isinstance(self.maximum_generations, int)
            or self.dimension < 1
            or self.multiplier < 1
            or self.size != max(5, self.dimension * self.multiplier)
            or self.maximum_generations < 1
        ):
            raise DirectTrfConstructionError("Invalid DE population topology")
        if (
            self.strategy != "best1bin"
            or self.initialization != "latinhypercube"
            or self.updating != "deferred"
        ):
            raise DirectTrfConstructionError("Unsupported DE backend topology")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-population",
                (
                    self.dimension,
                    self.multiplier,
                    self.size,
                    self.maximum_generations,
                    self.strategy,
                    self.initialization,
                    self.updating,
                ),
            ),
        )


def _validate_invocation_search_contract(
    root_problem_identity: str,
    root_problem: OptimizationProblem,
    coordinates: tuple[DeSearchCoordinate, ...],
    search_problem: OptimizationProblem,
) -> None:
    derivation = search_problem.derivation
    if (
        not root_problem_identity
        or not root_problem.acceptance_authority
        or root_problem.identity != root_problem_identity
        or not isinstance(derivation, DeSearchProblemDerivation)
        or derivation.root_problem_identity != root_problem_identity
    ):
        raise DirectTrfConstructionError(
            "DE search problem must retain its exact root derivation"
        )
    if not isinstance(coordinates, tuple) or not coordinates:
        raise DirectTrfConstructionError(
            "DE requires immutable explicitly selected coordinates"
        )
    if any(not isinstance(item, DeSearchCoordinate) for item in coordinates):
        raise DirectTrfConstructionError("DE search coordinates are malformed")
    selected_ids = tuple(item.param_id for item in coordinates)
    known_ids = {param_id for param_id, _value in search_problem.independent_items}
    if (
        selected_ids != search_problem.controlled_ids
        or selected_ids
        != tuple(sorted(selected_ids, key=lambda item: item.encode("utf-8")))
        or len(set(selected_ids)) != len(selected_ids)
        or set(selected_ids) - known_ids
    ):
        raise DirectTrfConstructionError(
            "DE coordinates must exactly match canonical selected ownership"
        )
    if {item.declaration_ordinal for item in coordinates} != set(
        range(len(coordinates))
    ):
        raise DirectTrfConstructionError(
            "DE coordinate declaration ordinals must be unique and complete"
        )
    specification_identity = _identity(
        "native-de-search-specification",
        tuple(item.identity for item in coordinates),
    )
    if derivation.search_specification_identity != specification_identity:
        raise DirectTrfConstructionError(
            "DE coordinate transforms or bounds differ from their parameter IDs"
        )
    root_indices = {
        param_id: index for index, param_id in enumerate(root_problem.controlled_ids)
    }
    if set(selected_ids) - set(root_indices):
        raise DirectTrfConstructionError(
            "DE search coordinates are not eligible root-controlled coordinates"
        )
    expected = _expected_search_projection(
        root_problem,
        selected_ids,
        specification_identity,
        root_indices,
    )
    _validate_exact_search_projection(
        root_problem,
        search_problem,
        selected_ids,
        expected,
    )
    for index, coordinate in enumerate(coordinates):
        if (
            coordinate.physical_lower < search_problem.lower_bounds[index]
            or coordinate.physical_upper > search_problem.upper_bounds[index]
        ):
            raise DirectTrfConstructionError(
                f"DE search range for {coordinate.param_id!r} exceeds physical bounds"
            )


def _expected_search_projection(
    root_problem: OptimizationProblem,
    selected_ids: tuple[str, ...],
    specification_identity: str,
    root_indices: dict[str, int],
) -> _ExpectedSearchProjection:
    selected = set(selected_ids)
    held_items = tuple(
        item for item in root_problem.independent_items if item[0] not in selected
    )
    root_starts = dict(root_problem.independent_items)
    start = tuple(root_starts[param_id] for param_id in selected_ids)
    return _ExpectedSearchProjection(
        held_items,
        start,
        tuple(
            root_problem.lower_bounds[root_indices[param_id]]
            for param_id in selected_ids
        ),
        tuple(
            root_problem.upper_bounds[root_indices[param_id]]
            for param_id in selected_ids
        ),
        DeSearchProblemDerivation(
            root_problem.identity,
            root_problem.affine_feasibility_identity,
            specification_identity,
            selected_ids,
            held_items,
            start,
        ),
    )


def _validate_exact_search_projection(
    root_problem: OptimizationProblem,
    search_problem: OptimizationProblem,
    selected_ids: tuple[str, ...],
    expected: _ExpectedSearchProjection,
) -> None:
    derivation = search_problem.derivation
    try:
        root_problem.validate_derived_problem(search_problem)
    except DirectTrfConstructionError as error:
        raise DirectTrfConstructionError(
            "DE search problem is not the exact canonical root projection"
        ) from error
    if (
        search_problem.controlled_ids != selected_ids
        or search_problem.start != expected.start
        or derivation != expected.derivation
    ):
        raise DirectTrfConstructionError(
            "DE search problem is not the exact canonical root projection"
        )


@dataclass(frozen=True, slots=True)
class DeSearchInvocation:
    """Product DE search policy with no local-fit, acceptance, or commit state."""

    root_problem_identity: str
    root_problem: OptimizationProblem = field(repr=False, compare=False)
    search_coordinates: tuple[DeSearchCoordinate, ...]
    search_problem: OptimizationProblem = field(repr=False, compare=False)
    root_seed: int
    population: DePopulation
    de_objective_request_budget: int
    mutation: float | tuple[float, float]
    recombination: float
    tol: float
    atol: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_invocation_search_contract(
            self.root_problem_identity,
            self.root_problem,
            self.search_coordinates,
            self.search_problem,
        )
        if (
            isinstance(self.root_seed, bool)
            or not isinstance(self.root_seed, int)
            or not 0 <= self.root_seed <= _UINT64_MAX
        ):
            raise DirectTrfConstructionError(
                "DE requires an explicit unsigned 64-bit root seed"
            )
        if not isinstance(self.population, DePopulation) or (
            self.population.dimension != len(self.search_coordinates)
        ):
            raise DirectTrfConstructionError(
                "DE population topology does not match selected coordinates"
            )
        de_budget = _positive_integer_budget(
            self.de_objective_request_budget,
            name="DE objective-request budget",
        )
        if de_budget < self.population.size:
            raise DirectTrfConstructionError(
                "DE objective-request budget must cover the initial population"
            )
        mutation = _mutation(self.mutation)
        recombination = _probability(
            self.recombination,
            name="DE recombination probability",
        )
        tol = _non_negative(self.tol, name="DE relative tolerance")
        atol = _non_negative(self.atol, name="DE absolute tolerance")
        object.__setattr__(self, "de_objective_request_budget", de_budget)
        object.__setattr__(self, "mutation", mutation)
        object.__setattr__(self, "recombination", recombination)
        object.__setattr__(self, "tol", tol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-search-invocation",
                (
                    _DE_BACKEND_POLICY_VERSION,
                    self.root_problem_identity,
                    tuple(item.identity for item in self.search_coordinates),
                    self.search_problem.identity,
                    self.root_seed,
                    self.population.identity,
                    de_budget,
                    (
                        tuple(_float_token(value) for value in mutation)
                        if isinstance(mutation, tuple)
                        else _float_token(mutation)
                    ),
                    _float_token(recombination),
                    _float_token(tol),
                    _float_token(atol),
                ),
            ),
        )

    @classmethod
    def for_product_problem(
        cls,
        problem: OptimizationProblem,
        *,
        search_coordinates: Sequence[
            tuple[str, float, float, DeCoordinateSemantics | str]
        ],
        root_seed: object,
    ) -> DeSearchInvocation:
        """Compile the versioned, non-configurable product DE backend policy."""
        dimension = len(search_coordinates)
        population_size = max(
            5,
            dimension * _PRODUCT_DE_POPULATION_MULTIPLIER,
        )
        return cls._for_problem_with_policy(
            problem,
            search_coordinates=search_coordinates,
            root_seed=root_seed,
            de_objective_request_budget=(
                population_size * (_PRODUCT_DE_MAXIMUM_GENERATIONS + 1)
            ),
            population_multiplier=_PRODUCT_DE_POPULATION_MULTIPLIER,
            maximum_generations=_PRODUCT_DE_MAXIMUM_GENERATIONS,
            mutation=(0.5, 1.0),
            recombination=0.7,
            tol=0.01,
            atol=0.0,
        )

    @classmethod
    def _for_problem_with_policy(  # noqa: C901 - closed DE construction policy
        cls,
        problem: OptimizationProblem,
        *,
        search_coordinates: Sequence[
            tuple[str, float, float, DeCoordinateSemantics | str]
        ],
        root_seed: object,
        de_objective_request_budget: int,
        population_multiplier: int,
        maximum_generations: int,
        mutation: float | tuple[float, float],
        recombination: float,
        tol: float,
        atol: float,
    ) -> DeSearchInvocation:
        if not problem.acceptance_authority:
            raise DirectTrfConstructionError("DE requires a complete root problem")
        if (
            isinstance(root_seed, bool)
            or not isinstance(root_seed, int)
            or not 0 <= root_seed <= _UINT64_MAX
        ):
            raise DirectTrfConstructionError(
                "DE requires an explicit unsigned 64-bit root seed"
            )
        if (
            isinstance(population_multiplier, bool)
            or not isinstance(population_multiplier, int)
            or population_multiplier < 1
            or isinstance(maximum_generations, bool)
            or not isinstance(maximum_generations, int)
            or maximum_generations < 1
        ):
            raise DirectTrfConstructionError(
                "DE population multiplier and generation limit must be positive integers"
            )
        declared = tuple(
            DeSearchCoordinate(param_id, lower, upper, semantics, ordinal)
            for ordinal, (param_id, lower, upper, semantics) in enumerate(
                search_coordinates
            )
        )
        if not declared:
            raise DirectTrfConstructionError(
                "DE requires at least one explicitly selected coordinate"
            )
        selected_ids = tuple(item.param_id for item in declared)
        if len(set(selected_ids)) != len(selected_ids):
            raise DirectTrfConstructionError("Duplicate DE search coordinate ID")
        unknown = set(selected_ids) - set(problem.controlled_ids)
        if unknown:
            raise DirectTrfConstructionError(
                "DE search coordinates must target controlled independent IDs: "
                + ", ".join(sorted(unknown))
            )
        coordinates = tuple(
            sorted(declared, key=lambda item: item.param_id.encode("utf-8"))
        )
        root_indices = {
            param_id: index for index, param_id in enumerate(problem.controlled_ids)
        }
        root_starts = dict(zip(problem.controlled_ids, problem.start, strict=True))
        for item in coordinates:
            index = root_indices[item.param_id]
            if (
                item.physical_lower < problem.lower_bounds[index]
                or item.physical_upper > problem.upper_bounds[index]
            ):
                raise DirectTrfConstructionError(
                    f"DE search range for {item.param_id!r} exceeds physical bounds"
                )
        canonical_selected_ids = tuple(item.param_id for item in coordinates)
        selected = set(canonical_selected_ids)
        held_items = tuple(
            (param_id, value)
            for param_id, value in problem.independent_items
            if param_id not in selected
        )
        start = tuple(root_starts[param_id] for param_id in canonical_selected_ids)

        def derive_search_problem() -> OptimizationProblem:
            search_specification_identity = _identity(
                "native-de-search-specification",
                tuple(item.identity for item in coordinates),
            )
            derivation = DeSearchProblemDerivation(
                problem.identity,
                problem.affine_feasibility_identity,
                search_specification_identity,
                canonical_selected_ids,
                held_items,
                start,
            )
            return problem.derive_child(
                controlled_ids=canonical_selected_ids,
                start=start,
                derivation=derivation,
            )

        search_problem = derive_search_problem()
        feasible = search_problem.feasible_coordinates
        if feasible is not None and not feasible.supports_box_only_algorithms:
            raise DirectTrfConstructionError(
                "DE search over relaxation-feasibility coordinates is not yet "
                "qualified; use Direct TRF for this method section"
            )
        if feasible is not None:
            feasible_lower, feasible_upper = feasible.solver_bounds
            child_indices = {
                param_id: index
                for index, param_id in enumerate(search_problem.controlled_ids)
            }
            coordinates = tuple(
                DeSearchCoordinate(
                    item.param_id,
                    max(
                        item.physical_lower,
                        feasible_lower[child_indices[item.param_id]],
                    ),
                    min(
                        item.physical_upper,
                        feasible_upper[child_indices[item.param_id]],
                    ),
                    item.semantics,
                    item.declaration_ordinal,
                )
                for item in coordinates
            )
            search_problem = derive_search_problem()
        population = DePopulation(
            len(coordinates),
            population_multiplier,
            max(5, len(coordinates) * population_multiplier),
            maximum_generations,
        )
        return cls(
            problem.identity,
            problem,
            coordinates,
            search_problem,
            cast("int", root_seed),
            population,
            de_objective_request_budget,
            mutation,
            recombination,
            tol,
            atol,
        )


def _positive_integer_budget(value: int, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise DirectTrfConstructionError(f"{name} must be a positive integer")
    return value


def _mutation(value: float | tuple[float, float]) -> float | tuple[float, float]:
    if isinstance(value, tuple):
        if len(value) != 2:
            raise DirectTrfConstructionError(
                "DE mutation dithering requires exactly two values"
            )
        lower = _finite_binary64(value[0], name="DE mutation lower")
        upper = _finite_binary64(value[1], name="DE mutation upper")
        if not 0.0 <= lower < upper < 2.0:
            raise DirectTrfConstructionError(
                "DE mutation dithering must satisfy 0 <= lower < upper < 2"
            )
        return lower, upper
    scalar = _finite_binary64(value, name="DE mutation")
    if not 0.0 <= scalar < 2.0:
        raise DirectTrfConstructionError("DE mutation must satisfy 0 <= value < 2")
    return scalar


def _probability(value: float, *, name: str) -> float:
    result = _finite_binary64(value, name=name)
    if not 0.0 <= result <= 1.0:
        raise DirectTrfConstructionError(f"{name} must lie in [0, 1]")
    return result


def _non_negative(value: float, *, name: str) -> float:
    result = _finite_binary64(value, name=name)
    if result < 0.0:
        raise DirectTrfConstructionError(f"{name} must be non-negative")
    return result


_DE_RESTART_TERMINALS = frozenset(
    {
        "population_converged",
        "generation_limit",
        "budget_exhausted",
    }
)
_DE_CANDIDATE_ORDER_VERSION = "chi-square-root-vector-request-ordinal-v1"


@dataclass(frozen=True, slots=True)
class DeSearchCandidate:
    """One valid DE evaluation under the root-coordinate total order."""

    solver_vector: tuple[float, ...]
    selected_vector: tuple[float, ...]
    full_vector: tuple[float, ...]
    chi_square: float
    request_ordinal: int
    evaluation_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        solver = tuple(
            _finite_binary64(value, name=f"DE solver vector[{index}]")
            for index, value in enumerate(self.solver_vector)
        )
        selected = tuple(
            _finite_binary64(value, name=f"DE selected vector[{index}]")
            for index, value in enumerate(self.selected_vector)
        )
        full = tuple(
            _finite_binary64(value, name=f"DE full vector[{index}]")
            for index, value in enumerate(self.full_vector)
        )
        objective = _finite_binary64(self.chi_square, name="DE candidate chi-square")
        if (
            not solver
            or len(solver) != len(selected)
            or not full
            or objective < 0.0
            or isinstance(self.request_ordinal, bool)
            or self.request_ordinal < 1
            or not self.evaluation_identity
        ):
            raise ValueError("Invalid DE search candidate evidence")
        object.__setattr__(self, "solver_vector", solver)
        object.__setattr__(self, "selected_vector", selected)
        object.__setattr__(self, "full_vector", full)
        object.__setattr__(self, "chi_square", objective)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-search-candidate",
                (
                    _DE_CANDIDATE_ORDER_VERSION,
                    tuple(_float_token(value) for value in solver),
                    tuple(_float_token(value) for value in selected),
                    tuple(_float_token(value) for value in full),
                    _float_token(objective),
                    self.request_ordinal,
                    self.evaluation_identity,
                ),
            ),
        )

    def ordering_key(self) -> tuple[float, tuple[float, ...], int]:
        return self.chi_square, self.full_vector, self.request_ordinal


class DeSearchTerminal(StrEnum):
    """Closed terminal outcomes for the selected-coordinate DE stage."""

    POPULATION_CONVERGED = "population_converged"
    GENERATION_LIMIT = "generation_limit"
    BUDGET_EXHAUSTED = "budget_exhausted"
    NO_VALID_CANDIDATE = "no_valid_candidate"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    BACKEND_FAILURE = "backend_failure"
    IMPLEMENTATION_FAILURE = "implementation_failure"


@dataclass(frozen=True, slots=True)
class DeBackendEvidence:
    """Detached closed evidence normalized from SciPy's DE result."""

    success: bool
    message: str
    generation_count: int
    backend_evaluation_count: int
    population_size: int
    finite_population_energies: int
    returned_solver_vector: tuple[float, ...]
    returned_objective: float | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-backend-evidence",
                (
                    self.success,
                    self.message,
                    self.generation_count,
                    self.backend_evaluation_count,
                    self.population_size,
                    self.finite_population_energies,
                    tuple(_float_token(value) for value in self.returned_solver_vector),
                    None
                    if self.returned_objective is None
                    else _float_token(self.returned_objective),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class DeSearchExecution:
    """Immutable accounting and candidate evidence for one DE attempt."""

    occurrence_identity: str = field(compare=False)
    problem_identity: str
    workflow_invocation_identity: str
    terminal: DeSearchTerminal
    counters: AttemptCounters
    population_identity: str
    candidate_ordering_policy: str
    valid_candidate_count: int
    rejected_trial_count: int
    best_candidate: DeSearchCandidate | None
    backend: DeBackendEvidence | None
    failure: TerminalFailure | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.valid_candidate_count < 0 or self.rejected_trial_count < 0:
            raise ValueError("DE candidate counts cannot be negative")
        if (self.best_candidate is None) != (self.valid_candidate_count == 0):
            raise ValueError("DE best-candidate evidence disagrees with its count")
        if self.restart_eligible and self.best_candidate is None:
            raise ValueError("Restart-eligible DE outcome lacks a candidate")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-search-execution",
                (
                    self.problem_identity,
                    self.workflow_invocation_identity,
                    self.terminal.value,
                    (
                        self.counters.solver_requests_received,
                        self.counters.objective_requests_accepted,
                        self.counters.objective_evaluations_completed,
                    ),
                    self.population_identity,
                    self.candidate_ordering_policy,
                    self.valid_candidate_count,
                    self.rejected_trial_count,
                    None
                    if self.best_candidate is None
                    else self.best_candidate.identity,
                    None if self.backend is None else self.backend.identity,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )

    @property
    def restart_eligible(self) -> bool:
        return (
            self.terminal.value in _DE_RESTART_TERMINALS
            and self.best_candidate is not None
        )


class _DeAttemptStop(RuntimeError):
    def __init__(self, terminal: DeSearchTerminal, failure: TerminalFailure) -> None:
        self.terminal = terminal
        self.failure = failure
        super().__init__(failure.message or failure.category)


class _LiveDeAttempt:
    def __init__(
        self,
        root_problem: OptimizationProblem,
        invocation: DeSearchInvocation,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
        cancellation: CancellationToken,
    ) -> None:
        self.root_problem = root_problem
        self.invocation = invocation
        self.parameterization = parameterization
        self.evaluator = engine.new_evaluator()
        self.cancellation = cancellation
        self.received = 0
        self.accepted = 0
        self.completed = 0
        self.rejected = 0
        self.candidates: list[DeSearchCandidate] = []
        self.best: DeSearchCandidate | None = None

    @property
    def counters(self) -> AttemptCounters:
        return AttemptCounters(self.received, self.accepted, self.completed)

    def map_solver_vector(
        self, solver_vector: Sequence[float] | Array
    ) -> tuple[float, ...]:
        values = tuple(float(value) for value in solver_vector)
        if len(values) != len(self.invocation.search_coordinates):
            raise ValueError("DE solver vector has the wrong dimension")
        return tuple(
            coordinate.to_physical(value)
            for coordinate, value in zip(
                self.invocation.search_coordinates,
                values,
                strict=True,
            )
        )

    def objective(self, solver_vector: Array) -> float:
        self.received += 1
        if self.cancellation.is_cancelled:
            raise _DeAttemptStop(
                DeSearchTerminal.CANCELLED,
                TerminalFailure("cancelled", "Cancellation observed before DE request"),
            )
        if self.accepted >= self.invocation.de_objective_request_budget:
            raise _DeAttemptStop(
                DeSearchTerminal.BUDGET_EXHAUSTED,
                TerminalFailure(
                    "de_objective_budget_exhausted",
                    "DE objective-request budget exhausted",
                ),
            )
        self.accepted += 1
        try:
            selected = self.map_solver_vector(solver_vector)
            lifecycle = self.invocation.search_problem.lifecycle_frame(
                selected,
                self.parameterization,
            )
            frame = EvaluationFrame.from_lifecycle_frame(
                self.parameterization,
                lifecycle,
            )
        except Exception as error:
            raise _DeAttemptStop(
                DeSearchTerminal.IMPLEMENTATION_FAILURE,
                TerminalFailure(
                    "de_candidate_contract_failure",
                    f"{type(error).__name__}: {error}",
                ),
            ) from error
        evaluated = self.evaluator.evaluate(frame)
        self.completed += 1
        if isinstance(evaluated, EvaluationFailure):
            if evaluated.validity == "INVALID_TRIAL":
                self.rejected += 1
                return math.inf
            raise _DeAttemptStop(
                DeSearchTerminal.IMPLEMENTATION_FAILURE,
                TerminalFailure(evaluated.category, evaluated.message, evaluated),
            )
        try:
            objective = canonical_chi_square(evaluated.residuals)
        except (TypeError, ValueError, ObjectiveScalarizationError):
            self.rejected += 1
            return math.inf
        try:
            full_vector = tuple(
                evaluated.resolved_values[param_id]
                for param_id in self.root_problem.controlled_ids
            )
            candidate = DeSearchCandidate(
                tuple(float(value) for value in solver_vector),
                selected,
                full_vector,
                objective,
                self.received,
                evaluated.identity,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise _DeAttemptStop(
                DeSearchTerminal.IMPLEMENTATION_FAILURE,
                TerminalFailure(
                    "de_candidate_evidence_failure",
                    f"{type(error).__name__}: {error}",
                ),
            ) from error
        self.candidates.append(candidate)
        if self.best is None or candidate.ordering_key() < self.best.ordering_key():
            self.best = candidate
        if self.cancellation.is_cancelled:
            raise _DeAttemptStop(
                DeSearchTerminal.CANCELLED,
                TerminalFailure("cancelled", "Cancellation observed after DE request"),
            )
        return objective


def _validate_de_context(
    problem: OptimizationProblem,
    invocation: DeSearchInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> None:
    if (
        not problem.acceptance_authority
        or invocation.root_problem_identity != problem.identity
        or engine.plan.identity != problem.evaluation_plan_identity
        or engine.plan.identity != invocation.search_problem.evaluation_plan_identity
    ):
        raise DirectTrfConstructionError(
            "DE invocation, root problem, and evaluator context differ"
        )
    problem.validate_parameterization(parameterization)
    invocation.search_problem.validate_parameterization(parameterization)


def _normalize_de_backend(
    result: object,
    *,
    dimension: int,
    population_size: int,
) -> DeBackendEvidence:
    backend = cast("_DeBackendResult", result)
    success = backend.success
    if not isinstance(success, (bool, np.bool_)):
        raise TypeError("DE backend success must be Boolean")
    message = backend.message
    if not isinstance(message, str):
        raise TypeError("DE backend message must be text")
    generation_count = backend.nit
    backend_count = backend.nfev
    if (
        isinstance(generation_count, bool)
        or not isinstance(generation_count, (int, np.integer))
        or int(generation_count) < 0
        or isinstance(backend_count, bool)
        or not isinstance(backend_count, (int, np.integer))
        or int(backend_count) < 0
    ):
        raise ValueError("DE backend counters must be non-negative integers")
    returned = np.asarray(backend.x, dtype=np.float64)
    population = np.asarray(backend.population, dtype=np.float64)
    energies = np.asarray(backend.population_energies, dtype=np.float64)
    if (
        returned.shape != (dimension,)
        or not np.all(np.isfinite(returned))
        or population.shape != (population_size, dimension)
        or energies.shape != (population_size,)
    ):
        raise ValueError("DE backend returned malformed population evidence")
    raw_objective = float(cast("Real", backend.fun))
    objective = raw_objective if math.isfinite(raw_objective) else None
    return DeBackendEvidence(
        bool(success),
        message,
        int(generation_count),
        int(backend_count),
        population_size,
        int(np.count_nonzero(np.isfinite(energies))),
        tuple(float(value) for value in returned),
        objective,
    )


class _DeBackendResult(Protocol):
    success: object
    message: object
    nit: object
    nfev: object
    x: object
    population: object
    population_energies: object
    fun: object


def _search_execution(
    live: _LiveDeAttempt,
    terminal: DeSearchTerminal,
    *,
    backend: DeBackendEvidence | None = None,
    failure: TerminalFailure | None = None,
) -> DeSearchExecution:
    effective_terminal = terminal
    effective_failure = failure
    if (
        terminal
        in {
            DeSearchTerminal.POPULATION_CONVERGED,
            DeSearchTerminal.GENERATION_LIMIT,
            DeSearchTerminal.BUDGET_EXHAUSTED,
        }
        and live.best is None
    ):
        effective_terminal = DeSearchTerminal.NO_VALID_CANDIDATE
        effective_failure = TerminalFailure(
            "de_no_valid_candidate",
            "DE completed without a valid finite objective evaluation",
        )
    return DeSearchExecution(
        uuid4().hex,
        live.invocation.search_problem.identity,
        live.invocation.identity,
        effective_terminal,
        live.counters,
        live.invocation.population.identity,
        _DE_CANDIDATE_ORDER_VERSION,
        len(live.candidates),
        live.rejected,
        live.best,
        backend,
        effective_failure,
    )


def _invoke_de_backend(
    live: _LiveDeAttempt,
    invocation: DeSearchInvocation,
    solver_start: Array,
) -> object:
    return differential_evolution(
        live.objective,
        Bounds(
            np.asarray(
                tuple(item.solver_bounds[0] for item in invocation.search_coordinates),
                dtype=np.float64,
            ),
            np.asarray(
                tuple(item.solver_bounds[1] for item in invocation.search_coordinates),
                dtype=np.float64,
            ),
        ),
        strategy="best1bin",
        maxiter=invocation.population.maximum_generations,
        popsize=invocation.population.multiplier,
        tol=invocation.tol,
        mutation=invocation.mutation,
        recombination=invocation.recombination,
        rng=np.random.Generator(np.random.PCG64(invocation.root_seed)),
        callback=None,
        disp=False,
        polish=False,
        init="latinhypercube",
        atol=invocation.atol,
        updating="deferred",
        workers=1,
        constraints=(),
        x0=solver_start,
        integrality=None,
        vectorized=False,
    )


def _finish_de_backend(
    live: _LiveDeAttempt,
    result: object,
) -> DeSearchExecution:
    try:
        backend = _normalize_de_backend(
            result,
            dimension=len(live.invocation.search_coordinates),
            population_size=live.invocation.population.size,
        )
    except Exception as error:  # noqa: BLE001 - mutable backend boundary
        return _search_execution(
            live,
            DeSearchTerminal.IMPLEMENTATION_FAILURE,
            failure=TerminalFailure(
                "malformed_de_backend_result",
                f"{type(error).__name__}: {error}",
            ),
        )
    if backend.success:
        terminal = DeSearchTerminal.POPULATION_CONVERGED
        failure = None
    elif "maximum number of iterations" in backend.message.lower():
        terminal = DeSearchTerminal.GENERATION_LIMIT
        failure = TerminalFailure("de_generation_limit", backend.message)
    else:
        terminal = DeSearchTerminal.BACKEND_FAILURE
        failure = TerminalFailure("de_backend_failure", backend.message)
    return _search_execution(
        live,
        terminal,
        backend=backend,
        failure=failure,
    )


def _unstarted_search_execution(
    invocation: DeSearchInvocation,
    terminal: DeSearchTerminal,
    failure: TerminalFailure,
) -> DeSearchExecution:
    """Freeze typed zero-counter evidence when DE cannot acquire an evaluator."""
    return DeSearchExecution(
        uuid4().hex,
        invocation.search_problem.identity,
        invocation.identity,
        terminal,
        AttemptCounters(0, 0, 0),
        invocation.population.identity,
        _DE_CANDIDATE_ORDER_VERSION,
        0,
        0,
        None,
        None,
        failure,
    )


def _start_live_de_attempt(
    problem: OptimizationProblem,
    invocation: DeSearchInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    cancellation: CancellationToken,
) -> _LiveDeAttempt | DeSearchExecution:
    if cancellation.is_cancelled:
        return _unstarted_search_execution(
            invocation,
            DeSearchTerminal.CANCELLED,
            TerminalFailure("cancelled", "Cancellation before DE evaluator binding"),
        )
    try:
        live = _LiveDeAttempt(
            problem,
            invocation,
            parameterization,
            engine,
            cancellation,
        )
    except KeyboardInterrupt:
        return _unstarted_search_execution(
            invocation,
            DeSearchTerminal.INTERRUPTED,
            TerminalFailure(
                "interrupted",
                "KeyboardInterrupt during DE evaluator binding",
            ),
        )
    except Exception as error:  # noqa: BLE001 - binding failure is typed evidence
        return _unstarted_search_execution(
            invocation,
            DeSearchTerminal.IMPLEMENTATION_FAILURE,
            TerminalFailure(
                "de_evaluator_binding_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    return live


def _run_de_search(
    problem: OptimizationProblem,
    invocation: DeSearchInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    cancellation: CancellationToken,
) -> DeSearchExecution:
    started = _start_live_de_attempt(
        problem,
        invocation,
        parameterization,
        engine,
        cancellation,
    )
    if isinstance(started, DeSearchExecution):
        return started
    live = started
    if cancellation.is_cancelled:
        return _search_execution(
            live,
            DeSearchTerminal.CANCELLED,
            failure=TerminalFailure(
                "cancelled", "Cancellation after evaluator binding"
            ),
        )
    solver_start = np.asarray(
        tuple(
            coordinate.solver_initial for coordinate in invocation.search_coordinates
        ),
        dtype=np.float64,
    )
    try:
        result = _invoke_de_backend(live, invocation, solver_start)
    except _DeAttemptStop as stop:
        return _search_execution(
            live,
            stop.terminal,
            failure=stop.failure,
        )
    except KeyboardInterrupt:
        return _search_execution(
            live,
            DeSearchTerminal.INTERRUPTED,
            failure=TerminalFailure("interrupted", "KeyboardInterrupt during DE"),
        )
    except Exception as error:  # noqa: BLE001 - undeclared backend errors fail closed
        return _search_execution(
            live,
            DeSearchTerminal.IMPLEMENTATION_FAILURE,
            failure=TerminalFailure(
                "unexpected_de_backend_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    if cancellation.is_cancelled:
        return _search_execution(
            live,
            DeSearchTerminal.CANCELLED,
            failure=TerminalFailure("cancelled", "Cancellation after DE backend"),
        )
    return _finish_de_backend(live, result)


def execute_de_search(
    problem: OptimizationProblem,
    invocation: DeSearchInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> DeSearchExecution:
    """Run only selected-coordinate DE, without fit acceptance or commit authority."""
    _validate_de_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    return _run_de_search(problem, invocation, parameterization, engine, token)
