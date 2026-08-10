"""Selected-coordinate differential evolution to native TRF qualification (#597).

This isolated native seam is intentionally not wired into production dispatch.
It owns an immutable selected-coordinate DE plan rooted in one complete native
problem.  Execution and the authoritative TRF transition are added in vertical
slices below the same seam.
"""

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
    EvaluationResult,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    AttemptCounters,
    CancellationToken,
    CandidateMaterialization,
    CommitReceipt,
    DePolishProblemDerivation,
    DeSearchProblemDerivation,
    DirectTrfCandidateOutcome,
    DirectTrfCandidateTerminal,
    DirectTrfConstructionError,
    DirectTrfInterrupted,
    DirectTrfInvocation,
    LiveFitCommitAuthority,
    MaterializationTerminal,
    MaterializedDirectTrfCandidate,
    ObjectiveScalarizationError,
    OptimizationProblem,
    TerminalFailure,
    _accept_materialized_fit_for_derived_workflow,
    _grant_derived_fit_commit_authority,
    _materialize_derived_direct_trf_candidate_for_root,
    canonical_chi_square,
    commit_accepted_fit,
    execute_direct_trf_candidate,
)
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.values import AnalysisValues
from chemex.typing import Array

_SCHEMA_VERSION = 1
_DE_WORKFLOW_VERSION = "native-selected-coordinate-de-direct-trf-v1"
_DE_BACKEND_POLICY_VERSION = "scipy-de-best1bin-lhs-deferred-pcg64-v1"
_UINT64_MAX = 2**64 - 1


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

    def to_solver(self, physical_value: float) -> float:
        value = _finite_binary64(physical_value, name=f"DE start {self.param_id!r}")
        if not self.physical_lower <= value <= self.physical_upper:
            raise DirectTrfConstructionError(
                f"DE start for {self.param_id!r} is outside its search range"
            )
        if self.semantics is DeCoordinateSemantics.LOG:
            if value <= 0.0:
                raise DirectTrfConstructionError(
                    f"Logarithmic DE start for {self.param_id!r} must be positive"
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
        coordinate.to_solver(search_problem.start[index])


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
    if (
        search_problem.evaluation_plan_identity != root_problem.evaluation_plan_identity
        or search_problem.parameterization_identity
        != root_problem.parameterization_identity
        or search_problem.evaluator_parameterization_identity
        != root_problem.evaluator_parameterization_identity
        or search_problem.constraint_program_identity
        != root_problem.constraint_program_identity
        or search_problem.configuration_identity != root_problem.configuration_identity
        or search_problem.source_snapshot != root_problem.source_snapshot
        or search_problem.independent_items != root_problem.independent_items
        or search_problem.controlled_ids != selected_ids
        or search_problem.held_items != expected.held_items
        or search_problem.start != expected.start
        or search_problem.lower_bounds != expected.lower_bounds
        or search_problem.upper_bounds != expected.upper_bounds
        or search_problem.commit_scope != root_problem.commit_scope
        or search_problem.scalarization_version != root_problem.scalarization_version
        or derivation != expected.derivation
    ):
        raise DirectTrfConstructionError(
            "DE search problem is not the exact canonical root projection"
        )


@dataclass(frozen=True, slots=True)
class DeDirectTrfInvocation:
    """Immutable selected-coordinate DE search and one full TRF polish policy."""

    root_problem_identity: str
    root_problem: OptimizationProblem = field(repr=False, compare=False)
    search_coordinates: tuple[DeSearchCoordinate, ...]
    search_problem: OptimizationProblem = field(repr=False, compare=False)
    root_seed: int
    population: DePopulation
    de_objective_request_budget: int
    polish_objective_request_budget: int
    mutation: float | tuple[float, float]
    recombination: float
    tol: float
    atol: float
    polish_x_scale: tuple[float, ...]
    polish_ftol: float | None
    polish_xtol: float | None
    polish_gtol: float | None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        coordinates = self.search_coordinates
        _validate_invocation_search_contract(
            self.root_problem_identity,
            self.root_problem,
            coordinates,
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
            self.population.dimension != len(coordinates)
        ):
            raise DirectTrfConstructionError(
                "DE population topology does not match selected coordinates"
            )
        de_budget = _positive_integer_budget(
            self.de_objective_request_budget,
            name="DE objective-request budget",
        )
        polish_budget = _positive_integer_budget(
            self.polish_objective_request_budget,
            name="TRF polish objective-request budget",
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
        polish = DirectTrfInvocation(
            self.root_problem_identity,
            polish_budget,
            self.polish_x_scale,
            self.polish_ftol,
            self.polish_xtol,
            self.polish_gtol,
        )
        if len(polish.x_scale) != len(self.root_problem.controlled_ids):
            raise DirectTrfConstructionError(
                "TRF polish scale must match the complete root coordinate dimension"
            )
        object.__setattr__(self, "de_objective_request_budget", de_budget)
        object.__setattr__(self, "polish_objective_request_budget", polish_budget)
        object.__setattr__(self, "mutation", mutation)
        object.__setattr__(self, "recombination", recombination)
        object.__setattr__(self, "tol", tol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(self, "polish_x_scale", polish.x_scale)
        object.__setattr__(self, "polish_ftol", polish.ftol)
        object.__setattr__(self, "polish_xtol", polish.xtol)
        object.__setattr__(self, "polish_gtol", polish.gtol)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-direct-trf-invocation",
                (
                    _DE_WORKFLOW_VERSION,
                    _DE_BACKEND_POLICY_VERSION,
                    self.root_problem_identity,
                    tuple(item.identity for item in self.search_coordinates),
                    self.search_problem.identity,
                    self.root_seed,
                    self.population.identity,
                    de_budget,
                    polish_budget,
                    (
                        tuple(_float_token(value) for value in mutation)
                        if isinstance(mutation, tuple)
                        else _float_token(mutation)
                    ),
                    _float_token(recombination),
                    _float_token(tol),
                    _float_token(atol),
                    tuple(_float_token(value) for value in polish.x_scale),
                    tuple(
                        None if value is None else _float_token(value)
                        for value in (
                            polish.ftol,
                            polish.xtol,
                            polish.gtol,
                        )
                    ),
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        """Return a JSON-safe declarative record with no live runtime state."""
        mutation: float | list[float]
        if isinstance(self.mutation, tuple):
            mutation = list(self.mutation)
        else:
            mutation = self.mutation
        return {
            "schema_version": _SCHEMA_VERSION,
            "workflow_version": _DE_WORKFLOW_VERSION,
            "backend_policy_version": _DE_BACKEND_POLICY_VERSION,
            "identity": self.identity,
            "root_problem_identity": self.root_problem_identity,
            "search_problem_identity": self.search_problem.identity,
            "search_problem": {
                "controlled_ids": list(self.search_problem.controlled_ids),
                "captured_held_items": [
                    [param_id, value]
                    for param_id, value in self.search_problem.held_items
                ],
                "physical_start": list(self.search_problem.start),
                "physical_lower_bound_tokens": [
                    _float_token(value) for value in self.search_problem.lower_bounds
                ],
                "physical_upper_bound_tokens": [
                    _float_token(value) for value in self.search_problem.upper_bounds
                ],
                "solver_x0": [
                    coordinate.to_solver(value)
                    for coordinate, value in zip(
                        self.search_coordinates,
                        self.search_problem.start,
                        strict=True,
                    )
                ],
            },
            "search_coordinates": [
                {
                    "identity": item.identity,
                    "param_id": item.param_id,
                    "physical_lower": item.physical_lower,
                    "physical_upper": item.physical_upper,
                    "semantics": cast("DeCoordinateSemantics", item.semantics).value,
                    "declaration_ordinal": item.declaration_ordinal,
                }
                for item in self.search_coordinates
            ],
            "root_seed": self.root_seed,
            "population": {
                "identity": self.population.identity,
                "dimension": self.population.dimension,
                "multiplier": self.population.multiplier,
                "size": self.population.size,
                "maximum_generations": self.population.maximum_generations,
                "strategy": self.population.strategy,
                "initialization": self.population.initialization,
                "updating": self.population.updating,
            },
            "budgets": {
                "de_objective_requests": self.de_objective_request_budget,
                "polish_objective_requests": self.polish_objective_request_budget,
            },
            "candidate_ordering_policy": _DE_CANDIDATE_ORDER_VERSION,
            "mutation": mutation,
            "recombination": self.recombination,
            "tol": self.tol,
            "atol": self.atol,
            "polish": {
                "x_scale": list(self.polish_x_scale),
                "ftol": self.polish_ftol,
                "xtol": self.polish_xtol,
                "gtol": self.polish_gtol,
            },
            "authority": {
                "de_acceptance": False,
                "de_commit": False,
                "eligible_transition_count": 1,
            },
        }

    @classmethod
    def for_problem(
        cls,
        problem: OptimizationProblem,
        *,
        search_coordinates: Sequence[
            tuple[str, float, float, DeCoordinateSemantics | str]
        ],
        root_seed: object,
        de_objective_request_budget: int,
        polish_objective_request_budget: int,
        population_multiplier: int = 15,
        mutation: float | tuple[float, float] = (0.5, 1.0),
        recombination: float = 0.7,
        tol: float = 0.01,
        atol: float = 0.0,
        maximum_generations: int = 1000,
        polish_x_scale: float | Sequence[float] | None = None,
        polish_ftol: float | None = 1.0e-8,
        polish_xtol: float | None = 1.0e-8,
        polish_gtol: float | None = 1.0e-8,
    ) -> DeDirectTrfInvocation:
        """Compile one canonical bounded selected-coordinate search plan."""
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
            item.to_solver(root_starts[item.param_id])
        canonical_selected_ids = tuple(item.param_id for item in coordinates)
        selected = set(canonical_selected_ids)
        held_items = tuple(
            (param_id, value)
            for param_id, value in problem.independent_items
            if param_id not in selected
        )
        start = tuple(root_starts[param_id] for param_id in canonical_selected_ids)
        search_specification_identity = _identity(
            "native-de-search-specification",
            tuple(item.identity for item in coordinates),
        )
        derivation = DeSearchProblemDerivation(
            problem.identity,
            search_specification_identity,
            canonical_selected_ids,
            held_items,
            start,
        )
        search_problem = OptimizationProblem(
            problem.evaluation_plan_identity,
            problem.parameterization_identity,
            problem.evaluator_parameterization_identity,
            problem.constraint_program_identity,
            problem.configuration_identity,
            problem.source_snapshot,
            problem.independent_items,
            canonical_selected_ids,
            held_items,
            start,
            tuple(
                problem.lower_bounds[root_indices[item]]
                for item in canonical_selected_ids
            ),
            tuple(
                problem.upper_bounds[root_indices[item]]
                for item in canonical_selected_ids
            ),
            problem.commit_scope,
            derivation,
            problem.scalarization_version,
        )
        population = DePopulation(
            len(coordinates),
            population_multiplier,
            max(5, len(coordinates) * population_multiplier),
            maximum_generations,
        )
        de_budget = _positive_integer_budget(
            de_objective_request_budget,
            name="DE objective-request budget",
        )
        polish_budget = _positive_integer_budget(
            polish_objective_request_budget,
            name="TRF polish objective-request budget",
        )
        normalized_mutation = _mutation(mutation)
        normalized_recombination = _probability(
            recombination,
            name="DE recombination probability",
        )
        normalized_tol = _non_negative(tol, name="DE relative tolerance")
        normalized_atol = _non_negative(atol, name="DE absolute tolerance")
        polish_template = DirectTrfInvocation.for_problem(
            problem,
            objective_request_budget=polish_budget,
            x_scale=polish_x_scale,
            ftol=polish_ftol,
            xtol=polish_xtol,
            gtol=polish_gtol,
        )
        return cls(
            problem.identity,
            problem,
            coordinates,
            search_problem,
            cast("int", root_seed),
            population,
            de_budget,
            polish_budget,
            normalized_mutation,
            normalized_recombination,
            normalized_tol,
            normalized_atol,
            polish_template.x_scale,
            polish_template.ftol,
            polish_template.xtol,
            polish_template.gtol,
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
    PREFLIGHT_INVALID = "preflight_invalid"
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
    preflight_evaluation_identity: str | None
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
                    self.preflight_evaluation_identity,
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
        invocation: DeDirectTrfInvocation,
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
    invocation: DeDirectTrfInvocation,
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


def _de_preflight(live: _LiveDeAttempt) -> EvaluationResult | TerminalFailure:
    try:
        lifecycle = live.invocation.search_problem.lifecycle_frame(
            live.invocation.search_problem.start,
            live.parameterization,
        )
        frame = EvaluationFrame.from_lifecycle_frame(live.parameterization, lifecycle)
        evaluated = live.evaluator.evaluate(frame)
    except Exception as error:  # noqa: BLE001 - preflight fails closed
        return TerminalFailure(
            "de_preflight_frame_failure",
            f"{type(error).__name__}: {error}",
        )
    if isinstance(evaluated, EvaluationFailure):
        return TerminalFailure(evaluated.category, evaluated.message, evaluated)
    try:
        canonical_chi_square(evaluated.residuals)
    except (TypeError, ValueError, ObjectiveScalarizationError) as error:
        return TerminalFailure(
            "de_preflight_scalarization_failure",
            f"{type(error).__name__}: {error}",
        )
    return evaluated


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
    preflight_identity: str | None,
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
        preflight_identity,
        len(live.candidates),
        live.rejected,
        live.best,
        backend,
        effective_failure,
    )


def _invoke_de_backend(
    live: _LiveDeAttempt,
    invocation: DeDirectTrfInvocation,
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
    preflight_identity: str,
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
            preflight_identity,
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
        preflight_identity,
        backend=backend,
        failure=failure,
    )


def _unstarted_search_execution(
    invocation: DeDirectTrfInvocation,
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
        None,
        0,
        0,
        None,
        None,
        failure,
    )


def _start_live_de_attempt(
    problem: OptimizationProblem,
    invocation: DeDirectTrfInvocation,
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
    invocation: DeDirectTrfInvocation,
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
    try:
        preflight = _de_preflight(live)
    except KeyboardInterrupt:
        return _search_execution(
            live,
            DeSearchTerminal.INTERRUPTED,
            None,
            failure=TerminalFailure(
                "interrupted",
                "KeyboardInterrupt during DE preflight",
            ),
        )
    if isinstance(preflight, TerminalFailure):
        return _search_execution(
            live,
            DeSearchTerminal.PREFLIGHT_INVALID,
            None,
            failure=preflight,
        )
    if cancellation.is_cancelled:
        return _search_execution(
            live,
            DeSearchTerminal.CANCELLED,
            preflight.identity,
            failure=TerminalFailure("cancelled", "Cancellation after DE preflight"),
        )
    solver_start = np.asarray(
        tuple(
            coordinate.to_solver(value)
            for coordinate, value in zip(
                invocation.search_coordinates,
                invocation.search_problem.start,
                strict=True,
            )
        ),
        dtype=np.float64,
    )
    try:
        result = _invoke_de_backend(live, invocation, solver_start)
    except _DeAttemptStop as stop:
        return _search_execution(
            live,
            stop.terminal,
            preflight.identity,
            failure=stop.failure,
        )
    except KeyboardInterrupt:
        return _search_execution(
            live,
            DeSearchTerminal.INTERRUPTED,
            preflight.identity,
            failure=TerminalFailure("interrupted", "KeyboardInterrupt during DE"),
        )
    except Exception as error:  # noqa: BLE001 - undeclared backend errors fail closed
        return _search_execution(
            live,
            DeSearchTerminal.IMPLEMENTATION_FAILURE,
            preflight.identity,
            failure=TerminalFailure(
                "unexpected_de_backend_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    if cancellation.is_cancelled:
        return _search_execution(
            live,
            DeSearchTerminal.CANCELLED,
            preflight.identity,
            failure=TerminalFailure("cancelled", "Cancellation after DE backend"),
        )
    return _finish_de_backend(live, result, preflight.identity)


def _derive_polish_problem(
    problem: OptimizationProblem,
    invocation: DeDirectTrfInvocation,
    search: DeSearchExecution,
) -> tuple[OptimizationProblem, DirectTrfInvocation]:
    candidate = search.best_candidate
    if not search.restart_eligible or candidate is None:
        raise DirectTrfConstructionError(
            "Only a restart-eligible DE outcome may create a TRF polish"
        )
    derivation = DePolishProblemDerivation(
        problem.identity,
        invocation.identity,
        invocation.search_problem.identity,
        search.identity,
        candidate.identity,
        problem.controlled_ids,
        candidate.full_vector,
    )
    child = OptimizationProblem(
        problem.evaluation_plan_identity,
        problem.parameterization_identity,
        problem.evaluator_parameterization_identity,
        problem.constraint_program_identity,
        problem.configuration_identity,
        problem.source_snapshot,
        problem.independent_items,
        problem.controlled_ids,
        problem.held_items,
        candidate.full_vector,
        problem.lower_bounds,
        problem.upper_bounds,
        problem.commit_scope,
        derivation,
        problem.scalarization_version,
    )
    polish = DirectTrfInvocation.for_problem(
        child,
        objective_request_budget=invocation.polish_objective_request_budget,
        x_scale=invocation.polish_x_scale,
        ftol=invocation.polish_ftol,
        xtol=invocation.polish_xtol,
        gtol=invocation.polish_gtol,
    )
    return child, polish


def _validate_de_transition_lineage(
    problem: OptimizationProblem,
    invocation: DeDirectTrfInvocation,
    search: DeSearchExecution,
) -> None:
    """Validate DE-owned candidate lineage before any TRF transfer exists."""
    candidate = search.best_candidate
    if candidate is None:
        raise DirectTrfConstructionError(
            "Restart-eligible DE search lacks a transition candidate"
        )
    selected_ids = invocation.search_problem.controlled_ids
    root_indices = {
        param_id: index for index, param_id in enumerate(problem.controlled_ids)
    }
    expected_selected = tuple(
        candidate.full_vector[root_indices[param_id]] for param_id in selected_ids
    )
    unselected_ids = tuple(
        param_id
        for param_id in problem.controlled_ids
        if param_id not in set(selected_ids)
    )
    root_start = dict(zip(problem.controlled_ids, problem.start, strict=True))
    expected_solver = tuple(
        coordinate.to_solver(value)
        for coordinate, value in zip(
            invocation.search_coordinates,
            candidate.selected_vector,
            strict=True,
        )
    )
    if (
        search.problem_identity != invocation.search_problem.identity
        or search.workflow_invocation_identity != invocation.identity
        or search.population_identity != invocation.population.identity
        or search.candidate_ordering_policy != _DE_CANDIDATE_ORDER_VERSION
        or len(candidate.full_vector) != len(problem.controlled_ids)
        or len(candidate.selected_vector) != len(selected_ids)
        or candidate.selected_vector != expected_selected
        or candidate.solver_vector != expected_solver
        or any(
            candidate.full_vector[root_indices[param_id]] != root_start[param_id]
            for param_id in unselected_ids
        )
        or candidate.request_ordinal > search.counters.objective_requests_accepted
        or search.valid_candidate_count
        > search.counters.objective_evaluations_completed
    ):
        raise DirectTrfConstructionError(
            "Eligible DE candidate lineage is incompatible with its root problem"
        )


def _validate_de_polish_transition_lineage(
    problem: OptimizationProblem,
    invocation: DeDirectTrfInvocation,
    search: DeSearchExecution,
    polish_problem: OptimizationProblem,
    polish_invocation: DirectTrfInvocation,
    polish: DirectTrfCandidateOutcome,
) -> None:
    """Bind successful Direct TRF evidence to this exact DE transition."""
    search_candidate = search.best_candidate
    polished_candidate = polish.candidate
    materialization = polish.materialization
    derivation = polish_problem.derivation
    expected_derivation = None
    if search_candidate is not None:
        expected_derivation = DePolishProblemDerivation(
            problem.identity,
            invocation.identity,
            invocation.search_problem.identity,
            search.identity,
            search_candidate.identity,
            problem.controlled_ids,
            search_candidate.full_vector,
        )
    if (
        search_candidate is None
        or polish.terminal is not DirectTrfCandidateTerminal.SUCCESS
        or polished_candidate is None
        or materialization is None
        or not isinstance(derivation, DePolishProblemDerivation)
        or derivation != expected_derivation
        or polish_problem.evaluation_plan_identity != problem.evaluation_plan_identity
        or polish_problem.parameterization_identity != problem.parameterization_identity
        or polish_problem.evaluator_parameterization_identity
        != problem.evaluator_parameterization_identity
        or polish_problem.constraint_program_identity
        != problem.constraint_program_identity
        or polish_problem.configuration_identity != problem.configuration_identity
        or polish_problem.source_snapshot != problem.source_snapshot
        or polish_problem.independent_items != problem.independent_items
        or polish_problem.controlled_ids != problem.controlled_ids
        or polish_problem.held_items != problem.held_items
        or polish_problem.start != search_candidate.full_vector
        or polish_problem.lower_bounds != problem.lower_bounds
        or polish_problem.upper_bounds != problem.upper_bounds
        or polish_problem.commit_scope != problem.commit_scope
        or polish_problem.scalarization_version != problem.scalarization_version
        or polish_invocation.problem_identity != polish_problem.identity
        or polish_invocation.objective_request_budget
        != invocation.polish_objective_request_budget
        or polish_invocation.x_scale != invocation.polish_x_scale
        or polish_invocation.ftol != invocation.polish_ftol
        or polish_invocation.xtol != invocation.polish_xtol
        or polish_invocation.gtol != invocation.polish_gtol
        or polish.execution.problem_identity != polish_problem.identity
        or polish.execution.invocation_identity != polish_invocation.identity
        or polished_candidate.problem_identity != polish_problem.identity
        or polished_candidate.invocation_identity != polish_invocation.identity
        or polished_candidate.execution_identity != polish.execution.identity
        or polish.materialization is not polished_candidate.materialization
        or materialization.problem_identity != polish_problem.identity
        or materialization.invocation_identity != polish_invocation.identity
        or materialization.execution_identity != polish.execution.identity
        or materialization.evaluation_identity
        != polished_candidate.evaluation_result.identity
    ):
        raise DirectTrfConstructionError(
            "Successful TRF candidate lacks exact DE polish transition lineage"
        )


@dataclass(frozen=True, slots=True)
class DeTrfTransitionAccounting:
    """Explicit non-fungible accounting owned by the DE→TRF transition."""

    de_budget: int
    de_counters: AttemptCounters
    polish_budget: int
    polish_counters: AttemptCounters | None
    search_to_polish_transfers: int
    root_materializations: int
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.search_to_polish_transfers not in {
            0,
            1,
        } or self.root_materializations not in {
            0,
            1,
        }:
            raise ValueError("DE→TRF transition counts must be zero or one")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-trf-transition-accounting",
                (
                    self.de_budget,
                    (
                        self.de_counters.solver_requests_received,
                        self.de_counters.objective_requests_accepted,
                        self.de_counters.objective_evaluations_completed,
                    ),
                    self.polish_budget,
                    None
                    if self.polish_counters is None
                    else (
                        self.polish_counters.solver_requests_received,
                        self.polish_counters.objective_requests_accepted,
                        self.polish_counters.objective_evaluations_completed,
                    ),
                    self.search_to_polish_transfers,
                    self.root_materializations,
                ),
            ),
        )


class DeDirectTrfTerminal(StrEnum):
    """Closed workflow outcomes; only ACCEPTED can expose fit authority."""

    ACCEPTED = "accepted"
    SEARCH_UNSUCCESSFUL = "search_unsuccessful"
    POLISH_UNSUCCESSFUL = "polish_unsuccessful"
    MATERIALIZATION_FAILURE = "materialization_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    EXECUTION_FAILURE = "execution_failure"


@dataclass(frozen=True, slots=True)
class DePolishProvenance:
    """Exact selected-restart and one-polish lineage for accepted evidence."""

    workflow_invocation_identity: str
    root_problem_identity: str
    search_execution_identity: str
    search_candidate_identity: str
    search_terminal: DeSearchTerminal
    polish_problem_identity: str
    polish_invocation_identity: str
    polish_execution_identity: str
    polished_candidate_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-de-polish-provenance",
                (
                    self.workflow_invocation_identity,
                    self.root_problem_identity,
                    self.search_execution_identity,
                    self.search_candidate_identity,
                    self.search_terminal.value,
                    self.polish_problem_identity,
                    self.polish_invocation_identity,
                    self.polish_execution_identity,
                    self.polished_candidate_identity,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class AcceptedDeDirectTrfResult(AcceptedFitResult):
    """Accepted root evidence retaining exact DE search and TRF polish lineage."""

    de_polish_provenance: DePolishProvenance
    workflow_invocation: DeDirectTrfInvocation = field(repr=False, compare=False)
    search_execution: DeSearchExecution = field(repr=False, compare=False)
    polish_problem: OptimizationProblem = field(repr=False, compare=False)
    polish_invocation: DirectTrfInvocation = field(repr=False, compare=False)
    polish_outcome: DirectTrfCandidateOutcome = field(repr=False, compare=False)
    fresh_candidate: MaterializedDirectTrfCandidate = field(
        repr=False,
        compare=False,
    )

    @classmethod
    def from_accepted(
        cls,
        accepted: AcceptedFitResult,
        provenance: DePolishProvenance,
        workflow_invocation: DeDirectTrfInvocation,
        search_execution: DeSearchExecution,
        polish_problem: OptimizationProblem,
        polish_invocation: DirectTrfInvocation,
        polish_outcome: DirectTrfCandidateOutcome,
        fresh_candidate: MaterializedDirectTrfCandidate,
    ) -> AcceptedDeDirectTrfResult:
        return cls(
            accepted.occurrence_identity,
            accepted.problem_identity,
            accepted.invocation_identity,
            accepted.execution_identity,
            accepted.materialization_identity,
            accepted.parameterization_identity,
            accepted.evaluator_parameterization_identity,
            accepted.source_occurrence_identity,
            accepted.source_revision,
            accepted.controlled_ids,
            accepted.vector,
            accepted.chi_square,
            accepted.evaluation_result,
            accepted.commit_scope,
            accepted.commit_items,
            accepted.origin_context_identity,
            provenance,
            workflow_invocation,
            search_execution,
            polish_problem,
            polish_invocation,
            polish_outcome,
            fresh_candidate,
        )


@dataclass(frozen=True, slots=True)
class DeDirectTrfOutcome:
    """Complete DE→TRF occurrence with separated stage and authority evidence."""

    terminal: DeDirectTrfTerminal
    search: DeSearchExecution
    accounting: DeTrfTransitionAccounting
    polish_problem: OptimizationProblem | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    polish_invocation: DirectTrfInvocation | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    polish: DirectTrfCandidateOutcome | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    root_materialization: CandidateMaterialization | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    accepted_result: AcceptedDeDirectTrfResult | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    commit_authority: LiveFitCommitAuthority | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    failure: TerminalFailure | None = None

    def __post_init__(self) -> None:
        accepted = self.terminal is DeDirectTrfTerminal.ACCEPTED
        if accepted != (
            self.accepted_result is not None and self.commit_authority is not None
        ):
            raise ValueError("Only accepted DE→TRF exposes fit commit authority")
        if self.polish is None and (
            self.polish_problem is not None or self.polish_invocation is not None
        ):
            raise ValueError("Unstarted TRF polish cannot expose partial topology")


class DeDirectTrfInterrupted(KeyboardInterrupt):
    """Propagated interruption carrying the frozen DE→TRF workflow outcome."""

    def __init__(self, outcome: DeDirectTrfOutcome) -> None:
        self.outcome = outcome
        super().__init__("Native DE→TRF was interrupted")


def _transition_accounting(
    invocation: DeDirectTrfInvocation,
    search: DeSearchExecution,
    polish: DirectTrfCandidateOutcome | None,
    *,
    transferred: bool,
    materialized: bool,
) -> DeTrfTransitionAccounting:
    return DeTrfTransitionAccounting(
        invocation.de_objective_request_budget,
        search.counters,
        invocation.polish_objective_request_budget,
        None if polish is None else polish.execution.counters,
        int(transferred),
        int(materialized),
    )


def _gate_de_to_trf_transition(
    invocation: DeDirectTrfInvocation,
    search: DeSearchExecution,
    token: CancellationToken,
) -> DeDirectTrfOutcome | None:
    try:
        cancelled_before_polish = token.is_cancelled
    except KeyboardInterrupt as error:
        outcome = DeDirectTrfOutcome(
            DeDirectTrfTerminal.INTERRUPTED,
            search,
            _transition_accounting(
                invocation,
                search,
                None,
                transferred=False,
                materialized=False,
            ),
            failure=TerminalFailure(
                "interrupted",
                "KeyboardInterrupt at DE to TRF transition gate",
            ),
        )
        raise DeDirectTrfInterrupted(outcome) from error
    if not cancelled_before_polish:
        return None
    return DeDirectTrfOutcome(
        DeDirectTrfTerminal.CANCELLED,
        search,
        _transition_accounting(
            invocation,
            search,
            None,
            transferred=False,
            materialized=False,
        ),
        failure=TerminalFailure(
            "cancelled",
            "Cancellation observed before DE to TRF transfer",
        ),
    )


def execute_de_direct_trf(
    problem: OptimizationProblem,
    invocation: DeDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
) -> DeDirectTrfOutcome:
    """Run one selected-coordinate DE stage and at most one native TRF polish."""
    _validate_de_context(problem, invocation, parameterization, engine)
    token = CancellationToken() if cancellation is None else cancellation
    search = _run_de_search(problem, invocation, parameterization, engine, token)
    if not search.restart_eligible:
        terminal = {
            DeSearchTerminal.CANCELLED: DeDirectTrfTerminal.CANCELLED,
            DeSearchTerminal.INTERRUPTED: DeDirectTrfTerminal.INTERRUPTED,
            DeSearchTerminal.BACKEND_FAILURE: DeDirectTrfTerminal.EXECUTION_FAILURE,
            DeSearchTerminal.IMPLEMENTATION_FAILURE: (
                DeDirectTrfTerminal.EXECUTION_FAILURE
            ),
        }.get(search.terminal, DeDirectTrfTerminal.SEARCH_UNSUCCESSFUL)
        outcome = DeDirectTrfOutcome(
            terminal,
            search,
            _transition_accounting(
                invocation,
                search,
                None,
                transferred=False,
                materialized=False,
            ),
            failure=search.failure,
        )
        if terminal is DeDirectTrfTerminal.INTERRUPTED:
            raise DeDirectTrfInterrupted(outcome)
        return outcome
    try:
        _validate_de_transition_lineage(problem, invocation, search)
    except Exception as error:  # noqa: BLE001 - transition lineage fails closed
        return DeDirectTrfOutcome(
            DeDirectTrfTerminal.EXECUTION_FAILURE,
            search,
            _transition_accounting(
                invocation,
                search,
                None,
                transferred=False,
                materialized=False,
            ),
            failure=TerminalFailure(
                "de_search_lineage_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    transition_failure = _gate_de_to_trf_transition(invocation, search, token)
    if transition_failure is not None:
        return transition_failure
    polish_problem, polish_invocation = _derive_polish_problem(
        problem,
        invocation,
        search,
    )
    try:
        polish = execute_direct_trf_candidate(
            polish_problem,
            polish_invocation,
            parameterization,
            engine,
            cancellation=token,
        )
    except DirectTrfInterrupted as error:
        polish = DirectTrfCandidateOutcome(
            DirectTrfCandidateTerminal.CANCELLED,
            error.execution,
            error.materialization,
        )
        outcome = DeDirectTrfOutcome(
            DeDirectTrfTerminal.INTERRUPTED,
            search,
            _transition_accounting(
                invocation,
                search,
                polish,
                transferred=True,
                materialized=False,
            ),
            polish_problem,
            polish_invocation,
            polish,
            failure=TerminalFailure(
                "interrupted", "KeyboardInterrupt during TRF polish"
            ),
        )
        raise DeDirectTrfInterrupted(outcome) from error
    if polish.terminal is not DirectTrfCandidateTerminal.SUCCESS:
        terminal = {
            DirectTrfCandidateTerminal.CANCELLED: DeDirectTrfTerminal.CANCELLED,
        }.get(polish.terminal, DeDirectTrfTerminal.POLISH_UNSUCCESSFUL)
        return DeDirectTrfOutcome(
            terminal,
            search,
            _transition_accounting(
                invocation,
                search,
                polish,
                transferred=True,
                materialized=False,
            ),
            polish_problem,
            polish_invocation,
            polish,
            failure=polish.execution.failure,
        )
    try:
        _validate_de_polish_transition_lineage(
            problem,
            invocation,
            search,
            polish_problem,
            polish_invocation,
            polish,
        )
    except Exception as error:  # noqa: BLE001 - polish lineage fails closed
        return DeDirectTrfOutcome(
            DeDirectTrfTerminal.EXECUTION_FAILURE,
            search,
            _transition_accounting(
                invocation,
                search,
                polish,
                transferred=True,
                materialized=False,
            ),
            polish_problem,
            polish_invocation,
            polish,
            failure=TerminalFailure(
                "de_polish_lineage_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    polished_candidate = cast("MaterializedDirectTrfCandidate", polish.candidate)
    fresh = _materialize_derived_direct_trf_candidate_for_root(
        problem,
        polished_candidate,
        parameterization,
        engine,
        cancellation=token,
    )
    if isinstance(fresh, CandidateMaterialization):
        terminal = {
            MaterializationTerminal.CANCELLED: DeDirectTrfTerminal.CANCELLED,
            MaterializationTerminal.INTERRUPTED: DeDirectTrfTerminal.INTERRUPTED,
        }.get(fresh.terminal, DeDirectTrfTerminal.MATERIALIZATION_FAILURE)
        materialization = fresh
        outcome = DeDirectTrfOutcome(
            terminal,
            search,
            _transition_accounting(
                invocation,
                search,
                polish,
                transferred=True,
                materialized=True,
            ),
            polish_problem,
            polish_invocation,
            polish,
            materialization,
            failure=materialization.failure,
        )
        if terminal is DeDirectTrfTerminal.INTERRUPTED:
            raise DeDirectTrfInterrupted(outcome)
        return outcome
    search_candidate = cast("DeSearchCandidate", search.best_candidate)
    provenance = DePolishProvenance(
        invocation.identity,
        problem.identity,
        search.identity,
        search_candidate.identity,
        search.terminal,
        polish_problem.identity,
        polish_invocation.identity,
        polish.execution.identity,
        polished_candidate.identity,
    )
    accepted_base = _accept_materialized_fit_for_derived_workflow(
        problem=problem,
        invocation_identity=polish_invocation.identity,
        execution_identity=polish.execution.identity,
        materialization=fresh.materialization,
        vector=fresh.vector,
        chi_square=fresh.chi_square,
        evaluation_result=fresh.evaluation_result,
        authority_context_identity=provenance.identity,
    )
    accepted = AcceptedDeDirectTrfResult.from_accepted(
        accepted_base,
        provenance,
        invocation,
        search,
        polish_problem,
        polish_invocation,
        polish,
        fresh,
    )
    authority = _grant_derived_fit_commit_authority(accepted, problem)
    return DeDirectTrfOutcome(
        DeDirectTrfTerminal.ACCEPTED,
        search,
        _transition_accounting(
            invocation,
            search,
            polish,
            transferred=True,
            materialized=True,
        ),
        polish_problem,
        polish_invocation,
        polish,
        fresh.materialization,
        accepted,
        authority,
    )


def validate_de_commit_lineage(
    accepted: AcceptedDeDirectTrfResult,
    problem: OptimizationProblem,
) -> None:
    """Fail closed unless accepted evidence retains exact DE→TRF lineage."""
    provenance = accepted.de_polish_provenance
    invocation = accepted.workflow_invocation
    search = accepted.search_execution
    polish_problem = accepted.polish_problem
    polish_invocation = accepted.polish_invocation
    polish = accepted.polish_outcome
    fresh = accepted.fresh_candidate
    search_candidate = search.best_candidate
    polished_candidate = polish.candidate
    derivation = polish_problem.derivation
    if (
        accepted.problem_identity != problem.identity
        or accepted.origin_context_identity != provenance.identity
        or provenance.workflow_invocation_identity != invocation.identity
        or provenance.root_problem_identity != problem.identity
        or provenance.search_execution_identity != search.identity
        or search_candidate is None
        or provenance.search_candidate_identity != search_candidate.identity
        or provenance.search_terminal is not search.terminal
        or not search.restart_eligible
        or not isinstance(derivation, DePolishProblemDerivation)
        or derivation.root_problem_identity != problem.identity
        or derivation.workflow_invocation_identity != invocation.identity
        or derivation.search_execution_identity != search.identity
        or derivation.search_candidate_identity != search_candidate.identity
        or provenance.polish_problem_identity != polish_problem.identity
        or provenance.polish_invocation_identity != polish_invocation.identity
        or provenance.polish_execution_identity != polish.execution.identity
        or polish.terminal is not DirectTrfCandidateTerminal.SUCCESS
        or polished_candidate is None
        or provenance.polished_candidate_identity != polished_candidate.identity
        or fresh.problem_identity != problem.identity
        or accepted.invocation_identity != polish_invocation.identity
        or accepted.execution_identity != polish.execution.identity
        or accepted.materialization_identity != fresh.materialization.identity
        or accepted.vector != fresh.vector
        or accepted.chi_square != fresh.chi_square
        or accepted.evaluation_result is not fresh.evaluation_result
    ):
        raise DirectTrfConstructionError(
            "Accepted DE→TRF result lacks exact canonical polish authority"
        )


def commit_de_accepted_fit(
    accepted: AcceptedDeDirectTrfResult,
    authority: LiveFitCommitAuthority,
    *,
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    analysis_values: AnalysisValues,
) -> CommitReceipt:
    """Commit one exact accepted DE→TRF result through the native boundary."""
    if not isinstance(accepted, AcceptedDeDirectTrfResult):
        raise DirectTrfConstructionError(
            "DE→TRF commit requires accepted workflow evidence"
        )
    validate_de_commit_lineage(accepted, problem)
    return commit_accepted_fit(
        accepted,
        authority,
        problem=problem,
        parameterization=parameterization,
        analysis_values=analysis_values,
    )
