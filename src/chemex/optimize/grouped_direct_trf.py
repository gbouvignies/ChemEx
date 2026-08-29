"""Exact native fit-component decomposition for grouped Direct TRF.

The production deterministic dispatcher composes this closed decomposition with
the native Direct TRF executor and one aggregate commit.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np

from chemex.evaluation.native import EvaluationEngine, EvaluationPlan, EvaluationResult
from chemex.optimize import direct_trf as direct_trf_owner
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    CancellationToken,
    CandidateMaterialization,
    ComponentProblemDerivation,
    DirectTrfCandidateTerminal,
    DirectTrfConstructionError,
    DirectTrfInterrupted,
    DirectTrfInvocation,
    DirectTrfTerminal,
    FinalResidualJacobianEvidence,
    LiveFitCommitAuthority,
    MaterializationTerminal,
    MaterializedDirectTrfCandidate,
    ObjectiveScalarizationError,
    OptimizationProblem,
    ResidualJacobianSource,
    TerminalFailure,
    accept_materialized_fit,
    canonical_chi_square,
    execute_direct_trf_candidate,
)
from chemex.optimize.progress import (
    ContextualProgressObserver,
    FitProgressContext,
    ProgressObserver,
)
from chemex.parameters.parameterization import ActiveParameterization

_SCHEMA_VERSION = 1
_DECOMPOSITION_POLICY = "profile-constraint-connectivity-v1"
_PROJECTION_POLICY = "complete-profile-root-order-v1"


class _GroupedValidationCancelled(Exception):
    """Internal control flow for cancellation during context recompilation."""


def _check_grouped_cancellation(cancellation: CancellationToken | None) -> None:
    if cancellation is not None and cancellation.is_cancelled:
        raise _GroupedValidationCancelled


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "record": record, "schema_version": _SCHEMA_VERSION},
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class FitComponent:
    """One deterministic exact child derived from a complete root problem."""

    controlled_ids: tuple[str, ...]
    root_profile_indices: tuple[int, ...]
    problem: OptimizationProblem = field(repr=False, compare=False)
    identity: str


class _ControlledDependencyResolver:
    def __init__(
        self,
        controlled_ids: tuple[str, ...],
        parameterization: ActiveParameterization,
    ) -> None:
        self._controlled = frozenset(controlled_ids)
        self._constraints = {
            constraint.target_id: constraint
            for constraint in parameterization.program.constraints
        }
        self._cache: dict[str, frozenset[str]] = {}
        self._path_cache: dict[
            str,
            tuple[tuple[str, tuple[str, ...]], ...],
        ] = {}

    def resolve(self, param_id: str) -> frozenset[str]:
        if param_id in self._cache:
            return self._cache[param_id]
        if param_id in self._controlled:
            result = frozenset((param_id,))
        elif constraint := self._constraints.get(param_id):
            result = frozenset(
                dependency
                for source in constraint.dependencies
                for dependency in self.resolve(source)
            )
        else:
            result = frozenset()
        self._cache[param_id] = result
        return result

    def paths(self, param_id: str) -> tuple[tuple[str, tuple[str, ...]], ...]:
        """Return every declared path from one resolved input to a control."""
        if param_id in self._path_cache:
            return self._path_cache[param_id]
        if param_id in self._controlled:
            result = ((param_id, (param_id,)),)
        elif constraint := self._constraints.get(param_id):
            result = tuple(
                (controlled_id, (param_id, *path))
                for dependency in constraint.dependencies
                for controlled_id, path in self.paths(dependency)
            )
        else:
            result = ()
        self._path_cache[param_id] = result
        return result


def _profile_dependencies(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> tuple[frozenset[str], ...]:
    resolver = _ControlledDependencyResolver(problem.controlled_ids, parameterization)
    return tuple(
        (
            frozenset(
                controlled_id
                for param_id in profile.param_ids
                for controlled_id in resolver.resolve(param_id)
            )
            if profile.retained_observation_indices
            else frozenset()
        )
        for profile in engine.plan.profiles
    )


def _profile_dependency_paths(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> tuple[tuple[tuple[str, str, tuple[str, ...]], ...], ...]:
    resolver = _ControlledDependencyResolver(problem.controlled_ids, parameterization)
    return tuple(
        (
            tuple(
                (param_id, controlled_id, path)
                for param_id in profile.param_ids
                for controlled_id, path in resolver.paths(param_id)
            )
            if profile.retained_observation_indices
            else ()
        )
        for profile in engine.plan.profiles
    )


def _ordered_component_controls(
    controlled_ids: tuple[str, ...],
    profile_dependencies: tuple[frozenset[str], ...],
    feasibility_dependencies: tuple[frozenset[str], ...] = (),
) -> tuple[tuple[str, ...], ...]:
    parent = {param_id: param_id for param_id in controlled_ids}

    def find(param_id: str) -> str:
        while parent[param_id] != param_id:
            parent[param_id] = parent[parent[param_id]]
            param_id = parent[param_id]
        return param_id

    for dependencies in (*profile_dependencies, *feasibility_dependencies):
        ordered = tuple(
            param_id for param_id in controlled_ids if param_id in dependencies
        )
        for param_id in ordered[1:]:
            left_root = find(ordered[0])
            right_root = find(param_id)
            if left_root != right_root:
                parent[right_root] = left_root
    observed = set().union(*profile_dependencies)
    missing = set(controlled_ids) - observed
    if missing:
        raise DirectTrfConstructionError(
            "Controlled parameters lack a retained-objective dependency: "
            + ", ".join(sorted(missing))
        )
    components = {
        tuple(
            candidate for candidate in controlled_ids if find(candidate) == find(root)
        )
        for root in controlled_ids
    }
    return tuple(
        sorted(
            components,
            key=lambda ids: tuple(param_id.encode("utf-8") for param_id in ids),
        )
    )


def _build_component(
    problem: OptimizationProblem,
    engine: EvaluationEngine,
    component_ids: tuple[str, ...],
    profile_dependencies: tuple[frozenset[str], ...],
) -> FitComponent:
    component_set = set(component_ids)
    profile_indices = tuple(
        index
        for index, dependencies in enumerate(profile_dependencies)
        if dependencies and dependencies.issubset(component_set)
    )
    child_plan = engine.project_plan(profile_indices)
    bounds = {
        param_id: (start, lower, upper)
        for param_id, start, lower, upper in zip(
            problem.controlled_ids,
            problem.start,
            problem.lower_bounds,
            problem.upper_bounds,
            strict=True,
        )
    }
    held_items = tuple(
        item for item in problem.independent_items if item[0] not in component_set
    )
    start = tuple(bounds[param_id][0] for param_id in component_ids)
    lower = tuple(bounds[param_id][1] for param_id in component_ids)
    upper = tuple(bounds[param_id][2] for param_id in component_ids)
    component_identity = _component_identity(
        problem,
        engine,
        component_ids,
        profile_indices,
        lower,
        upper,
    )
    derivation = ComponentProblemDerivation(
        problem.identity,
        problem.affine_feasibility_identity,
        component_identity,
        _PROJECTION_POLICY,
        child_plan.identity,
        component_ids,
        held_items,
    )
    child_problem = problem.derive_grouped_component(
        controlled_ids=component_ids,
        start=start,
        derivation=derivation,
    )
    return FitComponent(
        component_ids,
        profile_indices,
        child_problem,
        component_identity,
    )


def _component_identity(
    problem: OptimizationProblem,
    engine: EvaluationEngine,
    component_ids: tuple[str, ...],
    profile_indices: tuple[int, ...],
    lower_bounds: tuple[float, ...],
    upper_bounds: tuple[float, ...],
) -> str:
    """Derive the canonical identity of one exact grouped root projection."""
    return _identity(
        "native-fit-component",
        (
            _DECOMPOSITION_POLICY,
            _PROJECTION_POLICY,
            problem.evaluator_parameterization_identity,
            problem.constraint_program_identity,
            component_ids,
            tuple(
                engine.plan.profiles[index].source_identity for index in profile_indices
            ),
            engine.plan.identity,
            profile_indices,
            tuple(float(value).hex() for value in lower_bounds),
            tuple(float(value).hex() for value in upper_bounds),
        ),
    )


@dataclass(frozen=True, slots=True)
class FitPartitionProof:
    """Verifiable exact partition of current Profile objectives and box domain."""

    root_plan_identity: str
    constraint_program_identity: str
    controlled_ids: tuple[str, ...]
    profile_records: tuple[
        tuple[
            str,
            tuple[str, ...],
            tuple[str, ...],
            tuple[tuple[str, str, tuple[str, ...]], ...],
        ],
        ...,
    ]
    component_records: tuple[
        tuple[
            str,
            tuple[str, ...],
            tuple[int, ...],
            tuple[tuple[str, str, str], ...],
        ],
        ...,
    ]
    constant_profile_indices: tuple[int, ...]
    non_objective_profile_indices: tuple[int, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        assigned_profiles = (
            [
                index
                for _component_id, _controlled, indices, _bounds in self.component_records
                for index in indices
            ]
            + list(self.constant_profile_indices)
            + list(self.non_objective_profile_indices)
        )
        assigned_controls = [
            param_id
            for _component_id, controlled, _indices, _bounds in self.component_records
            for param_id in controlled
        ]
        if (
            sorted(assigned_profiles) != list(range(len(self.profile_records)))
            or len(assigned_profiles) != len(set(assigned_profiles))
            or set(assigned_controls) != set(self.controlled_ids)
            or len(assigned_controls) != len(set(assigned_controls))
        ):
            raise DirectTrfConstructionError("Fit partition is incomplete or overlaps")
        for _component_id, controlled, indices, bounds in self.component_records:
            controlled_set = set(controlled)
            if tuple(item[0] for item in bounds) != controlled:
                raise DirectTrfConstructionError(
                    "Fit partition bounds do not match component coordinates"
                )
            if any(
                not self.profile_records[index][2]
                or not set(self.profile_records[index][2]).issubset(controlled_set)
                for index in indices
            ):
                raise DirectTrfConstructionError(
                    "Fit partition component does not cover Profile dependencies"
                )
        if any(
            self.profile_records[index][2] for index in self.constant_profile_indices
        ):
            raise DirectTrfConstructionError(
                "Constant objective partition contains a controlled dependency"
            )
        if any(
            self.profile_records[index][2] or self.profile_records[index][3]
            for index in self.non_objective_profile_indices
        ):
            raise DirectTrfConstructionError(
                "Non-objective Profile partition contains a retained dependency"
            )
        for _profile_id, param_ids, controlled_ids, paths in self.profile_records:
            path_controls = {controlled_id for _local, controlled_id, _path in paths}
            if path_controls != set(controlled_ids) or any(
                local_id not in param_ids
                or not path
                or path[0] != local_id
                or path[-1] != controlled_id
                for local_id, controlled_id, path in paths
            ):
                raise DirectTrfConstructionError(
                    "Fit partition dependency provenance is incomplete"
                )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-fit-partition-proof",
                (
                    _DECOMPOSITION_POLICY,
                    _PROJECTION_POLICY,
                    "independent-box-domain-v1",
                    self.root_plan_identity,
                    self.constraint_program_identity,
                    self.controlled_ids,
                    self.profile_records,
                    self.component_records,
                    self.constant_profile_indices,
                    self.non_objective_profile_indices,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class FitDecomposition:
    """Maximal exact Profile/objective factorization of one native root fit."""

    root_problem_identity: str
    root_plan_identity: str
    components: tuple[FitComponent, ...]
    constant_profile_indices: tuple[int, ...]
    non_objective_profile_indices: tuple[int, ...]
    partition_proof: FitPartitionProof
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        expected_component_records = tuple(
            (
                component.identity,
                component.controlled_ids,
                component.root_profile_indices,
                tuple(
                    (param_id, float(lower).hex(), float(upper).hex())
                    for param_id, lower, upper in zip(
                        component.controlled_ids,
                        component.problem.lower_bounds,
                        component.problem.upper_bounds,
                        strict=True,
                    )
                ),
            )
            for component in self.components
        )
        if (
            self.partition_proof.root_plan_identity != self.root_plan_identity
            or self.partition_proof.component_records != expected_component_records
            or self.partition_proof.constant_profile_indices
            != self.constant_profile_indices
            or self.partition_proof.non_objective_profile_indices
            != self.non_objective_profile_indices
        ):
            raise DirectTrfConstructionError(
                "Fit decomposition differs from its exact partition proof"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-fit-decomposition",
                (
                    _DECOMPOSITION_POLICY,
                    _PROJECTION_POLICY,
                    self.root_plan_identity,
                    self.partition_proof.constraint_program_identity,
                    tuple(component.identity for component in self.components),
                    self.constant_profile_indices,
                    self.non_objective_profile_indices,
                    self.partition_proof.identity,
                ),
            ),
        )

    @classmethod
    def from_root(
        cls,
        problem: OptimizationProblem,
        parameterization: ActiveParameterization,
        engine: EvaluationEngine,
        *,
        cancellation: CancellationToken | None = None,
    ) -> FitDecomposition:
        """Derive exact components from semantic Profile and constraint paths."""
        if not problem.acceptance_authority:
            raise DirectTrfConstructionError(
                "Grouped decomposition requires one complete root problem"
            )
        problem.validate_parameterization(parameterization)
        if engine.plan.identity != problem.evaluation_plan_identity:
            raise DirectTrfConstructionError(
                "Root evaluator belongs to another problem"
            )

        profile_dependencies = _profile_dependencies(
            problem,
            parameterization,
            engine,
        )
        profile_paths = _profile_dependency_paths(
            problem,
            parameterization,
            engine,
        )
        components_list: list[FitComponent] = []
        feasible = problem.feasible_coordinates
        feasibility_dependencies = (
            () if feasible is None else feasible.controlled_domain_groups
        )
        for component_ids in _ordered_component_controls(
            problem.controlled_ids,
            profile_dependencies,
            feasibility_dependencies,
        ):
            _check_grouped_cancellation(cancellation)
            component = _build_component(
                problem,
                engine,
                component_ids,
                profile_dependencies,
            )
            _check_grouped_cancellation(cancellation)
            components_list.append(component)
        components = tuple(components_list)

        constant_indices = tuple(
            index
            for index, dependencies in enumerate(profile_dependencies)
            if not dependencies
            and engine.plan.profiles[index].retained_observation_indices
        )
        non_objective_indices = tuple(
            index
            for index, profile in enumerate(engine.plan.profiles)
            if not profile.retained_observation_indices
        )
        partition_proof = FitPartitionProof(
            engine.plan.identity,
            problem.constraint_program_identity,
            problem.controlled_ids,
            tuple(
                (
                    profile.identity,
                    profile.param_ids,
                    tuple(
                        param_id
                        for param_id in problem.controlled_ids
                        if param_id in profile_dependencies[index]
                    ),
                    profile_paths[index],
                )
                for index, profile in enumerate(engine.plan.profiles)
            ),
            tuple(
                (
                    component.identity,
                    component.controlled_ids,
                    component.root_profile_indices,
                    tuple(
                        (
                            param_id,
                            float(lower).hex(),
                            float(upper).hex(),
                        )
                        for param_id, lower, upper in zip(
                            component.controlled_ids,
                            component.problem.lower_bounds,
                            component.problem.upper_bounds,
                            strict=True,
                        )
                    ),
                )
                for component in components
            ),
            constant_indices,
            non_objective_indices,
        )
        return cls(
            problem.identity,
            engine.plan.identity,
            components,
            constant_indices,
            non_objective_indices,
            partition_proof,
        )


@dataclass(frozen=True, slots=True)
class GroupedDirectTrfInvocation:
    """Explicit immutable per-component Direct TRF budget allocation."""

    root_problem_identity: str
    decomposition_identity: str
    component_invocations: tuple[DirectTrfInvocation, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-grouped-direct-trf-invocation",
                (
                    self.root_problem_identity,
                    self.decomposition_identity,
                    tuple(item.identity for item in self.component_invocations),
                ),
            ),
        )

    @classmethod
    def for_decomposition(
        cls,
        decomposition: FitDecomposition,
        *,
        objective_request_budgets: Sequence[int],
    ) -> GroupedDirectTrfInvocation:
        """Bind one declared objective-request budget to every component."""
        budgets = tuple(objective_request_budgets)
        if len(budgets) != len(decomposition.components):
            raise DirectTrfConstructionError(
                "Grouped Direct TRF requires one explicit budget per component"
            )
        invocations = tuple(
            DirectTrfInvocation.for_problem(
                component.problem,
                objective_request_budget=budget,
            )
            for component, budget in zip(
                decomposition.components,
                budgets,
                strict=True,
            )
        )
        return cls(
            decomposition.root_problem_identity,
            decomposition.identity,
            invocations,
        )


class ComponentDisposition(StrEnum):
    """Closed lifecycle disposition for every required fit component."""

    SUCCEEDED = "succeeded"
    FAILED = "failed"
    EXECUTION_FAILURE = "execution_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    NOT_STARTED = "not_started"


@dataclass(frozen=True, slots=True)
class FitComponentOutcome:
    """Non-authoritative component evidence in canonical decomposition order."""

    component_identity: str
    controlled_ids: tuple[str, ...]
    disposition: ComponentDisposition
    execution_identity: str | None = None
    candidate: MaterializedDirectTrfCandidate | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    evaluation_plan: EvaluationPlan | None = field(
        default=None,
        repr=False,
        compare=False,
        kw_only=True,
    )
    final_residual_jacobian: FinalResidualJacobianEvidence | None = field(
        default=None,
        repr=False,
        compare=False,
        kw_only=True,
    )
    failure: TerminalFailure | None = None
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if (self.disposition is ComponentDisposition.SUCCEEDED) != (
            self.candidate is not None and self.evaluation_plan is not None
        ):
            raise ValueError(
                "Only a successful component may expose a candidate and plan"
            )
        if self.final_residual_jacobian is not None and (
            self.candidate is None
            or self.final_residual_jacobian.controlled_ids != self.controlled_ids
            or self.final_residual_jacobian.final_vector != self.candidate.vector
            or self.final_residual_jacobian.final_residuals
            != tuple(
                float(value) for value in self.candidate.evaluation_result.residuals
            )
        ):
            raise ValueError("Component final Jacobian differs from its candidate")
        if self.disposition is ComponentDisposition.NOT_STARTED and (
            self.execution_identity is not None or self.failure is not None
        ):
            raise ValueError("An unstarted component cannot expose attempt evidence")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-fit-component-outcome",
                (
                    self.component_identity,
                    self.controlled_ids,
                    self.disposition.value,
                    self.execution_identity,
                    None if self.candidate is None else self.candidate.identity,
                    None
                    if self.evaluation_plan is None
                    else self.evaluation_plan.identity,
                    None
                    if self.final_residual_jacobian is None
                    else self.final_residual_jacobian.identity,
                    None if self.failure is None else self.failure.identity,
                ),
            ),
        )


class GroupedDirectTrfTerminal(StrEnum):
    """Closed authoritative outcomes of grouped native Direct TRF."""

    ACCEPTED = "accepted"
    COMPONENT_FAILURE = "component_failure"
    EXECUTION_FAILURE = "execution_failure"
    DECOMPOSITION_VALIDATION_FAILURE = "decomposition_validation_failure"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"


@dataclass(frozen=True, slots=True)
class GroupedDirectTrfOutcome:
    """Complete grouped occurrence; only ACCEPTED may carry root authority."""

    terminal: GroupedDirectTrfTerminal
    components: tuple[FitComponentOutcome, ...]
    accepted_result: AcceptedFitResult | None = field(
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
        if (self.terminal is GroupedDirectTrfTerminal.ACCEPTED) != (
            self.accepted_result is not None and self.commit_authority is not None
        ):
            raise ValueError("Only accepted grouped execution has root fit authority")
        if self.terminal is not GroupedDirectTrfTerminal.ACCEPTED and (
            self.accepted_result is not None or self.commit_authority is not None
        ):
            raise ValueError("Unaccepted grouped execution cannot expose fit authority")


def _grouped_context_validation_identity(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> str:
    return _identity(
        "native-grouped-context-validation",
        (
            problem.identity,
            decomposition.identity,
            invocation.identity,
            parameterization.identity,
            engine.plan.identity,
        ),
    )


def _validate_grouped_context(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    cancellation: CancellationToken,
) -> str:
    _check_grouped_cancellation(cancellation)
    problem.validate_parameterization(parameterization)
    _check_grouped_cancellation(cancellation)
    if engine.plan.identity != problem.evaluation_plan_identity:
        raise DirectTrfConstructionError("Root evaluator belongs to another problem")
    # Recheck the cheap semantic partition records without rebuilding child
    # evaluators or their feasibility charts. Component execution validates each
    # projected plan against the already sealed child problem before evaluation.
    profile_dependencies = _profile_dependencies(
        problem,
        parameterization,
        engine,
    )
    _check_grouped_cancellation(cancellation)
    profile_paths = _profile_dependency_paths(
        problem,
        parameterization,
        engine,
    )
    feasible = problem.feasible_coordinates
    feasibility_dependencies = (
        () if feasible is None else feasible.controlled_domain_groups
    )
    expected_component_controls = _ordered_component_controls(
        problem.controlled_ids,
        profile_dependencies,
        feasibility_dependencies,
    )
    expected_component_profiles = tuple(
        tuple(
            index
            for index, dependencies in enumerate(profile_dependencies)
            if dependencies and dependencies.issubset(set(component_ids))
        )
        for component_ids in expected_component_controls
    )
    expected_constant_indices = tuple(
        index
        for index, dependencies in enumerate(profile_dependencies)
        if not dependencies and engine.plan.profiles[index].retained_observation_indices
    )
    expected_non_objective_indices = tuple(
        index
        for index, profile in enumerate(engine.plan.profiles)
        if not profile.retained_observation_indices
    )
    expected_profile_records = tuple(
        (
            profile.identity,
            profile.param_ids,
            tuple(
                param_id
                for param_id in problem.controlled_ids
                if param_id in profile_dependencies[index]
            ),
            profile_paths[index],
        )
        for index, profile in enumerate(engine.plan.profiles)
    )
    proof = decomposition.partition_proof
    if (
        decomposition.root_problem_identity != problem.identity
        or decomposition.root_plan_identity != engine.plan.identity
        or proof.root_plan_identity != engine.plan.identity
        or proof.constraint_program_identity != problem.constraint_program_identity
        or proof.controlled_ids != problem.controlled_ids
        or proof.profile_records != expected_profile_records
        or tuple(component.controlled_ids for component in decomposition.components)
        != expected_component_controls
        or tuple(
            component.root_profile_indices for component in decomposition.components
        )
        != expected_component_profiles
        or decomposition.constant_profile_indices != expected_constant_indices
        or decomposition.non_objective_profile_indices != expected_non_objective_indices
        or invocation.root_problem_identity != problem.identity
        or invocation.decomposition_identity != decomposition.identity
        or len(invocation.component_invocations) != len(decomposition.components)
    ):
        raise DirectTrfConstructionError(
            "Grouped Direct TRF context is not rooted in one problem"
        )
    for component, component_invocation in zip(
        decomposition.components,
        invocation.component_invocations,
        strict=True,
    ):
        _check_grouped_cancellation(cancellation)
        problem.validate_derived_problem(component.problem)
        derivation = component.problem.derivation
        root_indices = {
            param_id: index for index, param_id in enumerate(problem.controlled_ids)
        }
        expected_start = tuple(
            problem.start[root_indices[param_id]]
            for param_id in component.controlled_ids
        )
        lower_bounds = tuple(
            problem.lower_bounds[root_indices[param_id]]
            for param_id in component.controlled_ids
        )
        upper_bounds = tuple(
            problem.upper_bounds[root_indices[param_id]]
            for param_id in component.controlled_ids
        )
        expected_identity = _component_identity(
            problem,
            engine,
            component.controlled_ids,
            component.root_profile_indices,
            lower_bounds,
            upper_bounds,
        )
        expected_feasible = problem.project_grouped_feasibility(
            controlled_ids=component.controlled_ids,
            start=expected_start,
            lower_bounds=lower_bounds,
            upper_bounds=upper_bounds,
        )
        actual_feasible = component.problem.feasible_coordinates
        if (
            not isinstance(derivation, ComponentProblemDerivation)
            or component.identity != expected_identity
            or derivation.component_identity != component.identity
            or derivation.projection_policy != _PROJECTION_POLICY
            or component.problem.start != expected_start
            or derivation.projected_plan_identity
            != component.problem.evaluation_plan_identity
            or (expected_feasible is None) != (actual_feasible is None)
            or (
                expected_feasible is not None
                and actual_feasible is not None
                and expected_feasible.identity != actual_feasible.identity
            )
            or component_invocation.problem_identity != component.problem.identity
        ):
            raise DirectTrfConstructionError(
                "Grouped Direct TRF context is not rooted in one problem"
            )
    return _grouped_context_validation_identity(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
    )


def _not_started(component: FitComponent) -> FitComponentOutcome:
    return FitComponentOutcome(
        component.identity,
        component.controlled_ids,
        ComponentDisposition.NOT_STARTED,
    )


def execute_direct_trf_components(
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
    progress_observer: ContextualProgressObserver | None = None,
    grid_seed_ordinal: int | None = None,
    grid_seed_total: int | None = None,
) -> tuple[FitComponentOutcome, ...]:
    """Execute canonical components without acceptance, commit, or publication."""
    if (
        invocation.root_problem_identity != decomposition.root_problem_identity
        or invocation.decomposition_identity != decomposition.identity
        or len(invocation.component_invocations) != len(decomposition.components)
    ):
        raise DirectTrfConstructionError(
            "Component invocation belongs to another decomposition"
        )
    token = CancellationToken() if cancellation is None else cancellation
    outcomes: list[FitComponentOutcome] = []
    stopped = token.is_cancelled
    component_total = len(decomposition.components)
    for component_ordinal, (component, component_invocation) in enumerate(
        zip(
            decomposition.components,
            invocation.component_invocations,
            strict=True,
        ),
        start=1,
    ):
        if stopped:
            outcomes.append(_not_started(component))
            continue
        try:
            child_engine = engine.project_profiles(component.root_profile_indices)
            if child_engine.plan.identity != component.problem.evaluation_plan_identity:
                raise DirectTrfConstructionError(
                    "Component evaluator projection differs from its derivation"
                )
            context = FitProgressContext(
                component_ordinal,
                component_total,
                component.controlled_ids,
                grid_seed_ordinal,
                grid_seed_total,
            )
            local_observer: ProgressObserver | None = (
                None
                if progress_observer is None
                else lambda event, context=context: progress_observer(context, event)
            )
            result = execute_direct_trf_candidate(
                component.problem,
                component_invocation,
                parameterization,
                child_engine,
                cancellation=token,
                progress_observer=local_observer,
            )
        except DirectTrfInterrupted as interrupted:
            interruption_failure = (
                interrupted.materialization.failure
                if interrupted.materialization is not None
                else interrupted.execution.failure
            )
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    ComponentDisposition.INTERRUPTED,
                    interrupted.execution.identity,
                    failure=interruption_failure,
                )
            )
            stopped = True
            continue
        except KeyboardInterrupt:
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    ComponentDisposition.INTERRUPTED,
                    failure=TerminalFailure(
                        "interrupted",
                        "KeyboardInterrupt during component projection",
                    ),
                )
            )
            stopped = True
            continue
        except Exception as error:  # noqa: BLE001 - projection failures stop safely
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    ComponentDisposition.EXECUTION_FAILURE,
                    failure=TerminalFailure(
                        "component_projection_failure",
                        f"{type(error).__name__}: {error}",
                    ),
                )
            )
            stopped = True
            continue
        if result.terminal is DirectTrfCandidateTerminal.SUCCESS:
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    ComponentDisposition.SUCCEEDED,
                    result.execution.identity,
                    result.candidate,
                    evaluation_plan=child_engine.plan,
                    final_residual_jacobian=(
                        None
                        if result.execution.backend is None
                        else result.execution.backend.final_residual_jacobian
                    ),
                )
            )
        elif result.terminal is DirectTrfCandidateTerminal.CANCELLED:
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    ComponentDisposition.CANCELLED,
                    result.execution.identity,
                    failure=result.execution.failure,
                )
            )
            stopped = True
        else:
            execution_failure = (
                result.terminal is DirectTrfCandidateTerminal.MATERIALIZATION_FAILURE
                or result.execution.terminal
                in {
                    DirectTrfTerminal.PREFLIGHT_INVALID,
                    DirectTrfTerminal.BACKEND_FAILURE,
                    DirectTrfTerminal.IMPLEMENTATION_FAILURE,
                }
            )
            outcomes.append(
                FitComponentOutcome(
                    component.identity,
                    component.controlled_ids,
                    (
                        ComponentDisposition.EXECUTION_FAILURE
                        if execution_failure
                        else ComponentDisposition.FAILED
                    ),
                    result.execution.identity,
                    failure=(
                        result.materialization.failure
                        if result.materialization is not None
                        and result.materialization.failure is not None
                        else result.execution.failure
                    ),
                )
            )
            stopped = execution_failure
    return tuple(outcomes)


type _ComponentProjection = tuple[EvaluationPlan, EvaluationResult]


def _materialized_component_projection(
    child_plan: EvaluationPlan,
    component: FitComponent,
    invocation: DirectTrfInvocation,
    outcome: FitComponentOutcome,
) -> _ComponentProjection | TerminalFailure:
    candidate = outcome.candidate
    if (
        outcome.controlled_ids != component.controlled_ids
        or candidate is None
        or len(candidate.vector) != len(component.controlled_ids)
    ):
        return TerminalFailure(
            "aggregate_assignment_failure",
            "A successful component lacks its exact controlled assignment",
        )
    if any(
        not lower <= value <= upper
        for value, lower, upper in zip(
            candidate.vector,
            component.problem.lower_bounds,
            component.problem.upper_bounds,
            strict=True,
        )
    ):
        return TerminalFailure(
            "aggregate_feasibility_failure",
            "Combined component assignment violates the root feasible domain",
        )
    child_result = candidate.evaluation_result
    expected_compatibility_identity = child_plan.compatibility_identity
    try:
        child_chi_square = canonical_chi_square(child_result.residuals)
    except (TypeError, ValueError, ObjectiveScalarizationError):
        return TerminalFailure(
            "decomposition_projection_mismatch",
            "Fresh root objective differs from component evidence",
        )
    if (
        child_result.plan_identity != child_plan.identity
        or child_plan.identity != component.problem.evaluation_plan_identity
        or child_result.parameterization_identity
        != component.problem.evaluator_parameterization_identity
        or tuple(child_result.resolved_values) != component.problem.commit_scope
        or candidate.problem_identity != component.problem.identity
        or candidate.invocation_identity != invocation.identity
        or candidate.execution_identity != outcome.execution_identity
        or candidate.chi_square != child_chi_square
        or candidate.materialization.terminal is not MaterializationTerminal.SUCCESS
        or candidate.materialization.problem_identity != component.problem.identity
        or candidate.materialization.invocation_identity != invocation.identity
        or candidate.materialization.execution_identity != outcome.execution_identity
        or candidate.materialization.candidate.vector != candidate.vector
        or candidate.materialization.candidate.chi_square != candidate.chi_square
        or candidate.materialization.evaluation_identity != child_result.identity
        or candidate.materialization.evaluator_compatibility_identity
        != expected_compatibility_identity
        or child_result.evaluator_compatibility_identity
        != expected_compatibility_identity
        or candidate.vector
        != tuple(
            child_result.resolved_values[param_id]
            for param_id in component.controlled_ids
        )
    ):
        return TerminalFailure(
            "decomposition_projection_mismatch",
            "Fresh root objective differs from component evidence",
        )
    return child_plan, child_result


def _root_projection_matches(
    root_result: EvaluationResult,
    root_engine: EvaluationEngine,
    component: FitComponent,
    projection: _ComponentProjection,
) -> bool:
    child_plan, child_result = projection
    for child_index, root_index in enumerate(component.root_profile_indices):
        root_descriptor = root_engine.plan.profiles[root_index]
        child_descriptor = child_plan.profiles[child_index]
        root_profile = root_result.profiles[root_index]
        child_profile = child_result.profiles[child_index]
        root_start = root_descriptor.observation_offset
        root_stop = root_start + root_descriptor.observation_count
        child_start = child_descriptor.observation_offset
        child_stop = child_start + child_descriptor.observation_count
        root_residual_start = root_profile.residual_offset
        root_residual_stop = root_residual_start + root_profile.residual_count
        child_residual_start = child_profile.residual_offset
        child_residual_stop = child_residual_start + child_profile.residual_count
        if (
            root_descriptor.source_identity != child_descriptor.source_identity
            or root_profile.normalization_factor != child_profile.normalization_factor
            or root_profile.retained_observation_indices
            != child_profile.retained_observation_indices
            or root_profile.kernel_identity != child_profile.kernel_identity
            or not np.array_equal(
                root_result.unscaled_calculations[root_start:root_stop],
                child_result.unscaled_calculations[child_start:child_stop],
            )
            or not np.array_equal(
                root_result.normalized_calculations[root_start:root_stop],
                child_result.normalized_calculations[child_start:child_stop],
            )
            or not np.array_equal(
                root_result.residuals[root_residual_start:root_residual_stop],
                child_result.residuals[child_residual_start:child_residual_stop],
            )
        ):
            return False
    return True


def _projection_mismatch(
    outcomes: tuple[FitComponentOutcome, ...],
) -> GroupedDirectTrfOutcome:
    return GroupedDirectTrfOutcome(
        GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
        outcomes,
        failure=TerminalFailure(
            "decomposition_projection_mismatch",
            "Fresh root objective differs from component evidence",
        ),
    )


def _validate_component_projections(
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    outcomes: tuple[FitComponentOutcome, ...],
    token: CancellationToken,
) -> tuple[_ComponentProjection, ...] | GroupedDirectTrfOutcome:
    projections: list[_ComponentProjection] = []
    try:
        for component, component_invocation, outcome in zip(
            decomposition.components,
            invocation.component_invocations,
            outcomes,
            strict=True,
        ):
            if token.is_cancelled:
                return GroupedDirectTrfOutcome(
                    GroupedDirectTrfTerminal.CANCELLED,
                    outcomes,
                )
            child_plan = outcome.evaluation_plan
            if child_plan is None:
                return GroupedDirectTrfOutcome(
                    GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
                    outcomes,
                    failure=TerminalFailure(
                        "decomposition_projection_mismatch",
                        "Successful component lacks its execution-time plan",
                    ),
                )
            projection = _materialized_component_projection(
                child_plan,
                component,
                component_invocation,
                outcome,
            )
            if token.is_cancelled:
                return GroupedDirectTrfOutcome(
                    GroupedDirectTrfTerminal.CANCELLED,
                    outcomes,
                )
            if isinstance(projection, TerminalFailure):
                return GroupedDirectTrfOutcome(
                    GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
                    outcomes,
                    failure=projection,
                )
            projections.append(projection)
    except KeyboardInterrupt:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.INTERRUPTED,
            outcomes,
        )
    except Exception as error:  # noqa: BLE001 - validation failures fail closed
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=TerminalFailure(
                "aggregate_validation_exception",
                f"{type(error).__name__}: {error}",
            ),
        )
    return tuple(projections)


def _component_gate(
    outcomes: tuple[FitComponentOutcome, ...],
    cancellation: CancellationToken | None,
) -> GroupedDirectTrfOutcome | None:
    dispositions = tuple(item.disposition for item in outcomes)
    if ComponentDisposition.INTERRUPTED in dispositions:
        return GroupedDirectTrfOutcome(GroupedDirectTrfTerminal.INTERRUPTED, outcomes)
    if ComponentDisposition.CANCELLED in dispositions or (
        cancellation is not None and cancellation.is_cancelled
    ):
        return GroupedDirectTrfOutcome(GroupedDirectTrfTerminal.CANCELLED, outcomes)
    if ComponentDisposition.EXECUTION_FAILURE in dispositions:
        superseding = next(
            (item.failure for item in outcomes if item.failure is not None),
            TerminalFailure(
                "component_execution_failure",
                "A component invalidated grouped execution",
            ),
        )
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.EXECUTION_FAILURE,
            outcomes,
            failure=superseding,
        )
    if all(item is ComponentDisposition.SUCCEEDED for item in dispositions):
        return None
    primary = next(
        (item.failure for item in outcomes if item.failure is not None),
        TerminalFailure("component_failure", "A required component failed"),
    )
    return GroupedDirectTrfOutcome(
        GroupedDirectTrfTerminal.COMPONENT_FAILURE,
        outcomes,
        failure=primary,
    )


def _aggregate_vector(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    outcomes: tuple[FitComponentOutcome, ...],
) -> tuple[float, ...] | TerminalFailure:
    assignments: dict[str, float] = {}
    for component, outcome in zip(decomposition.components, outcomes, strict=True):
        candidate = outcome.candidate
        if (
            outcome.controlled_ids != component.controlled_ids
            or candidate is None
            or len(candidate.vector) != len(component.controlled_ids)
        ):
            return TerminalFailure(
                "aggregate_assignment_failure",
                "A successful component lacks its exact controlled assignment",
            )
        for param_id, value in zip(
            component.controlled_ids,
            candidate.vector,
            strict=True,
        ):
            if param_id in assignments or param_id not in problem.controlled_ids:
                return TerminalFailure(
                    "aggregate_assignment_failure",
                    "Component assignments are duplicate or out of root scope",
                )
            assignments[param_id] = value
    if set(assignments) != set(problem.controlled_ids):
        return TerminalFailure(
            "aggregate_assignment_failure",
            "Component assignments do not cover the complete root vector",
        )
    vector = tuple(assignments[param_id] for param_id in problem.controlled_ids)
    if any(
        not lower <= value <= upper
        for value, lower, upper in zip(
            vector,
            problem.lower_bounds,
            problem.upper_bounds,
            strict=True,
        )
    ):
        return TerminalFailure(
            "aggregate_feasibility_failure",
            "Combined component assignment violates the root feasible domain",
        )
    return vector


def _compose_partitioned_final_jacobian(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    engine: EvaluationEngine,
    outcomes: tuple[FitComponentOutcome, ...],
    aggregate: EvaluationResult,
) -> FinalResidualJacobianEvidence | None:
    """Compose exact independent component Jacobians in canonical root order."""
    proof = decomposition.partition_proof
    if (
        not decomposition.components
        or proof.controlled_ids != problem.controlled_ids
        or proof.root_plan_identity != engine.plan.identity
        or len(outcomes) != len(decomposition.components)
    ):
        return None
    row_ranges: list[tuple[int, int]] = []
    offset = 0
    for profile in engine.plan.profiles:
        count = len(profile.retained_observation_indices)
        row_ranges.append((offset, offset + count))
        offset += count
    if offset != aggregate.residuals.size:
        return None
    root_columns = {
        param_id: index for index, param_id in enumerate(problem.controlled_ids)
    }
    matrix = np.zeros((offset, len(problem.controlled_ids)), dtype=np.float64)
    source_identities: list[str] = [proof.identity]
    for component, outcome in zip(
        decomposition.components,
        outcomes,
        strict=True,
    ):
        retained = outcome.final_residual_jacobian
        candidate = outcome.candidate
        if (
            retained is None
            or candidate is None
            or retained.controlled_ids != component.controlled_ids
            or retained.final_vector != candidate.vector
            or retained.final_residuals
            != tuple(float(value) for value in candidate.evaluation_result.residuals)
        ):
            return None
        retained_matrix = np.asarray(retained.matrix, dtype=np.float64)
        child_offset = 0
        for root_profile_index in component.root_profile_indices:
            root_start, root_stop = row_ranges[root_profile_index]
            count = root_stop - root_start
            child_stop = child_offset + count
            if child_stop > retained.shape[0]:
                return None
            for child_column, param_id in enumerate(component.controlled_ids):
                matrix[root_start:root_stop, root_columns[param_id]] = retained_matrix[
                    child_offset:child_stop,
                    child_column,
                ]
            child_offset = child_stop
        if child_offset != retained.shape[0]:
            return None
        source_identities.append(retained.identity)
    return FinalResidualJacobianEvidence(
        ResidualJacobianSource.FIT_PARTITION_COMPOSITION,
        problem.controlled_ids,
        tuple(
            aggregate.resolved_values[param_id] for param_id in problem.controlled_ids
        ),
        tuple(float(value) for value in aggregate.residuals),
        matrix.shape,
        tuple(tuple(float(value) for value in row) for row in matrix),
        tuple(source_identities),
    )


def _accept_grouped_aggregate_unchecked(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    engine: EvaluationEngine,
    outcomes: tuple[FitComponentOutcome, ...],
    projections: tuple[_ComponentProjection, ...],
    token: CancellationToken,
    materialized: MaterializedDirectTrfCandidate,
) -> GroupedDirectTrfOutcome:
    aggregate = materialized.evaluation_result
    for component, projection in zip(
        decomposition.components,
        projections,
        strict=True,
    ):
        if token.is_cancelled:
            return GroupedDirectTrfOutcome(
                GroupedDirectTrfTerminal.CANCELLED,
                outcomes,
            )
        if not _root_projection_matches(
            aggregate,
            engine,
            component,
            projection,
        ):
            return _projection_mismatch(outcomes)
        if token.is_cancelled:
            return GroupedDirectTrfOutcome(
                GroupedDirectTrfTerminal.CANCELLED,
                outcomes,
            )
    if token.is_cancelled:
        return GroupedDirectTrfOutcome(GroupedDirectTrfTerminal.CANCELLED, outcomes)
    final_residual_jacobian = _compose_partitioned_final_jacobian(
        problem,
        decomposition,
        engine,
        outcomes,
        aggregate,
    )
    if token.is_cancelled:
        return GroupedDirectTrfOutcome(GroupedDirectTrfTerminal.CANCELLED, outcomes)
    accepted, authority = accept_materialized_fit(
        problem=problem,
        invocation_identity=invocation.identity,
        execution_identity=materialized.execution_identity,
        materialization=materialized.materialization,
        vector=materialized.vector,
        chi_square=materialized.chi_square,
        evaluation_result=aggregate,
        final_residual_jacobian=final_residual_jacobian,
    )
    return GroupedDirectTrfOutcome(
        GroupedDirectTrfTerminal.ACCEPTED,
        outcomes,
        accepted,
        authority,
    )


def _accept_grouped_aggregate(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    engine: EvaluationEngine,
    outcomes: tuple[FitComponentOutcome, ...],
    projections: tuple[_ComponentProjection, ...],
    token: CancellationToken,
    materialized: MaterializedDirectTrfCandidate,
) -> GroupedDirectTrfOutcome:
    try:
        return _accept_grouped_aggregate_unchecked(
            problem,
            decomposition,
            invocation,
            engine,
            outcomes,
            projections,
            token,
            materialized,
        )
    except KeyboardInterrupt:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.INTERRUPTED,
            outcomes,
        )
    except Exception as error:  # noqa: BLE001 - validation failures fail closed
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=TerminalFailure(
                "aggregate_validation_exception",
                f"{type(error).__name__}: {error}",
            ),
        )


def _grouped_materialization_failure(
    failure: direct_trf_owner.RootMaterializationFailure,
    outcomes: tuple[FitComponentOutcome, ...],
) -> GroupedDirectTrfOutcome:
    if failure.terminal is MaterializationTerminal.CANCELLED:
        return GroupedDirectTrfOutcome(GroupedDirectTrfTerminal.CANCELLED, outcomes)
    if failure.terminal is MaterializationTerminal.INTERRUPTED:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.INTERRUPTED,
            outcomes,
        )
    detail = failure.failure
    if detail.category in {
        "materialization_binding_failure",
        "materialization_exception",
    }:
        detail = TerminalFailure(
            "aggregate_materialization_exception",
            detail.message,
            detail.evaluation_failure,
        )
    return GroupedDirectTrfOutcome(
        GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
        outcomes,
        failure=detail,
    )


def materialize_grouped_direct_trf(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    components: Sequence[FitComponentOutcome],
    *,
    cancellation: CancellationToken | None = None,
) -> GroupedDirectTrfOutcome:
    """Freshly validate one standalone aggregate before granting authority."""
    outcomes = tuple(components)
    token = CancellationToken() if cancellation is None else cancellation
    try:
        validation_identity = _validate_grouped_context(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            token,
        )
    except _GroupedValidationCancelled:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.CANCELLED,
            outcomes,
        )
    except KeyboardInterrupt:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.INTERRUPTED,
            outcomes,
        )
    except Exception as error:  # noqa: BLE001 - proof/context mismatch fails closed
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=TerminalFailure(
                "decomposition_context_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    return _materialize_grouped_direct_trf_validated(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        outcomes,
        validation_identity,
        cancellation=token,
    )


def _materialize_grouped_direct_trf_validated(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    components: Sequence[FitComponentOutcome],
    validation_identity: str,
    *,
    cancellation: CancellationToken | None = None,
) -> GroupedDirectTrfOutcome:
    """Materialize one aggregate carrying exact prior context validation."""
    outcomes = tuple(components)
    token = CancellationToken() if cancellation is None else cancellation
    if validation_identity != _grouped_context_validation_identity(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
    ):
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=TerminalFailure(
                "decomposition_context_failure",
                "Grouped context validation belongs to another invocation",
            ),
        )
    if tuple(item.component_identity for item in outcomes) != tuple(
        component.identity for component in decomposition.components
    ):
        raise DirectTrfConstructionError(
            "Component outcomes do not match canonical decomposition order"
        )
    gated = _component_gate(outcomes, token)
    if gated is not None:
        return gated
    projections_or_outcome = _validate_component_projections(
        decomposition,
        invocation,
        outcomes,
        token,
    )
    if isinstance(projections_or_outcome, GroupedDirectTrfOutcome):
        return projections_or_outcome
    projections = projections_or_outcome
    try:
        _check_grouped_cancellation(token)
        vector_or_failure = _aggregate_vector(problem, decomposition, outcomes)
        _check_grouped_cancellation(token)
    except _GroupedValidationCancelled:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.CANCELLED,
            outcomes,
        )
    if isinstance(vector_or_failure, TerminalFailure):
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=vector_or_failure,
        )
    vector = vector_or_failure
    materialized = direct_trf_owner.materialize_root_candidate(
        problem,
        parameterization,
        engine,
        vector=vector,
        invocation_identity=invocation.identity,
        execution_identity=lambda aggregate: _identity(
            "native-grouped-direct-trf-execution",
            (
                invocation.identity,
                tuple(outcome.identity for outcome in outcomes),
                tuple(float(value).hex() for value in vector),
                aggregate.identity,
            ),
        ),
        cancellation=token,
    )
    if isinstance(materialized, direct_trf_owner.RootMaterializationFailure):
        return _grouped_materialization_failure(materialized, outcomes)
    if isinstance(materialized, CandidateMaterialization):
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.DECOMPOSITION_VALIDATION_FAILURE,
            outcomes,
            failure=materialized.failure,
        )
    return _accept_grouped_aggregate(
        problem,
        decomposition,
        invocation,
        engine,
        outcomes,
        projections,
        token,
        materialized,
    )


def execute_grouped_direct_trf(
    problem: OptimizationProblem,
    decomposition: FitDecomposition,
    invocation: GroupedDirectTrfInvocation,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
    *,
    cancellation: CancellationToken | None = None,
    progress_observer: ContextualProgressObserver | None = None,
) -> GroupedDirectTrfOutcome:
    """Execute all exact components, then freshly validate the root aggregate."""
    token = CancellationToken() if cancellation is None else cancellation
    try:
        validation_identity = _validate_grouped_context(
            problem,
            decomposition,
            invocation,
            parameterization,
            engine,
            token,
        )
    except _GroupedValidationCancelled:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.CANCELLED,
            tuple(_not_started(component) for component in decomposition.components),
        )
    except KeyboardInterrupt:
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.INTERRUPTED,
            tuple(_not_started(component) for component in decomposition.components),
        )
    except Exception as error:  # noqa: BLE001 - invalidation gives every disposition
        return GroupedDirectTrfOutcome(
            GroupedDirectTrfTerminal.EXECUTION_FAILURE,
            tuple(_not_started(component) for component in decomposition.components),
            failure=TerminalFailure(
                "grouped_context_failure",
                f"{type(error).__name__}: {error}",
            ),
        )
    components = execute_direct_trf_components(
        decomposition,
        invocation,
        parameterization,
        engine,
        cancellation=token,
        progress_observer=progress_observer,
    )
    return _materialize_grouped_direct_trf_validated(
        problem,
        decomposition,
        invocation,
        parameterization,
        engine,
        components,
        validation_identity,
        cancellation=token,
    )
