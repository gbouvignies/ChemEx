"""Compile Method language into value-independent fitting instructions."""

from __future__ import annotations

from dataclasses import dataclass, field

from chemex.configuration.method_plan import (
    GridSearch,
    MethodPlan,
    ProfilePreparation,
    ProfileSelection,
    StatisticsPlan,
)
from chemex.configuration.method_validation import (
    ResolvedGridAxis,
    project_grid_axes,
    project_preparation,
    resolve_method_plan,
)
from chemex.containers.experiment import Experiment, project_profile_selection
from chemex.containers.experiments import Experiments
from chemex.containers.profile import Profile
from chemex.parameters.parameterization import (
    SealedParameterModel,
    StaticParameterization,
    compile_static_parameterization,
)


@dataclass(frozen=True, slots=True)
class ProfileBinding:
    experiment: Experiment = field(compare=False, repr=False)
    selected: tuple[Profile, ...] = field(compare=False, repr=False)
    filtered: tuple[Profile, ...] = field(compare=False, repr=False)

    def activate(self) -> None:
        """Apply compiled membership without evaluating a selection."""
        self.experiment.profiles = list(self.selected)
        self.experiment.filtered_profiles = list(self.filtered)


@dataclass(frozen=True, slots=True)
class GridSearchInstruction:
    axes: tuple[ResolvedGridAxis, ...]


@dataclass(frozen=True, slots=True)
class ProfilePreparationInstruction:
    hold_ids: tuple[str, ...]
    mirror_ids: tuple[str, ...]


type SearchInstruction = GridSearchInstruction | ProfilePreparationInstruction | None


@dataclass(frozen=True, slots=True)
class FitStep:
    name: str
    ordinal: int
    bindings: tuple[ProfileBinding, ...]
    profile_count: int
    parameterization: StaticParameterization
    search: SearchInstruction
    statistics: StatisticsPlan | None


@dataclass(frozen=True, slots=True)
class SkippedStep:
    name: str
    ordinal: int
    reason: str = "no_profiles"


type ExecutableStep = FitStep | SkippedStep


@dataclass(frozen=True, slots=True)
class ExecutableMethodPlan:
    steps: tuple[ExecutableStep, ...]


def _select_profiles(
    previous: tuple[ProfileBinding, ...], selection: ProfileSelection
) -> tuple[ProfileBinding, ...]:
    selected_bindings: list[ProfileBinding] = []
    for binding in previous:
        selected, filtered = project_profile_selection(
            binding.selected, binding.filtered, selection
        )
        selected_bindings.append(ProfileBinding(binding.experiment, selected, filtered))
    return tuple(selected_bindings)


def compile_method_plan(
    plan: MethodPlan,
    model: SealedParameterModel,
    experiments: Experiments,
) -> ExecutableMethodPlan:
    """Resolve every Method Step before any numerical value or fit is read."""
    resolved = resolve_method_plan(plan, model)
    bindings = tuple(
        ProfileBinding(
            experiment,
            tuple(experiment.profiles),
            tuple(experiment.filtered_profiles),
        )
        for experiment in experiments
    )
    steps: list[ExecutableStep] = []
    for ordinal, (step, meaning) in enumerate(
        zip(plan.steps, resolved, strict=True), 1
    ):
        bindings = _select_profiles(bindings, step.selection)
        selected = tuple(
            profile for binding in bindings for profile in binding.selected
        )
        if not selected:
            steps.append(SkippedStep(step.name, ordinal))
            continue
        required_ids = set().union(*(profile.param_ids for profile in selected))
        parameterization = compile_static_parameterization(
            model, meaning.roles, meaning.constraints, required_ids
        )
        search: SearchInstruction = None
        if isinstance(step.search, GridSearch):
            axes = project_grid_axes(
                meaning.grid_axes,
                model,
                active_scope_ids=parameterization.program.scope_ids,
                final_fit_ids=parameterization.fit_ids,
            )
            search = GridSearchInstruction(axes)
        elif isinstance(step.search, ProfilePreparation):
            hold, mirror = project_preparation(
                step.search, meaning.preparation, model, parameterization.fit_ids
            )
            search = ProfilePreparationInstruction(hold, mirror)
        steps.append(
            FitStep(
                step.name,
                ordinal,
                bindings,
                len(selected),
                parameterization,
                search,
                step.statistics,
            )
        )
    return ExecutableMethodPlan(tuple(steps))
