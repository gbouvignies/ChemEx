"""Registries for functions used by scientific constraint expressions."""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, ClassVar, Literal

from chemex.configuration.conditions import Conditions
from chemex.parameters.setting import ParamLocalSetting

SettingsType = dict[str, ParamLocalSetting]
SettingMakerType = Callable[[Conditions], SettingsType]

type FunctionArgumentDomain = Literal["unbounded", "nonnegative", "positive"]
_MINIMUM_POSITIVE_FLOAT = math.ulp(0.0)


@dataclass(frozen=True, slots=True)
class NumericalFunctionLinearization:
    """Model-owned numerical derivative policy for one function component."""

    function_id: str
    component: str | None
    argument_scales: tuple[float, ...]
    argument_domains: tuple[FunctionArgumentDomain, ...]
    output_scale: float = 1.0
    normalized_population_components: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class AnalyticFunctionLinearization:
    """Model-owned analytic derivatives for one scientific function component."""

    function_id: str
    component: str | None
    implementation_identity: str
    partials: tuple[Callable[..., float], ...]


type FunctionLinearization = (
    NumericalFunctionLinearization | AnalyticFunctionLinearization
)


def population_linearizations(
    argument_scales: tuple[float, ...],
    argument_domains: tuple[FunctionArgumentDomain, ...],
    *components: str,
) -> tuple[NumericalFunctionLinearization, ...]:
    """Build coupled derivative policies for one normalized population vector."""
    population_components = tuple(components)
    if len(population_components) < 2 or len(set(population_components)) != len(
        population_components
    ):
        raise ValueError("Normalized populations require distinct components")
    # Let each population's sampled magnitude set its roundoff scale.  A unit
    # floor would overwhelm informative subnormal components.
    return tuple(
        NumericalFunctionLinearization(
            function_id="populations",
            component=component,
            argument_scales=argument_scales,
            argument_domains=argument_domains,
            output_scale=_MINIMUM_POSITIVE_FLOAT,
            normalized_population_components=population_components,
        )
        for component in population_components
    )


class Registry:
    """Registry for storing functions used by parameter constraints."""

    user_function_registry: ClassVar[dict[str, Any]] = {}

    def register(self, name: str, user_functions: dict[str, Any]) -> None:
        """Register a new set of user functions."""
        self.user_function_registry[name] = user_functions

    def get(self, name: str) -> dict[str, Any]:
        if name not in self.user_function_registry:
            return {}
        return self.user_function_registry[name]


user_function_registry = Registry()


class LinearizationRegistry:
    """Registry for explicitly approved model-owned function derivatives."""

    _items: ClassVar[
        dict[
            str,
            tuple[FunctionLinearization, ...],
        ]
    ] = {}

    def register(
        self,
        name: str,
        capabilities: tuple[FunctionLinearization, ...],
    ) -> None:
        self._items[name] = capabilities

    def get(
        self,
        name: str,
    ) -> tuple[FunctionLinearization, ...]:
        return self._items.get(name, ())


function_linearization_registry = LinearizationRegistry()
