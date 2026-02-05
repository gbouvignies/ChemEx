"""Registry for B1 distribution plugins."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, ClassVar

from chemex.nmr.constants import Distribution

# Type alias for distribution generator functions
DistributionGenerator = Callable[..., Distribution]


class DistributionRegistry:
    """Registry for B1 distribution generators."""

    _generators: ClassVar[dict[str, DistributionGenerator]] = {}
    _configs: ClassVar[dict[str, type[Any]]] = {}

    @classmethod
    def register(
        cls,
        name: str,
        generator: DistributionGenerator,
        config_class: type[Any],
    ) -> None:
        """Register a distribution generator and its config class.

        Parameters
        ----------
        name : str
            Distribution type name (e.g., "gaussian", "hermite")
        generator : DistributionGenerator
            Function that generates the distribution
        config_class : type
            Pydantic model class for configuration

        """
        cls._generators[name] = generator
        cls._configs[name] = config_class

    @classmethod
    def get_generator(cls, name: str) -> DistributionGenerator:
        """Get a distribution generator by name.

        Parameters
        ----------
        name : str
            Distribution type name

        Returns
        -------
        DistributionGenerator
            The generator function

        Raises
        ------
        ValueError
            If distribution type is not registered

        """
        try:
            return cls._generators[name]
        except KeyError:
            valid_types = ", ".join(cls._generators.keys())
            msg = f"Unknown B1 distribution type: {name}. Valid types: {valid_types}"
            raise ValueError(msg) from None

    @classmethod
    def get_all_generators(cls) -> dict[str, DistributionGenerator]:
        """Get all registered generators."""
        return cls._generators.copy()

    @classmethod
    def get_all_config_classes(cls) -> list[type[Any]]:
        """Get all registered config classes for discriminated union."""
        return list(cls._configs.values())

    @classmethod
    def get_config_class(cls, name: str) -> type[Any]:
        """Get a registered config class by distribution name.

        Parameters
        ----------
        name : str
            Distribution type name (e.g., "gaussian")

        Returns
        -------
        type[Any]
            The Pydantic configuration class

        Raises
        ------
        KeyError
            If no config class is registered under that name

        """
        return cls._configs[name]


# Singleton instance
registry = DistributionRegistry()


def get_b1_distribution(
    distribution_type: str,
    value: float,
    **kwargs: float | str,
) -> Distribution:
    """Get a B1 distribution using the registry.

    Parameters
    ----------
    distribution_type : str
        Type of distribution from the registry
    value : float
        Nominal B1 field value
    **kwargs
        Additional parameters for the distribution function

    Returns
    -------
    Distribution
        B1 values and weights for integration

    Raises
    ------
    ValueError
        If distribution_type is not in registry

    """
    generator = registry.get_generator(distribution_type)
    return generator(value, **kwargs)
