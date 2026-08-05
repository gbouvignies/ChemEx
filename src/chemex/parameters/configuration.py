"""Immutable, sealed per-analysis parameter configuration.

Parameter TOML is applied once to `ParameterDefinitions` defaults using the
documented precedence rules (`chemex.configuration.parameters.DefaultSetting`),
producing an immutable per-analysis `ParameterConfiguration`. For each stable
ID it owns the effective initial value, the effective lower and upper bounds,
the supported parameter-file search metadata (`brute_step`), and field-level
override provenance.

Method roles, active constraints, and runtime `FIT`/`FIX`/`CONSTRAINTS`
selections are not configuration fields; they belong to a later, method-scoped
Active Parameterization built from the sealed configuration.

`build_parameter_configuration` is the single one-way legacy adapter for this
layer: the model-free auxiliary-initialization precedence and constraint
evaluation used to resolve numeric values remain the legacy-authoritative
`ParameterCatalog`'s responsibility (see `chemex.parameters.database`); this
adapter only seals the already-resolved effective values into an immutable,
deterministically ordered native snapshot, while independently tracking
parameter-file override provenance.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, replace

from chemex.configuration.parameters import DefaultListType
from chemex.parameters.database import ParameterCatalog, ParameterIndex
from chemex.parameters.definitions import (
    ParameterDefinitions,
    SealedParameterModelError,
)

__all__ = [
    "ParameterConfiguration",
    "ParameterConfigurationBuilder",
    "ParameterConfigurationEntry",
    "build_parameter_configuration",
]


@dataclass(frozen=True, slots=True)
class ParameterConfigurationEntry:
    """The immutable, per-analysis configuration of a single parameter.

    Attributes:
        param_id: The stable ID (`ParamName.id_`) identifying the parameter.
        value: The effective initial value.
        min: The effective lower bound.
        max: The effective upper bound.
        brute_step: The effective grid/brute-force step, if any.
        overridden: Whether a parameter-file default entry matched and
            overrode this parameter's definition defaults.

    """

    param_id: str
    value: float | None
    min: float
    max: float
    brute_step: float | None
    overridden: bool


class ParameterConfiguration(Mapping[str, ParameterConfigurationEntry]):
    """A sealed, deterministically ordered set of parameter configurations.

    Instances are produced exclusively by
    `ParameterConfigurationBuilder.seal` and are immutable: iteration order is
    the canonical ascending order of stable parameter IDs.
    """

    def __init__(self, entries: Mapping[str, ParameterConfigurationEntry]) -> None:
        self._entries = dict(sorted(entries.items()))

    def __setattr__(self, name: str, value: object) -> None:
        if "_entries" in self.__dict__:
            msg = f"Cannot set {name!r}: ParameterConfiguration is sealed"
            raise SealedParameterModelError(msg)
        super().__setattr__(name, value)

    def __getitem__(self, param_id: str) -> ParameterConfigurationEntry:
        return self._entries[param_id]

    def __iter__(self) -> Iterator[str]:
        return iter(self._entries)

    def __len__(self) -> int:
        return len(self._entries)

    def __repr__(self) -> str:
        return f"ParameterConfiguration({len(self)} parameters)"

    @property
    def ids(self) -> tuple[str, ...]:
        """The stable parameter IDs in canonical ascending order."""
        return tuple(self._entries)


class ParameterConfigurationBuilder:
    """Builds an immutable `ParameterConfiguration` from sealed definitions.

    Initial entries mirror each definition's default value/bounds. Applying
    parameter-file defaults follows the exact precedence used by the
    legacy-authoritative `ParameterCatalog.set_defaults`: later matching
    entries never override an ID already claimed by an earlier (higher
    priority) one, and an explicit `min`/`max` override also carries the
    matching `brute_step`.
    """

    def __init__(self, definitions: ParameterDefinitions) -> None:
        self._index = ParameterIndex()
        for definition in definitions.values():
            self._index.add(definition.param_name)
        self._entries: dict[str, ParameterConfigurationEntry] = {
            param_id: ParameterConfigurationEntry(
                param_id=param_id,
                value=definition.value,
                min=definition.min,
                max=definition.max,
                brute_step=None,
                overridden=False,
            )
            for param_id, definition in definitions.items()
        }
        self._sealed = False

    def apply_defaults(self, defaults: DefaultListType) -> None:
        """Apply parameter-file defaults using the legacy matching precedence."""
        if self._sealed:
            msg = "Cannot apply defaults: ParameterConfigurationBuilder is sealed"
            raise SealedParameterModelError(msg)

        id_pool = set(self._entries)
        for name_to_set, setting in reversed(defaults):
            matching_ids = self._index.get_matching_ids(name_to_set) & id_pool
            id_pool -= matching_ids
            for matching_id in matching_ids:
                current = self._entries[matching_id]
                self._entries[matching_id] = ParameterConfigurationEntry(
                    param_id=matching_id,
                    value=setting.value,
                    min=setting.min if setting.min is not None else current.min,
                    max=setting.max if setting.max is not None else current.max,
                    brute_step=(
                        setting.brute_step
                        if setting.min is not None
                        else current.brute_step
                    ),
                    overridden=True,
                )

    def overwrite_values(self, values: Mapping[str, float | None]) -> None:
        """Overwrite effective initial values without touching bounds/provenance.

        Used by `build_parameter_configuration` to seal in the exact numeric
        values already resolved by the legacy-authoritative parameter catalog
        (including cross-catalog model-free auxiliary initialization), without
        reimplementing constraint evaluation natively.
        """
        if self._sealed:
            msg = "Cannot overwrite values: ParameterConfigurationBuilder is sealed"
            raise SealedParameterModelError(msg)

        for param_id, value in values.items():
            if value is None or param_id not in self._entries:
                continue
            self._entries[param_id] = replace(self._entries[param_id], value=value)

    def seal(self) -> ParameterConfiguration:
        """Seal the collected configuration into an immutable `ParameterConfiguration`."""
        if self._sealed:
            msg = "ParameterConfigurationBuilder is already sealed"
            raise SealedParameterModelError(msg)
        self._sealed = True
        return ParameterConfiguration(self._entries)


def build_parameter_configuration(
    definitions: ParameterDefinitions,
    catalog: ParameterCatalog,
    defaults: DefaultListType,
) -> ParameterConfiguration:
    """Adapt a legacy-authoritative catalog into a sealed `ParameterConfiguration`.

    Args:
        definitions: The sealed native definitions for the same catalog
            (ordinary or model-free).
        catalog: The legacy-authoritative catalog, already defaulted via
            `ParameterCatalog.set_defaults`/`ParameterStore.set_defaults`.
        defaults: The same parameter-file defaults passed to
            `ParameterStore.set_defaults`, used to independently compute
            field-level override provenance and effective bounds.

    Returns:
        A sealed `ParameterConfiguration` whose effective values are read
        directly from `catalog` (so they exactly match, by construction, the
        legacy-authoritative fitted/simulated values), while bounds,
        brute-force step, and override provenance are computed natively.

    """
    builder = ParameterConfigurationBuilder(definitions)
    builder.apply_defaults(defaults)

    resolved = catalog.get_parameters(definitions.ids)
    builder.overwrite_values(
        {param_id: setting.value for param_id, setting in resolved.items()},
    )

    return builder.seal()
