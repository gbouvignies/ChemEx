"""Canonical, sealed ChemEx parameter definitions.

A stable parameter ID identifies exactly one canonical, immutable
:class:`ParameterDefinition`: scientific identity and name, spin-system and
condition scope, base value and domain, and only the default constraint
metadata that genuinely belongs to the definition (for example the built-in
proton-exchange or model-free relations produced by
:mod:`chemex.parameters.spins` and :mod:`chemex.models.factory`).

Method roles, active constraints selected by ``FIT``/``FIX``/``CONSTRAINTS``,
current values, optimizer metadata, estimates, and uncertainty are not
definition fields and are never compared for conflicts here; they belong to a
later, method-scoped Active Parameterization.

Every experiment/profile contribution for a given stable ID must agree on the
canonicalized definition. There is no first-contribution-wins behavior:
:class:`ParameterDefinitionsBuilder` raises :class:`ParameterDefinitionConflictError`
the moment two contributions disagree, and :meth:`ParameterDefinitionsBuilder.seal`
returns an immutable, deterministically ordered :class:`ParameterDefinitions`
that cannot be mutated afterward.

:class:`ParameterDefinitionsCollector` is the single one-way adapter that
observes the legacy-authoritative :class:`chemex.parameters.factory.ParameterFactory`
while it constructs experiments, without changing its behavior: attaching a
collector is opt-in and has no effect on the legacy parameter catalogs.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from chemex.parameters.name import ParamName

if TYPE_CHECKING:
    from chemex.parameters.setting import Parameters, ParamSetting

__all__ = [
    "ParameterDefinition",
    "ParameterDefinitionConflictError",
    "ParameterDefinitions",
    "ParameterDefinitionsBuilder",
    "ParameterDefinitionsCollector",
    "SealedParameterModelError",
]

_COMPARABLE_FIELDS = ("value", "min", "max", "vary", "expr")


class SealedParameterModelError(RuntimeError):
    """Raised when a sealed parameter aggregate is mutated after sealing."""


class ParameterDefinitionConflictError(ValueError):
    """Raised when two contributions disagree on a parameter's canonical definition.

    Attributes:
        param_id: The stable parameter ID in conflict.
        field_name: The first canonical field found to disagree.
        values: The two disagreeing values, `(existing, incoming)`.
        contributors: Labels identifying the two contributors, `(existing, incoming)`.

    """

    def __init__(
        self,
        param_id: str,
        field_name: str,
        values: tuple[object, object],
        contributors: tuple[str, str],
    ) -> None:
        msg = (
            f"Conflicting {field_name!r} for parameter {param_id!r}: "
            f"{contributors[0]!r} defines {values[0]!r}, "
            f"{contributors[1]!r} defines {values[1]!r}"
        )
        super().__init__(msg)
        self.param_id = param_id
        self.field_name = field_name
        self.values = values
        self.contributors = contributors


@dataclass(frozen=True, slots=True)
class ParameterDefinition:
    """The immutable scientific identity of a single canonical parameter.

    Attributes:
        param_id: The stable ID (`ParamName.id_`) identifying the parameter.
        param_name: The canonical scientific name, spin-system, and
            condition scope.
        value: The default scientific value from the model/spin-system
            construction.
        min: The default lower bound.
        max: The default upper bound.
        vary: The default vary status coming from the model/spin-system
            construction, not a runtime `FIT`/`FIX` selection.
        expr: The default constraint expression coming from the
            model/spin-system construction, not a runtime `CONSTRAINTS`
            selection.

    """

    param_id: str
    param_name: ParamName
    value: float | None
    min: float
    max: float
    vary: bool
    expr: str

    @classmethod
    def from_setting(cls, param_id: str, setting: ParamSetting) -> ParameterDefinition:
        """Build a `ParameterDefinition` snapshot from a legacy `ParamSetting`."""
        return cls(
            param_id=param_id,
            param_name=setting.param_name,
            value=setting.value,
            min=setting.min,
            max=setting.max,
            vary=setting.vary,
            expr=setting.expr,
        )

    def _canonical_values(self) -> tuple[object, ...]:
        return tuple(getattr(self, name) for name in _COMPARABLE_FIELDS)


def _first_differing_field(
    existing: ParameterDefinition,
    incoming: ParameterDefinition,
) -> str:
    for name in _COMPARABLE_FIELDS:
        if getattr(existing, name) != getattr(incoming, name):
            return name
    return _COMPARABLE_FIELDS[-1]


class ParameterDefinitions(Mapping[str, ParameterDefinition]):
    """A sealed, deterministically ordered set of parameter definitions.

    Instances are produced exclusively by :meth:`ParameterDefinitionsBuilder.seal`
    and are immutable: iteration order is the canonical ascending order of
    stable parameter IDs, independent of contribution order.
    """

    def __init__(self, definitions: Mapping[str, ParameterDefinition]) -> None:
        self._definitions = dict(sorted(definitions.items()))

    def __setattr__(self, name: str, value: object) -> None:
        if "_definitions" in self.__dict__:
            msg = f"Cannot set {name!r}: ParameterDefinitions is sealed"
            raise SealedParameterModelError(msg)
        super().__setattr__(name, value)

    def __getitem__(self, param_id: str) -> ParameterDefinition:
        return self._definitions[param_id]

    def __iter__(self) -> Iterator[str]:
        return iter(self._definitions)

    def __len__(self) -> int:
        return len(self._definitions)

    def __repr__(self) -> str:
        return f"ParameterDefinitions({len(self)} parameters)"

    @property
    def ids(self) -> tuple[str, ...]:
        """The stable parameter IDs in canonical ascending order."""
        return tuple(self._definitions)


class ParameterDefinitionsBuilder:
    """Collects parameter-definition contributions and seals them once complete.

    Every contribution for a given stable ID must canonically agree with any
    prior contribution for that ID; a disagreement raises
    `ParameterDefinitionConflictError` immediately rather than silently keeping
    the first contribution. Once sealed, the builder can no longer be mutated.
    """

    def __init__(self) -> None:
        self._definitions: dict[str, ParameterDefinition] = {}
        self._contributors: dict[str, str] = {}
        self._sealed = False

    def add(self, param_id: str, setting: ParamSetting, *, source: str) -> None:
        """Add one parameter contribution, validating it against any existing one."""
        if self._sealed:
            msg = "Cannot add contributions: ParameterDefinitionsBuilder is sealed"
            raise SealedParameterModelError(msg)

        definition = ParameterDefinition.from_setting(param_id, setting)
        existing = self._definitions.get(param_id)
        if existing is None:
            self._definitions[param_id] = definition
            self._contributors[param_id] = source
            return

        if existing._canonical_values() != definition._canonical_values():
            field_name = _first_differing_field(existing, definition)
            raise ParameterDefinitionConflictError(
                param_id,
                field_name,
                (getattr(existing, field_name), getattr(definition, field_name)),
                (self._contributors[param_id], source),
            )

    def add_many(self, parameters: Parameters, *, source: str) -> None:
        """Add every parameter in `parameters`, validating each against prior ones."""
        for param_id, setting in parameters.items():
            self.add(param_id, setting, source=source)

    def seal(self) -> ParameterDefinitions:
        """Seal the collected contributions into an immutable `ParameterDefinitions`."""
        if self._sealed:
            msg = "ParameterDefinitionsBuilder is already sealed"
            raise SealedParameterModelError(msg)
        self._sealed = True
        return ParameterDefinitions(self._definitions)


@dataclass
class ParameterDefinitionsCollector:
    """Observes the legacy-authoritative `ParameterFactory` without changing it.

    Attach a fresh collector to `ParameterFactory.definitions` to build sealed
    native `ParameterDefinitions` for the ordinary and model-free catalogs
    while the legacy factory constructs experiments in the usual way. The
    legacy catalogs remain authoritative; this collector only mirrors their
    per-ID scientific definitions into an independently validated, sealed
    snapshot.
    """

    ordinary: ParameterDefinitionsBuilder = field(
        default_factory=ParameterDefinitionsBuilder,
    )
    model_free: ParameterDefinitionsBuilder = field(
        default_factory=ParameterDefinitionsBuilder,
    )

    def collect(
        self,
        parameters: Parameters,
        parameters_mf: Parameters,
        *,
        source: str,
    ) -> None:
        """Record one profile's ordinary and model-free contributions."""
        self.ordinary.add_many(parameters, source=source)
        self.model_free.add_many(parameters_mf, source=source)

    def seal(self) -> tuple[ParameterDefinitions, ParameterDefinitions]:
        """Seal both catalogs, returning `(ordinary, model_free)` definitions."""
        return self.ordinary.seal(), self.model_free.seal()
