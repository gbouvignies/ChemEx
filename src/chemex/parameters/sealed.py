"""Canonical immutable parameter definitions and per-analysis configuration.

These ChemEx-native objects are constructed beside the legacy mutable catalog.
The legacy execution path remains authoritative at this migration checkpoint.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import overload

from chemex.configuration.conditions import Conditions
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem

_CONDITION_FIELDS = ("d2o", "h_larmor_frq", "l_total", "p_total", "temperature")
_DEFINITION_FIELDS = (
    "name",
    "spin_system_name",
    "condition_entries",
    "default_value",
    "lower_bound",
    "upper_bound",
)


def extract_condition_entries(
    conditions: Conditions,
) -> tuple[tuple[str, float], ...]:
    """Extract scoped numeric condition fields into an immutable tuple."""
    return tuple(
        (field_name, value)
        for field_name in _CONDITION_FIELDS
        if (value := getattr(conditions, field_name, None)) is not None
    )


@dataclass(frozen=True, slots=True)
class ParamDefinition:
    """Immutable scientific identity and physical metadata of a parameter."""

    param_id: str
    name: str
    spin_system_name: str
    condition_entries: tuple[tuple[str, float], ...]
    default_value: float | None
    lower_bound: float
    upper_bound: float

    def __post_init__(self) -> None:
        entries = tuple(
            (str(name), float(value)) for name, value in self.condition_entries
        )
        object.__setattr__(self, "condition_entries", entries)


@dataclass(frozen=True, slots=True)
class DefinitionContribution:
    """One definition plus stable context identifying its construction source."""

    definition: ParamDefinition
    contributor: str


@dataclass(frozen=True, slots=True)
class DefinitionFieldConflict:
    """All contributed values for one disagreeing definition field."""

    field: str
    values: tuple[object, ...]


class DefinitionConflictError(ValueError):
    """Raised when contributions for one stable ID disagree."""

    def __init__(
        self,
        param_id: str,
        conflicts: Sequence[DefinitionFieldConflict],
        contributors: Sequence[str],
    ) -> None:
        self.param_id = param_id
        self.conflicts = tuple(conflicts)
        self.conflicting_fields = tuple(conflict.field for conflict in self.conflicts)
        self.contributors = tuple(contributors)
        # Retain the first-field convenience attributes used by early #582 callers.
        self.field = self.conflicts[0].field
        self.values = self.conflicts[0].values
        fields = "; ".join(
            f"{conflict.field}={conflict.values!r}" for conflict in self.conflicts
        )
        sources = ", ".join(repr(contributor) for contributor in self.contributors)
        super().__init__(
            f"Conflicting definitions for parameter {param_id!r}: {fields}; "
            f"contributors=({sources})"
        )


@dataclass(frozen=True, slots=True)
class ParamConfig:
    """Immutable per-analysis effective state after construction precedence."""

    param_id: str
    effective_value: float | None
    lower_bound: float
    upper_bound: float


def _float_token(value: float | None) -> str | None:
    return None if value is None else float(value).hex()


def _fingerprint(kind: str, records: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "schema": 1, "records": records},
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _definition_identity(definitions: Sequence[ParamDefinition]) -> str:
    records = [
        (
            definition.param_id,
            definition.name,
            definition.spin_system_name,
            tuple(
                (name, _float_token(value))
                for name, value in definition.condition_entries
            ),
            _float_token(definition.default_value),
            _float_token(definition.lower_bound),
            _float_token(definition.upper_bound),
        )
        for definition in definitions
    ]
    return _fingerprint("parameter-definitions", records)


def _configuration_identity(configs: Sequence[ParamConfig]) -> str:
    records = [
        (
            config.param_id,
            _float_token(config.effective_value),
            _float_token(config.lower_bound),
            _float_token(config.upper_bound),
        )
        for config in configs
    ]
    return _fingerprint("parameter-configuration", records)


@dataclass(frozen=True, slots=True)
class SealedDefinitions:
    """Deeply immutable, canonically ordered parameter definitions."""

    _definitions: tuple[ParamDefinition, ...]
    _index: Mapping[str, int]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        definitions = tuple(self._definitions)
        index_data = {
            definition.param_id: position
            for position, definition in enumerate(definitions)
        }
        if len(index_data) != len(definitions):
            msg = "Sealed definitions contain duplicate parameter IDs"
            raise ValueError(msg)
        index = MappingProxyType(index_data)
        object.__setattr__(self, "_definitions", definitions)
        object.__setattr__(self, "_index", index)
        object.__setattr__(self, "identity", _definition_identity(definitions))

    def __len__(self) -> int:
        return len(self._definitions)

    def __iter__(self) -> Iterator[ParamDefinition]:
        return iter(self._definitions)

    def __contains__(self, param_id: object) -> bool:
        return param_id in self._index

    @overload
    def __getitem__(self, key: str) -> ParamDefinition: ...

    @overload
    def __getitem__(self, key: int) -> ParamDefinition: ...

    def __getitem__(self, key: str | int) -> ParamDefinition:
        if isinstance(key, int):
            return self._definitions[key]
        return self._definitions[self._index[key]]


@dataclass(frozen=True, slots=True)
class SealedConfiguration:
    """Deeply immutable configuration in definition order."""

    _configs: tuple[ParamConfig, ...]
    _index: Mapping[str, int]
    definitions_identity: str = ""
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        configs = tuple(self._configs)
        index_data = {
            config.param_id: position for position, config in enumerate(configs)
        }
        if len(index_data) != len(configs):
            msg = "Sealed configuration contains duplicate parameter IDs"
            raise ValueError(msg)
        index = MappingProxyType(index_data)
        object.__setattr__(self, "_configs", configs)
        object.__setattr__(self, "_index", index)
        object.__setattr__(
            self,
            "identity",
            _fingerprint(
                "parameter-configuration",
                (self.definitions_identity, _configuration_identity(configs)),
            ),
        )

    def __len__(self) -> int:
        return len(self._configs)

    def __iter__(self) -> Iterator[ParamConfig]:
        return iter(self._configs)

    def __contains__(self, param_id: object) -> bool:
        return param_id in self._index

    @overload
    def __getitem__(self, key: str) -> ParamConfig: ...

    @overload
    def __getitem__(self, key: int) -> ParamConfig: ...

    def __getitem__(self, key: str | int) -> ParamConfig:
        if isinstance(key, int):
            return self._configs[key]
        return self._configs[self._index[key]]


type ContributionLike = (
    DefinitionContribution | ParamDefinition | tuple[ParamDefinition, str]
)


def _coerce_contribution(contribution: ContributionLike) -> DefinitionContribution:
    if isinstance(contribution, DefinitionContribution):
        return contribution
    if isinstance(contribution, ParamDefinition):
        return DefinitionContribution(contribution, "unspecified contributor")
    definition, contributor = contribution
    return DefinitionContribution(definition, contributor)


def _equal(left: object, right: object) -> bool:
    if (
        isinstance(left, float)
        and isinstance(right, float)
        and math.isnan(left)
        and math.isnan(right)
    ):
        return True
    return left == right


def _definition_sort_key(definition: ParamDefinition) -> ParamName:
    """Reproduce the legacy-authoritative ``ParamName`` ordering exactly."""
    return ParamName(
        definition.name,
        SpinSystem.from_name(definition.spin_system_name),
        Conditions.model_construct(None, **dict(definition.condition_entries)),
    )


class InvalidConstructionError(ValueError):
    """Base error for invalid initial values or physical bounds."""

    artifact = "parameter state"

    def __init__(
        self,
        param_id: str,
        field: str,
        value: object,
        detail: str,
    ) -> None:
        self.param_id = param_id
        self.field = field
        self.value = value
        self.detail = detail
        super().__init__(
            f"Invalid {self.artifact} for parameter {param_id!r}: "
            f"{field}={value!r}; {detail}"
        )


class InvalidDefinitionError(InvalidConstructionError):
    """Raised when a definition default or physical bound cannot be sealed."""

    artifact = "definition"


class InvalidConfigurationError(InvalidConstructionError):
    """Raised when an effective value or physical bound cannot be sealed."""

    artifact = "configuration"


def _validate_initial_state(
    param_id: str,
    value: float | None,
    lower_bound: float,
    upper_bound: float,
    *,
    value_field: str,
    error_type: type[InvalidConstructionError],
) -> None:
    for field_name, field_value in (
        (value_field, value),
        ("lower_bound", lower_bound),
        ("upper_bound", upper_bound),
    ):
        if field_value is not None and math.isnan(field_value):
            raise error_type(
                param_id,
                field_name,
                field_value,
                "NaN is not a valid configured value",
            )
    if lower_bound > upper_bound:
        raise error_type(
            param_id,
            "lower_bound",
            lower_bound,
            f"must not exceed upper_bound={upper_bound!r}",
        )
    if value is not None and not lower_bound <= value <= upper_bound:
        raise error_type(
            param_id,
            value_field,
            value,
            (
                "must be within configured physical bounds "
                f"[{lower_bound!r}, {upper_bound!r}]"
            ),
        )


def canonicalize_definitions(
    contributions: Mapping[str, Sequence[ContributionLike]],
) -> SealedDefinitions:
    """Validate contributions and return one canonical definition per stable ID."""
    canonical: list[ParamDefinition] = []
    for param_id, raw_contributions in contributions.items():
        normalized = tuple(
            sorted(
                (_coerce_contribution(item) for item in raw_contributions),
                key=lambda item: (item.contributor, repr(item.definition)),
            )
        )
        if not normalized:
            continue
        if any(item.definition.param_id != param_id for item in normalized):
            msg = f"Contribution key {param_id!r} does not match its parameter ID"
            raise ValueError(msg)

        reference = normalized[0].definition
        conflicts = tuple(
            DefinitionFieldConflict(
                field_name,
                tuple(getattr(item.definition, field_name) for item in normalized),
            )
            for field_name in _DEFINITION_FIELDS
            if any(
                not _equal(
                    getattr(reference, field_name),
                    getattr(item.definition, field_name),
                )
                for item in normalized[1:]
            )
        )
        if conflicts:
            raise DefinitionConflictError(
                param_id,
                conflicts,
                tuple(item.contributor for item in normalized),
            )
        _validate_initial_state(
            param_id,
            reference.default_value,
            reference.lower_bound,
            reference.upper_bound,
            value_field="default_value",
            error_type=InvalidDefinitionError,
        )
        canonical.append(reference)

    definitions = tuple(sorted(canonical, key=_definition_sort_key))
    index = {
        definition.param_id: position for position, definition in enumerate(definitions)
    }
    return SealedDefinitions(_definitions=definitions, _index=index)


class ConfigurationMismatchError(ValueError):
    """Raised when definitions and the active legacy catalog do not match."""

    def __init__(self, missing: Sequence[str], unexpected: Sequence[str]) -> None:
        self.missing = tuple(missing)
        self.unexpected = tuple(unexpected)
        super().__init__(
            "Parameter configuration does not match sealed definitions: "
            f"missing={self.missing!r}, unexpected={self.unexpected!r}"
        )


def build_sealed_configuration(
    definitions: SealedDefinitions,
    legacy_parameters: Mapping[str, ParamSetting],
) -> SealedConfiguration:
    """Build configuration from the fully initialized active legacy catalog."""
    definition_ids = {definition.param_id for definition in definitions}
    legacy_ids = set(legacy_parameters)
    if definition_ids != legacy_ids:
        raise ConfigurationMismatchError(
            sorted(definition_ids - legacy_ids),
            sorted(legacy_ids - definition_ids),
        )

    configs: list[ParamConfig] = []
    for definition in definitions:
        parameter = legacy_parameters[definition.param_id]
        _validate_initial_state(
            definition.param_id,
            parameter.value,
            parameter.min,
            parameter.max,
            value_field="effective_value",
            error_type=InvalidConfigurationError,
        )
        configs.append(
            ParamConfig(
                param_id=definition.param_id,
                effective_value=parameter.value,
                lower_bound=parameter.min,
                upper_bound=parameter.max,
            )
        )

    configs_tuple = tuple(configs)
    index = {config.param_id: position for position, config in enumerate(configs_tuple)}
    return SealedConfiguration(
        _configs=configs_tuple,
        _index=index,
        definitions_identity=definitions.identity,
    )
