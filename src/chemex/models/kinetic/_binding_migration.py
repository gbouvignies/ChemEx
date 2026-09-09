"""Migration diagnostics for the two directional-rate Category-C models."""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Literal

from chemex.exceptions import ChemExError
from chemex.parameters.name import ParamName

if TYPE_CHECKING:
    from chemex.configuration.parameters import DefaultType


@dataclass(frozen=True, slots=True)
class BindingRateMigration:
    """One model's complete breaking directional-rate migration policy."""

    model_name: str
    forward_name: str
    reverse_name: str
    equilibrium_name: str
    exchange_name: str
    equilibrium_lower_bound: float
    equilibrium_upper_bound: float = 1.0e6
    exchange_lower_bound: float = 0.0
    exchange_upper_bound: float = 1.0e6

    @property
    def equilibrium_allows_zero(self) -> bool:
        return self.equilibrium_lower_bound == 0.0

    @property
    def legacy_names(self) -> frozenset[str]:
        return frozenset((self.forward_name, self.reverse_name))

    @property
    def replacement_names(self) -> frozenset[str]:
        return frozenset((self.equilibrium_name, self.exchange_name))


_MIN_POSITIVE_FLOAT = math.ulp(0.0)
_LOG_MIN_POSITIVE_FLOAT = math.log(_MIN_POSITIVE_FLOAT)
_LOG_MAX_FLOAT = math.log(sys.float_info.max)

# The declaration is deliberately close to the top: a maintainer can read each
# old pair, replacement pair, and zero policy without following parser code.
_MIGRATIONS = {
    "3st_binding_cs": BindingRateMigration(
        model_name="3st_binding_cs",
        forward_name="KAB",
        reverse_name="KBA",
        equilibrium_name="KEQ_AB",
        exchange_name="KEX_AB",
        equilibrium_lower_bound=_MIN_POSITIVE_FLOAT,
    ),
    "3st_binding_if": BindingRateMigration(
        model_name="3st_binding_if",
        forward_name="KBC",
        reverse_name="KCB",
        equilibrium_name="KEQ_BC",
        exchange_name="KEX_BC",
        equilibrium_lower_bound=0.0,
    ),
}


class LegacyBindingParameterError(ChemExError, ValueError):
    """A Category-C input uses directional rates that are now derived outputs."""

    def __init__(self, model_name: str, explanation: str) -> None:
        super().__init__(explanation)
        self.model_name = model_name
        self.explanation = explanation


type _UnrepresentableKind = Literal["below_minimum", "above_maximum"]


@dataclass(frozen=True, slots=True)
class _ConvertedLegacyValue:
    value: float | None
    unrepresentable: _UnrepresentableKind | None = None


@dataclass(frozen=True, slots=True)
class _LegacyRatePairConversion:
    equilibrium: _ConvertedLegacyValue
    exchange: _ConvertedLegacyValue


def binding_rate_migration(model_name: str) -> BindingRateMigration | None:
    return _MIGRATIONS.get(model_name)


def _convert_legacy_rate_pair(
    forward: float,
    reverse: float,
) -> _LegacyRatePairConversion:
    """Classify exact positive legacy relationships before formatting advice."""
    log_ratio = math.log(forward) - math.log(reverse)
    if log_ratio < _LOG_MIN_POSITIVE_FLOAT:
        equilibrium = _ConvertedLegacyValue(None, "below_minimum")
    elif log_ratio > _LOG_MAX_FLOAT:
        equilibrium = _ConvertedLegacyValue(None, "above_maximum")
    else:
        ratio = forward / reverse
        if ratio == 0.0:
            equilibrium = _ConvertedLegacyValue(None, "below_minimum")
        elif not math.isfinite(ratio):
            equilibrium = _ConvertedLegacyValue(None, "above_maximum")
        else:
            equilibrium = _ConvertedLegacyValue(ratio)

    exact_sum = Fraction.from_float(forward) + Fraction.from_float(reverse)
    exchange = (
        _ConvertedLegacyValue(None, "above_maximum")
        if exact_sum > Fraction.from_float(sys.float_info.max)
        else _ConvertedLegacyValue(float(exact_sum))
    )
    return _LegacyRatePairConversion(equilibrium=equilibrium, exchange=exchange)


def _scope_label(scope: ParamName | None) -> str:
    return "[GLOBAL]" if scope is None or not scope else str(scope)


def _bound_override_guidance(
    migration: BindingRateMigration,
    conversion: _LegacyRatePairConversion,
) -> str:
    guidance: list[str] = []
    replacements = (
        (
            migration.equilibrium_name,
            conversion.equilibrium.value,
            migration.equilibrium_lower_bound,
            migration.equilibrium_upper_bound,
        ),
        (
            migration.exchange_name,
            conversion.exchange.value,
            migration.exchange_lower_bound,
            migration.exchange_upper_bound,
        ),
    )
    for name, value, lower, upper in replacements:
        if value is not None and not lower <= value <= upper:
            guidance.append(
                f"{name} = {value!r} lies outside the new default bound "
                f"[{lower!r}, {upper!r}]. Override it explicitly in the parameter "
                f"file with {name} = [{value!r}, {lower!r}, {value!r}]."
            )
    return " ".join(guidance)


def _conversion_guidance(
    migration: BindingRateMigration,
    conversion: _LegacyRatePairConversion,
) -> str:
    equilibrium = conversion.equilibrium
    if equilibrium.value is not None:
        equilibrium_text = f"{migration.equilibrium_name} = {equilibrium.value!r}"
    elif equilibrium.unrepresentable == "below_minimum":
        equilibrium_text = (
            f"the equivalent {migration.equilibrium_name} is below minimum positive "
            "binary64 representability and requires manual model/parameter "
            "reconsideration"
        )
    else:
        equilibrium_text = (
            f"the equivalent {migration.equilibrium_name} exceeds maximum finite "
            "binary64 representability and requires manual model/parameter "
            "reconsideration"
        )

    exchange = conversion.exchange
    exchange_text = (
        f"{migration.exchange_name} = {exchange.value!r}"
        if exchange.value is not None
        else (
            f"the equivalent {migration.exchange_name} exceeds maximum finite "
            "binary64 representability and requires manual model/parameter "
            "reconsideration"
        )
    )
    bounds = _bound_override_guidance(migration, conversion)
    suffix = f" {bounds}" if bounds else ""
    return f"Use replacement values {equilibrium_text}; {exchange_text}.{suffix}"


def legacy_binding_migration_message(
    migration: BindingRateMigration,
    supplied_names: set[str],
    values: Mapping[str, float] | None = None,
    *,
    scope: ParamName | None = None,
) -> str:
    old = migration.legacy_names & supplied_names
    new = migration.replacement_names & supplied_names
    prefix = (
        f"Scope {_scope_label(scope)}: model '{migration.model_name}' no longer "
        f"accepts {migration.forward_name}/{migration.reverse_name} as independent "
        "inputs; the directional rates are derived outputs. "
    )
    formulas = (
        f"{migration.equilibrium_name} = {migration.forward_name} / "
        f"{migration.reverse_name}; {migration.exchange_name} = "
        f"{migration.forward_name} + {migration.reverse_name}."
    )
    if new:
        return (
            prefix + "The old and new parameterizations cannot be combined. " + formulas
        )
    if len(old) == 1:
        missing = next(iter(migration.legacy_names - old))
        return (
            prefix
            + f"The legacy pair is incomplete: {missing} is also needed to reconstruct "
            f"{migration.equilibrium_name}/{migration.exchange_name}. " + formulas
        )
    if values is None:
        return prefix + "Replace the legacy pair using: " + formulas

    forward = values[migration.forward_name]
    reverse = values[migration.reverse_name]
    if not math.isfinite(forward) or not math.isfinite(reverse):
        return prefix + "Legacy directional rates must be finite. " + formulas
    if forward < 0.0 or reverse < 0.0:
        return prefix + "Legacy directional rates must be non-negative. " + formulas
    if forward == 0.0 and reverse == 0.0:
        return (
            prefix
            + f"{migration.exchange_name} = 0; {migration.equilibrium_name} cannot be "
            "inferred from (0,0), so choose the equilibrium ratio explicitly."
        )
    if reverse == 0.0:
        return (
            prefix
            + f"A finite {migration.equilibrium_name} cannot be reconstructed when "
            f"{migration.forward_name} > 0 and {migration.reverse_name} = 0."
        )
    if forward == 0.0:
        if not migration.equilibrium_allows_zero:
            return (
                prefix + f"The legacy pair gives {migration.equilibrium_name} = 0, but "
                f"{migration.model_name} requires finite positive "
                f"{migration.equilibrium_name}; choose it explicitly."
            )
        conversion = _LegacyRatePairConversion(
            equilibrium=_ConvertedLegacyValue(0.0),
            exchange=_ConvertedLegacyValue(reverse),
        )
    else:
        conversion = _convert_legacy_rate_pair(forward, reverse)
    return prefix + _conversion_guidance(migration, conversion)


def validate_legacy_binding_defaults(
    model_name: str,
    defaults: Sequence[DefaultType],
) -> None:
    migration = binding_rate_migration(model_name)
    if migration is None:
        return

    grouped: dict[str, tuple[ParamName, dict[str, float]]] = {}
    for name, setting in defaults:
        if name.name not in migration.legacy_names | migration.replacement_names:
            continue
        scope = ParamName("", name.spin_system, name.conditions)
        _existing_scope, supplied = grouped.setdefault(scope.id_, (scope, {}))
        supplied[name.name] = setting.value

    messages: list[str] = []
    for scope, supplied in grouped.values():
        supplied_names = set(supplied)
        if not migration.legacy_names & supplied_names:
            continue
        values = supplied if migration.legacy_names <= supplied_names else None
        messages.append(
            legacy_binding_migration_message(
                migration,
                supplied_names,
                values,
                scope=scope,
            )
        )
    if messages:
        raise LegacyBindingParameterError(model_name, "\n".join(messages))
