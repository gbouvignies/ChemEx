from __future__ import annotations

import math
import sys
from dataclasses import dataclass

from scipy.optimize import brentq

MIN_POSITIVE_FLOAT = math.nextafter(0.0, 1.0)
LOG_MIN_POSITIVE_FLOAT = math.log(MIN_POSITIVE_FLOAT)

_MASS_BALANCE_TOLERANCE = 64.0 * sys.float_info.epsilon
_LOG_BALANCE_TOLERANCE = 4096.0 * sys.float_info.epsilon
_LOG_MAX_FLOAT = math.log(sys.float_info.max)


@dataclass(frozen=True, slots=True)
class OligomerizationEquilibrium:
    """Normalized species fractions and their retained logarithms."""

    monomer_fraction: float
    oligomer_fractions: tuple[float, ...]
    log_monomer_fraction: float
    log_oligomer_fractions: tuple[float, ...]
    log_tagged_fractions: tuple[float, ...]


def validate_oligomerization_kd(kd: float) -> float:
    """Return a KD only when it has finite reversible-kinetic semantics."""
    if not math.isfinite(kd) or kd <= 0.0:
        msg = "Oligomerization KD must be finite and strictly positive"
        raise ValueError(msg)
    return kd


def log_equilibrium_coefficient(
    p_total: float,
    p_total_power: int,
    *kds: float,
) -> float:
    """Build an equilibrium coefficient logarithm without eager products."""
    if not math.isfinite(p_total) or p_total <= 0.0:
        msg = "Oligomerization P_total must be finite and strictly positive"
        raise ValueError(msg)
    if p_total_power < 0:
        msg = "Oligomerization concentration solver received an invalid power"
        raise ValueError(msg)
    validated_kds = tuple(validate_oligomerization_kd(kd) for kd in kds)
    return p_total_power * math.log(p_total) - math.fsum(
        math.log(kd) for kd in validated_kds
    )


def _logsumexp(values: tuple[float, ...]) -> float:
    maximum = max(values)
    return maximum + math.log(math.fsum(math.exp(value - maximum) for value in values))


def solve_oligomerization_fractions(
    terms: tuple[tuple[int, float], ...],
) -> OligomerizationEquilibrium:
    """Solve normalized oligomerization mass balance from log coefficients."""
    if not terms:
        msg = "Oligomerization concentration solver requires at least one term"
        raise ValueError(msg)

    weighted_terms: list[tuple[int, float, float]] = []
    for stoichiometry, log_coefficient in terms:
        if (
            stoichiometry < 2
            or not math.isfinite(log_coefficient)
            or isinstance(stoichiometry, bool)
        ):
            msg = "Oligomerization concentration solver received an invalid term"
            raise ValueError(msg)
        weighted_terms.append(
            (
                stoichiometry,
                log_coefficient,
                math.log(stoichiometry) + log_coefficient,
            )
        )

    def log_mass_balance(log_monomer_fraction: float) -> float:
        return _logsumexp(
            (
                log_monomer_fraction,
                *(
                    log_weighted_coefficient + stoichiometry * log_monomer_fraction
                    for stoichiometry, _, log_weighted_coefficient in weighted_terms
                ),
            ),
        )

    # At y=0 the monomer term alone is one, so the balance is non-negative.
    # This finite lower endpoint makes every one of N positive mass terms at
    # most exp(-1) / N, giving a deterministic strict negative bracket.
    intercept_maximum = max(
        0.0,
        *(log_weighted for _, _, log_weighted in weighted_terms),
    )
    term_count = len(weighted_terms) + 1
    lower = -(intercept_maximum + math.log(term_count) + 1.0)
    log_monomer_fraction, result = brentq(
        log_mass_balance,
        lower,
        0.0,
        xtol=MIN_POSITIVE_FLOAT,
        rtol=4.0 * sys.float_info.epsilon,
        maxiter=256,
        full_output=True,
        disp=False,
    )
    if not result.converged:
        msg = "Oligomerization concentration solver did not converge"
        raise RuntimeError(msg)
    if not math.isfinite(log_monomer_fraction) or log_monomer_fraction > 0.0:
        msg = (
            "Oligomerization concentration solver returned an invalid monomer fraction"
        )
        raise RuntimeError(msg)

    raw_log_oligomer_fractions = tuple(
        log_coefficient + stoichiometry * log_monomer_fraction
        for stoichiometry, log_coefficient, _ in weighted_terms
    )
    raw_log_tagged_fractions = (
        log_monomer_fraction,
        *(
            math.log(stoichiometry) + log_fraction
            for (stoichiometry, _, _), log_fraction in zip(
                weighted_terms,
                raw_log_oligomer_fractions,
                strict=True,
            )
        ),
    )
    log_normalizer = _logsumexp(raw_log_tagged_fractions)
    if abs(log_normalizer) > _LOG_BALANCE_TOLERANCE:
        msg = "Oligomerization concentration solver violated mass conservation"
        raise RuntimeError(msg)
    log_tagged_fractions = tuple(
        log_fraction - log_normalizer for log_fraction in raw_log_tagged_fractions
    )
    tagged_fractions = [math.exp(log_fraction) for log_fraction in log_tagged_fractions]

    # Exponentiation can miss one last bit of unit mass. Correct only the
    # largest tagged component, deterministically, at floating-point scale.
    normalized_mass = math.fsum(tagged_fractions)
    correction = 1.0 - normalized_mass
    largest_index = max(range(len(tagged_fractions)), key=tagged_fractions.__getitem__)
    corrected_largest = tagged_fractions[largest_index] + correction
    if (
        abs(correction) > _MASS_BALANCE_TOLERANCE
        or corrected_largest < 0.0
        or not math.isfinite(corrected_largest)
    ):
        msg = "Oligomerization concentration solver violated mass conservation"
        raise RuntimeError(msg)
    tagged_fractions[largest_index] = corrected_largest

    monomer_fraction = tagged_fractions[0]
    oligomer_fractions = tuple(
        tagged_fraction / stoichiometry
        for tagged_fraction, (stoichiometry, _, _) in zip(
            tagged_fractions[1:],
            weighted_terms,
            strict=True,
        )
    )
    log_oligomer_fractions = tuple(
        log_tagged_fraction - math.log(stoichiometry)
        for log_tagged_fraction, (stoichiometry, _, _) in zip(
            log_tagged_fractions[1:],
            weighted_terms,
            strict=True,
        )
    )

    return OligomerizationEquilibrium(
        monomer_fraction=monomer_fraction,
        oligomer_fractions=oligomer_fractions,
        log_monomer_fraction=log_tagged_fractions[0],
        log_oligomer_fractions=log_oligomer_fractions,
        log_tagged_fractions=log_tagged_fractions,
    )


def _validate_positive_rate_log(log_rate: float, *, description: str) -> None:
    if not math.isfinite(log_rate) or log_rate > _LOG_MAX_FLOAT:
        msg = f"{description} exceeds the maximum finite binary64 value"
        raise ValueError(msg)
    if log_rate < LOG_MIN_POSITIVE_FLOAT:
        msg = f"Positive {description} is below binary64 representability"
        raise ValueError(msg)


def _positive_rate_from_log(log_rate: float, *, description: str) -> float:
    _validate_positive_rate_log(log_rate, description=description)
    rate = math.exp(log_rate)
    if rate == 0.0:
        msg = f"Positive {description} is below binary64 representability"
        raise ValueError(msg)
    if not math.isfinite(rate):
        msg = f"{description} exceeds the maximum finite binary64 value"
        raise ValueError(msg)
    return rate


def concentration_from_log_fraction(p_total: float, log_fraction: float) -> float:
    """Reconstruct a concentration without first rounding its fraction."""
    if not math.isfinite(p_total) or p_total <= 0.0:
        msg = "Oligomerization P_total must be finite and strictly positive"
        raise ValueError(msg)
    if not math.isfinite(log_fraction):
        msg = "Oligomerization species fraction log must be finite"
        raise ValueError(msg)

    log_concentration = math.log(p_total) + log_fraction
    if log_concentration < LOG_MIN_POSITIVE_FLOAT:
        return 0.0
    if log_concentration > _LOG_MAX_FLOAT:
        msg = "Oligomerization concentration exceeds the maximum finite binary64 value"
        raise ValueError(msg)
    return math.exp(log_concentration)


def concentrations_from_log_fractions(
    p_total: float,
    species: tuple[tuple[int, float], ...],
) -> tuple[float, ...]:
    """Reconstruct species concentrations with a last-bit mass correction."""
    concentrations = [
        concentration_from_log_fraction(p_total, log_fraction)
        for _, log_fraction in species
    ]
    tagged_mass_fractions = [
        stoichiometry * (concentration / p_total)
        for concentration, (stoichiometry, _) in zip(
            concentrations,
            species,
            strict=True,
        )
    ]
    correction = 1.0 - math.fsum(tagged_mass_fractions)
    largest_index = max(
        range(len(concentrations)),
        key=tagged_mass_fractions.__getitem__,
    )
    largest_stoichiometry = species[largest_index][0]
    corrected_largest = (
        concentrations[largest_index] + correction * p_total / largest_stoichiometry
    )
    if (
        abs(correction) > _LOG_BALANCE_TOLERANCE
        or corrected_largest < 0.0
        or not math.isfinite(corrected_largest)
    ):
        msg = "Oligomerization concentration reconstruction violated mass conservation"
        raise RuntimeError(msg)
    concentrations[largest_index] = corrected_largest
    return tuple(concentrations)


def scale_reversible_rate(rate: float, factor: float) -> float:
    """Scale a non-negative rate while rejecting structural under/overflow."""
    if not math.isfinite(rate) or rate < 0.0:
        msg = "Oligomerization reverse rate must be finite and non-negative"
        raise ValueError(msg)
    if not math.isfinite(factor) or factor <= 0.0:
        msg = "Oligomerization rate scale must be finite and strictly positive"
        raise ValueError(msg)
    if rate == 0.0:
        return 0.0
    log_scaled_rate = math.log(rate) + math.log(factor)
    _validate_positive_rate_log(
        log_scaled_rate,
        description="tagged reverse rate",
    )
    scaled_rate = rate * factor
    if scaled_rate == 0.0:
        msg = "Positive tagged reverse rate is below binary64 representability"
        raise ValueError(msg)
    if not math.isfinite(scaled_rate):
        msg = "Tagged reverse rate exceeds the maximum finite binary64 value"
        raise ValueError(msg)
    return scaled_rate


def detailed_balance_forward_rate(
    reverse_rate: float,
    log_source_fraction: float,
    log_destination_fraction: float,
) -> float:
    """Construct one tagged forward rate without an eager population ratio."""
    if not math.isfinite(reverse_rate) or reverse_rate < 0.0:
        msg = "Oligomerization reverse rate must be finite and non-negative"
        raise ValueError(msg)
    if not math.isfinite(log_source_fraction) or not math.isfinite(
        log_destination_fraction
    ):
        msg = "Oligomerization tagged population logs must be finite"
        raise ValueError(msg)
    if reverse_rate == 0.0:
        return 0.0
    return _positive_rate_from_log(
        math.log(reverse_rate) + log_destination_fraction - log_source_fraction,
        description="tagged forward rate",
    )
