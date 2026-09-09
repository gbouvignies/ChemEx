"""Validated finite-pool equilibrium authority for one-ligand binding models."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

MIN_POSITIVE_FLOAT = math.nextafter(0.0, 1.0)

_LOG_TWO = math.log(2.0)
_LOG_FOUR = math.log(4.0)
_LOG_MIN_POSITIVE_FLOAT = math.log(MIN_POSITIVE_FLOAT)
_LOG_MIN_NORMAL_FLOAT = math.log(sys.float_info.min)
_LOG_MAX_FLOAT = math.log(sys.float_info.max)
_BALANCE_TOLERANCE = 64.0 * sys.float_info.epsilon
_LOG_EQUILIBRIUM_TOLERANCE = 4096.0 * sys.float_info.epsilon


@dataclass(frozen=True, slots=True)
class BindingEquilibrium:
    """Physical species composition for one finite ligand pool."""

    protein_free: float
    ligand_free: float
    bound_total: float
    free_proteins: tuple[float, ...]
    free_ligands: tuple[float, ...]
    complexes: tuple[float, ...]
    populations: tuple[float, ...]
    log_populations: tuple[float, ...]
    free_protein_log_fractions: tuple[float, ...]
    free_ligand_log_fractions: tuple[float, ...]
    complex_log_fractions: tuple[float, ...]
    edge_log_ratios: tuple[float, ...]
    log_apparent_kd: float


def validate_binding_totals(p_total: float, l_total: float) -> None:
    """Validate the finite-pool concentration domain used by binding models."""
    if not math.isfinite(p_total) or p_total <= 0.0:
        msg = "Binding P_total must be finite and strictly positive"
        raise ValueError(msg)
    if not math.isfinite(l_total) or l_total < 0.0:
        msg = "Binding L_total must be finite and non-negative"
        raise ValueError(msg)


def detailed_balance_rate(
    reverse_rate: float,
    log_source_population: float,
    log_destination_population: float,
) -> float:
    """Construct a tagged forward rate from the canonical equilibrium ratio."""
    if not math.isfinite(reverse_rate) or reverse_rate < 0.0:
        msg = "Binding reverse rate must be finite and non-negative"
        raise ValueError(msg)
    if reverse_rate == 0.0 or log_destination_population == -math.inf:
        return 0.0
    if not math.isfinite(log_destination_population):
        msg = "Binding destination population log is invalid"
        raise ValueError(msg)
    if log_source_population == -math.inf:
        msg = "Tagged binding forward rate exceeds binary64 representability"
        raise ValueError(msg)
    if not math.isfinite(log_source_population):
        msg = "Binding source population log is invalid"
        raise ValueError(msg)
    log_rate = (
        math.log(reverse_rate) + log_destination_population - log_source_population
    )
    if log_rate > _LOG_MAX_FLOAT:
        msg = "Tagged binding forward rate exceeds binary64 representability"
        raise ValueError(msg)
    if log_rate < _LOG_MIN_POSITIVE_FLOAT:
        msg = "Positive tagged binding forward rate is below binary64 representability"
        raise ValueError(msg)
    rate = math.exp(log_rate)
    if not math.isfinite(rate) or rate <= 0.0:
        msg = "Tagged binding forward rate is outside binary64 representability"
        raise ValueError(msg)
    return rate


def split_exchange_rate(
    total_rate: float,
    log_source_weight: float,
    log_destination_weight: float,
) -> tuple[float, float]:
    """Split a total exchange scale by an equilibrium ratio, including zero."""
    if not math.isfinite(total_rate) or total_rate < 0.0:
        msg = "Binding exchange-rate scale must be finite and non-negative"
        raise ValueError(msg)
    if total_rate == 0.0:
        return 0.0, 0.0
    weights = (log_source_weight, log_destination_weight)
    if any(math.isnan(weight) or weight == math.inf for weight in weights) or all(
        weight == -math.inf for weight in weights
    ):
        msg = "Binding exchange-rate weights are invalid"
        raise ValueError(msg)
    log_source_fraction, log_destination_fraction = _normalized_log_weights(weights)
    log_total_rate = math.log(total_rate)
    forward = (
        0.0
        if log_destination_fraction == -math.inf
        else _checked_positive_log_value(
            log_total_rate + log_destination_fraction,
            description="Binding forward exchange rate",
        )
    )
    reverse = (
        0.0
        if log_source_fraction == -math.inf
        else _checked_positive_log_value(
            log_total_rate + log_source_fraction,
            description="Binding reverse exchange rate",
        )
    )
    rates = [forward, reverse]
    correction = total_rate - math.fsum(rates)
    largest = max(range(len(rates)), key=rates.__getitem__)
    rates[largest] += correction
    if rates[largest] <= 0.0 or not math.isfinite(rates[largest]):
        msg = "Binding exchange rates could not preserve their total scale"
        raise RuntimeError(msg)
    return rates[0], rates[1]


def log_binding_weight(kd: float) -> float:
    """Return log(1 / KD) without constructing a possibly overflowing inverse."""
    if not math.isfinite(kd) or kd <= 0.0:
        msg = "Binding KD must be finite and strictly positive"
        raise ValueError(msg)
    return -math.log(kd)


def log_equilibrium_ratio(ratio: float) -> float:
    """Return a log equilibrium weight, preserving an exact-zero state removal."""
    if not math.isfinite(ratio) or ratio < 0.0:
        msg = "Binding equilibrium ratio must be finite and non-negative"
        raise ValueError(msg)
    return -math.inf if ratio == 0.0 else math.log(ratio)


def report_only_positive_log_value(log_value: float) -> float:
    """Return a positive report value or a non-finite typed-value sentinel."""
    if math.isnan(log_value) or log_value < _LOG_MIN_POSITIVE_FLOAT:
        return math.nan
    if log_value > _LOG_MAX_FLOAT:
        return math.inf
    return math.exp(log_value)


def _logsumexp(values: tuple[float, ...]) -> float:
    maximum = max(values)
    if maximum == -math.inf:
        return -math.inf
    return maximum + math.log(math.fsum(math.exp(value - maximum) for value in values))


def _validate_log_weights(
    weights: tuple[float, ...],
    *,
    description: str,
) -> None:
    if not weights or any(
        math.isnan(weight) or weight == math.inf for weight in weights
    ):
        msg = f"Binding {description} weights must be finite or exact-zero"
        raise ValueError(msg)
    if all(weight == -math.inf for weight in weights):
        msg = f"Binding {description} weights must retain at least one species"
        raise ValueError(msg)


def _normalized_weights(log_weights: tuple[float, ...]) -> tuple[float, ...]:
    maximum = max(log_weights)
    raw_weights = [
        0.0 if log_weight == -math.inf else math.exp(log_weight - maximum)
        for log_weight in log_weights
    ]
    raw_total = math.fsum(raw_weights)
    weights = [weight / raw_total for weight in raw_weights]
    correction = 1.0 - math.fsum(weights)
    largest = max(range(len(weights)), key=weights.__getitem__)
    weights[largest] += correction
    if (
        abs(correction) > _BALANCE_TOLERANCE
        or weights[largest] < 0.0
        or not math.isfinite(weights[largest])
    ):
        msg = "Binding equilibrium weights could not be normalized"
        raise RuntimeError(msg)
    return tuple(weights)


def _normalized_log_weights(log_weights: tuple[float, ...]) -> tuple[float, ...]:
    log_total = _logsumexp(log_weights)
    return tuple(
        -math.inf if log_weight == -math.inf else log_weight - log_total
        for log_weight in log_weights
    )


def _checked_positive_log_value(log_value: float, *, description: str) -> float:
    if log_value > _LOG_MAX_FLOAT:
        msg = f"{description} exceeds binary64 representability"
        raise ValueError(msg)
    if log_value < _LOG_MIN_POSITIVE_FLOAT:
        msg = f"Positive {description.lower()} is below binary64 representability"
        raise ValueError(msg)
    value = math.exp(log_value)
    if not math.isfinite(value) or value <= 0.0:
        msg = f"{description} is outside binary64 representability"
        raise ValueError(msg)
    return value


def _distribute(total: float, fractions: tuple[float, ...]) -> tuple[float, ...]:
    if total == 0.0:
        return tuple(0.0 for _ in fractions)
    concentrations = [total * fraction for fraction in fractions]
    correction = total - math.fsum(concentrations)
    largest = max(range(len(concentrations)), key=concentrations.__getitem__)
    concentrations[largest] += correction
    # A subnormal macro pool may be representable even when none of its exact
    # species shares are.  Preserve the pool mass in one concrete value while
    # the accompanying log fractions retain the equilibrium composition.
    correction_tolerance = max(
        _BALANCE_TOLERANCE * total,
        len(fractions) * MIN_POSITIVE_FLOAT,
    )
    if (
        abs(correction) > correction_tolerance
        or concentrations[largest] < 0.0
        or not math.isfinite(concentrations[largest])
    ):
        msg = "Binding species concentrations could not be normalized"
        raise RuntimeError(msg)
    return tuple(concentrations)


def _finite_pool_totals(
    p_total: float,
    l_total: float,
    log_apparent_kd: float,
) -> tuple[float, float, float, float, float, float]:
    """Solve P + L <-> PL without cancellation in either occupancy limit."""
    log_protein = math.log(p_total)
    log_ligand = math.log(l_total)
    log_sum = _logsumexp((log_protein, log_ligand, log_apparent_kd))
    difference = abs(p_total - l_total)
    log_difference = -math.inf if difference == 0.0 else math.log(difference)
    log_discriminant = _logsumexp(
        (
            2.0 * log_difference,
            _LOG_TWO + log_apparent_kd + _logsumexp((log_protein, log_ligand)),
            2.0 * log_apparent_kd,
        ),
    )
    log_sqrt_discriminant = 0.5 * log_discriminant
    log_bound = (
        _LOG_TWO
        + log_protein
        + log_ligand
        - _logsumexp((log_sum, log_sqrt_discriminant))
    )

    limiting = min(p_total, l_total)
    log_limiting = math.log(limiting)
    if log_bound - log_limiting <= -_LOG_TWO:
        bound = math.exp(log_bound)
        protein_free = p_total - bound
        ligand_free = l_total - bound
        log_protein_free = math.log(protein_free)
        log_ligand_free = math.log(ligand_free)
    else:
        log_linear = _logsumexp((log_difference, log_apparent_kd))
        log_sqrt = 0.5 * _logsumexp(
            (
                2.0 * log_linear,
                _LOG_FOUR + log_apparent_kd + log_limiting,
            ),
        )
        log_free_limiting = (
            _LOG_TWO
            + log_apparent_kd
            + log_limiting
            - _logsumexp((log_linear, log_sqrt))
        )
        free_limiting = (
            _equal_pool_free(limiting, log_apparent_kd)
            if difference == 0.0
            else math.exp(log_free_limiting)
        )
        if difference == 0.0 and free_limiting > 0.0:
            log_free_limiting = math.log(free_limiting)
        bound = limiting - free_limiting
        if p_total <= l_total:
            protein_free = free_limiting
            ligand_free = difference + free_limiting
            log_protein_free = log_free_limiting
            log_ligand_free = _logsumexp((log_difference, log_free_limiting))
        else:
            ligand_free = free_limiting
            protein_free = difference + free_limiting
            log_ligand_free = log_free_limiting
            log_protein_free = _logsumexp((log_difference, log_free_limiting))

    return (
        protein_free,
        ligand_free,
        bound,
        log_protein_free,
        log_ligand_free,
        log_bound,
    )


def _equal_pool_free(total: float, log_apparent_kd: float) -> float:
    """Evaluate the equal-pool free root without log subtraction error."""
    if not _LOG_MIN_POSITIVE_FLOAT <= log_apparent_kd <= _LOG_MAX_FLOAT:
        log_total = math.log(total)
        log_free = 0.5 * (log_apparent_kd + log_total)
        return math.exp(log_free)
    kd = math.exp(log_apparent_kd)
    scale = max(kd, total)
    sqrt_kd = math.sqrt(kd / scale)
    sqrt_sum = math.sqrt(kd / scale + 4.0 * (total / scale))
    fraction = 2.0 * sqrt_kd / (sqrt_kd + sqrt_sum)
    return total * fraction


def _validate_concentration_domain(equilibrium: BindingEquilibrium) -> None:
    concentrations = (
        equilibrium.protein_free,
        equilibrium.ligand_free,
        equilibrium.bound_total,
        *equilibrium.free_proteins,
        *equilibrium.free_ligands,
        *equilibrium.complexes,
        *equilibrium.populations,
    )
    if any(not math.isfinite(value) or value < 0.0 for value in concentrations):
        msg = "Binding equilibrium produced an invalid concentration"
        raise RuntimeError(msg)


def _validate_mass_balance(
    equilibrium: BindingEquilibrium,
    p_total: float,
    l_total: float,
) -> None:
    if not math.isclose(
        math.fsum((equilibrium.protein_free, equilibrium.bound_total)),
        p_total,
        rel_tol=_BALANCE_TOLERANCE,
        abs_tol=MIN_POSITIVE_FLOAT,
    ):
        msg = "Binding equilibrium violated protein mass conservation"
        raise RuntimeError(msg)
    if not math.isclose(
        math.fsum((equilibrium.ligand_free, equilibrium.bound_total)),
        l_total,
        rel_tol=_BALANCE_TOLERANCE,
        abs_tol=MIN_POSITIVE_FLOAT,
    ):
        msg = "Binding equilibrium violated ligand mass conservation"
        raise RuntimeError(msg)
    if (
        not math.isclose(
            math.fsum(equilibrium.free_proteins),
            equilibrium.protein_free,
            rel_tol=_BALANCE_TOLERANCE,
            abs_tol=MIN_POSITIVE_FLOAT,
        )
        or not math.isclose(
            math.fsum(equilibrium.free_ligands),
            equilibrium.ligand_free,
            rel_tol=_BALANCE_TOLERANCE,
            abs_tol=MIN_POSITIVE_FLOAT,
        )
        or not math.isclose(
            math.fsum(equilibrium.complexes),
            equilibrium.bound_total,
            rel_tol=_BALANCE_TOLERANCE,
            abs_tol=MIN_POSITIVE_FLOAT,
        )
    ):
        msg = "Binding equilibrium species do not conserve their pools"
        raise RuntimeError(msg)


def _validate_populations(
    equilibrium: BindingEquilibrium,
    p_total: float,
) -> None:
    if not math.isclose(
        math.fsum(equilibrium.populations),
        1.0,
        rel_tol=_BALANCE_TOLERANCE,
        abs_tol=_BALANCE_TOLERANCE,
    ):
        msg = "Binding equilibrium populations are not normalized"
        raise RuntimeError(msg)
    expected_populations = tuple(
        0.0 if log_population == -math.inf else math.exp(log_population)
        for log_population in equilibrium.log_populations
    )
    if any(
        not math.isclose(
            actual,
            expected,
            rel_tol=_BALANCE_TOLERANCE,
            abs_tol=_BALANCE_TOLERANCE,
        )
        for actual, expected in zip(
            equilibrium.populations,
            expected_populations,
            strict=True,
        )
    ):
        msg = "Binding equilibrium populations do not match their log authority"
        raise RuntimeError(msg)
    log_p_total = math.log(p_total)
    species = (*equilibrium.free_proteins, *equilibrium.complexes)
    for concentration, log_population in zip(
        species,
        equilibrium.log_populations,
        strict=True,
    ):
        log_concentration = log_population + log_p_total
        if concentration == 0.0:
            consistent = (
                log_population == -math.inf
                or log_concentration <= _LOG_MIN_POSITIVE_FLOAT
            )
        elif (
            concentration < sys.float_info.min
            and log_concentration < _LOG_MIN_NORMAL_FLOAT
        ):
            consistent = True
        else:
            consistent = abs(math.log(concentration) - log_concentration) <= (
                _LOG_EQUILIBRIUM_TOLERANCE * max(1.0, abs(log_concentration))
            )
        if not consistent:
            msg = "Binding equilibrium populations disagree with species concentrations"
            raise RuntimeError(msg)


def _validate_species_distribution(
    concentrations: tuple[float, ...],
    total: float,
    log_weights: tuple[float, ...],
    fractions: tuple[float, ...],
    *,
    description: str,
) -> None:
    for concentration, log_weight, fraction in zip(
        concentrations,
        log_weights,
        fractions,
        strict=True,
    ):
        if log_weight == -math.inf and concentration != 0.0:
            msg = f"Binding equilibrium retained a zero-weight {description}"
            raise RuntimeError(msg)
        expected = total * fraction
        if not math.isclose(
            concentration,
            expected,
            rel_tol=_BALANCE_TOLERANCE,
            abs_tol=4.0 * MIN_POSITIVE_FLOAT,
        ):
            msg = f"Binding equilibrium violated a {description} relation"
            raise RuntimeError(msg)


def _validate_equilibrium(
    equilibrium: BindingEquilibrium,
    *,
    protein_total: float,
    ligand_total: float,
    log_protein_free: float,
    log_ligand_free: float,
    log_bound: float,
    free_protein_log_weights: tuple[float, ...],
    free_ligand_log_weights: tuple[float, ...],
    complex_log_weights: tuple[float, ...],
    free_protein_fractions: tuple[float, ...],
    free_ligand_fractions: tuple[float, ...],
    complex_fractions: tuple[float, ...],
) -> None:
    _validate_concentration_domain(equilibrium)
    _validate_mass_balance(equilibrium, protein_total, ligand_total)
    _validate_populations(equilibrium, protein_total)
    _validate_species_distribution(
        equilibrium.free_proteins,
        equilibrium.protein_free,
        free_protein_log_weights,
        free_protein_fractions,
        description="free-protein species",
    )
    if ligand_total > 0.0:
        log_residual = (
            log_protein_free + log_ligand_free - log_bound - equilibrium.log_apparent_kd
        )
        if abs(log_residual) > _LOG_EQUILIBRIUM_TOLERANCE * max(
            1.0,
            abs(equilibrium.log_apparent_kd),
        ):
            msg = "Binding equilibrium violated the association law"
            raise RuntimeError(msg)
        # The aggregate association law above and these normalized scientific
        # weights together imply every individual species relation.  Compare
        # against the fractions used to construct the result so the validator
        # does not reject representable subnormal values merely because a
        # second log/exp reconstruction rounds differently.
        _validate_species_distribution(
            equilibrium.free_ligands,
            equilibrium.ligand_free,
            free_ligand_log_weights,
            free_ligand_fractions,
            description="free-ligand species",
        )
        _validate_species_distribution(
            equilibrium.complexes,
            equilibrium.bound_total,
            complex_log_weights,
            complex_fractions,
            description="complex",
        )


def solve_binding_equilibrium(
    p_total: float,
    l_total: float,
    *,
    free_protein_log_weights: tuple[float, ...] = (0.0,),
    free_ligand_log_weights: tuple[float, ...],
    complex_log_weights: tuple[float, ...],
    edge_log_ratios: tuple[float, ...] = (),
) -> BindingEquilibrium:
    """Solve one finite ligand pool and distribute its free and bound species."""
    validate_binding_totals(p_total, l_total)
    _validate_log_weights(free_protein_log_weights, description="free-protein")
    _validate_log_weights(free_ligand_log_weights, description="free-ligand")
    _validate_log_weights(complex_log_weights, description="complex")
    if any(math.isnan(ratio) or ratio == math.inf for ratio in edge_log_ratios):
        msg = "Binding edge equilibrium ratios must be finite or exact-zero"
        raise ValueError(msg)

    protein_fractions = _normalized_weights(free_protein_log_weights)
    free_fractions = _normalized_weights(free_ligand_log_weights)
    complex_fractions = _normalized_weights(complex_log_weights)
    protein_log_fractions = _normalized_log_weights(free_protein_log_weights)
    free_log_fractions = _normalized_log_weights(free_ligand_log_weights)
    complex_log_fractions = _normalized_log_weights(complex_log_weights)
    log_protein_weight = _logsumexp(free_protein_log_weights)
    log_free_weight = _logsumexp(free_ligand_log_weights)
    log_complex_weight = _logsumexp(complex_log_weights)
    log_apparent_kd = log_protein_weight + log_free_weight - log_complex_weight

    if l_total == 0.0:
        protein_free = p_total
        ligand_free = bound_total = 0.0
        log_protein_free = math.log(p_total)
        log_ligand_free = log_bound = -math.inf
    else:
        (
            protein_free,
            ligand_free,
            bound_total,
            log_protein_free,
            log_ligand_free,
            log_bound,
        ) = _finite_pool_totals(p_total, l_total, log_apparent_kd)
    free_proteins = _distribute(protein_free, protein_fractions)
    free_ligands = _distribute(ligand_free, free_fractions)
    complexes = _distribute(bound_total, complex_fractions)
    log_p_total = math.log(p_total)
    raw_log_populations = (
        *(
            -math.inf
            if log_weight == -math.inf
            else log_protein_free + log_weight - log_protein_weight - log_p_total
            for log_weight in free_protein_log_weights
        ),
        *(
            -math.inf
            if log_weight == -math.inf
            else log_bound + log_weight - log_complex_weight - log_p_total
            for log_weight in complex_log_weights
        ),
    )
    populations = _normalized_weights(raw_log_populations)
    log_population_total = _logsumexp(raw_log_populations)
    log_populations = tuple(
        -math.inf
        if log_population == -math.inf
        else log_population - log_population_total
        for log_population in raw_log_populations
    )
    equilibrium = BindingEquilibrium(
        protein_free=protein_free,
        ligand_free=ligand_free,
        bound_total=bound_total,
        free_proteins=free_proteins,
        free_ligands=free_ligands,
        complexes=complexes,
        populations=populations,
        log_populations=log_populations,
        free_protein_log_fractions=protein_log_fractions,
        free_ligand_log_fractions=free_log_fractions,
        complex_log_fractions=complex_log_fractions,
        edge_log_ratios=edge_log_ratios,
        log_apparent_kd=log_apparent_kd,
    )
    _validate_equilibrium(
        equilibrium,
        protein_total=p_total,
        ligand_total=l_total,
        log_protein_free=log_protein_free,
        log_ligand_free=log_ligand_free,
        log_bound=log_bound,
        free_protein_log_weights=free_protein_log_weights,
        free_ligand_log_weights=free_ligand_log_weights,
        complex_log_weights=complex_log_weights,
        free_protein_fractions=protein_fractions,
        free_ligand_fractions=free_fractions,
        complex_fractions=complex_fractions,
    )
    return equilibrium
