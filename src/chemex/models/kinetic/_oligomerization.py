from __future__ import annotations

import math
import sys

from scipy.optimize import brentq

_MASS_BALANCE_TOLERANCE = 64.0 * sys.float_info.epsilon


def solve_oligomerization_fractions(
    terms: tuple[tuple[int, float], ...],
) -> tuple[float, tuple[float, ...]]:
    """Solve a positive normalized oligomerization mass balance."""
    if not terms:
        msg = "Oligomerization concentration solver requires at least one term"
        raise RuntimeError(msg)

    weighted_terms: list[tuple[int, float, float]] = []
    for stoichiometry, coefficient in terms:
        weighted_coefficient = stoichiometry * coefficient
        if (
            stoichiometry < 2
            or not math.isfinite(coefficient)
            or coefficient <= 0.0
            or not math.isfinite(weighted_coefficient)
        ):
            msg = "Oligomerization concentration solver received an invalid term"
            raise RuntimeError(msg)
        weighted_terms.append((stoichiometry, coefficient, weighted_coefficient))

    def mass_balance(monomer_fraction: float) -> float:
        contributions = [
            weighted_coefficient * monomer_fraction**stoichiometry
            for stoichiometry, _, weighted_coefficient in weighted_terms
        ]
        return math.fsum((monomer_fraction, *contributions)) - 1.0

    monomer_fraction, result = brentq(
        mass_balance,
        0.0,
        1.0,
        xtol=sys.float_info.min,
        rtol=4.0 * sys.float_info.epsilon,
        maxiter=256,
        full_output=True,
        disp=False,
    )
    if not result.converged:
        msg = "Oligomerization concentration solver did not converge"
        raise RuntimeError(msg)
    if (
        not math.isfinite(monomer_fraction)
        or monomer_fraction < 0.0
        or monomer_fraction > 1.0
    ):
        msg = (
            "Oligomerization concentration solver returned an invalid monomer fraction"
        )
        raise RuntimeError(msg)

    oligomer_fractions = tuple(
        coefficient * monomer_fraction**stoichiometry
        for stoichiometry, coefficient, _ in weighted_terms
    )
    if any(
        not math.isfinite(fraction) or fraction < 0.0 for fraction in oligomer_fractions
    ):
        msg = (
            "Oligomerization concentration solver returned an invalid oligomer fraction"
        )
        raise RuntimeError(msg)

    normalized_mass = math.fsum(
        (
            monomer_fraction,
            *(
                stoichiometry * fraction
                for (stoichiometry, _, _), fraction in zip(
                    weighted_terms,
                    oligomer_fractions,
                    strict=True,
                )
            ),
        ),
    )
    if abs(normalized_mass - 1.0) > _MASS_BALANCE_TOLERANCE:
        msg = "Oligomerization concentration solver violated mass conservation"
        raise RuntimeError(msg)

    return monomer_fraction, oligomer_fractions
