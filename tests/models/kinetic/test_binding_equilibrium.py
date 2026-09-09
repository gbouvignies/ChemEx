"""Independent numerical qualifications for finite-pool binding equilibrium."""

from __future__ import annotations

import math
import sys
from decimal import Decimal, localcontext

import pytest

from chemex.models.kinetic._binding import (
    log_binding_weight,
    log_equilibrium_ratio,
    solve_binding_equilibrium,
    split_exchange_rate,
)


def _decimal_one_site(
    p_total: float,
    l_total: float,
    kd: float,
) -> tuple[float, float, float]:
    """Return an independent stable high-precision physical solution."""
    values = tuple(Decimal.from_float(value) for value in (p_total, l_total, kd))
    exponent_span = max(value.adjusted() for value in values) - min(
        value.adjusted() for value in values
    )
    with localcontext() as context:
        context.prec = exponent_span + 100
        protein, ligand, dissociation = values
        total = protein + ligand + dissociation
        discriminant = (total * total - Decimal(4) * protein * ligand).sqrt()
        bound = Decimal(2) * protein * ligand / (total + discriminant)
        return float(protein - bound), float(ligand - bound), float(bound)


@pytest.mark.parametrize(
    ("p_total", "l_total", "kd"),
    (
        (1.0e-3, 1.0e-3, 1.0e-300),
        (1.0e-100, 1.0e-100, 1.0e-300),
        (1.0e-9, 2.0e-3, 1.0e-50),
        (2.0e-3, 1.0e-9, 1.0e-50),
        (1.0e-12, 2.0e-12, 3.0e-12),
        (1.0e-3, 2.0e-3, 1.0e6),
    ),
)
def test_one_site_binding_matches_high_precision_physical_root(
    p_total: float,
    l_total: float,
    kd: float,
) -> None:
    equilibrium = solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(log_binding_weight(kd),),
    )
    expected = _decimal_one_site(p_total, l_total, kd)

    assert equilibrium.protein_free == pytest.approx(expected[0], rel=2.0e-14, abs=0.0)
    assert equilibrium.ligand_free == pytest.approx(expected[1], rel=2.0e-14, abs=0.0)
    assert equilibrium.bound_total == pytest.approx(expected[2], rel=2.0e-14, abs=0.0)
    assert math.fsum((equilibrium.protein_free, equilibrium.bound_total)) == (
        pytest.approx(p_total, rel=2.0e-15, abs=0.0)
    )
    assert math.fsum((equilibrium.ligand_free, equilibrium.bound_total)) == (
        pytest.approx(l_total, rel=2.0e-15, abs=0.0)
    )
    assert math.fsum(equilibrium.populations) == pytest.approx(1.0, rel=2.0e-15)


@pytest.mark.parametrize(
    ("p_total", "l_total", "kd"),
    (
        (1.0e-320, 1.0e308, 1.0e308),
        (1.0e-3, math.nextafter(1.0e-3, math.inf), 1.0e-300),
    ),
    ids=("extreme-asymmetric-totals", "adjacent-nearly-equal-totals"),
)
def test_finite_pool_binding_preserves_decisive_input_information(
    p_total: float,
    l_total: float,
    kd: float,
) -> None:
    expected = _decimal_one_site(p_total, l_total, kd)
    equilibrium = solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(log_binding_weight(kd),),
    )

    assert equilibrium.protein_free == pytest.approx(
        expected[0], rel=2.0e-13, abs=math.nextafter(0.0, 1.0)
    )
    assert equilibrium.ligand_free == pytest.approx(
        expected[1], rel=2.0e-13, abs=math.nextafter(0.0, 1.0)
    )
    assert equilibrium.bound_total == pytest.approx(
        expected[2], rel=2.0e-13, abs=math.nextafter(0.0, 1.0)
    )


def test_population_survives_dimensional_complex_underflow() -> None:
    minimum = math.nextafter(0.0, 1.0)
    equilibrium = solve_binding_equilibrium(
        1.0e-3,
        2.0e-3,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(
            log_binding_weight(minimum),
            log_binding_weight(1.0),
        ),
    )

    assert equilibrium.complexes[1] == 0.0
    assert equilibrium.populations[2] == minimum
    assert math.fsum(equilibrium.populations) == 1.0
    assert equilibrium.log_populations[2] > -math.inf


@pytest.mark.parametrize(
    ("total_rate", "ratio", "expected"),
    (
        (10.0, 0.0, (0.0, 10.0)),
        (10.0, 1.0, (5.0, 5.0)),
        (1.0, math.nextafter(0.0, 1.0), (math.nextafter(0.0, 1.0), 1.0)),
        (1.0, sys.float_info.max, (1.0, 1.0 / sys.float_info.max)),
    ),
)
def test_exchange_rate_split_preserves_boundaries_and_ratio(
    total_rate: float,
    ratio: float,
    expected: tuple[float, float],
) -> None:
    forward, reverse = split_exchange_rate(
        total_rate,
        0.0,
        log_equilibrium_ratio(ratio),
    )

    assert forward == pytest.approx(expected[0], rel=2.0e-15, abs=0.0)
    assert reverse == pytest.approx(expected[1], rel=3.0e-14, abs=0.0)
    assert forward + reverse == pytest.approx(total_rate, rel=2.0e-15)
    if ratio > 0.0:
        assert forward / reverse == pytest.approx(ratio, rel=3.0e-14, abs=0.0)


@pytest.mark.parametrize(
    ("total_rate", "ratio"),
    (
        (0.5, math.nextafter(0.0, 1.0)),
        (math.nextafter(0.0, 1.0), 1.0),
    ),
)
def test_exchange_rate_split_rejects_positive_unrepresentable_direction(
    total_rate: float,
    ratio: float,
) -> None:
    with pytest.raises(ValueError, match="below binary64 representability"):
        split_exchange_rate(
            total_rate,
            0.0,
            log_equilibrium_ratio(ratio),
        )


def test_alternative_bound_modes_use_normalized_equilibrium_weights() -> None:
    kd_ab = 4.0e-4
    kd_ac = 1.3e-3
    equilibrium = solve_binding_equilibrium(
        8.0e-4,
        2.3e-3,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(
            log_binding_weight(kd_ab),
            log_binding_weight(kd_ac),
        ),
    )
    pl1, pl2 = equilibrium.complexes

    assert pl1 / pl2 == pytest.approx(kd_ac / kd_ab, rel=2.0e-15)
    assert kd_ab * pl1 == pytest.approx(
        equilibrium.protein_free * equilibrium.ligand_free,
        rel=2.0e-14,
    )
    assert kd_ac * pl2 == pytest.approx(
        equilibrium.protein_free * equilibrium.ligand_free,
        rel=2.0e-14,
    )


def test_free_protein_conformers_contribute_to_apparent_equilibrium() -> None:
    equilibrium = solve_binding_equilibrium(
        8.0e-4,
        2.3e-3,
        free_protein_log_weights=(0.0, math.log(2.5)),
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(math.log(3.5 / 4.0e-4),),
    )
    free_a, free_b = equilibrium.free_proteins

    assert free_b / free_a == pytest.approx(2.5, rel=2.0e-15)
    assert equilibrium.log_apparent_kd == pytest.approx(math.log(4.0e-4))
    assert equilibrium.populations == pytest.approx(
        tuple(
            concentration / 8.0e-4
            for concentration in (*equilibrium.free_proteins, *equilibrium.complexes)
        ),
        rel=2.0e-15,
        abs=0.0,
    )


@pytest.mark.parametrize(
    ("free_weights", "complex_weights"),
    (
        ((1.0,), (1.0 / 4.0e-4, 1.0 / 1.3e-3)),
        ((1.0, 2.5), (1.0 / 4.0e-4, 2.5 / 1.3e-3)),
        (
            (1.0, 2.5),
            (1.0 / 4.0e-4, 2.5 / 1.3e-3, 2.5 * 1.7 / 1.3e-3),
        ),
        (
            (1.0,),
            (
                1.0 / 3.1e-3,
                2.5 / 3.1e-3,
                2.5 * 1.7 / 3.1e-3,
            ),
        ),
    ),
    ids=(
        "alternative-complexes",
        "partner-conformers",
        "partner-conformers-three-complexes",
        "three-bound-conformers",
    ),
)
def test_multi_species_binding_matches_independent_effective_kd_matrix(
    free_weights: tuple[float, ...],
    complex_weights: tuple[float, ...],
) -> None:
    p_total = 8.0e-4
    l_total = 2.3e-3
    with localcontext() as context:
        context.prec = 100
        free_total = sum(Decimal(str(value)) for value in free_weights)
        complex_total = sum(Decimal(str(value)) for value in complex_weights)
        kd_app = float(free_total / complex_total)
    expected_p, expected_l, expected_bound = _decimal_one_site(
        p_total,
        l_total,
        kd_app,
    )
    equilibrium = solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=tuple(math.log(value) for value in free_weights),
        complex_log_weights=tuple(math.log(value) for value in complex_weights),
    )

    assert equilibrium.protein_free == pytest.approx(expected_p, rel=3.0e-14)
    assert equilibrium.ligand_free == pytest.approx(expected_l, rel=3.0e-14)
    assert equilibrium.bound_total == pytest.approx(expected_bound, rel=3.0e-14)
    for concentration, weight in zip(
        equilibrium.free_ligands,
        free_weights,
        strict=True,
    ):
        assert concentration == pytest.approx(
            expected_l * weight / math.fsum(free_weights),
            rel=3.0e-14,
        )
    for concentration, weight in zip(
        equilibrium.complexes,
        complex_weights,
        strict=True,
    ):
        assert concentration == pytest.approx(
            expected_bound * weight / math.fsum(complex_weights),
            rel=3.0e-14,
        )


def test_zero_equilibrium_weight_removes_only_that_species() -> None:
    equilibrium = solve_binding_equilibrium(
        8.0e-4,
        2.3e-3,
        free_ligand_log_weights=(0.0, log_equilibrium_ratio(0.0)),
        complex_log_weights=(
            log_binding_weight(4.0e-4),
            log_equilibrium_ratio(0.0) + log_binding_weight(1.3e-3),
        ),
    )

    assert equilibrium.free_ligands[1] == 0.0
    assert equilibrium.complexes[1] == 0.0
    assert equilibrium.populations[2] == 0.0
    assert equilibrium.populations[0] + equilibrium.populations[1] == (
        pytest.approx(1.0)
    )


def test_zero_ligand_preserves_the_unbound_protein_population() -> None:
    equilibrium = solve_binding_equilibrium(
        8.0e-4,
        0.0,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(log_binding_weight(1.0e-3),),
    )

    assert equilibrium.protein_free == 8.0e-4
    assert equilibrium.ligand_free == 0.0
    assert equilibrium.bound_total == 0.0
    assert equilibrium.free_proteins == (8.0e-4,)
    assert equilibrium.free_ligands == (0.0,)
    assert equilibrium.complexes == (0.0,)
    assert equilibrium.populations == (1.0, 0.0)


def test_zero_ligand_distributes_a_free_protein_conformational_equilibrium() -> None:
    equilibrium = solve_binding_equilibrium(
        8.0e-4,
        0.0,
        free_protein_log_weights=(0.0, math.log(3.0)),
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(math.log(4.0 / 1.0e-3),),
    )

    assert equilibrium.free_proteins == pytest.approx((2.0e-4, 6.0e-4))
    assert equilibrium.populations == pytest.approx((0.25, 0.75, 0.0))
    assert equilibrium.log_populations == pytest.approx(
        (math.log(0.25), math.log(0.75), -math.inf)
    )


@pytest.mark.parametrize("kd", (0.0, -1.0, math.inf, -math.inf, math.nan))
def test_binding_kd_must_be_finite_and_strictly_positive(kd: float) -> None:
    with pytest.raises(ValueError, match="KD must be finite and strictly positive"):
        log_binding_weight(kd)


@pytest.mark.parametrize("p_total", (0.0, -1.0, math.inf, math.nan))
def test_binding_protein_total_must_be_finite_and_strictly_positive(
    p_total: float,
) -> None:
    with pytest.raises(
        ValueError,
        match="P_total must be finite and strictly positive",
    ):
        solve_binding_equilibrium(
            p_total,
            1.0e-3,
            free_ligand_log_weights=(0.0,),
            complex_log_weights=(log_binding_weight(1.0e-3),),
        )


@pytest.mark.parametrize("l_total", (-1.0, math.inf, math.nan))
def test_binding_ligand_total_must_be_finite_and_non_negative(
    l_total: float,
) -> None:
    with pytest.raises(ValueError, match="L_total must be finite and non-negative"):
        solve_binding_equilibrium(
            1.0e-3,
            l_total,
            free_ligand_log_weights=(0.0,),
            complex_log_weights=(log_binding_weight(1.0e-3),),
        )
