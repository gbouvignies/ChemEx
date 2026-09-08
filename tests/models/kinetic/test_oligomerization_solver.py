from __future__ import annotations

import math
import sys
from collections.abc import Callable
from types import ModuleType, SimpleNamespace

import pytest

from chemex.models.kinetic import _oligomerization
from chemex.models.kinetic import settings_2st_monomer_dimer as dimer_model
from chemex.models.kinetic import settings_2st_monomer_tetramer as tetramer_model
from chemex.models.kinetic import settings_2st_monomer_trimer as trimer_model
from chemex.models.kinetic import (
    settings_3st_monomer_dimer_tetramer as dimer_tetramer_model,
)
from chemex.models.kinetic import (
    settings_3st_monomer_dimer_trimer as dimer_trimer_model,
)
from chemex.models.kinetic.settings_2st_monomer_dimer import (
    calculate_concentrations as calculate_dimer,
)
from chemex.models.kinetic.settings_2st_monomer_tetramer import (
    calculate_concentrations as calculate_tetramer,
)
from chemex.models.kinetic.settings_2st_monomer_trimer import (
    calculate_concentrations as calculate_trimer,
)
from chemex.models.kinetic.settings_3st_monomer_dimer_tetramer import (
    calculate_concentrations as calculate_dimer_tetramer,
)
from chemex.models.kinetic.settings_3st_monomer_dimer_trimer import (
    calculate_concentrations as calculate_dimer_trimer,
)

ConcentrationSolver = Callable[..., dict[str, float]]

MASS_TOLERANCE = 64.0 * sys.float_info.epsilon
EQUILIBRIUM_TOLERANCE = 4096.0 * sys.float_info.epsilon

DIRECT_MODELS = (dimer_model, trimer_model, tetramer_model)
SEQUENTIAL_MODELS = (dimer_trimer_model, dimer_tetramer_model)
UNSUPPORTED_P_TOTALS = (("zero-p-total", 0.0), ("non-finite-p-total", math.nan))
UNSUPPORTED_KDS = (
    ("zero-kd", 0.0),
    ("negative-kd", -1.0),
    ("non-finite-kd", math.inf),
)
POSITIVE_KD_STRESS = (
    1e-31,
    1e-32,
    1e-33,
    1e-50,
    1e-100,
    sys.float_info.min,
    math.nextafter(0.0, 1.0),
)


@pytest.mark.parametrize(
    ("solver", "arguments", "expected"),
    [
        (
            calculate_dimer,
            (1e-3, 1e-20),
            {
                "monomer": 2.2360679749997897e-12,
                "dimer": 4.999999988819660e-4,
            },
        ),
        (
            calculate_trimer,
            (1e-3, 1e-16),
            {
                "monomer": 3.217952700630540e-7,
                "trimer": 3.332260682433123e-4,
            },
        ),
        (
            calculate_tetramer,
            (1e-3, 1e-19),
            {
                "monomer": 2.2348176281143595e-6,
                "tetramer": 2.494412955929714e-4,
            },
        ),
        (
            calculate_dimer_trimer,
            (1e-3, 1e-5, 1e-3),
            {
                "monomer": 6.524662842584409e-5,
                "dimer": 4.257122520940166e-4,
                "trimer": 2.777628912870757e-5,
            },
        ),
        (
            calculate_dimer_tetramer,
            (1e-3, 1e-4, 1e-4),
            {
                "monomer": 1.1227497354437556e-4,
                "dimer": 1.2605669684390233e-4,
                "tetramer": 1.5890290819195495e-4,
            },
        ),
        (
            calculate_dimer_tetramer,
            (1e-3, 1e-8, 1e-8),
            {
                "monomer": 1.2564002054890984e-7,
                "dimer": 1.5785414763530487e-6,
                "tetramer": 2.4917931925668625e-4,
            },
        ),
    ],
    ids=[
        "dimer-failed-hybr",
        "trimer-failed-hybr",
        "tetramer-failed-hybr",
        "dimer-trimer-negative-branch",
        "dimer-tetramer-negative-branch",
        "dimer-tetramer-successful-wrong-branch",
    ],
)
def test_historical_hybr_failures_use_physical_solution(
    solver: ConcentrationSolver,
    arguments: tuple[float, ...],
    expected: dict[str, float],
) -> None:
    assert solver(*arguments) == pytest.approx(expected, rel=2e-14, abs=0.0)


@pytest.mark.parametrize(
    ("solver", "oligomer", "stoichiometry", "p_total", "kd"),
    [
        (solver, oligomer, stoichiometry, p_total, kd)
        for solver, oligomer, stoichiometry in (
            (calculate_dimer, "dimer", 2),
            (calculate_trimer, "trimer", 3),
            (calculate_tetramer, "tetramer", 4),
        )
        for p_total, kd in (
            (1e-9, 1.0),
            (1e-6, 1e-6),
            (1e-3, 1e-12),
            *((1e-100, kd) for kd in POSITIVE_KD_STRESS),
        )
    ],
)
def test_direct_positive_domain_stress(
    solver: ConcentrationSolver,
    oligomer: str,
    stoichiometry: int,
    p_total: float,
    kd: float,
) -> None:
    concentrations = solver(p_total, kd)
    monomer = concentrations["monomer"]
    oligomer_concentration = concentrations[oligomer]

    _assert_physical_concentrations(concentrations, p_total)
    _assert_normalized_mass(
        p_total,
        ((1, monomer), (stoichiometry, oligomer_concentration)),
    )
    _assert_scale_aware_close(
        kd * oligomer_concentration,
        monomer**stoichiometry,
    )


@pytest.mark.parametrize(
    ("solver", "higher_oligomer", "p_total", "kd1", "kd2"),
    [
        (solver, higher_oligomer, p_total, kd1, kd2)
        for solver, higher_oligomer in (
            (calculate_dimer_trimer, "trimer"),
            (calculate_dimer_tetramer, "tetramer"),
        )
        for p_total, kd1, kd2 in (
            (1e-9, 1.0, 1.0),
            (1e-3, 1e-5, 1e-3),
            *((1e-100, kd, 1e-6) for kd in POSITIVE_KD_STRESS),
            *((1e-100, 1e-6, kd) for kd in POSITIVE_KD_STRESS),
            *((1e-100, kd, kd) for kd in POSITIVE_KD_STRESS),
        )
    ],
)
def test_sequential_positive_domain_stress(
    solver: ConcentrationSolver,
    higher_oligomer: str,
    p_total: float,
    kd1: float,
    kd2: float,
) -> None:
    concentrations = solver(p_total, kd1, kd2)
    monomer = concentrations["monomer"]
    dimer = concentrations["dimer"]
    higher = concentrations[higher_oligomer]
    higher_stoichiometry = 3 if higher_oligomer == "trimer" else 4

    _assert_physical_concentrations(concentrations, p_total)
    _assert_normalized_mass(
        p_total,
        ((1, monomer), (2, dimer), (higher_stoichiometry, higher)),
    )
    _assert_scale_aware_close(kd1 * dimer, monomer**2)
    higher_equilibrium = monomer * dimer if higher_oligomer == "trimer" else dimer**2
    _assert_scale_aware_close(kd2 * higher, higher_equilibrium)


def test_helper_rejects_non_converged_solver_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        _oligomerization,
        "brentq",
        lambda *_args, **_kwargs: (0.5, SimpleNamespace(converged=False)),
    )

    with pytest.raises(RuntimeError, match="did not converge"):
        _oligomerization.solve_oligomerization_fractions(((2, 0.0),))


@pytest.mark.parametrize("root", [math.nan, math.inf, 0.1])
def test_helper_rejects_invalid_root(
    monkeypatch: pytest.MonkeyPatch,
    root: float,
) -> None:
    monkeypatch.setattr(
        _oligomerization,
        "brentq",
        lambda *_args, **_kwargs: (root, SimpleNamespace(converged=True)),
    )

    with pytest.raises(RuntimeError, match="invalid monomer fraction"):
        _oligomerization.solve_oligomerization_fractions(((2, 0.0),))


def test_helper_rejects_reconstructed_state_that_violates_mass(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        _oligomerization,
        "brentq",
        lambda *_args, **_kwargs: (0.0, SimpleNamespace(converged=True)),
    )

    with pytest.raises(RuntimeError, match="violated mass conservation"):
        _oligomerization.solve_oligomerization_fractions(((2, 0.0),))


@pytest.mark.parametrize(
    ("model", "arguments"),
    [
        pytest.param(model, (p_total, 1e-3), id=f"{model.NAME}-{case}")
        for model in DIRECT_MODELS
        for case, p_total in UNSUPPORTED_P_TOTALS
    ]
    + [
        pytest.param(model, (p_total, 1e-3, 1e-3), id=f"{model.NAME}-{case}")
        for model in SEQUENTIAL_MODELS
        for case, p_total in UNSUPPORTED_P_TOTALS
    ],
)
def test_unsupported_p_total_retains_legacy_hybr_policy(
    monkeypatch: pytest.MonkeyPatch,
    model: ModuleType,
    arguments: tuple[float, ...],
) -> None:
    legacy_call_count = 0
    sentinel = tuple(float(index) for index in range(len(arguments)))

    def legacy_root(*_args: object, **_kwargs: object) -> dict[str, tuple[float, ...]]:
        nonlocal legacy_call_count
        legacy_call_count += 1
        return {"x": sentinel}

    calculator = model.calculate_concentrations
    calculator.cache_clear()
    try:
        monkeypatch.setattr(model, "root", legacy_root)
        calculator(*arguments)
        assert legacy_call_count == 1
    finally:
        calculator.cache_clear()


@pytest.mark.parametrize(
    ("model", "arguments"),
    [
        pytest.param(model, (1e-3, kd), id=f"{model.NAME}-{case}")
        for model in DIRECT_MODELS
        for case, kd in UNSUPPORTED_KDS
    ]
    + [
        pytest.param(model, (1e-3, kd, 1e-3), id=f"{model.NAME}-kd1-{case}")
        for model in SEQUENTIAL_MODELS
        for case, kd in UNSUPPORTED_KDS
    ]
    + [
        pytest.param(model, (1e-3, 1e-3, kd), id=f"{model.NAME}-kd2-{case}")
        for model in SEQUENTIAL_MODELS
        for case, kd in UNSUPPORTED_KDS
    ],
)
def test_every_oligomerization_kd_position_rejects_invalid_domain(
    model: ModuleType,
    arguments: tuple[float, ...],
) -> None:
    model.calculate_concentrations.cache_clear()
    try:
        with pytest.raises(ValueError, match="KD must be finite and strictly positive"):
            model.calculate_concentrations(*arguments)
    finally:
        model.calculate_concentrations.cache_clear()


@pytest.mark.parametrize("kd", POSITIVE_KD_STRESS)
@pytest.mark.parametrize("model", (*DIRECT_MODELS, *SEQUENTIAL_MODELS))
def test_positive_kd_never_dispatches_to_legacy_hybr(
    monkeypatch: pytest.MonkeyPatch,
    model: ModuleType,
    kd: float,
) -> None:
    def fail_legacy(*_args: object, **_kwargs: object) -> None:
        pytest.fail("positive KD dispatched to legacy HYBR")

    arguments = (1e-100, kd) if model in DIRECT_MODELS else (1e-100, kd, kd)
    model.calculate_concentrations.cache_clear()
    try:
        monkeypatch.setattr(model, "root", fail_legacy)
        model.calculate_concentrations(*arguments)
    finally:
        model.calculate_concentrations.cache_clear()


@pytest.mark.parametrize("p_total", [-1.0, math.inf, math.nan])
@pytest.mark.parametrize("model", (*DIRECT_MODELS, *SEQUENTIAL_MODELS))
def test_rate_functions_retain_negative_and_non_finite_total_rejection(
    model: ModuleType,
    p_total: float,
) -> None:
    arguments = (
        (p_total, 1e-3, 1.0)
        if model in DIRECT_MODELS
        else (p_total, 1e-3, 1e-3, 1.0, 1.0)
    )
    model.calculate_rates.cache_clear()
    try:
        with pytest.raises(
            ValueError, match="P_total must be finite and strictly positive"
        ):
            model.calculate_rates(*arguments)
    finally:
        model.calculate_rates.cache_clear()


def _assert_physical_concentrations(
    concentrations: dict[str, float],
    p_total: float,
) -> None:
    assert all(math.isfinite(value) for value in concentrations.values())
    assert all(value >= 0.0 for value in concentrations.values())
    assert 0.0 <= concentrations["monomer"] <= p_total


def _assert_normalized_mass(
    p_total: float,
    species: tuple[tuple[int, float], ...],
) -> None:
    normalized_mass = math.fsum(
        stoichiometry * concentration / p_total
        for stoichiometry, concentration in species
    )
    assert abs(normalized_mass - 1.0) <= MASS_TOLERANCE


def _assert_scale_aware_close(left: float, right: float) -> None:
    scale = max(abs(left), abs(right), sys.float_info.min)
    assert abs(left - right) <= EQUILIBRIUM_TOLERANCE * scale


@pytest.mark.parametrize("kd", [1e-33, 1e-50])
def test_direct_dimer_uses_literal_sub_floor_kd(kd: float) -> None:
    concentrations = calculate_dimer(1e-3, kd)

    assert concentrations["monomer"] ** 2 / concentrations["dimer"] == (
        pytest.approx(kd, rel=1e-12, abs=0.0)
    )


def test_sub_floor_values_are_not_treated_as_legacy_effective_kd() -> None:
    baseline = calculate_dimer(1e-3, 1e-32)
    below = calculate_dimer(1e-3, 1e-33)
    far_below = calculate_dimer(1e-3, 1e-50)

    assert below["monomer"] / baseline["monomer"] == pytest.approx(
        math.sqrt(0.1),
        rel=1e-12,
    )
    assert far_below["monomer"] / baseline["monomer"] == pytest.approx(
        1e-9,
        rel=1e-12,
    )
    assert dimer_model.calculate_rates(1e-3, 1e-33, 1.0)["kab"] != (
        pytest.approx(
            dimer_model.calculate_rates(1e-3, 1e-32, 1.0)["kab"],
            rel=1e-3,
        )
    )


@pytest.mark.parametrize(
    ("model", "stoichiometry"),
    ((dimer_model, 2), (trimer_model, 3), (tetramer_model, 4)),
)
@pytest.mark.parametrize("kd", POSITIVE_KD_STRESS)
def test_direct_stress_equilibrium_is_satisfied_in_log_space(
    model: ModuleType,
    stoichiometry: int,
    kd: float,
) -> None:
    p_total = 1e-100
    equilibrium = model._calculate_equilibrium(p_total, kd)
    log_monomer = math.log(p_total) + equilibrium.log_monomer_fraction
    log_oligomer = math.log(p_total) + equilibrium.log_oligomer_fractions[0]

    assert math.log(kd) + log_oligomer == pytest.approx(
        stoichiometry * log_monomer,
        abs=2e-12,
    )


@pytest.mark.parametrize(
    ("model", "higher_stoichiometry"),
    ((dimer_trimer_model, 3), (dimer_tetramer_model, 4)),
)
@pytest.mark.parametrize("vary", ("kd1", "kd2", "both"))
@pytest.mark.parametrize("kd", POSITIVE_KD_STRESS)
def test_sequential_stress_equilibria_are_satisfied_in_log_space(
    model: ModuleType,
    higher_stoichiometry: int,
    vary: str,
    kd: float,
) -> None:
    p_total = 1e-100
    kd1 = kd if vary in {"kd1", "both"} else 1e-6
    kd2 = kd if vary in {"kd2", "both"} else 1e-6
    equilibrium = model._calculate_equilibrium(p_total, kd1, kd2)
    log_p_total = math.log(p_total)
    log_monomer = log_p_total + equilibrium.log_monomer_fraction
    log_dimer = log_p_total + equilibrium.log_oligomer_fractions[0]
    log_higher = log_p_total + equilibrium.log_oligomer_fractions[1]

    assert math.log(kd1) + log_dimer == pytest.approx(
        2.0 * log_monomer,
        abs=2e-12,
    )
    right = log_monomer + log_dimer if higher_stoichiometry == 3 else 2.0 * log_dimer
    assert math.log(kd2) + log_higher == pytest.approx(right, abs=2e-12)


def test_minimum_subnormal_kd_does_not_fail_on_coefficient_overflow() -> None:
    kd = math.nextafter(0.0, 1.0)
    concentrations = calculate_tetramer(1e-3, kd)

    _assert_physical_concentrations(concentrations, 1e-3)
    _assert_normalized_mass(
        1e-3,
        ((1, concentrations["monomer"]), (4, concentrations["tetramer"])),
    )


@pytest.mark.parametrize("kd", [0.0, -1.0, math.inf, math.nan])
def test_shared_helper_rejects_non_positive_or_non_finite_kd(kd: float) -> None:
    with pytest.raises(ValueError, match="KD must be finite and strictly positive"):
        _oligomerization.validate_oligomerization_kd(kd)


def test_log_equilibrium_retains_underflowed_fraction_authority() -> None:
    log_coefficient = _oligomerization.log_equilibrium_coefficient(
        1e-100,
        3,
        1e-31,
    )
    equilibrium = _oligomerization.solve_oligomerization_fractions(
        ((4, log_coefficient),)
    )

    assert equilibrium.oligomer_fractions[0] > 0.0
    assert 1e-100 * equilibrium.oligomer_fractions[0] == 0.0
    assert math.isfinite(equilibrium.log_oligomer_fractions[0])
    assert (
        math.log(1e-31) + equilibrium.log_oligomer_fractions[0] + math.log(1e-100)
    ) == pytest.approx(
        4.0 * (math.log(1e-100) + equilibrium.log_monomer_fraction),
        abs=2e-12,
    )


def test_concentration_reconstruction_keeps_representable_absolute_value() -> None:
    log_fraction = -750.0

    assert math.exp(log_fraction) == 0.0
    concentration = _oligomerization.concentration_from_log_fraction(
        1e100,
        log_fraction,
    )

    assert concentration > 0.0
    assert concentration == pytest.approx(
        math.exp(math.log(1e100) + log_fraction),
        rel=2e-13,
    )


def test_concentration_reconstruction_keeps_representable_fraction() -> None:
    assert _oligomerization.concentration_from_log_fraction(
        1e3,
        math.log(0.25),
    ) == pytest.approx(250.0)


def test_concentration_reconstruction_returns_zero_only_below_representability() -> (
    None
):
    assert _oligomerization.concentration_from_log_fraction(1e-100, -750.0) == 0.0


def test_dimer_trimer_reconstructs_audited_subnormal_trimer_concentration() -> None:
    concentrations = calculate_dimer_trimer(
        1e100,
        math.nextafter(0.0, 1.0),
        sys.float_info.max,
    )

    # Rounded binary64 value of the independently evaluated high-precision
    # equilibrium concentration 4.3715130080390679...e-321.
    assert concentrations["trimer"] == float.fromhex("0x0.0000000000375p-1022")


def test_koff_over_kd_may_overflow_when_direct_final_rate_is_finite() -> None:
    minimum = math.nextafter(0.0, 1.0)

    assert math.isinf(1.0 / minimum)
    rate = dimer_model.calculate_rates(minimum, minimum, 1.0)["kab"]
    assert math.isfinite(rate) and rate > 0.0


def test_koff_over_kd_may_overflow_when_sequential_final_rates_are_finite() -> None:
    minimum = math.nextafter(0.0, 1.0)

    assert math.isinf(1.0 / minimum)
    rates = dimer_trimer_model.calculate_rates(
        minimum,
        minimum,
        minimum,
        1.0,
        1.0,
    )
    assert all(math.isfinite(rate) and rate > 0.0 for rate in rates.values())


def test_population_ratio_may_overflow_when_log_product_rate_is_finite() -> None:
    reverse_rate = 1e-308
    log_source = -710.0
    log_destination = 0.0

    with pytest.raises(OverflowError):
        math.exp(log_destination - log_source)
    forward_rate = _oligomerization.detailed_balance_forward_rate(
        reverse_rate,
        log_source,
        log_destination,
    )

    assert forward_rate == pytest.approx(
        math.exp(math.log(reverse_rate) + 710.0),
        rel=2e-13,
    )


def test_actual_positive_forward_rate_overflow_is_rejected() -> None:
    with pytest.raises(ValueError, match="exceeds the maximum finite binary64"):
        _oligomerization.detailed_balance_forward_rate(
            sys.float_info.max,
            -1.0,
            0.0,
        )


def test_actual_positive_forward_rate_underflow_is_rejected() -> None:
    with pytest.raises(ValueError, match="below binary64 representability"):
        _oligomerization.detailed_balance_forward_rate(
            math.nextafter(0.0, 1.0),
            0.0,
            -1.0,
        )


@pytest.mark.parametrize("log_rate", [-745.5, -744.5])
def test_positive_forward_rate_below_exact_minimum_is_rejected_before_rounding(
    log_rate: float,
) -> None:
    assert log_rate < _oligomerization.LOG_MIN_POSITIVE_FLOAT

    with pytest.raises(ValueError, match="below binary64 representability"):
        _oligomerization.detailed_balance_forward_rate(1.0, 0.0, log_rate)


@pytest.mark.parametrize(
    "rate",
    [
        math.nextafter(0.0, 1.0),
        math.nextafter(math.nextafter(0.0, 1.0), math.inf),
    ],
)
def test_positive_forward_rate_at_or_above_exact_minimum_is_accepted(
    rate: float,
) -> None:
    assert _oligomerization.detailed_balance_forward_rate(rate, 0.0, 0.0) == rate


def test_scaled_positive_rate_uses_exact_minimum_boundary() -> None:
    minimum = math.nextafter(0.0, 1.0)
    next_above = math.nextafter(minimum, math.inf)

    for factor in (0.25, 0.75):
        with pytest.raises(ValueError, match="below binary64 representability"):
            _oligomerization.scale_reversible_rate(minimum, factor)
    assert _oligomerization.scale_reversible_rate(minimum, 1.0) == minimum
    assert _oligomerization.scale_reversible_rate(next_above, 1.0) == next_above


def test_exact_zero_reverse_rate_returns_exact_zero_forward_rate() -> None:
    assert _oligomerization.detailed_balance_forward_rate(0.0, -1000.0, 0.0) == 0.0
