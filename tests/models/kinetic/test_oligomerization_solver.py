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
EQUILIBRIUM_TOLERANCE = 256.0 * sys.float_info.epsilon

DIRECT_MODELS = (dimer_model, trimer_model, tetramer_model)
SEQUENTIAL_MODELS = (dimer_trimer_model, dimer_tetramer_model)
UNSUPPORTED_P_TOTALS = (("zero-p-total", 0.0), ("non-finite-p-total", math.nan))
UNSUPPORTED_KDS = (
    ("zero-kd", 0.0),
    ("sub-threshold-kd", 0.5e-32),
    ("non-finite-kd", math.inf),
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
            (1e-1, 1e-32),
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
            (1e-6, 1e-32, 1e-3),
            (1e-3, 1e-5, 1e-3),
            (1e-2, 1e-3, 1e-32),
            (1e-1, 1e-32, 1e-32),
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
        _oligomerization.solve_oligomerization_fractions(((2, 1.0),))


@pytest.mark.parametrize("root", [math.nan, math.inf, -0.1, 1.1])
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
        _oligomerization.solve_oligomerization_fractions(((2, 1.0),))


def test_helper_rejects_reconstructed_state_that_violates_mass(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        _oligomerization,
        "brentq",
        lambda *_args, **_kwargs: (0.0, SimpleNamespace(converged=True)),
    )

    with pytest.raises(RuntimeError, match="violated mass conservation"):
        _oligomerization.solve_oligomerization_fractions(((2, 1.0),))


@pytest.mark.parametrize(
    ("model", "arguments"),
    [
        pytest.param(model, (p_total, 1e-3), id=f"{model.NAME}-{case}")
        for model in DIRECT_MODELS
        for case, p_total in UNSUPPORTED_P_TOTALS
    ]
    + [
        pytest.param(model, (1e-3, kd), id=f"{model.NAME}-{case}")
        for model in DIRECT_MODELS
        for case, kd in UNSUPPORTED_KDS
    ]
    + [
        pytest.param(model, (p_total, 1e-3, 1e-3), id=f"{model.NAME}-{case}")
        for model in SEQUENTIAL_MODELS
        for case, p_total in UNSUPPORTED_P_TOTALS
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
def test_unsupported_domain_routes_to_legacy_hybr(
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
