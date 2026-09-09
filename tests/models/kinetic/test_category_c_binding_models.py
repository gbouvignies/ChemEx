"""Equilibrium-authoritative qualifications for Category-C binding models."""

from __future__ import annotations

import math
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm

from chemex.chemex import run
from chemex.cli import build_parser
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import DefaultSetting, read_defaults
from chemex.experiments.builder import build_experiments
from chemex.models.kinetic import settings_3st_binding_cs as cs
from chemex.models.kinetic import settings_3st_binding_if as induced_fit
from chemex.models.kinetic._binding_migration import (
    LegacyBindingParameterError,
    validate_legacy_binding_defaults,
)
from chemex.optimize import uncertainty as uncertainty_module
from chemex.optimize.native_deterministic import run_native_deterministic
from chemex.parameters.parameterization import (
    ConstraintDomainError,
    NonFiniteParameterValueError,
    ParameterRole,
)
from chemex.parameters.sealed import (
    InvalidConfigurationError,
    parameter_name_from_definition,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession
from tests.models.kinetic.test_binding_models import _prepare_model

BINDING_EXAMPLE = Path(__file__).parents[3] / "examples/Combinations/2stBinding"
BINDING_CPMG = BINDING_EXAMPLE / "Experiments/cpmg_13p.toml"
BINDING_PARAMETERS = BINDING_EXAMPLE / "Parameters/params.toml"


def _write_execution_fixture(
    directory: Path,
    model_name: str,
) -> tuple[Path, tuple[Path, ...], Path, Path]:
    if model_name == "3st_binding_cs":
        kinetic_values = "KOFF_BC = 80.0\nKEQ_AB = 3.0\nKEX_AB = 250.0"
        fitted_kinetics = '"KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"'
    else:
        kinetic_values = "KOFF_AB = 80.0\nKEQ_BC = 3.0\nKEX_BC = 250.0"
        fitted_kinetics = '"KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"'
    parameters = directory / "parameters.toml"
    parameters.write_text(
        f"""[GLOBAL]
R1_A = 1.5
R2_A = 4.7
KD_APP = [2e-4, 1e-8, 1.0]
{kinetic_values}
""",
        encoding="utf-8",
    )
    method = directory / "method.toml"
    method.write_text(
        f"""FORMAT_VERSION = 2

[STEP]
INCLUDE = ["486"]
ROLES = [
  {{ FIX = ["R2_A", "R2_B", "R2_C", "DW_AB", "DW_AC"] }},
  {{ FIT = [{fitted_kinetics}] }},
]
STATISTICS = {{ MC = {{ REPLICATES = 1, SEED = 731 }}, BS = {{ REPLICATES = 1, SEED = 732 }}, MCMC = {{ STEPS = 2, BURN = 0, SEED = 733 }} }}
""",
        encoding="utf-8",
    )
    return BINDING_CPMG, (BINDING_PARAMETERS, parameters), method, directory / "Output"


def _parameter_output_line(output: str, name: str) -> str:
    matches = [
        line
        for line in output.splitlines()
        if line.startswith((f"{name} ", f'"{name},'))
    ]
    assert len(matches) == 1
    return matches[0]


def _output_value(line: str) -> float:
    return float(line.split("=", 1)[1].split("#", 1)[0])


def _output_unavailable_reason(line: str) -> str:
    marker = "error unavailable: "
    assert marker in line
    return line.split(marker, 1)[1].split(";", 1)[0]


@pytest.mark.parametrize("kex_ab", (300.0, 1.0e-10, 1.0e-50, 0.0))
@pytest.mark.parametrize("koff_bc", (75.0, 1.0e-50, 0.0))
def test_cs_equilibrium_is_independent_of_every_kinetic_scale(
    kex_ab: float,
    koff_bc: float,
) -> None:
    concentrations = cs.calculate_concentrations(5.0e-4, 8.0e-4, 8.0e-4, 2.0)
    populations = cs.calculate_populations(5.0e-4, 8.0e-4, 8.0e-4, 2.0)
    rates = cs.calculate_rates(
        5.0e-4,
        8.0e-4,
        8.0e-4,
        2.0,
        kex_ab,
        koff_bc,
    )

    assert concentrations == pytest.approx(
        {
            "a": 9.605091023733683e-05,
            "b": 1.9210182047467366e-04,
            "c": 2.118472692879895e-04,
            "l": 5.881527307120105e-04,
        },
        rel=3.0e-14,
    )
    assert populations == pytest.approx(
        {
            "pa": 0.19210182047467367,
            "pb": 0.38420364094934734,
            "pc": 0.42369453857597894,
        },
        rel=3.0e-14,
    )
    assert rates["kab"] + rates["kba"] == pytest.approx(kex_ab)
    if kex_ab == 0.0:
        assert rates["kab"] == rates["kba"] == 0.0
    if koff_bc == 0.0:
        assert rates["kbc"] == rates["kcb"] == 0.0
    assert populations["pa"] * rates["kab"] == pytest.approx(
        populations["pb"] * rates["kba"], rel=3.0e-14, abs=math.ulp(0.0)
    )
    assert populations["pb"] * rates["kbc"] == pytest.approx(
        populations["pc"] * rates["kcb"], rel=3.0e-14, abs=math.ulp(0.0)
    )


@pytest.mark.parametrize("kex_bc", (300.0, 1.0e-10, 1.0e-50, 0.0))
@pytest.mark.parametrize("koff_ab", (75.0, 1.0e-50, 0.0))
def test_if_equilibrium_is_independent_of_every_kinetic_scale(
    kex_bc: float,
    koff_ab: float,
) -> None:
    concentrations = induced_fit.calculate_concentrations(5.0e-4, 8.0e-4, 8.0e-4, 2.0)
    populations = induced_fit.calculate_populations(5.0e-4, 8.0e-4, 8.0e-4, 2.0)
    rates = induced_fit.calculate_rates(
        5.0e-4,
        8.0e-4,
        8.0e-4,
        2.0,
        kex_bc,
        koff_ab,
    )

    assert concentrations == pytest.approx(
        {
            "a": 2.881527307120105e-04,
            "b": 7.061575642932983e-05,
            "c": 1.4123151285865966e-04,
            "l": 5.881527307120105e-04,
        },
        rel=3.0e-14,
    )
    assert populations == pytest.approx(
        {
            "pa": 0.5763054614240212,
            "pb": 0.14123151285865965,
            "pc": 0.2824630257173193,
        },
        rel=3.0e-14,
    )
    assert rates["kbc"] + rates["kcb"] == pytest.approx(kex_bc)
    if kex_bc == 0.0:
        assert rates["kbc"] == rates["kcb"] == 0.0
    if koff_ab == 0.0:
        assert rates["kab"] == rates["kba"] == 0.0
    assert populations["pa"] * rates["kab"] == pytest.approx(
        populations["pb"] * rates["kba"], rel=3.0e-14, abs=math.ulp(0.0)
    )
    assert populations["pb"] * rates["kbc"] == pytest.approx(
        populations["pc"] * rates["kcb"], rel=3.0e-14, abs=math.ulp(0.0)
    )


@pytest.mark.parametrize(
    ("model_name", "defaults", "message"),
    (
        (
            "3st_binding_cs",
            {"KAB": 30.0, "KBA": 270.0},
            "KEQ_AB = 0.1111111111111111; KEX_AB = 300.0",
        ),
        (
            "3st_binding_if",
            {"KBC": 30.0, "KCB": 270.0},
            "KEQ_BC = 0.1111111111111111; KEX_BC = 300.0",
        ),
        (
            "3st_binding_cs",
            {"KAB": 30.0},
            "pair is incomplete",
        ),
        (
            "3st_binding_if",
            {"KBC": 0.0, "KCB": 0.0},
            "KEX_BC = 0; KEQ_BC cannot be inferred",
        ),
        (
            "3st_binding_if",
            {"KBC": 30.0, "KCB": 0.0},
            "finite KEQ_BC cannot be reconstructed",
        ),
        (
            "3st_binding_if",
            {"KBC": 0.0, "KCB": 30.0},
            "KEQ_BC = 0.0; KEX_BC = 30.0",
        ),
        (
            "3st_binding_cs",
            {"KAB": 0.0, "KBA": 30.0},
            "3st_binding_cs requires finite positive KEQ_AB",
        ),
        (
            "3st_binding_cs",
            {"KAB": 30.0, "KBA": 270.0, "KEQ_AB": 0.1},
            "old and new parameterizations cannot be combined",
        ),
    ),
)
def test_legacy_directional_parameter_defaults_are_rejected_with_migration(
    model_name: str,
    defaults: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(LegacyBindingParameterError, match=message):
        _prepare_model(model_name, 1.0e-3, 2.0e-3, defaults)


@pytest.mark.parametrize(
    ("model_name", "rates", "expected", "forbidden"),
    (
        (
            "3st_binding_if",
            (math.ulp(0.0), sys.float_info.max),
            "equivalent KEQ_BC is below minimum positive binary64",
            "KEQ_BC = 0.0",
        ),
        (
            "3st_binding_if",
            (sys.float_info.max, math.ulp(0.0)),
            "equivalent KEQ_BC exceeds maximum finite binary64",
            "KEQ_BC = inf",
        ),
        (
            "3st_binding_if",
            (sys.float_info.max, sys.float_info.max),
            "equivalent KEX_BC exceeds maximum finite binary64",
            "KEX_BC = inf",
        ),
        (
            "3st_binding_if",
            (sys.float_info.min, 1.0),
            f"KEQ_BC = {sys.float_info.min!r}",
            "outside binary64 representability",
        ),
    ),
)
def test_legacy_migration_never_invents_unrepresentable_values(
    model_name: str,
    rates: tuple[float, float],
    expected: str,
    forbidden: str,
) -> None:
    with pytest.raises(LegacyBindingParameterError) as error:
        _prepare_model(
            model_name,
            1.0e-3,
            2.0e-3,
            {"KBC": rates[0], "KCB": rates[1]},
        )

    assert expected in str(error.value)
    assert forbidden not in str(error.value)


@pytest.mark.parametrize(
    ("rates", "expected"),
    (
        (
            (999_999.0, 0.5),
            "KEQ_BC = [1999998.0, 0.0, 1999998.0]",
        ),
        (
            (600_000.0, 600_000.0),
            "KEX_BC = [1200000.0, 0.0, 1200000.0]",
        ),
    ),
)
def test_legacy_migration_explains_explicit_default_bound_override(
    rates: tuple[float, float],
    expected: str,
) -> None:
    with pytest.raises(
        LegacyBindingParameterError, match="outside the new default"
    ) as error:
        _prepare_model(
            "3st_binding_if",
            1.0e-3,
            2.0e-3,
            {"KBC": rates[0], "KCB": rates[1]},
        )

    assert expected in str(error.value)


def test_legacy_rate_pairs_are_diagnosed_per_normalized_condition_scope(
    tmp_path: Path,
) -> None:
    scoped = tmp_path / "scoped.toml"
    scoped.write_text(
        """[GLOBAL]
"kab, T->25.0C" = 30.0
"kba, T->35.0C" = 270.0
""",
        encoding="utf-8",
    )

    with pytest.raises(LegacyBindingParameterError) as error:
        validate_legacy_binding_defaults(
            "3st_binding_cs",
            read_defaults([scoped]),
        )

    message = str(error.value)
    assert "T->25.0C" in message and "KBA is also needed" in message
    assert "T->35.0C" in message and "KAB is also needed" in message
    assert "KEQ_AB = 0.1111111111111111" not in message


def test_complete_legacy_pairs_are_converted_independently_per_temperature(
    tmp_path: Path,
) -> None:
    scoped = tmp_path / "serialized-constrained.toml"
    scoped.write_text(
        """[GLOBAL]
"KBC, T->25.0C" = 30.0
"KCB, T->25.0C" = 270.0
"KBC, T->35.0C" = 50.0
"KCB, T->35.0C" = 50.0
""",
        encoding="utf-8",
    )

    with pytest.raises(LegacyBindingParameterError) as error:
        validate_legacy_binding_defaults(
            "3st_binding_if",
            read_defaults([scoped]),
        )

    message = str(error.value)
    assert "T->25.0C" in message and "KEQ_BC = 0.1111111111111111" in message
    assert "T->35.0C" in message and "KEQ_BC = 1.0" in message


@pytest.mark.parametrize(
    ("model_name", "method_entry", "legacy_name"),
    (
        ("3st_binding_cs", 'FIT = ["KAB"]', "KAB"),
        ("3st_binding_cs", 'FIX = ["KBA"]', "KBA"),
        ("3st_binding_if", 'CONSTRAINTS = ["[KBC] = [KEX_BC]"]', "KBC"),
        ("3st_binding_if", 'CONSTRAINTS = ["[KEQ_BC] = [KCB]"]', "KCB"),
        ("3st_binding_if", 'GRID = ["[KCB] = lin(1.0, 2.0, 2)"]', "KCB"),
        (
            "3st_binding_cs",
            '[STEP.SEARCH.DE]\nSEED = 19\nCOORDINATES = ["[KAB] = lin(1, 2)"]',
            "KAB",
        ),
    ),
)
def test_legacy_directional_method_roles_are_rejected_with_migration(
    tmp_path,
    model_name: str,
    method_entry: str,
    legacy_name: str,
) -> None:
    session, _ = _prepare_model(model_name, 1.0e-3, 2.0e-3, {})
    assert session.try_build_analysis_values()
    method_path = tmp_path / "method.toml"
    version = "FORMAT_VERSION = 2\n" if "SEARCH.DE" in method_entry else ""
    method_path.write_text(
        f"{version}[STEP]\n{method_entry}\n",
        encoding="utf-8",
    )
    plan = read_method_plan([method_path])

    with pytest.raises(ValueError, match=f"{legacy_name} is a derived output"):
        session.validate_method_plan(plan)


def test_stale_constrained_parameter_file_is_rejected_before_defaults_apply(
    tmp_path: Path,
) -> None:
    stale = tmp_path / "constrained.toml"
    stale.write_text("[GLOBAL]\nKAB = 30.0\nKBA = 270.0\n", encoding="utf-8")
    session, _ = _prepare_model("3st_binding_cs", 1.0e-3, 2.0e-3, {})

    with pytest.raises(
        LegacyBindingParameterError,
        match=r"KEQ_AB = 0\.1111111111111111; KEX_AB = 300\.0",
    ):
        session.parameters.set_defaults(read_defaults([stale]))


def test_old_and_new_method_syntax_has_no_precedence(tmp_path: Path) -> None:
    session, _ = _prepare_model("3st_binding_if", 1.0e-3, 2.0e-3, {})
    assert session.try_build_analysis_values()
    method_path = tmp_path / "method.toml"
    method_path.write_text(
        '[STEP]\nFIT = ["KBC", "KEQ_BC"]\n',
        encoding="utf-8",
    )
    plan = read_method_plan([method_path])

    with pytest.raises(
        ValueError, match="old and new parameterizations cannot be combined"
    ):
        session.validate_method_plan(plan)


def test_legacy_method_selectors_are_diagnosed_per_condition_scope(
    tmp_path: Path,
) -> None:
    session, _ = _prepare_model("3st_binding_cs", 1.0e-3, 2.0e-3, {})
    assert session.try_build_analysis_values()
    method_path = tmp_path / "method.toml"
    method_path.write_text(
        '[STEP]\nFIT = ["KAB, T->25.0C", "KBA, T->35.0C"]\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError) as error:
        session.validate_method_plan(read_method_plan([method_path]))

    message = str(error.value)
    assert "T->25.0C" in message and "KBA is also needed" in message
    assert "T->35.0C" in message and "KAB is also needed" in message


ORDINARY_CASES = (
    ("balanced", 1.0e-3, 2.0e-3, 1.0e-3, 100.0, 100.0, 75.0),
    ("biased", 1.0e-3, 2.0e-3, 2.0e-4, 900.0, 10.0, 33.0),
    ("slow", 5.0e-4, 8.0e-4, 8.0e-4, 1.0e-4, 2.0e-4, 0.03),
    ("fast", 5.0e-4, 8.0e-4, 8.0e-4, 4.0e5, 2.0e5, 3.0e5),
    ("weak", 2.0e-4, 3.0e-4, 0.7, 40.0, 160.0, 12.0),
    ("strong", 2.0e-4, 3.0e-4, 1.0e-12, 40.0, 160.0, 12.0),
    ("depletion", 1.0e-3, 8.0e-4, 1.0e-7, 30.0, 270.0, 45.0),
)


def _decimal_macrostate(
    p_total: float,
    l_total: float,
    kd_app: float,
) -> tuple[float, float, float]:
    with localcontext() as context:
        context.prec = 100
        protein, ligand, kd = (
            Decimal.from_float(value) for value in (p_total, l_total, kd_app)
        )
        total = protein + ligand + kd
        discriminant = (total * total - Decimal(4) * protein * ligand).sqrt()
        bound = Decimal(2) * protein * ligand / (total + discriminant)
        return float(protein - bound), float(ligand - bound), float(bound)


@pytest.mark.parametrize(
    ("_case", "p_total", "l_total", "kd_app", "legacy_kab", "legacy_kba", "koff_bc"),
    ORDINARY_CASES,
)
def test_cs_positive_legacy_parameterization_round_trips_to_equilibrium_authority(
    _case: str,
    p_total: float,
    l_total: float,
    kd_app: float,
    legacy_kab: float,
    legacy_kba: float,
    koff_bc: float,
) -> None:
    keq_ab = legacy_kab / legacy_kba
    kex_ab = legacy_kab + legacy_kba
    expected_unbound, expected_ligand, expected_bound = _decimal_macrostate(
        p_total, l_total, kd_app
    )
    expected_a = expected_unbound / (1.0 + keq_ab)
    expected_b = expected_unbound - expected_a
    kd_bc = kd_app * keq_ab / (1.0 + keq_ab)
    expected_kbc = koff_bc * expected_ligand / kd_bc

    concentrations = cs.calculate_concentrations(p_total, l_total, kd_app, keq_ab)
    populations = cs.calculate_populations(p_total, l_total, kd_app, keq_ab)
    rates = cs.calculate_rates(p_total, l_total, kd_app, keq_ab, kex_ab, koff_bc)

    assert concentrations == pytest.approx(
        {"a": expected_a, "b": expected_b, "c": expected_bound, "l": expected_ligand},
        rel=4.0e-14,
        abs=math.ulp(0.0),
    )
    assert populations == pytest.approx(
        {
            "pa": expected_a / p_total,
            "pb": expected_b / p_total,
            "pc": expected_bound / p_total,
        },
        rel=4.0e-14,
        abs=math.ulp(0.0),
    )
    assert rates == pytest.approx(
        {
            "kab": legacy_kab,
            "kba": legacy_kba,
            "kbc": expected_kbc,
            "kcb": koff_bc,
        },
        rel=5.0e-14,
        abs=math.ulp(0.0),
    )
    assert cs.calculate_intrinsic_kd(kd_app, keq_ab) == pytest.approx(kd_bc)
    assert cs.calculate_kon(koff_bc, kd_app, keq_ab) == pytest.approx(koff_bc / kd_bc)


@pytest.mark.parametrize(
    ("_case", "p_total", "l_total", "kd_app", "legacy_kbc", "legacy_kcb", "koff_ab"),
    ORDINARY_CASES,
)
def test_if_positive_legacy_parameterization_round_trips_to_equilibrium_authority(
    _case: str,
    p_total: float,
    l_total: float,
    kd_app: float,
    legacy_kbc: float,
    legacy_kcb: float,
    koff_ab: float,
) -> None:
    keq_bc = legacy_kbc / legacy_kcb
    kex_bc = legacy_kbc + legacy_kcb
    expected_free, expected_ligand, expected_bound = _decimal_macrostate(
        p_total, l_total, kd_app
    )
    expected_b = expected_bound / (1.0 + keq_bc)
    expected_c = expected_bound - expected_b
    kd_ab = kd_app * (1.0 + keq_bc)
    expected_kab = koff_ab * expected_ligand / kd_ab

    concentrations = induced_fit.calculate_concentrations(
        p_total, l_total, kd_app, keq_bc
    )
    populations = induced_fit.calculate_populations(p_total, l_total, kd_app, keq_bc)
    rates = induced_fit.calculate_rates(
        p_total, l_total, kd_app, keq_bc, kex_bc, koff_ab
    )

    assert concentrations == pytest.approx(
        {"a": expected_free, "b": expected_b, "c": expected_c, "l": expected_ligand},
        rel=4.0e-14,
        abs=math.ulp(0.0),
    )
    assert populations == pytest.approx(
        {
            "pa": expected_free / p_total,
            "pb": expected_b / p_total,
            "pc": expected_c / p_total,
        },
        rel=4.0e-14,
        abs=math.ulp(0.0),
    )
    assert rates == pytest.approx(
        {
            "kab": expected_kab,
            "kba": koff_ab,
            "kbc": legacy_kbc,
            "kcb": legacy_kcb,
        },
        rel=5.0e-14,
        abs=math.ulp(0.0),
    )
    assert induced_fit.calculate_intrinsic_kd(kd_app, keq_bc) == pytest.approx(kd_ab)
    assert induced_fit.calculate_kon(koff_ab, kd_app, keq_bc) == pytest.approx(
        koff_ab / kd_ab
    )


@pytest.mark.parametrize(
    ("module", "expected"),
    (
        (
            cs,
            (0.8331367642181426, 0.13577169691101548, 0.031091538870841925),
        ),
        (
            induced_fit,
            (0.9060016771131519, 0.08542529780633287, 0.008573025080515152),
        ),
    ),
)
def test_balanced_legacy_rates_preserve_a_frozen_propagated_observable(
    module,
    expected: tuple[float, float, float],
) -> None:
    rates = module.calculate_rates(1.0e-3, 2.0e-3, 1.0e-3, 1.0, 200.0, 75.0)
    kab, kba, kbc, kcb = (rates[name] for name in ("kab", "kba", "kbc", "kcb"))
    generator = np.array(
        (
            (-kab, kba, 0.0),
            (kab, -kba - kbc, kcb),
            (0.0, kbc, -kcb),
        )
    )
    propagated = expm(0.002 * generator) @ np.array((1.0, 0.0, 0.0))

    assert propagated == pytest.approx(expected, rel=2.0e-15, abs=0.0)


def test_cs_zero_ligand_retains_the_apo_equilibrium() -> None:
    concentrations = cs.calculate_concentrations(8.0e-4, 0.0, 4.0e-4, 3.0)
    populations = cs.calculate_populations(8.0e-4, 0.0, 4.0e-4, 3.0)
    rates = cs.calculate_rates(8.0e-4, 0.0, 4.0e-4, 3.0, 700.0, 80.0)

    assert concentrations == pytest.approx(
        {"a": 2.0e-4, "b": 6.0e-4, "c": 0.0, "l": 0.0}
    )
    assert populations == pytest.approx({"pa": 0.25, "pb": 0.75, "pc": 0.0})
    assert rates == pytest.approx({"kab": 525.0, "kba": 175.0, "kbc": 0.0, "kcb": 80.0})


def test_cs_minimum_positive_keq_retains_a_positive_forward_rate() -> None:
    minimum_positive = math.ulp(0.0)
    rates = cs.calculate_conformational_rates(1.0, minimum_positive)

    assert rates["kab"] == minimum_positive
    assert rates["kba"] > 0.0
    assert math.fsum(rates.values()) == 1.0


def test_if_zero_ligand_is_entirely_free_and_zero_keq_removes_c() -> None:
    zero_ligand = induced_fit.calculate_populations(8.0e-4, 0.0, 4.0e-4, 3.0)
    zero_ligand_rates = induced_fit.calculate_rates(
        8.0e-4, 0.0, 4.0e-4, 3.0, 700.0, 80.0
    )
    populations = induced_fit.calculate_populations(8.0e-4, 2.3e-3, 4.0e-4, 0.0)
    rates = induced_fit.calculate_rates(8.0e-4, 2.3e-3, 4.0e-4, 0.0, 700.0, 80.0)

    assert zero_ligand == {"pa": 1.0, "pb": 0.0, "pc": 0.0}
    assert zero_ligand_rates["kab"] == 0.0
    assert populations["pc"] == 0.0
    assert rates["kbc"] == 0.0
    assert rates["kcb"] == 700.0
    assert induced_fit.calculate_intrinsic_kd(4.0e-4, 0.0) == pytest.approx(4.0e-4)


@pytest.mark.parametrize("model_name", ("3st_binding_cs", "3st_binding_if"))
def test_category_c_public_independent_parameters_exclude_directional_rates(
    model_name: str,
) -> None:
    session, local_ids = _prepare_model(model_name, 1.0e-3, 2.0e-3, {})
    assert session.try_build_analysis_values()
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    independent = {
        parameter_model.definitions[param_id].name
        for param_id, declaration in parameter_model.declarations.items()
        if declaration.requires_independent
    }
    expected = (
        {"KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"}
        if model_name == "3st_binding_cs"
        else {"KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"}
    )

    assert expected <= independent
    assert {"KAB", "KBA", "KBC", "KCB"}.isdisjoint(independent - expected)
    resolved = session.resolve_current_values(
        {local_ids[name] for name in ("kab", "kba", "kbc", "kcb", "pa", "pb", "pc")}
    )
    assert all(math.isfinite(value) for value in resolved.values())


@pytest.mark.parametrize(
    ("module", "keq"),
    ((cs, math.ulp(0.0)), (cs, 1.0e6), (induced_fit, 0.0), (induced_fit, 1.0e6)),
)
@pytest.mark.parametrize(
    ("p_total", "l_total", "kd_app"),
    (
        (1.0e-300, 2.0e-300, math.ulp(0.0)),
        (1.0e-3, 1.0e-3, 1.0e-100),
        (1.0e-9, 2.0e-3, 1.0),
        (2.0e-3, 1.0e-9, 1.0e-31),
    ),
)
def test_category_c_extreme_equilibria_are_nonnegative_and_conservative(
    module,
    keq: float,
    p_total: float,
    l_total: float,
    kd_app: float,
) -> None:
    concentrations = module.calculate_concentrations(p_total, l_total, kd_app, keq)
    populations = module.calculate_populations(p_total, l_total, kd_app, keq)

    assert all(
        math.isfinite(value) and value >= 0.0 for value in concentrations.values()
    )
    assert math.fsum(populations.values()) == pytest.approx(1.0, rel=3.0e-14)
    assert math.fsum(concentrations[name] for name in ("a", "b", "c")) == (
        pytest.approx(p_total, rel=3.0e-14, abs=math.ulp(0.0))
    )
    bound_names = ("c",) if module is cs else ("b", "c")
    assert concentrations["l"] + math.fsum(
        concentrations[name] for name in bound_names
    ) == (pytest.approx(l_total, rel=3.0e-14, abs=math.ulp(0.0)))


def test_cs_zero_keq_is_rejected_at_metadata_and_runtime() -> None:
    session, _ = _prepare_model("3st_binding_cs", 1.0e-3, 2.0e-3, {"keq_ab": 0.0})
    assert not session.try_build_analysis_values()
    assert isinstance(
        session.parameter_factory.native_construction_error, InvalidConfigurationError
    )

    overridden, _ = _prepare_model(
        "3st_binding_cs",
        1.0e-3,
        2.0e-3,
        {"keq_ab": DefaultSetting(0.0, min=0.0, max=1.0e6)},
    )
    assert not overridden.try_build_analysis_values()
    assert isinstance(
        overridden.parameter_factory.native_construction_error, ConstraintDomainError
    )


@pytest.mark.parametrize("model_name", ("3st_binding_cs", "3st_binding_if"))
def test_category_c_zero_kd_and_zero_protein_are_rejected(model_name: str) -> None:
    zero_kd, _ = _prepare_model(model_name, 1.0e-3, 2.0e-3, {"kd_app": 0.0})
    assert not zero_kd.try_build_analysis_values()
    assert isinstance(
        zero_kd.parameter_factory.native_construction_error, InvalidConfigurationError
    )

    zero_protein, _ = _prepare_model(model_name, 0.0, 2.0e-3, {})
    assert not zero_protein.try_build_analysis_values()
    assert isinstance(
        zero_protein.parameter_factory.native_construction_error, ConstraintDomainError
    )


def test_unrepresentable_cs_intrinsic_reports_do_not_block_tagged_dynamics() -> None:
    session, local_ids = _prepare_model(
        "3st_binding_cs",
        1.0e-300,
        2.0e-300,
        {
            "kd_app": math.ulp(0.0),
            "keq_ab": math.ulp(0.0),
            "kex_ab": 0.0,
            "koff_bc": 0.0,
        },
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    dynamics = {
        local_ids[name] for name in ("pa", "pb", "pc", "kab", "kba", "kbc", "kcb")
    }
    assert all(
        math.isfinite(value)
        for value in session.resolve_current_values(dynamics).values()
    )
    reports = session.resolve_report_only_values()
    assert local_ids["kd_bc"] not in reports
    assert reports[local_ids["kon_bc"]] == 0.0
    with pytest.raises(NonFiniteParameterValueError):
        session.resolve_current_values({local_ids["kd_bc"]})


def test_unrepresentable_if_kon_does_not_block_finite_tagged_dynamics() -> None:
    session, local_ids = _prepare_model(
        "3st_binding_if",
        1.0e-300,
        2.0e-300,
        {
            "kd_app": math.ulp(0.0),
            "keq_bc": 0.0,
            "kex_bc": 0.0,
            "koff_ab": 1.0,
        },
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    rates = session.resolve_current_values(
        {local_ids[name] for name in ("kab", "kba", "kbc", "kcb")}
    )
    assert all(math.isfinite(value) for value in rates.values())
    reports = session.resolve_report_only_values()
    assert local_ids["kon_ab"] not in reports
    assert reports[local_ids["kd_ab"]] == math.ulp(0.0)
    with pytest.raises(NonFiniteParameterValueError):
        session.resolve_current_values({local_ids["kon_ab"]})


def test_unrepresentable_if_intrinsic_kd_is_report_only() -> None:
    maximum = sys.float_info.max
    session, local_ids = _prepare_model(
        "3st_binding_if",
        1.0e-3,
        2.0e-3,
        {
            "kd_app": DefaultSetting(maximum, min=math.ulp(0.0), max=maximum),
            "keq_bc": DefaultSetting(maximum, min=0.0, max=maximum),
            "kex_bc": 0.0,
            "koff_ab": 0.0,
        },
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    populations = session.resolve_current_values(
        {local_ids[name] for name in ("pa", "pb", "pc")}
    )
    assert all(math.isfinite(value) for value in populations.values())
    reports = session.resolve_report_only_values()
    assert local_ids["kd_ab"] not in reports
    assert reports[local_ids["kon_ab"]] == 0.0
    with pytest.raises(NonFiniteParameterValueError):
        session.resolve_current_values({local_ids["kd_ab"]})


@pytest.mark.parametrize("model_name", ("3st_binding_cs", "3st_binding_if"))
@pytest.mark.filterwarnings("ignore:invalid value encountered in divide:RuntimeWarning")
def test_category_c_parameterization_executes_fit_mcmc_and_resampling(
    tmp_path: Path,
    model_name: str,
) -> None:
    experiment, parameter_files, method, output = _write_execution_fixture(
        tmp_path,
        model_name,
    )
    args = build_parser().parse_args(
        [
            "fit",
            "-e",
            str(experiment),
            "-p",
            *(str(path) for path in parameter_files),
            "-m",
            str(method),
            "-d",
            model_name,
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )

    session = AnalysisSession.create()
    run(args, session=session)

    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "KD_APP" in fitted
    assert all(name not in fitted for name in ("KAB", "KBA", "KBC", "KCB"))
    assert all(name in constrained for name in ("KAB", "KBA", "KBC", "KCB"))
    report_only_names = (
        ("KD_BC", "KON_BC") if model_name == "3st_binding_cs" else ("KD_AB", "KON_AB")
    )
    directional_name = "KAB" if model_name == "3st_binding_cs" else "KBC"
    for name in (*report_only_names, directional_name, "PA"):
        line = _parameter_output_line(constrained, name)
        assert math.isfinite(_output_value(line))
        assert "# ±" in line
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    kex_name = "KEX_AB" if model_name == "3st_binding_cs" else "KEX_BC"
    [kex_id] = [
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name == kex_name
    ]
    fitted_kex = session.analysis_values.snapshot()[kex_id]
    assert math.isfinite(fitted_kex)
    if model_name == "3st_binding_if":
        # This fitted boundary coordinate varies slightly across platforms. The
        # invariant is positive near-zero exchange; state-0 tests cover exact zero.
        assert 0.0 < fitted_kex <= 1.0e-9
    else:
        assert fitted_kex >= 0.0
    assert (output / "Statistics" / "MonteCarlo" / "summary.toml").is_file()
    assert (output / "Statistics" / "Bootstrap" / "summary.toml").is_file()
    assert (output / "Statistics" / "MCMC" / "summary.toml").is_file()
    restart = output / "run_info" / "restart.toml"
    restart_defaults = read_defaults([restart])
    validate_legacy_binding_defaults(model_name, restart_defaults)
    restart_names = {name.name for name, _setting in restart_defaults}
    new_names = (
        {"KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"}
        if model_name == "3st_binding_cs"
        else {"KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"}
    )
    legacy_names = {"KAB", "KBA"} if model_name == "3st_binding_cs" else {"KBC", "KCB"}
    assert new_names <= restart_names
    assert restart_names.isdisjoint(legacy_names)

    reloaded_output = tmp_path / "Reloaded"
    reload_args = build_parser().parse_args(
        [
            "simulate",
            "-e",
            str(experiment),
            "-p",
            str(restart),
            "-d",
            model_name,
            "-o",
            str(reloaded_output),
            "--plot",
            "nothing",
        ]
    )
    run(reload_args, session=AnalysisSession.create())
    assert (reloaded_output / "Parameters" / "constrained.toml").is_file()

    with pytest.raises(LegacyBindingParameterError, match="derived outputs"):
        validate_legacy_binding_defaults(
            model_name,
            read_defaults([output / "Parameters" / "constrained.toml"]),
        )


@pytest.mark.parametrize("model_name", ("3st_binding_cs", "3st_binding_if"))
def test_category_c_no_variable_output_uses_reportable_parameter_set(
    tmp_path: Path,
    model_name: str,
) -> None:
    experiment_path, parameter_files, _method, output = _write_execution_fixture(
        tmp_path,
        model_name,
    )
    session = AnalysisSession.create()
    session.set_model(model_name)
    experiments = build_experiments(
        [experiment_path],
        Selection(include=[SpinSystem.from_name("486N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults(parameter_files))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    baseline = session.compile_parameterization(Method(), experiments.param_ids)
    fixed_names = sorted(
        {
            parameter_name_from_definition(parameter_model.definitions[param_id]).name
            for param_id in baseline.independent_ids
        }
    )
    parameterization = session.compile_parameterization(
        Method(fix=fixed_names),
        experiments.param_ids,
    )
    assert all(
        parameterization.role(param_id) is not ParameterRole.FIT
        for param_id in parameterization.independent_ids
    )

    fit = run_native_deterministic(
        experiments,
        output,
        "nothing",
        session=session,
        parameterization=parameterization,
    )

    assert fit is None
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    report_only_names = (
        ("KD_BC", "KON_BC") if model_name == "3st_binding_cs" else ("KD_AB", "KON_AB")
    )
    for name in (*report_only_names, "PA"):
        line = _parameter_output_line(constrained, name)
        assert math.isfinite(_output_value(line))
        assert "# ±" not in line
        assert "error unavailable" not in line


def test_category_c_report_only_output_explains_unavailable_uncertainty(
    tmp_path: Path,
) -> None:
    experiment, parameter_files, method, output = _write_execution_fixture(
        tmp_path,
        "3st_binding_cs",
    )
    method.write_text(
        method.read_text(encoding="utf-8").split("STATISTICS", 1)[0],
        encoding="utf-8",
    )
    original_svd = uncertainty_module.svd

    def rank_deficient_svd(*args, **kwargs):
        left, singular, right = original_svd(*args, **kwargs)
        singular = np.array(singular, copy=True)
        singular[-1] = 0.0
        return left, singular, right

    args = build_parser().parse_args(
        [
            "fit",
            "-e",
            str(experiment),
            "-p",
            *(str(path) for path in parameter_files),
            "-m",
            str(method),
            "-d",
            "3st_binding_cs",
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )

    with patch("chemex.optimize.uncertainty.svd", side_effect=rank_deficient_svd):
        run(args, session=AnalysisSession.create())

    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    for name in ("KD_BC", "KON_BC"):
        line = _parameter_output_line(constrained, name)
        assert math.isfinite(_output_value(line))
        assert _output_unavailable_reason(line) == (
            "constrained propagation unavailable"
        )
