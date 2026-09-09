"""Scientific invariants for oligomerization kinetic models."""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.models.kinetic import (
    settings_3st_monomer_dimer_tetramer as dimer_tetramer_model,
)
from chemex.models.kinetic import (
    settings_3st_monomer_dimer_trimer as dimer_trimer_model,
)
from chemex.nmr.basis import Basis
from chemex.optimize.deterministic_uncertainty import (
    compile_model_constraint_linearization_capabilities,
)
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import ConstraintDomainError
from chemex.parameters.sealed import InvalidConfigurationError
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

# SciPy's concentration roots are stable well inside 1e-9 relative for these
# ordinary positive fixtures across the supported Python/platform matrix.
RELATIVE_TOLERANCE = 1.0e-9
ABSOLUTE_TOLERANCE = 1.0e-14
EQUILIBRIUM_ABSOLUTE_TOLERANCE = 1.0e-24


def _prepare_model(
    model_name: str,
    p_total: float,
    defaults: Mapping[str, float | DefaultSetting],
) -> tuple[AnalysisSession, dict[str, str], dict[str, str]]:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=p_total,
    )
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    settings = model_factory.create(model_name, conditions)
    local_ids = {
        name: setting.name_setting.get_param_name(spin_system, conditions).id_
        for name, setting in settings.items()
    }
    config = SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )

    name_map = session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=spin_system,
    )
    session.parameters.set_defaults(
        [
            (
                ParamName.from_section(name),
                value if isinstance(value, DefaultSetting) else DefaultSetting(value),
            )
            for name, value in defaults.items()
        ],
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error,
    )
    return session, name_map, local_ids


def _construct_and_resolve(
    model_name: str,
    p_total: float,
    defaults: Mapping[str, float | DefaultSetting],
) -> tuple[dict[str, str], dict[str, str], dict[str, float]]:
    session, name_map, local_ids = _prepare_model(model_name, p_total, defaults)
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error,
    )
    resolved = session.resolve_current_values(set(local_ids.values()))

    return (
        name_map,
        local_ids,
        {name: resolved[param_id] for name, param_id in local_ids.items()},
    )


@pytest.mark.parametrize(
    ("model_name", "population_names"),
    (
        ("2st_monomer_dimer", ("pa", "pb")),
        ("2st_monomer_trimer", ("pa", "pb")),
        ("2st_monomer_tetramer", ("pa", "pb")),
        ("3st_monomer_dimer_trimer", ("pa", "pb", "pc")),
        ("3st_monomer_dimer_tetramer", ("pa", "pb", "pc")),
    ),
)
def test_oligomer_population_outputs_have_production_uncertainty_capabilities(
    model_name: str,
    population_names: tuple[str, ...],
) -> None:
    session, _name_map, local_ids = _prepare_model(model_name, 1.0e-3, {})
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    output_ids = tuple(local_ids[name] for name in population_names)
    parameterization = session.compile_parameterization(Method(), set(output_ids))

    compiled = compile_model_constraint_linearization_capabilities(
        parameterization,
        output_ids,
    )

    assert compiled.output_scope == output_ids
    population_capabilities = tuple(
        capability
        for capability in compiled.capabilities
        if capability.component in population_names
    )
    assert {capability.component for capability in population_capabilities} == set(
        population_names
    )
    assert all(
        capability.normalized_population_components == population_names
        for capability in population_capabilities
    )


@pytest.mark.parametrize(
    ("model_name", "oligomer", "stoichiometry", "p_total", "kd", "koff"),
    [
        ("2st_monomer_dimer", "dimer", 2, 1.3e-3, 2.7e-3, 83.0),
        ("2st_monomer_trimer", "trimer", 3, 1.1e-3, 7.3e-7, 127.0),
        ("2st_monomer_tetramer", "tetramer", 4, 1.4e-3, 4.2e-10, 91.0),
    ],
)
def test_direct_oligomerization_resolves_tagged_chemical_equilibrium(
    model_name: str,
    oligomer: str,
    stoichiometry: int,
    p_total: float,
    kd: float,
    koff: float,
) -> None:
    name_map, local_ids, values = _construct_and_resolve(
        model_name,
        p_total,
        {"kd": kd, "koff": koff},
    )
    oligomer_concentration = values[f"c_{oligomer}"]
    chemical_pa = values["c_monomer"] / p_total
    chemical_pb = stoichiometry * oligomer_concentration / p_total

    assert {"kab", "kba", "pa", "pb"} <= name_map.keys()
    assert len(local_ids) == len(set(local_ids.values()))
    assert values["c_monomer"] >= 0.0
    assert oligomer_concentration >= 0.0
    assert values["c_monomer"] + stoichiometry * oligomer_concentration == (
        pytest.approx(
            p_total,
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )
    )
    assert values["pa"] + values["pb"] == pytest.approx(
        1.0,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["pa"] == pytest.approx(
        chemical_pa,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["pb"] == pytest.approx(
        chemical_pb,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["kd"] * oligomer_concentration == pytest.approx(
        values["c_monomer"] ** stoichiometry,
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kab"] == pytest.approx(
        stoichiometry
        * values["koff"]
        / values["kd"]
        * values["c_monomer"] ** (stoichiometry - 1),
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["kba"] == pytest.approx(koff)
    assert values["pa"] * values["kab"] == pytest.approx(
        values["pb"] * values["kba"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_direct_dimer_resolved_oracle_is_unchanged() -> None:
    _, _, values = _construct_and_resolve(
        "2st_monomer_dimer",
        1.3e-3,
        {"kd": 2.7e-3, "koff": 83.0},
    )

    assert values == pytest.approx(
        {
            "koff": 83.0,
            "kd": 2.7e-3,
            "c_monomer": 0.0008118170701199257,
            "c_dimer": 0.0002440914649400371,
            "kab": 49.91171616292876,
            "kba": 83.0,
            "pa": 0.6244746693230198,
            "pb": 0.37552533067698024,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_sequential_dimer_trimer_resolved_oracle_and_tagged_equilibrium() -> None:
    p_total = 1.25e-3
    name_map, local_ids, values = _construct_and_resolve(
        "3st_monomer_dimer_trimer",
        p_total,
        {"kd1": 1.7e-3, "kd2": 6.3e-4, "koff1": 73.0, "koff2": 149.0},
    )
    chemical_populations = {
        "pa": values["c_monomer"] / p_total,
        "pb": 2.0 * values["c_dimer"] / p_total,
        "pc": 3.0 * values["c_trimer"] / p_total,
    }

    assert {"kab", "kba", "kac", "kca", "kbc", "kcb"} <= name_map.keys()
    assert len(local_ids) == len(set(local_ids.values()))
    assert all(values[name] >= 0.0 for name in ("c_monomer", "c_dimer", "c_trimer"))
    assert (
        values["c_monomer"] + 2.0 * values["c_dimer"] + 3.0 * values["c_trimer"]
    ) == pytest.approx(
        p_total,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert sum(values[name] for name in ("pa", "pb", "pc")) == pytest.approx(
        1.0,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    for name, chemical_population in chemical_populations.items():
        assert values[name] == pytest.approx(
            chemical_population,
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )

    assert values["kd1"] * values["c_dimer"] == pytest.approx(
        values["c_monomer"] ** 2,
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kd2"] * values["c_trimer"] == pytest.approx(
        values["c_monomer"] * values["c_dimer"],
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kab"] == pytest.approx(
        2.0 * values["koff1"] / values["kd1"] * values["c_monomer"],
    )
    assert values["kba"] == pytest.approx(values["koff1"])
    assert values["kac"] == pytest.approx(
        values["koff2"] / values["kd2"] * values["c_dimer"]
    )
    assert values["kca"] == pytest.approx(values["koff2"] / 3.0)
    assert values["kbc"] == pytest.approx(
        values["koff2"] / values["kd2"] * values["c_monomer"]
    )
    assert values["kcb"] == pytest.approx(2.0 * values["koff2"] / 3.0)
    for left_population, forward_rate, right_population, reverse_rate in (
        ("pa", "kab", "pb", "kba"),
        ("pa", "kac", "pc", "kca"),
        ("pb", "kbc", "pc", "kcb"),
    ):
        assert values[left_population] * values[forward_rate] == pytest.approx(
            values[right_population] * values[reverse_rate],
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )

    assert values == pytest.approx(
        {
            "kd1": 1.7e-3,
            "kd2": 6.3e-4,
            "koff1": 73.0,
            "koff2": 149.0,
            "c_monomer": 0.0005239890527646985,
            "c_dimer": 0.0001615085455406397,
            "c_trimer": 0.00013433128538467405,
            "kab": 45.001412766850585,
            "kba": 73.0,
            "kac": 38.198052834214785,
            "kca": 49.666666666666664,
            "kbc": 123.92756962212711,
            "kcb": 99.33333333333333,
            "pa": 0.4191912422143663,
            "pb": 0.25841367286552414,
            "pc": 0.3223950849201095,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_sequential_dimer_tetramer_constructs_distinct_tagged_equilibrium() -> None:
    p_total = 1.2e-3
    name_map, local_ids, values = _construct_and_resolve(
        "3st_monomer_dimer_tetramer",
        p_total,
        {"kd1": 1.9e-3, "kd2": 4.1e-4, "koff1": 79.0, "koff2": 163.0},
    )
    chemical_populations = {
        "pa": values["c_monomer"] / p_total,
        "pb": 2.0 * values["c_dimer"] / p_total,
        "pc": 4.0 * values["c_tetramer"] / p_total,
    }
    distinct_names = {
        "c_monomer",
        "c_dimer",
        "c_tetramer",
        "kab",
        "kba",
        "kbc",
        "kcb",
    }

    assert {"kab", "kba", "kbc", "kcb", "pa", "pb", "pc"} <= name_map.keys()
    assert {"kac", "kca"}.isdisjoint(name_map)
    assert distinct_names <= local_ids.keys()
    assert len(local_ids) == len(set(local_ids.values()))
    assert len({local_ids[name] for name in distinct_names}) == len(distinct_names)
    for name in distinct_names:
        assert local_ids[name].startswith(f"__{name.upper()}__")

    assert all(values[name] >= 0.0 for name in ("c_monomer", "c_dimer", "c_tetramer"))
    assert (
        values["c_monomer"] + 2.0 * values["c_dimer"] + 4.0 * values["c_tetramer"]
    ) == pytest.approx(
        p_total,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert sum(values[name] for name in ("pa", "pb", "pc")) == pytest.approx(
        1.0,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    for name, chemical_population in chemical_populations.items():
        assert values[name] == pytest.approx(
            chemical_population,
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )

    assert values["kd1"] * values["c_dimer"] == pytest.approx(
        values["c_monomer"] ** 2,
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kd2"] * values["c_tetramer"] == pytest.approx(
        values["c_dimer"] ** 2,
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kab"] == pytest.approx(
        2.0 * values["koff1"] / values["kd1"] * values["c_monomer"],
    )
    assert values["kba"] == pytest.approx(values["koff1"])
    assert values["kbc"] == pytest.approx(
        2.0 * values["koff2"] / values["kd2"] * values["c_dimer"],
    )
    assert values["kcb"] == pytest.approx(values["koff2"])
    for left_population, forward_rate, right_population, reverse_rate in (
        ("pa", "kab", "pb", "kba"),
        ("pb", "kbc", "pc", "kcb"),
    ):
        assert values[left_population] * values[forward_rate] == pytest.approx(
            values[right_population] * values[reverse_rate],
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )

    assert values == pytest.approx(
        {
            "kd1": 1.9e-3,
            "kd2": 4.1e-4,
            "koff1": 79.0,
            "koff2": 163.0,
            "c_monomer": 0.0005706471522611405,
            "c_dimer": 0.00017138851179980075,
            "c_tetramer": 7.164395603481447e-5,
            "kab": 47.45381581961063,
            "kba": 79.0,
            "kbc": 136.27476791886596,
            "kcb": 163.0,
            "pa": 0.47553929356807656,
            "pb": 0.28564751964515517,
            "pc": 0.23881318678676824,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


@pytest.mark.parametrize(
    "model_name",
    (
        "2st_monomer_dimer",
        "2st_monomer_trimer",
        "2st_monomer_tetramer",
        "3st_monomer_dimer_trimer",
        "3st_monomer_dimer_tetramer",
    ),
)
def test_oligomerization_models_do_not_construct_eager_kon_settings(
    model_name: str,
) -> None:
    settings = model_factory.create(
        model_name,
        Conditions(h_larmor_frq=600.0, temperature=25.0, p_total=1e-3),
    )

    assert {"kon", "kon1", "kon2"}.isdisjoint(settings)


@pytest.mark.parametrize(
    ("model_name", "defaults", "concentration_names", "expected_rates"),
    (
        (
            "2st_monomer_dimer",
            {"kd": 1e-6, "koff": 90.0},
            ("c_monomer", "c_dimer"),
            {"kab": 0.0, "kba": 90.0},
        ),
        (
            "2st_monomer_trimer",
            {"kd": 1e-6, "koff": 90.0},
            ("c_monomer", "c_trimer"),
            {"kab": 0.0, "kba": 90.0},
        ),
        (
            "2st_monomer_tetramer",
            {"kd": 1e-6, "koff": 90.0},
            ("c_monomer", "c_tetramer"),
            {"kab": 0.0, "kba": 90.0},
        ),
        (
            "3st_monomer_dimer_trimer",
            {"kd1": 1e-6, "kd2": 2e-6, "koff1": 90.0, "koff2": 150.0},
            ("c_monomer", "c_dimer", "c_trimer"),
            {
                "kab": 0.0,
                "kba": 90.0,
                "kac": 0.0,
                "kca": 50.0,
                "kbc": 0.0,
                "kcb": 100.0,
            },
        ),
        (
            "3st_monomer_dimer_tetramer",
            {"kd1": 1e-6, "kd2": 2e-6, "koff1": 90.0, "koff2": 150.0},
            ("c_monomer", "c_dimer", "c_tetramer"),
            {"kab": 0.0, "kba": 90.0, "kbc": 0.0, "kcb": 150.0},
        ),
    ),
)
def test_zero_total_full_model_retains_zero_association_and_state_a(
    model_name: str,
    defaults: dict[str, float],
    concentration_names: tuple[str, ...],
    expected_rates: dict[str, float],
) -> None:
    _, _, values = _construct_and_resolve(model_name, 0.0, defaults)

    assert {name: values[name] for name in concentration_names} == dict.fromkeys(
        concentration_names, 0.0
    )
    assert {name: values[name] for name in expected_rates} == pytest.approx(
        expected_rates
    )
    assert values["pa"] == 1.0
    assert values["pb"] == 0.0
    if "pc" in values:
        assert values["pc"] == 0.0


KD_STRESS = (
    1e-31,
    1e-32,
    1e-33,
    1e-50,
    1e-100,
    sys.float_info.min,
    math.nextafter(0.0, 1.0),
)


@pytest.mark.parametrize(
    ("model_name", "oligomer", "stoichiometry"),
    (
        ("2st_monomer_dimer", "dimer", 2),
        ("2st_monomer_trimer", "trimer", 3),
        ("2st_monomer_tetramer", "tetramer", 4),
    ),
)
@pytest.mark.parametrize("kd", KD_STRESS)
def test_direct_positive_kd_stress_resolves_finite_detailed_balance(
    model_name: str,
    oligomer: str,
    stoichiometry: int,
    kd: float,
) -> None:
    p_total = 1e-100
    _, _, values = _construct_and_resolve(
        model_name,
        p_total,
        {"kd": kd, "koff": 1.0},
    )

    assert values["kd"] == kd
    assert all(
        math.isfinite(values[name]) and values[name] >= 0.0
        for name in ("c_monomer", f"c_{oligomer}", "kab", "kba", "pa", "pb")
    )
    assert values["c_monomer"] + stoichiometry * values[f"c_{oligomer}"] == (
        pytest.approx(p_total, rel=2e-12, abs=math.nextafter(0.0, 1.0))
    )
    assert values["pa"] * values["kab"] == pytest.approx(
        values["pb"] * values["kba"],
        rel=2e-12,
        abs=math.nextafter(0.0, 1.0),
    )


@pytest.mark.parametrize(
    ("model_name", "higher", "higher_stoichiometry"),
    (
        ("3st_monomer_dimer_trimer", "trimer", 3),
        ("3st_monomer_dimer_tetramer", "tetramer", 4),
    ),
)
@pytest.mark.parametrize("vary", ("kd1", "kd2", "both"))
@pytest.mark.parametrize("kd", KD_STRESS)
def test_sequential_positive_kd_stress_resolves_finite_detailed_balance(
    model_name: str,
    higher: str,
    higher_stoichiometry: int,
    vary: str,
    kd: float,
) -> None:
    p_total = 1e-100
    kd1 = kd if vary in {"kd1", "both"} else 1e-6
    kd2 = kd if vary in {"kd2", "both"} else 1e-6
    model = (
        dimer_trimer_model
        if model_name == "3st_monomer_dimer_trimer"
        else dimer_tetramer_model
    )
    log_pa, log_pb, log_pc = model._calculate_equilibrium(
        p_total,
        kd1,
        kd2,
    ).log_tagged_fractions
    koff1 = math.exp(min(0.0, log_pa - log_pb))
    koff2 = math.exp(min(0.0, log_pa - log_pc, log_pb - log_pc))
    _, _, values = _construct_and_resolve(
        model_name,
        p_total,
        {"kd1": kd1, "kd2": kd2, "koff1": koff1, "koff2": koff2},
    )

    assert values["kd1"] == kd1
    assert values["kd2"] == kd2
    names = {
        "c_monomer",
        "c_dimer",
        f"c_{higher}",
        "kab",
        "kba",
        "kbc",
        "kcb",
        "pa",
        "pb",
        "pc",
    }
    if higher == "trimer":
        names |= {"kac", "kca"}
    assert all(math.isfinite(values[name]) and values[name] >= 0.0 for name in names)
    assert (
        values["c_monomer"]
        + 2.0 * values["c_dimer"]
        + higher_stoichiometry * values[f"c_{higher}"]
    ) == pytest.approx(
        p_total,
        rel=2e-12,
        abs=math.nextafter(0.0, 1.0),
    )
    edges = [("pa", "kab", "pb", "kba"), ("pb", "kbc", "pc", "kcb")]
    if higher == "trimer":
        edges.append(("pa", "kac", "pc", "kca"))
    for source, forward, destination, reverse in edges:
        assert values[source] * values[forward] == pytest.approx(
            values[destination] * values[reverse],
            rel=2e-10,
            abs=ABSOLUTE_TOLERANCE,
        )


@pytest.mark.parametrize(
    ("model_name", "defaults"),
    (
        ("2st_monomer_dimer", {"kd": 0.0}),
        ("2st_monomer_trimer", {"kd": 0.0}),
        ("2st_monomer_tetramer", {"kd": 0.0}),
        ("3st_monomer_dimer_trimer", {"kd1": 0.0}),
        ("3st_monomer_dimer_trimer", {"kd2": 0.0}),
        ("3st_monomer_dimer_tetramer", {"kd1": 0.0}),
        ("3st_monomer_dimer_tetramer", {"kd2": 0.0}),
    ),
)
def test_zero_kd_is_rejected_by_default_parameter_bounds(
    model_name: str,
    defaults: dict[str, float],
) -> None:
    session, _, _ = _prepare_model(model_name, 1e-3, defaults)

    assert not session.try_build_analysis_values()
    assert isinstance(
        session.parameter_factory.native_construction_error,
        InvalidConfigurationError,
    )


@pytest.mark.parametrize(
    ("model_name", "kd_name"),
    (
        ("2st_monomer_dimer", "kd"),
        ("2st_monomer_trimer", "kd"),
        ("2st_monomer_tetramer", "kd"),
        ("3st_monomer_dimer_trimer", "kd1"),
        ("3st_monomer_dimer_trimer", "kd2"),
        ("3st_monomer_dimer_tetramer", "kd1"),
        ("3st_monomer_dimer_tetramer", "kd2"),
    ),
)
def test_zero_kd_override_reaches_runtime_domain_authority(
    model_name: str,
    kd_name: str,
) -> None:
    session, _, _ = _prepare_model(
        model_name,
        1e-3,
        {kd_name: DefaultSetting(0.0, min=0.0, max=1.0)},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "KD must be finite and strictly positive" in str(error.__cause__)


@pytest.mark.parametrize(
    ("model_name", "p_total", "defaults", "message"),
    (
        (
            "2st_monomer_dimer",
            sys.float_info.max,
            {"kd": math.nextafter(0.0, 1.0), "koff": 1.0},
            "exceeds the maximum finite binary64 value",
        ),
        (
            "2st_monomer_tetramer",
            1e-100,
            {"kd": 1e-31, "koff": 1e-100},
            "below binary64 representability",
        ),
    ),
)
def test_unrepresentable_actual_forward_rate_is_constraint_domain_error(
    model_name: str,
    p_total: float,
    defaults: dict[str, float],
    message: str,
) -> None:
    session, _, _ = _prepare_model(model_name, p_total, defaults)

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert message in str(error.__cause__)


def test_monomer_trimer_rejects_upward_rounded_subminimum_forward_rate() -> None:
    p_total = math.exp((-744.5 - math.log(3.0)) / 2.0)
    session, _, _ = _prepare_model(
        "2st_monomer_trimer",
        p_total,
        {"kd": 1.0, "koff": 1.0},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "below binary64 representability" in str(error.__cause__)


@pytest.mark.parametrize(
    ("model_name", "oligomer", "stoichiometry"),
    (
        ("2st_monomer_dimer", "dimer", 2.0),
        ("2st_monomer_trimer", "trimer", 3.0),
        ("2st_monomer_tetramer", "tetramer", 4.0),
    ),
)
def test_direct_zero_koff_freezes_exchange_without_erasing_equilibrium_species(
    model_name: str,
    oligomer: str,
    stoichiometry: float,
) -> None:
    p_total = 1.0e-3
    _, _, values = _construct_and_resolve(
        model_name,
        p_total,
        {"kd": 1.0e-6, "koff": 0.0},
    )

    assert values["kab"] == 0.0
    assert values["kba"] == 0.0
    assert values["pa"] == pytest.approx(values["c_monomer"] / p_total)
    assert values["pb"] == pytest.approx(
        stoichiometry * values[f"c_{oligomer}"] / p_total
    )
    assert values["pa"] > 0.0
    assert values["pb"] > 0.0


@pytest.mark.parametrize(
    ("model_name", "higher", "higher_stoichiometry"),
    (
        ("3st_monomer_dimer_trimer", "trimer", 3.0),
        ("3st_monomer_dimer_tetramer", "tetramer", 4.0),
    ),
)
@pytest.mark.parametrize(
    ("koff1", "koff2"),
    ((0.0, 0.0), (0.0, 150.0), (90.0, 0.0), (90.0, 150.0)),
)
def test_sequential_populations_are_invariant_to_every_koff_disconnection(
    model_name: str,
    higher: str,
    higher_stoichiometry: float,
    koff1: float,
    koff2: float,
) -> None:
    p_total = 1.0e-3
    _, _, values = _construct_and_resolve(
        model_name,
        p_total,
        {"kd1": 1.0e-6, "kd2": 2.0e-6, "koff1": koff1, "koff2": koff2},
    )

    assert values["pa"] == pytest.approx(values["c_monomer"] / p_total)
    assert values["pb"] == pytest.approx(2.0 * values["c_dimer"] / p_total)
    assert values["pc"] == pytest.approx(
        higher_stoichiometry * values[f"c_{higher}"] / p_total
    )
    assert values["pa"] > 0.0
    assert values["pb"] > 0.0
    assert values["pc"] > 0.0
