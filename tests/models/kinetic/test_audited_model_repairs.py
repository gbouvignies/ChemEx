"""Regression tests for definite kinetic-model audit findings."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.nmr.basis import Basis
from chemex.parameters.name import ParamName
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

# The concentration helpers use SciPy's nonlinear root solver. Repeated runs agree
# below 1e-12 relative; 1e-10 leaves headroom for supported-platform solver noise.
# The 1e-14 absolute tolerance applies only near exact-zero limits; relative
# tolerance dominates the order-one algebraic H/D population comparisons.
RELATIVE_TOLERANCE = 1.0e-10
ABSOLUTE_TOLERANCE = 1.0e-14
PARTNER_CONDITIONS = Conditions(
    h_larmor_frq=600.0,
    temperature=25.0,
    p_total=8.0e-4,
    l_total=2.3e-3,
)


def _construct_and_resolve(
    model_name: str,
    conditions: Conditions,
    defaults: Mapping[str, float],
    resolved_names: set[str],
) -> tuple[dict[str, str], dict[str, float]]:
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    config = SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    settings = model_factory.create(model_name, conditions)
    local_ids = {
        name: setting.name_setting.get_param_name(spin_system, conditions).id_
        for name, setting in settings.items()
    }

    name_map = session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=spin_system,
    )
    session.parameters.set_defaults(
        [
            (ParamName.from_section(name), DefaultSetting(value))
            for name, value in defaults.items()
        ],
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error,
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error,
    )
    resolved = session.resolve_current_values(
        {local_ids[name] for name in resolved_names},
    )

    return name_map, {name: resolved[local_ids[name]] for name in resolved_names}


def test_2st_hd_exact_zero_kdh_preserves_equilibrium_composition() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        d2o=0.25,
    )

    _, values = _construct_and_resolve(
        "2st_hd",
        conditions,
        {"kdh": 0.0, "phi": 1.2},
        {"kab", "kba", "pa", "pb"},
    )

    assert values == pytest.approx(
        {
            "kab": 0.0,
            "kba": 0.0,
            "pa": 5.0 / 7.0,
            "pb": 2.0 / 7.0,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_2st_hd_populations_are_invariant_to_kdh_scale() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        d2o=0.25,
    )
    expected_populations = {"pa": 5.0 / 7.0, "pb": 2.0 / 7.0}

    for kdh in (0.0, 1.0e-20, 1.0, 1.0e6):
        _, values = _construct_and_resolve(
            "2st_hd",
            conditions,
            {"kdh": kdh, "phi": 1.2},
            {"kab", "kba", "pa", "pb"},
        )

        assert {name: values[name] for name in ("pa", "pb")} == pytest.approx(
            expected_populations,
            rel=RELATIVE_TOLERANCE,
            abs=ABSOLUTE_TOLERANCE,
        )
        if kdh > 0.0:
            rate_sum = values["kab"] + values["kba"]
            assert values["pa"] == pytest.approx(
                values["kba"] / rate_sum,
                rel=RELATIVE_TOLERANCE,
                abs=ABSOLUTE_TOLERANCE,
            )
            assert values["pb"] == pytest.approx(
                values["kab"] / rate_sum,
                rel=RELATIVE_TOLERANCE,
                abs=ABSOLUTE_TOLERANCE,
            )


def test_2st_hd_preserves_ordinary_positive_rate_behavior() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        d2o=0.2,
    )

    _, values = _construct_and_resolve(
        "2st_hd",
        conditions,
        {"kdh": 10.0, "phi": 1.25},
        {"kab", "kba", "pa", "pb"},
    )

    assert values == pytest.approx(
        {
            "kab": 2.5,
            "kba": 8.0,
            "pa": 8.0 / 10.5,
            "pb": 2.5 / 10.5,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


@pytest.mark.parametrize(
    ("d2o", "phi", "expected_populations"),
    [
        (0.0, 0.75, {"pa": 1.0, "pb": 0.0}),
        (1.0, 1.5, {"pa": 0.0, "pb": 1.0}),
    ],
)
def test_2st_hd_exact_zero_kdh_preserves_solvent_endpoints(
    d2o: float,
    phi: float,
    expected_populations: dict[str, float],
) -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        d2o=0.25,
    )

    _, values = _construct_and_resolve(
        "2st_hd",
        conditions,
        {"d2o": d2o, "kdh": 0.0, "phi": phi},
        {"pa", "pb"},
    )

    assert values == pytest.approx(
        expected_populations,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


@pytest.mark.parametrize(
    ("d2o", "phi"),
    [
        (0.0, 0.75),
        (1.0e-12, 1.5),
        (0.25, 1.2),
        (0.5, 1.0),
        (1.0 - 1.0e-12, 0.75),
        (1.0, 1.5),
    ],
)
def test_2st_hd_populations_satisfy_equilibrium_invariants(
    d2o: float,
    phi: float,
) -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        d2o=0.25,
    )

    _, values = _construct_and_resolve(
        "2st_hd",
        conditions,
        {"d2o": d2o, "kdh": 0.0, "phi": phi},
        {"pa", "pb"},
    )
    denominator = 1.0 + d2o * (phi - 1.0)
    expected = {
        "pa": (1.0 - d2o) / denominator,
        "pb": d2o * phi / denominator,
    }

    assert all(isfinite(values[name]) for name in ("pa", "pb"))
    assert values["pa"] >= 0.0
    assert values["pb"] >= 0.0
    assert values["pa"] + values["pb"] == pytest.approx(
        1.0,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values == pytest.approx(
        expected,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_4st_hd_constructs_the_conformer_hd_square_topology() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
        d2o=0.2,
    )
    active_rates = {"kab", "kba", "kac", "kca", "kbd", "kdb", "kcd", "kdc"}
    absent_rates = {"kad", "kda", "kbc", "kcb"}
    resolved_names = active_rates | {"pa", "pb", "pc", "pd"}

    name_map, values = _construct_and_resolve(
        "4st_hd",
        conditions,
        {
            "pop_b": 0.3,
            "kex_ab": 400.0,
            "kdh_a": 10.0,
            "kdh_b": 20.0,
            "phi_a": 1.25,
        },
        resolved_names,
    )

    assert active_rates <= name_map.keys()
    assert absent_rates.isdisjoint(name_map)
    assert values == pytest.approx(
        {
            "kab": 120.0,
            "kba": 280.0,
            "kac": 2.5,
            "kca": 8.0,
            "kbd": 5.0,
            "kdb": 16.0,
            "kcd": 120.0,
            "kdc": 280.0,
            "pa": 0.5333333333333333,
            "pb": 0.22857142857142856,
            "pc": 0.16666666666666666,
            "pd": 0.07142857142857142,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert sum(values[name] for name in ("pa", "pb", "pc", "pd")) == (
        pytest.approx(1.0, rel=RELATIVE_TOLERANCE, abs=ABSOLUTE_TOLERANCE)
    )


def test_4st_hd_phi_b_has_repaired_constraint_metadata() -> None:
    settings = model_factory.create(
        "4st_hd",
        Conditions(
            h_larmor_frq=600.0,
            temperature=25.0,
            p_total=1.0e-3,
            l_total=2.0e-3,
            d2o=0.2,
        ),
    )

    phi_b = settings["phi_b"]
    assert phi_b.value == 1.1
    assert phi_b.min == 0.75
    assert phi_b.max == 1.50
    assert phi_b.expr == "{phi_a}"


def test_4st_hd_constructs_distinct_population_definitions_for_two_residues() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
        d2o=0.2,
    )
    session = AnalysisSession.create()
    session.set_model("4st_hd")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )

    name_maps = [
        session.parameter_factory.create_parameters(
            config,  # ty: ignore[invalid-argument-type]
            basis=basis,
            spin_system=SpinSystem.from_name(spin_system),
        )
        for spin_system in ("G23N-HN", "A24N-HN")
    ]
    population_ids = [
        {name_map[name] for name in ("pa", "pb", "pc", "pd")} for name_map in name_maps
    ]

    assert population_ids[0].isdisjoint(population_ids[1])
    session.parameters.set_defaults([])
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error,
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error,
    )
    resolved = session.resolve_current_values(population_ids[0] | population_ids[1])
    assert population_ids[0] | population_ids[1] <= resolved.keys()
    for name_map in name_maps:
        assert sum(resolved[name_map[name]] for name in ("pa", "pb", "pc", "pd")) == (
            pytest.approx(
                1.0,
                rel=RELATIVE_TOLERANCE,
                abs=ABSOLUTE_TOLERANCE,
            )
        )


def test_4st_hd_separates_population_definitions_by_d2o_condition() -> None:
    conditions = [
        Conditions(
            h_larmor_frq=600.0,
            temperature=25.0,
            p_total=1.0e-3,
            l_total=2.0e-3,
            d2o=d2o,
        )
        for d2o in (0.2, 0.8)
    ]
    session = AnalysisSession.create()
    session.set_model("4st_hd")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    name_maps = []

    for condition in conditions:
        config = SimpleNamespace(
            conditions=condition,
            to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
        )
        name_maps.append(
            session.parameter_factory.create_parameters(
                config,  # ty: ignore[invalid-argument-type]
                basis=basis,
                spin_system=spin_system,
            ),
        )

    population_names = [
        {
            name: session.parameters.get_parameters([name_map[name]])[
                name_map[name]
            ].param_name
            for name in ("pa", "pb", "pc", "pd")
        }
        for name_map in name_maps
    ]
    for name in ("pa", "pb", "pc", "pd"):
        assert population_names[0][name] != population_names[1][name]

    population_ids = [
        {name_map[name] for name in ("pa", "pb", "pc", "pd")} for name_map in name_maps
    ]
    all_population_ids = population_ids[0] | population_ids[1]
    assert population_ids[0].isdisjoint(population_ids[1])

    session.parameters.set_defaults([])
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error,
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error,
    )
    resolved = session.resolve_current_values(all_population_ids)
    assert all_population_ids <= resolved.keys()
    for name_map in name_maps:
        assert sum(resolved[name_map[name]] for name in ("pa", "pb", "pc", "pd")) == (
            pytest.approx(
                1.0,
                rel=RELATIVE_TOLERANCE,
                abs=ABSOLUTE_TOLERANCE,
            )
        )
    assert resolved[name_maps[0]["pa"]] > resolved[name_maps[1]["pa"]]
    assert resolved[name_maps[0]["pc"]] < resolved[name_maps[1]["pc"]]


def test_4st_hd_zero_conformational_exchange_preserves_both_hd_pairs() -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
        d2o=0.2,
    )

    _, values = _construct_and_resolve(
        "4st_hd",
        conditions,
        {
            "pop_b": 0.3,
            "kex_ab": 0.0,
            "kdh_a": 10.0,
            "kdh_b": 20.0,
            "phi_a": 1.25,
        },
        {"kab", "kba", "kcd", "kdc", "pa", "pb", "pc", "pd"},
    )

    assert values == pytest.approx(
        {
            "kab": 0.0,
            "kba": 0.0,
            "kcd": 0.0,
            "kdc": 0.0,
            "pa": 0.5333333333333333,
            "pb": 0.22857142857142856,
            "pc": 0.16666666666666666,
            "pd": 0.07142857142857142,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


@pytest.mark.parametrize("d2o", [1.0e-6, 1.0 - 1.0e-6])
def test_4st_hd_populations_follow_d2o_at_solvent_limits(d2o: float) -> None:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
        d2o=d2o,
    )

    _, values = _construct_and_resolve(
        "4st_hd",
        conditions,
        {"pop_b": 0.3, "phi_a": 1.25},
        {"pa", "pb", "pc", "pd"},
    )
    denominator = 1.0 + d2o * (1.25 - 1.0)

    assert values == pytest.approx(
        {
            "pa": 0.7 * (1.0 - d2o) / denominator,
            "pb": 0.3 * (1.0 - d2o) / denominator,
            "pc": 0.7 * d2o * 1.25 / denominator,
            "pd": 0.3 * d2o * 1.25 / denominator,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


@pytest.mark.parametrize(
    "model_name",
    ["3st_binding_partner_2st", "4st_binding_partner_2st"],
)
def test_partner_binding_models_have_unique_derived_parameter_identities(
    model_name: str,
) -> None:
    conditions = PARTNER_CONDITIONS
    settings = model_factory.create(model_name, conditions)
    spin_system = SpinSystem.from_name("G23N-HN")
    identities = [
        setting.name_setting.get_param_name(spin_system, conditions).id_
        for setting in settings.values()
    ]

    assert len(identities) == len(set(identities))


def test_3st_partner_binding_uses_distinct_asymmetric_free_species() -> None:
    conditions = PARTNER_CONDITIONS
    resolved_names = {
        "l1_free",
        "l2_free",
        "pl1",
        "pl2",
        "kab",
        "kac",
        "kbc",
        "kcb",
        "pa",
        "pb",
        "pc",
    }

    _, values = _construct_and_resolve(
        "3st_binding_partner_2st",
        conditions,
        {
            "koff_ab": 80.0,
            "kd_ab": 4.0e-4,
            "koff_ac": 130.0,
            "kd_ac": 1.3e-3,
            "keq": 2.5,
            "kex_bc": 700.0,
        },
        resolved_names,
    )

    assert values == pytest.approx(
        {
            "l1_free": 4.997662866500894e-4,
            "l2_free": 1.2494157166252236e-3,
            "pl1": 3.11331911192151e-4,
            "pl2": 2.3948608553253585e-4,
            "kab": 99.95325733001788,
            "kac": 124.94157166252236,
            "kbc": 304.347826087037,
            "kcb": 395.65217391296306,
            "pa": 0.3114775040940694,
            "pb": 0.3891648889902396,
            "pc": 0.2993576069156911,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["l1_free"] != pytest.approx(values["l2_free"])
    p_free = conditions.p_total - values["pl1"] - values["pl2"]
    assert conditions.l_total == pytest.approx(
        values["l1_free"] + values["l2_free"] + values["pl1"] + values["pl2"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert 4.0e-4 * values["pl1"] == pytest.approx(
        p_free * values["l1_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert 1.3e-3 * values["pl2"] == pytest.approx(
        p_free * values["l2_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["l2_free"] == pytest.approx(
        2.5 * values["l1_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )


def test_4st_partner_binding_uses_all_intended_free_and_bound_species() -> None:
    conditions = PARTNER_CONDITIONS
    resolved_names = {
        "p_free",
        "l1_free",
        "l2_free",
        "pl1",
        "pl2",
        "pl3",
        "kab",
        "kac",
        "kbc",
        "kcb",
        "kcd",
        "kdc",
        "pa",
        "pb",
        "pc",
        "pd",
    }

    _, values = _construct_and_resolve(
        "4st_binding_partner_2st",
        conditions,
        {
            "koff_ab": 80.0,
            "kd_ab": 4.0e-4,
            "koff_ac": 130.0,
            "kd_ac": 1.3e-3,
            "keq_l": 2.5,
            "keq_pl": 1.7,
            "kex_bc": 700.0,
            "kex_cd": 900.0,
        },
        resolved_names,
    )

    assert values == pytest.approx(
        {
            "p_free": 1.7119763210448377e-4,
            "l1_free": 4.774850377441382e-4,
            "l2_free": 1.1937125943603454e-3,
            "pl1": 2.043607695642459e-4,
            "pl2": 1.5720059197454455e-4,
            "pl3": 2.6724100635672573e-4,
            "kab": 95.49700754882764,
            "kac": 119.37125943603455,
            "kbc": 304.34782608919727,
            "kcb": 395.6521739108027,
            "kcd": 566.6666666666666,
            "kdc": 333.3333333333333,
            "pa": 0.2139970401306047,
            "pb": 0.25545096195530735,
            "pc": 0.19650073996818068,
            "pd": 0.33405125794590716,
        },
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["l1_free"] != pytest.approx(values["l2_free"])
    assert sum(values[name] for name in ("pa", "pb", "pc", "pd")) == (
        pytest.approx(1.0, rel=RELATIVE_TOLERANCE, abs=ABSOLUTE_TOLERANCE)
    )
    assert conditions.p_total == pytest.approx(
        values["p_free"] + values["pl1"] + values["pl2"] + values["pl3"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert conditions.l_total == pytest.approx(
        values["l1_free"]
        + values["l2_free"]
        + values["pl1"]
        + values["pl2"]
        + values["pl3"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert 4.0e-4 * values["pl1"] == pytest.approx(
        values["p_free"] * values["l1_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert 1.3e-3 * values["pl2"] == pytest.approx(
        values["p_free"] * values["l2_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["l2_free"] == pytest.approx(
        2.5 * values["l1_free"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
    assert values["pl3"] == pytest.approx(
        1.7 * values["pl2"],
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    )
