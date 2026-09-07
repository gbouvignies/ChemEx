"""Scientific invariants for oligomerization kinetic models."""

from __future__ import annotations

from collections.abc import Mapping
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.nmr.basis import Basis
from chemex.parameters.name import ParamName
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

# SciPy's concentration roots are stable well inside 1e-9 relative for these
# ordinary positive fixtures across the supported Python/platform matrix.
RELATIVE_TOLERANCE = 1.0e-9
ABSOLUTE_TOLERANCE = 1.0e-14
EQUILIBRIUM_ABSOLUTE_TOLERANCE = 1.0e-24


def _construct_and_resolve(
    model_name: str,
    p_total: float,
    defaults: Mapping[str, float],
) -> tuple[dict[str, str], dict[str, str], dict[str, float]]:
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
    resolved = session.resolve_current_values(set(local_ids.values()))

    return (
        name_map,
        local_ids,
        {name: resolved[param_id] for name, param_id in local_ids.items()},
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
    assert values["kon"] == pytest.approx(values["koff"] / values["kd"])
    assert values["kd"] * oligomer_concentration == pytest.approx(
        values["c_monomer"] ** stoichiometry,
        rel=RELATIVE_TOLERANCE,
        abs=EQUILIBRIUM_ABSOLUTE_TOLERANCE,
    )
    assert values["kab"] == pytest.approx(
        stoichiometry * values["kon"] * values["c_monomer"] ** (stoichiometry - 1),
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
            "kon": 30740.74074074074,
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

    assert values["kon1"] == pytest.approx(values["koff1"] / values["kd1"])
    assert values["kon2"] == pytest.approx(values["koff2"] / values["kd2"])
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
        2.0 * values["kon1"] * values["c_monomer"],
    )
    assert values["kba"] == pytest.approx(values["koff1"])
    assert values["kac"] == pytest.approx(values["kon2"] * values["c_dimer"])
    assert values["kca"] == pytest.approx(values["koff2"] / 3.0)
    assert values["kbc"] == pytest.approx(values["kon2"] * values["c_monomer"])
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
            "kon1": 42941.17647058824,
            "kon2": 236507.9365079365,
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

    assert values["kon1"] == pytest.approx(values["koff1"] / values["kd1"])
    assert values["kon2"] == pytest.approx(values["koff2"] / values["kd2"])
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
        2.0 * values["kon1"] * values["c_monomer"],
    )
    assert values["kba"] == pytest.approx(values["koff1"])
    assert values["kbc"] == pytest.approx(
        2.0 * values["kon2"] * values["c_dimer"],
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
            "kon1": 41578.94736842105,
            "kon2": 397560.97560975613,
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
