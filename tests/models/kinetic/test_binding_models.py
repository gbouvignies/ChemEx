"""Association-equilibrium authority across finite-pool binding models."""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.nmr.basis import Basis
from chemex.optimize.deterministic_uncertainty import (
    compile_model_constraint_linearization_capabilities,
)
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ConstraintDomainError,
    NonFiniteParameterValueError,
)
from chemex.parameters.sealed import InvalidConfigurationError
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.parameters import write_parameters
from chemex.runtime import AnalysisSession


def _prepare_model(
    model_name: str,
    p_total: float,
    l_total: float,
    defaults: Mapping[str, float | DefaultSetting],
) -> tuple[AnalysisSession, dict[str, str]]:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=25.0,
        p_total=p_total,
        l_total=l_total,
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
    session.parameter_factory.create_parameters(
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
    return session, local_ids


def _construct_and_resolve(
    model_name: str,
    p_total: float,
    l_total: float,
    defaults: Mapping[str, float | DefaultSetting],
) -> dict[str, float]:
    session, local_ids = _prepare_model(model_name, p_total, l_total, defaults)
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error,
    )
    resolved = session.resolve_current_values(set(local_ids.values()))
    return {name: resolved[param_id] for name, param_id in local_ids.items()}


@pytest.mark.parametrize(
    ("model_name", "population_names", "defaults"),
    (
        ("2st_binding", ("pa", "pb"), {}),
        ("3st_double_binding", ("pa", "pb", "pc"), {}),
        ("3st_binding_partner_2st", ("pa", "pb", "pc"), {}),
        ("4st_binding_partner_2st", ("pa", "pb", "pc", "pd"), {}),
        (
            "4st_binding_3_bound_states",
            ("pa", "pb", "pc", "pd", "kd_eff"),
            {},
        ),
    ),
)
def test_binding_population_outputs_have_production_uncertainty_capabilities(
    model_name: str,
    population_names: tuple[str, ...],
    defaults: Mapping[str, float],
) -> None:
    session, local_ids = _prepare_model(
        model_name,
        1.0e-3,
        2.0e-3,
        defaults,
    )
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
    expected_components = tuple(
        component for component in population_names if component != "kd_eff"
    )
    population_capabilities = tuple(
        capability
        for capability in compiled.capabilities
        if capability.component in expected_components
    )
    assert {capability.component for capability in population_capabilities} == set(
        expected_components
    )
    assert all(
        capability.normalized_population_components == expected_components
        for capability in population_capabilities
    )


@pytest.mark.parametrize(
    ("p_total", "l_total", "expected_p_free"),
    (
        (1.0e-100, 1.0e-100, 1.0e-200),
        (1.0e-3, 1.0e-3, math.sqrt(1.0e-303)),
    ),
)
def test_2st_binding_resolves_audited_strong_binding_cases(
    p_total: float,
    l_total: float,
    expected_p_free: float,
) -> None:
    values = _construct_and_resolve(
        "2st_binding",
        p_total,
        l_total,
        {"kd": 1.0e-300, "koff": 0.0},
    )

    assert values["p_free"] == pytest.approx(expected_p_free, rel=2.0e-14, abs=0.0)
    assert values["l_free"] == pytest.approx(expected_p_free, rel=2.0e-14, abs=0.0)
    assert values["pl"] > 0.999999999999 * p_total
    assert values["p_free"] + values["pl"] == pytest.approx(
        p_total, rel=2.0e-15, abs=0.0
    )
    assert values["l_free"] + values["pl"] == pytest.approx(
        l_total, rel=2.0e-15, abs=0.0
    )
    assert values["pa"] == pytest.approx(
        values["p_free"] / p_total, rel=5.0e-14, abs=0.0
    )
    assert values["pb"] == pytest.approx(values["pl"] / p_total, rel=2.0e-14, abs=0.0)
    assert values["kab"] == 0.0
    assert values["kba"] == 0.0


@pytest.mark.parametrize(
    "kd",
    (1.0e-31, 1.0e-32, 1.0e-50, 1.0e-100, 1.0e-300),
)
def test_2st_binding_preserves_literal_positive_kd_below_legacy_floor(
    kd: float,
) -> None:
    values = _construct_and_resolve(
        "2st_binding",
        1.0e-3,
        2.0e-3,
        {"kd": kd, "koff": kd},
    )

    assert values["kd"] == kd
    assert values["kon"] == pytest.approx(1.0)
    assert values["kab"] == pytest.approx(values["l_free"])
    assert values["pa"] * values["kab"] == pytest.approx(
        values["pb"] * values["kba"],
        rel=2.0e-12,
        abs=math.nextafter(0.0, 1.0),
    )


def test_2st_binding_zero_ligand_is_unbound() -> None:
    values = _construct_and_resolve(
        "2st_binding",
        1.0e-3,
        0.0,
        {"kd": 1.0e-6, "koff": 100.0},
    )

    assert values["p_free"] == 1.0e-3
    assert values["l_free"] == 0.0
    assert values["pl"] == 0.0
    assert values["kab"] == 0.0
    assert values["kba"] == 100.0
    assert values["pa"] == 1.0
    assert values["pb"] == 0.0


def test_unrepresentable_public_kon_does_not_block_finite_tagged_dynamics() -> None:
    minimum = math.nextafter(0.0, 1.0)
    session, local_ids = _prepare_model(
        "2st_binding",
        1.0e-300,
        1.0e-300,
        {"kd": minimum, "koff": 1.0},
    )

    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    dynamics = {local_ids[name] for name in ("pa", "pb", "kab", "kba")}
    resolved = session.resolve_current_values(dynamics)
    assert resolved[local_ids["kab"]] == pytest.approx(4.5e11, rel=2.0e-2)
    assert resolved[local_ids["kba"]] == 1.0
    assert resolved[local_ids["pa"]] + resolved[local_ids["pb"]] == 1.0

    with pytest.raises(NonFiniteParameterValueError):
        session.resolve_current_values({local_ids["kon"]})

    assert session.resolve_report_only_values() == {}


def test_supplied_unrepresentable_report_only_kon_remains_deferred() -> None:
    minimum = math.nextafter(0.0, 1.0)
    session, local_ids = _prepare_model(
        "2st_binding",
        1.0e-300,
        1.0e-300,
        {"kd": minimum, "koff": 1.0, "kon": 1.0},
    )

    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    snapshot = session.analysis_values.snapshot()
    assert local_ids["kon"] not in snapshot
    dynamics = {local_ids[name] for name in ("pa", "pb", "kab", "kba")}
    parameterization = session.compile_parameterization(Method(), dynamics)
    assert local_ids["kon"] not in parameterization.scope_ids
    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))
    assert all(math.isfinite(resolved[param_id]) for param_id in dynamics)
    assert session.resolve_report_only_values() == {}

    with pytest.raises(NonFiniteParameterValueError):
        session.resolve_current_values({local_ids["kon"]})


def test_representable_public_kon_remains_serialized_at_report_time(
    tmp_path: Path,
) -> None:
    session, local_ids = _prepare_model(
        "2st_binding",
        1.0e-3,
        2.0e-3,
        {"kd": 2.0e-6, "koff": 40.0, "kon": 1.0},
    )
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    assert local_ids["kon"] not in session.analysis_values.snapshot()
    dynamics = {local_ids[name] for name in ("pa", "pb", "kab", "kba")}
    parameterization = session.compile_parameterization(Method(), dynamics)
    assert local_ids["kon"] not in parameterization.scope_ids
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    assert resolved[local_ids["kab"]] == pytest.approx(
        resolved[local_ids["pb"]]
        * resolved[local_ids["kba"]]
        / resolved[local_ids["pa"]]
    )
    assert session.resolve_report_only_values()[local_ids["kon"]] == 2.0e7

    write_parameters(
        tmp_path,
        parameter_model=parameter_model,
        parameter_values=resolved,
        parameterization=parameterization,
        report_only_values=session.resolve_report_only_values(),
    )

    constrained = (tmp_path / "Parameters" / "constrained.toml").read_text()
    assert "KON" in constrained
    assert "2.00000e+07" in constrained
    assert "[KOFF, T->25.0C] / [KD, T->25.0C]" in constrained


def test_2st_binding_zero_kd_is_rejected_by_default_bounds() -> None:
    session, _ = _prepare_model(
        "2st_binding",
        1.0e-3,
        2.0e-3,
        {"kd": 0.0},
    )

    assert not session.try_build_analysis_values()
    assert isinstance(
        session.parameter_factory.native_construction_error,
        InvalidConfigurationError,
    )


def test_2st_binding_zero_kd_override_reaches_runtime_domain_authority() -> None:
    session, _ = _prepare_model(
        "2st_binding",
        1.0e-3,
        2.0e-3,
        {"kd": DefaultSetting(0.0, min=0.0, max=1.0)},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "KD must be finite and strictly positive" in str(error.__cause__)


def test_category_a_binding_rejects_zero_protein_total() -> None:
    session, _ = _prepare_model(
        "2st_binding",
        0.0,
        2.0e-3,
        {"kd": 1.0e-6, "koff": 100.0},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "P_total must be finite and strictly positive" in str(error.__cause__)


@pytest.mark.parametrize(
    ("koff_ab", "koff_ac"),
    ((0.0, 0.0), (0.0, 130.0), (80.0, 0.0), (80.0, 130.0)),
)
def test_double_binding_populations_are_invariant_to_branch_disconnection(
    koff_ab: float,
    koff_ac: float,
) -> None:
    p_total = 8.0e-4
    values = _construct_and_resolve(
        "3st_double_binding",
        p_total,
        2.3e-3,
        {
            "koff_ab": koff_ab,
            "kd_ab": 4.0e-4,
            "koff_ac": koff_ac,
            "kd_ac": 1.3e-3,
        },
    )

    assert values["pa"] + values["pb"] + values["pc"] == pytest.approx(1.0)
    assert values["pa"] == pytest.approx(values["pfree"] / p_total)
    assert all(values[name] > 0.0 for name in ("pa", "pb", "pc"))
    assert values["kd_ab"] * values["pb"] == pytest.approx(
        values["pfree"] * values["l_free"] / p_total,
        rel=2.0e-14,
    )
    assert values["kd_ac"] * values["pc"] == pytest.approx(
        values["pfree"] * values["l_free"] / p_total,
        rel=2.0e-14,
    )
    assert values["pa"] * values["kab"] == pytest.approx(values["pb"] * values["kba"])
    assert values["pa"] * values["kac"] == pytest.approx(values["pc"] * values["kca"])


def test_double_binding_strong_limit_retains_both_alternative_complexes() -> None:
    values = _construct_and_resolve(
        "3st_double_binding",
        1.0e-3,
        1.0e-3,
        {
            "koff_ab": 0.0,
            "kd_ab": 1.0e-300,
            "koff_ac": 0.0,
            "kd_ac": 1.0e-300,
        },
    )

    assert values["pfree"] == pytest.approx(math.sqrt(5.0e-304), rel=5.0e-14, abs=0.0)
    assert values["l_free"] == pytest.approx(values["pfree"], rel=5.0e-14, abs=0.0)
    assert values["pb"] == pytest.approx(0.5, rel=2.0e-14)
    assert values["pc"] == pytest.approx(0.5, rel=2.0e-14)
    assert values["pa"] > 0.0


@pytest.mark.parametrize("kd_name", ("kd_ab", "kd_ac"))
def test_double_binding_zero_kd_is_rejected_at_metadata_and_runtime(
    kd_name: str,
) -> None:
    session, _ = _prepare_model(
        "3st_double_binding",
        1.0e-3,
        2.0e-3,
        {kd_name: 0.0},
    )
    assert not session.try_build_analysis_values()
    assert isinstance(
        session.parameter_factory.native_construction_error,
        InvalidConfigurationError,
    )

    overridden, _ = _prepare_model(
        "3st_double_binding",
        1.0e-3,
        2.0e-3,
        {kd_name: DefaultSetting(0.0, min=0.0, max=1.0)},
    )
    assert not overridden.try_build_analysis_values()
    error = overridden.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "KD must be finite and strictly positive" in str(error.__cause__)


@pytest.mark.parametrize(
    ("koff_ab", "koff_ac", "kex_bc"),
    (
        (0.0, 0.0, 0.0),
        (80.0, 0.0, 0.0),
        (0.0, 130.0, 0.0),
        (0.0, 0.0, 700.0),
        (80.0, 130.0, 0.0),
        (80.0, 0.0, 700.0),
        (0.0, 130.0, 700.0),
        (80.0, 130.0, 700.0),
    ),
)
def test_partner_binding_populations_are_invariant_to_every_kinetic_scale(
    koff_ab: float,
    koff_ac: float,
    kex_bc: float,
) -> None:
    p_total = 8.0e-4
    values = _construct_and_resolve(
        "3st_binding_partner_2st",
        p_total,
        2.3e-3,
        {
            "koff_ab": koff_ab,
            "kd_ab": 4.0e-4,
            "koff_ac": koff_ac,
            "kd_ac": 1.3e-3,
            "keq": 2.5,
            "kex_bc": kex_bc,
        },
    )

    assert values["pa"] + values["pb"] + values["pc"] == pytest.approx(1.0)
    assert values["pa"] == pytest.approx(1.0 - values["pb"] - values["pc"])
    assert values["pb"] == pytest.approx(values["pl1"] / p_total)
    assert values["pc"] == pytest.approx(values["pl2"] / p_total)
    assert values["l2_free"] == pytest.approx(2.5 * values["l1_free"])
    p_free = p_total * values["pa"]
    assert 4.0e-4 * values["pl1"] == pytest.approx(
        p_free * values["l1_free"], rel=2.0e-14
    )
    assert 1.3e-3 * values["pl2"] == pytest.approx(
        p_free * values["l2_free"], rel=2.0e-14
    )
    assert values["pa"] * values["kab"] == pytest.approx(values["pb"] * values["kba"])
    assert values["pa"] * values["kac"] == pytest.approx(values["pc"] * values["kca"])
    assert values["pb"] * values["kbc"] == pytest.approx(values["pc"] * values["kcb"])


def test_partner_binding_zero_keq_removes_second_conformer_without_a_floor() -> None:
    values = _construct_and_resolve(
        "3st_binding_partner_2st",
        8.0e-4,
        2.3e-3,
        {
            "koff_ab": 80.0,
            "kd_ab": 4.0e-4,
            "koff_ac": 130.0,
            "kd_ac": 1.3e-3,
            "keq": 0.0,
            "kex_bc": 700.0,
        },
    )

    assert values["l2_free"] == 0.0
    assert values["pl2"] == 0.0
    assert values["pc"] == 0.0
    assert values["kac"] == 0.0
    assert values["kbc"] == 0.0
    assert values["kcb"] == 700.0
    assert values["pa"] + values["pb"] == pytest.approx(1.0)


def test_partner_binding_zero_ligand_retains_parameter_defined_rate_splitting() -> None:
    values = _construct_and_resolve(
        "3st_binding_partner_2st",
        8.0e-4,
        0.0,
        {
            "koff_ab": 80.0,
            "kd_ab": 4.0e-4,
            "koff_ac": 130.0,
            "kd_ac": 1.3e-3,
            "keq": 2.5,
            "kex_bc": 700.0,
        },
    )

    assert values["l1_free"] == 0.0
    assert values["l2_free"] == 0.0
    assert values["pl1"] == 0.0
    assert values["pl2"] == 0.0
    assert (values["pa"], values["pb"], values["pc"]) == (1.0, 0.0, 0.0)
    assert values["kab"] == 0.0
    assert values["kac"] == 0.0
    assert values["kbc"] + values["kcb"] == pytest.approx(700.0)
    assert values["kbc"] / values["kcb"] == pytest.approx(2.5 * 4.0e-4 / 1.3e-3)


@pytest.mark.parametrize("edge_mask", range(16))
def test_four_state_partner_populations_ignore_all_edge_scale_patterns(
    edge_mask: int,
) -> None:
    p_total = 8.0e-4
    scales = (80.0, 130.0, 700.0, 900.0)
    koff_ab, koff_ac, kex_bc, kex_cd = (
        scale if edge_mask & (1 << index) else 0.0 for index, scale in enumerate(scales)
    )
    values = _construct_and_resolve(
        "4st_binding_partner_2st",
        p_total,
        2.3e-3,
        {
            "koff_ab": koff_ab,
            "kd_ab": 4.0e-4,
            "koff_ac": koff_ac,
            "kd_ac": 1.3e-3,
            "keq_l": 2.5,
            "keq_pl": 1.7,
            "kex_bc": kex_bc,
            "kex_cd": kex_cd,
        },
    )

    assert sum(values[name] for name in ("pa", "pb", "pc", "pd")) == (
        pytest.approx(1.0)
    )
    assert values["pa"] == pytest.approx(values["p_free"] / p_total)
    assert values["pb"] == pytest.approx(values["pl1"] / p_total)
    assert values["pc"] == pytest.approx(values["pl2"] / p_total)
    assert values["pd"] == pytest.approx(values["pl3"] / p_total)
    assert values["l2_free"] == pytest.approx(2.5 * values["l1_free"])
    assert values["pl3"] == pytest.approx(1.7 * values["pl2"])
    assert 4.0e-4 * values["pl1"] == pytest.approx(
        values["p_free"] * values["l1_free"], rel=2.0e-14
    )
    assert 1.3e-3 * values["pl2"] == pytest.approx(
        values["p_free"] * values["l2_free"], rel=2.0e-14
    )


def test_four_state_partner_strong_binding_is_non_negative_and_conservative() -> None:
    p_total = 1.0e-3
    l_total = 1.0e-3
    values = _construct_and_resolve(
        "4st_binding_partner_2st",
        p_total,
        l_total,
        {
            "koff_ab": 0.0,
            "kd_ab": 1.0e-300,
            "koff_ac": 0.0,
            "kd_ac": 1.0e-300,
            "keq_l": 2.5,
            "keq_pl": 1.7,
            "kex_bc": 0.0,
            "kex_cd": 0.0,
        },
    )

    concentrations = (
        values["p_free"],
        values["l1_free"],
        values["l2_free"],
        values["pl1"],
        values["pl2"],
        values["pl3"],
    )
    assert all(math.isfinite(value) and value >= 0.0 for value in concentrations)
    assert values["p_free"] + values["pl1"] + values["pl2"] + values["pl3"] == (
        pytest.approx(p_total, rel=2.0e-15, abs=0.0)
    )
    assert (
        values["l1_free"]
        + values["l2_free"]
        + values["pl1"]
        + values["pl2"]
        + values["pl3"]
    ) == pytest.approx(l_total, rel=2.0e-15, abs=0.0)


@pytest.mark.parametrize(("keq_l", "keq_pl"), ((0.0, 1.7), (2.5, 0.0)))
def test_four_state_partner_zero_equilibrium_ratio_removes_exact_state_weights(
    keq_l: float,
    keq_pl: float,
) -> None:
    values = _construct_and_resolve(
        "4st_binding_partner_2st",
        8.0e-4,
        2.3e-3,
        {
            "koff_ab": 80.0,
            "kd_ab": 4.0e-4,
            "koff_ac": 130.0,
            "kd_ac": 1.3e-3,
            "keq_l": keq_l,
            "keq_pl": keq_pl,
            "kex_bc": 700.0,
            "kex_cd": 900.0,
        },
    )

    if keq_l == 0.0:
        assert values["l2_free"] == 0.0
        assert values["pl2"] == 0.0
        assert values["pl3"] == 0.0
        assert values["kbc"] == 0.0
        assert values["kcb"] == 700.0
    if keq_pl == 0.0:
        assert values["pl3"] == 0.0
        assert values["kcd"] == 0.0
        assert values["kdc"] == 900.0


@pytest.mark.parametrize("edge_mask", range(8))
def test_three_bound_state_populations_ignore_all_kinetic_scale_patterns(
    edge_mask: int,
) -> None:
    p_total = 8.0e-4
    scales = (80.0, 700.0, 900.0)
    koff_ab, kex_bc, kex_cd = (
        scale if edge_mask & (1 << index) else 0.0 for index, scale in enumerate(scales)
    )
    values = _construct_and_resolve(
        "4st_binding_3_bound_states",
        p_total,
        2.3e-3,
        {
            "kd_app": 4.0e-4,
            "koff_ab": koff_ab,
            "kex_bc": kex_bc,
            "keq_bc": 2.5,
            "kex_cd": kex_cd,
            "keq_cd": 1.7,
        },
    )

    populations = tuple(values[name] for name in ("pa", "pb", "pc", "pd"))
    assert sum(populations) == pytest.approx(1.0)
    assert values["pa"] == pytest.approx(values["c_p"] / p_total)
    assert values["pb"] == pytest.approx(values["c_pl1"] / p_total)
    assert values["pc"] == pytest.approx(values["c_pl2"] / p_total)
    assert values["pd"] == pytest.approx(values["c_pl3"] / p_total)
    assert values["c_pl2"] == pytest.approx(2.5 * values["c_pl1"])
    assert values["c_pl3"] == pytest.approx(1.7 * values["c_pl2"])
    assert values["kd_ab"] == pytest.approx(4.0e-4 * (1.0 + 2.5 + 2.5 * 1.7))
    assert values["kd_ab"] * values["c_pl1"] == pytest.approx(
        values["c_p"] * values["c_l"], rel=2.0e-14
    )
    assert values["kd_eff"] == values["kd_app"]


def test_three_bound_state_strong_binding_is_non_negative_and_conservative() -> None:
    p_total = 1.0e-3
    l_total = 1.0e-3
    values = _construct_and_resolve(
        "4st_binding_3_bound_states",
        p_total,
        l_total,
        {
            "kd_app": 1.0e-300,
            "koff_ab": 0.0,
            "kex_bc": 0.0,
            "keq_bc": 2.5,
            "kex_cd": 0.0,
            "keq_cd": 1.7,
        },
    )

    concentrations = tuple(
        values[name] for name in ("c_p", "c_l", "c_pl1", "c_pl2", "c_pl3")
    )
    assert all(math.isfinite(value) and value >= 0.0 for value in concentrations)
    assert values["c_p"] + values["c_pl"] == pytest.approx(
        p_total, rel=2.0e-15, abs=0.0
    )
    assert values["c_l"] + values["c_pl"] == pytest.approx(
        l_total, rel=2.0e-15, abs=0.0
    )
    assert values["kd_eff"] == 1.0e-300


@pytest.mark.parametrize(("keq_bc", "keq_cd"), ((0.0, 1.7), (2.5, 0.0)))
def test_three_bound_state_zero_equilibrium_ratios_remove_exact_states(
    keq_bc: float,
    keq_cd: float,
) -> None:
    values = _construct_and_resolve(
        "4st_binding_3_bound_states",
        8.0e-4,
        2.3e-3,
        {
            "kd_app": 4.0e-4,
            "koff_ab": 80.0,
            "kex_bc": 700.0,
            "keq_bc": keq_bc,
            "kex_cd": 900.0,
            "keq_cd": keq_cd,
        },
    )

    if keq_bc == 0.0:
        assert values["c_pl2"] == 0.0
        assert values["c_pl3"] == 0.0
        assert values["kbc"] == 0.0
        assert values["kcb"] == 700.0
    if keq_cd == 0.0:
        assert values["c_pl3"] == 0.0
        assert values["kcd"] == 0.0
        assert values["kdc"] == 900.0


def test_three_bound_state_zero_kd_app_is_rejected_at_metadata_and_runtime() -> None:
    session, _ = _prepare_model(
        "4st_binding_3_bound_states",
        1.0e-3,
        2.0e-3,
        {"kd_app": 0.0},
    )
    assert not session.try_build_analysis_values()
    assert isinstance(
        session.parameter_factory.native_construction_error,
        InvalidConfigurationError,
    )

    overridden, _ = _prepare_model(
        "4st_binding_3_bound_states",
        1.0e-3,
        2.0e-3,
        {"kd_app": DefaultSetting(0.0, min=0.0, max=1.0)},
    )
    assert not overridden.try_build_analysis_values()
    error = overridden.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "KD must be finite and strictly positive" in str(error.__cause__)


CATEGORY_A_BINDING_MODELS = (
    "2st_binding",
    "3st_double_binding",
    "3st_binding_partner_2st",
    "4st_binding_partner_2st",
    "4st_binding_3_bound_states",
)


@pytest.mark.parametrize("model_name", CATEGORY_A_BINDING_MODELS)
def test_every_category_a_binding_model_rejects_zero_protein_total(
    model_name: str,
) -> None:
    defaults = _zero_rate_defaults(
        model_name,
        "kd_app" if model_name == "4st_binding_3_bound_states" else "kd",
        1.0e-6,
    )
    if model_name != "2st_binding" and "kd" in defaults:
        defaults.pop("kd")
    session, _ = _prepare_model(model_name, 0.0, 2.0e-3, defaults)

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "P_total must be finite and strictly positive" in str(error.__cause__)


@pytest.mark.parametrize(
    ("model_name", "defaults"),
    (
        ("2st_binding", {"kd": 1.0e-6, "koff": 100.0}),
        (
            "3st_double_binding",
            {
                "kd_ab": 1.0e-6,
                "kd_ac": 2.0e-6,
                "koff_ab": 100.0,
                "koff_ac": 120.0,
            },
        ),
        (
            "3st_binding_partner_2st",
            {
                "kd_ab": 1.0e-6,
                "kd_ac": 2.0e-6,
                "koff_ab": 100.0,
                "koff_ac": 120.0,
                "keq": 2.5,
                "kex_bc": 700.0,
            },
        ),
        (
            "4st_binding_partner_2st",
            {
                "kd_ab": 1.0e-6,
                "kd_ac": 2.0e-6,
                "koff_ab": 100.0,
                "koff_ac": 120.0,
                "keq_l": 2.5,
                "keq_pl": 1.7,
                "kex_bc": 700.0,
                "kex_cd": 900.0,
            },
        ),
        (
            "4st_binding_3_bound_states",
            {
                "kd_app": 1.0e-6,
                "koff_ab": 100.0,
                "keq_bc": 2.5,
                "keq_cd": 1.7,
                "kex_bc": 700.0,
                "kex_cd": 900.0,
            },
        ),
    ),
)
def test_every_category_a_binding_model_accepts_zero_ligand(
    model_name: str,
    defaults: dict[str, float],
) -> None:
    values = _construct_and_resolve(model_name, 1.0e-3, 0.0, defaults)
    populations = tuple(
        values[name] for name in ("pa", "pb", "pc", "pd") if name in values
    )

    assert populations[0] == 1.0
    assert all(population == 0.0 for population in populations[1:])


KD_CASES = (
    ("2st_binding", "kd"),
    ("3st_double_binding", "kd_ab"),
    ("3st_double_binding", "kd_ac"),
    ("3st_binding_partner_2st", "kd_ab"),
    ("3st_binding_partner_2st", "kd_ac"),
    ("4st_binding_partner_2st", "kd_ab"),
    ("4st_binding_partner_2st", "kd_ac"),
    ("4st_binding_3_bound_states", "kd_app"),
)


def _zero_rate_defaults(model_name: str, kd_name: str, kd: float) -> dict[str, float]:
    defaults: dict[str, float] = {kd_name: kd}
    if model_name == "2st_binding":
        return {**defaults, "koff": 0.0}
    if model_name == "3st_double_binding":
        return {
            "kd_ab": 1.0e-6,
            "kd_ac": 2.0e-6,
            "koff_ab": 0.0,
            "koff_ac": 0.0,
            **defaults,
        }
    if model_name == "3st_binding_partner_2st":
        return {
            "kd_ab": 1.0e-6,
            "kd_ac": 2.0e-6,
            "koff_ab": 0.0,
            "koff_ac": 0.0,
            "keq": 2.5,
            "kex_bc": 0.0,
            **defaults,
        }
    if model_name == "4st_binding_partner_2st":
        return {
            "kd_ab": 1.0e-6,
            "kd_ac": 2.0e-6,
            "koff_ab": 0.0,
            "koff_ac": 0.0,
            "keq_l": 2.5,
            "keq_pl": 1.7,
            "kex_bc": 0.0,
            "kex_cd": 0.0,
            **defaults,
        }
    return {
        "kd_app": kd,
        "koff_ab": 0.0,
        "keq_bc": 2.5,
        "keq_cd": 1.7,
        "kex_bc": 0.0,
        "kex_cd": 0.0,
    }


@pytest.mark.parametrize(("model_name", "kd_name"), KD_CASES)
@pytest.mark.parametrize(
    "kd",
    (1.0e-31, 1.0e-32, 1.0e-50, 1.0e-100, sys.float_info.min, math.nextafter(0.0, 1.0)),
)
def test_every_binding_kd_preserves_positive_values_below_legacy_floors(
    model_name: str,
    kd_name: str,
    kd: float,
) -> None:
    values = _construct_and_resolve(
        model_name,
        1.0e-3,
        2.0e-3,
        _zero_rate_defaults(model_name, kd_name, kd),
    )

    assert values[kd_name] == kd
    assert all(math.isfinite(value) and value >= 0.0 for value in values.values())
    populations = tuple(
        values[name] for name in ("pa", "pb", "pc", "pd") if name in values
    )
    assert sum(populations) == pytest.approx(1.0, rel=2.0e-14)


@pytest.mark.parametrize(("model_name", "kd_name"), KD_CASES)
def test_every_binding_kd_rejects_zero_even_with_a_bounds_override(
    model_name: str,
    kd_name: str,
) -> None:
    defaults = _zero_rate_defaults(model_name, kd_name, 1.0e-6)
    defaults[kd_name] = DefaultSetting(0.0, min=0.0, max=1.0)  # ty: ignore[invalid-assignment]
    session, _ = _prepare_model(model_name, 1.0e-3, 2.0e-3, defaults)

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "KD must be finite and strictly positive" in str(error.__cause__)
