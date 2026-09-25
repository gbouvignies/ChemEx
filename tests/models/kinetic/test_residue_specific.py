from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession


def setup_module() -> None:
    register_kinetic_settings()


def test_2st_rs_modifier_preserves_settings_except_kinetic_scope() -> None:
    conditions = Conditions()

    base_settings = model_factory.create("2st", conditions)
    settings = model_factory.create_for_model(ModelSpec.from_name("2st.rs"), conditions)

    assert settings.keys() == base_settings.keys()

    for key, base in base_settings.items():
        modified = settings[key]
        assert modified.name_setting == replace(
            base.name_setting,
            spin_system_part="g",
        )
        assert modified.value == base.value
        assert modified.min == base.min
        assert modified.max == base.max
        assert modified.vary == base.vary
        assert modified.expr == base.expr
        assert modified.supports_estimation == base.supports_estimation
        assert modified.report_only == base.report_only
        assert modified.model_owned == base.model_owned


def test_three_state_models_become_residue_specific_with_rs_suffix() -> None:
    base_settings = model_factory.create("3st", Conditions())
    residue_specific_settings = model_factory.create_for_model(
        ModelSpec.from_name("3st.rs"),
        Conditions(),
    )

    assert base_settings["pb"].name_setting.spin_system_part == ""

    for key in ("pb", "pc", "kex_ab", "kex_ac", "kab", "kba", "kac", "kca", "pa"):
        assert residue_specific_settings[key].name_setting.spin_system_part == "g"


def test_hd_models_keep_d2o_global_with_rs_suffix() -> None:
    settings = model_factory.create_for_model(
        ModelSpec.from_name("2st_hd.rs"),
        Conditions(d2o=0.2),
    )

    assert settings["d2o"].name_setting.spin_system_part == ""

    for key in ("kdh", "phi", "kab", "kba", "pa", "pb"):
        assert settings[key].name_setting.spin_system_part == "g"


@pytest.mark.parametrize(
    "name",
    ("3st_binding_cs", "3st_monomer_dimer_trimer", "3st_eyring_linear"),
)
def test_rs_changes_only_eligible_kinetic_scope_across_families(name: str) -> None:
    conditions = Conditions(temperature=25.0, p_total=1e-3, l_total=2e-3)
    base = model_factory.create(name, conditions)
    residue_specific = model_factory.create_for_model(
        ModelSpec.from_name(f"{name}.rs"),
        conditions,
    )

    assert base.keys() == residue_specific.keys()
    for key, setting in base.items():
        modified = residue_specific[key]
        assert modified.name_setting.name == setting.name_setting.name
        assert (
            modified.name_setting.conditions_part
            == setting.name_setting.conditions_part
        )
        assert modified.name_setting.spin_system_part == "g"
        assert modified.value == setting.value
        assert modified.min == setting.min
        assert modified.max == setting.max
        assert modified.vary == setting.vary
        assert modified.expr == setting.expr
        assert modified.report_only == setting.report_only


def test_parameter_factory_uses_residue_specific_model_settings() -> None:
    config = SimpleNamespace(
        conditions=Conditions(),
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    spin_system = SpinSystem.from_name("G23N-HN")

    session_base = AnalysisSession.create()
    session_base.set_model("2st")
    base_basis = Basis(type="iz", spin_system="nh", model=session_base.model.spec)
    base_ids = session_base.parameter_factory.create_parameters(
        config,
        basis=base_basis,
        spin_system=spin_system,
    )
    base_pb = session_base.parameters.get_parameters([base_ids["pb"]])[base_ids["pb"]]

    session_rs = AnalysisSession.create()
    session_rs.set_model("2st.rs")
    rs_basis = Basis(type="iz", spin_system="nh", model=session_rs.model.spec)
    rs_ids = session_rs.parameter_factory.create_parameters(
        config,
        basis=rs_basis,
        spin_system=spin_system,
    )
    rs_pb = session_rs.parameters.get_parameters([rs_ids["pb"]])[rs_ids["pb"]]

    assert base_pb.param_name.spin_system.name == ""
    assert rs_pb.param_name.spin_system.name == "G23"
