from __future__ import annotations

from types import SimpleNamespace

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession


def setup_module() -> None:
    register_kinetic_settings()


def test_2st_residue_specific_matches_legacy_alias() -> None:
    conditions = Conditions()

    settings = model_factory.create_for_model(ModelSpec.from_name("2st.rs"), conditions)
    legacy_settings = model_factory.create("2st_rs", conditions)

    assert settings.keys() == legacy_settings.keys()

    for key in settings:
        assert settings[key].name_setting == legacy_settings[key].name_setting
        assert settings[key].value == legacy_settings[key].value
        assert settings[key].min == legacy_settings[key].min
        assert settings[key].max == legacy_settings[key].max
        assert settings[key].vary == legacy_settings[key].vary
        assert settings[key].expr == legacy_settings[key].expr


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
