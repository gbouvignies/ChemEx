from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.parameters import read_defaults
from chemex.models.factory import model_factory
from chemex.models.loader import register_kinetic_settings
from chemex.nmr.basis import Basis
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spins import create_base_param_settings
from chemex.runtime import AnalysisSession

FULL_CONDITIONS = Conditions(
    h_larmor_frq=800.0,
    temperature=25.0,
    p_total=1.0e-4,
    l_total=2.0e-4,
    d2o=0.1,
    label=("2h",),
)


def setup_module() -> None:
    register_kinetic_settings()


def _assert_finite_ordered_bounds(
    name: str,
    value: float | None,
    lower: float,
    upper: float,
) -> None:
    assert math.isfinite(lower), f"{name} has a nonfinite lower bound"
    assert math.isfinite(upper), f"{name} has a nonfinite upper bound"
    assert lower < upper, f"{name} has an empty default domain"
    assert value is not None, f"{name} has no independent default value"
    assert lower <= value <= upper, f"{name} defaults outside its domain"


def test_every_registered_kinetic_independent_has_a_finite_default_domain() -> None:
    for model_name in sorted(model_factory.set):
        settings = model_factory.create(model_name, FULL_CONDITIONS)
        for local_name, setting in settings.items():
            if setting.expr:
                continue
            _assert_finite_ordered_bounds(
                f"{model_name}:{local_name}",
                setting.value,
                setting.min,
                setting.max,
            )


def test_approved_kinetic_family_domains_apply_only_to_independent_settings() -> None:
    expected_domains = {
        "kex": (0.0, 1.0e6),
        "koff": (0.0, 1.0e6),
        "kdh": (0.0, 1.0e6),
        "kd": (0.0, 1.0),
        "keq": (0.0, 1.0e6),
    }
    direct_rates = {"kab", "kba", "kbc", "kcb"}
    approved_equilibria = {"keq", "keq_l", "keq_pl"}

    for model_name in sorted(model_factory.set):
        for local_name, setting in model_factory.create(
            model_name, FULL_CONDITIONS
        ).items():
            if setting.expr:
                continue
            if local_name in direct_rates:
                expected = (0.0, 1.0e6)
            elif local_name in approved_equilibria:
                expected = expected_domains["keq"]
            else:
                family = next(
                    (
                        prefix
                        for prefix in ("kex", "koff", "kdh", "kd")
                        if local_name.startswith(prefix)
                    ),
                    None,
                )
                if family is None:
                    continue
                expected = expected_domains[family]
            assert (setting.min, setting.max) == expected, (
                model_name,
                local_name,
            )

    derived_kab = model_factory.create("2st", FULL_CONDITIONS)["kab"]
    assert derived_kab.expr
    assert math.isinf(derived_kab.max)


@pytest.mark.parametrize("extension", ("", "dq", "tq"))
def test_every_standard_spin_setting_has_a_finite_default_domain(
    extension: str,
) -> None:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    basis = Basis(
        type="ixyzsxyz_eq",
        spin_system="nh",
        extension=extension,  # ty: ignore[invalid-argument-type]
        model=session.model.spec,
    )

    for state in session.model.states:
        for local_name, setting in create_base_param_settings(
            basis, state, FULL_CONDITIONS
        ).items():
            _assert_finite_ordered_bounds(
                f"{extension or 'sq'}:{local_name}",
                setting.value,
                setting.min,
                setting.max,
            )


@pytest.mark.parametrize("extension", ("", "dq", "tq"))
def test_approved_standard_spin_family_domains(extension: str) -> None:
    session = AnalysisSession.create()
    session.set_model("2st")
    basis = Basis(
        type="ixyzsxyz_eq",
        spin_system="nh",
        extension=extension,  # ty: ignore[invalid-argument-type]
        model=session.model.spec,
    )
    settings = create_base_param_settings(basis, "a", FULL_CONDITIONS)
    expected = {
        "tauc_a": (0.0, 1000.0),
        "s2_a": (0.0, 1.0),
        "khh_a": (0.0, 1.0e6),
        "r2_i_a": (0.0, 1000.0),
        "r2_s_a": (0.0, 1000.0),
        "r1_i_a": (0.0, 1000.0),
        "r1_s_a": (0.0, 1000.0),
        "r2a_i_a": (0.0, 1000.0),
        "r2a_s_a": (0.0, 1000.0),
        "r2mq_is_a": (0.0, 1000.0),
        "r1a_is_a": (0.0, 1000.0),
        "etaxy_i_a": (-1000.0, 1000.0),
        "etaxy_s_a": (-1000.0, 1000.0),
        "etaz_i_a": (-1000.0, 1000.0),
        "etaz_s_a": (-1000.0, 1000.0),
        "sigma_is_a": (-1000.0, 1000.0),
        "mu_is_a": (-1000.0, 1000.0),
        "cs_i_a": (-100.0, 300.0),
        "cs_s_a": (-100.0, 300.0),
        "j_is_a": (-1000.0, 1000.0),
        "d_a": (0.0, 1000.0),
    }

    assert {
        local_name: (settings[local_name].min, settings[local_name].max)
        for local_name in expected
    } == expected


@pytest.mark.parametrize(
    ("model_name", "extension", "expected"),
    (
        ("2st", "", {"KEX_AB": (0.0, 1.0e6), "R2_A": (0.0, 1000.0)}),
        ("2st.rs", "", {"KEX_AB": (0.0, 1.0e6)}),
        (
            "2st.mf",
            "",
            {
                "TAUC_A": (0.0, 1000.0),
                "S2_A": (0.0, 1.0),
                "KHH_A": (0.0, 1.0e6),
            },
        ),
        (
            "2st.tc",
            "",
            {"CS_A": (-100.0, 300.0), "DWM_AB": (-1.0, 1.0)},
        ),
        (
            "2st",
            "dq",
            {
                "R1DQ_A": (0.0, 1000.0),
                "R2DQ_A": (0.0, 1000.0),
                "R1ADQ_A": (0.0, 1000.0),
                "R2ADQ_A": (0.0, 1000.0),
                "ETAXYDQ_A": (-1000.0, 1000.0),
                "ETAZDQ_A": (-1000.0, 1000.0),
            },
        ),
        (
            "2st",
            "tq",
            {
                "R1TQ_A": (0.0, 1000.0),
                "R2TQ_A": (0.0, 1000.0),
                "R1ATQ_A": (0.0, 1000.0),
                "R2ATQ_A": (0.0, 1000.0),
                "ETAXYTQ_A": (-1000.0, 1000.0),
                "ETAZTQ_A": (-1000.0, 1000.0),
            },
        ),
        (
            "3st_binding_cs",
            "",
            {
                "KAB": (0.0, 1.0e6),
                "KBA": (0.0, 1.0e6),
                "KOFF_BC": (0.0, 1.0e6),
                "KD_APP": (0.0, 1.0),
            },
        ),
    ),
)
def test_representative_parameterizations_preserve_family_bounds(
    model_name: str,
    extension: str,
    expected: dict[str, tuple[float, float]],
) -> None:
    session = AnalysisSession.create()
    session.set_model(model_name)
    config = SimpleNamespace(
        conditions=FULL_CONDITIONS,
        to_be_fitted=SimpleNamespace(
            rates=[
                "r1",
                "r2",
                "r1a",
                "r2a",
                "r2mq",
                "etaxy",
                "etaz",
                "sigma",
                "mu",
                "d",
            ],
            model_free=["tauc", "s2", "khh"],
        ),
    )
    basis = Basis(
        type="ixyzsxyz_eq",
        spin_system="nh",
        extension=extension,  # ty: ignore[invalid-argument-type]
        model=session.model.spec,
    )
    session.parameter_factory.create_parameters(
        config,
        basis=basis,
        spin_system=SpinSystem.from_name("G23N-HN"),
    )
    session.parameters.set_defaults([])
    session.parameter_factory.seal_definitions()
    session.parameter_factory.seal_configuration()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None

    for name, bounds in expected.items():
        matches = [definition for definition in definitions if definition.name == name]
        assert matches, name
        assert {
            (definition.lower_bound, definition.upper_bound) for definition in matches
        } == {bounds}


@pytest.mark.parametrize(
    ("parameter_value", "expected"),
    (
        ("25.0", (25.0, 0.0, 1000.0)),
        ("[25.0, 2.0, 50.0]", (25.0, 2.0, 50.0)),
    ),
)
def test_user_defaults_inherit_or_explicitly_override_standard_bounds(
    tmp_path: Path,
    parameter_value: str,
    expected: tuple[float, float, float],
) -> None:
    session = AnalysisSession.create()
    session.set_model("2st")
    config = SimpleNamespace(
        conditions=FULL_CONDITIONS,
        to_be_fitted=SimpleNamespace(rates=["r2"], model_free=[]),
    )
    basis = Basis(
        type="ixy",
        spin_system="nh",
        model=session.model.spec,
    )
    ids = session.parameter_factory.create_parameters(
        config,
        basis=basis,
        spin_system=SpinSystem.from_name("G23N-HN"),
    )
    param_id = ids["r2_i_a"]
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        f"[R2_A]\nG23N = {parameter_value}\n",
        encoding="utf-8",
    )
    session.parameters.set_defaults(read_defaults([parameters]))
    session.parameter_factory.seal_definitions()
    session.parameter_factory.seal_configuration()
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    resolved = configuration[param_id]

    assert (
        resolved.effective_value,
        resolved.lower_bound,
        resolved.upper_bound,
    ) == expected
