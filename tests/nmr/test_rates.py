from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from chemex import _entrypoint
from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.exceptions import ChemExError
from chemex.experiments.builder import build_experiments
from chemex.models.model import ModelSpec
from chemex.nmr import rates as nmr_rates
from chemex.nmr.basis import Basis
from chemex.nmr.rates import get_model_free_expressions, rate_functions
from chemex.parameters.spins import create_base_param_settings
from chemex.runtime import AnalysisSession
from chemex.terminal import format_known_failure

_EXAMPLES_ROOT = Path(__file__).resolve().parents[2] / "examples" / "Experiments"

_SWAPPED_COMPONENTS = {
    "r2_i": "r2_s",
    "r2_s": "r2_i",
    "r1_i": "r1_s",
    "r1_s": "r1_i",
    "r2a_i": "r2a_s",
    "r2a_s": "r2a_i",
    "r2mq_is": "r2mq_is",
    "r1a_is": "r1a_is",
    "etaxy_i": "etaxy_s",
    "etaxy_s": "etaxy_i",
    "etaz_i": "etaz_s",
    "etaz_s": "etaz_i",
    "sigma_is": "sigma_is",
    "mu_is": "mu_is",
}


def _model_free_expressions(
    spin_system: str,
    *,
    deuterated: bool = False,
) -> dict[str, str]:
    basis = Basis(
        type="ixyzsxyz",
        spin_system=spin_system,
        model=ModelSpec(),
    )
    label = ("2h",) if deuterated else ()
    return get_model_free_expressions(
        basis,
        Conditions(h_larmor_frq=799.708, label=label),
    )


@pytest.mark.parametrize(
    ("canonical_orientation", "reversed_orientation"),
    (("nh", "hn"), ("ch", "hc")),
)
@pytest.mark.parametrize("state", ("a", "b"))
def test_reversed_orientation_expressions_use_one_physical_identity(
    canonical_orientation: str,
    reversed_orientation: str,
    state: str,
) -> None:
    canonical = _model_free_expressions(canonical_orientation)
    reversed_ = _model_free_expressions(reversed_orientation)

    for component, reversed_component in _SWAPPED_COMPONENTS.items():
        assert (
            canonical[f"{component}_{state}"]
            == reversed_[f"{reversed_component}_{state}"]
        )


def test_nh_hn_expressions_canonicalize_physical_n_and_h_components() -> None:
    nh = _model_free_expressions("nh")
    hn = _model_free_expressions("hn")

    assert (
        nh["r1_i_a"]
        == hn["r1_s_a"]
        == ("nh(799.708, {tauc_a}, {s2_a}, {khh_a})['r1_i']")
    )
    assert (
        nh["r1_s_a"]
        == hn["r1_i_a"]
        == ("nh(799.708, {tauc_a}, {s2_a}, {khh_a})['r1_s']")
    )


def test_diffusion_coefficient_is_non_negative() -> None:
    basis = Basis(type="ixyzsxyz_diff", spin_system="ch", model=ModelSpec())

    settings = create_base_param_settings(basis, "a", Conditions())

    assert settings["d_a"].min == 0.0


def test_relaxation_psd_specs_match_the_actual_basis_transition_support() -> None:
    basis = Basis(
        type="ixyzsxyz_eq",
        spin_system="nh",
        model=ModelSpec(),
    )
    non_equilibrium = tuple(
        index
        for index, component in enumerate(basis.components_states)
        if not component.startswith(("ie_", "se_"))
    )

    for specification in basis.relaxation_psd_specs:
        diagonal_supports = []
        for name in specification.diagonal_names:
            matrix = basis.matrices[name]
            diagonal_supports.append(
                {index for index in non_equilibrium if matrix[index, index] < 0.0}
            )
        for row, column, name in specification.off_diagonal_names:
            matrix = basis.matrices[name]
            observed = {
                (left, right)
                for left in non_equilibrium
                for right in non_equilibrium
                if left != right and matrix[left, right] != 0.0
            }
            expected_support = {
                (left, right)
                for left in diagonal_supports[row]
                for right in diagonal_supports[column]
            } | {
                (right, left)
                for left in diagonal_supports[row]
                for right in diagonal_supports[column]
            }
            assert observed
            assert observed <= expected_support


@pytest.mark.parametrize("component", ("r2mq_is", "r1a_is", "sigma_is", "mu_is"))
def test_reversed_orientation_preserves_symmetric_is_component(
    component: str,
) -> None:
    nh = _model_free_expressions("nh")
    hn = _model_free_expressions("hn")

    assert nh[f"{component}_a"] == hn[f"{component}_a"]


@pytest.mark.parametrize(
    ("canonical_orientation", "reversed_orientation"),
    (("nh", "hn"), ("ch", "hc")),
)
def test_deuterated_reversed_orientation_expressions_are_canonical(
    canonical_orientation: str,
    reversed_orientation: str,
) -> None:
    canonical = _model_free_expressions(canonical_orientation, deuterated=True)
    reversed_ = _model_free_expressions(reversed_orientation, deuterated=True)

    for component, reversed_component in _SWAPPED_COMPONENTS.items():
        assert canonical[f"{component}_a"] == reversed_[f"{reversed_component}_a"]
        assert f"{canonical_orientation}_d(" in canonical[f"{component}_a"]


@pytest.mark.parametrize(
    ("canonical_orientation", "reversed_orientation", "arguments"),
    (
        ("nh", "hn", (799.708, 8.7, 0.83, 2.4)),
        ("nh_d", "hn_d", (799.708, 8.7, 0.83, 2.4)),
        ("ch", "hc", (799.708, 8.7, 0.83)),
        ("ch_d", "hc_d", (799.708, 8.7, 0.83)),
    ),
)
def test_reversed_rate_functions_are_numerically_equivalent(
    canonical_orientation: str,
    reversed_orientation: str,
    arguments: tuple[float, ...],
) -> None:
    canonical = rate_functions[canonical_orientation](*arguments)
    reversed_ = rate_functions[reversed_orientation](*arguments)

    for component, reversed_component in _SWAPPED_COMPONENTS.items():
        np.testing.assert_allclose(
            canonical[component],
            reversed_[reversed_component],
            rtol=1e-13,
            atol=0.0,
        )


def test_khh_remains_attached_to_the_physical_proton() -> None:
    nh_without_exchange = rate_functions["nh"](799.708, 8.7, 0.83)
    nh_with_exchange = rate_functions["nh"](799.708, 8.7, 0.83, 2.4)
    hn_without_exchange = rate_functions["hn"](799.708, 8.7, 0.83)
    hn_with_exchange = rate_functions["hn"](799.708, 8.7, 0.83, 2.4)

    assert nh_with_exchange["r1_s"] - nh_without_exchange["r1_s"] == pytest.approx(
        2.4,
    )
    assert hn_with_exchange["r1_i"] - hn_without_exchange["r1_i"] == pytest.approx(
        2.4,
    )
    assert nh_with_exchange["r1_i"] == nh_without_exchange["r1_i"]
    assert hn_with_exchange["r1_s"] == hn_without_exchange["r1_s"]


@pytest.mark.parametrize(("tauc", "s2"), ((0.0, 0.9), (4.0, 0.0)))
def test_model_free_exact_zero_limits_produce_zero_rates(
    tauc: float,
    s2: float,
) -> None:
    rates = rate_functions["nh"](799.708, tauc, s2, 0.0)
    assert all(rate == 0.0 for rate in rates.values())


def test_cn_model_free_aliases_are_absent_and_explicitly_rejected() -> None:
    assert "cn" not in rate_functions
    assert "cn_d" not in rate_functions
    assert not nmr_rates.supports_model_free("cn")
    for supported in ("nh", "hn", "ch", "hc"):
        assert nmr_rates.supports_model_free(supported)

    basis = Basis(
        type="ixyzsz",
        spin_system="cn",
        model=ModelSpec(name="2st.mf", states="ab", model_free=True),
    )
    with pytest.raises(nmr_rates.UnsupportedModelFreeError) as exc_info:
        get_model_free_expressions(basis, Conditions(h_larmor_frq=800.0))

    error = exc_info.value
    assert isinstance(error, ChemExError)
    assert error.spin_system == "cn"
    assert error.model_name == "2st.mf"
    formatted = format_known_failure(error)
    assert "Model-free (.mf) relaxation is not supported" in formatted
    assert "Model: 2st.mf" in formatted
    assert "Spin System: cn" in formatted


def test_tauc_docstrings_specify_nanoseconds() -> None:
    for obj in (
        nmr_rates._calculate_jw,
        nmr_rates.RatesIS.__call__,
        nmr_rates.RateNH.__call__,
    ):
        assert obj.__doc__ is not None
        assert "nanoseconds" in obj.__doc__
        assert "Correlation time in seconds" not in obj.__doc__


def test_cli_entrypoint_reports_known_failure_for_cn_model_free(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    exp_dir = _EXAMPLES_ROOT / "CPMG_13CO_AP"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "chemex",
            "simulate",
            "-e",
            str(exp_dir / "Experiments" / "800MHz.toml"),
            "-p",
            str(exp_dir / "Parameters" / "parameters.toml"),
            "-d",
            "2st.mf",
            "-o",
            str(tmp_path / "out_mf"),
            "--plot",
            "nothing",
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        _entrypoint.main()
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    combined = f"{captured.out}\n{captured.err}"
    assert "Model-free (.mf) relaxation is not supported" in combined
    assert "unexpected internal error" not in combined.lower()


def test_cpmg_13co_ap_rejects_mf_and_honors_r1a_in_non_mf() -> None:
    exp_dir = _EXAMPLES_ROOT / "CPMG_13CO_AP"
    exp_files = [
        exp_dir / "Experiments" / "500MHz.toml",
        exp_dir / "Experiments" / "800MHz.toml",
    ]
    par_files = [exp_dir / "Parameters" / "parameters.toml"]

    # 1. Explicit rejection of 2st.mf
    session_mf = AnalysisSession.create()
    session_mf.set_model("2st.mf")
    with pytest.raises(nmr_rates.UnsupportedModelFreeError):
        build_experiments(
            exp_files,
            Selection(None, None),
            session=session_mf,
        )

    # 2. Non-.mf 2st succeeds, honors R1A_A and R1_A, and never seeds RateCH
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        exp_files,
        Selection(None, None),
        session=session,
    )
    assert session.parameter_factory.sealed_definitions is not None
    session.parameters.set_defaults(read_defaults(par_files))
    assert session.try_build_analysis_values()

    profile_800 = next(
        profile
        for exp in experiments
        if exp.filename.name == "800MHz.toml"
        for profile in exp
    )
    required_ids = set(profile_800.name_map.values())
    parameterization = session.compile_parameterization(Method(), required_ids)
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )

    r1_id = profile_800.name_map["r1_i_a"]
    r1a_id = profile_800.name_map["r1a_is_a"]
    r2_id = profile_800.name_map["r2_i_a"]
    r2a_id = profile_800.name_map["r2a_i_a"]

    assert resolved[r1_id] == pytest.approx(1.00)
    assert resolved[r1a_id] == pytest.approx(3.30)
    assert resolved[r2a_id] == pytest.approx(resolved[r2_id] + 3.30 - 1.00)


def test_non_mf_cn_without_defaults_uses_spin_defaults_in_mixed_session() -> None:
    co_file = _EXAMPLES_ROOT / "CPMG_13CO_AP" / "Experiments" / "800MHz.toml"
    nh_file = _EXAMPLES_ROOT / "CPMG_15N_IP_0013" / "Experiments" / "600MHz.toml"
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [co_file, nh_file],
        Selection(None, None),
        session=session,
    )
    assert session.parameter_factory.sealed_definitions is not None
    session.parameters.set_defaults([])
    assert session.try_build_analysis_values()

    exp_list = list(experiments)
    co_profile = exp_list[0].profiles[0]
    nh_profile = exp_list[1].profiles[0]
    required_ids = set(co_profile.name_map.values()) | set(nh_profile.name_map.values())
    parameterization = session.compile_parameterization(Method(), required_ids)
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )

    # Carbonyl profile keeps base spin defaults rather than aliphatic RateCH(800, 4.0, 0.9)
    assert resolved[co_profile.name_map["r1_i_a"]] == pytest.approx(1.5)
    assert resolved[co_profile.name_map["r1a_is_a"]] == pytest.approx(4.0)
    assert resolved[co_profile.name_map["r2_i_a"]] == pytest.approx(10.0)
    assert resolved[co_profile.name_map["r2a_i_a"]] == pytest.approx(12.5)

    # Co-loaded 15N profile still receives its model-free bootstrap defaults
    expected_nh = rate_functions["nh_d"](600.25, 4.0, 0.9, 0.0)
    assert resolved[nh_profile.name_map["r1_i_a"]] == pytest.approx(
        float(expected_nh["r1_i"])
    )


def test_cpmg_ch3_13c_h2c_honors_r1a_a_in_non_mf() -> None:
    exp_dir = _EXAMPLES_ROOT / "CPMG_CH3_13C_H2C"
    exp_files = [exp_dir / "Experiments" / "600MHz.toml"]
    par_files = [exp_dir / "Parameters" / "parameters.toml"]

    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        exp_files,
        Selection(None, None),
        session=session,
    )
    assert session.parameter_factory.sealed_definitions is not None
    session.parameters.set_defaults(read_defaults(par_files))
    assert session.try_build_analysis_values()

    profile = next(iter(experiments)).profiles[0]
    parameterization = session.compile_parameterization(
        Method(),
        set(profile.name_map.values()),
    )
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )

    assert resolved[profile.name_map["r1_i_a"]] == pytest.approx(1.5)
    assert resolved[profile.name_map["r1a_is_a"]] == pytest.approx(3.0)
    assert resolved[profile.name_map["r2a_i_a"]] == pytest.approx(
        resolved[profile.name_map["r2_i_a"]] + 3.0 - 1.5
    )


@pytest.mark.parametrize(
    ("canonical", "reversed_"),
    (("nh", "hn"), ("ch", "hc")),
)
def test_non_mf_antiphase_expressions_preserve_is_exchange_symmetry(
    canonical: str,
    reversed_: str,
) -> None:
    cond = Conditions(h_larmor_frq=800.0)
    can_settings = create_base_param_settings(
        Basis(type="ixyzsxyz", spin_system=canonical, model=ModelSpec()),
        "a",
        cond,
    )
    rev_settings = create_base_param_settings(
        Basis(type="ixyzsxyz", spin_system=reversed_, model=ModelSpec()),
        "a",
        cond,
    )

    assert can_settings["r1a_is_a"].expr == ""
    assert can_settings["r2a_i_a"].expr == "{r2_i_a} + {r1a_is_a} - {r1_i_a}"
    assert can_settings["r2a_s_a"].expr == "{r2_s_a} - {r1_i_a}"

    assert rev_settings["r1a_is_a"].expr == "{r1_i_a} + {r1_s_a}"
    assert rev_settings["r2a_i_a"].expr == "{r2_i_a} - {r1_s_a}"
    assert rev_settings["r2a_s_a"].expr == "{r2_s_a} + {r1a_is_a} - {r1_s_a}"


@pytest.mark.parametrize("tauc", (1.0, 4.0, 10.0, 20.0))
@pytest.mark.parametrize("deuterated", (False, True))
def test_non_mf_antiphase_r2a_expressions_track_redfield_physics(
    tauc: float,
    deuterated: bool,
) -> None:
    suffix = "_d" if deuterated else ""
    for spin_system in ("nh", "ch"):
        rf = rate_functions[f"{spin_system}{suffix}"](800.0, tauc, 0.9)
        approx_r2a_i = float(rf["r2_i"] + rf["r1a_is"] - rf["r1_i"])
        assert rf["r2a_i"] > rf["r2_i"]
        assert approx_r2a_i == pytest.approx(float(rf["r2a_i"]), rel=0.02, abs=0.05)

    for spin_system in ("hn", "hc"):
        rf = rate_functions[f"{spin_system}{suffix}"](800.0, tauc, 0.9)
        approx_r2a_i = float(rf["r2_i"] - rf["r1_s"])
        approx_r1a_is = float(rf["r1_i"] + rf["r1_s"])
        approx_r2a_s = float(rf["r2_s"] + approx_r1a_is - rf["r1_s"])
        old_r2a_s = float(rf["r2_s"] - rf["r1_i"])
        assert rf["r2a_s"] > rf["r2_s"]
        assert abs(approx_r2a_s - float(rf["r2a_s"])) < abs(
            old_r2a_s - float(rf["r2a_s"])
        )
        if spin_system == "hc":
            assert rf["r2a_i"] < rf["r2_i"]
            assert abs(approx_r2a_i - float(rf["r2a_i"])) < abs(
                float(rf["r2_i"] + rf["r1_s"]) - float(rf["r2a_i"])
            )
            if tauc >= 4.0:
                assert approx_r2a_i == pytest.approx(
                    float(rf["r2a_i"]),
                    rel=0.02,
                    abs=0.05,
                )
                assert approx_r2a_s == pytest.approx(
                    float(rf["r2a_s"]),
                    rel=0.02,
                    abs=0.05,
                )
        elif tauc >= 4.0:
            assert rf["r2a_i"] < rf["r2_i"]
            assert abs(approx_r2a_i - float(rf["r2a_i"])) < abs(
                float(rf["r2_i"] + rf["r1_s"]) - float(rf["r2a_i"])
            )
