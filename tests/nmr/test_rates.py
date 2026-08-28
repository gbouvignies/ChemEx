from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.rates import get_model_free_expressions, rate_functions
from chemex.parameters.spins import create_base_param_settings

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
