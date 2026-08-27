from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.rates import get_model_free_expressions, rate_functions

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
