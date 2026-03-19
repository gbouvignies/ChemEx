from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.detection import build_detection_vector
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.parameters.spin_system import SpinSystem


def test_build_detection_vector_supports_addition_and_subtraction() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())

    vector = build_detection_vector("[2izsz_a] - [iz_a] + [iz_b]", basis.vectors)
    expected = (
        basis.vectors["2izsz_a"] - basis.vectors["iz_a"] + basis.vectors["iz_b"]
    )

    np.testing.assert_allclose(vector, expected)


def test_build_detection_vector_rejects_invalid_syntax() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())

    with pytest.raises(ValueError, match="Invalid detection expression"):
        build_detection_vector("[iz_a] [2izsz_a]", basis.vectors)


def test_build_detection_vector_rejects_unknown_components() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())

    with pytest.raises(ValueError, match="Unknown detection component"):
        build_detection_vector("[unknown_component]", basis.vectors)


def test_liouvillian_detect_uses_parsed_detection_expression() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())
    liouvillian = LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.detection = "[2izsz_a] - [iz_a]"

    magnetization = 3.0 * basis.vectors["2izsz_a"] + 1.5 * basis.vectors["iz_a"]

    assert liouvillian.detect(magnetization) == pytest.approx(1.5)


def test_liouvillian_detection_state_does_not_change_on_parse_failure() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())
    liouvillian = LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.detection = "[iz_a]"

    with pytest.raises(ValueError, match="Unknown detection component"):
        liouvillian.detection = "[unknown_component]"

    assert liouvillian.detection == "[iz_a]"
    assert liouvillian.detect(basis.vectors["iz_a"]) == pytest.approx(1.0)
