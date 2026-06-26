from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr._engine.detection import (
    build_detection_vector,
    build_state_detection_expression,
)
from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.nmr.basis import Basis
from chemex.parameters.spin_system import SpinSystem


def test_build_detection_vector_supports_addition_and_subtraction() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())

    vector = build_detection_vector("[2izsz_a] - [iz_a] + [iz_b]", basis.vectors)
    expected = (
        basis.vectors["2izsz_a"] - basis.vectors["iz_a"] + basis.vectors["iz_b"]
    )

    np.testing.assert_allclose(vector, expected)


def test_build_state_detection_expression_preserves_component_signs() -> None:
    expression = build_state_detection_expression(
        "[2izsz] - [iz]",
        ("a", "c"),
    )

    assert expression == "[2izsz_a] - [iz_a] + [2izsz_c] - [iz_c]"


def test_multi_state_detection_sums_selected_final_components() -> None:
    model = ModelSpec(name="3st", states="abc")
    basis = Basis(type="ixyz", spin_system="nh", model=model)
    expression = build_state_detection_expression("[iz]", ("a", "c"))
    detection_vector = build_detection_vector(expression, basis.vectors).transpose()
    magnetization = (
        2.0 * basis.vectors["iz_a"]
        + 5.0 * basis.vectors["iz_b"]
        + 7.0 * basis.vectors["iz_c"]
    )

    detected = (detection_vector @ magnetization).item()

    assert detected == pytest.approx(9.0)


def test_explicit_all_state_detection_matches_legacy_unsuffixed_component() -> None:
    model = ModelSpec(name="3st", states="abc")
    basis = Basis(type="ixyz", spin_system="nh", model=model)
    expression = build_state_detection_expression("[iz]", model.states)

    explicit = build_detection_vector(expression, basis.vectors)
    legacy = build_detection_vector("[iz]", basis.vectors)

    np.testing.assert_allclose(explicit, legacy)


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
    liouvillian = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.detection = "[2izsz_a] - [iz_a]"

    magnetization = 3.0 * basis.vectors["2izsz_a"] + 1.5 * basis.vectors["iz_a"]

    assert liouvillian.detect(magnetization) == pytest.approx(1.5)


def test_liouvillian_detection_state_does_not_change_on_parse_failure() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())
    liouvillian = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.detection = "[iz_a]"

    with pytest.raises(ValueError, match="Unknown detection component"):
        liouvillian.detection = "[unknown_component]"

    assert liouvillian.detection == "[iz_a]"
    assert liouvillian.detect(basis.vectors["iz_a"]) == pytest.approx(1.0)
