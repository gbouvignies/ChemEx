"""Public topology and settings signatures for generic N-state models."""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_nst import (
    complete_edges,
    fork_edges,
    linear_edges,
    register,
)
from chemex.parameters.setting import ParamLocalSetting

type Edge = tuple[str, str]

EXPECTED_COMPLETE: dict[str, tuple[Edge, ...]] = {
    "abc": (("a", "b"), ("a", "c"), ("b", "c")),
    "abcd": (
        ("a", "b"),
        ("a", "c"),
        ("a", "d"),
        ("b", "c"),
        ("b", "d"),
        ("c", "d"),
    ),
    "abcde": (
        ("a", "b"),
        ("a", "c"),
        ("a", "d"),
        ("a", "e"),
        ("b", "c"),
        ("b", "d"),
        ("b", "e"),
        ("c", "d"),
        ("c", "e"),
        ("d", "e"),
    ),
    "abcdef": (
        ("a", "b"),
        ("a", "c"),
        ("a", "d"),
        ("a", "e"),
        ("a", "f"),
        ("b", "c"),
        ("b", "d"),
        ("b", "e"),
        ("b", "f"),
        ("c", "d"),
        ("c", "e"),
        ("c", "f"),
        ("d", "e"),
        ("d", "f"),
        ("e", "f"),
    ),
}

EXPECTED_LINEAR: dict[str, tuple[Edge, ...]] = {
    "abc": (("a", "b"), ("b", "c")),
    "abcd": (("a", "b"), ("b", "c"), ("c", "d")),
    "abcde": (("a", "b"), ("b", "c"), ("c", "d"), ("d", "e")),
    "abcdef": (
        ("a", "b"),
        ("b", "c"),
        ("c", "d"),
        ("d", "e"),
        ("e", "f"),
    ),
}

EXPECTED_FORK: dict[str, tuple[Edge, ...]] = {
    "abc": (("a", "b"), ("a", "c")),
    "abcd": (("a", "b"), ("a", "c"), ("a", "d")),
    "abcde": (("a", "b"), ("a", "c"), ("a", "d"), ("a", "e")),
    "abcdef": (
        ("a", "b"),
        ("a", "c"),
        ("a", "d"),
        ("a", "e"),
        ("a", "f"),
    ),
}


@pytest.mark.parametrize(
    ("builder", "expected"),
    (
        (complete_edges, EXPECTED_COMPLETE),
        (linear_edges, EXPECTED_LINEAR),
        (fork_edges, EXPECTED_FORK),
    ),
)
def test_topology_edge_builders_are_explicit(
    builder: Callable[[str], tuple[Edge, ...]],
    expected: dict[str, tuple[Edge, ...]],
) -> None:
    for states, edges in expected.items():
        assert builder(states) == edges


MODEL_TO_EDGES = {
    **{f"{len(states)}st": edges for states, edges in EXPECTED_COMPLETE.items()},
    "3st_triangle": EXPECTED_COMPLETE["abc"],
    **{f"{len(states)}st_linear": edges for states, edges in EXPECTED_LINEAR.items()},
    **{f"{len(states)}st_fork": edges for states, edges in EXPECTED_FORK.items()},
}


def _setting_signature(setting: ParamLocalSetting) -> tuple[object, ...]:
    name = setting.name_setting
    return (
        name.name,
        name.spin_system_part,
        name.conditions_part,
        name.allow_residue_specific,
        setting.value,
        setting.min,
        setting.max,
        setting.vary,
        setting.supports_estimation,
        setting.report_only,
        setting.expr,
    )


@pytest.mark.parametrize(("model_name", "edges"), MODEL_TO_EDGES.items())
def test_public_model_settings_match_declared_topology(
    model_name: str,
    edges: tuple[Edge, ...],
) -> None:
    register()
    settings = model_factory.create(model_name, Conditions())
    state_count = int(model_name[0])
    states = "abcdef"[:state_count]
    population_names = {f"p{state}" for state in states[1:]}
    kex_names = {f"kex_{left}{right}" for left, right in edges}
    rate_names = {
        f"k{source}{target}"
        for left, right in edges
        for source, target in ((left, right), (right, left))
    }

    assert population_names <= settings.keys()
    assert {name for name in settings if name.startswith("kex_")} == kex_names
    assert {
        name
        for name in settings
        if name.startswith("k") and not name.startswith("kex_")
    } == rate_names
    assert set(settings) == population_names | {"pa"} | kex_names | rate_names

    for name in population_names:
        setting = settings[name]
        assert (setting.value, setting.min, setting.max, setting.vary) == (
            0.02,
            0.0,
            1.0,
            True,
        )
        assert setting.name_setting.conditions_part == (
            "temperature",
            "p_total",
            "l_total",
        )
    for name in kex_names:
        setting = settings[name]
        assert (setting.value, setting.min, setting.max, setting.vary) == (
            200.0,
            0.0,
            1.0e6,
            True,
        )
        assert setting.name_setting.conditions_part == (
            "temperature",
            "p_total",
            "l_total",
        )
    assert "population_complement" in settings["pa"].expr
    assert (settings["pa"].min, settings["pa"].max) == (
        (-math.inf, math.inf) if state_count == 3 else (0.0, 1.0)
    )
    assert all("pair_rates" in settings[name].expr for name in rate_names)


def test_three_state_triangle_is_an_exact_signature_alias() -> None:
    register()
    canonical = model_factory.create("3st", Conditions())
    compatibility = model_factory.create("3st_triangle", Conditions())

    assert tuple(canonical) == tuple(compatibility)
    assert {
        name: _setting_signature(setting) for name, setting in canonical.items()
    } == {name: _setting_signature(setting) for name, setting in compatibility.items()}


@pytest.mark.parametrize(
    ("model_name", "edges"),
    (
        ("3st", EXPECTED_COMPLETE["abc"]),
        ("3st_triangle", EXPECTED_COMPLETE["abc"]),
        ("3st_linear", EXPECTED_LINEAR["abc"]),
        ("3st_fork", EXPECTED_FORK["abc"]),
    ),
)
def test_three_state_settings_preserve_frozen_base_metadata(
    model_name: str,
    edges: tuple[Edge, ...],
) -> None:
    """Freeze the public metadata shipped at the audited base commit."""
    register()
    settings = model_factory.create(model_name, Conditions())

    for name in ("pb", "pc"):
        setting = settings[name]
        assert _setting_signature(setting)[:-1] == (
            name,
            "",
            ("temperature", "p_total", "l_total"),
            True,
            0.02,
            0.0,
            1.0,
            True,
            False,
            False,
        )
    pa = settings["pa"]
    assert _setting_signature(pa)[:-1] == (
        "pa",
        "",
        ("temperature", "p_total", "l_total"),
        True,
        None,
        -math.inf,
        math.inf,
        False,
        False,
        False,
    )
    assert pa.dependencies == {"pb", "pc"}

    for left, right in edges:
        kex = settings[f"kex_{left}{right}"]
        assert _setting_signature(kex)[:-1] == (
            f"kex_{left}{right}",
            "",
            ("temperature", "p_total", "l_total"),
            True,
            200.0,
            0.0,
            1.0e6,
            True,
            False,
            False,
        )
        for source, target in ((left, right), (right, left)):
            rate = settings[f"k{source}{target}"]
            assert _setting_signature(rate)[:-1] == (
                f"k{source}{target}",
                "",
                ("temperature", "p_total", "l_total"),
                True,
                None,
                -math.inf,
                math.inf,
                False,
                False,
                False,
            )
            assert rate.dependencies == {
                f"kex_{left}{right}",
                f"p{left}",
                f"p{right}",
            }


def test_all_audited_public_model_names_are_registered() -> None:
    register()

    assert MODEL_TO_EDGES.keys() <= model_factory.set
    assert not {"4st_complete", "5st_complete", "6st_complete"} & model_factory.set
