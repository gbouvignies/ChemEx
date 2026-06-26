from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemex.configuration.experiment import (
    DetectionSettings,
    RelaxationSettings,
)
from chemex.experiments.catalog.shift_15n_sq import Shift15NSqSettings
from chemex.models.model import ModelSpec


@pytest.mark.parametrize("value", ["a", ["a"]])
def test_observed_state_normalizes_single_state_inputs(value: object) -> None:
    settings = DetectionSettings.model_validate({"observed_state": value})

    assert settings.observed_states == ("a",)
    assert settings.primary_state == "a"


def test_observed_state_normalizes_multiple_states() -> None:
    settings = DetectionSettings.model_validate(
        {"observed_state": ["a", "c"]},
    )

    assert settings.observed_states == ("a", "c")
    assert settings.primary_state == "a"
    assert settings.model_dump()["observed_state"] == ("a", "c")


def test_string_and_single_item_list_have_equivalent_detection() -> None:
    string_settings = RelaxationSettings.model_validate({"observed_state": "a"})
    list_settings = RelaxationSettings.model_validate({"observed_state": ["a"]})

    assert string_settings.get_detection_expression("[iz]") == "[iz_a]"
    assert (
        list_settings.get_detection_expression("[iz]")
        == string_settings.get_detection_expression("[iz]")
    )


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ([], "must contain at least one state"),
        (["a", "a"], "contains duplicate states: 'a'"),
        (["a", 1], "Every entry in 'observed_state' must be a string"),
        (1, "must be a state string or a non-empty list"),
    ],
)
def test_observed_state_rejects_invalid_values(
    value: object,
    message: str,
) -> None:
    with pytest.raises(ValidationError, match=message):
        DetectionSettings.model_validate({"observed_state": value})


def test_observed_state_rejects_unknown_state() -> None:
    with pytest.raises(
        ValidationError,
        match="Unknown observed state: 'z'",
    ):
        DetectionSettings.model_validate({"observed_state": ["a", "z"]})


def test_observed_state_rejects_state_missing_from_model() -> None:
    model = ModelSpec(name="2st", states="ab")

    with pytest.raises(
        ValidationError,
        match=r"Unknown observed state: 'c'.*Available states: 'a', 'b'",
    ):
        DetectionSettings.model_validate(
            {"observed_state": ["a", "c"]},
            context={"model": model},
        )


def test_observed_state_supports_states_from_larger_models() -> None:
    model = ModelSpec(name="6st", states="abcdef")

    settings = DetectionSettings.model_validate(
        {"observed_state": ["a", "f"]},
        context={"model": model},
    )

    assert settings.observed_states == ("a", "f")


def test_non_detection_experiment_rejects_observed_state_list() -> None:
    with pytest.raises(ValidationError, match="must be a state string"):
        Shift15NSqSettings.model_validate(
            {
                "name": "shift_15n_sq",
                "observed_state": ["a", "c"],
            },
        )


def test_detect_all_states_reports_replacement_syntax() -> None:
    model = ModelSpec(name="3st", states="abc")

    with pytest.raises(
        ValidationError,
        match=r'observed_state = \["a", "b", "c"\]',
    ):
        RelaxationSettings.model_validate(
            {"detect_all_states": True},
            context={"model": model},
        )


def test_disabled_detect_all_states_reports_removal() -> None:
    with pytest.raises(
        ValidationError,
        match=r"remove it.*single-state detection",
    ):
        RelaxationSettings.model_validate({"detect_all_states": False})


def test_enabled_cs_evolution_prior_reports_replacement_syntax() -> None:
    with pytest.raises(
        ValidationError,
        match=r"set 'start_state' to the observed state or states",
    ):
        RelaxationSettings.model_validate({"cs_evolution_prior": True})


def test_disabled_cs_evolution_prior_reports_removal() -> None:
    with pytest.raises(
        ValidationError,
        match=r"remove it.*equilibrium preparation",
    ):
        RelaxationSettings.model_validate({"cs_evolution_prior": False})


def test_start_state_omitted_uses_equilibrium_components() -> None:
    settings = RelaxationSettings.model_validate(
        {"observed_state": ["a", "c"]},
    )

    assert settings.start_states == ()
    assert settings.get_start_terms("iz") == ["iz"]


def test_start_state_accepts_single_state() -> None:
    settings = RelaxationSettings.model_validate(
        {
            "observed_state": ["a", "c"],
            "start_state": "b",
        },
    )

    assert settings.start_states == ("b",)
    assert settings.get_start_terms("iz") == ["iz_b"]


def test_start_state_accepts_multiple_states() -> None:
    settings = RelaxationSettings.model_validate(
        {
            "observed_state": "a",
            "start_state": ["a", "c"],
        },
    )

    assert settings.start_states == ("a", "c")
    assert settings.get_start_terms("2izsz", "-iz") == [
        "2izsz_a",
        "-iz_a",
        "2izsz_c",
        "-iz_c",
    ]


def test_empty_start_state_selects_equilibrium_preparation() -> None:
    settings = RelaxationSettings.model_validate({"start_state": []})

    assert settings.start_states == ()
    assert settings.get_start_terms("iz") == ["iz"]


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (["a", "a"], "'start_state' contains duplicate states"),
        (["a", "z"], "Unknown start state: 'z'"),
    ],
)
def test_start_state_rejects_invalid_values(
    value: object,
    message: str,
) -> None:
    with pytest.raises(ValidationError, match=message):
        RelaxationSettings.model_validate({"start_state": value})
