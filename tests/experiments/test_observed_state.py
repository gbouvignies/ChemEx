from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemex.experiments.catalog.cest_1hn_ap import Cest1HnApSettings
from chemex.experiments.catalog.cest_1hn_ip_ap import Cest1HnIpApSettings
from chemex.experiments.catalog.cest_ch3_1h_ip_ap import CestCh31HIpApSettings
from chemex.experiments.catalog.coscest_1hn_ip_ap import CosCest1HnIpApSettings
from chemex.experiments.catalog.cpmg_1hn_ap import Cpmg1HnApSettings
from chemex.experiments.catalog.cpmg_1hn_ap_0013 import Cpmg1HnAp0013Settings
from chemex.experiments.catalog.cpmg_15n_ip import Cpmg15NIpSettings
from chemex.experiments.catalog.cpmg_ch3_1h_sq import CpmgCh31HSqSettings
from chemex.experiments.catalog.dcest_15n import DCest15NSettings


def _settings(
    observed_state: object,
    *,
    start_state: object = None,
) -> Cpmg15NIpSettings:
    return Cpmg15NIpSettings.model_validate(
        {
            "name": "cpmg_15n_ip",
            "time_t2": 0.04,
            "carrier": 118.0,
            "pw90": 40.0e-6,
            "observed_state": observed_state,
            "start_state": start_state,
        },
    )


def test_cpmg_observed_state_string_remains_unchanged() -> None:
    assert _settings("a").detection == "[iz_a]"


def test_cpmg_observed_state_list_sums_selected_states() -> None:
    assert _settings(["a", "c"]).detection == "[iz_a] + [iz_c]"


def test_cpmg_start_state_list_selects_starting_components() -> None:
    assert _settings("a", start_state=["a", "c"]).start_terms == ["iz_a", "iz_c"]


@pytest.mark.parametrize(
    ("settings_cls", "data"),
    [
        (
            Cest1HnIpApSettings,
            {
                "name": "cest_1hn_ip_ap",
                "time_t1": 0.2,
                "carrier": 8.0,
                "d1": 0.5,
            },
        ),
        (
            CestCh31HIpApSettings,
            {
                "name": "cest_ch3_1h_ip_ap",
                "time_t1": 0.2,
                "carrier": 0.8,
                "d1": 0.5,
            },
        ),
        (
            CosCest1HnIpApSettings,
            {
                "name": "coscest_1hn_ip_ap",
                "time_t1": 0.2,
                "carrier": 8.0,
                "sw": 1000.0,
                "cos_n": 2,
                "d1": 0.5,
            },
        ),
    ],
)
def test_start_state_is_applied_to_identity_component_preparation(
    settings_cls,
    data: dict[str, object],
) -> None:
    settings = settings_cls.model_validate(
        {**data, "start_state": ["a", "c"]},
    )

    assert settings.start_terms == ["ie_a", "ie_c"]


@pytest.mark.parametrize(
    ("settings_cls", "data"),
    [
        (
            Cest1HnApSettings,
            {
                "name": "cest_1hn_ap",
                "time_t1": 0.2,
                "carrier": 8.0,
                "b1_frq": 25.0,
            },
        ),
        (
            Cpmg1HnApSettings,
            {
                "name": "cpmg_1hn_ap",
                "time_t2": 0.04,
                "carrier": 8.0,
                "pw90": 10.0e-6,
            },
        ),
        (
            Cpmg1HnAp0013Settings,
            {
                "name": "cpmg_1hn_ap_0013",
                "time_t2": 0.04,
                "carrier": 8.0,
                "pw90": 10.0e-6,
                "ncyc_max": 20,
            },
        ),
        (
            CpmgCh31HSqSettings,
            {
                "name": "cpmg_ch3_1h_sq",
                "time_t2": 0.04,
                "carrier": 0.8,
                "pw90": 10.0e-6,
                "ncyc_max": 20,
            },
        ),
    ],
)
def test_legacy_non_equilibrium_defaults_are_preserved(
    settings_cls,
    data: dict[str, object],
) -> None:
    settings = settings_cls.model_validate(data)

    assert settings.start_states == ("a",)
    assert all(term.endswith("_a") for term in settings.start_terms)


def test_legacy_non_equilibrium_default_can_be_overridden_to_equilibrium() -> None:
    settings = Cpmg1HnApSettings.model_validate(
        {
            "name": "cpmg_1hn_ap",
            "time_t2": 0.04,
            "carrier": 8.0,
            "pw90": 10.0e-6,
            "start_state": [],
        },
    )

    assert settings.start_states == ()
    assert settings.start_terms == ["2izsz"]


def test_disabled_cs_evolution_prior_reports_equilibrium_override() -> None:
    with pytest.raises(
        ValidationError,
        match=r"start_state = \[\].*equilibrium preparation",
    ):
        Cpmg1HnApSettings.model_validate(
            {
                "name": "cpmg_1hn_ap",
                "time_t2": 0.04,
                "carrier": 8.0,
                "pw90": 10.0e-6,
                "cs_evolution_prior": False,
            },
        )


def test_legacy_non_equilibrium_default_tracks_observed_states() -> None:
    settings = Cpmg1HnApSettings.model_validate(
        {
            "name": "cpmg_1hn_ap",
            "time_t2": 0.04,
            "carrier": 8.0,
            "pw90": 10.0e-6,
            "observed_state": ["a", "b"],
        },
    )

    assert settings.start_states == ("a", "b")
    assert settings.start_terms == ["2izsz_a", "2izsz_b"]


def test_dcest_hd_equilibrium_override_is_not_treated_as_omitted() -> None:
    data = {
        "name": "dcest_15n",
        "time_t1": 0.2,
        "sw": 1000.0,
        "carrier": 118.0,
        "b1_eff": 25.0,
        "pw90": 40.0e-6,
        "model_name": "2st_hd",
    }

    default = DCest15NSettings.model_validate(data)
    equilibrium = DCest15NSettings.model_validate({**data, "start_state": []})

    assert default.start_terms == ["iz_a"]
    assert equilibrium.start_terms == ["iz"]
