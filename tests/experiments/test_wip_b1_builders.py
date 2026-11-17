from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest
from pydantic import BaseModel

from chemex.experiments.catalog.wip import (
    cest_15n_test,
    relaxation_15n_r1rho,
    relaxation_15n_r1rho_eig,
)
from chemex.nmr.b1 import DephasingB1Model, GaussianB1Model
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem


@pytest.mark.parametrize(
    ("config_cls", "build_spectrometer", "experiment"),
    [
        (
            cest_15n_test.Cest15NTestConfig,
            cest_15n_test.build_spectrometer,
            {
                "name": "wip.cest_15n_test",
                "time_t1": 0.2,
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": 0.15,
                "b1_inh_res": 9,
            },
        ),
        (
            relaxation_15n_r1rho.Relaxation15NR1RhoConfig,
            relaxation_15n_r1rho.build_spectrometer,
            {
                "name": "wip.relaxation_15n_r1rho",
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": 0.15,
                "b1_inh_res": 9,
            },
        ),
        (
            relaxation_15n_r1rho_eig.Relaxation15NR1RhoConfig,
            relaxation_15n_r1rho_eig.build_spectrometer,
            {
                "name": "wip.relaxation_15n_r1rho_eig",
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": 0.15,
                "b1_inh_res": 9,
            },
        ),
    ],
)
def test_wip_builders_use_b1_inhomogeneity_api(
    monkeypatch: pytest.MonkeyPatch,
    config_cls: type[BaseModel],
    build_spectrometer: Callable[[Any, SpinSystem], object],
    experiment: dict[str, object],
) -> None:
    calls: list[tuple[float, GaussianB1Model | None]] = []

    def fake_set_b1_i_inhomogeneity(
        _self: Spectrometer,
        nominal: float,
        distribution: GaussianB1Model | None = None,
    ) -> None:
        calls.append((nominal, distribution))

    monkeypatch.setattr(
        Spectrometer,
        "set_b1_i_inhomogeneity",
        fake_set_b1_i_inhomogeneity,
    )
    config = config_cls.model_validate(
        {
            "experiment": experiment,
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )

    build_spectrometer(config, SpinSystem(name="G23N-HN"))

    assert len(calls) == 1
    nominal, distribution = calls[0]
    assert nominal == experiment["b1_frq"]
    assert isinstance(distribution, GaussianB1Model)
    assert distribution.scale == experiment["b1_inh_scale"]
    assert distribution.res == experiment["b1_inh_res"]


@pytest.mark.parametrize(
    ("config_cls", "build_spectrometer", "experiment", "expected_type", "checks"),
    [
        (
            cest_15n_test.Cest15NTestConfig,
            cest_15n_test.build_spectrometer,
            {
                "name": "wip.cest_15n_test",
                "time_t1": 0.2,
                "carrier": 118.0,
                "b1_frq": 25.0,
            },
            GaussianB1Model,
            {"scale": 0.1, "res": 11},
        ),
        (
            relaxation_15n_r1rho.Relaxation15NR1RhoConfig,
            relaxation_15n_r1rho.build_spectrometer,
            {
                "name": "wip.relaxation_15n_r1rho",
                "carrier": 118.0,
                "b1_frq": 25.0,
            },
            DephasingB1Model,
            {},
        ),
        (
            relaxation_15n_r1rho_eig.Relaxation15NR1RhoConfig,
            relaxation_15n_r1rho_eig.build_spectrometer,
            {
                "name": "wip.relaxation_15n_r1rho_eig",
                "carrier": 118.0,
                "b1_frq": 25.0,
            },
            DephasingB1Model,
            {},
        ),
    ],
)
def test_wip_builders_normalize_default_legacy_b1_settings(
    monkeypatch: pytest.MonkeyPatch,
    config_cls: type[BaseModel],
    build_spectrometer: Callable[[Any, SpinSystem], object],
    experiment: dict[str, object],
    expected_type: type[GaussianB1Model | DephasingB1Model],
    checks: dict[str, object],
) -> None:
    calls: list[tuple[float, GaussianB1Model | DephasingB1Model | None]] = []

    def fake_set_b1_i_inhomogeneity(
        _self: Spectrometer,
        nominal: float,
        distribution: GaussianB1Model | DephasingB1Model | None = None,
    ) -> None:
        calls.append((nominal, distribution))

    monkeypatch.setattr(
        Spectrometer,
        "set_b1_i_inhomogeneity",
        fake_set_b1_i_inhomogeneity,
    )
    config = config_cls.model_validate(
        {
            "experiment": experiment,
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )

    build_spectrometer(config, SpinSystem(name="G23N-HN"))

    assert len(calls) == 1
    nominal, distribution = calls[0]
    assert nominal == experiment["b1_frq"]
    assert isinstance(distribution, expected_type)
    for name, value in checks.items():
        assert getattr(distribution, name) == value


def test_wip_cest_15n_test_auto_res_uses_time_t1() -> None:
    def _build_res(time_t1: float) -> int:
        config = cest_15n_test.Cest15NTestConfig.model_validate(
            {
                "experiment": {
                    "name": "wip.cest_15n_test",
                    "time_t1": time_t1,
                    "carrier": 118.0,
                    "b1_frq": 25.0,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                },
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            }
        )
        spectrometer = cest_15n_test.build_spectrometer(
            config, SpinSystem(name="G23N-HN")
        )
        return spectrometer.b1_i_distribution.values.size

    assert _build_res(1.5) > _build_res(0.1)
