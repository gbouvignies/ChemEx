"""Scientific reference and regression tests for B-05 catalog physics, schemas, and CPMG guards."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Selection
from chemex.containers.data import Data
from chemex.experiments import experiment_types
from chemex.experiments.catalog import (
    cest_1hn_ap,
    cest_1hn_ip_ap,
    cest_15n_cw,
    cest_15n_tr,
    cest_ch3_1h_ip_ap,
    coscest_1hn_ip_ap,
    coscest_13c,
    cpmg_1hn_ap,
    cpmg_1hn_ap_0013,
    cpmg_13c_ip,
    cpmg_13co_ap,
    cpmg_15n_ip,
    cpmg_15n_ip_0013,
    cpmg_15n_tr,
    cpmg_15n_tr_0013,
    cpmg_ch3_1h_dq,
    cpmg_ch3_1h_sq,
    cpmg_ch3_1h_tq,
    cpmg_ch3_1h_tq_diff,
    cpmg_ch3_13c_h2c,
    cpmg_ch3_13c_h2c_0013,
    cpmg_chd2_1h_ap,
    cpmg_hn_dq_zq,
    dcest_13c,
    dcest_15n,
)
from chemex.experiments.experiment_types import (
    CpmgDutyCycleError,
    validate_cpmg_tau_cps,
)
from chemex.experiments.loader import register_experiments
from chemex.models.model import ModelSpec
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession, ensure_plugins_registered


def _two_state_model() -> ModelSpec:
    ensure_plugins_registered()
    return ModelSpec.from_name("2st")


def _make_cest_data(offsets: list[float]) -> Data:
    arr = np.array(offsets, dtype=float)
    return Data(
        exp=np.ones_like(arr),
        err=np.full_like(arr, 0.01),
        metadata=arr,
    )


def _two_state_hn_params() -> dict[str, float]:
    return {
        "pa": 0.8,
        "pb": 0.2,
        "kab": 30.0,
        "kba": 120.0,
        "cs_i_a": 8.4,
        "cs_i_b": 7.6,
        "cs_s_a": 120.0,
        "cs_s_b": 117.0,
        "j_is_a": -93.0,
        "j_is_b": -93.0,
        "r1_i_a": 2.0,
        "r1_i_b": 2.0,
        "r1_s_a": 1.5,
        "r1_s_b": 1.5,
        "r1a_is_a": 3.2,
        "r1a_is_b": 3.2,
        "r2_i_a": 15.0,
        "r2_i_b": 15.0,
        "r2a_i_a": 16.5,
        "r2a_i_b": 16.5,
        "etaz_i_a": 0.6,
        "etaz_i_b": 0.6,
        "etaxy_i_a": 1.5,
        "etaxy_i_b": 1.5,
    }


def _two_state_c_params() -> dict[str, float]:
    return {
        "pa": 0.8,
        "pb": 0.2,
        "kab": 25.0,
        "kba": 100.0,
        "cs_i_a": 55.0,
        "cs_i_b": 58.0,
        "r1_i_a": 1.8,
        "r1_i_b": 1.8,
        "r2_i_a": 18.0,
        "r2_i_b": 18.0,
    }


def _two_state_n_params() -> dict[str, float]:
    return {
        "pa": 0.85,
        "pb": 0.15,
        "kab": 15.0,
        "kba": 85.0,
        "cs_i_a": 120.0,
        "cs_i_b": 116.5,
        "r1_i_a": 1.4,
        "r1_i_b": 1.4,
        "r2_i_a": 12.0,
        "r2_i_b": 12.0,
    }


@pytest.mark.parametrize(
    ("module", "spin_name", "exp_base"),
    [
        (
            cest_1hn_ip_ap,
            "G23HN-N",
            {
                "name": "cest_1hn_ip_ap",
                "time_t1": 0.3,
                "carrier": 8.0,
                "b1_frq": 25.0,
                "b1_distribution": {"type": "homogeneous"},
                "d1": 0.6,
            },
        ),
        (
            coscest_1hn_ip_ap,
            "G23HN-N",
            {
                "name": "coscest_1hn_ip_ap",
                "time_t1": 0.3,
                "carrier": 8.0,
                "b1_frq": 25.0,
                "b1_distribution": {"type": "homogeneous"},
                "cos_n": 2,
                "sw": 400.0,
                "d1": 0.6,
            },
        ),
        (
            cest_ch3_1h_ip_ap,
            "L10HD1-CD1",
            {
                "name": "cest_ch3_1h_ip_ap",
                "time_t1": 0.3,
                "carrier": 0.8,
                "b1_frq": 25.0,
                "b1_distribution": {"type": "homogeneous"},
                "d1": 0.6,
                "taua": 0.0019,
            },
        ),
    ],
)
def test_ip_ap_cest_preserves_thermal_bath_under_start_state(
    module,
    spin_name: str,
    exp_base: dict[str, object],
) -> None:
    """IP-AP 1H CEST keeps the full identity bath populated when start_state selects longitudinal states."""
    model = _two_state_model()
    params = _two_state_hn_params()
    if spin_name.startswith("L10"):
        params = {
            **params,
            "cs_i_a": 0.9,
            "cs_i_b": 0.5,
            "cs_s_a": 24.0,
            "cs_s_b": 22.0,
            "j_is_a": 125.0,
            "j_is_b": 125.0,
        }
    data = _make_cest_data([-20000.0, -300.0, 0.0, 240.0])

    def run_calc(
        start_state: object,
        *,
        time_t1: float | None = None,
        b1_frq: float | None = None,
        d1: float | None = None,
    ) -> np.ndarray:
        exp_cfg = dict(exp_base)
        if start_state is not None:
            exp_cfg["start_state"] = start_state
        if time_t1 is not None:
            exp_cfg["time_t1"] = time_t1
        if b1_frq is not None:
            exp_cfg["b1_frq"] = b1_frq
        if d1 is not None:
            exp_cfg["d1"] = d1
        config = module.EXPERIMENT_TYPE.config_type.model_validate(
            {
                "experiment": exp_cfg,
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            },
            context={"model": model},
        )
        config.model = model
        config.experiment.model_name = model.name
        calc = module.EXPERIMENT_TYPE.create_profile_calculation(
            config,
            SpinSystem(name=spin_name),
        )
        calc.spectrometer.update(params)
        return calc.pulse_sequence.calculate(calc.spectrometer, data)

    # 1. Selecting all states ["a", "b"] must match start_state=None identically.
    calc_default = run_calc(None)
    calc_all = run_calc(["a", "b"])
    np.testing.assert_allclose(calc_all, calc_default, rtol=1e-13, atol=1e-15)

    # 2. Long-T1 zero-RF limit (b1_frq = 0, time_t1 = 25 s >> 1/R1):
    # Because the thermal identity operator sum_s I_{e,s} = 1 is conserved across
    # all states, start_state="a" and start_state="b" must relax to the exact same
    # Boltzmann steady state as start_state=None (not p_a or p_b times equilibrium).
    eq_default = run_calc(None, time_t1=25.0, b1_frq=0.0)
    eq_state_a = run_calc("a", time_t1=25.0, b1_frq=0.0)
    eq_state_b = run_calc("b", time_t1=25.0, b1_frq=0.0)
    assert np.abs(eq_default[0]) > 0.2
    np.testing.assert_allclose(eq_state_a, eq_default, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(eq_state_b, eq_default, rtol=1e-6, atol=1e-9)

    # Also under continuous RF (b1_frq = 25 Hz), long-T1 steady state is independent of start_state.
    rf_ss_default = run_calc(None, time_t1=25.0, b1_frq=25.0)
    rf_ss_a = run_calc("a", time_t1=25.0, b1_frq=25.0)
    rf_ss_b = run_calc("b", time_t1=25.0, b1_frq=25.0)
    np.testing.assert_allclose(rf_ss_a, rf_ss_default, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(rf_ss_b, rf_ss_default, rtol=1e-6, atol=1e-9)

    # 3. Exact affine thermal-bath superposition identity at finite time_t1 = 0.3 s:
    # start("a") + start("b") == [ie, iz_a, 0] + [ie, 0, iz_b]
    #                         == [ie, iz_a, iz_b] + [ie, 0, 0] == start(None) + start(d1=0).
    calc_a = run_calc("a")
    calc_b = run_calc("b")
    calc_d1_zero = run_calc(None, d1=0.0)
    np.testing.assert_allclose(
        calc_a + calc_b,
        calc_default + calc_d1_zero,
        rtol=1e-13,
        atol=1e-15,
    )

    # 4. Short-T1 limit (time_t1 = 1e-8 s): on the INEPT reference plane (-20000 Hz),
    # start_state="a" retains d1-recovered iz_a while start_state="b" suppresses it
    # (leaving only tiny b->a transfer during 2*taua).
    short_a = run_calc("a", time_t1=1e-8, b1_frq=0.0)
    short_b = run_calc("b", time_t1=1e-8, b1_frq=0.0)
    assert np.abs(short_a[0]) > 0.2
    assert np.abs(short_b[0]) < 0.05 * np.abs(short_a[0])


def test_coscest_13c_honors_start_state() -> None:
    """coscest_13c uses start_terms rather than unconditional get_equilibrium()."""
    model = _two_state_model()
    params = _two_state_c_params()
    data = _make_cest_data([-10000.0, -450.0, 0.0, 450.0])

    def run_calc(start_state: object, *, time_t1: float = 1e-6) -> np.ndarray:
        exp_cfg: dict[str, object] = {
            "name": "coscest_13c",
            "time_t1": time_t1,
            "carrier": 56.0,
            "b1_frq": 25.0,
            "b1_distribution": {"type": "homogeneous"},
            "cos_n": 2,
            "sw": 500.0,
        }
        if start_state is not None:
            exp_cfg["start_state"] = start_state
        config = coscest_13c.CosCest13CConfig.model_validate(
            {
                "experiment": exp_cfg,
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            },
            context={"model": model},
        )
        config.model = model
        config.experiment.model_name = model.name
        calc = coscest_13c.create_profile_calculation(config, SpinSystem(name="A10CA"))
        calc.spectrometer.update(params)
        return calc.pulse_sequence.calculate(calc.spectrometer, data)

    np.testing.assert_allclose(run_calc(["a", "b"]), run_calc(None), rtol=1e-13)
    assert np.all(np.abs(run_calc("a")) > 0.7)
    assert np.all(np.abs(run_calc("b")) < 1e-3)


@pytest.mark.parametrize(
    ("module", "spin_name", "params_fn", "exp_cfg"),
    [
        (
            dcest_15n,
            "G23N-HN",
            _two_state_n_params,
            {
                "name": "dcest_15n",
                "time_t1": 0.25,
                "time_equil": 0.05,
                "carrier": 118.0,
                "b1_eff": 20.0,
                "pw90": 40.0e-6,
                "b1_distribution": {"type": "homogeneous"},
                "sw": 800.0,
            },
        ),
        (
            dcest_13c,
            "A10CA",
            _two_state_c_params,
            {
                "name": "dcest_13c",
                "time_t1": 0.25,
                "time_equil": 0.05,
                "carrier": 56.0,
                "b1_eff": 25.0,
                "pw90": 15.0e-6,
                "b1_distribution": {"type": "homogeneous"},
                "sw": 800.0,
            },
        ),
        (
            coscest_13c,
            "A10CA",
            _two_state_c_params,
            {
                "name": "coscest_13c",
                "time_t1": 0.25,
                "time_equil": 0.05,
                "carrier": 56.0,
                "b1_frq": 25.0,
                "b1_distribution": {"type": "homogeneous"},
                "cos_n": 2,
                "sw": 800.0,
            },
        ),
    ],
)
def test_equilibration_delay_resets_offset_i_to_zero(
    module,
    spin_name: str,
    params_fn,
    exp_cfg: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """dcest_15n, dcest_13c, and coscest_13c evaluate d_eq with offset_i == 0.0."""
    model = _two_state_model()
    config = module.EXPERIMENT_TYPE.config_type.model_validate(
        {
            "experiment": exp_cfg,
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        },
        context={"model": model},
    )
    config.model = model
    config.experiment.model_name = model.name
    calc = module.EXPERIMENT_TYPE.create_profile_calculation(
        config,
        SpinSystem(name=spin_name),
    )
    calc.spectrometer.update(params_fn())
    data = _make_cest_data([-1500.0, -350.0, 250.0, 1850.0])

    observed_equil_offsets: list[float] = []
    orig_delays = calc.spectrometer.delays

    def tracked_delays(times):
        if times == 0.05:
            observed_equil_offsets.append(calc.spectrometer.offset_i)
        return orig_delays(times)

    monkeypatch.setattr(calc.spectrometer, "delays", tracked_delays)

    # Dirty offset_i prior to first calculation and also check repeated calculation.
    calc.spectrometer.offset_i = 4321.0
    first = calc.pulse_sequence.calculate(calc.spectrometer, data)
    second = calc.pulse_sequence.calculate(calc.spectrometer, data)

    assert observed_equil_offsets == [0.0, 0.0]
    np.testing.assert_allclose(second, first, rtol=0.0, atol=1e-15)


def _settings_cls(module):
    return module.EXPERIMENT_TYPE.config_type.model_fields["experiment"].annotation


@pytest.mark.parametrize(
    ("module", "carrier_fields"),
    [
        (cest_15n_cw, ("carrier", "carrier_dec")),
        (cest_15n_tr, ("carrier",)),
        (cest_1hn_ap, ("carrier",)),
        (cest_1hn_ip_ap, ("carrier",)),
        (cest_ch3_1h_ip_ap, ("carrier",)),
        (coscest_13c, ("carrier",)),
        (coscest_1hn_ip_ap, ("carrier",)),
        (cpmg_13c_ip, ("carrier",)),
        (cpmg_13co_ap, ("carrier",)),
        (cpmg_15n_ip, ("carrier",)),
        (cpmg_15n_ip_0013, ("carrier",)),
        (cpmg_15n_tr, ("carrier",)),
        (cpmg_15n_tr_0013, ("carrier",)),
        (cpmg_1hn_ap, ("carrier",)),
        (cpmg_1hn_ap_0013, ("carrier",)),
        (cpmg_ch3_13c_h2c, ("carrier",)),
        (cpmg_ch3_13c_h2c_0013, ("carrier",)),
        (cpmg_ch3_1h_dq, ("carrier",)),
        (cpmg_ch3_1h_sq, ("carrier",)),
        (cpmg_ch3_1h_tq, ("carrier",)),
        (cpmg_ch3_1h_tq_diff, ("carrier",)),
        (cpmg_chd2_1h_ap, ("carrier",)),
        (cpmg_hn_dq_zq, ("carrier_h", "carrier_n")),
    ],
)
def test_catalog_carrier_fields_use_chemical_shift_ppm_schema(
    module,
    carrier_fields: tuple[str, ...],
) -> None:
    """Every catalog carrier setting reports ppm chemical-shift units in JSON schema."""
    settings_cls = _settings_cls(module)
    properties = settings_cls.model_json_schema()["properties"]
    for field_name in carrier_fields:
        desc = properties[field_name]["description"]
        assert "ppm" in desc.lower(), f"{settings_cls.__name__}.{field_name}: {desc!r}"
        assert "hz" not in desc.lower(), (
            f"{settings_cls.__name__}.{field_name}: {desc!r}"
        )


@pytest.mark.parametrize(
    ("module", "exp_cfg", "bad_ncyc"),
    [
        (
            cpmg_15n_ip,
            {"name": "cpmg_15n_ip", "time_t2": 0.01, "carrier": 118.0, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_15n_ip_0013,
            {
                "name": "cpmg_15n_ip_0013",
                "time_t2": 0.01,
                "carrier": 118.0,
                "pw90": 50e-6,
                "ncyc_max": 100,
            },
            100.0,
        ),
        (
            cpmg_15n_tr,
            {"name": "cpmg_15n_tr", "time_t2": 0.01, "carrier": 118.0, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_15n_tr_0013,
            {
                "name": "cpmg_15n_tr_0013",
                "time_t2": 0.01,
                "carrier": 118.0,
                "pw90": 50e-6,
                "ncyc_max": 100,
            },
            100.0,
        ),
        (
            cpmg_13c_ip,
            {"name": "cpmg_13c_ip", "time_t2": 0.01, "carrier": 55.0, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_13co_ap,
            {
                "name": "cpmg_13co_ap",
                "time_t2": 0.01,
                "carrier": 176.0,
                "pw90": 50e-6,
                "pw180_ref": 100e-6,
            },
            100.0,
        ),
        (
            cpmg_1hn_ap,
            {"name": "cpmg_1hn_ap", "time_t2": 0.01, "carrier": 8.0, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_1hn_ap_0013,
            {
                "name": "cpmg_1hn_ap_0013",
                "time_t2": 0.01,
                "carrier": 8.0,
                "pw90": 50e-6,
                "ncyc_max": 100,
            },
            100.0,
        ),
        (
            cpmg_ch3_13c_h2c,
            {
                "name": "cpmg_ch3_13c_h2c",
                "time_t2": 0.01,
                "carrier": 20.0,
                "pw90": 50e-6,
                "tauc": 0.001,
            },
            100.0,
        ),
        (
            cpmg_ch3_13c_h2c_0013,
            {
                "name": "cpmg_ch3_13c_h2c_0013",
                "time_t2": 0.01,
                "carrier": 20.0,
                "pw90": 50e-6,
                "tauc": 0.001,
                "ncyc_max": 100,
            },
            100.0,
        ),
        (
            cpmg_ch3_1h_dq,
            {"name": "cpmg_ch3_1h_dq", "time_t2": 0.01, "carrier": 0.8, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_ch3_1h_sq,
            {
                "name": "cpmg_ch3_1h_sq",
                "time_t2": 0.01,
                "carrier": 0.8,
                "pw90": 50e-6,
                "ncyc_max": 100,
            },
            100.0,
        ),
        (
            cpmg_ch3_1h_tq,
            {"name": "cpmg_ch3_1h_tq", "time_t2": 0.01, "carrier": 0.8, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_ch3_1h_tq_diff,
            {
                "name": "cpmg_ch3_1h_tq_diff",
                "time_t2": 0.01,
                "carrier": 0.8,
                "pw90": 50e-6,
                "delta": 1.0e-3,
                "gradient": 10.0,
            },
            100.0,
        ),
        (
            cpmg_chd2_1h_ap,
            {"name": "cpmg_chd2_1h_ap", "time_t2": 0.01, "carrier": 0.8, "pw90": 50e-6},
            100.0,
        ),
        (
            cpmg_hn_dq_zq,
            {
                "name": "cpmg_hn_dq_zq",
                "time_t2": 0.01,
                "carrier_h": 8.0,
                "carrier_n": 118.0,
                "pw90_h": 50e-6,
                "pw90_n": 50e-6,
                "dq_flg": True,
            },
            40.0,
        ),
    ],
)
def test_cpmg_negative_interpulse_delay_raises_duty_cycle_error(
    module,
    exp_cfg: dict[str, object],
    bad_ncyc: float,
) -> None:
    """CPMG sequences reject negative inter-pulse delays with CpmgDutyCycleError."""
    settings = module.EXPERIMENT_TYPE.config_type.model_fields[
        "experiment"
    ].annotation.model_validate(exp_cfg)
    seq_cls = next(
        obj
        for name, obj in vars(module).items()
        if name.endswith("Sequence") and isinstance(obj, type)
    )
    sequence = seq_cls(settings)
    delay_fn = getattr(sequence, "_get_delays", None) or sequence._get_tau_cps

    with pytest.raises(CpmgDutyCycleError, match="negative delay"):
        delay_fn(np.array([0.0, bad_ncyc]))


def test_cpmg_build_rejects_duty_cycle_overflow_with_file_and_profile_context(
    tmp_path: Path,
) -> None:
    """experiment_types.build() fails fast with CpmgDutyCycleError carrying filename and spin system."""
    register_experiments()
    data_file = tmp_path / "10N-HN.out"
    data_file.write_text("0 1.0 0.02\n100 0.8 0.02\n", encoding="utf-8")
    exp_file = tmp_path / "bad_cpmg.toml"
    exp_file.write_text(
        """[experiment]
name = "cpmg_15n_ip"
time_t2 = 0.01
carrier = 118.0
pw90 = 50e-6

[conditions]
h_larmor_frq = 600.0

[data]
path = "."

[data.profiles]
10N-HN = "10N-HN.out"
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()
    session.set_model("2st")
    source = experiment_types.open(exp_file)

    with pytest.raises(CpmgDutyCycleError) as err_info:
        experiment_types.build(
            source,
            selection=Selection(include=None, exclude=None),
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

    assert err_info.value.filename == exp_file
    assert err_info.value.spin_system == "10N-H"
    assert err_info.value.experiment_type == "cpmg_15n_ip"
    assert err_info.value.tau_cp < 0.0


def test_cpmg_15n_ip_preserves_intentional_negative_t_neg_phase_correction() -> None:
    """Analytical negative delay t_neg = -2*pw90/pi in cpmg_15n_ip remains permitted when tau_cp >= 0."""
    model = _two_state_model()
    config = cpmg_15n_ip.Cpmg15NIpConfig.model_validate(
        {
            "experiment": {
                "name": "cpmg_15n_ip",
                "time_t2": 0.04,
                "carrier": 118.0,
                "pw90": 40.0e-6,
            },
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        },
        context={"model": model},
    )
    config.model = model
    config.experiment.model_name = model.name
    calc = cpmg_15n_ip.create_profile_calculation(config, SpinSystem(name="G23N-HN"))
    calc.spectrometer.update(_two_state_n_params())
    ncycs = np.array([0.0, 2.0, 10.0])
    data = Data(exp=np.ones(3), err=np.full(3, 0.01), metadata=ncycs)

    tau_cps, all_delays = calc.pulse_sequence._get_delays(ncycs)
    assert all(v > 0.0 for v in tau_cps.values())
    assert any(d < 0.0 for d in all_delays)  # t_neg = -2 * pw90 / pi

    intensities = calc.pulse_sequence.calculate(calc.spectrometer, data)
    assert intensities.shape == (3,)
    assert np.all(np.isfinite(intensities))


def test_cpmg_hn_dq_zq_uses_max_pw90_and_handles_reference_only_planes() -> None:
    """cpmg_hn_dq_zq subtracts 7/3 * max(pw90_n, pw90_h) and supports reference-only ncycs."""
    model = _two_state_model()
    pw90_n = 15.0e-6
    pw90_h = 48.0e-6  # pw90_h > pw90_n: composite pulse duration is 14/3 * pw90_h
    time_t2 = 0.04
    config = cpmg_hn_dq_zq.CpmgHNDqZqConfig.model_validate(
        {
            "experiment": {
                "name": "cpmg_hn_dq_zq",
                "time_t2": time_t2,
                "carrier_h": 8.2,
                "carrier_n": 118.0,
                "pw90_h": pw90_h,
                "pw90_n": pw90_n,
                "dq_flg": True,
            },
            "conditions": Conditions(h_larmor_frq=600.0).model_dump(),
            "data": {},
        },
        context={"model": model},
    )
    config.model = model
    config.experiment.model_name = model.name
    calc = cpmg_hn_dq_zq.create_profile_calculation(config, SpinSystem(name="G23N-HN"))
    calc.spectrometer.update(
        {
            **_two_state_n_params(),
            "cs_s_a": 8.3,
            "cs_s_b": 7.9,
            "r2mq_is_a": 14.0,
            "r2mq_is_b": 14.0,
            "mu_is_a": 2.0,
            "mu_is_b": 2.0,
        }
    )

    # 1. Verify exact total CPMG element timing: 4*ncyc*tau_cp + 2*ncyc*(14/3 * max(pw90)) == time_t2
    tau_cps = calc.pulse_sequence._get_tau_cps(np.array([0.0, 4.0, 12.0]))
    composite_pulse_duration = (14.0 / 3.0) * max(pw90_n, pw90_h)
    for ncyc, tau_cp in tau_cps.items():
        total_time = 4.0 * ncyc * tau_cp + 2.0 * ncyc * composite_pulse_duration
        np.testing.assert_allclose(total_time, time_t2, rtol=1e-14, atol=0.0)

    # 2. Verify ncyc that fits pw90_n alone (15 us) but overflows pw90_h (48 us) is rejected!
    # Maximum ncyc for pw90_h = 48 us at time_t2 = 0.04 is 0.04 / (4 * 7/3 * 48e-6) ≈ 89.28.
    with pytest.raises(CpmgDutyCycleError):
        calc.pulse_sequence._get_tau_cps(np.array([0.0, 92.0]))

    # 3. Reference-only ncycs = [0.0, 0.0] succeeds without zip(*...) ValueError.
    ref_data = Data(
        exp=np.ones(2),
        err=np.full(2, 0.01),
        metadata=np.array([0.0, 0.0]),
    )
    ref_calc = calc.pulse_sequence.calculate(calc.spectrometer, ref_data)
    assert ref_calc.shape == (2,)
    np.testing.assert_allclose(ref_calc, [0.85, 0.85], rtol=1e-12)


def test_validate_cpmg_tau_cps_clamps_sub_ulp_negative_tau_and_delta() -> None:
    """Sub-ULP negative floating-point noise in tau_cps and deltas is clamped to 0.0."""
    deltas = {0.0: -5e-17, 4.0: 1.0e-3}
    validated = validate_cpmg_tau_cps(
        "cpmg_15n_ip_0013",
        {4.0: -2e-17},
        deltas=deltas,
    )
    assert validated == {4.0: 0.0}
    assert deltas[0.0] == 0.0
    assert deltas[4.0] == 1.0e-3


def test_cpmg_ch3_1h_sq_validates_positive_ncyc_max_and_cycle_cap() -> None:
    """cpmg_ch3_1h_sq rejects non-positive ncyc_max and data ncyc > ncyc_max."""
    with pytest.raises(ValueError, match="greater than 0"):
        cpmg_ch3_1h_sq.CpmgCh31HSqSettings.model_validate(
            {
                "name": "cpmg_ch3_1h_sq",
                "time_t2": 0.04,
                "carrier": 0.8,
                "pw90": 10e-6,
                "ncyc_max": 0,
            }
        )

    settings = cpmg_ch3_1h_sq.CpmgCh31HSqSettings.model_validate(
        {
            "name": "cpmg_ch3_1h_sq",
            "time_t2": 0.04,
            "carrier": 0.8,
            "pw90": 10e-6,
            "ncyc_max": 10,
        }
    )
    seq = cpmg_ch3_1h_sq.CpmgCh31HSqSequence(settings)
    with pytest.raises(CpmgDutyCycleError, match="delta"):
        seq._get_delays(np.array([0.0, 12.0]))
