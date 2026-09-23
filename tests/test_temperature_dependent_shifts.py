"""Focused scientific regressions for reference-centered ``.tc`` shifts."""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.method_plan import MethodFormatError
from chemex.configuration.methods import Method, read_method_plan
from chemex.configuration.parameters import DefaultSetting
from chemex.containers.data import Data
from chemex.experiments.catalog.cest_15n import Cest15NConfig
from chemex.experiments.catalog.cest_15n import (
    create_profile_calculation as create_cest_profile_calculation,
)
from chemex.experiments.catalog.cpmg_15n_ip import Cpmg15NIpConfig
from chemex.experiments.catalog.cpmg_15n_ip import (
    create_profile_calculation as create_cpmg_profile_calculation,
)
from chemex.experiments.catalog.shift_15n_sq import (
    Shift15NSqSequence,
    Shift15NSqSettings,
)
from chemex.nmr.basis import Basis
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.database import TemperatureReferenceConfigurationError
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ModelDerivationOverrideError,
    ParameterRole,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.run_info import serialize_parameter_file
from chemex.runtime import AnalysisSession
from chemex.typing import Array

FULL_CONDITIONS = {
    "h_larmor_frq": 800.0,
    "p_total": 1.0e-4,
    "l_total": 2.0e-4,
    "d2o": 0.1,
    "label": ("2h",),
}


def _defaults(**values: float) -> list[tuple[ParamName, DefaultSetting]]:
    return [
        (ParamName.from_section(name), DefaultSetting(value))
        for name, value in values.items()
    ]


def _build_shift_session(
    temperatures: tuple[float, ...],
    defaults: list[tuple[ParamName, DefaultSetting]],
    *,
    model_name: str = "2st.tc",
    basis_type: str = "ixy",
    fitted: tuple[str, ...] = ("r2",),
) -> tuple[AnalysisSession, dict[float, dict[str, str]], dict[str, float]]:
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(
        type=basis_type,  # ty: ignore[invalid-argument-type]
        spin_system="nh",
        model=session.model.spec,
    )
    parameter_ids: dict[float, dict[str, str]] = {}
    for temperature in temperatures:
        config = SimpleNamespace(
            conditions=Conditions(temperature=temperature, **FULL_CONDITIONS),
            to_be_fitted=SimpleNamespace(
                rates=list(fitted),
                model_free=list(fitted),
            ),
        )
        parameter_ids[temperature] = session.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=SpinSystem.from_name("G23N-HN"),
        )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error
    )
    session.parameters.set_defaults(defaults)
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    required_ids = {
        param_id
        for mapping in parameter_ids.values()
        for local_name, param_id in mapping.items()
        if local_name.startswith("cs_")
    }
    return session, parameter_ids, dict(session.resolve_current_values(required_ids))


def _resolved_lines(
    parameter_ids: dict[float, dict[str, str]],
    resolved: dict[str, float],
) -> Array:
    return np.asarray(
        [
            (
                resolved[parameter_ids[temperature]["cs_i_a"]],
                resolved[parameter_ids[temperature]["cs_i_b"]]
                - resolved[parameter_ids[temperature]["cs_i_a"]],
                resolved[parameter_ids[temperature]["cs_i_b"]],
            )
            for temperature in parameter_ids
        ]
    )


def test_reference_centered_cs_and_dw_lines_are_independent() -> None:
    temperatures = (0.0, 10.0, 25.0, 40.0)
    _, parameter_ids, resolved = _build_shift_session(
        temperatures,
        _defaults(TREF=25.0, CS0_A=8.0, CS1_A=-0.01, DW0_AB=1.5, DW1_AB=0.02),
    )
    np.testing.assert_allclose(
        _resolved_lines(parameter_ids, resolved),
        np.asarray(
            [
                (8.25, 1.0, 9.25),
                (8.15, 1.2, 9.35),
                (8.0, 1.5, 9.5),
                (7.85, 1.8, 9.65),
            ]
        ),
        rtol=0.0,
        atol=1.0e-12,
    )


def test_zero_slopes_match_the_temperature_independent_model() -> None:
    temperatures = (0.0, 25.0, 40.0)
    _, tc_ids, tc_values = _build_shift_session(
        temperatures,
        _defaults(CS0_A=8.0, CS1_A=0.0, DW0_AB=1.5, DW1_AB=0.0),
    )
    _, ordinary_ids, ordinary_values = _build_shift_session(
        temperatures,
        _defaults(CS_A=8.0, DW_AB=1.5),
        model_name="2st",
    )
    np.testing.assert_allclose(
        _resolved_lines(tc_ids, tc_values)[:, (0, 2)],
        _resolved_lines(ordinary_ids, ordinary_values)[:, (0, 2)],
        rtol=0.0,
        atol=1.0e-12,
    )


def test_reference_change_transformation_preserves_resolved_shifts() -> None:
    temperatures = (5.0, 17.5, 40.0)
    cs0, cs1, dw0, dw1 = 118.3, -0.012, 2.1, 0.018
    _, ids_25, values_25 = _build_shift_session(
        temperatures,
        _defaults(TREF=25.0, CS0_A=cs0, CS1_A=cs1, DW0_AB=dw0, DW1_AB=dw1),
    )
    new_reference = 10.0
    _, ids_10, values_10 = _build_shift_session(
        temperatures,
        _defaults(
            TREF=new_reference,
            CS0_A=cs0 + cs1 * (new_reference - 25.0),
            CS1_A=cs1,
            DW0_AB=dw0 + dw1 * (new_reference - 25.0),
            DW1_AB=dw1,
        ),
    )
    np.testing.assert_allclose(
        _resolved_lines(ids_25, values_25),
        _resolved_lines(ids_10, values_10),
        rtol=0.0,
        atol=1.0e-12,
    )


def test_coefficients_share_across_conditions_but_not_nuclei_or_residues() -> None:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    basis = Basis(type="ixyzsxyz", spin_system="nh", model=session.model.spec)

    def create(temperature: float, field: float, residue: int) -> dict[str, str]:
        config = SimpleNamespace(
            conditions=Conditions(
                temperature=temperature,
                h_larmor_frq=field,
                p_total=1.0e-4,
                l_total=2.0e-4,
                d2o=0.1,
                label=("2h",),
            ),
            to_be_fitted=SimpleNamespace(rates=["r2"], model_free=[]),
        )
        return session.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=SpinSystem.from_name(f"G{residue}N-HN"),
        )

    low = create(15.0, 600.0, 23)
    high = create(35.0, 900.0, 23)
    other_field = create(15.0, 900.0, 23)
    other_residue = create(15.0, 600.0, 24)
    assert low["cs_i_a"] != high["cs_i_a"]
    assert low["cs_i_a"] == other_field["cs_i_a"]
    assert low["cs_i_a"] != other_residue["cs_i_a"]
    assert session.parameter_factory.try_seal_definitions()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    assert sum(item.name == "TREF" for item in definitions) == 1
    for name in ("CS0_A", "CS1_A", "DW0_AB", "DW1_AB"):
        spin_systems = {
            item.spin_system_name for item in definitions if item.name == name
        }
        assert spin_systems == {"G23N", "G23H", "G24N", "G24H"}


def test_multistate_direct_shift_observables_use_a_centered_composition() -> None:
    temperature = 30.0
    _, parameter_ids, resolved = _build_shift_session(
        (temperature,),
        _defaults(
            TREF=20.0,
            CS0_A=8.0,
            CS1_A=0.01,
            DW0_AB=1.0,
            DW1_AB=-0.02,
            DW0_AC=2.0,
            DW1_AC=0.03,
        ),
        model_name="3st.tc",
    )
    model_session = AnalysisSession.create()
    model_session.set_model("3st.tc")
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem.from_name("G23N-HN"),
        Basis(type="ixy", spin_system="nh", model=model_session.model.spec),
        Conditions(h_larmor_frq=800.0, temperature=temperature),
    )
    ids = parameter_ids[temperature]
    spectrometer.update(
        {
            "r2_i_a": 4.0,
            "r2_i_b": 4.0,
            "r2_i_c": 4.0,
            "cs_i_a": resolved[ids["cs_i_a"]],
            "cs_i_b": resolved[ids["cs_i_b"]],
            "cs_i_c": resolved[ids["cs_i_c"]],
            "kab": 0.0,
            "kba": 0.0,
            "kac": 0.0,
            "kca": 0.0,
            "kbc": 0.0,
            "kcb": 0.0,
            "pa": 1.0,
            "pb": 0.0,
            "pc": 0.0,
        }
    )
    data = Data(
        exp=np.asarray([0.0]), err=np.asarray([1.0]), metadata=np.asarray([0.0])
    )
    observed = np.asarray(
        [
            Shift15NSqSequence(
                Shift15NSqSettings(name="shift_15n_sq", observed_state=state)
            ).calculate(spectrometer, data)[0]
            for state in ("b", "c")
        ]
    )
    np.testing.assert_allclose(observed, np.asarray([8.9, 10.4]), atol=1.0e-12)


@pytest.mark.parametrize("model_name", ("2st_eyring.tc", "2st.mf.tc"))
def test_tc_composes_with_representative_model_extensions(model_name: str) -> None:
    _, ids, resolved = _build_shift_session(
        (20.0, 30.0),
        _defaults(CS0_A=8.0, CS1_A=-0.01),
        model_name=model_name,
    )
    assert tuple(
        resolved[ids[temperature]["cs_i_a"]] for temperature in (20.0, 30.0)
    ) == pytest.approx((8.05, 7.95), abs=1.0e-12)


def test_direct_shift_authority_expands_cs_a_to_both_coefficients() -> None:
    session, parameter_ids, _ = _build_shift_session(
        (15.0, 35.0),
        _defaults(CS0_A=118.0, CS1_A=-0.01, DW0_AB=2.0, DW1_AB=0.01),
        fitted=("cs_i_a",),
    )
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    by_name = {item.name: item.param_id for item in model.definitions}
    required = {
        param_id for mapping in parameter_ids.values() for param_id in mapping.values()
    }
    parameterization = session.compile_parameterization(
        Method(fix=("DW0_AB", "DW1_AB")), required
    )
    assert parameterization.role(by_name["CS0_A"]) is ParameterRole.FIT
    assert parameterization.role(by_name["CS1_A"]) is ParameterRole.FIT
    assert parameterization.role(by_name["TREF"]) is ParameterRole.FIX
    assert parameterization.role(by_name["CS_A"]) is ParameterRole.DERIVED


@pytest.mark.parametrize(
    "method",
    (
        Method(fit=("TREF",)),
        Method(fix=("TREF",)),
        Method(constraints=("[TREF] = 25.0",)),
        Method(constraints=("[CS0_A] = [TREF]",)),
    ),
)
def test_tref_is_protected_from_method_operations(method: Method) -> None:
    session, parameter_ids, _ = _build_shift_session(
        (15.0, 35.0), _defaults(TREF=20.0, CS0_A=118.0, CS1_A=-0.01)
    )
    with pytest.raises(ModelDerivationOverrideError):
        session.compile_parameterization(method, set(parameter_ids[15.0].values()))


@pytest.mark.parametrize(
    ("selector", "guidance"),
    (
        ("CS_A", "use CS0_A and CS1_A"),
        ("DW_AB", "use DW0_AB and DW1_AB"),
        ("CS_B", "use CS0_A, CS1_A, DW0_AB, and DW1_AB"),
    ),
)
def test_derived_shift_overrides_name_independent_coordinates(
    selector: str,
    guidance: str,
) -> None:
    session, parameter_ids, _ = _build_shift_session(
        (15.0, 35.0),
        _defaults(CS0_A=118.0, CS1_A=-0.01, DW0_AB=2.0, DW1_AB=0.02),
    )
    required = {
        param_id for mapping in parameter_ids.values() for param_id in mapping.values()
    }
    with pytest.raises(ModelDerivationOverrideError, match=guidance):
        session.compile_parameterization(Method(fit=(selector,)), required)


def test_tref_configuration_must_be_unscoped_scalar() -> None:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=Conditions(h_larmor_frq=800.0, temperature=25.0),
        to_be_fitted=SimpleNamespace(rates=["r2"], model_free=[]),
    )
    session.parameter_factory.create_parameters(
        config, basis=basis, spin_system=SpinSystem.from_name("G23N-HN")
    )
    with pytest.raises(TemperatureReferenceConfigurationError, match="unqualified"):
        session.parameters.set_defaults(
            [(ParamName("TREF", SpinSystem.from_name("G23N")), DefaultSetting(25.0))]
        )
    with pytest.raises(TemperatureReferenceConfigurationError, match="bounds"):
        session.parameters.set_defaults(
            [(ParamName("TREF"), DefaultSetting(25.0, 0.0, 50.0))]
        )
    with pytest.raises(TemperatureReferenceConfigurationError, match="greater"):
        session.parameters.set_defaults(_defaults(TREF=-273.15))


def test_canonical_temperature_coefficients_support_constraints_and_grid(
    tmp_path: Path,
) -> None:
    session, parameter_ids, _ = _build_shift_session(
        (15.0, 25.0, 35.0),
        _defaults(CS0_A=118.0, CS1_A=-0.01, DW0_AB=2.0, DW1_AB=0.02),
        fitted=("cs_i_a",),
    )
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    by_name = {item.name: item.param_id for item in model.definitions}
    required = {
        param_id for mapping in parameter_ids.values() for param_id in mapping.values()
    }
    constrained = session.compile_parameterization(
        Method(constraints=("[CS1_A] = -0.5 * [DW1_AB]",)), required
    )
    values = constrained.resolve(
        constrained.frame_from_snapshot(session.analysis_values.snapshot())
    )
    assert values[by_name["CS1_A"]] == pytest.approx(-0.5 * values[by_name["DW1_AB"]])

    method = tmp_path / "method.toml"
    method.write_text(
        """FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = ["[DW1_AB] = values(-0.01, 0.0, 0.01)"]
""",
        encoding="utf-8",
    )
    session.validate_method_plan(read_method_plan((method,)))

    protected = tmp_path / "protected.toml"
    protected.write_text(
        """FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = ["[TREF] = values(15.0, 20.0)"]
""",
        encoding="utf-8",
    )
    with pytest.raises(MethodFormatError, match="model-owned"):
        session.validate_method_plan(read_method_plan((protected,)))


def _calculate_cpmg(
    *, temperature: float, tref: float, cs0: float, cs1: float, dw0: float, dw1: float
) -> Array:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    config = Cpmg15NIpConfig.model_validate(
        {
            "model": session.model.spec,
            "experiment": {
                "name": "cpmg_15n_ip",
                "time_t2": 0.03,
                "carrier": 118.0,
                "pw90": 40.0e-6,
                "time_equil": 0.002,
            },
            "conditions": {
                "h_larmor_frq": 800.0,
                "temperature": temperature,
                "label": ["2h"],
            },
            "data": {},
        },
        context={"model": session.model.spec},
    )
    calculation = create_cpmg_profile_calculation(
        config, SpinSystem.from_name("G23N-HN")
    )
    parameter_ids = session.parameter_factory.create_parameters(
        config,
        basis=calculation.spectrometer.basis,
        spin_system=calculation.spectrometer.spin_system,
    )
    assert session.parameter_factory.try_seal_definitions()
    session.parameters.set_defaults(
        _defaults(
            TREF=tref,
            CS0_A=cs0,
            CS1_A=cs1,
            DW0_AB=dw0,
            DW1_AB=dw1,
            PB=0.1,
            KEX_AB=500.0,
        )
    )
    assert session.try_build_analysis_values()
    resolved = session.resolve_current_values(set(parameter_ids.values()))
    calculation.spectrometer.update(
        {name: resolved[param_id] for name, param_id in parameter_ids.items()}
    )
    data = Data(
        exp=np.ones(4),
        err=np.ones(4),
        metadata=np.asarray([0.0, -1.0, 2.0, 8.0]),
    )
    return calculation.pulse_sequence.calculate(calculation.spectrometer, data)


def test_reference_transformation_preserves_cpmg_observables() -> None:
    original = _calculate_cpmg(
        temperature=40.0, tref=25.0, cs0=118.0, cs1=-0.02, dw0=2.0, dw1=0.01
    )
    transformed = _calculate_cpmg(
        temperature=40.0, tref=10.0, cs0=118.3, cs1=-0.02, dw0=1.85, dw1=0.01
    )
    np.testing.assert_allclose(original, transformed, rtol=2.0e-13, atol=1.0e-14)
    np.testing.assert_allclose(
        original,
        np.asarray([0.89636109, 0.20194106, 0.22946011, 0.50416439]),
        rtol=2.0e-7,
        atol=1.0e-9,
    )


def _calculate_cest(
    *, temperature: float, tref: float, cs0: float, cs1: float, dw0: float, dw1: float
) -> Array:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    config = Cest15NConfig.model_validate(
        {
            "model": session.model.spec,
            "experiment": {
                "name": "cest_15n",
                "time_t1": 0.25,
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": 0.0,
                "b1_inh_res": 1,
            },
            "conditions": {
                "h_larmor_frq": 800.0,
                "temperature": temperature,
                "label": ["2h"],
            },
            "data": {},
        },
        context={"model": session.model.spec},
    )
    calculation = create_cest_profile_calculation(
        config, SpinSystem.from_name("G23N-HN")
    )
    parameter_ids = session.parameter_factory.create_parameters(
        config,
        basis=calculation.spectrometer.basis,
        spin_system=calculation.spectrometer.spin_system,
    )
    assert session.parameter_factory.try_seal_definitions()
    session.parameters.set_defaults(
        _defaults(
            TREF=tref,
            CS0_A=cs0,
            CS1_A=cs1,
            DW0_AB=dw0,
            DW1_AB=dw1,
            PB=0.1,
            KEX_AB=500.0,
        )
    )
    assert session.try_build_analysis_values()
    resolved = session.resolve_current_values(set(parameter_ids.values()))
    calculation.spectrometer.update(
        {name: resolved[param_id] for name, param_id in parameter_ids.items()}
    )
    data = Data(
        exp=np.ones(5),
        err=np.ones(5),
        metadata=np.asarray([2.0e4, -200.0, 0.0, 200.0, 400.0]),
    )
    return calculation.pulse_sequence.calculate(calculation.spectrometer, data)


def test_reference_transformation_preserves_cest_observables() -> None:
    original = _calculate_cest(
        temperature=40.0, tref=25.0, cs0=118.0, cs1=-0.02, dw0=2.0, dw1=0.01
    )
    transformed = _calculate_cest(
        temperature=40.0, tref=10.0, cs0=118.3, cs1=-0.02, dw0=1.85, dw1=0.01
    )
    without_slopes = _calculate_cest(
        temperature=40.0, tref=25.0, cs0=118.0, cs1=0.0, dw0=2.0, dw1=0.0
    )
    np.testing.assert_allclose(original, transformed, rtol=2.0e-13, atol=1.0e-14)
    assert not np.allclose(original[1:], without_slopes[1:], rtol=1.0e-6, atol=1.0e-8)


def test_longitudinal_basis_does_not_create_unused_shift_coefficients() -> None:
    session = AnalysisSession.create()
    session.set_model("2st.tc")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=Conditions(h_larmor_frq=800.0, temperature=25.0),
        to_be_fitted=SimpleNamespace(rates=["r1"], model_free=[]),
    )
    session.parameter_factory.create_parameters(
        config, basis=basis, spin_system=SpinSystem.from_name("G23N-HN")
    )
    assert session.parameter_factory.try_seal_definitions()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    assert not {item.name for item in definitions} & {
        "TREF",
        "CS0_A",
        "CS1_A",
        "DW0_AB",
        "DW1_AB",
    }


def test_restart_serialization_keeps_tref_and_canonical_coordinates(
    tmp_path: Path,
) -> None:
    session, _parameter_ids, _resolved = _build_shift_session(
        (15.0, 25.0, 35.0),
        _defaults(TREF=20.0, CS0_A=118.2, CS1_A=0.015, DW0_AB=2.0, DW1_AB=0.012),
    )
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    text = serialize_parameter_file(
        model, session.analysis_values.snapshot(), state_kind="restart"
    )
    assert '"TREF" = 20.0' in text
    assert all(name in text for name in ("CS0_A", "CS1_A", "DW0_AB", "DW1_AB"))
    restart = tmp_path / "restart.toml"
    restart.write_text(text, encoding="utf-8")
    with restart.open("rb") as file:
        assert tomllib.load(file)["GLOBAL"]["TREF"] == 20.0


def _write_shift_fit_case(root: Path) -> tuple[list[Path], Path, Path]:
    temperatures = (10.0, 20.0, 30.0)
    cs_values = (118.05, 118.2, 118.35)
    state_b_values = (119.93, 120.2, 120.47)
    experiments: list[Path] = []
    for temperature, cs_a, cs_b in zip(
        temperatures, cs_values, state_b_values, strict=True
    ):
        for state, value in (("a", cs_a), ("b", cs_b)):
            data = root / f"shift_{state}_{temperature:g}.txt"
            data.write_text(
                f"23N-HN {value:.8f} 0.001\n24N-HN {value + 2.0:.8f} 0.001\n",
                encoding="utf-8",
            )
            experiment = root / f"shift_{state}_{temperature:g}.toml"
            experiment.write_text(
                f"""[experiment]
name = "shift_15n_sq"
observed_state = "{state}"

[conditions]
h_larmor_frq = 800.0
temperature = {temperature}

[data]
path = "{data.name}"
""",
                encoding="utf-8",
            )
            experiments.append(experiment)
    parameters = root / "parameters.toml"
    parameters.write_text(
        """[GLOBAL]
TREF = 20.0
PB = 0.1
KEX_AB = 0.0
DW0_AB = 1.5
DW1_AB = 0.0

[CS0_A]
23N = 118.0
24N = 120.0

[CS1_A]
23N = 0.0
24N = 0.0
""",
        encoding="utf-8",
    )
    method = root / "method.toml"
    method.write_text(
        """FORMAT_VERSION = 2

[FIT]
ROLES = [
  { FIX = ["PB", "KEX_AB", "R2_A"] },
  { FIT = ["CS0_A", "CS1_A", "DW0_AB", "DW1_AB"] },
]
""",
        encoding="utf-8",
    )
    return experiments, parameters, method


def test_representative_multi_temperature_fit_recovers_planted_coefficients(
    tmp_path: Path,
) -> None:
    experiments, parameters, method = _write_shift_fit_case(tmp_path)
    output = tmp_path / "Output"
    subprocess.run(  # noqa: S603 - fixed project executable and test-owned inputs
        [
            sys.executable,
            "-m",
            "chemex",
            "fit",
            "-e",
            *(str(path) for path in experiments),
            "-p",
            str(parameters),
            "-m",
            str(method),
            "-d",
            "2st.tc",
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    with (output / "Parameters" / "fitted.toml").open("rb") as file:
        fitted = tomllib.load(file)
    with (output / "Parameters" / "fixed.toml").open("rb") as file:
        fixed = tomllib.load(file)
    assert fixed["GLOBAL"]["TREF"] == pytest.approx(20.0)
    assert fitted["CS0_A"]["23N"] == pytest.approx(118.2, abs=2.0e-4)
    assert fitted["CS1_A"]["23N"] == pytest.approx(0.015, abs=2.0e-5)
    assert fitted["DW0_AB"]["23N"] == pytest.approx(2.0, abs=2.0e-4)
    assert fitted["DW1_AB"]["23N"] == pytest.approx(0.012, abs=2.0e-5)
