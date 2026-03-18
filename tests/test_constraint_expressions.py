from __future__ import annotations

import pytest

import chemex.parameters.database as database_module
from chemex.parameters.database import ParameterCatalog
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem


def make_param(name: str, spin_system_name: str = "") -> ParamSetting:
    spin_system = (
        SpinSystem.from_name(spin_system_name) if spin_system_name else SpinSystem()
    )
    return ParamSetting(ParamName(name, spin_system), value=1.0, vary=False)


def make_catalog(*parameters: ParamSetting) -> ParameterCatalog:
    catalog = ParameterCatalog()
    catalog.add_multiple({parameter.id_: parameter for parameter in parameters})
    return catalog


def get_expression(catalog: ParameterCatalog, param_id: str) -> str:
    return catalog.get_parameters([param_id])[param_id].expr


def test_constraint_reference_uses_non_self_match_when_available() -> None:
    global_pb = make_param("PB")
    local_pb = make_param("PB", "L55HD2")
    catalog = make_catalog(global_pb, local_pb)

    catalog.set_expressions(["PB, NUC->L55HD2 = [PB]"])

    assert get_expression(catalog, local_pb.id_) == global_pb.id_
    params = catalog.build_lmfit_params(model_name="2st")
    assert params[local_pb.id_].expr == global_pb.id_


def test_constraint_reference_exits_on_missing_match(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    global_pb = make_param("PB")
    catalog = make_catalog(global_pb)
    recorded: dict[str, str] = {}

    monkeypatch.setattr(
        database_module,
        "print_error_constraints",
        lambda expression, detail=None: recorded.update(
            {"expression": expression, "detail": detail or ""}
        ),
    )

    with pytest.raises(SystemExit):
        catalog.set_expressions(["PB = [BOGUS]"])

    assert recorded == {
        "expression": "PB = [BOGUS]",
        "detail": 'No parameter matches reference "[BOGUS]"',
    }
    assert get_expression(catalog, global_pb.id_) == ""


def test_constraint_reference_exits_on_missing_separator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    global_pb = make_param("PB")
    catalog = make_catalog(global_pb)
    recorded: dict[str, str] = {}

    monkeypatch.setattr(
        database_module,
        "print_error_constraints",
        lambda expression, detail=None: recorded.update(
            {"expression": expression, "detail": detail or ""}
        ),
    )

    with pytest.raises(SystemExit):
        catalog.set_expressions(["PB [BOGUS]"])

    assert recorded == {
        "expression": "PB [BOGUS]",
        "detail": "Expected exactly one '=' in the constraint expression",
    }
    assert get_expression(catalog, global_pb.id_) == ""


def test_constraint_reference_exits_on_self_only_match(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    local_pb = make_param("PB", "L55HD2")
    catalog = make_catalog(local_pb)
    recorded: dict[str, str] = {}

    monkeypatch.setattr(
        database_module,
        "print_error_constraints",
        lambda expression, detail=None: recorded.update(
            {"expression": expression, "detail": detail or ""}
        ),
    )

    with pytest.raises(SystemExit):
        catalog.set_expressions(["PB, NUC->L55HD2 = [PB]"])

    assert recorded == {
        "expression": "PB, NUC->L55HD2 = [PB]",
        "detail": 'Reference "[PB]" resolves only to the constrained parameter',
    }
    assert get_expression(catalog, local_pb.id_) == ""
