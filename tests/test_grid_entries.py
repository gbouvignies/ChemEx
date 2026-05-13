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
    return ParamSetting(ParamName(name, spin_system), value=1.0, vary=True)


def make_catalog(*parameters: ParamSetting) -> ParameterCatalog:
    catalog = ParameterCatalog()
    catalog.add_multiple({parameter.id_: parameter for parameter in parameters})
    return catalog


def get_setting(catalog: ParameterCatalog, param_id: str) -> ParamSetting:
    return catalog.get_parameters([param_id])[param_id]


def test_parse_grid_preserves_specific_entry_precedence() -> None:
    global_pb = make_param("PB")
    local_pb = make_param("PB", "10N-H")
    catalog = make_catalog(global_pb, local_pb)

    grid = catalog.parse_grid(
        ["PB = lin(1, 2, 2)", "PB, NUC->10N-H = lin(3, 4, 2)"],
    )

    assert grid[global_pb.id_].tolist() == [1.0, 2.0]
    assert grid[local_pb.id_].tolist() == [3.0, 4.0]
    assert get_setting(catalog, global_pb.id_).vary is False
    assert get_setting(catalog, local_pb.id_).vary is False


def test_parse_grid_last_matching_entry_wins() -> None:
    global_pb = make_param("PB")
    local_pb = make_param("PB", "10N-H")
    catalog = make_catalog(global_pb, local_pb)

    grid = catalog.parse_grid(
        ["PB, NUC->10N-H = lin(3, 4, 2)", "PB = lin(1, 2, 2)"],
    )

    assert grid[global_pb.id_].tolist() == [1.0, 2.0]
    assert grid[local_pb.id_].tolist() == [1.0, 2.0]


def test_parse_grid_exits_on_missing_separator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    catalog = make_catalog(make_param("PB"))
    recorded: dict[str, str] = {}

    monkeypatch.setattr(
        database_module,
        "print_error_grid_settings",
        lambda entry, detail=None: recorded.update(
            {"entry": entry, "detail": detail or ""}
        ),
    )

    with pytest.raises(SystemExit):
        catalog.parse_grid(["PB lin(1, 2, 2)"])

    assert recorded == {
        "entry": "PB lin(1, 2, 2)",
        "detail": "Expected exactly one '=' in the grid entry",
    }


def test_parse_grid_exits_on_invalid_definition(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    catalog = make_catalog(make_param("PB"))
    recorded: dict[str, str] = {}

    monkeypatch.setattr(
        database_module,
        "print_error_grid_settings",
        lambda entry, detail=None: recorded.update(
            {"entry": entry, "detail": detail or ""}
        ),
    )

    with pytest.raises(SystemExit):
        catalog.parse_grid(["PB = bad(1, 2, 2)"])

    assert recorded == {
        "entry": "PB = bad(1, 2, 2)",
        "detail": "Unsupported grid definition",
    }
