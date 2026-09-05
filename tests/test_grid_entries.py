import pytest

from chemex.parameters.database import GridExpressionError, ParameterCatalog
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


def test_parse_grid_preserves_precedence_without_mutating_parameter_roles() -> None:
    global_pb = make_param("PB")
    local_pb = make_param("PB", "10N-H")
    catalog = make_catalog(global_pb, local_pb)

    grid = catalog.parse_grid(
        ["PB = lin(1, 2, 2)", "PB, NUC->10N-H = lin(3, 4, 2)"],
    )

    assert grid[global_pb.id_].tolist() == [1.0, 2.0]
    assert grid[local_pb.id_].tolist() == [3.0, 4.0]
    assert get_setting(catalog, global_pb.id_).vary is True
    assert get_setting(catalog, local_pb.id_).vary is True


def test_parse_grid_last_matching_entry_wins() -> None:
    global_pb = make_param("PB")
    local_pb = make_param("PB", "10N-H")
    catalog = make_catalog(global_pb, local_pb)

    grid = catalog.parse_grid(
        ["PB, NUC->10N-H = lin(3, 4, 2)", "PB = lin(1, 2, 2)"],
    )

    assert grid[global_pb.id_].tolist() == [1.0, 2.0]
    assert grid[local_pb.id_].tolist() == [1.0, 2.0]


def test_parse_grid_raises_typed_error_on_missing_separator() -> None:
    catalog = make_catalog(make_param("PB"))

    with pytest.raises(GridExpressionError) as error_info:
        catalog.parse_grid(["PB lin(1, 2, 2)"])

    assert error_info.value.entry == "PB lin(1, 2, 2)"
    assert error_info.value.detail == "Expected exactly one '=' in the grid entry"


def test_parse_grid_raises_typed_error_on_invalid_definition() -> None:
    catalog = make_catalog(make_param("PB"))

    with pytest.raises(GridExpressionError) as error_info:
        catalog.parse_grid(["PB = bad(1, 2, 2)"])

    assert error_info.value.entry == "PB = bad(1, 2, 2)"
    assert error_info.value.detail == "Unsupported grid definition"
