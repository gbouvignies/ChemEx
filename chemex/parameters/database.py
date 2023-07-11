from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from lmfit import Parameters as ParametersLF

from chemex.messages import (
    print_error_constraints,
    print_error_grid_settings,
    print_status_changes,
    print_warning_negative_jch,
    print_warning_positive_jnh,
)
from chemex.models.model import model
from chemex.nmr.rates import rate_functions
from chemex.parameters.name import ParamName
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable, Sequence

    from chemex.configuration.methods import Method
    from chemex.configuration.parameters import DefaultListType
    from chemex.parameters.setting import Parameters, ParamSetting
    from chemex.typing import ArrayFloat


_PARAM_NAME = r"\[(.+?)\]"
_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_LINEAR = rf"lin[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_GEOMETRIC = rf"log[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_LIST = rf"([(](({_FLOAT})(,|[)]$))+)"
_GRID_DEFINTION = (_LINEAR, _GEOMETRIC, _LIST)


class ParameterIndex:
    def __init__(self) -> None:
        self._index: defaultdict[Hashable, set[str]] = defaultdict(set)

    def add(self, param_name: ParamName) -> None:
        for search_key in param_name.search_keys:
            self._index[search_key].add(param_name.id)

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        search_keys = param_name.search_keys
        return set[str].intersection(
            *(self._index.get(search_key, set()) for search_key in search_keys)
        )


def _convert_grid_expression_to_values(grid_expression: str) -> ArrayFloat:
    if match := re.match(_LINEAR, grid_expression):
        return np.linspace(
            float(match.group("start")),
            float(match.group("end")),
            int(match.group("num")),
        )
    if match := re.match(_GEOMETRIC, grid_expression):
        return np.geomspace(
            float(match.group("start")),
            float(match.group("end")),
            int(match.group("num")),
        )
    return np.fromstring(grid_expression.strip("(){}[]"), sep=",")


@dataclass
class ParameterCatalog:
    _parameters: Parameters = field(default_factory=dict)
    _index: ParameterIndex = field(default_factory=ParameterIndex)

    def _add(self, parameter: ParamSetting) -> None:
        if parameter.id not in self._parameters:
            self._parameters[parameter.id] = parameter
            self._index.add(parameter.param_name)

        if parameter.vary:
            self._parameters[parameter.id].vary = True
            self._parameters[parameter.id].expr = ""

    def add_multiple(self, parameters: Parameters) -> None:
        for parameter in parameters.values():
            self._add(parameter)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        relevant_ids: set[str] = set()

        pool_ids = set(self._parameters) & set(param_ids)

        while pool_ids:
            for param_id in pool_ids.copy():
                relevant_ids.add(param_id)
                pool_ids.discard(param_id)
                pool_ids.update(self._parameters[param_id].dependencies)
            pool_ids -= relevant_ids

        return {
            param_id: parameter
            for param_id, parameter in self._parameters.items()
            if param_id in relevant_ids
        }

    def build_lmfit_params(
        self, param_ids: Iterable[str] | None = None
    ) -> ParametersLF:
        if param_ids is None:
            param_ids = set(self._parameters)

        parameters = self.get_parameters(param_ids)

        parameter_args = (parameter.args for parameter in parameters.values())

        usersyms = rate_functions | user_function_registry.get(model.name)
        lmfit_params = ParametersLF(usersyms=usersyms)
        lmfit_params.add_many(*parameter_args)
        lmfit_params.update_constraints()

        for id, lf_param in lmfit_params.items():
            lf_param.stderr = parameters[id].stderr

        return lmfit_params

    def update_from_lmfit_params(self, parameters: ParametersLF) -> None:
        for param_id, parameter in parameters.items():
            self._parameters[param_id].value = parameter.value
            self._parameters[param_id].stderr = parameter.stderr

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        return self._index.get_matching_ids(param_name)

    def get_value(self, id_: str) -> float | None:
        return self._parameters[id_].value

    def set_values(self, par_values: dict[str, float]) -> None:
        for id_, value in par_values.items():
            if id_ in self._parameters:
                self._parameters[id_].value = value

    def set_defaults(self, defaults: DefaultListType) -> None:
        id_pool = set(self._parameters)
        for name_to_set, setting in reversed(defaults):
            matching_ids = self.get_matching_ids(name_to_set) & id_pool
            id_pool -= matching_ids
            for matching_id in matching_ids:
                self._parameters[matching_id].set(setting)

    def _count_per_section(self, param_ids: set[str]) -> Counter[str]:
        return Counter(
            self._parameters[param_id].param_name.section for param_id in param_ids
        )

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        parameters = self._parameters

        ids_modified: set[str] = set()
        ids_pool = {
            param_id
            for param_id, setting in parameters.items()
            if setting.vary != vary or setting.expr
        }

        for section_name in reversed(section_names):
            param_name = ParamName.from_section(section_name)
            matching_ids = self.get_matching_ids(param_name) & ids_pool
            for matching_id in matching_ids:
                parameters[matching_id].vary = vary
                if vary:
                    parameters[matching_id].expr = ""
            ids_modified.update(matching_ids)

        return self._count_per_section(ids_modified)

    def fix_all(self) -> None:
        for parameter in self._parameters.values():
            parameter.vary = False
            parameter.expr = ""

    def _get_ids_left(self, expression: str) -> set[str]:
        return self.get_matching_ids(ParamName.from_section(expression))

    def _get_ids_right(self, expression: str) -> dict[str, set[str]]:
        ids_right: dict[str, set[str]] = {}
        for match in re.finditer(_PARAM_NAME, expression):
            param_name = ParamName.from_section(match.group(1))
            ids_right[match.group(0)] = self.get_matching_ids(param_name)
        return ids_right

    def _set_expression(self, expression: str, ids_pool: set[str]) -> set[str]:
        left, right, *something_else = expression.split("=")

        if something_else:
            print_error_constraints(expression)
            sys.exit()

        ids_left = self._get_ids_left(left) & ids_pool
        ids_right_dict = self._get_ids_right(right)

        for id_left in ids_left:
            new_expression = right.strip()
            for section_name, ids_right in ids_right_dict.items():
                name = self._parameters[id_left].param_name
                id_replace = name.get_closest_id(ids_right)
                new_expression = new_expression.replace(section_name, id_replace)
            self._parameters[id_left].expr = new_expression

        return ids_left

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        ids_modified: set[str] = set()
        ids_pool = set(self._parameters)

        for expression in reversed(expression_list):
            ids_changed = self._set_expression(expression, ids_pool)
            ids_pool -= ids_changed
            ids_modified.update(ids_changed)

        return self._count_per_section(ids_modified)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, ArrayFloat]:
        ids_pool = set(self._parameters)

        grid_values: dict[str, ArrayFloat] = {}

        for entry in reversed(grid_entries):
            name, expression, *something_else = entry.replace(" ", "").split("=")

            valid_definition = any(
                re.match(regex, expression) for regex in _GRID_DEFINTION
            )

            if something_else or not valid_definition:
                print_error_grid_settings(entry)
                sys.exit()

            ids_left = self.get_matching_ids(ParamName.from_section(name))
            values = _convert_grid_expression_to_values(expression)

            grid_values |= {param_id: values for param_id in ids_left & ids_pool}

            self.set_vary([name], False)

        return grid_values

    def check_params(self):
        """Check whether the J couplings have the right sign."""
        positive_jnh = False
        negative_jch = False
        for setting in self._parameters.values():
            param_name = setting.param_name
            if not param_name.name.startswith("J_") or setting.value is None:
                continue
            atoms = {atom.name for atom in param_name.spin_system.atoms.values()}
            if atoms == {"N", "H"} and setting.value > 0:
                positive_jnh = True
            if atoms == {"C", "H"} and setting.value < 0:
                negative_jch = True
        if positive_jnh:
            print_warning_positive_jnh()
        if negative_jch:
            print_warning_negative_jch()

    def sort(self) -> None:
        sorted_items = sorted(
            self._parameters.items(), key=lambda item: item[1].param_name
        )
        self._parameters = dict(sorted_items)


@dataclass
class ParamManager:
    _database: ParameterCatalog
    _database_mf: ParameterCatalog

    @property
    def database(self) -> ParameterCatalog:
        return self._database_mf if model.model_free else self._database

    def add_multiple(self, parameters: Parameters) -> None:
        self._database.add_multiple(parameters)

    def add_multiple_mf(self, parameters: Parameters) -> None:
        self._database_mf.add_multiple(parameters)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        return self.database.get_parameters(param_ids)

    def get_value(self, param_id: str) -> float | None:
        return self.database.get_value(param_id)

    def build_lmfit_params(self, param_ids: Iterable[str]) -> ParametersLF:
        return self.database.build_lmfit_params(param_ids)

    def sort(self) -> None:
        self.database.sort()

    def update_from_parameters(self, parameters: ParametersLF) -> None:
        self.database.update_from_lmfit_params(parameters)

    def set_values(self, par_values: dict[str, float]) -> None:
        return self.database.set_values(par_values)

    def set_defaults(self, defaults: DefaultListType) -> None:
        self._database_mf.set_defaults(defaults)

        if model.model_free:
            return

        params_mf = self._database_mf.build_lmfit_params()
        self.database.set_values(params_mf.valuesdict())

        self.database.set_defaults(defaults)

        self.database.check_params()

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        return self.database.set_vary(section_names, vary)

    def fix_all(self) -> None:
        self.database.fix_all()

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        return self.database.set_expressions(expression_list)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, ArrayFloat]:
        return self.database.parse_grid(grid_entries)


_parameter_catalog = ParameterCatalog()
_parameter_catalog_mf = ParameterCatalog()
_manager = ParamManager(_parameter_catalog, _parameter_catalog_mf)

set_param_vary = _manager.set_vary
set_param_expressions = _manager.set_expressions
add_parameters = _manager.add_multiple
add_parameters_mf = _manager.add_multiple_mf
get_parameters = _manager.get_parameters
build_lmfit_params = _manager.build_lmfit_params
update_from_parameters = _manager.update_from_parameters
parse_grid = _manager.parse_grid
set_param_values = _manager.set_values
set_param_defaults = _manager.set_defaults
sort_parameters = _manager.sort
fix_all_parameters = _manager.fix_all


def set_parameter_status(method: Method):
    """Set whether or not to vary a fitting parameter or to use a mathematical
    expression.
    """
    matches_con = set_param_expressions(method.constraints)
    matches_fix = set_param_vary(method.fix, vary=False)
    matches_fit = set_param_vary(method.fit, vary=True)

    print_status_changes(matches_fit, matches_fix, matches_con)
