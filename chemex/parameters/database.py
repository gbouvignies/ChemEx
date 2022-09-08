from __future__ import annotations

from collections import Counter
from collections import defaultdict
from collections.abc import Hashable
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from re import compile
from typing import DefaultDict

import numpy as np
from lmfit import Parameters as ParametersLF

from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultListType
from chemex.messages import print_status_changes
from chemex.model import model
from chemex.nmr.rates import rate_functions
from chemex.parameters.name import ParamName
from chemex.parameters.setting import Parameters
from chemex.parameters.setting import ParamSetting

_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_RE_PARAM_NAME = compile(r"\[(.+?)\]")
_RE_GRID_DEFINITION = compile(
    rf"(lin[(]{_FLOAT},{_FLOAT},\d+[)]$)|"
    rf"(log[(]{_FLOAT},{_FLOAT},\d+[)]$)|"
    rf"([(](({_FLOAT})(,|[)]$))+)"
)


class ParameterIndex:
    def __init__(self) -> None:
        self._index: DefaultDict[Hashable, set[str]] = defaultdict(set)

    def add(self, param_name: ParamName) -> None:
        for search_key in param_name.search_keys:
            self._index[search_key].add(param_name.id)

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        search_keys = param_name.search_keys
        return set.intersection(
            *(self._index.get(search_key, set()) for search_key in search_keys)
        )


def _convert_grid_expression_to_values(grid_expression: str) -> np.ndarray:
    to_be_evaluated = grid_expression.replace("lin", "np.linspace")
    to_be_evaluated = to_be_evaluated.replace("log", "np.geomspace")
    return np.asarray(eval(to_be_evaluated))


@dataclass
class ParameterCatalog:
    __parameters: Parameters = field(default_factory=dict)
    __index: ParameterIndex = ParameterIndex()

    def _add(self, parameter: ParamSetting) -> None:

        if parameter.id not in self.__parameters:
            self.__parameters[parameter.id] = parameter
            self.__index.add(parameter.param_name)

        if parameter.vary:
            self.__parameters[parameter.id].vary = True
            self.__parameters[parameter.id].expr = ""

    def add_multiple(self, parameters: Parameters) -> None:
        for parameter in parameters.values():
            self._add(parameter)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:

        relevant_ids = set()

        pool_ids = set(self.__parameters) & set(param_ids)

        while pool_ids:
            for param_id in pool_ids.copy():
                relevant_ids.add(param_id)
                pool_ids.discard(param_id)
                pool_ids.update(self.__parameters[param_id].dependencies)
            pool_ids -= relevant_ids

        return {
            param_id: parameter
            for param_id, parameter in self.__parameters.items()
            if param_id in relevant_ids
        }

    def build_lmfit_params(
        self, param_ids: Iterable[str] | None = None
    ) -> ParametersLF:

        if param_ids is None:
            param_ids = set(self.__parameters)

        parameters = self.get_parameters(param_ids)

        parameter_args = (parameter.args for parameter in parameters.values())

        lmfit_params = ParametersLF(usersyms=rate_functions)
        lmfit_params.add_many(*parameter_args)
        lmfit_params.update_constraints()

        for id, lf_param in lmfit_params.items():
            lf_param.stderr = parameters[id].stderr

        return lmfit_params

    def update_from_lmfit_params(self, parameters: ParametersLF) -> None:
        for param_id, parameter in parameters.items():
            self.__parameters[param_id].value = parameter.value
            self.__parameters[param_id].stderr = parameter.stderr

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        return self.__index.get_matching_ids(param_name)

    def get_value(self, id_: str) -> float | None:
        return self.__parameters[id_].value

    def set_values(self, par_values: dict[str, float]) -> None:
        for id_, value in par_values.items():
            if id_ in self.__parameters:
                self.__parameters[id_].value = value

    def set_defaults(self, defaults: DefaultListType) -> None:
        id_pool = set(self.__parameters)
        for name_to_set, setting in reversed(defaults):
            matching_ids = self.get_matching_ids(name_to_set) & id_pool
            id_pool -= matching_ids
            for matching_id in matching_ids:
                self.__parameters[matching_id].set(setting)

    def _count_per_section(self, param_ids: set[str]) -> Counter[str]:
        return Counter(
            self.__parameters[param_id].param_name.section for param_id in param_ids
        )

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:

        parameters = self.__parameters

        ids_modified = set()
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
        for parameter in self.__parameters.values():
            parameter.vary = False

    def _get_ids_left(self, expression: str) -> set[str]:
        return self.get_matching_ids(ParamName.from_section(expression))

    def _get_ids_right(self, expression: str) -> dict[str, set[str]]:
        ids_right = {}
        for match in _RE_PARAM_NAME.finditer(expression):
            param_name = ParamName.from_section(match.group(1))
            ids_right[match.group(0)] = self.get_matching_ids(param_name)
        return ids_right

    def _set_expression(self, expression: str, ids_pool: set[str]) -> set[str]:

        left, right, *something_else = expression.split("=")

        if something_else:
            print(
                f'\nError reading constraints:\n  -> "{expression}"\n\nProgram aborted\n'
            )
            exit()

        ids_left = self._get_ids_left(left) & ids_pool
        ids_right_dict = self._get_ids_right(right)

        for id_left in ids_left:
            new_expression = right.strip()
            for section_name, ids_right in ids_right_dict.items():
                name = self.__parameters[id_left].param_name
                id_replace = name.get_closest_id(ids_right)
                new_expression = new_expression.replace(section_name, id_replace)
            self.__parameters[id_left].expr = new_expression

        return ids_left

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:

        ids_modified = set()
        ids_pool = set(self.__parameters)

        for expression in reversed(expression_list):
            ids_changed = self._set_expression(expression, ids_pool)
            ids_pool -= ids_changed
            ids_modified.update(ids_changed)

        return self._count_per_section(ids_modified)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, np.ndarray]:

        ids_pool = set(self.__parameters)

        grid_values: dict[str, np.ndarray] = {}

        for entry in reversed(grid_entries):
            name, expression, *something_else = entry.replace(" ", "").split("=")

            if something_else or not _RE_GRID_DEFINITION.match(expression):
                print(
                    f'\nError reading grid settings:\n  -> "{entry}"\n\nProgram aborted\n'
                )
                exit()

            ids_left = self.get_matching_ids(ParamName.from_section(name))
            values = _convert_grid_expression_to_values(expression)

            grid_values |= {param_id: values for param_id in ids_left & ids_pool}

            self.set_vary([name], False)

        return grid_values

    def check_params(self):
        """Check whether the J couplings have the right sign"""
        messages = []
        for setting in self.__parameters.values():
            param_name = setting.param_name
            if not param_name.name.startswith("J_") or setting.value is None:
                continue
            atoms = {atom.name for atom in param_name.spin_system.atoms.values()}
            if atoms == {"N", "H"} and setting.value > 0.0:
                messages.append(
                    "Warning: Some 1J(NH) scalar couplings are set with positive values.\n"
                    "This can cause the TROSY and anti-TROSY components to be switched in\n"
                    "some experiments."
                )
            if atoms == {"C", "H"} and setting.value < 0.0:
                messages.append(
                    "Warning: Some 1J(CH) scalar couplings are set with negative values.\n"
                    "This can cause the TROSY and anti-TROSY components to be switched in\n"
                    "some experiments."
                )
        if messages:
            print("")
            print("\n\n".join(np.unique(messages)))

    def sort(self) -> None:
        sorted_items = sorted(
            self.__parameters.items(), key=lambda item: item[1].param_name
        )
        self.__parameters = dict(sorted_items)


@dataclass
class ParamManager:
    __database: ParameterCatalog
    __database_mf: ParameterCatalog

    @property
    def database(self) -> ParameterCatalog:
        return self.__database_mf if model.model_free else self.__database

    def add_multiple(self, parameters: Parameters) -> None:
        self.__database.add_multiple(parameters)

    def add_multiple_mf(self, parameters: Parameters) -> None:
        self.__database_mf.add_multiple(parameters)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        return self.database.get_parameters(param_ids)

    def get_value(self, param_id: str) -> float | None:
        return self.database.get_value(param_id)

    def build_lmfit_params(self, param_ids: Iterable[str]) -> ParametersLF:
        return self.database.build_lmfit_params(param_ids)

    def sort(self):
        self.database.sort()

    def update_from_parameters(self, parameters: ParametersLF) -> None:
        self.database.update_from_lmfit_params(parameters)

    def set_values(self, par_values: dict[str, float]) -> None:
        return self.database.set_values(par_values)

    def set_defaults(self, defaults: DefaultListType) -> None:

        self.__database_mf.set_defaults(defaults)

        if model.model_free:
            return

        params_mf = self.__database_mf.build_lmfit_params()
        self.__database_mf.set_values(params_mf.valuesdict())

        self.database.set_defaults(defaults)

        params = self.database.build_lmfit_params()
        self.database.set_values(params.valuesdict())

        self.database.check_params()

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        return self.database.set_vary(section_names, vary)

    def fix_all(self) -> None:
        self.database.fix_all()

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        return self.database.set_expressions(expression_list)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, np.ndarray]:
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
    expression."""

    matches_con = set_param_expressions(method.constraints)
    matches_fix = set_param_vary(method.fix, vary=False)
    matches_fit = set_param_vary(method.fit, vary=True)

    print_status_changes(matches_fit, matches_fix, matches_con)
