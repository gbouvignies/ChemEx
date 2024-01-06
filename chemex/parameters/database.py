"""This module offers a framework for managing NMR analysis parameters.

It includes classes and methods for indexing, updating, and retrieving parameters,
handling parameter constraints, and integration with the lmfit library for
optimization. Additionally, it supports grid definition for parameter exploration and
ensures physical validity of J couplings.
"""
from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
from lmfit import Parameters as ParametersLF

from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultListType
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
from chemex.parameters.setting import Parameters, ParamSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

_PARAM_NAME = r"\[(.+?)\]"
_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_LINEAR = rf"lin[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_GEOMETRIC = rf"log[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_LIST = rf"([(](({_FLOAT})(,|[)]$))+)"
_GRID_DEFINTION = (_LINEAR, _GEOMETRIC, _LIST)


class ParameterIndex:
    """Maintains an index of parameter names for efficient searching and retrieval.

    This class creates an index structure that allows for quick lookup of parameter
    IDs based on parameter names, which is particularly useful in large datasets
    with numerous parameters.

    Attributes:
        _index (defaultdict[Hashable, set[str]]): Maps search keys to parameter IDs.
    """

    def __init__(self) -> None:
        """Initializes an empty parameter index."""
        self._index: defaultdict[Hashable, set[str]] = defaultdict(set)

    def add(self, param_name: ParamName) -> None:
        """Adds a parameter name to the index.

        Args:
            param_name (ParamName): The parameter name to be indexed.
        """
        for search_key in param_name.search_keys:
            self._index[search_key].add(param_name.id_)

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        """Retrieves parameter IDs matching a given parameter name.

        Args:
            param_name (ParamName): Parameter name to match in the index.

        Returns:
            set[str]: Set of matching parameter IDs.
        """
        search_keys = param_name.search_keys
        return set[str].intersection(
            *(self._index.get(search_key, set()) for search_key in search_keys),
        )


def _convert_grid_expression_to_values(grid_expression: str) -> ArrayFloat:
    """Converts grid expression to floating-point values.

    Parses grid expressions used for parameter exploration and converts them into
    an array of floating-point values for further processing.

    Args:
        grid_expression (str): Grid expression to be converted.

    Returns:
        ArrayFloat: Array of floating-point values from the grid expression.
    """
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
    """Manages a collection of parameter settings for NMR data analysis.

    This class is responsible for storing, updating, and retrieving various
    parameter settings, including their values, expressions, and dependencies.
    It also interfaces with the lmfit library for optimization purposes.

    Attributes:
        _parameters (Parameters): The parameters managed by the catalog.
        _index (ParameterIndex): An index for efficient parameter lookup.
    """

    _parameters: Parameters = field(default_factory=dict)
    _index: ParameterIndex = field(default_factory=ParameterIndex)

    def _add(self, parameter: ParamSetting) -> None:
        """Adds a single parameter to the catalog.

        Args:
            parameter (ParamSetting): The parameter to be added.
        """
        if parameter.id_ not in self._parameters:
            self._parameters[parameter.id_] = parameter
            self._index.add(parameter.param_name)

        if parameter.vary:
            self._parameters[parameter.id_].vary = True
            self._parameters[parameter.id_].expr = ""

    def add_multiple(self, parameters: Parameters) -> None:
        """Adds multiple parameters to the catalog.

        Args:
            parameters (Parameters): Dictionary of parameters to be added.
        """
        for parameter in parameters.values():
            self._add(parameter)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        """Retrieves a subset of parameters based on their IDs.

        Args:
            param_ids (Iterable[str]): IDs of the parameters to retrieve.

        Returns:
            Parameters: Subset of parameters matching the provided IDs.
        """
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
        self,
        param_ids: Iterable[str] | None = None,
    ) -> ParametersLF:
        """Constructs lmfit Parameters from the catalog's parameters.

        This method converts the internal parameter representation to the format
        required by the lmfit library for optimization.

        Args:
            param_ids (Iterable[str] | None): Specific parameter IDs to include.
                                              If None, includes all parameters.

        Returns:
            ParametersLF: lmfit-compatible parameter collection.
        """
        if param_ids is None:
            param_ids = set(self._parameters)

        parameters = self.get_parameters(param_ids)

        parameter_args = (parameter.args for parameter in parameters.values())

        usersyms = rate_functions | user_function_registry.get(model.name)
        lmfit_params = ParametersLF(usersyms=usersyms)
        lmfit_params.add_many(*parameter_args)
        lmfit_params.update_constraints()

        for id_, lf_param in lmfit_params.items():
            lf_param.stderr = parameters[id_].stderr

        return lmfit_params

    def update_from_lmfit_params(self, parameters: ParametersLF) -> None:
        """Updates catalog parameters from lmfit Parameters.

        Args:
            parameters (ParametersLF): lmfit Parameters to update from.
        """
        for param_id, parameter in parameters.items():
            self._parameters[param_id].value = parameter.value
            self._parameters[param_id].stderr = parameter.stderr

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        """Finds parameter IDs matching a specified parameter name.

        Args:
            param_name (ParamName): Name of the parameter to match.

        Returns:
            set[str]: Set of IDs matching the parameter name.
        """
        return self._index.get_matching_ids(param_name)

    def get_value(self, id_: str) -> float | None:
        """Retrieves the value of a parameter by its ID.

        Args:
            id_ (str): The ID of the parameter.

        Returns:
            float | None: The value of the parameter, or None if not found.
        """
        return self._parameters[id_].value

    def set_values(self, par_values: dict[str, float]) -> None:
        """Sets the values of specified parameters.

        Args:
            par_values (dict[str, float]): Dictionary of parameter IDs and values.
        """
        for id_, value in par_values.items():
            if id_ in self._parameters:
                self._parameters[id_].value = value

    def set_defaults(self, defaults: DefaultListType) -> None:
        """Sets default settings for a group of parameters.

        Args:
            defaults (DefaultListType): Defaults to apply.
        """
        id_pool = set(self._parameters)
        for name_to_set, setting in reversed(defaults):
            matching_ids = self.get_matching_ids(name_to_set) & id_pool
            id_pool -= matching_ids
            for matching_id in matching_ids:
                self._parameters[matching_id].set(setting)

    def _count_per_section(self, param_ids: set[str]) -> Counter[str]:
        """Counts parameters per section for a given set of IDs.

        Args:
            param_ids (set[str]): Set of parameter IDs to count.

        Returns:
            Counter[str]: Counts of parameters per section.
        """
        return Counter(
            self._parameters[param_id].param_name.section for param_id in param_ids
        )

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        """Sets the variability of parameters by section name.

        Args:
            section_names (Sequence[str]): Names of the sections to update.
            vary (bool): Whether to set parameters as variable.

        Returns:
            Counter[str]: Counts of parameters updated per section.
        """
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
        """Fixes all parameters, preventing them from varying during fitting."""
        for parameter in self._parameters.values():
            parameter.vary = False
            parameter.expr = ""

    def _get_ids_left(self, expression: str) -> set[str]:
        """Extracts parameter IDs from the left side of an expression.

        Args:
            expression (str): The expression to parse.

        Returns:
            set[str]: Set of parameter IDs extracted.
        """
        return self.get_matching_ids(ParamName.from_section(expression))

    def _get_ids_right(self, expression: str) -> dict[str, set[str]]:
        """Extracts parameter IDs from the right side of an expression.

        Args:
            expression (str): The expression to parse.

        Returns:
            dict[str, set[str]]: Mapping of section names to parameter IDs.
        """
        ids_right: dict[str, set[str]] = {}
        for match in re.finditer(_PARAM_NAME, expression):
            param_name = ParamName.from_section(match.group(1))
            ids_right[match.group(0)] = self.get_matching_ids(param_name)
        return ids_right

    def _set_expression(self, expression: str, ids_pool: set[str]) -> set[str]:
        """Sets expressions for parameters based on an input expression.

        Args:
            expression (str): The expression to apply.
            ids_pool (set[str]): Set of parameter IDs to consider.

        Returns:
            set[str]: Set of parameter IDs updated.
        """
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
        """Applies a list of expressions to parameters.

        Args:
            expression_list (Sequence[str]): Expressions to apply.

        Returns:
            Counter[str]: Counts of parameters updated per section.
        """
        ids_modified: set[str] = set()
        ids_pool = set(self._parameters)

        for expression in reversed(expression_list):
            ids_changed = self._set_expression(expression, ids_pool)
            ids_pool -= ids_changed
            ids_modified.update(ids_changed)

        return self._count_per_section(ids_modified)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, ArrayFloat]:
        """Parses grid definitions and sets up parameters accordingly.

        Args:
            grid_entries (list[str]): List of grid definitions.

        Returns:
            dict[str, ArrayFloat]: Mapping of parameter IDs to grid values.
        """
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

            self.set_vary([name], vary=False)

        return grid_values

    def check_params(self) -> None:
        """Checks and warns about the physical validity of J couplings."""
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
        """Sorts the parameters in the catalog by their names."""
        sorted_items = sorted(
            self._parameters.items(),
            key=lambda item: item[1].param_name,
        )
        self._parameters = dict(sorted_items)


@dataclass
class ParamManager:
    """Manages a collection of ParameterCatalogs for different fitting models.

    This class acts as an interface to multiple ParameterCatalogs, allowing
    for easy management of parameters across different fitting models used
    in NMR data analysis.

    Attributes:
        _database (ParameterCatalog): The primary parameter catalog.
        _database_mf (ParameterCatalog): The model-free parameter catalog.
    """

    _database: ParameterCatalog
    _database_mf: ParameterCatalog

    @property
    def database(self) -> ParameterCatalog:
        """Returns the appropriate parameter catalog based on the model.

        Returns:
            ParameterCatalog: The active parameter catalog.
        """
        return self._database_mf if model.model_free else self._database

    def add_multiple(self, parameters: Parameters) -> None:
        """Adds multiple parameters to the primary catalog.

        Args:
            parameters (Parameters): Parameters to be added.
        """
        self._database.add_multiple(parameters)

    def add_multiple_mf(self, parameters: Parameters) -> None:
        """Adds multiple parameters to the model-free catalog.

        Args:
            parameters (Parameters): Parameters to be added.
        """
        self._database_mf.add_multiple(parameters)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        """Retrieves parameters from the active catalog by IDs.

        Args:
            param_ids (Iterable[str]): IDs of the parameters to retrieve.

        Returns:
            Parameters: Retrieved parameters.
        """
        return self.database.get_parameters(param_ids)

    def get_value(self, param_id: str) -> float | None:
        """Gets the value of a parameter by its ID from the active catalog.

        Args:
            param_id (str): ID of the parameter.

        Returns:
            float | None: The value of the parameter or None if not found.
        """
        return self.database.get_value(param_id)

    def build_lmfit_params(self, param_ids: Iterable[str]) -> ParametersLF:
        """Builds lmfit Parameters from the active catalog.

        Args:
            param_ids (Iterable[str]): IDs of parameters to include.

        Returns:
            ParametersLF: lmfit-compatible parameters.
        """
        return self.database.build_lmfit_params(param_ids)

    def sort(self) -> None:
        """Sorts parameters in the active catalog."""
        self.database.sort()

    def update_from_parameters(self, parameters: ParametersLF) -> None:
        """Updates the active catalog from lmfit Parameters.

        Args:
            parameters (ParametersLF): lmfit Parameters to update from.
        """
        self.database.update_from_lmfit_params(parameters)

    def set_values(self, par_values: dict[str, float]) -> None:
        """Sets the values of specific parameters in the active catalog.

        Args:
            par_values (dict[str, float]): Parameter values to set.
        """
        return self.database.set_values(par_values)

    def set_defaults(self, defaults: DefaultListType) -> None:
        """Sets defaults for parameters in both catalogs.

        Args:
            defaults (DefaultListType): Default settings to apply.
        """
        self._database_mf.set_defaults(defaults)

        if model.model_free:
            return

        params_mf = self._database_mf.build_lmfit_params()
        self.database.set_values(params_mf.valuesdict())

        self.database.set_defaults(defaults)

        self.database.check_params()

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        """Sets variability of parameters in the active catalog by section name.

        Args:
            section_names (Sequence[str]): Section names to update.
            vary (bool): Whether parameters should vary.

        Returns:
            Counter[str]: Counts of updated parameters per section.
        """
        return self.database.set_vary(section_names, vary)

    def fix_all(self) -> None:
        """Fixes all parameters in the active catalog, preventing them from varying."""
        self.database.fix_all()

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        """Applies expressions to parameters in the active catalog.

        Args:
            expression_list (Sequence[str]): Expressions to apply.

        Returns:
            Counter[str]: Counts of updated parameters per section.
        """
        return self.database.set_expressions(expression_list)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, ArrayFloat]:
        """Parses grid definitions and sets up parameters in the active catalog.

        Args:
            grid_entries (list[str]): Grid definitions to parse.

        Returns:
            dict[str, ArrayFloat]: Parameters mapped to grid values.
        """
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


def set_parameter_status(method: Method) -> None:
    """Sets parameter status based on a given method configuration.

    This function configures parameters to vary, be fixed, or follow certain
    expressions based on the provided method configuration.

    Args:
        method (Method): Method configuration to apply.
    """
    matches_con = set_param_expressions(method.constraints)
    matches_fix = set_param_vary(method.fix, vary=False)
    matches_fit = set_param_vary(method.fit, vary=True)

    print_status_changes(matches_fit, matches_fix, matches_con)
