"""Index, configure, and expose ChemEx analysis parameters."""

import re
import sys
from collections import Counter, defaultdict
from collections.abc import Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np

from chemex.configuration.parameters import DefaultListType
from chemex.messages import (
    print_error_constraints,
    print_error_grid_settings,
    print_warning_negative_jch,
    print_warning_positive_jnh,
)
from chemex.parameters.name import ParamName
from chemex.parameters.setting import Parameters, ParamSetting
from chemex.typing import Array

_PARAM_NAME = r"\[(.+?)\]"
_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_LINEAR = rf"lin[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_GEOMETRIC = rf"log[(](?P<start>{_FLOAT}),(?P<end>{_FLOAT}),(?P<num>\d+)[)]$"
_LIST = rf"([(](({_FLOAT})(,|[)]$))+)"
_GRID_DEFINTION = (_LINEAR, _GEOMETRIC, _LIST)


class ConstraintExpressionError(ValueError):
    def __init__(self, expression: str, detail: str) -> None:
        super().__init__(detail)
        self.expression = expression
        self.detail = detail


class GridExpressionError(ValueError):
    def __init__(self, entry: str, detail: str) -> None:
        super().__init__(detail)
        self.entry = entry
        self.detail = detail


class ModelReader(Protocol):
    @property
    def name(self) -> str: ...

    @property
    def model_free(self) -> bool: ...

    @property
    def identity(self) -> str: ...


class ParameterIndex:
    """Maintains an index of parameter names for efficient searching and retrieval.

    This class creates an index structure that allows for quick lookup of parameter
    IDs based on parameter names, which is particularly useful in large datasets
    with numerous parameters.

    Attributes:
        _index (defaultdict[Hashable, set[str]]): Maps search keys to parameter IDs.

    """

    def __init__(self) -> None:
        """Initialize an empty parameter index."""
        self._index: defaultdict[Hashable, set[str]] = defaultdict(set)

    def add(self, param_name: ParamName) -> None:
        """Add a parameter name to the index.

        Args:
            param_name (ParamName): The parameter name to be indexed.

        """
        for search_key in param_name.search_keys:
            self._index[search_key].add(param_name.id_)

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        """Retrieve parameter IDs matching a given parameter name.

        Args:
            param_name (ParamName): Parameter name to match in the index.

        Returns:
            set[str]: Set of matching parameter IDs.

        """
        search_keys = param_name.search_keys
        return set[str].intersection(
            *(self._index.get(search_key, set()) for search_key in search_keys),
        )

    def clear(self) -> None:
        """Remove all indexed parameter names."""
        self._index.clear()


def _convert_grid_expression_to_values(grid_expression: str) -> Array:
    """Convert grid expression to floating-point values.

    Parses grid expressions used for parameter exploration and converts them into
    an array of floating-point values for further processing.

    Args:
        grid_expression (str): Grid expression to be converted.

    Returns:
        Array: Array of floating-point values from the grid expression.

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
    """Manage parameter settings, expressions, and dependencies.

    Attributes:
        _parameters (Parameters): The parameters managed by the catalog.
        _index (ParameterIndex): An index for efficient parameter lookup.

    """

    _parameters: Parameters = field(default_factory=dict)
    _index: ParameterIndex = field(default_factory=ParameterIndex)

    def _add(self, parameter: ParamSetting) -> None:
        """Add a single parameter to the catalog.

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
        """Add multiple parameters to the catalog.

        Args:
            parameters (Parameters): Dictionary of parameters to be added.

        """
        for parameter in parameters.values():
            self._add(parameter)

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        """Retrieve a subset of parameters based on their IDs.

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

    def get_matching_ids(self, param_name: ParamName) -> set[str]:
        """Find parameter IDs matching a specified parameter name.

        Args:
            param_name (ParamName): Name of the parameter to match.

        Returns:
            set[str]: Set of IDs matching the parameter name.

        """
        return self._index.get_matching_ids(param_name)

    def _set_values(self, par_values: Mapping[str, float]) -> None:
        """Apply construction-time auxiliary values before configuration sealing."""
        for id_, value in par_values.items():
            if id_ in self._parameters:
                self._parameters[id_].value = value

    def set_defaults(self, defaults: DefaultListType) -> None:
        """Set default settings for a group of parameters.

        Args:
            defaults (DefaultListType): Defaults to apply.

        """
        id_pool = set(self._parameters)
        for name_to_set, setting in reversed(defaults):
            matching_ids = self.get_matching_ids(name_to_set) & id_pool
            id_pool -= matching_ids
            for matching_id in matching_ids:
                self._parameters[matching_id].set(setting)
                self._parameters[matching_id].expr = ""

    def _count_per_section(self, param_ids: set[str]) -> Counter[str]:
        """Count parameters per section for a given set of IDs.

        Args:
            param_ids (set[str]): Set of parameter IDs to count.

        Returns:
            Counter[str]: Counts of parameters per section.

        """
        return Counter(
            self._parameters[param_id].param_name.section for param_id in param_ids
        )

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        """Set the variability of parameters by section name.

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
        """Fix all parameters, preventing them from varying during fitting."""
        for parameter in self._parameters.values():
            parameter.vary = False

    def _get_ids_left(self, expression: str) -> set[str]:
        """Extract parameter IDs from the left side of an expression.

        Args:
            expression (str): The expression to parse.

        Returns:
            set[str]: Set of parameter IDs extracted.

        """
        return self.get_matching_ids(ParamName.from_section(expression))

    def _get_ids_right(self, expression: str) -> dict[str, set[str]]:
        """Extract parameter IDs from the right side of an expression.

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

    def _split_constraint_expression(self, expression: str) -> tuple[str, str]:
        separator_count = expression.count("=")
        if separator_count != 1:
            detail = "Expected exactly one '=' in the constraint expression"
            raise ConstraintExpressionError(expression, detail)

        left, _, right = expression.partition("=")
        return left, right.strip()

    def _resolve_constraint_reference(
        self,
        expression: str,
        left_id: str,
        reference: str,
        candidate_ids: set[str],
    ) -> str:
        if not candidate_ids:
            detail = f'No parameter matches reference "{reference}"'
            raise ConstraintExpressionError(expression, detail)

        non_self_candidate_ids = candidate_ids - {left_id}
        if not non_self_candidate_ids:
            detail = (
                f'Reference "{reference}" resolves only to the constrained parameter'
            )
            raise ConstraintExpressionError(expression, detail)

        param_name = self._parameters[left_id].param_name
        return param_name.get_closest_id(non_self_candidate_ids)

    def _set_expression(self, expression: str, ids_pool: set[str]) -> set[str]:
        """Set expressions for parameters based on an input expression.

        Args:
            expression (str): The expression to apply.
            ids_pool (set[str]): Set of parameter IDs to consider.

        Returns:
            set[str]: Set of parameter IDs updated.

        """
        left, right = self._split_constraint_expression(expression)

        ids_left = self._get_ids_left(left) & ids_pool
        ids_right_dict = self._get_ids_right(right)

        for id_left in ids_left:
            new_expression = right.strip()
            for section_name, ids_right in ids_right_dict.items():
                id_replace = self._resolve_constraint_reference(
                    expression,
                    id_left,
                    section_name,
                    ids_right,
                )
                new_expression = new_expression.replace(section_name, id_replace)
            self._parameters[id_left].expr = new_expression

        return ids_left

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        """Apply a list of expressions to parameters.

        Args:
            expression_list (Sequence[str]): Expressions to apply.

        Returns:
            Counter[str]: Counts of parameters updated per section.

        """
        ids_modified: set[str] = set()
        ids_pool = set(self._parameters)

        try:
            for expression in reversed(expression_list):
                ids_changed = self._set_expression(expression, ids_pool)
                ids_pool -= ids_changed
                ids_modified.update(ids_changed)
        except ConstraintExpressionError as error:
            print_error_constraints(error.expression, error.detail)
            sys.exit(1)

        return self._count_per_section(ids_modified)

    def _split_grid_entry(self, entry: str) -> tuple[str, str]:
        compact_entry = entry.replace(" ", "")

        separator_count = compact_entry.count("=")
        if separator_count != 1:
            detail = "Expected exactly one '=' in the grid entry"
            raise GridExpressionError(entry, detail)

        name, _, expression = compact_entry.partition("=")
        if not any(re.fullmatch(regex, expression) for regex in _GRID_DEFINTION):
            detail = "Unsupported grid definition"
            raise GridExpressionError(entry, detail)

        return name, expression

    def parse_grid(self, grid_entries: list[str]) -> dict[str, Array]:
        """Parse grid definitions and sets up parameters accordingly.

        Args:
            grid_entries (list[str]): List of grid definitions.

        Returns:
            dict[str, Array]: Mapping of parameter IDs to grid values.

        """
        ids_pool = set(self._parameters)

        grid_values: dict[str, Array] = {}

        try:
            for entry in reversed(grid_entries):
                name, expression = self._split_grid_entry(entry)
                ids_selected = self.get_matching_ids(ParamName.from_section(name))
                ids_changed = ids_selected & ids_pool
                values = _convert_grid_expression_to_values(expression)

                grid_values.update(dict.fromkeys(ids_changed, values))
                ids_pool -= ids_changed

        except GridExpressionError as error:
            print_error_grid_settings(error.entry, error.detail)
            sys.exit(1)

        return grid_values

    def check_params(self) -> None:
        """Check and warns about the physical validity of J couplings."""
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

    def clear(self) -> None:
        """Remove all registered parameters and rebuild the search index."""
        self._parameters.clear()
        self._index = ParameterIndex()


@dataclass
class ParameterStore:
    """Manages a collection of ParameterCatalogs for different fitting models.

    This class acts as an interface to multiple ParameterCatalogs, allowing
    for easy management of parameters across different fitting models used
    in NMR data analysis.

    Attributes:
        _database (ParameterCatalog): The primary parameter catalog.
        _database_mf (ParameterCatalog): The model-free parameter catalog.

    """

    model: ModelReader
    _database: ParameterCatalog
    _database_mf: ParameterCatalog
    _defaults_applied: bool = field(default=False, init=False, repr=False)
    _defaults: DefaultListType = field(default_factory=list, init=False, repr=False)
    _configuration_locked: bool = field(default=False, init=False, repr=False)

    @property
    def database(self) -> ParameterCatalog:
        """Returns the appropriate parameter catalog based on the model.

        Returns:
            ParameterCatalog: The active parameter catalog.

        """
        return self._database_mf if self.model.model_free else self._database

    @property
    def defaults_applied(self) -> bool:
        """Whether default and ordinary model-free initialization completed."""
        return self._defaults_applied

    def lock_configuration(self) -> None:
        """Reject later definition and configuration mutations."""
        self._configuration_locked = True

    def _ensure_configuration_open(self) -> None:
        if self._configuration_locked:
            msg = "Parameter configuration is sealed"
            raise RuntimeError(msg)

    def add_multiple(self, parameters: Parameters) -> None:
        """Add multiple parameters to the primary catalog.

        Args:
            parameters (Parameters): Parameters to be added.

        """
        self._ensure_configuration_open()
        self._database.add_multiple(parameters)
        self._defaults_applied = False

    def add_multiple_mf(self, parameters: Parameters) -> None:
        """Add multiple parameters to the model-free catalog.

        Args:
            parameters (Parameters): Parameters to be added.

        """
        self._ensure_configuration_open()
        self._database_mf.add_multiple(parameters)
        self._defaults_applied = False

    def get_parameters(self, param_ids: Iterable[str]) -> Parameters:
        """Retrieve parameters from the active catalog by IDs.

        Args:
            param_ids (Iterable[str]): IDs of the parameters to retrieve.

        Returns:
            Parameters: Retrieved parameters.

        """
        return self.database.get_parameters(param_ids)

    def sort(self) -> None:
        """Sort parameters in the active catalog."""
        self.database.sort()

    def set_defaults(self, defaults: DefaultListType) -> None:
        """Set defaults for parameters in both catalogs.

        Args:
            defaults (DefaultListType): Default settings to apply.

        """
        self._ensure_configuration_open()
        if self._defaults_applied:
            msg = "Parameter defaults have already been applied"
            raise RuntimeError(msg)
        self._defaults_applied = False
        self._defaults = list(defaults)
        self._database_mf.set_defaults(defaults)

        if self.model.model_free:
            self._defaults_applied = True
            return

        self.database.set_defaults(defaults)

        self.database.check_params()
        self._defaults_applied = True

    def seed_from_model_free_values(self, values: Mapping[str, float]) -> None:
        """Apply native model-free values before explicit ordinary defaults."""
        self._ensure_configuration_open()
        if self.model.model_free:
            return
        if not self._defaults_applied:
            msg = "Parameter defaults must be parsed before model-free initialization"
            raise RuntimeError(msg)
        self._database._set_values(values)
        self._database.set_defaults(self._defaults)
        self._database.check_params()

    def set_vary(self, section_names: Sequence[str], vary: bool) -> Counter[str]:
        """Set variability of parameters in the active catalog by section name.

        Args:
            section_names (Sequence[str]): Section names to update.
            vary (bool): Whether parameters should vary.

        Returns:
            Counter[str]: Counts of updated parameters per section.

        """
        return self.database.set_vary(section_names, vary)

    def fix_all(self) -> None:
        """Fix all parameters in the active catalog, preventing them from varying."""
        self.database.fix_all()

    def set_expressions(self, expression_list: Sequence[str]) -> Counter[str]:
        """Apply expressions to parameters in the active catalog.

        Args:
            expression_list (Sequence[str]): Expressions to apply.

        Returns:
            Counter[str]: Counts of updated parameters per section.

        """
        return self.database.set_expressions(expression_list)

    def parse_grid(self, grid_entries: list[str]) -> dict[str, Array]:
        """Parse grid definitions and sets up parameters in the active catalog.

        Args:
            grid_entries (list[str]): Grid definitions to parse.

        Returns:
            dict[str, Array]: Parameters mapped to grid values.

        """
        return self.database.parse_grid(grid_entries)

    def reset(self) -> None:
        """Clear all parameter catalogs used by the current process."""
        self._database.clear()
        self._database_mf.clear()
        self._defaults_applied = False
        self._defaults.clear()
        self._configuration_locked = False


def create_parameter_store(model_reader: ModelReader) -> ParameterStore:
    return ParameterStore(model_reader, ParameterCatalog(), ParameterCatalog())
