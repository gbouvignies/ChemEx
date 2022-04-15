from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from operator import and_
from pathlib import Path

from chemex.containers.experiments import Experiments
from chemex.parameters import database
from chemex.parameters.name import ParamName


@dataclass
class Group:
    path: Path
    message: str
    experiments: Experiments


def ids_to_param_name(param_ids: set[str]) -> ParamName:
    parameters = database.get_parameters(param_ids)
    param_names = (parameter.param_name for parameter in parameters.values())
    return reduce(and_, param_names)


def group_ids(experiments: Experiments) -> list[tuple[ParamName, set[str]]]:

    parameters = database.get_parameters(experiments.param_ids)

    ids_vary: set[str] = {
        param_id
        for param_id, param in parameters.items()
        if param.vary and not param.expr
    }

    param_id_sets = (
        set(database.get_parameters(param_id_set))
        for param_id_set in experiments.param_id_sets
    )

    groups: dict[frozenset[str], set[str]] = {
        frozenset(param_id_set & ids_vary): set(param_id_set) & ids_vary
        for param_id_set in param_id_sets
    }

    grouped = True

    while grouped:
        grouped = False
        for (i1, u1), (i2, u2) in combinations(groups.items(), r=2):
            if intersection := frozenset(u1 & u2):
                union = groups.pop(i1) | groups.pop(i2)
                groups.setdefault(intersection, set()).update(union)
                grouped = True
                break

    return sorted((ids_to_param_name(union), union) for union in groups.values())


def create_groups(experiments: Experiments) -> list[Group]:
    """Create groups of datapoints that depend on disjoint sets of variables.

    For example, if the population of the minor state and the exchange
    rate are set to 'fix', chances are that the fit can be decomposed
    residue-specifically.

    """
    id_groups = group_ids(experiments)

    group_name_template = ""
    path_name_template = ""
    group_message_template = ""

    if len(id_groups) > 1:
        group_name_template = "{index}_{pname.folder}"
        path_name_template += "Groups/{group_name}"
        group_message_template = (
            f"Group [magenta]{{group_name}}[/] ({{index}}/{len(id_groups)})"
        )

    groups = []
    for index, (pname, param_ids) in enumerate(id_groups, start=1):
        group_name = group_name_template.format(index=index, pname=pname)
        group_path = Path(path_name_template.format(group_name=group_name))
        param_ids = set(database.get_parameters(param_ids))
        group_experiments = experiments.get_relevant_subset(param_ids)
        group_message = group_message_template.format(
            index=index, group_name=group_name
        )
        groups.append(Group(group_path, group_message, group_experiments))

    return groups
