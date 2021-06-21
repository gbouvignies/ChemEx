import functools as fn
import itertools as it
import operator
import pathlib


def create_groups(experiments, params):
    """Create groups of datapoints that depend on disjoint sets of variables.

    For example, if the population of the minor state and the exchange
    rate are set to 'fix', chances are that the fit can be decomposed
    residue-specifically.
    """
    params_exp = experiments.select_params(params)
    pname_groups = group_pnames(experiments, params_exp)
    number = len(pname_groups)
    if number > 1:
        path = pathlib.Path("Groups")
        message = f"-- Group {{index}}/{number} ({{name}}) --\n"
    else:
        path = pathlib.Path("")
        message = ""
    groups_dict = {}
    for index, (gname, pnames) in enumerate(pname_groups.items(), start=1):
        exps = experiments.get_relevant_subset(pnames)
        group = {
            "experiments": exps,
            "params": exps.select_params(params_exp),
            "path": path / gname.folder,
        }
        groups_dict[gname] = group
    groups = []
    for index, (name, group) in enumerate(sorted(groups_dict.items()), start=1):
        group["message"] = message.format(index=index, name=name)
        groups.append(group)
    return groups


def group_pnames(experiments, params):

    pnames_vary = {
        param.name for param in params.values() if param.vary and not param.expr
    }

    groups = {
        frozenset(pname_set & pnames_vary): frozenset(pname_set & pnames_vary)
        for pname_set in experiments.pname_sets
    }

    grouped = True

    while grouped:
        grouped = False
        for g1, g2 in it.combinations(groups, r=2):
            intersection = g1 & g2
            if intersection:
                union = groups.pop(g1) | groups.pop(g2)
                groups.setdefault(intersection, set()).update(union)
                grouped = True
                break

    return {fnames_to_pname(union, params): union for union in groups.values()}


def fnames_to_pname(fnames, params):
    pnames = (params[fname].user_data["pname"] for fname in fnames)
    return fn.reduce(operator.and_, pnames)
