import pathlib

import chemex.parameters.name as cpn


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
    for index, pnames in enumerate(pname_groups, start=1):
        exps = experiments.get_relevant_subset(pnames)
        name = exps.get_group_name() if number > 1 else cpn.ParamName()
        group = {
            "experiments": exps,
            "params": exps.select_params(params_exp),
            "path": path / name.folder,
        }
        groups_dict[name] = group
    groups = []
    for index, (name, group) in enumerate(sorted(groups_dict.items()), start=1):
        group["message"] = message.format(index=index, name=name)
        groups.append(group)
    return groups


def group_pnames(experiments, params):
    pnames_vary = {
        name for name, param in params.items() if param.vary and not param.expr
    }
    groups = []
    for pnames in experiments.pname_sets:
        varies = pnames & pnames_vary
        found = False
        for group in groups:
            if varies & group:
                group |= varies
                found = True
                break
        if not found and varies:
            groups.append(varies)
    return groups
