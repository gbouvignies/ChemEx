import lmfit as lf

import chemex.parameters.name as cpn


def make_params(settings, conditions, spin_system=None):
    fnames = cpn.get_fnames(settings, conditions, spin_system)
    params = _settings_to_params(settings, fnames)
    return fnames, params


def _settings_to_params(settings, fnames):
    parameter_list = []
    for name, setting in settings.items():
        fname = fnames[name]
        param = lf.Parameter(
            name=fname,
            value=setting.get("value"),
            min=setting.get("min"),
            max=setting.get("max"),
            vary=setting.get("vary"),
            expr=setting.get("expr", "").format_map(fnames),
        )
        parameter_list.append(param)
    params = lf.Parameters()
    params.add_many(*parameter_list)
    return params


def merge(params_list):
    # The symtable fix will be unnecessary in the next version of lmfit
    params_merged = {}
    symtable_merged = {}
    for params in params_list:
        unique_symbols = {
            key: params._asteval.symtable[key]
            for key in params._asteval.user_defined_symbols()
        }
        symtable_merged.update(unique_symbols)
        for name, param in params.items():
            if name in params_merged and params_merged[name].vary:
                continue
            params_merged[name] = param
    params = lf.Parameters(usersyms=symtable_merged)
    params.add_many(*params_merged.values())
    return params
