import lmfit as lf

import chemex.parameters.name as cpn


def settings_to_params(settings, fnames):
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
        param._delay_asteval = True
        parameter_list.append(param)
    params = lf.Parameters()
    params.add_many(*parameter_list)
    set_delay_eval(params, False)
    return params


def make_params(settings, conditions, spin_system=None):
    fnames = cpn._get_fnames(settings, conditions, spin_system)
    params = settings_to_params(settings, fnames)
    return fnames, params


def merge(params_list):
    params_merged = lf.Parameters()
    for params in params_list:
        set_delay_eval(params, delay=True)
        for name, param in params.items():
            if name in params_merged and params_merged[name].vary:
                continue
            params_merged[name] = param
    set_delay_eval(params_merged, delay=False)
    return params_merged


def set_delay_eval(params, delay=False):
    for param in params.values():
        param._delay_asteval = delay
