import chemex.parameters.helper as cph
import chemex.parameters.kinetics as cpk
import chemex.parameters.liouvillian as cpl


def create_params(basis, model, conditions, spin_system=None, constraints=None):
    fnames_k, params_k = cpk.create_params_k(model, conditions, spin_system)
    fnames_l, params_l = cpl.create_params_l(
        basis, model, conditions, spin_system, constraints
    )
    pnames = {**fnames_k, **fnames_l}
    params = cph.merge((params_k, params_l))
    return pnames, params
