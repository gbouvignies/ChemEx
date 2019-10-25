import pathlib as pl

import chemex.containers.experiment as cce
import chemex.helper as ch
import chemex.nmr.spin_system as cns
import chemex.parameters as cp
import chemex.parameters.settings as cps


def get_type():
    print(__name__)
    return __name__.split(".")[-1]


def read(
    config,
    pulse_seq_cls,
    propagator_cls,
    container_cls,
    rates_cls=None,
    fit_setting=None,
):
    filename = config["filename"]
    exp_type = config["experiment"]["name"]
    paths = _get_profile_paths(config)
    profiles = _read_profiles(
        paths, config, pulse_seq_cls, propagator_cls, container_cls, fit_setting
    )
    h_frq = config["conditions"]["h_larmor_frq"]
    rates = None
    if rates_cls is not None:
        rates = rates_cls(h_frq, config["spin_system"].get("rates"))
    experiment = cce.RelaxationExperiment(filename, exp_type, profiles, rates)
    experiment.estimate_noise(config["data"]["error"])
    experiment.merge_same_profiles()
    return experiment


def _get_profile_paths(config):
    path = ch.normalize_path(config["filename"].parent, pl.Path(config["data"]["path"]))
    paths = {
        path / profile[1]: cns.SpinSystem(profile[0])
        for profile in config["data"]["profiles"]
    }
    selection = config["selection"]
    if selection is not None:
        for path, spin_system in paths.copy().items():
            if not any(spin_system & name == name for name in selection):
                del paths[path]
    return paths


def _read_profiles(
    paths, config, pulse_seq_cls, propagator_cls, container_cls, fit_setting
):
    model = config["model"]
    conditions = config["conditions"]
    basis = config["spin_system"]["basis"]
    atoms = config["spin_system"]["atoms"]
    constraints = config["spin_system"].get("constraints")
    h_frq = conditions["h_larmor_frq"]
    propagator = propagator_cls(basis=basis, model=model, atoms=atoms, h_frq=h_frq)
    profiles = []
    for path, spin_system in paths.items():
        config["spin_system"]["spin_system"] = spin_system
        par_names, params_default = cp.create_params(
            basis=basis,
            model=model,
            conditions=conditions,
            spin_system=spin_system,
            constraints=constraints,
        )
        pulse_seq = pulse_seq_cls(config=config, propagator=propagator)
        profile = container_cls.from_file(
            path=path,
            config=config,
            pulse_seq=pulse_seq,
            par_names=par_names,
            params_default=params_default,
        )
        cps.set_status(profile.params_default, fit_setting, verbose=False)
        profiles.append(profile)
    return sorted(profiles)
