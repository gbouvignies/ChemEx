import pathlib as pl

import chemex.containers.experiment as cce
import chemex.helper as ch
import chemex.nmr.helper as cnn
import chemex.parameters as cp


def read(config, pulse_seq_cls, propagator_cls, container_cls, rates_cls=None):
    filename = config["filename"]
    exp_type = config["experiment"]["name"]
    paths = _get_profile_paths(config)
    profiles = _read_profiles(
        paths, config, pulse_seq_cls, propagator_cls, container_cls
    )
    h_frq = config["conditions"]["h_larmor_frq"]
    rates = None
    if rates_cls is not None:
        rates = rates_cls(h_frq, config["spin_system"].get("rates"))
    experiment = cce.RelaxationExperiment(filename, exp_type, profiles, rates)
    experiment.estimate_noise(config["data"]["error"])
    return experiment


def _get_profile_paths(config):
    path = ch.normalize_path(config["filename"].parent, pl.Path(config["data"]["path"]))
    paths = {
        cnn.SpinSystem(profile[0]): path / profile[1]
        for profile in config["data"]["profiles"]
    }
    return paths


def _read_profiles(paths, config, pulse_seq_cls, propagator_cls, container_cls):
    model = config["model"]
    conditions = config["conditions"]
    basis = config["spin_system"]["basis"]
    atoms = config["spin_system"]["atoms"]
    constraints = config["spin_system"].get("constraints")
    h_frq = conditions["h_larmor_frq"]
    profiles = {}
    for spin_system, path in paths.items():
        config["spin_system"]["spin_system"] = spin_system
        par_names, params_default = cp.create_params(
            basis=basis,
            model=model,
            conditions=conditions,
            spin_system=spin_system,
            constraints=constraints,
        )
        propagator = propagator_cls(basis=basis, model=model, atoms=atoms, h_frq=h_frq)
        pulse_seq = pulse_seq_cls(config=config, propagator=propagator)
        profiles[path] = container_cls.from_file(
            path=path,
            config=config,
            pulse_seq=pulse_seq,
            par_names=par_names,
            params_default=params_default,
        )
    return profiles
