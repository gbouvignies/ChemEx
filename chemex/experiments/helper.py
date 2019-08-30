import pathlib as pl

import jsonschema as js

import chemex.containers.experiment as cce
import chemex.helper as ch
import chemex.nmr.helper as cnn
import chemex.parameters as cp


def validate(config, schema):
    try:
        Validator(schema).validate(config)
    except js.ValidationError as err:
        print("Validation error: {0}".format(err))
        raise


def read(config, pulse_seq_cls, propagator_cls, container_cls):
    filename = config["filename"]
    exp_type = config["experiment"]["name"]
    paths = _get_profile_paths(config)
    profiles = _read_profiles(
        paths, config, pulse_seq_cls, propagator_cls, container_cls
    )
    experiment = cce.RelaxationExperiment(filename, exp_type, profiles)
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


def _extend_with_default(validator_class):
    """Taken from: https://python-jsonschema.readthedocs.io/en/stable/faq/"""
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property_, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(property_, subschema["default"])
        for error in validate_properties(validator, properties, instance, schema):
            yield error

    return js.validators.extend(validator_class, {"properties": set_defaults})


Validator = _extend_with_default(js.Draft7Validator)
