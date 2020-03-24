import pathlib as pl

import numpy as np

import chemex.containers.cest as ccce
import chemex.containers.cpmg as cccp
import chemex.containers.experiment as cce
import chemex.containers.relaxation as ccr
import chemex.containers.shift as ccs
import chemex.helper as ch
import chemex.nmr.propagator as cnp
import chemex.nmr.spin_system as cns
import chemex.parameters.helper as cph


def load_experiment(config, pulse_seq_cls):
    read = experiment_cls = profile_cls = schema = None
    for key, container in _CONTAINERS.items():
        if config["experiment"]["name"].startswith(key):
            read = container["read"]
            experiment_cls = container["experiment"]
            profile_cls = container["profile"]
            schema = container["schema"]
            break
    ch.validate(config, schema)
    profiles = read(config, pulse_seq_cls, profile_cls)
    experiment = experiment_cls(config=config, profiles=profiles)
    experiment.estimate_noise(config["data"]["error"])
    experiment.merge_same_profiles()
    return experiment


def _read_profiles(config, pulse_seq_cls, profile_cls):
    propagator = cnp.PropagatorIS.from_config(config)
    paths = _get_profile_paths(config)
    profiles = []
    for path, spin_system in paths.items():
        config["spin_system"] = spin_system
        pnames, params = cph.create_params(config, propagator,)
        pulse_seq = pulse_seq_cls(config, propagator)
        profile = profile_cls.from_file(path, config, pulse_seq, pnames, params,)
        profiles.append(profile)
    return sorted(profiles)


def _read_shifts(config, pulse_seq_cls, profile_cls):
    propagator = cnp.PropagatorIS.from_config(config)
    shifts = _get_shifts(config)
    profiles = []
    for spin_system, data in shifts.items():
        config["spin_system"] = spin_system
        pnames, params = cph.create_params(config, propagator)
        pulse_seq = pulse_seq_cls(config, propagator)
        profile = profile_cls(spin_system, data, pulse_seq, pnames, params)
        profiles.append(profile)
    return sorted(profiles)


def _get_profile_paths(config):
    path = ch.normalize_path(config["filename"].parent, pl.Path(config["data"]["path"]))
    include = config["selection"]["include"]
    exclude = config["selection"]["exclude"]
    paths = {}
    for name, filename in config["data"]["profiles"]:
        spin_system = cns.SpinSystem(name)
        included = include is None or spin_system.part_of(include)
        excluded = exclude is not None and spin_system.part_of(exclude)
        if included and not excluded:
            paths[path / filename] = spin_system
    return paths


def _get_shifts(config):
    path = ch.normalize_path(config["filename"].parent, pl.Path(config["data"]["path"]))
    include = config["selection"]["include"]
    exclude = config["selection"]["exclude"]
    data = np.loadtxt(
        path / config["data"]["shifts"],
        dtype=[("name", "U15"), ("shift", "f8"), ("error", "f8")],
    )
    shifts = {}
    for name, shift, error in data:
        spin_system = cns.SpinSystem(name)
        included = include is None or spin_system.part_of(include)
        excluded = exclude is not None and spin_system.part_of(exclude)
        if included and not excluded:
            shifts[spin_system] = {"shift": shift, "error": error}
    return shifts


_CONTAINERS = {
    "relaxation": {
        "experiment": cce.RelaxationExperiment,
        "profile": ccr.RelaxationProfile,
        "read": _read_profiles,
        "schema": ccr.RELAXATION_SCHEMA,
    },
    "cest": {
        "experiment": cce.RelaxationExperiment,
        "profile": ccce.CestProfile,
        "read": _read_profiles,
        "schema": ccce.CEST_SCHEMA,
    },
    "dcest": {
        "experiment": cce.RelaxationExperiment,
        "profile": ccce.CestProfile,
        "read": _read_profiles,
        "schema": ccce.CEST_SCHEMA,
    },
    "coscest": {
        "experiment": cce.RelaxationExperiment,
        "profile": ccce.CestProfile,
        "read": _read_profiles,
        "schema": ccce.CEST_SCHEMA,
    },
    "cpmg": {
        "experiment": cce.RelaxationExperiment,
        "profile": cccp.CpmgProfile,
        "read": _read_profiles,
        "schema": cccp.CPMG_SCHEMA,
    },
    "shift": {
        "experiment": cce.ShiftExperiment,
        "profile": ccs.ShiftProfile,
        "read": _read_shifts,
        "schema": ccs.SHIFT_SCHEMA,
    },
}
