import pathlib as pl

import numpy as np

import chemex.containers.experiment as cce
import chemex.helper as ch
import chemex.nmr.spin_system as cns
import chemex.parameters as cp
import chemex.parameters.settings as cps


def read(
    config,
    pulse_seq_cls,
    propagator_cls,
    container_cls,
    rates_cls=None,
    fit_setting=None,
):
    profiles = _read_profiles(
        config, pulse_seq_cls, propagator_cls, container_cls, fit_setting
    )
    rates = None
    if rates_cls is not None:
        rates = rates_cls(
            h_larmor_frq=config["conditions"]["h_larmor_frq"],
            spins=config["spin_system"].get("rates"),
        )
    experiment = cce.RelaxationExperiment(
        filename=config["filename"],
        name=config["experiment"]["name"],
        profiles=profiles,
        rates=rates,
    )
    experiment.estimate_noise(config["data"]["error"])
    experiment.merge_same_profiles()
    return experiment


def read_shift(
    config,
    pulse_seq_cls,
    propagator_cls,
    container_cls,
    rates_cls=None,
    fit_setting=None,
):
    profiles = _read_shifts(
        config, pulse_seq_cls, propagator_cls, container_cls, fit_setting
    )
    rates = None
    if rates_cls is not None:
        rates = rates_cls(
            h_larmor_frq=config["conditions"]["h_larmor_frq"],
            spins=config["spin_system"].get("rates"),
        )
    experiment = cce.ShiftExperiment(
        filename=config["filename"],
        name=config["experiment"]["name"],
        profiles=profiles,
        rates=rates,
    )
    experiment.estimate_noise(config["data"]["error"])
    experiment.merge_same_profiles()
    return experiment


def _read_profiles(config, pulse_seq_cls, propagator_cls, container_cls, fit_setting):
    propagator = propagator_cls.from_config(config)
    paths = _get_profile_paths(config)
    profiles = []
    for path, spin_system in paths.items():
        config["spin_system"]["spin_system"] = spin_system
        pnames, params_default = cp.create_params(
            basis=config["spin_system"]["basis"],
            model=config["model"],
            conditions=config["conditions"],
            spin_system=spin_system,
            constraints=config["spin_system"].get("constraints"),
        )
        pulse_seq = pulse_seq_cls(config=config, propagator=propagator)
        profile = container_cls.from_file(
            path=path,
            config=config,
            pulse_seq=pulse_seq,
            pnames=pnames,
            params_default=params_default,
        )
        cps.set_status(profile.params_default, fit_setting, verbose=False)
        profiles.append(profile)
    return sorted(profiles)


def _read_shifts(config, pulse_seq_cls, propagator_cls, container_cls, fit_setting):
    propagator = propagator_cls.from_config(config)
    shifts = _get_shifts(config)
    profiles = []
    for spin_system, data in shifts.items():
        config["spin_system"]["spin_system"] = spin_system
        pnames, params_default = cp.create_params(
            basis=config["spin_system"]["basis"],
            model=config["model"],
            conditions=config["conditions"],
            spin_system=spin_system,
            constraints=config["spin_system"].get("constraints"),
        )
        pulse_seq = pulse_seq_cls(config=config, propagator=propagator)
        profile = container_cls(
            name=spin_system,
            data=data,
            pulse_seq=pulse_seq,
            pnames=pnames,
            params_default=params_default,
        )
        cps.set_status(profile.params_default, fit_setting, verbose=False)
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
