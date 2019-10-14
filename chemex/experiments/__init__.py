import importlib as il
import importlib.resources as ir
import pathlib as pl
import pkgutil as pu
import sys

import chemex.experiments.configs as cec
import chemex.helper as ch
import chemex.parameters.kinetics as cpk


def read(filename, model):
    config = ch.read_toml(filename)
    config["filename"] = pl.Path(filename)
    config["model"] = model
    experiment_name = get_experiment_name(config)
    cpk.validates_conditions(config)
    module = grab(experiment_name)
    experiment = module.read(config)
    return experiment


def get_experiment_name(config):
    filename = config["filename"]
    if "experiment" not in config:
        sys.exit(
            f"\nerror: The experiment file '{filename}' has no section '[experiment]'."
        )
    elif "name" not in config["experiment"]:
        sys.exit(
            f"\nerror: The experiment file '{filename}' has no entry 'name' in the "
            f"section '[experiment]'."
        )
    return config["experiment"]["name"]


def grab(name):
    try:
        module = il.import_module(f"{__package__}.{name}")
    except ModuleNotFoundError:
        sys.exit(
            f"\nerror: '{name}' is not part of our experiment collection! "
            f"Run 'chemex info' to obtain the full list of the available experiments."
        )
    else:
        return module


def get_info():
    docs = {}
    for module in pu.iter_modules(__path__, __name__ + "."):
        if module.ispkg or "helper" in module.name:
            continue
        imported_module = il.import_module(module.name)
        exp_name = module.name.replace(__name__ + ".", "")
        docs[exp_name] = imported_module.__doc__
    return docs


def get_config():
    cfgs = {
        name.replace(".toml", ""): ir.read_text(cec, name)
        for name in ir.contents(cec)
        if name.endswith(".toml")
    }
    return cfgs
