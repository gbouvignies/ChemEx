from __future__ import annotations

import configparser
from argparse import Namespace

import matplotlib.pyplot as plt

import chemex.parameters.name as cpn
import chemex.parameters.spin_system as cns
from chemex.configuration.parameters import ParameterConfigurationError
from chemex.exceptions import ChemExError
from chemex.messages import print_making_plots, print_section


def plot_param(args: Namespace) -> None:
    """Plot values of a parameter versus residue number."""
    params = configparser.ConfigParser()

    if len(args.parameters) > 1:
        raise ChemExError(
            "Multiple parameter files were given. 'chemex plot_param' "
            "requires exactly one parameter file."
        )

    filename = args.parameters.pop()
    try:
        loaded = params.read(str(filename))
        if not loaded:
            raise FileNotFoundError(filename)
    except (OSError, configparser.Error) as error:
        raise ParameterConfigurationError(
            (filename,),
            "The fitted parameter file could not be parsed.",
        ) from error
    param_name = cpn.ParamName.from_section(args.parname)
    curves = {}

    print_making_plots()

    try:
        for section in params.sections():
            section_name = cpn.ParamName.from_section(section.strip('"'))
            if param_name.match(section_name):
                print_section(section)
                residues: list[int] = []
                values: list[float] = []
                errors: list[float] = []
                for key, entry in params.items(section):
                    residues.append(int(cns.SpinSystem(name=key).numbers["i"]))
                    split = entry.split()
                    values.append(float(split[0]))
                    try:
                        uncertainty = float(split[2].strip("±"))
                    except ValueError:
                        uncertainty = 0.0
                    errors.append(uncertainty)
                curves[section] = (residues, values, errors)
    except (IndexError, KeyError, ValueError, configparser.Error) as error:
        raise ParameterConfigurationError(
            (filename,),
            "The fitted parameter values could not be parsed.",
        ) from error

    _, axis = plt.subplots(figsize=(12, 5))
    axis.yaxis.grid(visible=True)

    for section, (residues, values, errors) in curves.items():
        axis.errorbar(
            residues,
            values,
            yerr=errors,
            label=section,
            fmt=".",
            barsabove=True,
        )

    plt.legend()
    plt.show()
