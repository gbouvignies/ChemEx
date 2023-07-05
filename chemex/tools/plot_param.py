from __future__ import annotations

import configparser
import sys
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

import chemex.parameters.name as cpn
import chemex.parameters.spin_system as cns
from chemex.messages import print_making_plots, print_section

if TYPE_CHECKING:
    from argparse import Namespace


def plot_param(args: Namespace):
    """Plot values of a parameter versus residue number."""
    params = configparser.ConfigParser()

    if len(args.parameters) > 1:
        sys.exit(
            "\nError: Multiple parameter files were given. 'chemex plot_param' "
            "should only be run with a single parameter file.\n"
        )

    params.read(str(args.parameters.pop()))
    param_name = cpn.ParamName.from_section(args.parname)
    curves = {}

    print_making_plots()

    for section in params.sections():
        section_name = cpn.ParamName.from_section(section.strip('"'))
        if param_name.match(section_name):
            print_section(section)
            residues: list[int] = []
            values: list[float] = []
            errors: list[float] = []
            for key, entry in params.items(section):
                residues.append(int(cns.SpinSystem(key).numbers["i"]))
                split = entry.split()
                values.append(float(split[0]))
                try:
                    error = float(split[2].strip("±"))
                except ValueError:
                    error = 0.0
                errors.append(error)
            curves[section] = (residues, values, errors)

    _, axis = plt.subplots(figsize=(12, 5))
    axis.yaxis.grid(True)

    for section, (residues, values, errors) in curves.items():
        axis.errorbar(
            residues, values, yerr=errors, label=section, fmt=".", barsabove=True
        )

    plt.legend()
    plt.show()
