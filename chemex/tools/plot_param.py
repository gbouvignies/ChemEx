import argparse
import configparser
import pathlib

from matplotlib import pyplot as plt

from chemex import parameters, peaks


def plot_param(args):

    params = configparser.ConfigParser()
    params.read(str(args.parameters))

    parname = args.parname

    curves = {}

    print("Plotting...")

    for section in params.sections():

        short_name = parameters.ParameterName().from_section(section).name

        if parname.lower() in short_name:
            print("".join(["  [", section, "]"]))
            points = []

            for key, value in params.items(section):

                res = int(peaks.Peak(key).numbers["i"])
                split = value.split()
                value = float(split[0])

                try:
                    error = float(split[2])
                except ValueError:
                    error = 0.0

                points.append((res, value, error))

            curves[section] = zip(*points)

    _, axis = plt.subplots(figsize=(12, 5))

    axis.yaxis.grid(True)

    for section, (res, vals, errors) in curves.items():
        axis.errorbar(res, vals, yerr=errors, label=section, fmt=".", barsabove=True)

    plt.legend()
    plt.show()

    return 0
