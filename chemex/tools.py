import argparse
import configparser
import pathlib

import contextlib
import re

import numpy as np
from matplotlib.widgets import Button
from scipy import interpolate

import matplotlib.pyplot as plt

from chemex import peaks, parameters


def plot_param():

    parser = argparse.ArgumentParser(
        description="Plot one selected parameter from a 'parameters.fit' file"
    )
    parser.add_argument("--file", "-f", type=pathlib.Path, required=True)
    parser.add_argument("--par", "-p", required=True)

    args = parser.parse_args()

    params = configparser.ConfigParser()
    params.read(str(args.file))
    parname = args.par
    curves = {}

    for section in params.sections():

        short_name = parameters.ParameterName().from_section(section).name

        if parname.lower() in short_name:
            print("".join(["[", section, "]"]))
            points = []

            for key, value in params.items(section):

                res = int(peaks.Peak(key).numbers["i"])
                split = value.split()
                value = float(split[0])

                try:
                    error = float(split[2])
                except ValueError:
                    error = 0.0
                    print(res, value, error)

                points.append((res, value, error))

            curves[section] = zip(*points)

    _, axis = plt.subplots(figsize=(12, 5))

    axis.yaxis.grid(True)

    for section, (res, vals, errors) in curves.items():
        axis.errorbar(res, vals, yerr=errors, label=section, fmt=".", barsabove=True)

    plt.legend()
    plt.show()

    return 0


def cest_pick():

    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", nargs="+", type=pathlib.Path)
    parser.add_argument("--carrier", type=float)
    parser.add_argument("--frq", type=float)
    parser.add_argument("--out", type=pathlib.Path, default="Out")
    args = parser.parse_args()

    callback = Buttons(args)

    axprevious = plt.axes([0.825, 0.1, 0.075, 0.075])
    axnext = plt.axes([0.9, 0.1, 0.075, 0.075])

    bprevious = Button(axprevious, "Previous")
    bprevious.on_clicked(callback.previous)

    bnext = Button(axnext, "Next")
    bnext.on_clicked(callback.next)

    plt.show()


def get_numbers(name):
    return [int(number) for number in re.findall("\d+", name)]


def prepare_profile(filename, frq, carrier):

    data = np.loadtxt(filename)
    mask = data[:, 0] > -10000.0
    norm = np.mean(data[~mask][:, 1])
    data = data[mask]
    data[:, 1:] /= norm

    offsets = data[:, 0] / frq + carrier

    profile = {"x": offsets, "y": data[:, 1], "e": data[:, 2]}

    spline = interpolate.CubicSpline(profile["x"], profile["y"])

    fine_offsets = np.linspace(min(offsets), max(offsets), 1000)

    profile_cs = {"x": fine_offsets, "y": spline(fine_offsets)}

    return profile, profile_cs


class Buttons(object):
    def __init__(self, args):

        self.out = args.out

        self.profile, self.profile_spline = None, None
        self.cs_a, self.cs_b = {}, {}
        self.positions = []
        self.cid = None

        self.fig, self.axis = plt.subplots(figsize=(12, 5))
        self.fig.subplots_adjust(left=0.07, bottom=0.1, right=0.8, top=0.9)

        self.filenames = {}

        for name in args.inp:
            short_name = name.stem.split(".")[0].split("-")[0]
            self.filenames[short_name] = name

        self.frq = args.frq

        self.carrier = args.carrier

        self.names = sorted(self.filenames, key=get_numbers)
        self.index = -1
        self.next(None)

    def previous(self, event):
        if self.index > 0:
            self.index -= 1
        else:
            self.index = len(self.names) - 1
        self.name = self.names[self.index]
        self.prepare_profiles(self.name)
        self.clear()
        self.plot_profiles_cs()

    def next(self, event):
        if self.index < len(self.names) - 1:
            self.index += 1
        else:
            self.index = 0
        self.name = self.names[self.index]
        self.prepare_profiles(self.name)
        self.clear()
        self.plot_profiles_cs()

    def clear(self):
        while self.positions:
            self.positions.pop().remove()
        self.axis.clear()

    def prepare_profiles(self, name):
        self.profile, self.profile_spline = prepare_profile(
            self.filenames[name], self.frq, self.carrier
        )

    def plot_profiles(self):

        self.axis.plot(self.profile["x"], self.profile["y"], ".", color="C0")

        self.axis.set_title(self.name)

        self.axis.set_xlabel(r"$^{15}$N (ppm)")
        self.axis.set_ylabel("$I/I_0$")

        self.fig.canvas.draw()

    def plot_profiles_cs(self, event=None):

        self.cid = self.fig.canvas.mpl_connect(
            "button_press_event", self.plot_positions
        )

        self.clear()
        self.plot_positions()
        self.axis.plot(
            self.profile_spline["x"],
            self.profile_spline["y"],
            "--",
            color="C0",
            alpha=0.66,
        )

        self.plot_profiles()

        self.fig.canvas.draw()

    def plot_positions(self, event=None):

        cs_a = self.cs_a.get(self.name)
        cs_b = self.cs_b.get(self.name)

        if event and event.inaxes == self.axis:
            if cs_a is None:
                self.cs_a[self.name] = event.xdata
            elif cs_b is None:
                self.cs_b[self.name] = event.xdata
            else:
                while self.positions:
                    self.positions.pop().remove()
                self.cs_a[self.name] = event.xdata
                self.cs_b[self.name] = None

        cs_a = self.cs_a.get(self.name)
        cs_b = self.cs_b.get(self.name)

        kwargs = {"color": "0.74", "linewidth": 1.0}

        if cs_a is not None:
            self.positions.append(self.axis.axvline(cs_a, **kwargs))
            text = "".join([r"$\varpi_a$ = ", "{:.3f} ppm".format(float(cs_a))])
            self.positions.append(self.fig.text(0.82, 0.8, text))

        if cs_b is not None:
            self.positions.append(self.axis.axvline(cs_b, linestyle="dashed", **kwargs))
            text = "".join([r"$\varpi_b$ = ", "{:.3f} ppm".format(float(cs_b))])
            self.positions.append(self.fig.text(0.82, 0.75, text))

        if cs_a is not None and cs_b is not None:
            text = "".join(
                [r"$\Delta\varpi_{ab}$ = ", "{:.3f} ppm".format(float(cs_b - cs_a))]
            )
            self.positions.append(self.fig.text(0.82, 0.7, text))
            self.save()

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def save(self):

        self.out.mkdir(parents=True, exist_ok=True)

        fname1 = self.out / "cs_a.txt"
        fname2 = self.out / "cs_b.txt"
        fname3 = self.out / "dw_ab.txt"

        with contextlib.ExitStack() as stack:
            file1 = stack.enter_context(fname1.open("w"))
            file2 = stack.enter_context(fname2.open("w"))
            file3 = stack.enter_context(fname3.open("w"))

            for name in self.names:
                cs_a = self.cs_a.get(name)
                cs_b = self.cs_b.get(name)

                if cs_a is None:
                    continue

                cs_a = float(cs_a)

                if cs_b is None:
                    cs_b = cs_a
                else:
                    cs_b = float(cs_b)

                dw_ab = cs_b - cs_a

                file1.write(f"{name:10s} {cs_a:8.3f}\n")
                file2.write(f"{name:10s} {cs_b:8.3f}\n")
                file3.write(f"{name:10s} {dw_ab:8.3f}\n")
