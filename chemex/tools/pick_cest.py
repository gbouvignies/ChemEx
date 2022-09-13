from __future__ import annotations

import contextlib
import sys
from argparse import Namespace

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button
from matplotlib.widgets import Cursor
from matplotlib.widgets import Slider
from scipy.interpolate import CubicSpline

from chemex import model
from chemex.configuration.methods import Selection
from chemex.containers.profile import Profile
from chemex.experiments.builder import build_experiments
from chemex.experiments.loader import register_experiments
from chemex.parameters.spin_system import Nucleus
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import create_plot_data_exp


class Curve:
    def __init__(self, profile: Profile, sw: float | None = None):
        data = create_plot_data_exp(profile)
        ppms = data.metadata
        xmin, xmax = min(ppms), max(ppms)
        self.xrange = xmax - xmin
        self.xcentre = 0.5 * (xmin + xmax)
        self.x = ppms
        self.y = data.exp
        _, indices = np.unique(ppms, return_index=True)
        x = ppms[indices]
        y = data.exp[indices]
        if sw is not None:
            y[0] = y[-1] = 0.5 * (y[0] + y[-1])
            self.spline = CubicSpline(x, y, bc_type="periodic", extrapolate="periodic")
        else:
            self.spline = CubicSpline(x, y)

    def get_xrange(self, sw=None):
        if sw is None:
            sw = 1.0
        return self.xcentre + sw * self.xrange * np.array([-0.5, 0.5])


class Buttons:

    KWARGS = {"color": "0.74", "linewidth": 1.0}
    LSTYLE = {"a": "-", "b": ":"}
    TEXT_Y = {"a": 0.8, "b": 0.75}

    def __init__(self, experiments, path, sw):

        self.data = {}
        names = set()

        for experiment in experiments:
            for profile in experiment:
                name = SpinSystem(profile.name.names["i"])
                names.add(name)
                self.data.setdefault(name, []).append(Curve(profile, sw))

        self.names = sorted(names)
        self.out = path

        self.name = None
        self.curves = None
        self.cs_a, self.cs_b = {}, {}
        self.lines = []
        self.cid = None
        self.sw = sw

        self.fig, self.axis = plt.subplots(figsize=(12, 5))
        self.fig.subplots_adjust(left=0.07, bottom=0.1, right=0.8, top=0.9)

        self.cursor = Cursor(self.axis, horizOn=False, useblit=True)

        self.index = -1
        self.next(None)

    def _clear_lines(self):
        while self.lines:
            self.lines.pop().remove()

    def _clear_axis(self):
        self._clear_lines()
        self.axis.clear()

    def _show_labels(self):
        self.axis.set_title(str(self.name).upper())
        atom = SpinSystem(str(self.name)).atoms["i"]
        if atom.nucleus == Nucleus.N15:
            self.axis.set_xlabel(r"$^{15}$N (ppm)")
        elif atom.nucleus == Nucleus.C13:
            self.axis.set_xlabel(r"$^{13}$C (ppm)")
        elif atom.nucleus == Nucleus.H1:
            self.axis.set_xlabel(r"$^{1}$H (ppm)")
        else:
            self.axis.set_xlabel(r"Chemical shift (ppm)")
        self.axis.set_ylabel("$I/I_0$")

    def _plot_profiles(self):
        if self.curves is None:
            return
        xranges = np.concatenate([curve.get_xrange(self.sw) for curve in self.curves])
        grid = np.linspace(min(xranges), max(xranges), 1000)

        for index, curve in enumerate(self.curves):
            color = f"C{index}"
            self.axis.plot(curve.x, curve.y, ".", color=color)
            self.axis.plot(grid, curve.spline(grid), "--", color=color, alpha=0.66)
        self.axis.invert_xaxis()
        self._show_labels()
        self.fig.canvas.draw_idle()

    def _get_click_position(self, event=None):
        return event.xdata if event and event.inaxes == self.axis else None

    def _add_line(self, position, state):
        text_ = rf"$\varpi_{state}$ = {position:.3f} ppm"
        text = self.fig.text(0.82, self.TEXT_Y[state], text_)
        line = self.axis.axvline(position, linestyle=self.LSTYLE[state], **self.KWARGS)
        self.lines.extend([line, text])

    def _add_text_dw(self, dw_ab):
        text_ = rf"$\Delta\varpi_{{ab}}$ = {dw_ab:.3f} ppm"
        text = self.fig.text(0.82, 0.7, text_)
        self.lines.append(text)

    def _save(self):

        self.out.mkdir(parents=True, exist_ok=True)

        fname1 = self.out / "cs_a.toml"
        fname2 = self.out / "dw_ab.toml"

        with contextlib.ExitStack() as stack:

            file1 = stack.enter_context(fname1.open("w"))
            file2 = stack.enter_context(fname2.open("w"))
            file1.write("[CS_A]\n")
            file2.write("[DW_AB]\n")

            for name in self.names:

                cs_a = self.cs_a.get(name)
                cs_b = self.cs_b.get(name)

                if cs_a is None:
                    continue

                file1.write(f"{str(name).upper():10s} = {cs_a:8.3f}\n")

                if cs_b is not None:
                    dw_ab = cs_b - cs_a
                    file2.write(f"{str(name).upper():10s} = {dw_ab:8.3f}\n")

    def _plot_lines(self, event=None):

        cs_a = self.cs_a.get(self.name)
        cs_b = self.cs_b.get(self.name)

        xdata = self._get_click_position(event)

        if xdata is not None:
            if cs_a is None:
                cs_a = xdata
            else:
                cs_b = xdata

        self._clear_lines()

        if cs_a is not None:
            self._add_line(cs_a, "a")

        if cs_b is not None:
            dw_ab = cs_b - cs_a
            self._add_line(cs_b, "b")
            self._add_text_dw(dw_ab)

        self.cs_a[self.name] = cs_a
        self.cs_b[self.name] = cs_b

        self._save()

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def _plot(self, event=None):
        self.cid = self.fig.canvas.mpl_connect("button_press_event", self._plot_lines)
        self._clear_axis()
        self._plot_lines()
        self._plot_profiles()
        self.fig.canvas.draw_idle()

    def _shift(self, step):
        self.index += step
        self.index %= len(self.names)
        self.name = self.names[self.index]
        self.curves = self.data[self.name]
        self._clear_axis()
        self._plot()

    def next(self, event):
        """Go to next residue."""
        self._shift(+1)

    def previous(self, event):
        """Go to previous residue."""
        self._shift(-1)

    def swap(self, event):
        """Swap peak peak positions for major/minor states."""
        name = self.name
        if self.cs_b[name] is not None:
            self.cs_a[name], self.cs_b[name] = self.cs_b[name], self.cs_a[name]
        self._plot_lines()

    def clear(self, event):
        name = self.name
        self.cs_a[name], self.cs_b[name] = None, None
        self._plot_lines()

    def set_sw(self, text):
        try:
            self.sw = float(text)
        except ValueError:
            pass
        self._clear_axis()
        self._plot()


def pick_cest(args: Namespace):
    """Pick peak positions in CEST profiles."""

    register_experiments()

    # Read experimental setup and data
    model.set_model("2st")

    # Read experimental setup and data
    no_selection = Selection(include=None, exclude=None)
    experiments = build_experiments(args.experiments, no_selection)

    sw = None

    for experiment in experiments:
        if not experiment.name.startswith(("cest", "dcest", "coscest")):
            sys.exit(
                f"\nError: '{experiment.name}' experiment not supported. "
                "The command 'chemex pick_cest' only works with CEST experiments.\n"
            )
        if experiment.name.startswith("dcest") or experiment.name.startswith("coscest"):
            sw = 3.0

    callback = Buttons(experiments, args.output, sw)

    if sw is not None:
        axsw = plt.axes((0.855, 0.425, 0.1, 0.02))
        sw_slider = Slider(
            ax=axsw, label="Scale", valmin=1.0, valmax=10.0, valfmt="%3.1f", valinit=sw
        )
        # bsw = Textbox(axsw, "Number of SWs:", initial="3.0", label_pad=0.15)
        sw_slider.on_changed(callback.set_sw)

    axprevious = plt.axes((0.825, 0.3, 0.075, 0.075))
    bprevious = Button(axprevious, "Previous")
    bprevious.on_clicked(callback.previous)

    axnext = plt.axes((0.9, 0.3, 0.075, 0.075))
    bnext = Button(axnext, "Next")
    bnext.on_clicked(callback.next)

    axswap = plt.axes((0.825, 0.2, 0.15, 0.075))
    bswap = Button(axswap, "Swap")
    bswap.on_clicked(callback.swap)

    axclear = plt.axes((0.825, 0.1, 0.15, 0.075))
    bclear = Button(axclear, "Clear")
    bclear.on_clicked(callback.clear)

    plt.show()
