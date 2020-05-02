import contextlib
import sys

import matplotlib.pyplot as plt
import matplotlib.widgets as mw
import numpy as np
import scipy.interpolate as si

import chemex.experiments as ce
import chemex.parameters.kinetics as cpk
from chemex.nmr import spin_system


def pick_cest(args):

    if len(args.experiments) > 1:
        sys.exit(
            "\nError: Multiple experiment files were given. 'chemex pick_cest' "
            "should only be run with a single experiment files.\n"
        )

    # Read experimental setup and data
    model = cpk.parse_model(name="2st")
    # Read experimental setup and data
    experiment = ce.read(filename=args.experiments.pop(), model=model)

    if not experiment.name.startswith("cest"):
        sys.exit(
            "\nError: The command 'chemex pick_cest' only works with CEST "
            "experiments.\n"
        )

    callback = Buttons(experiment, args.out_dir)

    axclear = plt.axes([0.825, 0.3, 0.15, 0.075])
    axprevious = plt.axes([0.825, 0.2, 0.075, 0.075])
    axnext = plt.axes([0.9, 0.2, 0.075, 0.075])
    axclose = plt.axes([0.825, 0.1, 0.15, 0.075])

    bclear = mw.Button(axclear, "Clear")
    bclear.on_clicked(callback.clear)

    bprevious = mw.Button(axprevious, "Previous")
    bprevious.on_clicked(callback.previous)

    bnext = mw.Button(axnext, "Next")
    bnext.on_clicked(callback.next)

    bclose = mw.Button(axclose, "Quit")
    bclose.on_clicked(callback.close)

    plt.show()


class Buttons:

    KWARGS = {"color": "0.74", "linewidth": 1.0}
    LSTYLE = {"a": "-", "b": ":"}
    TEXT_Y = {"a": 0.8, "b": 0.75}

    def __init__(self, data, path):

        self.data = data._profiles
        self.out = path

        self.names = [spin_system.SpinSystem(profile.name).names["i"] for profile in self.data]

        self.curve, self.curve_sp = None, None
        self.cs_a, self.cs_b = {}, {}
        self.lines = []
        self.cid = None

        self.fig, self.axis = plt.subplots(figsize=(12, 5))
        self.fig.subplots_adjust(left=0.07, bottom=0.1, right=0.8, top=0.9)

        self.index = -1
        self.next(None)

    def clear(self, event):
        name = self.name.names["i"]
        self.cs_a[name], self.cs_b[name] = None, None
        self._clear_lines()

    def previous(self, event):
        self._shift(-1)

    def next(self, event):
        self._shift(+1)

    def close(self, event):
        plt.close()

    def _shift(self, step):
        self.index += step
        self.index %= len(self.data)
        self._profile = self.data[self.index]
        self.name = spin_system.SpinSystem(self._profile.name)
        self._profile_to_curve()
        self._clear_axis()
        self._plot()

    def _plot(self, event=None):
        self.cid = self.fig.canvas.mpl_connect("button_press_event", self._plot_lines)
        self._clear_axis()
        self._plot_lines()
        self._plot_profile()
        self.fig.canvas.draw()

    def _plot_profile(self):

        x, y = self.curve["x"], self.curve["y"]
        x_sp, y_sp = self.curve_sp["x"], self.curve_sp["y"]

        self.axis.plot(x, y, ".", color="C0")
        self.axis.plot(x_sp, y_sp, "--", color="C0", alpha=0.66)

        self._show_labels()

        self.fig.canvas.draw()

    def _show_labels(self):
        self.axis.set_title(str(self.name).upper())
        atom = (spin_system.SpinSystem(self._profile.name).atoms["i"])
        if atom == "N":
            self.axis.set_xlabel(r"$^{15}$N (ppm)")
        elif atom == "C":
            self.axis.set_xlabel(r"$^{13}$C (ppm)")
        elif atom == "H":
            self.axis.set_xlabel(r"$^{1}$H (ppm)")
        else:
            self.axis.set_xlabel(r"Chemical shift (ppm)")
        self.axis.set_ylabel("$I/I_0$")

    def _profile_to_curve(self):
        profile = self._profile
        data_exp = profile._get_plot_data_exp()
        _, indices = np.unique(data_exp["ppms"], return_index=True)
        spline = si.CubicSpline(
            data_exp["ppms"][indices], data_exp["intensities"][indices]
        )
        fine_offsets = np.linspace(min(data_exp["ppms"]), max(data_exp["ppms"]), 1000)
        self.curve = {
            "x": data_exp["ppms"],
            "y": data_exp["intensities"],
            "e": data_exp["errors"],
        }
        self.curve_sp = {"x": fine_offsets, "y": spline(fine_offsets)}

    def _plot_lines(self, event=None):

        name = self.name.names["i"]

        cs_a = self.cs_a.get(name)
        cs_b = self.cs_b.get(name)

        xdata = self._get_click_position(event)

        if xdata is not None:
            if cs_a is None or cs_b is not None:
                cs_a = xdata
                cs_b = None
            else:
                cs_b = xdata

        self._clear_lines()

        if cs_a is not None:
            self._add_line(cs_a, "a")

        if cs_b is not None:
            dw_ab = cs_b - cs_a
            self._add_line(cs_b, "b")
            self._add_text_dw(dw_ab)

        self.cs_a[name] = cs_a
        self.cs_b[name] = cs_b

        self._save()

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def _get_click_position(self, event=None):
        result = None
        if event and event.inaxes == self.axis:
            result = event.xdata
        return result

    def _add_line(self, position, state):
        text_ = fr"$\varpi_{state}$ = {position:.3f} ppm"
        text = self.fig.text(0.82, self.TEXT_Y[state], text_)
        line = self.axis.axvline(position, linestyle=self.LSTYLE[state], **self.KWARGS)
        self.lines.extend([line, text])

    def _add_text_dw(self, dw_ab):
        text_ = fr"$\Delta\varpi_{{ab}}$ = {dw_ab:.3f} ppm"
        text = self.fig.text(0.82, 0.7, text_)
        self.lines.append(text)

    def _clear_axis(self):
        self._clear_lines()
        self.axis.clear()

    def _clear_lines(self):
        while self.lines:
            self.lines.pop().remove()

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
                elif cs_b is None:
                    cs_b = cs_a

                dw_ab = cs_b - cs_a

                file1.write("{:10s} = {:8.3f}\n".format(name.upper(), cs_a))
                file2.write("{:10s} = {:8.3f}\n".format(name.upper(), dw_ab))
