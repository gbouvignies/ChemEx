import contextlib
from pathlib import Path
from typing import Literal

import numpy as np
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.backend_bases import Event, LocationEvent
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

from chemex.containers.experiments import Experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus

from .curve import Curve

LSTYLE = {"a": "-", "b": ":"}
TEXT_Y = {"a": 0.8, "b": 0.75}
XLABELS = {
    Nucleus.N15: r"$^{15}$N (ppm)",
    Nucleus.C13: r"$^{13}$C (ppm)",
    Nucleus.H1: r"$^{1}$H (ppm)",
    Nucleus.X: r"Chemical shift (ppm)",
}


class Buttons:
    def __init__(
        self,
        figure: Figure,
        axis: Axes,
        experiments: Experiments,
        path: Path,
        sw: float | None,
    ) -> None:
        self.fig = figure
        self.axis = axis
        self.data: dict[SpinSystem, list[Curve]] = {}
        spin_systems: set[SpinSystem] = set()

        for experiment in experiments:
            for profile in experiment:
                spin_system = profile.spin_system
                spin_systems.add(spin_system)
                self.data.setdefault(spin_system, []).append(Curve(profile, sw))

        self.spin_systems = sorted(spin_systems)
        self.out = path

        self.spin_system = SpinSystem("")
        self.curves: list[Curve] = []
        self.cs_a: dict[SpinSystem, float | None] = {}
        self.cs_b: dict[SpinSystem, float | None] = {}
        self.artists: list[Artist] = []
        self.sw = sw

        self.cursor = Cursor(self.axis, horizOn=False, useblit=True)

        self.index = -1
        self.next()

    def _clear_artists(self) -> None:
        while self.artists:
            self.artists.pop().remove()

    def _clear_axis(self) -> None:
        self._clear_artists()
        self.axis.clear()

    def _show_labels(self) -> None:
        self.axis.set_title(str(self.spin_system))
        self.axis.set_xlabel(XLABELS[self.spin_system.nuclei["i"]])
        self.axis.set_ylabel("$I/I_0$")

    def _plot_profiles(self) -> None:
        if not self.curves:
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

    def _get_click_position(self, event: Event) -> float | None:
        if not isinstance(event, LocationEvent):
            return None
        if event.inaxes != self.axis:
            return None
        return event.xdata

    def _add_line(self, position: float, state: Literal["a", "b"]):
        text_ = rf"$\varpi_{state}$ = {position:.3f} ppm"
        text = self.fig.text(0.82, TEXT_Y[state], text_)
        line = self.axis.axvline(
            x=position,
            linestyle=LSTYLE[state],
            color="0.74",
            linewidth=1.0,
        )
        self.artists.extend([line, text])

    def _add_text_dw(self, dw_ab: float) -> None:
        text_ = rf"$\Delta\varpi_{{ab}}$ = {dw_ab:.3f} ppm"
        text = self.fig.text(0.82, 0.7, text_)
        self.artists.append(text)

    def _save(self):
        self.out.mkdir(parents=True, exist_ok=True)

        fname1 = self.out / "cs_a.toml"
        fname2 = self.out / "dw_ab.toml"

        with contextlib.ExitStack() as stack:
            file1 = stack.enter_context(fname1.open("w", encoding="utf-8"))
            file2 = stack.enter_context(fname2.open("w", encoding="utf-8"))
            file1.write("[CS_A]\n")
            file2.write("[DW_AB]\n")

            for name in self.spin_systems:
                cs_a = self.cs_a.get(name)
                cs_b = self.cs_b.get(name)

                if cs_a is None:
                    continue

                file1.write(f"{str(name).upper():10s} = {cs_a:8.3f}\n")

                if cs_b is not None:
                    dw_ab = cs_b - cs_a
                    file2.write(f"{str(name).upper():10s} = {dw_ab:8.3f}\n")

    def set_cs(self, event: Event) -> None:
        """Set the chemical shift."""
        xdata = self._get_click_position(event)
        if xdata is None:
            return

        key = self.spin_system

        if self.cs_a.get(key) is None or self.cs_b.get(key) is not None:
            self.cs_a[key] = xdata
            self.cs_b[key] = None
        else:
            self.cs_b[key] = xdata

        self._plot_lines()

    def _plot_lines(self) -> None:
        key = self.spin_system
        cs_a = self.cs_a.get(key)
        cs_b = self.cs_b.get(key)

        self._clear_artists()

        if cs_a is not None:
            self._add_line(cs_a, "a")

        if cs_b is not None:
            self._add_line(cs_b, "b")

        if cs_a is not None and cs_b is not None:
            self._add_text_dw(cs_b - cs_a)

        self._save()

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def _plot(self, event: Event | None = None) -> None:
        self._clear_axis()
        self._plot_lines()
        self._plot_profiles()
        self.fig.canvas.draw_idle()

    def _shift(self, step: int):
        self.index += step
        self.index %= len(self.spin_systems)
        self.spin_system = self.spin_systems[self.index]
        self.curves = self.data[self.spin_system]
        self._clear_axis()
        self._plot()

    def next(self, _event: Event | None = None) -> None:
        """Go to next residue."""
        self._shift(+1)

    def previous(self, _event: Event) -> None:
        """Go to previous residue."""
        self._shift(-1)

    def swap(self, event: Event) -> None:
        """Swap peak peak positions for major/minor states."""
        name = self.spin_system
        if self.cs_b[name] is not None:
            self.cs_a[name], self.cs_b[name] = self.cs_b[name], self.cs_a[name]
        self._plot_lines()

    def clear(self, event: Event) -> None:
        name = self.spin_system
        self.cs_a[name], self.cs_b[name] = None, None
        self._plot_lines()

    def set_sw(self, sw: float) -> None:
        with contextlib.suppress(ValueError):
            self.sw = sw

        self._clear_axis()
        self._plot()
