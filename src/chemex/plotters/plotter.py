from __future__ import annotations

from pathlib import Path
from typing import Protocol

from chemex.containers.profile import Profile


class Plotter(Protocol):
    def plot(self, output_stem: Path, profiles: list[Profile]) -> None: ...

    def plot_simulation(self, output_stem: Path, profiles: list[Profile]) -> None: ...
