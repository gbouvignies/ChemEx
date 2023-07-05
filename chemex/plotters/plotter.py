from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from pathlib import Path

    from chemex.containers.profile import Profile


class Plotter(Protocol):
    def plot(self, path: Path, profiles: list[Profile]) -> None:
        ...

    def plot_simulation(self, path: Path, profiles: list[Profile]) -> None:
        ...
