from __future__ import annotations

from pathlib import Path
from typing import Protocol

from chemex.containers.profile import Profile


class Plotter(Protocol):
    def plot(self, path: Path, profiles: list[Profile]) -> None:
        ...
