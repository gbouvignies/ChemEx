from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from dataclasses import field
from functools import cached_property
from operator import attrgetter
from typing import Protocol

import numpy as np
from cachetools import cachedmethod
from cachetools import LRUCache
from lmfit import Parameters as ParametersLF
from numpy.typing import NDArray

from chemex.containers.data import Data
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.data import Printer

NDArrayFloat = NDArray[np.float_]
NDArrayBool = NDArray[np.bool_]


class PulseSequence(Protocol):
    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:
        ...

    def is_reference(self, metadata: NDArrayFloat) -> NDArrayBool:
        ...


class Filterer(Protocol):
    def filter(self, data: Data) -> None:
        ...


def _cache_key(self, params: ParametersLF) -> tuple[float, ...]:
    return (
        *(params[param_id].value for param_id in self.param_ids),
        *self.data.metadata,
    )


@dataclass(order=True)
class Profile:
    data: Data = field(compare=False)
    spectrometer: Spectrometer = field(compare=False)
    pulse_sequence: PulseSequence = field(compare=False)
    name_map: dict[str, str] = field(compare=False)
    printer: Printer = field(compare=False)
    filterer: Filterer | None = field(compare=False, default=None)
    is_scaled: bool = field(compare=False, default=True)
    name: SpinSystem = field(compare=True, init=False)
    cache: LRUCache = field(compare=False, init=False)

    def __post_init__(self):
        self.name = self.spectrometer.liouvillian.spin_system
        self.data.refs = self.pulse_sequence.is_reference(self.data.metadata)
        self.cache = LRUCache(maxsize=5)

    @cached_property
    def param_ids(self) -> set[str]:
        return set(self.name_map.values())

    def _get_par_values(self, params: ParametersLF) -> dict[str, float]:
        return {
            local_name: params[param_id].value
            for local_name, param_id in self.name_map.items()
        }

    def update_spectrometer(self, params: ParametersLF) -> None:
        par_values = self._get_par_values(params)
        self.spectrometer.update(par_values)

    def calculate(self, params: ParametersLF) -> np.ndarray:
        self.update_spectrometer(params)
        self.data.calc = self.pulse_sequence.calculate(self.spectrometer, self.data)
        if self.is_scaled:
            self.data.scale_calc()
        return self.data.calc

    @cachedmethod(attrgetter("cache"), key=_cache_key)
    def residuals(self, params: ParametersLF) -> np.ndarray:
        residuals = (self.calculate(params) - self.data.exp) / self.data.err
        return residuals[self.data.mask]

    def filter(self, params: ParametersLF) -> None:
        if self.filterer is None:
            return
        self.update_spectrometer(params)
        self.filterer.filter(self.data)

    def set_noise(self, value: float):
        self.data.err[:] = value

    def monte_carlo(self) -> Profile:
        profile = deepcopy(self)
        profile.data = profile.data.monte_carlo()
        return profile

    def bootstrap(self) -> Profile:
        """Make a profile for bootstrap analysis."""
        profile = deepcopy(self)
        profile.data = profile.data.bootstrap()
        return profile

    def any_duplicate(self):
        return self.data.any_duplicate()

    def __add__(self: Profile, other: Profile) -> Profile:
        profile = deepcopy(self)
        profile.data = self.data + other.data
        return profile

    def __str__(self) -> str:
        return self.printer.print(str(self.name), self.data)
