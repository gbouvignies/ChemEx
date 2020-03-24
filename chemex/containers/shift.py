import copy
import functools as ft

import numpy as np


SHIFT_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "duplicates"],
                    "default": "file",
                },
                "path": {"type": "string", "default": "./"},
            },
            "required": ["path"],
        }
    },
}


@ft.total_ordering
class ShiftProfile:
    def __init__(self, name, data, pulse_seq, pnames, params):
        self.name = name
        self.data = data
        self._pulse_seq = pulse_seq
        self._pnames = pnames
        self.params = params

    def residuals(self, params):
        residuals = (self.calculate(params) - self.data["shift"]) / self.data["error"]
        return np.array([residuals])

    def calculate(self, params):
        par_values = self._get_parvals(params)
        return self._pulse_seq.calculate(par_values)

    def estimate_noise_variance(self, kind):
        return self.data.estimate_noise_variance(kind)

    def set_noise(self, value):
        self.data["error"] = value

    def print(self, params):
        shift = self.data["shift"]
        error = self.data["error"]
        value = self.calculate(params)
        return f"{str(self.name):10s} {shift: 17.8e} {error: 17.8e} {value: 17.8e}\n"

    def filter(self, params):
        pass

    def monte_carlo(self, params):
        shift_ref = self.calculate(params)
        profile = copy.copy(self)
        profile.data["shift"] = np.random.normal(shift_ref, profile.data["error"])
        return profile

    def get_cs_value(self, params):
        name = self._pulse_seq.cs_i_state
        fname = self._pnames[name]
        return params[fname]

    def _get_parvals(self, params):
        return tuple(
            (name1, params[name2].value) for name1, name2 in self._pnames.items()
        )

    def __add__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        data = {
            "shift": 0.5 * (self.data["shift"] + other.data["shift"]),
            "error": np.sqrt(self.data["error"] ** 2 + other.data["error"] ** 2),
        }
        return ShiftProfile(self.name, data, self._pulse_seq, self._pnames, self.params)

    def __eq__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __lt__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name < other.name