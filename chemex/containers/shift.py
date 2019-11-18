import copy
import functools as ft

import numpy as np

import chemex.parameters.name as cpn


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
    def __init__(self, name, data, pulse_seq, pnames, params_default):
        self.name = name
        self.data = data
        self.params_default = params_default
        self._pulse_seq = pulse_seq
        self._pnames = pnames

    def residuals(self, params):
        residuals = (self.calculate(params) - self.data["shift"]) / self.data["error"]
        return np.array([residuals])

    def calculate(self, params):
        par_values = self._get_parvals(params)
        calculated = self._pulse_seq.calculate(par_values)
        return calculated

    def estimate_noise_variance(self, kind):
        return self.data.estimate_noise_variance(kind)

    def set_noise(self, value):
        self.data["error"] = value

    def print(self, params):
        shift = self.data["shift"]
        error = self.data["error"]
        value = self.calculate(params)
        output = f"{str(self.name):10s} {shift: 17.8e} {error: 17.8e} {value: 17.8e}\n"
        return output

    def filter(self, params):
        pass

    def monte_carlo(self, params):
        shift_ref = self.calculate(params)
        profile = copy.copy(self)
        profile.data["shift"] = np.random.normal(shift_ref, profile.data["error"])
        return profile

    def set_params(self, params, rates):
        for name1, name2 in self._pnames.items():
            name = cpn.remove_state(name1)
            if name in rates:
                params[name2].value = rates[name]

    def get_cs_value(self, params):
        name = self._pulse_seq.cs_i_state
        fname = self._pnames[name]
        return params[fname]

    def _get_parvals(self, params):
        parvals = tuple(
            (name1, params[name2].value) for name1, name2 in self._pnames.items()
        )
        return parvals

    def __add__(self, other: "ShiftProfile"):
        data = {}
        data["shift"] = 0.5 * (self.data["shift"] + other.data["shift"])
        data["error"] = np.sqrt(self.data["error"] ** 2 + other.data["error"] ** 2)
        return ShiftProfile(
            self.name, data, self._pulse_seq, self._pnames, self.params_default
        )

    def __eq__(self, other: "ShiftProfile"):
        return self.name == other.name

    def __lt__(self, other: "ShiftProfile"):
        return self.name < other.name
