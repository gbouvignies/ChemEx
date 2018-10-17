"""TODO: module docstring."""
import abc
import copy
from functools import lru_cache

import numpy as np

from chemex import peaks
from chemex.spindynamics import util

EXP_DETAILS = {"name": {"type": str}}

CONDITIONS = {
    "h_larmor_frq": {"type": float},
    "temperature": {"type": float},
    "p_total": {"default": None, "type": float},
    "l_total": {"default": None, "type": float},
}


class BaseProfile(metaclass=abc.ABCMeta):
    """TODO: class docstring."""

    subclasses = []

    BOOLEAN_STATES = {
        "1": True,
        "yes": True,
        "true": True,
        "on": True,
        "0": False,
        "no": False,
        "false": False,
        "off": False,
    }

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses.append(cls)

    def __init__(self, name=None, data=None, exp_details=None, model=None):
        """TODO: method docstring."""

        if name is None:
            name = ""

        self.name = name
        self.peak = peaks.Peak(self.name)
        self.model = util.parse_model(model)
        self.experiment_name = self.check_exp_details(exp_details, EXP_DETAILS)["name"]
        self.conditions = self.check_exp_details(exp_details, CONDITIONS)
        self.data = data
        self.mask = None
        self.map_names = {}
        self.basis = None

        self.calculate_unscaled_profile = lru_cache(256)(
            self._calculate_unscaled_profile
        )

    def __len__(self):
        return self.data.size

    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and back-calculated
        values."""

        values = self.calculate_profile(params)
        residuals = (self.data["intensities"] - values) / self.data["errors"]

        return residuals[self.mask]

    def calculate_profile(self, params=None, **kwargs):
        """Calculate the CEST profile."""
        params_local = tuple(
            (name_s, params[name_l].value) for name_s, name_l in self.map_names.items()
        )

        values = self.calculate_unscaled_profile(params_local)

        try:
            scale = sum(
                values * self.data["intensities"] / self.data["errors"] ** 2
            ) / sum((values / self.data["errors"]) ** 2)

        except ZeroDivisionError:
            scale = 0.0

        if kwargs:
            values = self._calculate_unscaled_profile(params_local, **kwargs)

        return values * scale

    @abc.abstractmethod
    def _calculate_unscaled_profile(self, params_local, **kwargs):
        """Calculate the unscaled CEST profile."""
        pass

    @abc.abstractmethod
    def print_profile(self, params=None):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def filter_points(self, params=None):
        """TODO: method docstring."""
        pass

    def make_mc_profile(self, params):
        """Make a profile for MC analysis."""

        profile = copy.copy(self)
        profile.data["intensities"] = (
            self.calculate_profile(params)
            + np.random.randn(len(self.data["intensities"])) * self.data["errors"]
        )

        return profile

    @staticmethod
    def check_exp_details(exp_details=None, expected=None):
        """Check for experimental parameters and convert them to their appropriate
        type."""

        if expected is None:
            expected = {}

        details = {}
        missing = []

        for name, item in expected.items():
            default = item.get("default")
            value = exp_details.get(name, default)
            required = "default" not in item
            if value is None and required:
                missing.append("'" + name + "'")
            elif isinstance(value, str):
                value = value.split()
            if value is None:
                details[name] = None
            else:
                details[name] = np.array(value, dtype=item.get("type")).reshape(-1)

                if details[name].shape == (1,):
                    details[name] = details[name][0]

        if missing:
            exit(
                "\n  - Missing experimental parameter detected. Please set: "
                "{:s}".format(", ".join(missing))
            )

        return details

    def get_bool(self, value):

        if value.lower() not in self.BOOLEAN_STATES:
            raise ValueError("Not a boolean: %s" % value)

        return self.BOOLEAN_STATES[value.lower()]
