"""TODO: module docstring."""

import abc

import numpy as np

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

    def __init__(self, name=None, exp_details=None, model=None):
        """TODO: method docstring."""

        if name is None:
            name = ""

        self.profile_name = name
        self.model = util.parse_model(model)
        self.experiment_name = self.check_exp_details(exp_details, EXP_DETAILS)["name"]
        self.conditions = self.check_exp_details(exp_details, CONDITIONS)
        self.val = None
        self.err = None
        self.data = None
        self.mask = None

    def __len__(self):
        return self.val.size

    @abc.abstractmethod
    def print_profile(self, params=None):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def filter_points(self, params=None):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and back-calculated
        values."""
        pass

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
