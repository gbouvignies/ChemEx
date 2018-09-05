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
    def print_profile(self):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def calculate_profile(self, params):
        """TODO: method docstring."""
        pass

    def create_default_parameters(self):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def filter_points(self, params=None):
        """TODO: method docstring."""
        pass

    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and back-calculated
        values."""

        values = self.calculate_profile(params)
        residuals = (self.val - values) / self.err

        return residuals[self.mask]

    @staticmethod
    def check_par(
        parameters=None, name=None, convert=None, default=None, required=True
    ):
        """Check for experimental parameters and convert them to their appropriate
        type."""
        value = parameters.get(name, default)

        if isinstance(value, str):
            value = value.split()
            if len(value) == 1:
                value = value[0]

        if required and value is None:
            exit(
                "Missing experimental parameter detected. Please set: {:s}".format(name)
            )

        if convert is not None and value is not None:
            try:
                if isinstance(value, list):
                    value = [convert(_) for _ in value]
                else:
                    value = convert(value)
            except ValueError:
                exit(
                    "\nExperimental parameter of wrong type detected. Please make"
                    " sure that {:s} is a {:s}".format(name, str(convert))
                )

        return value

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
