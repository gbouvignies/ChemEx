"""TODO: module docstring."""
import abc
import copy
from functools import lru_cache

import numpy as np
from numpy.lib import recfunctions as rfn

from chemex import peaks
from chemex.spindynamics import basis
from chemex.spindynamics import default
from chemex.spindynamics import util


class BaseProfile(metaclass=abc.ABCMeta):
    """TODO: class docstring."""

    EXP_DETAILS = {"name": {"type": str}}
    CONDITIONS = {
        "h_larmor_frq": {"type": float},
        "temperature": {"type": float},
        "p_total": {"default": None, "type": float},
        "l_total": {"default": None, "type": float},
    }
    SPIN_SYSTEM = None
    CONSTRAINTS = None
    EQUILIBRIUM = False
    DTYPE = [("par", "f8"), ("intensity", "f8"), ("error", "f8")]

    def __init__(self, name=None, data=None, exp_details=None, model=None):
        """TODO: method docstring."""

        if name is None:
            name = ""

        self.conditions = self.check_exp_details(exp_details, self.CONDITIONS)
        self.exp_details = self.check_exp_details(exp_details, self.EXP_DETAILS)

        self.name = name
        self.peak = peaks.Peak(self.name)
        self.model = util.parse_model(model)
        self.experiment_name = self.exp_details["name"]
        self.data = data
        self.mask = np.ones(len(self.data), dtype=np.bool)

        # Set the Liouvillian
        self.liouv = basis.Liouvillian(
            system=self.SPIN_SYSTEM,
            state_nb=self.model.state_nb,
            atoms=self.peak.atoms,
            h_larmor_frq=self.conditions["h_larmor_frq"],
            equilibrium=self.EQUILIBRIUM,
        )

        self.ppms_i = self.liouv.ppms["i"]

        # Set the parameters this profile depends on
        self.map_names, self.params = default.create_params(
            basis=self.liouv,
            model=self.model,
            nuclei=self.peak.names,
            conditions=self.conditions,
            constraints=self.CONSTRAINTS,
        )

        self.calculate_unscaled_profile = lru_cache(256)(
            self._calculate_unscaled_profile
        )

    def __len__(self):
        return self.data.size

    @property
    @abc.abstractmethod
    def reference(self):
        pass

    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and back-calculated
        values."""

        values = self.calculate_profile(params)
        residuals = (self.data["intensity"] - values) / self.data["error"]

        return residuals[self.mask]

    def calculate_profile(self, params=None, **kwargs):
        """Calculate the CEST profile."""
        params_local = tuple(
            (name_s, params[name_l].value) for name_s, name_l in self.map_names.items()
        )

        values = self.calculate_unscaled_profile(params_local)

        try:
            scale = sum(
                values * self.data["intensity"] / self.data["error"] ** 2
            ) / sum((values / self.data["error"]) ** 2)

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
        profile.data["intensity"] = (
            self.calculate_profile(params)
            + np.random.randn(len(self.data["intensity"])) * self.data["error"]
        )

        return profile

    def make_bs_profile(self):
        """Make a profile for boostrap analysis."""

        indexes = np.array(range(len(self.data["intensity"])))
        pool1 = indexes[self.reference]
        pool2 = indexes[~self.reference]

        bs_indexes = []
        if pool1.size:
            bs_indexes.extend(np.random.choice(pool1, len(pool1)))
        bs_indexes.extend(np.random.choice(pool2, len(pool2)))

        bs_indexes = sorted(bs_indexes)

        profile = copy.copy(self)
        profile.data = self.data[bs_indexes]
        profile.calculate_unscaled_profile = lru_cache(256)(
            profile._calculate_unscaled_profile
        )

        return profile

    def normalize_profile(self, params=None):

        ndata = self.data.copy()

        scale = 1.0 / np.mean(ndata[self.reference]["intensity"])

        ndata["intensity"] *= scale
        ndata["error"] *= abs(scale)

        if params is not None:
            values = self.calculate_profile(params) * scale
            ndata = rfn.append_fields(ndata, "intensity_calc", values, usemask=False)

        return ndata, scale

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
                missing.append(f"'{name}'")
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

    @staticmethod
    def get_bool(value):

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

        if value.lower() not in BOOLEAN_STATES:
            raise ValueError("Not a boolean: %s" % value)

        return BOOLEAN_STATES[value.lower()]
