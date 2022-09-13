from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from typing import Optional
from typing import Union

from pydantic import Field
from pydantic import ValidationError
from pydantic.types import PositiveInt

from chemex.configuration.base import BaseModelLowerCase
from chemex.parameters.name import SpinSystem
from chemex.toml import read_toml

# Type definitions
AllType = Literal["*", "all", "ALL", "All", "ALl", "AlL"]
SelectionType = Union[list[SpinSystem], AllType, None]


class Statistics(BaseModelLowerCase):
    mc: Optional[PositiveInt] = None
    bs: Optional[PositiveInt] = None
    bsn: Optional[PositiveInt] = None


@dataclass
class Selection:
    include: SelectionType
    exclude: SelectionType


class Method(BaseModelLowerCase):
    class Config:
        anystr_lower = True

    fitmethod: str = "leastsq"
    include: SelectionType = None
    exclude: SelectionType = None
    fit: list[str] = Field(default_factory=list)
    fix: list[str] = Field(default_factory=list)
    constraints: list[str] = Field(default_factory=list)
    grid: list[str] = Field(default_factory=list)
    statistics: Optional[Statistics] = None

    @property
    def selection(self) -> Selection:
        return Selection(include=self.include, exclude=self.exclude)


FITMETHOD_ERROR_MESSAGE = """  - "FITMETHOD" must be one of:
        "leastsq": Levenberg-Marquardt (default)
        "least_squares": Least-Squares minimization, using Trust Region Reflective method
        "differential_evolution": differential evolution
        "brute": brute force method
        "basinhopping": basinhopping
        "ampgo": Adaptive Memory Programming for Global Optimization
        "nelder": Nelder-Mead
        "lbfgsb": L-BFGS-B
        "powell": Powell
        "cg": Conjugate-Gradient
        "newton": Newton-CG
        "cobyla": Cobyla
        "bfgs": BFGS
        "tnc": Truncated Newton
        "trust-ncg": Newton-CG trust-region
        "trust-exact": nearly exact trust-region
        "trust-krylov": Newton GLTR trust-region
        "trust-constr": trust-region for constrained optimization
        "dogleg": Dog-leg trust-region
        "slsqp": Sequential Linear Squares Programming
        "emcee": Maximum likelihood via Monte-Carlo Markov Chain
        "shgo": Simplicial Homology Global Optimization
        "dual_annealing": Dual Annealing optimization"""

INCLUDE_ERROR_MESSAGE = """  - "INCLUDE" must be a list of nuclei to be included or '*' to select all the profiles.
    Example: ['G3H', 'S4H', 'W5H'] or [3, 4, 5]"""

EXCLUDE_ERROR_MESSAGE = """  - "EXCLUDE" must be a list of nuclei to be excluded.
    Example: ['G3H', 'S4H', 'W5H'] or [3, 4, 5]"""

FIT_ERROR_MESSAGE = """  - "FIT" must be a list of parameters to be fitted.
    Example: ["PB", "KEX_AB"]"""

FIX_ERROR_MESSAGE = """  - "FIX" must be a list of parameters to be fixed.
    Example: ["PB", "KEX_AB"]"""

CONSTRAINTS_ERROR_MESSAGE = """  - "CONSTRAINTS" must be a list of expressions.
    Example: ["[R2_B] = [R2_A] / 2", "[R1_B] = [R1_A]"]"""

GRID_ERROR_MESSAGE = """  - "GRID" must be a list of parameters to be gridded with the associated grid definition.
    Example: GRID    = [
        "[KEX_AB] = log(100.0, 600.0, 10)",
        "[PB] = log(0.03, 0.15, 10)",
        "[DW_AB] = lin(0.0, 10.0, 5)",
    ]"""

STATISTICS_ERROR_MESSAGE = """  - "STATISTICS" must be a dictionary with keys 'MC', 'BS', 'BSN'
    Example: { "MC"=10 } or { "MC"=10, "BS"=10 }"""

METHOD_ERROR_MESSAGES = {
    "fitmethod": FITMETHOD_ERROR_MESSAGE,
    "include": INCLUDE_ERROR_MESSAGE,
    "exclude": EXCLUDE_ERROR_MESSAGE,
    "fit": FIT_ERROR_MESSAGE,
    "fix": FIX_ERROR_MESSAGE,
    "constraints": CONSTRAINTS_ERROR_MESSAGE,
    "grid": GRID_ERROR_MESSAGE,
    "statistics": STATISTICS_ERROR_MESSAGE,
}

Methods = dict[str, Method]


def read_methods(filenames: Iterable[Path]) -> Methods:
    methods = {}

    for filename in filenames:
        methods_dict = read_toml(filename)
        for section, settings in methods_dict.items():
            try:
                method = Method(**settings)
            except ValidationError as error:
                options = {option for err in error.errors() for option in err["loc"]}
                print(f"\nERROR found in '[{section}]' of '{filename}':")
                for option in options:
                    print(METHOD_ERROR_MESSAGES.get(str(option)))
                sys.exit()
            methods[section] = method
    return methods
