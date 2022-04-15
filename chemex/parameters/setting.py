from __future__ import annotations

import re
from dataclasses import dataclass
from dataclasses import field
from re import Pattern
from typing import Literal

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.configuration.parameters import DefaultSetting
from chemex.parameters.name import ParamName
from chemex.parameters.spin_system import SpinSystem

_RE_IDS = re.compile(r"__[a-zA-Z_][a-zA-Z0-9_]*")
_RE_NAMES = re.compile(r"\{([a-zA-Z0-9_]+)\}")


@dataclass(frozen=True)
class NameSetting:
    name: str
    spin_system_part: Literal["i", "s", "is", "g", ""] = ""
    conditions_part: tuple[str, ...] = field(default_factory=tuple)

    def get_param_name(self, spin_system: SpinSystem, conditions: Conditions):
        param_spin_system = spin_system.build_sub_spin_system(self.spin_system_part)
        param_conditions = conditions.select_conditions(self.conditions_part)
        return ParamName(self.name, param_spin_system, param_conditions)


class ExpressionSetting:
    def __init__(self, re_dependencies: Pattern[str]):
        self.re_dependencies = re_dependencies
        self.__expr: str = ""
        self.__dependencies: set[str] = set()

    @property
    def expr(self) -> str:
        return self.__expr

    @expr.setter
    def expr(self, expression: str) -> None:
        self.__expr = expression
        self.__dependencies = set(self.re_dependencies.findall(expression))

    @property
    def dependencies(self) -> set[str]:
        return self.__dependencies

    def __str__(self) -> str:
        return self.__expr


class ParamLocalSetting:
    def __init__(
        self,
        name_setting: NameSetting,
        value: float | None = None,
        min: float = -np.inf,
        max: float = np.inf,
        vary: bool = False,
        expr: str = "",
    ) -> None:
        self.__expr: ExpressionSetting = ExpressionSetting(_RE_NAMES)
        self.name_setting = name_setting
        self.value = value
        self.min = min
        self.max = max
        self.vary = vary
        self.expr = expr

    @property
    def expr(self) -> str:
        return str(self.__expr)

    @expr.setter
    def expr(self, value: str) -> None:
        self.__expr.expr = value

    @property
    def dependencies(self) -> set[str]:
        return self.__expr.dependencies


class ParamSetting:
    def __init__(
        self,
        param_name: ParamName,
        value: float | None = None,
        stderr: float | None = None,
        min: float = -np.inf,
        max: float = np.inf,
        vary: bool = False,
        expr: str = "",
        brute_step: float | None = None,
    ) -> None:

        self.__expr: ExpressionSetting = ExpressionSetting(_RE_IDS)
        self.param_name = param_name
        self.value = value
        self.stderr = stderr
        self.min = min
        self.max = max
        self.vary = vary
        self.expr = expr
        self.brute_step = brute_step

    @property
    def expr(self) -> str:
        return str(self.__expr)

    @expr.setter
    def expr(self, value: str) -> None:
        self.__expr.expr = value

    @property
    def dependencies(self) -> set[str]:
        return self.__expr.dependencies

    @property
    def id(self) -> str:
        return self.param_name.id

    @property
    def args(self) -> tuple:
        return self.id, self.value, self.vary, self.min, self.max, self.expr

    def set(self, default_setting: DefaultSetting) -> None:
        self.value = default_setting.value
        if default_setting.min is not None:
            self.min = default_setting.min
        if default_setting.max is not None:
            self.max = default_setting.max
        if default_setting.min is not None:
            self.brute_step = default_setting.brute_step


LocalSettings = dict[str, ParamLocalSetting]

Parameters = dict[str, ParamSetting]
