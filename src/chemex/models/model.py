from __future__ import annotations

import sys
from dataclasses import dataclass
from string import ascii_lowercase

from chemex.messages import print_model_error
from chemex.models.factory import model_factory


@dataclass
class _Model:
    _name: str = "2st"
    _states: str = "ab"
    _model_free: bool = False
    _temp_coef: bool = False

    @staticmethod
    def validate_model_name(name: str) -> str:
        if name not in model_factory.set:
            print_model_error(name)
            sys.exit()
        return name

    def set_model(self, name: str) -> None:
        kinetic_model_name, *ext = name.split(".")
        self._name = self.validate_model_name(kinetic_model_name)
        state_nb = int(kinetic_model_name[0])
        self._states = ascii_lowercase[:state_nb]
        self._model_free = "mf" in ext if ext else False
        self._temp_coef = "tc" in ext if ext else False

    @property
    def name(self) -> str:
        return self._name

    @property
    def states(self) -> str:
        return self._states

    @property
    def model_free(self) -> bool:
        return self._model_free

    @property
    def temp_coef(self) -> bool:
        return self._temp_coef


model = _Model()

set_model = model.set_model
