from __future__ import annotations

from dataclasses import dataclass
from string import ascii_lowercase

KINETIC_MODEL_LIST = (
    "1st",
    "2st",
    "3st",
    "4st",
    "2st_rs",
    "2st_hd",
    "2st_eyring",
    "3st_eyring",
    "2st_binding",
    "4st_hd",
)


@dataclass
class _Model:
    __name: str = "2st"
    __states: str = "ab"
    __model_free: bool = False

    @staticmethod
    def validate_model_name(name: str) -> str:
        if name not in KINETIC_MODEL_LIST:
            print("Warning: The 'model' option should either be:")
            for model_name in KINETIC_MODEL_LIST:
                print(f"    - '{model_name}'")
            print("Set to the default value: '2st'.")
            return "2st"
        return name

    def set_model(self, name: str) -> None:
        kinetic_model_name, *ext = name.split(".")
        self.__name = self.validate_model_name(kinetic_model_name)
        state_nb = int(kinetic_model_name[0])
        self.__states = ascii_lowercase[:state_nb]
        self.__model_free = ext[0] == "mf" if ext else False

    @property
    def name(self) -> str:
        return self.__name

    @property
    def states(self) -> str:
        return self.__states

    @property
    def model_free(self) -> bool:
        return self.__model_free


model = _Model()

set_model = model.set_model
