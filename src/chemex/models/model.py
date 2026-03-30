from __future__ import annotations

import sys
from dataclasses import dataclass, field
from string import ascii_lowercase

from chemex.messages import print_model_error, print_model_suffix_error

SUPPORTED_MODEL_EXTENSIONS = ("mf", "rs", "tc")


@dataclass(frozen=True, order=True)
class ModelSpec:
    name: str = "2st"
    states: str = "ab"
    model_free: bool = False
    temp_coef: bool = False
    residue_specific: bool = False

    @staticmethod
    def validate_model_name(name: str) -> str:
        from chemex.models.factory import model_factory

        if name not in model_factory.set:
            print_model_error(name)
            sys.exit()
        return name

    @staticmethod
    def validate_extensions(name: str, extensions: list[str]) -> set[str]:
        unknown_suffixes = set(extensions) - set(SUPPORTED_MODEL_EXTENSIONS)
        if unknown_suffixes:
            print_model_suffix_error(
                name,
                unknown_suffixes,
                SUPPORTED_MODEL_EXTENSIONS,
            )
            sys.exit()
        return set(extensions)

    @classmethod
    def from_name(cls, name: str) -> ModelSpec:
        kinetic_model_name, *extensions = name.split(".")
        validated_name = cls.validate_model_name(kinetic_model_name)
        validated_extensions = cls.validate_extensions(name, extensions)
        state_nb = int(validated_name[0])
        return cls(
            name=validated_name,
            states=ascii_lowercase[:state_nb],
            model_free="mf" in validated_extensions,
            temp_coef="tc" in validated_extensions,
            residue_specific="rs" in validated_extensions,
        )


_DEFAULT_MODEL = ModelSpec()


@dataclass
class ModelState:
    _spec: ModelSpec = field(default_factory=ModelSpec)

    def set_model(self, name: str) -> None:
        self._spec = ModelSpec.from_name(name)

    def reset(self) -> None:
        self._spec = _DEFAULT_MODEL

    @property
    def spec(self) -> ModelSpec:
        return self._spec

    @property
    def name(self) -> str:
        return self._spec.name

    @property
    def states(self) -> str:
        return self._spec.states

    @property
    def model_free(self) -> bool:
        return self._spec.model_free

    @property
    def temp_coef(self) -> bool:
        return self._spec.temp_coef

    @property
    def residue_specific(self) -> bool:
        return self._spec.residue_specific
