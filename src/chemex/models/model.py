from __future__ import annotations

from dataclasses import dataclass, field
from string import ascii_lowercase

from chemex.exceptions import ChemExError

SUPPORTED_MODEL_EXTENSIONS = ("mf", "rs", "tc")


class ModelSelectionError(ChemExError, ValueError):
    """The requested kinetics model or extension is not supported."""

    def __init__(
        self,
        name: str,
        explanation: str,
        *,
        available: tuple[str, ...],
    ) -> None:
        super().__init__(explanation)
        self.name = name
        self.explanation = explanation
        self.available = available


@dataclass(frozen=True, order=True)
class ModelSpec:
    name: str = "2st"
    states: str = "ab"
    model_free: bool = False
    temp_coef: bool = False
    residue_specific: bool = False

    @property
    def identity(self) -> str:
        """Return the exact semantic identity of the selected model variant."""
        return (
            f"{self.name}|states={self.states}|model_free={self.model_free}|"
            f"temp_coef={self.temp_coef}|residue_specific={self.residue_specific}"
        )

    @staticmethod
    def validate_model_name(name: str) -> str:
        from chemex.models.factory import model_factory

        if name not in model_factory.set:
            raise ModelSelectionError(
                name,
                f"The model {name!r} is not available.",
                available=tuple(sorted(model_factory.set)),
            )
        return name

    @staticmethod
    def validate_extensions(name: str, extensions: list[str]) -> set[str]:
        unknown_suffixes = set(extensions) - set(SUPPORTED_MODEL_EXTENSIONS)
        if unknown_suffixes:
            suffixes = ", ".join(f".{suffix}" for suffix in sorted(unknown_suffixes))
            raise ModelSelectionError(
                name,
                f"The model {name!r} uses unsupported suffixes: {suffixes}.",
                available=tuple(f".{suffix}" for suffix in SUPPORTED_MODEL_EXTENSIONS),
            )
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

    @property
    def identity(self) -> str:
        return self._spec.identity
