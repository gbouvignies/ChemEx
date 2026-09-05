import tomllib
from collections.abc import Iterable, MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemex.exceptions import ChemExError


@dataclass(eq=False)
class TomlReadError(ChemExError):
    """A user-supplied TOML file could not be read or parsed."""

    filename: Path
    error: OSError | UnicodeDecodeError | tomllib.TOMLDecodeError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


def _deep_update(target: dict, src: dict) -> dict:
    """Recursively update target dict with src dict (deep merge)."""
    for key, value in src.items():
        if key in target and isinstance(target[key], dict) and isinstance(value, dict):
            target[key] = _deep_update(target[key], value)
        else:
            target[key] = value
    return target


def load_toml(filename: Path) -> dict[str, Any]:
    """Load a TOML file without applying command-line error handling."""
    with filename.open(mode="rb") as file:
        return tomllib.load(file)


def read_toml(filename: Path) -> dict[str, Any]:
    """Read and parse the experiment configuration file with 'toml."""
    try:
        config = load_toml(filename)
    except OSError as error:
        raise TomlReadError(filename, error) from error
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise TomlReadError(filename, error) from error

    return config


def read_toml_multi(filenames: Iterable[Path]) -> MutableMapping[str, Any]:
    """Read and parse multiple experiment configuration files with 'toml'."""
    config: dict[str, Any] = {}
    for filename in filenames:
        config = _deep_update(config, read_toml(filename))
    return config


def normalize_path(working_dir: Path, filename: Path) -> Path:
    """Normalize the path of a filename relative to a specific directory."""
    path = filename
    if not path.is_absolute():
        path = working_dir / path
    return path.resolve()
