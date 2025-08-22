from __future__ import annotations

import sys
import tomllib
from collections.abc import Iterable, MutableMapping
from pathlib import Path
from typing import Any

from pydantic.v1.utils import deep_update

from chemex.messages import print_file_not_found, print_toml_error


def read_toml(filename: Path) -> dict[str, Any]:
    """Read and parse the experiment configuration file with 'toml."""
    with filename.open(mode="rb") as file:
        try:
            config = tomllib.load(file)
        except FileNotFoundError:
            print_file_not_found(filename)
            sys.exit(1)
        except (tomllib.TOMLDecodeError, TypeError) as error:
            print_toml_error(filename, error)
            sys.exit(1)

        return config


def read_toml_multi(filenames: Iterable[Path]) -> MutableMapping[str, Any]:
    """Read and parse multiple experiment configuration files with 'toml'."""
    config: dict[str, Any] = {}
    for filename in filenames:
        config = deep_update(config, read_toml(filename))
    return config


def normalize_path(working_dir: Path, filename: Path) -> Path:
    """Normalize the path of a filename relative to a specific directory."""
    path = filename
    if not path.is_absolute():
        path = working_dir / path
    return path.resolve()
