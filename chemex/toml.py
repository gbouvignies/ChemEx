from __future__ import annotations

import sys
from collections.abc import Iterable
from collections.abc import MutableMapping
from pathlib import Path
from typing import Any

import tomli as tomllib  # Will be part of standard library in Python 3.11

from chemex.messages import print_file_not_found
from chemex.messages import print_toml_error


def read_toml(filename: Path) -> MutableMapping[str, Any]:
    """Read and parse the experiment configuration file with 'toml."""

    with open(filename, "rb") as file:
        try:
            config = tomllib.load(file)
        except FileNotFoundError:
            print_file_not_found(filename)
            sys.exit(1)
        except (tomllib.TOMLDecodeError, TypeError) as error:
            print_toml_error(filename, error)
            sys.exit(1)

        else:
            return config

    # Recursive dictionary merge
    # Copyright (C) 2016 Paul Durivage <pauldurivage+github@gmail.com>
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.


def _dict_merge(
    dct: MutableMapping[str, Any], merge_dct: MutableMapping[str, Any]
) -> None:
    """Recursive dictionary merge.

    Inspired by :meth:``dict.update()``, instead of updating only top-level
    keys, dict_merge recurses down into dicts nested to an arbitrary depth,
    updating keys. The ``merge_dct`` is merged into ``dct``.

    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None

    """
    for key in merge_dct:
        if (
            key in dct
            and isinstance(dct[key], dict)
            and isinstance(merge_dct[key], MutableMapping)
        ):
            _dict_merge(dct[key], merge_dct[key])
        else:
            dct[key] = merge_dct[key]


def read_toml_multi(filenames: Iterable[Path]) -> MutableMapping[str, Any]:
    """Read and parse multiple experiment configuration files with 'toml'."""
    config: dict[str, Any] = {}
    for filename in filenames:
        _dict_merge(config, read_toml(filename))
    return config


def normalize_path(working_dir: Path, filename: Path) -> Path:
    """Normalize the path of a filename relative to a specific directory."""
    path = filename
    if not path.is_absolute():
        path = working_dir / path
    return path.resolve()
