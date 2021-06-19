"""The util module contains a variety of utility functions."""
import collections
import sys

import jsonschema as js
import tomlkit


def read_toml_multi(filenames=None, schema=None):
    """Read and parse the experiment configuration file with configparser."""
    config = {}
    for filename in filenames:
        dict_merge(config, read_toml(filename, schema))
    return config


def read_toml(filename=None, schema=None):
    """Read and parse the experiment configuration file with configparser."""
    try:
        config = tomlkit.parse(filename.read_text())
        config = _check_case(config)
        if schema is not None:
            Validator(schema).validate(config)
    except FileNotFoundError:
        sys.exit(f"\nERROR: The file '{filename}' is empty or does not exist!\n")
    except (TypeError, tomlkit.exceptions.TOMLKitError) as error:
        sys.exit(f"\nERROR: '{filename}': {error}\n")
    except js.ValidationError as error:
        sys.exit(f"ERROR: '{filename}': {error.message}")

    else:
        return config


def _check_case(config):
    return {
        section.lower(): {key.lower(): value for key, value in setting.items()}
        for section, setting in config.items()
    }


def validate(config, schema):
    """Validate the syntax of the configuration file."""
    try:
        Validator(schema).validate(config)
    except js.ValidationError as error:
        sys.exit(f"ERROR: {error.message}")


def normalize_path(working_dir, filename):
    """Normalize the path of a filename relative to a specific directory."""
    path = filename
    if not path.is_absolute():
        path = working_dir / path
    return path.resolve()


def header1(string_):
    """Print a formatted heading."""
    if string_:
        print("\n".join(["", "", string_, "=" * len(string_)]))


def header2(string_):
    """Print a formatted subheading."""
    if string_:
        print("\n".join(["", string_, "-" * len(string_)]))


def header3(string_):
    """Print a formatted subsubheading."""
    if string_:
        print(f"-- {string_} --\n")


# Taken from: https://python-jsonschema.readthedocs.io/en/stable/faq/
def _extend_with_default(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property_, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(property_, subschema["default"])
        yield from validate_properties(validator, properties, instance, schema)

    return js.validators.extend(validator_class, {"properties": set_defaults})


Validator = _extend_with_default(js.Draft7Validator)


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


def dict_merge(dct, merge_dct):
    """Recursive dictionary merge.

    Inspired by :meth:``dict.update()``, instead of updating only top-level
    keys, dict_merge recurses down into dicts nested to an arbitrary depth,
    updating keys. The ``merge_dct`` is merged into ``dct``.

    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None

    """
    for k, _ in merge_dct.items():
        if (
            k in dct
            and isinstance(dct[k], dict)
            and isinstance(merge_dct[k], collections.Mapping)
        ):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
