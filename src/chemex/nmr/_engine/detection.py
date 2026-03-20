from __future__ import annotations

import re
from collections.abc import Mapping

import numpy as np

from chemex.typing import Array

_FIRST_TERM = re.compile(r"\s*([+-]?)\s*\[([^\[\]\s]+)\]\s*")
_NEXT_TERM = re.compile(r"([+-])\s*\[([^\[\]\s]+)\]\s*")
_SIGNS = {"": 1.0, "+": 1.0, "-": -1.0}


def _get_component_vector(name: str, vectors: Mapping[str, Array]) -> Array:
    try:
        return vectors[name]
    except KeyError as error:
        msg = f"Unknown detection component: {name!r}"
        raise ValueError(msg) from error


def build_detection_vector(expression: str, vectors: Mapping[str, Array]) -> Array:
    """Build a detection vector from a simple linear expression.

    Supported syntax is a sum/subtraction of bracketed basis components, such as
    ``[iz_a]`` or ``[2izsz_a] - [iz_a]``.
    """
    expression_stripped = expression.strip()
    if not expression_stripped:
        msg = "Detection expression cannot be empty"
        raise ValueError(msg)

    try:
        sample_vector = next(iter(vectors.values()))
    except StopIteration as error:
        msg = "Cannot build detection vector without basis vectors"
        raise ValueError(msg) from error

    match = _FIRST_TERM.match(expression_stripped)
    if match is None:
        msg = f"Invalid detection expression: {expression!r}"
        raise ValueError(msg)

    sign_token, component_name = match.groups()
    sign = _SIGNS[sign_token]
    detection_vector = sign * _get_component_vector(component_name, vectors)
    position = match.end()

    while position < len(expression_stripped):
        match = _NEXT_TERM.match(expression_stripped, position)
        if match is None:
            msg = f"Invalid detection expression: {expression!r}"
            raise ValueError(msg)

        operator_token, component_name = match.groups()
        sign = _SIGNS[operator_token]
        detection_vector = detection_vector + sign * _get_component_vector(
            component_name,
            vectors,
        )
        position = match.end()

    return np.array(detection_vector, copy=True, dtype=sample_vector.dtype)
