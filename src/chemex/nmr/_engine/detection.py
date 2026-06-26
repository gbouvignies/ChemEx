from __future__ import annotations

import re
from collections.abc import Iterable, Mapping

import numpy as np

from chemex.typing import Array

_FIRST_TERM = re.compile(r"\s*([+-]?)\s*\[([^\[\]\s]+)\]\s*")
_NEXT_TERM = re.compile(r"([+-])\s*\[([^\[\]\s]+)\]\s*")
_SIGNS = {"": 1.0, "+": 1.0, "-": -1.0}


def _parse_detection_expression(expression: str) -> list[tuple[float, str]]:
    expression_stripped = expression.strip()
    if not expression_stripped:
        msg = "Detection expression cannot be empty"
        raise ValueError(msg)

    match = _FIRST_TERM.match(expression_stripped)
    if match is None:
        msg = f"Invalid detection expression: {expression!r}"
        raise ValueError(msg)

    sign_token, component_name = match.groups()
    terms = [(_SIGNS[sign_token], component_name)]
    position = match.end()

    while position < len(expression_stripped):
        match = _NEXT_TERM.match(expression_stripped, position)
        if match is None:
            msg = f"Invalid detection expression: {expression!r}"
            raise ValueError(msg)

        operator_token, component_name = match.groups()
        terms.append((_SIGNS[operator_token], component_name))
        position = match.end()

    return terms


def _format_detection_expression(terms: Iterable[tuple[float, str]]) -> str:
    formatted: list[str] = []
    for index, (sign, component_name) in enumerate(terms):
        if index == 0:
            prefix = "- " if sign < 0.0 else ""
        else:
            prefix = " - " if sign < 0.0 else " + "
        formatted.append(f"{prefix}[{component_name}]")
    return "".join(formatted)


def build_state_detection_expression(
    expression: str,
    states: Iterable[str],
) -> str:
    """Expand a state-independent detection expression over selected states."""
    state_tuple = tuple(states)
    if not state_tuple:
        msg = "At least one observed state is required"
        raise ValueError(msg)

    terms = _parse_detection_expression(expression)
    state_terms = (
        (sign, f"{component_name}_{state}")
        for state in state_tuple
        for sign, component_name in terms
    )
    return _format_detection_expression(state_terms)


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
    try:
        sample_vector = next(iter(vectors.values()))
    except StopIteration as error:
        msg = "Cannot build detection vector without basis vectors"
        raise ValueError(msg) from error

    terms = _parse_detection_expression(expression)
    sign, component_name = terms[0]
    detection_vector = sign * _get_component_vector(component_name, vectors)
    for sign, component_name in terms[1:]:
        detection_vector = detection_vector + sign * _get_component_vector(
            component_name,
            vectors,
        )

    return np.array(detection_vector, copy=True, dtype=sample_vector.dtype)
