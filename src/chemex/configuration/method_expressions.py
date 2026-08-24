from __future__ import annotations

import ast
import io
import math
import re
import tokenize
from collections.abc import Callable
from typing import cast

from chemex.configuration.method_plan import (
    BinaryExpression,
    Constraint,
    ConstraintExpression,
    DeCoordinate,
    DeRange,
    GridAxis,
    GridRange,
    GridValues,
    LiteralExpression,
    MethodFormatError,
    ParameterSelector,
    SearchScale,
    SelectorExpression,
    SourceRef,
    UnaryExpression,
)
from chemex.parameters.name import ParamName
from chemex.parameters.spin_system import SpinSystem

_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_QUALIFIER = re.compile(
    r"(?P<key>NUC|T|B0|\[P\]|\[L\]|D2O)\s*->\s*(?P<value>.+)\Z",
    re.IGNORECASE,
)
_NUMERIC_QUALIFIER = {
    "T": re.compile(rf"(?P<number>{_NUMBER})(?P<unit>C)?\Z", re.IGNORECASE),
    "B0": re.compile(
        rf"(?P<number>{_NUMBER})(?P<unit>MHz)?\Z",
        re.IGNORECASE,
    ),
    "[P]": re.compile(rf"(?P<number>{_NUMBER})(?P<unit>M)?\Z", re.IGNORECASE),
    "[L]": re.compile(rf"(?P<number>{_NUMBER})(?P<unit>M)?\Z", re.IGNORECASE),
    "D2O": re.compile(rf"(?P<number>{_NUMBER})\Z", re.IGNORECASE),
}
_PUBLIC_OPERATORS = {"+", "-", "*", "/", "(", ")"}
_FIELD_FOR_QUALIFIER = {
    "T": "temperature",
    "B0": "h_larmor_frq",
    "[P]": "p_total",
    "[L]": "l_total",
    "D2O": "d2o",
}
type _SelectorValues = dict[str, str | float | None]


def parse_legacy_selector(text: str, source: SourceRef) -> ParameterSelector:
    parsed = ParamName.from_section(text)
    conditions = parsed.conditions
    return ParameterSelector(
        parsed.name,
        str(parsed.spin_system) or None,
        conditions.temperature,
        conditions.h_larmor_frq,
        conditions.p_total,
        conditions.l_total,
        conditions.d2o,
        source,
    )


def _error(message: str, source: SourceRef, start: int, end: int) -> MethodFormatError:
    return MethodFormatError(message, _span_source(source, start, end))


def _span_source(source: SourceRef, start: int, end: int) -> SourceRef:
    offset = source.start or 0
    return SourceRef(
        source.filename,
        source.step,
        source.field,
        source.index,
        offset + start,
        offset + end,
    )


def _parse_nucleus(
    value_text: str,
    source: SourceRef,
    start: int,
    end: int,
) -> str:
    try:
        spin_system = SpinSystem.from_name_strict(value_text)
    except (TypeError, ValueError) as error:
        raise _error(
            f"Invalid NUC spin-system selector {value_text!r}", source, start, end
        ) from error
    return str(spin_system)


def _parse_numeric_qualifier(
    key: str,
    value_text: str,
    source: SourceRef,
    start: int,
    end: int,
) -> float:
    numeric = _NUMERIC_QUALIFIER[key].fullmatch(value_text)
    if numeric is None:
        raise _error(
            f"Invalid value or unit for selector qualifier {key}", source, start, end
        )
    value = float(numeric.group("number"))
    if not math.isfinite(value):
        raise _error(f"Selector qualifier {key} must be finite", source, start, end)
    if key in {"B0", "[P]", "[L]"} and value < 0.0:
        raise _error(
            f"Selector qualifier {key} must be nonnegative", source, start, end
        )
    if key == "D2O" and not 0.0 < value < 1.0:
        raise _error(
            "Selector qualifier D2O must be between 0 and 1", source, start, end
        )
    return round(value, 1) if key in {"T", "B0"} else value


def _parse_qualifier(
    raw_fragment: str,
    start: int,
    source: SourceRef,
    seen: set[str],
    values: _SelectorValues,
) -> None:
    fragment = raw_fragment.strip()
    content_start = start + len(raw_fragment) - len(raw_fragment.lstrip())
    end = content_start + len(fragment)
    if not fragment:
        raise _error(
            "Empty selector qualifier", source, start, start + len(raw_fragment)
        )
    match = _QUALIFIER.fullmatch(fragment)
    if match is None:
        raise _error(
            f"Unknown or malformed selector qualifier {fragment!r}",
            source,
            content_start,
            end,
        )
    key = match.group("key").upper()
    if key in seen:
        raise _error(f"Duplicate selector qualifier {key}", source, content_start, end)
    seen.add(key)
    value_text = match.group("value").strip()
    if key == "NUC":
        values["spin_system"] = _parse_nucleus(value_text, source, content_start, end)
    else:
        values[_FIELD_FOR_QUALIFIER[key]] = _parse_numeric_qualifier(
            key, value_text, source, content_start, end
        )


def parse_strict_selector(text: str, source: SourceRef) -> ParameterSelector:
    stripped = text.strip()
    leading = len(text) - len(text.lstrip())
    if not stripped:
        raise _error("Parameter selector cannot be empty", source, 0, len(text))
    fragments = stripped.split(",")
    name = fragments[0].strip()
    if not _NAME.fullmatch(name):
        raise _error(
            f"Invalid parameter-name token {name!r}",
            source,
            leading,
            leading + len(fragments[0]),
        )
    values: _SelectorValues = dict.fromkeys(
        ("spin_system", "temperature", "h_larmor_frq", "p_total", "l_total", "d2o")
    )
    cursor = len(fragments[0])
    seen: set[str] = set()
    for raw_fragment in fragments[1:]:
        _parse_qualifier(raw_fragment, leading + cursor + 1, source, seen, values)
        cursor += len(raw_fragment) + 1
    source_offset = source.start or 0
    selector_source = SourceRef(
        source.filename,
        source.step,
        source.field,
        source.index,
        source_offset + leading,
        source_offset + leading + len(stripped),
    )
    return ParameterSelector(
        name.upper(),
        cast("str | None", values["spin_system"]),
        cast("float | None", values["temperature"]),
        cast("float | None", values["h_larmor_frq"]),
        cast("float | None", values["p_total"]),
        cast("float | None", values["l_total"]),
        cast("float | None", values["d2o"]),
        selector_source,
    )


def _scan_bracketed_selectors(
    text: str, source: SourceRef
) -> tuple[tuple[int, int, str], ...]:
    selectors: list[tuple[int, int, str]] = []
    position = 0
    while position < len(text):
        if text[position] != "[":
            position += 1
            continue
        start = position
        depth = 1
        position += 1
        while position < len(text) and depth:
            if text[position] == "[":
                depth += 1
            elif text[position] == "]":
                depth -= 1
            position += 1
        if depth:
            raise _error("Unclosed parameter selector", source, start, len(text))
        selectors.append((start, position, text[start + 1 : position - 1]))
    return tuple(selectors)


def _compile_constraint_ast(
    node: ast.AST,
    references: dict[str, ParameterSelector],
    source: SourceRef,
) -> ConstraintExpression:
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        try:
            value = float(node.value)
        except (OverflowError, ValueError):
            value = math.inf
        if math.isfinite(value):
            return LiteralExpression(value)
    if isinstance(node, ast.Name) and node.id in references:
        return SelectorExpression(references[node.id])
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        operator = "+" if isinstance(node.op, ast.UAdd) else "-"
        return UnaryExpression(
            operator,
            _compile_constraint_ast(node.operand, references, source),
        )
    if isinstance(node, ast.BinOp) and isinstance(
        node.op,
        (ast.Add, ast.Sub, ast.Mult, ast.Div),
    ):
        if isinstance(node.op, ast.Add):
            operator = "+"
        elif isinstance(node.op, ast.Sub):
            operator = "-"
        elif isinstance(node.op, ast.Mult):
            operator = "*"
        else:
            operator = "/"
        return BinaryExpression(
            operator,
            _compile_constraint_ast(node.left, references, source),
            _compile_constraint_ast(node.right, references, source),
        )
    raise MethodFormatError(
        f"Constraint contains unsupported scalar syntax {type(node).__name__}",
        source,
    )


def _constraint_target(
    left: str,
    source: SourceRef,
    selector_parser: Callable[[str, SourceRef], ParameterSelector],
) -> ParameterSelector:
    selectors = _scan_bracketed_selectors(left, source)
    if (
        len(selectors) != 1
        or left[: selectors[0][0]].strip()
        or left[selectors[0][1] :].strip()
    ):
        raise MethodFormatError(
            "Constraint left side must be one bracketed selector",
            _span_source(source, 0, len(left)),
        )
    start, end, text = selectors[0]
    return selector_parser(
        text,
        SourceRef(
            source.filename,
            source.step,
            source.field,
            source.index,
            start + 1,
            end - 1,
        ),
    )


def _constraint_python_source(
    left: str,
    right: str,
    source: SourceRef,
    selector_parser: Callable[[str, SourceRef], ParameterSelector],
) -> tuple[str, dict[str, ParameterSelector]]:
    references: dict[str, ParameterSelector] = {}
    chunks: list[str] = []
    cursor = 0
    right_source = _span_source(source, len(left) + 1, len(left) + 1 + len(right))
    for index, (start, end, selector_text) in enumerate(
        _scan_bracketed_selectors(right, right_source)
    ):
        token = f"__ref_{index}"
        chunks.extend((right[cursor:start], token))
        references[token] = selector_parser(
            selector_text,
            SourceRef(
                source.filename,
                source.step,
                source.field,
                source.index,
                len(left) + start + 2,
                len(left) + end,
            ),
        )
        cursor = end
    chunks.append(right[cursor:])
    return "".join(chunks).strip(), references


def _parse_public_scalar_expression(
    python_source: str,
    references: dict[str, ParameterSelector],
    source: SourceRef,
) -> ConstraintExpression:
    try:
        for token in tokenize.generate_tokens(io.StringIO(python_source).readline):
            ignored = token.type in {
                tokenize.NEWLINE,
                tokenize.NL,
                tokenize.ENDMARKER,
                tokenize.ENCODING,
            }
            valid_number = token.type == tokenize.NUMBER and re.fullmatch(
                _NUMBER, token.string
            )
            valid_reference = token.type == tokenize.NAME and token.string in references
            valid_operator = (
                token.type == tokenize.OP and token.string in _PUBLIC_OPERATORS
            )
            if ignored or valid_number or valid_reference or valid_operator:
                continue
            raise MethodFormatError(
                f"Constraint contains unsupported token {token.string!r}",
                source,
            )
        parsed = ast.parse(python_source, mode="eval")
    except (IndentationError, SyntaxError, tokenize.TokenError, ValueError) as error:
        raise MethodFormatError(
            "Constraint is not a valid scalar expression", source
        ) from error
    return _compile_constraint_ast(parsed.body, references, source)


def parse_constraint(
    text: str,
    source: SourceRef,
    selector_parser: Callable[[str, SourceRef], ParameterSelector],
) -> Constraint:
    entry_source = _span_source(source, 0, len(text))
    if text.count("=") != 1:
        raise MethodFormatError("Constraint must contain exactly one '='", entry_source)
    left, right = text.split("=", maxsplit=1)
    if not right.strip():
        raise MethodFormatError("Constraint right side cannot be empty", entry_source)
    target = _constraint_target(left, source, selector_parser)
    python_source, references = _constraint_python_source(
        left, right, source, selector_parser
    )
    right_start = len(left) + 1 + len(right) - len(right.lstrip())
    expression_source = _span_source(source, right_start, len(text))
    return Constraint(
        target,
        _parse_public_scalar_expression(python_source, references, expression_source),
        entry_source,
    )


def _search_parts(
    text: str, source: SourceRef
) -> tuple[str, SourceRef, str, list[str]]:
    selectors = _scan_bracketed_selectors(text, source)
    if len(selectors) != 1 or selectors[0][0] != 0:
        raise MethodFormatError(
            "V2 search entries must start with one bracketed selector",
            source,
        )
    start, end, selector_text = selectors[0]
    match = re.fullmatch(
        r"=\s*(?P<function>[A-Za-z_][A-Za-z0-9_]*)\s*\((?P<arguments>.*)\)",
        text[end:].strip(),
    )
    if match is None:
        raise MethodFormatError("Malformed search expression", source)
    arguments = [part.strip() for part in match.group("arguments").split(",")]
    if not arguments or any(not argument for argument in arguments):
        raise MethodFormatError("Search arguments cannot be empty", source)
    selector_source = SourceRef(
        source.filename,
        source.step,
        source.field,
        source.index,
        (source.start or 0) + start + 1,
        (source.start or 0) + end - 1,
    )
    return selector_text, selector_source, match.group("function").lower(), arguments


def _search_numbers(arguments: list[str], source: SourceRef) -> tuple[float, ...]:
    if any(re.fullmatch(_NUMBER, argument) is None for argument in arguments):
        raise MethodFormatError(
            "Search arguments must be finite decimal literals",
            source,
        )
    try:
        values = tuple(float(argument) for argument in arguments)
    except (OverflowError, ValueError) as error:
        raise MethodFormatError("Search arguments must be finite", source) from error
    if not all(math.isfinite(value) for value in values):
        raise MethodFormatError("Search arguments must be finite", source)
    return values


def parse_strict_grid_axis(text: str, source: SourceRef) -> GridAxis:
    leading = len(text) - len(text.lstrip())
    text = text.strip()
    entry_source = _span_source(source, leading, leading + len(text))
    selector_text, selector_source, function, arguments = _search_parts(
        text, entry_source
    )
    selector = parse_strict_selector(selector_text, selector_source)
    values = _search_numbers(arguments, entry_source)
    if function == "values":
        if len(set(values)) != len(values):
            raise MethodFormatError(
                "GRID values cannot contain duplicates", entry_source
            )
        return GridAxis(selector, GridValues(values), entry_source)
    if function not in {"lin", "log"} or len(arguments) != 3:
        message = (
            "GRID accepts only lin(low, high, count), log(low, high, count), "
            "or values(...)"
        )
        raise MethodFormatError(message, entry_source)
    if re.fullmatch(r"\d+", arguments[2]) is None:
        raise MethodFormatError("GRID count must be an integer", entry_source)
    low, high = values[:2]
    try:
        count = int(arguments[2])
    except ValueError as error:
        raise MethodFormatError(
            "GRID count must be an integer", entry_source
        ) from error
    if low >= high:
        raise MethodFormatError("GRID range must satisfy low < high", entry_source)
    if count < 2:
        raise MethodFormatError("GRID count must be at least 2", entry_source)
    if function == "log" and low <= 0.0:
        raise MethodFormatError(
            "Logarithmic GRID endpoints must be positive", entry_source
        )
    scale = SearchScale.LINEAR if function == "lin" else SearchScale.LOGARITHMIC
    return GridAxis(selector, GridRange(scale, low, high, count), entry_source)


def parse_strict_de_coordinate(text: str, source: SourceRef) -> DeCoordinate:
    leading = len(text) - len(text.lstrip())
    text = text.strip()
    entry_source = _span_source(source, leading, leading + len(text))
    selector_text, selector_source, function, arguments = _search_parts(
        text, entry_source
    )
    if function not in {"lin", "log"} or len(arguments) != 2:
        raise MethodFormatError(
            "DE accepts only lin(low, high) or log(low, high)",
            entry_source,
        )
    low, high = _search_numbers(arguments, entry_source)
    if low >= high:
        raise MethodFormatError("DE range must satisfy low < high", entry_source)
    if function == "log" and low <= 0.0:
        raise MethodFormatError(
            "Logarithmic DE endpoints must be positive", entry_source
        )
    scale = SearchScale.LINEAR if function == "lin" else SearchScale.LOGARITHMIC
    return DeCoordinate(
        parse_strict_selector(selector_text, selector_source),
        DeRange(scale, low, high),
        entry_source,
    )


def parse_legacy_grid_axis(text: str, source: SourceRef) -> GridAxis:
    compact = text.replace(" ", "")
    if compact.count("=") != 1:
        raise MethodFormatError(
            "Legacy GRID entry must contain exactly one '='", source
        )
    selector_text, expression = compact.split("=", maxsplit=1)
    selector = parse_legacy_selector(selector_text, source)
    function_match = re.fullmatch(
        r"(?P<function>lin|log)\((?P<arguments>.*)\)",
        expression,
    )
    if function_match is not None:
        function = function_match.group("function")
        arguments = function_match.group("arguments").split(",")
        if len(arguments) != 3 or re.fullmatch(r"\d+", arguments[2]) is None:
            raise MethodFormatError("Unsupported legacy GRID definition", source)
        low, high, _ = _search_numbers(arguments, source)
        scale = SearchScale.LINEAR if function == "lin" else SearchScale.LOGARITHMIC
        try:
            count = int(arguments[2])
        except ValueError as error:
            raise MethodFormatError(
                "Unsupported legacy GRID definition", source
            ) from error
        return GridAxis(selector, GridRange(scale, low, high, count), source)
    tuple_match = re.fullmatch(r"\((?P<arguments>.*)\)", expression)
    if tuple_match is None:
        raise MethodFormatError("Unsupported legacy GRID definition", source)
    arguments = [part.strip() for part in tuple_match.group("arguments").split(",")]
    values = _search_numbers(arguments, source)
    return GridAxis(selector, GridValues(values), source)
