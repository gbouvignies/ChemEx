"""Deletion gate for Method interpretation outside compilation."""

from __future__ import annotations

import ast
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "src/chemex"
METHOD_LANGUAGE = {
    "MethodPlan",
    "StepPlan",
    "ProfileSelection",
    "ParameterSelector",
    "FitAction",
    "FixAction",
    "ConstrainAction",
    "GridSearch",
    "DeSearch",
}
OLD_INTERPRETERS = {
    "effective_role_actions",
    "compile_parameterization_from_actions",
    "select_profiles",
    "resolve_grid_axes",
    "resolve_de_coordinates",
}


def _tree(relative_path: str) -> ast.Module:
    return ast.parse((SOURCE / relative_path).read_text(encoding="utf-8"))


def _imported_names(tree: ast.Module) -> set[str]:
    direct = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
    }
    modules = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
        if alias.name.startswith("chemex.configuration.method_")
    }
    return direct | modules


def _called_names(tree: ast.Module) -> set[str]:
    return {
        node.func.attr if isinstance(node.func, ast.Attribute) else node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, (ast.Attribute, ast.Name))
    }


def test_execution_and_numerical_fit_cannot_interpret_method_language() -> None:
    # Only the input façade and compiler may import Method language in optimize.
    for path in (SOURCE / "optimize").rglob("*.py"):
        if path.name in {"fitting.py", "method_compiler.py"}:
            continue
        module = str(path.relative_to(SOURCE))
        tree = _tree(module)
        assert not METHOD_LANGUAGE & _imported_names(tree), module
        assert not {
            name
            for name in _imported_names(tree)
            if name.startswith("chemex.configuration.method_")
        }, module
        assert not OLD_INTERPRETERS & _called_names(tree), module

    # Any new production Method-language importer needs an explicit boundary
    # decision, rather than silently becoming another execution interpreter.
    allowed = {
        Path("chemex.py"),
        Path("optimize/fitting.py"),
        Path("optimize/method_compiler.py"),
        Path("runtime/session.py"),
        Path("parameters/parameterization.py"),
        Path("containers/experiment.py"),
        Path("containers/experiments.py"),
    }
    for path in SOURCE.rglob("*.py"):
        module = path.relative_to(SOURCE)
        if module.parts[0] == "configuration" or module in allowed:
            continue
        assert not METHOD_LANGUAGE & _imported_names(_tree(str(module))), module


def test_compilation_has_no_analysis_values_dependency_or_action_replay() -> None:
    compiler = _tree("optimize/method_compiler.py")
    assert not {
        node.attr for node in ast.walk(compiler) if isinstance(node, ast.Attribute)
    } & {"analysis_values", "effective_value", "snapshot"}
    assert "effective_role_actions" not in _called_names(compiler)

    parameterization = _tree("parameters/parameterization.py")
    obsolete_rules = {"_build_action_rules", "_role_for", "_match_selector"}
    assert not obsolete_rules & {
        node.name
        for node in ast.walk(parameterization)
        if isinstance(node, ast.FunctionDef)
    }
