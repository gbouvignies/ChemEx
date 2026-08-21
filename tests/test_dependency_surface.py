from __future__ import annotations

import re
import subprocess
import sys
import textwrap
import tomllib
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
REMOVED_RUNTIME_PACKAGES = {"asteval", "lmfit", "numdifftools"}


def _requirement_name(requirement: str) -> str:
    match = re.match(r"[A-Za-z0-9_.-]+", requirement)
    assert match is not None
    return match.group().lower()


def test_root_runtime_dependencies_exclude_removed_legacy_stack() -> None:
    project = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    dependency_names = {
        _requirement_name(requirement)
        for requirement in project["project"]["dependencies"]
    }

    assert dependency_names.isdisjoint(REMOVED_RUNTIME_PACKAGES)


def test_application_imports_when_removed_legacy_stack_is_unavailable() -> None:
    script = textwrap.dedent(
        """
        import importlib.abc
        import sys

        blocked = {"asteval", "lmfit", "numdifftools"}

        class BlockRemovedPackages(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.partition(".")[0] in blocked:
                    raise ModuleNotFoundError(fullname)
                return None

        sys.meta_path.insert(0, BlockRemovedPackages())

        import chemex.chemex
        from chemex.runtime import AnalysisSession

        AnalysisSession.create()
        """
    )
    result = subprocess.run(  # noqa: S603 - fixed interpreter and inline test script
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
