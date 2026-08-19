from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

_LAUNCHER = Path(__file__).parents[1] / "src/chemex/numerical_lanes/launcher.py"


def _wheel(
    path: Path,
    marker: str,
    *,
    preload_attestor: bool = False,
    path_hook: bool = False,
) -> Path:
    files = {
        "chemex/__init__.py": (
            "from . import numerical_lanes\n" if preload_attestor else ""
        ),
        "chemex/application_probe.py": f"MARKER = {marker!r}\n",
        "chemex/numerical_lanes/__init__.py": (
            "raise RuntimeError('application wheel attestor was imported')\n"
        ),
        "chemex-1.0.dist-info/METADATA": (
            "Metadata-Version: 2.1\nName: chemex\nVersion: 1.0\n"
        ),
        "chemex-1.0.dist-info/WHEEL": (
            "Wheel-Version: 1.0\n"
            "Generator: test_numerical_lane_launcher\n"
            "Root-Is-Purelib: true\n"
            "Tag: py3-none-any\n"
        ),
    }
    if path_hook:
        files["application-overlay.pth"] = "/untrusted\n"
    rows: list[tuple[str, str, str]] = []
    for name, content in files.items():
        data = content.encode("utf-8")
        digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=")
        rows.append((name, f"sha256={digest.decode('ascii')}", str(len(data))))
    record_name = "chemex-1.0.dist-info/RECORD"
    rows.append((record_name, "", ""))
    output = io.StringIO(newline="")
    csv.writer(output, lineterminator="\n").writerows(rows)
    files[record_name] = output.getvalue()
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _qualification_module(root: Path) -> None:
    module_root = root / "tests/qualification"
    module_root.mkdir(parents=True)
    (root / "tests/__init__.py").write_text("", encoding="ascii")
    (module_root / "probe.py").write_text(
        """
import json
import os
from pathlib import Path

from chemex.application_probe import MARKER
import chemex.numerical_lanes as numerical_lanes

print(json.dumps({
    "attestor": str(Path(numerical_lanes.__file__).resolve()),
    "implementation": os.environ["CHEMEX_IMPLEMENTATION_WHEEL_SHA256"],
    "marker": MARKER,
    "pythonpath": os.environ.get("PYTHONPATH"),
}, sort_keys=True))
""".lstrip(),
        encoding="ascii",
    )


def _run_launcher(
    wheel: Path, expected_sha256: str, qualification_root: Path, apps: Path
) -> subprocess.CompletedProcess[str]:
    overlay = qualification_root / "untrusted-overlay"
    overlay.mkdir(exist_ok=True)
    (overlay / "sitecustomize.py").write_text(
        f"from pathlib import Path\nPath({str(qualification_root / 'ambient-import-ran')!r}).touch()\n",
        encoding="ascii",
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(overlay)
    return subprocess.run(  # noqa: S603 - the repository interpreter is fixed
        [
            sys.executable,
            "-I",
            str(_LAUNCHER),
            "--implementation-wheel",
            str(wheel),
            "--implementation-sha256",
            expected_sha256,
            "--application-directory",
            str(apps),
            "--qualification-root",
            str(qualification_root),
            "--module",
            "tests.qualification.probe",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )


def test_launcher_separates_implementation_identity_from_lane_attestor(
    tmp_path: Path,
) -> None:
    (tmp_path / "wheel-a").mkdir()
    (tmp_path / "wheel-b").mkdir()
    wheel_a = _wheel(tmp_path / "wheel-a/chemex-1.0-py3-none-any.whl", "A")
    wheel_b = _wheel(tmp_path / "wheel-b/chemex-1.0-py3-none-any.whl", "B")
    qualification_root = tmp_path / "qualification"
    _qualification_module(qualification_root)

    result_a = _run_launcher(
        wheel_a, _sha256(wheel_a), qualification_root, tmp_path / "apps-a"
    )
    mismatch = _run_launcher(
        wheel_b, _sha256(wheel_a), qualification_root, tmp_path / "apps-mismatch"
    )
    result_b = _run_launcher(
        wheel_b, _sha256(wheel_b), qualification_root, tmp_path / "apps-b"
    )

    assert result_a.returncode == 0, result_a.stderr
    assert mismatch.returncode != 0
    assert "implementation wheel SHA-256 does not match" in mismatch.stderr
    assert result_b.returncode == 0, result_b.stderr
    record_a = json.loads(result_a.stdout)
    record_b = json.loads(result_b.stdout)
    assert record_a["marker"] == "A"
    assert record_b["marker"] == "B"
    assert record_a["implementation"] != record_b["implementation"]
    assert record_a["pythonpath"] is None
    assert record_b["pythonpath"] is None
    assert not (qualification_root / "ambient-import-ran").exists()
    assert record_a["attestor"] == record_b["attestor"]
    assert record_a["attestor"] == str((_LAUNCHER.parent / "__init__.py").resolve())


def test_launcher_rejects_application_preloading_its_attestor(tmp_path: Path) -> None:
    (tmp_path / "wheel").mkdir()
    wheel = _wheel(
        tmp_path / "wheel/chemex-1.0-py3-none-any.whl",
        "shadow",
        preload_attestor=True,
    )
    qualification_root = tmp_path / "qualification"
    _qualification_module(qualification_root)

    result = _run_launcher(wheel, _sha256(wheel), qualification_root, tmp_path / "apps")

    assert result.returncode != 0
    assert "application wheel attestor was imported" in result.stderr


def test_launcher_requires_isolation_and_rejects_path_hooks(tmp_path: Path) -> None:
    without_isolation = subprocess.run(  # noqa: S603
        [sys.executable, str(_LAUNCHER)],
        check=False,
        capture_output=True,
        text=True,
    )
    (tmp_path / "wheel").mkdir()
    wheel = _wheel(
        tmp_path / "wheel/chemex-1.0-py3-none-any.whl",
        "path-hook",
        path_hook=True,
    )
    qualification_root = tmp_path / "qualification"
    _qualification_module(qualification_root)

    path_hook = _run_launcher(
        wheel, _sha256(wheel), qualification_root, tmp_path / "apps"
    )

    assert without_isolation.returncode != 0
    assert "Python isolated mode (-I) is required" in without_isolation.stderr
    assert path_hook.returncode != 0
    assert "may not install Python path hooks" in path_hook.stderr
