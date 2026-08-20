"""Run exact ChemEx wheels under the immutable numerical-lane attestor."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import os
import runpy
import shutil
import subprocess
import sys
import tempfile
import venv
import zipfile
from pathlib import Path, PurePosixPath
from types import ModuleType

_IMPLEMENTATION_DIGEST_ENV = "CHEMEX_IMPLEMENTATION_WHEEL_SHA256"
_INSTALLED_STAGE_ENV = "CHEMEX_NUMERICAL_LANE_INSTALLED_STAGE"
_APPLICATION_PREFIX_ENV = "CHEMEX_NUMERICAL_LANE_APPLICATION_PREFIX"
_QUALIFICATION_PREFIX = "tests.qualification."


class LaneLaunchError(RuntimeError):
    """The independently supplied application failed closed validation."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise LaneLaunchError("implementation wheel is unavailable") from error
    return digest.hexdigest()


def _verify_wheel(path: Path, expected_sha256: str) -> Path:
    wheel = path.resolve()
    if wheel.suffix != ".whl":
        raise LaneLaunchError("implementation artifact must be a wheel")
    if (
        len(expected_sha256) != 64
        or expected_sha256 != expected_sha256.lower()
        or any(character not in "0123456789abcdef" for character in expected_sha256)
    ):
        raise LaneLaunchError("implementation SHA-256 is invalid")
    if _sha256(wheel) != expected_sha256:
        raise LaneLaunchError("implementation wheel SHA-256 does not match")
    try:
        with zipfile.ZipFile(wheel) as archive:
            top_level: set[str] = set()
            for name in archive.namelist():
                member = PurePosixPath(name)
                if (
                    not member.parts
                    or member.is_absolute()
                    or any(part in {".", ".."} for part in member.parts)
                ):
                    raise LaneLaunchError(
                        "implementation wheel contains a noncanonical path"
                    )
                top_level.add(member.parts[0])
                if ".data/" in name:
                    raise LaneLaunchError(
                        "implementation wheel uses an unsupported install layout"
                    )
                if len(member.parts) == 1 and (
                    member.suffix == ".pth"
                    or member.name in {"sitecustomize.py", "usercustomize.py"}
                ):
                    raise LaneLaunchError(
                        "implementation wheel may not install Python path hooks"
                    )
            metadata_roots = {
                name
                for name in top_level
                if name.startswith("chemex-") and name.endswith(".dist-info")
            }
            if len(metadata_roots) != 1 or top_level != {"chemex", *metadata_roots}:
                raise LaneLaunchError(
                    "implementation wheel may contain only the ChemEx package and metadata"
                )
    except zipfile.BadZipFile as error:
        raise LaneLaunchError("implementation wheel is invalid") from error
    return wheel


def _verify_installed_implementation(wheel: Path) -> None:
    try:
        distribution = importlib.metadata.distribution("chemex")
        with zipfile.ZipFile(wheel) as archive:
            wheel_names = {
                name
                for name in archive.namelist()
                if not name.endswith("/") and not name.endswith(".dist-info/RECORD")
            }
            for name in wheel_names:
                if ".data/" in name:
                    raise LaneLaunchError(
                        "implementation wheel uses an unsupported install layout"
                    )
                installed = Path(str(distribution.locate_file(name))).resolve()
                if not installed.is_file() or installed.read_bytes() != archive.read(
                    name
                ):
                    raise LaneLaunchError(
                        "installed implementation does not match the verified wheel"
                    )
    except (
        importlib.metadata.PackageNotFoundError,
        OSError,
        zipfile.BadZipFile,
    ) as error:
        raise LaneLaunchError(
            "installed implementation cannot be verified against the wheel"
        ) from error

    package_root = Path(str(distribution.locate_file("chemex"))).resolve()
    expected_package = {name for name in wheel_names if name.startswith("chemex/")}
    try:
        installed_package = {
            path.relative_to(package_root.parent).as_posix()
            for path in package_root.rglob("*")
            if path.is_file()
            and path.suffix != ".pyc"
            and "__pycache__" not in path.parts
        }
    except OSError as error:
        raise LaneLaunchError(
            "installed implementation source is unavailable"
        ) from error
    if installed_package != expected_package:
        raise LaneLaunchError(
            "installed implementation file set does not match the wheel"
        )
    package_spec = importlib.util.find_spec("chemex")
    expected_origin = package_root / "__init__.py"
    if (
        package_spec is None
        or package_spec.origin is None
        or Path(package_spec.origin).resolve() != expected_origin
    ):
        raise LaneLaunchError(
            "importable implementation does not match the verified wheel"
        )


def _install_and_reexec(args: argparse.Namespace, original_args: list[str]) -> None:
    wheel = _verify_wheel(args.implementation_wheel, args.implementation_sha256)
    application_directory = args.application_directory.resolve()
    application_directory.mkdir(parents=True, exist_ok=True)
    prefix = Path(
        tempfile.mkdtemp(prefix="chemex-application-", dir=application_directory)
    ).resolve()
    environment = {
        name: value
        for name, value in os.environ.items()
        if not name.startswith("PYTHON")
    }
    try:
        venv.EnvBuilder(with_pip=True, system_site_packages=True).create(prefix)
        python = prefix / "bin/python"
        subprocess.run(  # noqa: S603 - the executable belongs to the new venv
            [
                str(python),
                "-I",
                "-m",
                "pip",
                "install",
                "--disable-pip-version-check",
                "--no-deps",
                "--no-index",
                str(wheel),
            ],
            check=True,
            capture_output=True,
            env=environment,
            text=True,
        )
        subprocess.run(  # noqa: S603 - the executable belongs to the new venv
            [str(python), "-I", "-m", "pip", "check"],
            check=True,
            capture_output=True,
            env=environment,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        shutil.rmtree(prefix, ignore_errors=True)
        raise LaneLaunchError(
            "exact implementation wheel installation failed"
        ) from error

    environment[_INSTALLED_STAGE_ENV] = "1"
    environment[_APPLICATION_PREFIX_ENV] = str(prefix)
    environment[_IMPLEMENTATION_DIGEST_ENV] = args.implementation_sha256
    try:
        os.execve(  # noqa: S606 - replacing this process is the authority boundary
            str(python),
            [str(python), "-I", str(Path(__file__).resolve()), *original_args],
            environment,
        )
    except OSError:
        shutil.rmtree(prefix, ignore_errors=True)
        raise


def _load_lane_attestor() -> ModuleType:
    import chemex

    module_name = "chemex.numerical_lanes"
    if module_name in sys.modules or hasattr(chemex, "numerical_lanes"):
        raise LaneLaunchError("application imported its own numerical-lane attestor")
    lane_root = Path(__file__).resolve().parent
    attestor_path = lane_root / "__init__.py"
    spec = importlib.util.spec_from_file_location(
        module_name,
        attestor_path,
        submodule_search_locations=[str(lane_root)],
    )
    if spec is None or spec.loader is None:
        raise LaneLaunchError("lane-owned attestor is unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    chemex.__dict__["numerical_lanes"] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        chemex.__dict__.pop("numerical_lanes", None)
        raise
    return module


def _validated_application_prefix(args: argparse.Namespace) -> Path:
    expected_prefix = Path(os.environ.get(_APPLICATION_PREFIX_ENV, "")).resolve()
    application_directory = args.application_directory.resolve()
    if (
        expected_prefix == Path(".").resolve()
        or expected_prefix.parent != application_directory
        or not expected_prefix.name.startswith("chemex-application-")
        or Path(sys.prefix).resolve() != expected_prefix
    ):
        raise LaneLaunchError("application environment does not match launcher state")
    return expected_prefix


def _run_qualification(args: argparse.Namespace) -> None:
    prefix = _validated_application_prefix(args)
    try:
        _verify_wheel(args.implementation_wheel, args.implementation_sha256)
        _verify_installed_implementation(args.implementation_wheel)
        if os.environ.get(_IMPLEMENTATION_DIGEST_ENV) != args.implementation_sha256:
            raise LaneLaunchError(
                "implementation authority does not match launcher state"
            )
        if not args.module.startswith(_QUALIFICATION_PREFIX):
            raise LaneLaunchError("launcher accepts only qualification modules")
        qualification_root = args.qualification_root.resolve()
        if not qualification_root.is_dir():
            raise LaneLaunchError("qualification checkout is unavailable")
        if (qualification_root / "chemex").exists() or (
            qualification_root / "chemex.py"
        ).exists():
            raise LaneLaunchError("qualification checkout may not shadow ChemEx")

        _load_lane_attestor()
        os.chdir(qualification_root)
        sys.path.append(str(qualification_root))
        try:
            qualification_spec = importlib.util.find_spec(args.module)
        except (ImportError, AttributeError, ValueError) as error:
            raise LaneLaunchError("qualification module is unavailable") from error
        if (
            qualification_spec is None
            or qualification_spec.origin is None
            or not Path(qualification_spec.origin)
            .resolve()
            .is_relative_to(qualification_root)
        ):
            raise LaneLaunchError("qualification module is outside the exact checkout")
        arguments = list(args.arguments)
        if arguments[:1] == ["--"]:
            arguments.pop(0)
        sys.argv = [args.module, *arguments]
        runpy.run_module(args.module, run_name="__main__", alter_sys=True)
    finally:
        shutil.rmtree(prefix, ignore_errors=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--implementation-wheel", type=Path, required=True)
    parser.add_argument("--implementation-sha256", required=True)
    parser.add_argument("--application-directory", type=Path, required=True)
    parser.add_argument("--qualification-root", type=Path, required=True)
    parser.add_argument("--module", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    return parser


def main(argv: list[str] | None = None) -> int:
    if not sys.flags.isolated:
        print(
            "numerical lane launcher: Python isolated mode (-I) is required",
            file=sys.stderr,
        )
        return 1
    original_args = list(sys.argv[1:] if argv is None else argv)
    args = _parser().parse_args(original_args)
    try:
        if os.environ.get(_INSTALLED_STAGE_ENV) == "1":
            _run_qualification(args)
        else:
            _install_and_reexec(args, original_args)
    except (LaneLaunchError, OSError, subprocess.CalledProcessError) as error:
        print(f"numerical lane launcher: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
