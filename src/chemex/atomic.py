"""Small filesystem primitives for atomic product publication."""

from __future__ import annotations

import ctypes
import errno
import os
import sys
import tempfile
from pathlib import Path

_AT_FDCWD = -100
_RENAME_NOREPLACE = 1
_RENAME_EXCL = 4


def write_text_atomic(destination: Path, content: str) -> None:
    """Replace one text file atomically from a same-directory temporary file."""
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}-",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            output.write(content)
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def publish_directory_noreplace(staging: Path, destination: Path) -> None:
    """Atomically rename a staged directory without replacing a destination."""
    if sys.platform == "win32":
        try:
            os.rename(staging, destination)
        except FileExistsError as error:
            raise FileExistsError(
                errno.EEXIST,
                "Publication destination exists",
                destination,
            ) from error
        return

    source_bytes = os.fsencode(staging)
    destination_bytes = os.fsencode(destination)
    if sys.platform == "linux":
        function_name = "renameat2"
        argument_types = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        arguments = (
            _AT_FDCWD,
            source_bytes,
            _AT_FDCWD,
            destination_bytes,
            _RENAME_NOREPLACE,
        )
    elif sys.platform == "darwin":
        function_name = "renamex_np"
        argument_types = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
        arguments = (source_bytes, destination_bytes, _RENAME_EXCL)
    else:
        raise OSError(
            errno.ENOTSUP,
            "Atomic no-replace directory rename is unavailable on this platform",
            destination,
        )

    libc = ctypes.CDLL(None, use_errno=True)
    try:
        rename_noreplace = getattr(libc, function_name)
    except AttributeError as error:
        raise OSError(
            errno.ENOSYS,
            f"{function_name} is unavailable for atomic publication",
            destination,
        ) from error
    rename_noreplace.argtypes = argument_types
    rename_noreplace.restype = ctypes.c_int
    result = rename_noreplace(*arguments)
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        raise FileExistsError(
            errno.EEXIST,
            "Publication destination exists",
            destination,
        )
    raise OSError(error_number, os.strerror(error_number), destination)
