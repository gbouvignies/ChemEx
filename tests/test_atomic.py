from __future__ import annotations

from pathlib import Path

import pytest

from chemex.atomic import open_text_atomic, remove_paths_best_effort


def test_streaming_atomic_write_never_exposes_a_truncated_final_file(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "artifact.tsv"
    destination.write_text("previous complete artifact\n", encoding="utf-8")

    with pytest.raises(KeyboardInterrupt), open_text_atomic(destination) as output:
        output.write("header\n")
        output.write("first complete row\n")
        raise KeyboardInterrupt

    assert destination.read_text(encoding="utf-8") == "previous complete artifact\n"
    assert not tuple(tmp_path.glob(".artifact.tsv-*"))


def test_atomic_replace_failure_preserves_previous_complete_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "artifact.tsv"
    destination.write_text("previous complete artifact\n", encoding="utf-8")
    original_replace = Path.replace

    def fail_destination_replace(source: Path, target: Path) -> Path:
        if target == destination:
            raise OSError("replace refused")
        return original_replace(source, target)

    monkeypatch.setattr(Path, "replace", fail_destination_replace)
    with (
        pytest.raises(OSError, match="replace refused"),
        open_text_atomic(destination) as output,
    ):
        output.write("complete new artifact\n")

    assert destination.read_text(encoding="utf-8") == "previous complete artifact\n"
    assert not tuple(tmp_path.glob(".artifact.tsv-*"))


def test_atomic_cleanup_failure_is_noted_without_replacing_original_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "artifact.tsv"
    destination.write_text("previous complete artifact\n", encoding="utf-8")
    original_unlink = Path.unlink

    def refuse_temporary_cleanup(path: Path, *, missing_ok: bool = False) -> None:
        if path.name.startswith(".artifact.tsv-"):
            raise OSError("unlink refused")
        original_unlink(path, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", refuse_temporary_cleanup)
    with (
        pytest.raises(KeyboardInterrupt) as error_info,
        open_text_atomic(destination) as output,
    ):
        output.write("new partial artifact\n")
        raise KeyboardInterrupt

    assert destination.read_text(encoding="utf-8") == "previous complete artifact\n"
    assert tuple(tmp_path.glob(".artifact.tsv-*"))
    assert error_info.value.__notes__ == [
        "ChemEx could not remove an incomplete temporary artifact."
    ]


def test_best_effort_removal_attempts_every_path_and_preserves_original_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = tuple(tmp_path / f"artifact-{index}.tsv" for index in range(3))
    for path in paths:
        path.write_text("stale complete artifact\n", encoding="utf-8")
    original_unlink = Path.unlink
    attempted: list[Path] = []

    def fail_first_cleanup(path: Path, *, missing_ok: bool = False) -> None:
        attempted.append(path)
        if path == paths[0]:
            raise OSError("stale cleanup refused")
        original_unlink(path, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", fail_first_cleanup)
    interruption = KeyboardInterrupt()

    remove_paths_best_effort(paths, interruption)

    assert attempted == list(paths)
    assert paths[0].is_file()
    assert all(not path.exists() for path in paths[1:])
    assert interruption.__notes__ == ["ChemEx could not remove a stale artifact."]
