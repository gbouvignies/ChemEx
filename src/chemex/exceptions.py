"""Semantic exceptions shared below ChemEx's terminal presentation boundary."""

from pathlib import Path


class ChemExError(Exception):
    """Known ChemEx failure safe for concise terminal presentation."""


class InputFileReadError(ChemExError, OSError):
    """A known user input could not be captured from a concrete path."""

    def __init__(self, path: Path, error: OSError) -> None:
        OSError.__init__(self, *error.args)
        self.path = path
        self.error = error


class ArtifactPublicationError(ChemExError, OSError):
    """A named ChemEx artifact could not be published at a concrete path."""

    def __init__(self, operation: str, path: Path, error: OSError) -> None:
        OSError.__init__(self, *error.args)
        self.operation = operation
        self.path = path
        self.error = error
