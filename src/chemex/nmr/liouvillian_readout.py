from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from chemex.nmr.detection import build_detection_vector
from chemex.nmr.magnetization import detect_signal
from chemex.typing import Array


class LiouvillianReadout:
    """Detection expression state and scalar readout for a Liouvillian."""

    def __init__(self, vectors: Mapping[str, Array]) -> None:
        self._vectors = vectors
        self._detection = ""
        self._detect_vector: Array = np.array([])

    @property
    def detection(self) -> str:
        return self._detection

    @detection.setter
    def detection(self, value: str) -> None:
        detect_vector = build_detection_vector(value, self._vectors).transpose()
        self._detection = value
        self._detect_vector = detect_vector

    def detect(self, magnetization: Array, weights: Array) -> float:
        return detect_signal(self._detect_vector, magnetization, weights)
