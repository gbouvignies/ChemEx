from __future__ import annotations

import numpy as np

from chemex.nmr.is_liouvillian_engine import ISLiouvillianEngine
from chemex.nmr.liouvillian_views import reshape_single_liouvillian
from chemex.typing import Array

# A small value used for numerical stability
SMALL_VALUE = 1e-6


class SpectrometerAnalysis:
    """Small analysis surface built on top of a spectrometer engine."""

    def __init__(self, engine: ISLiouvillianEngine) -> None:
        self._engine = engine

    def calculate_shifts(self) -> Array:
        liouv = reshape_single_liouvillian(
            self._engine.l_free,
            self._engine.size,
            purpose="Shift eigenvalue calculation",
        )
        return np.linalg.eigvals(liouv).imag

    def calculate_r1rho(self) -> float:
        liouv = reshape_single_liouvillian(
            self._engine.l_free + self._engine.l_b1x_i,
            self._engine.size,
            purpose="R1rho eigenvalue calculation",
        )
        eigenvalues = np.linalg.eigvals(liouv)
        real_eigenvalues = eigenvalues[
            np.isclose(eigenvalues.imag, 0.0, atol=SMALL_VALUE)
        ].real
        if real_eigenvalues.size == 0:
            smallest_imag = float(np.min(np.abs(eigenvalues.imag)))
            msg = (
                "R1rho eigenvalue calculation did not find a nearly real eigenvalue "
                f"within atol={SMALL_VALUE}; smallest |imag| was {smallest_imag:.3e}."
            )
            raise ValueError(msg)
        return -float(np.max(real_eigenvalues))
