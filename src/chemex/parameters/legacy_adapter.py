"""The single one-way adapter from legacy lmfit values to native central values."""

from __future__ import annotations

from collections.abc import Mapping

from lmfit import Parameters as LmfitParameters

from chemex.parameters.values import AnalysisValues


class LegacyValuesAdapter:
    """Mirror accepted legacy mutations without affecting legacy authority."""

    def __init__(self, analysis_values: AnalysisValues) -> None:
        self._analysis_values = analysis_values
        self._failure: Exception | None = None

    @property
    def failure(self) -> Exception | None:
        """Return the first failure that disabled native value mirroring."""
        return self._failure

    def try_commit(self, parameters: LmfitParameters) -> bool:
        """Best-effort commit of one complete lmfit parameter scope."""
        if self._failure is not None or not self._analysis_values.is_initialized:
            return False
        try:
            scope = tuple(parameters)
            candidate = {
                param_id: float(parameters[param_id].value) for param_id in scope
            }
        except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
            self._failure = error
            return False
        return self.try_commit_values(candidate)

    def try_commit_values(self, candidate: Mapping[str, float]) -> bool:
        """Best-effort commit of one native mapping from a legacy mutation."""
        if self._failure is not None or not self._analysis_values.is_initialized:
            return False
        try:
            expected = self._analysis_values.snapshot()
            self._analysis_values.commit(
                candidate,
                expected=expected,
                scope=tuple(candidate),
            )
        except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
            self._failure = error
            return False
        return True

    def disable(self, error: Exception) -> None:
        """Disable mirroring after an earlier native candidate failure."""
        if self._failure is None:
            self._failure = error

    def reset(self) -> None:
        """Clear failure state for a full session rebuild."""
        self._failure = None
