"""Model-owned physical feasibility declarations for relaxation operators."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True, order=True)
class RelaxationPsdBlock:
    """One symmetric relaxation-rate block that must remain PSD."""

    domain_id: str
    state: str
    diagonal_ids: tuple[str, ...]
    off_diagonal_ids: tuple[tuple[int, int, str], ...]

    def __post_init__(self) -> None:
        size = len(self.diagonal_ids)
        if size not in {2, 3}:
            raise ValueError("Relaxation PSD block has invalid diagonal scope")
        if any(
            not 0 <= row < column < size or not param_id
            for row, column, param_id in self.off_diagonal_ids
        ):
            raise ValueError("Relaxation PSD block has invalid off-diagonal scope")


def relaxation_blocks_identity(blocks: tuple[RelaxationPsdBlock, ...]) -> str:
    """Return a stable identity for canonical relaxation-domain declarations."""
    records = tuple(
        (
            block.domain_id,
            block.state,
            block.diagonal_ids,
            block.off_diagonal_ids,
        )
        for block in blocks
    )
    encoded = json.dumps(records, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class SealedRelaxationDomains:
    """Canonical immutable relaxation-domain declarations."""

    blocks: tuple[RelaxationPsdBlock, ...] = ()
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        blocks = tuple(sorted(set(self.blocks)))
        object.__setattr__(self, "blocks", blocks)
        object.__setattr__(self, "identity", relaxation_blocks_identity(blocks))
