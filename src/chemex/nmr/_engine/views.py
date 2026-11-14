from __future__ import annotations

import numpy as np

from chemex.typing import Array


def reshape_single_liouvillian(
    liouvillian: Array,
    size: int,
    *,
    purpose: str,
) -> Array:
    expected_shape = (size, size)
    squeeze_axes = tuple(
        axis for axis, dim in enumerate(liouvillian.shape[:-2]) if dim == 1
    )
    squeezed = (
        np.squeeze(liouvillian, axis=squeeze_axes) if squeeze_axes else liouvillian
    )
    if squeezed.shape != expected_shape:
        msg = (
            f"{purpose} requires single-point B1 and Jeff distributions; "
            f"got Liouvillian shape {liouvillian.shape} "
            f"(squeezed to {squeezed.shape})."
        )
        raise ValueError(msg)
    return squeezed.reshape(expected_shape)
