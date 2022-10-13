from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy import linalg


@lru_cache(maxsize=100)
def pop_2st(kab: float = 0.0, kba: float = 0.0) -> dict[str, float]:
    mat = np.array([[kab, -kba], [1.0, 1.0]])
    vec = np.array([0.0, 1.0])
    res = linalg.solve(mat, vec)
    return dict(zip(("pa", "pb"), res))


@lru_cache(maxsize=100)
def pop_3st(
    kab: float = 0.0,
    kba: float = 0.0,
    kac: float = 0.0,
    kca: float = 0.0,
    kbc: float = 0.0,
    kcb: float = 0.0,
) -> dict[str, float]:
    if all(np.isclose([kab + kba, kac + kca], [0.0, 0.0])):
        pb, pc = pop_2st(kbc, kcb).values()
        return {"pa": 0.0, "pb": pb, "pc": pc}

    if all(np.isclose([kab + kba, kbc + kcb], [0.0, 0.0])):
        pa, pc = pop_2st(kac, kca).values()
        return {"pa": pa, "pb": 0.0, "pc": pc}

    if all(np.isclose([kac + kca, kbc + kcb], [0.0, 0.0])):
        pa, pb = pop_2st(kab, kba).values()
        return {"pa": pa, "pb": pb, "pc": 0.0}

    mat = np.array(
        [
            [-kab - kac, kba, kca],
            [kab, -kba - kbc, kcb],
            [1.0, 1.0, 1.0],
        ]
    )
    vec = np.array([0.0, 0.0, 1.0])
    res = linalg.solve(mat, vec)

    return dict(zip(("pa", "pb", "pc"), res))


@lru_cache(maxsize=100)
def pop_4st(
    kab: float = 0.0,
    kba: float = 0.0,
    kac: float = 0.0,
    kca: float = 0.0,
    kad: float = 0.0,
    kda: float = 0.0,
    kbc: float = 0.0,
    kcb: float = 0.0,
    kbd: float = 0.0,
    kdb: float = 0.0,
    kcd: float = 0.0,
    kdc: float = 0.0,
) -> dict[str, float]:
    if all(np.isclose([kab + kba, kac + kca, kad + kda], [0.0, 0.0, 0.0])):
        pb, pc, pd = pop_3st(kbc, kcb, kbd, kdb, kcd, kdc).values()
        return {"pa": 0.0, "pb": pb, "pc": pc, "pd": pd}

    if all(np.isclose([kab + kba, kbc + kcb, kbd + kdb], [0.0, 0.0, 0.0])):
        pa, pc, pd = pop_3st(kac, kca, kad, kda, kcd, kdc).values()
        return {"pa": pa, "pb": 0.0, "pc": pc, "pd": pd}

    if all(np.isclose([kac + kca, kbc + kcb, kcd + kdc], [0.0, 0.0, 0.0])):
        pa, pb, pd = pop_3st(kab, kba, kad, kda, kbd, kdb).values()
        return {"pa": pa, "pb": pb, "pc": 0.0, "pd": pd}

    if all(np.isclose([kad + kda, kbd + kdb, kcd + kdc], [0.0, 0.0, 0.0])):
        pa, pb, pc = pop_3st(kab, kba, kac, kca, kbc, kcb).values()
        return {"pa": pa, "pb": pb, "pc": pc, "pd": 0.0}

    mat = np.array(
        [
            [-kab - kac - kad, kba, kca, kda],
            [kab, -kba - kbc - kbd, kcb, kdb],
            [kac, kbc, -kca - kcb - kcd, kdc],
            [1.0, 1.0, 1.0, 1.0],
        ]
    )
    vec = np.array([0.0, 0.0, 0.0, 1.0])
    res = linalg.solve(mat, vec)

    return dict(zip(("pa", "pb", "pc", "pd"), res))
