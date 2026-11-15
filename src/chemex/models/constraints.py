from __future__ import annotations

import math
from functools import lru_cache


def _validate_rates(*rates: float) -> None:
    for rate in rates:
        if not math.isfinite(rate) or rate < 0.0:
            msg = (
                f"Kinetic rate constants must be finite and non-negative, got {rate!r}"
            )
            raise ValueError(msg)


def _normalize_arborescence_weights(
    state_trees: tuple[tuple[str, tuple[tuple[float, ...], ...]], ...],
) -> dict[str, float] | None:
    weights = {
        name: math.fsum(math.prod(tree) for tree in trees)
        for name, trees in state_trees
    }
    total = math.fsum(weights.values())
    if total > 0.0 and math.isfinite(total):
        return {name: weight / total for name, weight in weights.items()}

    log_trees = {
        name: [
            math.fsum(math.log(edge) for edge in tree)
            for tree in trees
            if all(edge > 0.0 for edge in tree)
        ]
        for name, trees in state_trees
    }
    max_log = max(
        (log_weight for entries in log_trees.values() for log_weight in entries),
        default=None,
    )
    if max_log is None:
        return None

    scaled_weights = {
        name: math.fsum(math.exp(log_weight - max_log) for log_weight in entries)
        for name, entries in log_trees.items()
    }
    scaled_total = math.fsum(scaled_weights.values())
    return {name: weight / scaled_total for name, weight in scaled_weights.items()}


def _four_state_root_arborescences(
    kur: float,
    kvr: float,
    kwr: float,
    kuv: float,
    kvu: float,
    kuw: float,
    kwu: float,
    kvw: float,
    kwv: float,
) -> tuple[tuple[float, float, float], ...]:
    """Return the 16 directed spanning arborescences from (u, v, w) into root r."""
    return (
        (kur, kvr, kwr),
        (kur, kvr, kwu),
        (kur, kvr, kwv),
        (kur, kvu, kwr),
        (kur, kvw, kwr),
        (kuv, kvr, kwr),
        (kuw, kvr, kwr),
        (kur, kvu, kwu),
        (kur, kvu, kwv),
        (kur, kvw, kwu),
        (kuv, kvr, kwv),
        (kuv, kvr, kwu),
        (kuw, kvr, kwv),
        (kuw, kvw, kwr),
        (kuw, kvu, kwr),
        (kuv, kvw, kwr),
    )


@lru_cache(maxsize=100)
def pop_2st(kab: float = 0.0, kba: float = 0.0) -> dict[str, float]:
    _validate_rates(kab, kba)
    scale = max(kab, kba)
    if scale == 0.0:
        return {"pa": 1.0, "pb": 0.0}

    kab_scaled = kab / scale
    kba_scaled = kba / scale
    total = kab_scaled + kba_scaled

    return {"pa": kba_scaled / total, "pb": kab_scaled / total}


@lru_cache(maxsize=100)
def pop_3st(
    kab: float = 0.0,
    kba: float = 0.0,
    kac: float = 0.0,
    kca: float = 0.0,
    kbc: float = 0.0,
    kcb: float = 0.0,
) -> dict[str, float]:
    _validate_rates(kab, kba, kac, kca, kbc, kcb)
    scale = max(kab, kba, kac, kca, kbc, kcb)
    if scale == 0.0:
        return {"pa": 1.0, "pb": 0.0, "pc": 0.0}

    if kab == 0.0 and kba == 0.0 and kac == 0.0 and kca == 0.0:
        pb, pc = pop_2st(kbc, kcb).values()
        return {"pa": 0.0, "pb": pb, "pc": pc}

    if kab == 0.0 and kba == 0.0 and kbc == 0.0 and kcb == 0.0:
        pa, pc = pop_2st(kac, kca).values()
        return {"pa": pa, "pb": 0.0, "pc": pc}

    if kac == 0.0 and kca == 0.0 and kbc == 0.0 and kcb == 0.0:
        pa, pb = pop_2st(kab, kba).values()
        return {"pa": pa, "pb": pb, "pc": 0.0}

    # Directed spanning-tree weights give the stationary distribution without
    # mistaking a slow positive pathway for a structurally absent zero pathway.
    kab, kba, kac, kca, kbc, kcb = (
        rate / scale for rate in (kab, kba, kac, kca, kbc, kcb)
    )
    populations = _normalize_arborescence_weights(
        (
            ("pa", ((kba, kca), (kba, kcb), (kbc, kca))),
            ("pb", ((kab, kca), (kab, kcb), (kac, kcb))),
            ("pc", ((kac, kba), (kac, kbc), (kab, kbc))),
        )
    )
    if populations is None:
        msg = "Kinetic rates do not determine a unique stationary distribution"
        raise ValueError(msg)

    return populations


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
    _validate_rates(kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc)
    scale = max(kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc)
    if scale == 0.0:
        return {"pa": 1.0, "pb": 0.0, "pc": 0.0, "pd": 0.0}

    if (
        kab == 0.0
        and kba == 0.0
        and kac == 0.0
        and kca == 0.0
        and kad == 0.0
        and kda == 0.0
    ):
        pb, pc, pd = pop_3st(kbc, kcb, kbd, kdb, kcd, kdc).values()
        return {"pa": 0.0, "pb": pb, "pc": pc, "pd": pd}

    if (
        kab == 0.0
        and kba == 0.0
        and kbc == 0.0
        and kcb == 0.0
        and kbd == 0.0
        and kdb == 0.0
    ):
        pa, pc, pd = pop_3st(kac, kca, kad, kda, kcd, kdc).values()
        return {"pa": pa, "pb": 0.0, "pc": pc, "pd": pd}

    if (
        kac == 0.0
        and kca == 0.0
        and kbc == 0.0
        and kcb == 0.0
        and kcd == 0.0
        and kdc == 0.0
    ):
        pa, pb, pd = pop_3st(kab, kba, kad, kda, kbd, kdb).values()
        return {"pa": pa, "pb": pb, "pc": 0.0, "pd": pd}

    if (
        kad == 0.0
        and kda == 0.0
        and kbd == 0.0
        and kdb == 0.0
        and kcd == 0.0
        and kdc == 0.0
    ):
        pa, pb, pc = pop_3st(kab, kba, kac, kca, kbc, kcb).values()
        return {"pa": pa, "pb": pb, "pc": pc, "pd": 0.0}

    kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc = (
        rate / scale
        for rate in (kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc)
    )
    populations = _normalize_arborescence_weights(
        (
            (
                "pa",
                _four_state_root_arborescences(
                    kba, kca, kda, kbc, kcb, kbd, kdb, kcd, kdc
                ),
            ),
            (
                "pb",
                _four_state_root_arborescences(
                    kab, kcb, kdb, kac, kca, kad, kda, kcd, kdc
                ),
            ),
            (
                "pc",
                _four_state_root_arborescences(
                    kac, kbc, kdc, kab, kba, kad, kda, kbd, kdb
                ),
            ),
            (
                "pd",
                _four_state_root_arborescences(
                    kad, kbd, kcd, kab, kba, kac, kca, kbc, kcb
                ),
            ),
        )
    )
    if populations is None:
        msg = "Kinetic rates do not determine a unique stationary distribution"
        raise ValueError(msg)

    return populations
