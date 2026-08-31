from __future__ import annotations
import numpy as np
from numba import njit

@njit(cache=True)
def _vp_distance_numba(a: np.ndarray, b: np.ndarray, q_per_ms: float) -> float:
    n = a.size
    m = b.size
    prev = np.empty(m + 1, dtype=np.float64)
    curr = np.empty(m + 1, dtype=np.float64)
    for j in range(m + 1):
        prev[j] = float(j)
    for i in range(1, n + 1):
        curr[0] = float(i)
        ai = a[i - 1]
        for j in range(1, m + 1):
            delete = prev[j] + 1.0
            insert = curr[j - 1] + 1.0
            move = prev[j - 1] + q_per_ms * abs(ai - b[j - 1])
            best = delete
            if insert < best:
                best = insert
            if move < best:
                best = move
            curr[j] = best
        tmp = prev
        prev = curr
        curr = tmp
    return prev[m]


def victor_purpura(a, b, tau_ms: float = 10.0, normalize: bool = True) -> float:
    a = np.ascontiguousarray(np.sort(np.asarray(a, dtype=np.float64)))
    b = np.ascontiguousarray(np.sort(np.asarray(b, dtype=np.float64)))
    if tau_ms <= 0:
        raise ValueError("tau_ms must be > 0")
    q = 2.0 / float(tau_ms)
    d = float(_vp_distance_numba(a, b, q))
    if normalize:
        d /= max(a.size, b.size, 1)
    return d
