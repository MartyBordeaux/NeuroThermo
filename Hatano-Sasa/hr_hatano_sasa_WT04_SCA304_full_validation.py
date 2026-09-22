#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Полная воспроизводимая проверка Хатано-Саса для зафиксированной пары
WT_04 <-> SCA3_04 в стохастической модели Hindmarsh-Rose.

Этот файл консолидирует расчеты, использованные для финальной проверки
медленного протокола (50 мс на каждый шаг пути):

1) детерминированная проверка реальной динамики x(t) на 31 точке пути;
2) стационарные стохастические выборки rho_ss(x,y,z|p) на 31 точке;
3) прямое оценивание соседних отношений плотностей rho_p/rho_{p+dp}
   методом k ближайших соседей (k=200, зафиксировано после независимой QC);
4) локальная hold-out проверка отношения плотностей;
5) сравнение скоростей протокола 10, 25 и 50 мс при одинаковом N=1000;
6) две независимые серии по N=5000 для 50-мс протокола;
7) объединенный результат N=10000, бутстрэп-ДИ, эффективный размер выборки
   и концентрация экспоненциальных весов;
8) таблицы сходимости и диагностические графики.

Научное ограничение
-------------------
Y является безразмерным функционалом Хатано-Саса. Код НЕ интерпретирует Y
как энергию, работу или теплоту в физических единицах.

Зафиксированные endpoints
-------------------------
Параметры b,r,s,kappa_I соответствуют строкам WT_04 и SCA3_04 из
cell_fit_summary.csv. Значения J являются зафиксированными входами q=0.75,
которые использовались в уже выполненной проверке:

WT_04:   J = 0.9508790471690817 pA/pF
SCA3_04: J = 5.2405950924786335 pA/pF

Ключевые численные настройки заморожены и не подбираются по результату:
- p: 31 положение от 0 до 1;
- dt stochastic = 0.025 ms;
- D = (0.0025, 0.01, 0.00025);
- burn-in = 2400 ms;
- stationary sample = 6000 ms, retain every 0.5 ms;
- 5 фиксированных seed для стационарных выборок;
- kNN ratio: 10000+10000 обучающих точек на соседнюю пару, k=200;
- медленный протокол: dwell = 50 ms;
- replicate A: seeds 9100000... / 9200000...;
- replicate B: seeds 9310000... / 9410000...;
- N=5000 на направление в каждой реплике.

Зависимости
-----------
python -m pip install numpy scipy pandas matplotlib numba

Пример запуска
--------------
python hr_hatano_sasa_WT04_SCA304_full_validation.py \
    --outdir hs_WT04_SCA304_full_validation \
    --cell-fit-summary cell_fit_summary.csv

Повторный запуск использует стационарный кэш. Чтобы пересоздать его:
python hr_hatano_sasa_WT04_SCA304_full_validation.py \
    --outdir hs_WT04_SCA304_full_validation --rebuild-cache

Для быстрой проверки только детерминированной динамики:
python hr_hatano_sasa_WT04_SCA304_full_validation.py \
    --outdir hs_WT04_SCA304_full_validation --preflight-only
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numba import njit
from scipy.spatial import cKDTree
from scipy.special import logsumexp


# =============================================================================
# 1. Замороженный бенчмарк
# =============================================================================

@dataclass(frozen=True)
class Endpoint:
    cell_id: str
    b: float
    r: float
    s: float
    kappa_I: float
    J: float


WT = Endpoint(
    cell_id="WT_04",
    b=3.07526787577599,
    r=0.03667512217506929,
    s=0.519708187871091,
    kappa_I=0.61900642338991,
    J=0.9508790471690817,
)

SCA3 = Endpoint(
    cell_id="SCA3_04",
    b=3.1841739616986,
    r=0.06852728277066979,
    s=1.64104809448976,
    kappa_I=0.22679725511273896,
    J=5.2405950924786335,
)

# HR fixed parameters
a = 1.0
c = 1.0
d = 5.0
xR = -1.6
X0 = np.array([-1.6, -10.0, 1.0], dtype=np.float64)

# Stochastic benchmark
DT_STOCH = 0.025  # ms
D_X, D_Y, D_Z = 0.0025, 0.01, 0.00025
BURN_MS = 2400.0
SAMPLE_MS = 6000.0
RETAIN_MS = 0.5
DENSITY_SEEDS = np.array(
    [20260803, 20261812, 20262821, 20263830, 20264839], dtype=np.int64
)
HOLDOUT_SEEDS = np.array(
    [30360803, 30361812, 30362821, 30363830, 30364839], dtype=np.int64
)

# Path benchmark
N_PATH = 31
P_GRID = np.linspace(0.0, 1.0, N_PATH)

# kNN density-ratio benchmark
KNN_K = 200
KNN_ALPHA = 1.0
PAIR_TRAIN_PER_CLASS = 10_000
PAIR_RNG_SEED = 8811
HOLDOUT_QC_PER_CLASS = 10_000

# Protocol benchmark
DWELL_COMPARE = (10.0, 25.0, 50.0)
N_MATCHED = 1000
DWELL_FINAL = 50.0
N_REPLICATE = 5000
CONVERGENCE_N = (100, 300, 1000, 3000, 5000)

# Trial seed blocks. Do not change if exact replication is desired.
SEEDS = {
    "A_forward": 9_100_000,
    "A_reverse": 9_200_000,
    "B_forward": 9_310_000,
    "B_reverse": 9_410_000,
}

# Deterministic preflight benchmark
DT_DET = 0.05  # ms
PRE_MS = 500.0
STIM_MS = 1000.0
SPIKE_THRESHOLD = 0.0

# Bootstrap benchmark
N_BOOT = 5000
BOOT_SEEDS = {
    "forward_A": 71_230_001,
    "forward_B": 71_230_002,
    "forward_pooled": 71_230_003,
    "reverse_A": 81_230_001,
    "reverse_B": 81_230_002,
    "reverse_pooled": 81_230_003,
}


# =============================================================================
# 2. HR dynamics
# =============================================================================

def hr_rhs(state: np.ndarray, b: float, r: float, s: float, kappa: float, J: float) -> np.ndarray:
    x, y, z = state
    return np.array(
        [
            y - a * x**3 + b * x**2 - z + kappa * J,
            c - d * x**2 - y,
            r * (s * (x - xR) - z),
        ],
        dtype=np.float64,
    )


def rk4_step(state: np.ndarray, dt: float, b: float, r: float, s: float, kappa: float, J: float) -> np.ndarray:
    k1 = hr_rhs(state, b, r, s, kappa, J)
    k2 = hr_rhs(state + 0.5 * dt * k1, b, r, s, kappa, J)
    k3 = hr_rhs(state + 0.5 * dt * k2, b, r, s, kappa, J)
    k4 = hr_rhs(state + dt * k3, b, r, s, kappa, J)
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def deterministic_preflight_row(p: float, pars: np.ndarray) -> Dict[str, float]:
    b, r, s, kappa, J = map(float, pars)
    state = X0.copy()

    # Published fitting convention: pre-simulation with no injected current.
    for _ in range(int(round(PRE_MS / DT_DET))):
        state = rk4_step(state, DT_DET, b, r, s, kappa, 0.0)

    n_steps = int(round(STIM_MS / DT_DET))
    xs = np.empty(n_steps + 1, dtype=np.float64)
    xs[0] = state[0]
    for i in range(n_steps):
        state = rk4_step(state, DT_DET, b, r, s, kappa, J)
        xs[i + 1] = state[0]

    # Upward x=0 crossings with linear interpolation.
    idx = np.flatnonzero((xs[:-1] < SPIKE_THRESHOLD) & (xs[1:] >= SPIKE_THRESHOLD))
    spike_t = []
    for ii in idx:
        x0v = xs[ii]
        x1v = xs[ii + 1]
        frac = (SPIKE_THRESHOLD - x0v) / (x1v - x0v)
        spike_t.append((ii + frac) * DT_DET)
    spike_t = np.asarray(spike_t, dtype=np.float64)

    n_spikes = int(spike_t.size)
    rate_hz = n_spikes / (STIM_MS / 1000.0)
    mean_isi = float(np.mean(np.diff(spike_t))) if n_spikes >= 2 else np.nan

    return {
        "p": p,
        "b": b,
        "r": r,
        "s": s,
        "kappa_I": kappa,
        "J": J,
        "kappaJ": kappa * J,
        "n_spikes": n_spikes,
        "rate_Hz": rate_hz,
        "mean_ISI_ms": mean_isi,
        "x_min": float(xs.min()),
        "x_max": float(xs.max()),
    }


@njit(cache=True)
def stationary_samples(
    b: float, r: float, s: float, kappa: float, J: float,
    dt: float, burn_steps: int, sample_steps: int, retain_every: int,
    seed: int,
) -> np.ndarray:
    np.random.seed(seed)
    x, y, z = -1.6, -10.0, 1.0

    sx = math.sqrt(2.0 * D_X * dt)
    sy = math.sqrt(2.0 * D_Y * dt)
    sz = math.sqrt(2.0 * D_Z * dt)

    for _ in range(burn_steps):
        fx = y - x*x*x + b*x*x - z + kappa*J
        fy = 1.0 - 5.0*x*x - y
        fz = r*(s*(x + 1.6) - z)
        x += fx*dt + sx*np.random.normal()
        y += fy*dt + sy*np.random.normal()
        z += fz*dt + sz*np.random.normal()

    n_keep = sample_steps // retain_every
    out = np.empty((n_keep, 3), dtype=np.float64)
    q = 0
    for i in range(sample_steps):
        fx = y - x*x*x + b*x*x - z + kappa*J
        fy = 1.0 - 5.0*x*x - y
        fz = r*(s*(x + 1.6) - z)
        x += fx*dt + sx*np.random.normal()
        y += fy*dt + sy*np.random.normal()
        z += fz*dt + sz*np.random.normal()
        if (i + 1) % retain_every == 0:
            out[q, 0] = x
            out[q, 1] = y
            out[q, 2] = z
            q += 1
    return out


@njit(cache=True)
def protocol_jump_states(
    param_path: np.ndarray,
    start_pool: np.ndarray,
    order: np.ndarray,
    dt: float,
    dwell_steps: int,
    trial_seeds: np.ndarray,
) -> np.ndarray:
    """State X_k is recorded immediately before each parameter jump k->k+1."""
    n_trials = len(trial_seeds)
    n_jumps = len(order) - 1
    states = np.empty((n_trials, n_jumps, 3), dtype=np.float64)

    sx = math.sqrt(2.0 * D_X * dt)
    sy = math.sqrt(2.0 * D_Y * dt)
    sz = math.sqrt(2.0 * D_Z * dt)

    for t in range(n_trials):
        np.random.seed(int(trial_seeds[t]))
        q = np.random.randint(0, start_pool.shape[0])
        x, y, z = start_pool[q, 0], start_pool[q, 1], start_pool[q, 2]

        for m in range(n_jumps):
            states[t, m, 0] = x
            states[t, m, 1] = y
            states[t, m, 2] = z

            k1 = int(order[m + 1])
            b = param_path[k1, 0]
            r = param_path[k1, 1]
            s = param_path[k1, 2]
            kappa = param_path[k1, 3]
            J = param_path[k1, 4]

            for _ in range(dwell_steps):
                fx = y - x*x*x + b*x*x - z + kappa*J
                fy = 1.0 - 5.0*x*x - y
                fz = r*(s*(x + 1.6) - z)
                x += fx*dt + sx*np.random.normal()
                y += fy*dt + sy*np.random.normal()
                z += fz*dt + sz*np.random.normal()

    return states


# =============================================================================
# 3. Path, validation and data loading
# =============================================================================

def build_path() -> np.ndarray:
    p = P_GRID
    out = np.empty((N_PATH, 5), dtype=np.float64)
    # frozen published convention: b,s,J linear; r,kappa logarithmic
    out[:, 0] = (1-p)*WT.b + p*SCA3.b
    out[:, 2] = (1-p)*WT.s + p*SCA3.s
    out[:, 4] = (1-p)*WT.J + p*SCA3.J
    out[:, 1] = np.exp((1-p)*np.log(WT.r) + p*np.log(SCA3.r))
    out[:, 3] = np.exp((1-p)*np.log(WT.kappa_I) + p*np.log(SCA3.kappa_I))
    return out


def verify_cell_fit_summary(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    needed = {"cell_id", "b", "r", "s", "kappa_I"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"В {csv_path} отсутствуют столбцы: {sorted(missing)}")

    rows = []
    for ep in (WT, SCA3):
        hit = df[df["cell_id"] == ep.cell_id]
        if len(hit) != 1:
            raise ValueError(f"Ожидалась ровно одна строка {ep.cell_id}, найдено {len(hit)}")
        row = hit.iloc[0]
        for key in ("b", "r", "s", "kappa_I"):
            observed = float(row[key])
            expected = float(getattr(ep, key))
            if not np.isclose(observed, expected, rtol=1e-10, atol=1e-12):
                raise ValueError(
                    f"Замороженный benchmark не совпадает: {ep.cell_id}.{key}: "
                    f"CSV={observed:.17g}, frozen={expected:.17g}"
                )
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# 4. Stationary cache and kNN ratio estimation
# =============================================================================

def build_or_load_stationary_cache(cache_dir: Path, path: np.ndarray, rebuild: bool):
    cache_dir.mkdir(parents=True, exist_ok=True)
    train_file = cache_dir / "stationary_train.npy"
    hold_file = cache_dir / "stationary_holdout.npy"
    meta_file = cache_dir / "stationary_meta.json"

    if (not rebuild) and train_file.exists() and hold_file.exists() and meta_file.exists():
        train = np.load(train_file, mmap_mode=None)
        hold = np.load(hold_file, mmap_mode=None)
        if train.shape[0] != N_PATH or hold.shape[0] != N_PATH:
            raise RuntimeError("Кэш имеет несовместимую форму; используйте --rebuild-cache")
        return train, hold

    burn_steps = int(round(BURN_MS / DT_STOCH))
    sample_steps = int(round(SAMPLE_MS / DT_STOCH))
    retain_every = int(round(RETAIN_MS / DT_STOCH))
    n_keep = sample_steps // retain_every
    n_per_p = len(DENSITY_SEEDS) * n_keep

    train = np.empty((N_PATH, n_per_p, 3), dtype=np.float64)
    hold = np.empty_like(train)

    print("[stationary] Генерация стационарных train/holdout выборок...", flush=True)
    for j in range(N_PATH):
        b, r, s, kappa, J = path[j]
        a_train = []
        a_hold = []
        for seed in DENSITY_SEEDS:
            a_train.append(
                stationary_samples(
                    b, r, s, kappa, J,
                    DT_STOCH, burn_steps, sample_steps, retain_every,
                    int(seed + 100003*j),
                )
            )
        for seed in HOLDOUT_SEEDS:
            a_hold.append(
                stationary_samples(
                    b, r, s, kappa, J,
                    DT_STOCH, burn_steps, sample_steps, retain_every,
                    int(seed + 100003*j),
                )
            )
        train[j] = np.vstack(a_train)
        hold[j] = np.vstack(a_hold)
        if j % 5 == 0 or j == N_PATH - 1:
            print(f"  p-index {j:02d}/{N_PATH-1}", flush=True)

    np.save(train_file, train)
    np.save(hold_file, hold)
    with open(meta_file, "w", encoding="utf-8") as f:
        json.dump(
            {
                "shape_train": list(train.shape),
                "shape_holdout": list(hold.shape),
                "dt_ms": DT_STOCH,
                "burn_ms": BURN_MS,
                "sample_ms": SAMPLE_MS,
                "retain_ms": RETAIN_MS,
                "density_seeds": DENSITY_SEEDS.tolist(),
                "holdout_seeds": HOLDOUT_SEEDS.tolist(),
            },
            f,
            indent=2,
        )
    return train, hold


def standardization_from_train(train: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # Exactly the decimation used in the final cached run.
    dec = np.vstack([train[j, ::10, :] for j in range(N_PATH)])
    mu = dec.mean(axis=0)
    sd = dec.std(axis=0)
    if np.any(sd <= 0):
        raise RuntimeError("Нулевая дисперсия при стандартизации")
    return mu, sd


def make_pair_training_sets(
    train: np.ndarray,
    mu: np.ndarray,
    sd: np.ndarray,
    cache_dir: Path,
    rebuild: bool,
) -> List[np.ndarray]:
    """Create the exact balanced 10000+10000 adjacent-pair samples."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    files = [cache_dir / f"pair_{j:02d}.npy" for j in range(N_PATH - 1)]
    if (not rebuild) and all(x.exists() for x in files):
        return [np.load(x) for x in files]

    rng = np.random.default_rng(PAIR_RNG_SEED)
    pairs = []
    for j in range(N_PATH - 1):
        t0 = (train[j] - mu) / sd
        t1 = (train[j+1] - mu) / sd
        i0 = rng.choice(len(t0), PAIR_TRAIN_PER_CLASS, replace=False)
        i1 = rng.choice(len(t1), PAIR_TRAIN_PER_CLASS, replace=False)
        X = np.vstack([t0[i0], t1[i1]])
        np.save(files[j], X)
        pairs.append(X)
    return pairs


def labels_for_pair() -> np.ndarray:
    return np.r_[
        np.zeros(PAIR_TRAIN_PER_CLASS, dtype=np.int8),
        np.ones(PAIR_TRAIN_PER_CLASS, dtype=np.int8),
    ]


def log_ratio_from_indices(labels: np.ndarray, idx: np.ndarray, k: int = KNN_K, alpha: float = KNN_ALPHA) -> np.ndarray:
    lab = labels[idx]
    n0 = np.sum(lab == 0, axis=1)
    n1 = k - n0
    return np.log(n0 + alpha) - np.log(n1 + alpha)


def knn_local_qc(
    pairs: List[np.ndarray],
    hold: np.ndarray,
    mu: np.ndarray,
    sd: np.ndarray,
    outdir: Path,
) -> pd.DataFrame:
    """Independent one-step identity check; does not tune k."""
    labels = labels_for_pair()
    rng = np.random.default_rng(18_811)  # QC-only; does not affect frozen pair construction.
    rows = []

    for j, Xtr in enumerate(pairs):
        tree = cKDTree(Xtr)
        h0 = (hold[j] - mu) / sd
        h1 = (hold[j+1] - mu) / sd
        q0 = rng.choice(len(h0), HOLDOUT_QC_PER_CLASS, replace=False)
        q1 = rng.choice(len(h1), HOLDOUT_QC_PER_CLASS, replace=False)
        idx0 = tree.query(h0[q0], k=KNN_K, workers=-1)[1]
        idx1 = tree.query(h1[q1], k=KNN_K, workers=-1)[1]
        lr0 = log_ratio_from_indices(labels, idx0)
        lr1 = log_ratio_from_indices(labels, idx1)

        # For exact ratios:
        # E_{rho_j}[rho_{j+1}/rho_j] = E exp(-log rho_j/rho_j1) = 1
        # E_{rho_j1}[rho_j/rho_j1] = E exp(+log rho_j/rho_j1) = 1
        Af = float(np.mean(np.exp(-lr0)))
        Ar = float(np.mean(np.exp(+lr1)))
        rows.append(
            {
                "pair": j,
                "p0": P_GRID[j],
                "p1": P_GRID[j+1],
                "A_forward": Af,
                "logA_forward": math.log(Af),
                "A_reverse": Ar,
                "logA_reverse": math.log(Ar),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "knn_step_identity.csv", index=False)
    return df


# =============================================================================
# 5. Protocol Y and statistics
# =============================================================================

def hs_stats(Y: np.ndarray) -> Dict[str, float]:
    Y = np.asarray(Y, dtype=np.float64)
    if np.any(~np.isfinite(Y)):
        Y = Y[np.isfinite(Y)]
    if len(Y) == 0:
        raise ValueError("Нет конечных значений Y")

    logw = -Y
    lse = logsumexp(logw)
    logmean = float(lse - math.log(len(Y)))
    ess = float(math.exp(2.0*lse - logsumexp(2.0*logw)))
    w = np.exp(logw - lse)
    ws = np.sort(w)[::-1]

    return {
        "N": int(len(Y)),
        "logmean": logmean,
        "meanexp": float(math.exp(logmean)),
        "meanY": float(np.mean(Y)),
        "medianY": float(np.median(Y)),
        "ess": ess,
        "ess_fraction": float(ess / len(Y)),
        "top1": float(ws[0]),
        "top10": float(ws[:min(10, len(ws))].sum()),
        "top100": float(ws[:min(100, len(ws))].sum()),
    }


def bootstrap_logmean(Y: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Y = np.asarray(Y, dtype=np.float64)
    n = len(Y)
    out = np.empty(n_boot, dtype=np.float64)
    # Looping avoids allocating a n_boot x N index matrix.
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yy = Y[idx]
        out[b] = logsumexp(-yy) - math.log(n)
    return out


def evaluate_protocol_states(
    state_sets: Dict[str, np.ndarray],
    pairs: List[np.ndarray],
    mu: np.ndarray,
    sd: np.ndarray,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Convert protocol states into forward/reverse Y for each named state set.

    state_sets keys are like 'd10_A', 'd25_A', 'd50_A', 'd50_B'.
    Each value is a tuple (XF, XR).
    """
    labels = labels_for_pair()
    out = {}
    work = {}
    for name, (XF, XR) in state_sets.items():
        work[name] = [np.zeros(XF.shape[0]), np.zeros(XR.shape[0])]

    for j, Xtr in enumerate(pairs):
        tree = cKDTree(Xtr)
        m = (N_PATH - 2) - j
        for name, (XF, XR) in state_sets.items():
            qf = (XF[:, j, :] - mu) / sd
            qr = (XR[:, m, :] - mu) / sd
            idxf = tree.query(qf, k=KNN_K, workers=-1)[1]
            idxr = tree.query(qr, k=KNN_K, workers=-1)[1]
            work[name][0] += log_ratio_from_indices(labels, idxf)
            work[name][1] += -log_ratio_from_indices(labels, idxr)
        if j % 5 == 0 or j == N_PATH - 2:
            print(f"  ratio pair {j:02d}/{N_PATH-2}", flush=True)

    for name, yy in work.items():
        out[name] = (yy[0], yy[1])
    return out


def simulate_named_protocol(
    path: np.ndarray,
    start0: np.ndarray,
    start1: np.ndarray,
    dwell_ms: float,
    n: int,
    forward_seed_base: int,
    reverse_seed_base: int,
) -> Tuple[np.ndarray, np.ndarray]:
    steps = int(round(dwell_ms / DT_STOCH))
    order_f = np.arange(N_PATH, dtype=np.int64)
    order_r = order_f[::-1].copy()
    sf = np.arange(forward_seed_base, forward_seed_base + n, dtype=np.int64)
    sr = np.arange(reverse_seed_base, reverse_seed_base + n, dtype=np.int64)

    XF = protocol_jump_states(path, start0, order_f, DT_STOCH, steps, sf)
    XR = protocol_jump_states(path, start1, order_r, DT_STOCH, steps, sr)
    return XF, XR


# =============================================================================
# 6. Output helpers
# =============================================================================

def save_convergence_rows(Y: np.ndarray, dwell: float, direction: str, replicate: str) -> List[dict]:
    rows = []
    for n in CONVERGENCE_N:
        if n <= len(Y):
            s = hs_stats(Y[:n])
            rows.append({"dwell_ms": dwell, "direction": direction, "replicate": replicate, **s})
    return rows


def make_plots(outdir: Path, conv: pd.DataFrame, validation: pd.DataFrame, Yfp: np.ndarray, Yrp: np.ndarray):
    # 1) convergence
    plt.figure(figsize=(7.5, 4.8))
    for (direction, rep), g in conv.groupby(["direction", "replicate"]):
        label = f"{direction}, {rep}"
        plt.plot(g["N"], g["logmean"], marker="o", label=label)
    plt.axhline(0.0, linestyle="--", linewidth=1.0)
    plt.xscale("log")
    plt.xlabel("Число траекторий")
    plt.ylabel("ln <exp(-Y)>; теория = 0")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "convergence_50ms_two_replicates.png", dpi=180)
    plt.close()

    # 2) pooled Y distributions
    plt.figure(figsize=(7.5, 4.8))
    plt.hist(Yfp, bins=60, density=True, alpha=0.55, label="WT->SCA3")
    plt.hist(Yrp, bins=60, density=True, alpha=0.55, label="SCA3->WT")
    plt.xlabel("Y")
    plt.ylabel("Плотность")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "Y_pooled_50ms_N10000.png", dpi=180)
    plt.close()

    # 3) central pooled result with CI as point plot
    pooled = validation[validation["replicate"] == "pooled_A+B"].copy()
    if len(pooled) == 2:
        x = np.arange(2)
        y = pooled["logmean"].to_numpy()
        lo = pooled["ci95_low"].to_numpy()
        hi = pooled["ci95_high"].to_numpy()
        err = np.vstack([y-lo, hi-y])
        plt.figure(figsize=(6.0, 4.5))
        plt.errorbar(x, y, yerr=err, fmt="o", capsize=5)
        plt.axhline(0.0, linestyle="--", linewidth=1.0)
        plt.xticks(x, pooled["direction"].tolist())
        plt.ylabel("ln <exp(-Y)>; теория = 0")
        plt.tight_layout()
        plt.savefig(outdir / "pooled_50ms_HS_with_CI.png", dpi=180)
        plt.close()


# =============================================================================
# 7. Main
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Полная замороженная проверка Хатано-Саса WT_04 <-> SCA3_04"
    )
    ap.add_argument("--outdir", type=Path, default=Path("hs_WT04_SCA304_full_validation"))
    ap.add_argument(
        "--cell-fit-summary",
        type=Path,
        default=None,
        help="Опционально: cell_fit_summary.csv для строгой проверки b,r,s,kappa_I",
    )
    ap.add_argument("--rebuild-cache", action="store_true", help="Пересоздать стационарные выборки и pair-кэш")
    ap.add_argument("--preflight-only", action="store_true", help="Выполнить только детерминированный preflight")
    args = ap.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = outdir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    if args.cell_fit_summary is not None:
        verified = verify_cell_fit_summary(args.cell_fit_summary)
        verified.to_csv(outdir / "verified_cell_fit_rows.csv", index=False)
        print("[QC] cell_fit_summary совпадает с замороженными b,r,s,kappa_I.", flush=True)

    path = build_path()
    pd.DataFrame(
        path, columns=["b", "r", "s", "kappa_I", "J"]
    ).assign(p=P_GRID).to_csv(outdir / "parameter_path.csv", index=False)

    # Save exact benchmark metadata.
    metadata = {
        "WT": asdict(WT),
        "SCA3": asdict(SCA3),
        "fixed": {"a": a, "c": c, "d": d, "xR": xR, "x0": X0.tolist()},
        "stochastic": {
            "dt_ms": DT_STOCH,
            "D": [D_X, D_Y, D_Z],
            "burn_ms": BURN_MS,
            "sample_ms": SAMPLE_MS,
            "retain_ms": RETAIN_MS,
            "density_seeds": DENSITY_SEEDS.tolist(),
            "holdout_seeds": HOLDOUT_SEEDS.tolist(),
        },
        "ratio": {
            "method": "balanced kNN class-count density ratio",
            "k": KNN_K,
            "alpha": KNN_ALPHA,
            "train_per_class_per_adjacent_pair": PAIR_TRAIN_PER_CLASS,
            "pair_rng_seed": PAIR_RNG_SEED,
        },
        "protocol": {
            "n_path": N_PATH,
            "dwell_compare_ms": list(DWELL_COMPARE),
            "dwell_final_ms": DWELL_FINAL,
            "N_matched": N_MATCHED,
            "N_replicate": N_REPLICATE,
            "trial_seed_blocks": SEEDS,
        },
        "bootstrap": {"n_boot": N_BOOT, "seeds": BOOT_SEEDS},
    }
    with open(outdir / "frozen_benchmark.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    # -------------------------------------------------------------------------
    # Deterministic preflight across the full path
    # -------------------------------------------------------------------------
    print("[preflight] Проверка x(t) на 31 положениях пути...", flush=True)
    preflight = pd.DataFrame(
        [deterministic_preflight_row(float(P_GRID[j]), path[j]) for j in range(N_PATH)]
    )
    preflight.to_csv(outdir / "preflight_path.csv", index=False)

    if preflight["n_spikes"].min() <= 0:
        raise RuntimeError("Preflight failed: на пути найдено состояние без спайков")

    print(
        f"  WT: {preflight.iloc[0].rate_Hz:.1f} Hz, ISI={preflight.iloc[0].mean_ISI_ms:.3f} ms; "
        f"SCA3: {preflight.iloc[-1].rate_Hz:.1f} Hz, ISI={preflight.iloc[-1].mean_ISI_ms:.3f} ms",
        flush=True,
    )

    if args.preflight_only:
        print(f"Готово: {outdir / 'preflight_path.csv'}")
        return

    # -------------------------------------------------------------------------
    # Stationary train/holdout distributions and fixed kNN ratio models
    # -------------------------------------------------------------------------
    train, hold = build_or_load_stationary_cache(cache_dir, path, args.rebuild_cache)
    mu, sd = standardization_from_train(train)
    np.save(cache_dir / "standardization_mu.npy", mu)
    np.save(cache_dir / "standardization_sd.npy", sd)

    pairs = make_pair_training_sets(train, mu, sd, cache_dir, args.rebuild_cache)
    qc = knn_local_qc(pairs, hold, mu, sd, outdir)
    qc_summary = {
        "sum_logA_forward": float(qc["logA_forward"].sum()),
        "sum_logA_reverse": float(qc["logA_reverse"].sum()),
        "max_abs_logA_forward": float(qc["logA_forward"].abs().max()),
        "max_abs_logA_reverse": float(qc["logA_reverse"].abs().max()),
    }
    with open(outdir / "knn_qc_summary.json", "w", encoding="utf-8") as f:
        json.dump(qc_summary, f, indent=2)

    # Start pools from the exact stationary endpoint ensembles.
    start0 = train[0]
    start1 = train[-1]

    # -------------------------------------------------------------------------
    # Protocol states: matched N=1000 for 10/25 ms, full A/B at 50 ms.
    # Replicate A 50ms contains the matched N=1000 as its prefix.
    # -------------------------------------------------------------------------
    state_sets: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    print("[protocol] dwell=10 ms, replicate A, N=1000", flush=True)
    state_sets["d10_A"] = simulate_named_protocol(
        path, start0, start1, 10.0, N_MATCHED, SEEDS["A_forward"], SEEDS["A_reverse"]
    )
    print("[protocol] dwell=25 ms, replicate A, N=1000", flush=True)
    state_sets["d25_A"] = simulate_named_protocol(
        path, start0, start1, 25.0, N_MATCHED, SEEDS["A_forward"], SEEDS["A_reverse"]
    )
    print("[protocol] dwell=50 ms, replicate A, N=5000", flush=True)
    state_sets["d50_A"] = simulate_named_protocol(
        path, start0, start1, 50.0, N_REPLICATE, SEEDS["A_forward"], SEEDS["A_reverse"]
    )
    print("[protocol] dwell=50 ms, replicate B, N=5000", flush=True)
    state_sets["d50_B"] = simulate_named_protocol(
        path, start0, start1, 50.0, N_REPLICATE, SEEDS["B_forward"], SEEDS["B_reverse"]
    )

    print("[ratio] Оценка Y прямым kNN-отношением плотностей...", flush=True)
    Y = evaluate_protocol_states(state_sets, pairs, mu, sd)

    # Save raw Y arrays.
    for name, (Yf, Yr) in Y.items():
        np.save(outdir / f"Y_forward_{name}.npy", Yf)
        np.save(outdir / f"Y_reverse_{name}.npy", Yr)

    # -------------------------------------------------------------------------
    # Matched-N slow-protocol comparison
    # -------------------------------------------------------------------------
    matched_rows = []
    for key, dwell in [("d10_A", 10.0), ("d25_A", 25.0), ("d50_A", 50.0)]:
        Yf, Yr = Y[key]
        for direction, yy in [("forward", Yf[:N_MATCHED]), ("reverse", Yr[:N_MATCHED])]:
            matched_rows.append({"dwell_ms": dwell, "direction": direction, **hs_stats(yy)})
    matched = pd.DataFrame(matched_rows)
    matched.to_csv(outdir / "slow_protocol_matchedN_comparison.csv", index=False)

    # -------------------------------------------------------------------------
    # 50-ms: two independent N=5000 replicates + pooled N=10000
    # -------------------------------------------------------------------------
    YfA, YrA = Y["d50_A"]
    YfB, YrB = Y["d50_B"]
    YfP = np.concatenate([YfA, YfB])
    YrP = np.concatenate([YrA, YrB])
    np.save(outdir / "Y_forward_pooled_N10000.npy", YfP)
    np.save(outdir / "Y_reverse_pooled_N10000.npy", YrP)

    validation_rows = []
    boot_output = {}
    configs = [
        ("forward", "A", YfA, BOOT_SEEDS["forward_A"]),
        ("forward", "B", YfB, BOOT_SEEDS["forward_B"]),
        ("forward", "pooled_A+B", YfP, BOOT_SEEDS["forward_pooled"]),
        ("reverse", "A", YrA, BOOT_SEEDS["reverse_A"]),
        ("reverse", "B", YrB, BOOT_SEEDS["reverse_B"]),
        ("reverse", "pooled_A+B", YrP, BOOT_SEEDS["reverse_pooled"]),
    ]
    for direction, rep, yy, bs in configs:
        st = hs_stats(yy)
        boot = bootstrap_logmean(yy, N_BOOT, bs)
        lo, hi = np.quantile(boot, [0.025, 0.975])
        validation_rows.append(
            {
                "direction": direction,
                "replicate": rep,
                **st,
                "ci95_low": float(lo),
                "ci95_high": float(hi),
            }
        )
        boot_output[(direction, rep)] = boot
        pd.DataFrame({"bootstrap_logmean": boot}).to_csv(
            outdir / f"bootstrap_{direction}_{rep.replace('+','plus')}.csv", index=False
        )

    validation = pd.DataFrame(validation_rows)
    # Useful column order.
    cols = [
        "direction", "replicate", "N", "logmean", "meanexp", "ci95_low", "ci95_high",
        "meanY", "medianY", "ess", "ess_fraction", "top1", "top10", "top100"
    ]
    validation = validation[cols]
    validation.to_csv(outdir / "validation_summary.csv", index=False)

    # Convergence for A and B.
    conv_rows = []
    for direction, rep, yy in [
        ("forward", "A", YfA), ("reverse", "A", YrA),
        ("forward", "B", YfB), ("reverse", "B", YrB),
    ]:
        conv_rows += save_convergence_rows(yy, 50.0, direction, rep)
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(outdir / "convergence_50ms_two_replicates.csv", index=False)

    # Weight concentration table for pooled results.
    wc = validation[validation["replicate"] == "pooled_A+B"][
        ["direction", "N", "ess", "ess_fraction", "top1", "top10", "top100"]
    ].copy()
    wc.to_csv(outdir / "weight_concentration_pooled.csv", index=False)

    make_plots(outdir, conv, validation, YfP, YrP)

    # Console summary.
    print("\n=== SLOW PROTOCOL, MATCHED N=1000 ===")
    print(matched[["dwell_ms", "direction", "logmean", "meanexp", "ess_fraction"]].to_string(index=False))
    print("\n=== FINAL 50 ms VALIDATION ===")
    print(validation.to_string(index=False))
    print("\n=== kNN LOCAL QC ===")
    print(json.dumps(qc_summary, indent=2))
    print(f"\nРезультаты сохранены в: {outdir.resolve()}")


if __name__ == "__main__":
    main()