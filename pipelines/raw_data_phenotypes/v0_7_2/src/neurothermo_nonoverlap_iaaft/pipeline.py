from __future__ import annotations

import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import __version__
from .config import REQUIRED_UPSTREAM
from .metrics import (
    iaaft_surrogate,
    ordinal_pi_lags,
    ordinal_predictive_information,
    shuffled_surrogate,
    surrogate_fidelity,
)
from .raw import resolve_source_path, stationary_trace
from .statistics import (
    cell_auc,
    exact_cell_difference,
    exact_lag_family_tests,
    label_masks,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_dump(path: Path, payload: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)


def _read_domain(config: Mapping[str, Any], upstream: Path) -> pd.DataFrame:
    sweeps = pd.read_csv(upstream / "sweep_features.csv")
    if "animal_id" in sweeps.columns:
        sweeps = sweeps.drop(columns=["animal_id"])
    currents = [int(x) for x in config["analysis"]["currents_pA"]]
    groups = [str(x) for x in config["analysis"]["groups"]]
    domain = sweeps[sweeps.current_pA.isin(currents) & sweeps.group.isin(groups)].copy()
    domain["cell_id"] = domain.cell_id.astype(str)
    domain["current_pA"] = domain.current_pA.astype(int)
    return domain.sort_values(["group", "cell_id", "current_pA"]).reset_index(drop=True)


def _lag_samples(config: Mapping[str, Any]) -> Dict[int, int]:
    dt = float(config["surrogates"]["resample_dt_ms"])
    output = {}
    for lag_ms in config["surrogates"]["code_pair_lags_ms"]:
        lag_ms = int(lag_ms)
        samples = int(round(float(lag_ms) / dt))
        if not np.isclose(samples * dt, float(lag_ms), rtol=0.0, atol=1e-9):
            raise ValueError("Every code-pair lag must be an integer multiple of resample_dt_ms.")
        output[lag_ms] = samples
    return output


def validate_inputs(config: Mapping[str, Any], upstream: Path, raw_root: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    missing = [str(upstream / x) for x in REQUIRED_UPSTREAM if not (upstream / x).is_file()]
    rows.append({"check": "required_v031_files", "passed": not missing,
                 "detail": "none" if not missing else "; ".join(missing)})
    if missing:
        return pd.DataFrame(rows)
    domain = _read_domain(config, upstream)
    required = {
        "group", "cell_id", "source_path", "sweep_index", "current_pA",
        "stim_start_s", "stim_end_s", "stationary_samples", "predictive_information_nats",
    }
    absent = sorted(required - set(domain.columns))
    rows.append({"check": "required_columns", "passed": not absent,
                 "detail": "none" if not absent else ", ".join(absent)})
    if absent:
        return pd.DataFrame(rows)

    groups = [str(x) for x in config["analysis"]["groups"]]
    currents = [int(x) for x in config["analysis"]["currents_pA"]]
    expected = config["analysis"]["expected_cells"]
    counts = domain[["group", "cell_id"]].drop_duplicates().groupby("group").size().to_dict()
    rows.append({
        "check": "frozen_cell_counts",
        "passed": all(int(counts.get(group, 0)) == int(expected[group]) for group in groups),
        "detail": json.dumps(counts, sort_keys=True),
    })
    duplicates = int(domain.duplicated(["cell_id", "current_pA"]).sum())
    rows.append({"check": "one_row_per_cell_current", "passed": duplicates == 0,
                 "detail": "duplicates={}".format(duplicates)})
    per_cell = domain.groupby("cell_id").current_pA.nunique()
    incomplete = per_cell[per_cell != len(currents)].to_dict()
    rows.append({"check": "complete_current_grid", "passed": not incomplete,
                 "detail": "none" if not incomplete else json.dumps(incomplete, sort_keys=True)})
    rows.append({
        "check": "expected_total_rows",
        "passed": len(domain) == sum(int(expected[x]) for x in groups) * len(currents),
        "detail": "rows={}".format(len(domain)),
    })
    missing_raw = []
    resolved = set()
    for _, row in domain[["group", "source_path"]].drop_duplicates().iterrows():
        try:
            resolved.add(str(resolve_source_path(str(row.source_path), str(row.group), raw_root)))
        except FileNotFoundError as exc:
            missing_raw.append(str(exc))
    rows.append({"check": "raw_trace_files", "passed": not missing_raw,
                 "detail": "resolved_unique_files={}".format(len(resolved)) if not missing_raw else " | ".join(missing_raw)})

    lag_map = _lag_samples(config)
    order = int(config["surrogates"]["permutation_order"])
    delay_samples = int(round(float(config["surrogates"]["permutation_delay_ms"]) /
                              float(config["surrogates"]["resample_dt_ms"])))
    ordinal_width = (order - 1) * delay_samples + 1
    primary_lag = int(config["surrogates"]["primary_code_pair_lag_ms"])
    rows.append({
        "check": "primary_ordinal_windows_do_not_overlap",
        "passed": primary_lag in lag_map and lag_map.get(primary_lag, 0) >= ordinal_width,
        "detail": "ordinal_width_samples={}; primary_code_lag_samples={}".format(
            ordinal_width, lag_map.get(primary_lag)
        ),
    })
    rows.append({
        "check": "all_lags_unique_and_ordered",
        "passed": list(lag_map) == sorted(set(lag_map)),
        "detail": json.dumps(lag_map, sort_keys=True),
    })
    rows.append({"check": "cell_level_only", "passed": "animal_id" not in domain.columns,
                 "detail": "animal identifiers removed; no animal-level inference"})
    return pd.DataFrame(rows)


def _seed(master: int, cell_id: str, current_pA: int, kind: str) -> int:
    key = "{}|{}|{}|{}".format(master, cell_id, current_pA, kind).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "little")


def _surrogate_worker(payload):
    (
        cell_id, current_pA, values, order, delay, code_lags, n_shuffle, n_iaaft,
        master_seed, max_iterations, improvement_tolerance, patience, maximum_acf_lag,
    ) = payload
    shuffle_rng = np.random.default_rng(_seed(master_seed, cell_id, current_pA, "shuffle_v072"))
    iaaft_rng = np.random.default_rng(_seed(master_seed, cell_id, current_pA, "iaaft_v072"))
    shuffle_pi = np.empty((n_shuffle, len(code_lags)), dtype=float)
    iaaft_pi = np.empty((n_iaaft, len(code_lags)), dtype=float)
    diagnostics = []
    for index in range(n_shuffle):
        surrogate = shuffled_surrogate(values, shuffle_rng)
        shuffle_pi[index] = ordinal_pi_lags(surrogate, order, delay, code_lags)
    for index in range(n_iaaft):
        surrogate, convergence = iaaft_surrogate(
            values, iaaft_rng, max_iterations=max_iterations,
            improvement_tolerance=improvement_tolerance, patience=patience,
        )
        iaaft_pi[index] = ordinal_pi_lags(surrogate, order, delay, code_lags)
        fidelity = surrogate_fidelity(values, surrogate, maximum_acf_lag)
        fidelity.update(convergence)
        diagnostics.append(fidelity)
    return shuffle_pi, iaaft_pi, diagnostics


def _surrogate_statistics(observed: float, values: np.ndarray, prefix: str) -> Dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    suffixes = (
        "_median_nats", "_mean_nats", "_sd_nats", "_q025_nats", "_q25_nats",
        "_q75_nats", "_q975_nats", "_centered_PI_nats", "_z",
        "_p_lower", "_p_upper", "_p_two_sided",
    )
    if len(finite) == 0:
        return {prefix + suffix: np.nan for suffix in suffixes}
    median = float(np.median(finite))
    mean = float(np.mean(finite))
    sd = float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan
    quantiles = np.quantile(finite, [0.025, 0.25, 0.75, 0.975])
    centered_distance = abs(float(observed) - median)
    return {
        prefix + "_median_nats": median,
        prefix + "_mean_nats": mean,
        prefix + "_sd_nats": sd,
        prefix + "_q025_nats": float(quantiles[0]),
        prefix + "_q25_nats": float(quantiles[1]),
        prefix + "_q75_nats": float(quantiles[2]),
        prefix + "_q975_nats": float(quantiles[3]),
        prefix + "_centered_PI_nats": float(observed - median),
        prefix + "_z": float((observed - mean) / sd) if np.isfinite(sd) and sd > 1e-15 else np.nan,
        prefix + "_p_lower": float((1 + np.sum(finite <= observed)) / (1 + len(finite))),
        prefix + "_p_upper": float((1 + np.sum(finite >= observed)) / (1 + len(finite))),
        prefix + "_p_two_sided": float((1 + np.sum(np.abs(finite - median) >= centered_distance)) / (1 + len(finite))),
    }


def _extract_and_compute(config, domain: pd.DataFrame, raw_root: Path):
    trace_rows = []
    traces = []
    order = int(config["surrogates"]["permutation_order"])
    dt_ms = float(config["surrogates"]["resample_dt_ms"])
    delay = max(1, int(round(float(config["surrogates"]["permutation_delay_ms"]) / dt_ms)))
    lag_map = _lag_samples(config)
    lag_ms_values = list(lag_map)
    code_lags = [lag_map[x] for x in lag_ms_values]
    for _, row in domain.iterrows():
        values, metadata = stationary_trace(row, raw_root, config)
        adjacent = ordinal_predictive_information(values, order, delay, code_lag=1)
        observed_lags = ordinal_pi_lags(values, order, delay, code_lags)
        base = {
            "group": str(row.group), "cell_id": str(row.cell_id),
            "current_pA": int(row.current_pA), "sweep_index": int(row.sweep_index),
            "source_path_upstream": str(row.source_path),
            "predictive_information_upstream_nats": float(row.predictive_information_nats),
            "adjacent_PI_recomputed_nats": float(adjacent),
            "absolute_reproduction_error_nats": abs(adjacent - float(row.predictive_information_nats)),
            "upstream_stationary_samples": int(row.stationary_samples),
        }
        for lag_ms, value in zip(lag_ms_values, observed_lags):
            base["observed_PI_lag_{}ms_nats".format(lag_ms)] = float(value)
        base.update(metadata)
        trace_rows.append(base)
        traces.append(values)

    surrogate_config = config["surrogates"]
    n_shuffle = int(surrogate_config["n_shuffle"])
    n_iaaft = int(surrogate_config["n_iaaft"])
    maximum_acf_lag = max(1, int(round(float(surrogate_config["fidelity_acf_max_lag_ms"]) / dt_ms)))
    payloads = [
        (
            str(row.cell_id), int(row.current_pA), values, order, delay, code_lags,
            n_shuffle, n_iaaft, int(config["analysis"]["random_seed"]),
            int(surrogate_config["iaaft_max_iterations"]),
            float(surrogate_config["iaaft_improvement_tolerance"]),
            int(surrogate_config["iaaft_patience"]), maximum_acf_lag,
        )
        for (_, row), values in zip(domain.iterrows(), traces)
    ]
    workers = int(config.get("runtime", {}).get("workers", 2))
    if workers <= 1:
        results = [_surrogate_worker(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_surrogate_worker, payloads, chunksize=1))

    summary_rows = []
    long_rows = []
    arrays = {}
    fidelity_rows = []
    for base, (shuffle_pi, iaaft_pi, diagnostics) in zip(trace_rows, results):
        summary = dict(base)
        key = (str(base["cell_id"]), int(base["current_pA"]))
        arrays[key] = {"shuffle": shuffle_pi, "iaaft": iaaft_pi}
        for lag_index, lag_ms in enumerate(lag_ms_values):
            observed = float(base["observed_PI_lag_{}ms_nats".format(lag_ms)])
            summary.update(_surrogate_statistics(
                observed, shuffle_pi[:, lag_index], "shuffle_lag_{}ms".format(lag_ms)
            ))
            summary.update(_surrogate_statistics(
                observed, iaaft_pi[:, lag_index], "iaaft_lag_{}ms".format(lag_ms)
            ))
        summary["n_shuffle"] = int(n_shuffle)
        summary["n_iaaft"] = int(n_iaaft)
        sweep_diag = pd.DataFrame(diagnostics)
        for column in ["spectral_amplitude_nrmse", "log_psd_rmse", "acf_rmse", "iterations"]:
            values = sweep_diag[column].to_numpy(float)
            summary["iaaft_{}_median".format(column)] = float(np.median(values))
            summary["iaaft_{}_q95".format(column)] = float(np.quantile(values, .95))
            summary["iaaft_{}_max".format(column)] = float(np.max(values))
        summary["iaaft_converged_fraction"] = float(sweep_diag.converged.astype(bool).mean())
        summary["iaaft_amplitude_exact_fraction"] = float(sweep_diag.amplitude_distribution_exact.astype(bool).mean())
        summary_rows.append(summary)

        for kind, matrix in (("shuffle", shuffle_pi), ("iaaft", iaaft_pi)):
            for index in range(len(matrix)):
                record = {
                    "group": base["group"], "cell_id": base["cell_id"],
                    "current_pA": int(base["current_pA"]), "surrogate_type": kind,
                    "surrogate_index": int(index),
                }
                for lag_index, lag_ms in enumerate(lag_ms_values):
                    record["PI_lag_{}ms_nats".format(lag_ms)] = float(matrix[index, lag_index])
                if kind == "iaaft":
                    record.update(diagnostics[index])
                    fidelity = dict(record)
                    for lag_ms in lag_ms_values:
                        fidelity.pop("PI_lag_{}ms_nats".format(lag_ms), None)
                    fidelity_rows.append(fidelity)
                long_rows.append(record)
    return (
        pd.DataFrame(summary_rows), pd.DataFrame(long_rows), arrays,
        pd.DataFrame(fidelity_rows), lag_ms_values,
    )


def _ordered_cells(domain: pd.DataFrame, groups: Sequence[str]):
    order = domain[["group", "cell_id"]].drop_duplicates().sort_values(
        ["group", "cell_id"],
        key=lambda column: column.map({groups[0]: 0, groups[1]: 1}) if column.name == "group" else column,
    )
    cell_ids = order.cell_id.astype(str).tolist()
    by_cell = order.set_index("cell_id").group.to_dict()
    labels = np.asarray([by_cell[cell] == groups[1] for cell in cell_ids], dtype=bool)
    return cell_ids, labels


def _fidelity_summary(config, fidelity: pd.DataFrame):
    by_sweep_rows = []
    for keys, sub in fidelity.groupby(["group", "cell_id", "current_pA"], sort=True):
        row = {"group": keys[0], "cell_id": keys[1], "current_pA": int(keys[2]), "n_iaaft": int(len(sub))}
        for column in ["spectral_amplitude_nrmse", "log_psd_rmse", "acf_rmse", "iterations"]:
            values = sub[column].to_numpy(float)
            row[column + "_median"] = float(np.median(values))
            row[column + "_q95"] = float(np.quantile(values, .95))
            row[column + "_max"] = float(np.max(values))
        row["converged_fraction"] = float(sub.converged.astype(bool).mean())
        row["amplitude_exact_fraction"] = float(sub.amplitude_distribution_exact.astype(bool).mean())
        by_sweep_rows.append(row)
    by_sweep = pd.DataFrame(by_sweep_rows)

    gates = config["surrogates"]["fidelity_gates"]
    spectral_q95 = float(fidelity.spectral_amplitude_nrmse.quantile(.95))
    acf_q95 = float(fidelity.acf_rmse.quantile(.95))
    converged_fraction = float(fidelity.converged.astype(bool).mean())
    amplitude_exact_fraction = float(fidelity.amplitude_distribution_exact.astype(bool).mean())
    maximum_amplitude_error = float(fidelity.amplitude_sorted_max_abs_error.max())
    overall = pd.DataFrame([
        {"check": "IAAFT_spectral_amplitude_NRMSE_q95", "value": spectral_q95,
         "threshold": float(gates["maximum_spectral_amplitude_nrmse_q95"]),
         "comparison": "<=", "passed": spectral_q95 <= float(gates["maximum_spectral_amplitude_nrmse_q95"])},
        {"check": "IAAFT_ACF_RMSE_q95", "value": acf_q95,
         "threshold": float(gates["maximum_acf_rmse_q95"]),
         "comparison": "<=", "passed": acf_q95 <= float(gates["maximum_acf_rmse_q95"])},
        {"check": "IAAFT_converged_fraction", "value": converged_fraction,
         "threshold": float(gates["minimum_converged_fraction"]),
         "comparison": ">=", "passed": converged_fraction >= float(gates["minimum_converged_fraction"])},
        {"check": "IAAFT_amplitude_exact_fraction", "value": amplitude_exact_fraction,
         "threshold": 1.0, "comparison": "==", "passed": amplitude_exact_fraction == 1.0},
        {"check": "IAAFT_maximum_sorted_amplitude_error", "value": maximum_amplitude_error,
         "threshold": 0.0, "comparison": "==", "passed": maximum_amplitude_error == 0.0},
    ])
    return by_sweep, overall


def _families(lag_ms_values: Sequence[int]):
    return {
        "shuffle": {int(lag): "shuffle_lag_{}ms_centered_PI_nats".format(lag) for lag in lag_ms_values},
        "iaaft": {int(lag): "iaaft_lag_{}ms_centered_PI_nats".format(lag) for lag in lag_ms_values},
    }


def _count_stability(config, summary, arrays, cell_ids, labels, currents, masks, lag_ms_values):
    rows = []
    requested = sorted(set(int(x) for x in config["surrogates"]["count_sensitivity"]))
    observed_columns = ["observed_PI_lag_{}ms_nats".format(lag) for lag in lag_ms_values]
    for kind, maximum in (("shuffle", int(config["surrogates"]["n_shuffle"])),
                          ("iaaft", int(config["surrogates"]["n_iaaft"]))):
        for count in [x for x in requested if x <= maximum]:
            for lag_index, (lag_ms, observed_column) in enumerate(zip(lag_ms_values, observed_columns)):
                work = summary[["group", "cell_id", "current_pA", observed_column]].copy()
                medians = []
                for _, row in work.iterrows():
                    values = arrays[(str(row.cell_id), int(row.current_pA))][kind][:count, lag_index]
                    medians.append(float(np.median(values)))
                metric = "centered"
                work[metric] = work[observed_column].to_numpy(float) - np.asarray(medians)
                matrix = work.pivot(index="cell_id", columns="current_pA", values=metric).reindex(
                    index=cell_ids, columns=currents
                ).to_numpy(float)
                auc = cell_auc(matrix, currents)
                observed, p_expected, p_two, _ = exact_cell_difference(auc, masks, labels, -1.0)
                rows.append({
                    "surrogate_type": kind, "lag_ms": int(lag_ms), "n_surrogates": count,
                    "observed_raw_AUC_difference_SCA3_minus_WT": observed,
                    "observed_expected_direction_AUC_difference_WT_minus_SCA3": -observed,
                    "exact_p_AUC_expected_direction": p_expected,
                    "exact_p_AUC_two_sided": p_two,
                })
    return pd.DataFrame(rows)


def _group_curves(summary, lag_ms_values):
    rows = []
    metrics = []
    for lag in lag_ms_values:
        metrics.extend([
            "observed_PI_lag_{}ms_nats".format(lag),
            "shuffle_lag_{}ms_centered_PI_nats".format(lag),
            "iaaft_lag_{}ms_centered_PI_nats".format(lag),
        ])
    for (group, current), sub in summary.groupby(["group", "current_pA"], sort=True):
        for metric in metrics:
            values = sub[metric].dropna()
            rows.append({
                "group": group, "current_pA": int(current), "metric": metric,
                "n": int(len(values)), "median": float(values.median()),
                "q25": float(values.quantile(.25)), "q75": float(values.quantile(.75)),
            })
    return pd.DataFrame(rows)


def _plot_band(ax, frame, metric, title, groups):
    colors = {groups[0]: "#2878b5", groups[1]: "#be2334"}
    for group in groups:
        sub = frame[(frame.group == group) & (frame.metric == metric)].sort_values("current_pA")
        x = sub.current_pA.to_numpy(float)
        median = sub["median"].to_numpy(float)
        q25 = sub.q25.to_numpy(float)
        q75 = sub.q75.to_numpy(float)
        ax.plot(x, median, marker="o", color=colors[group], label=group)
        ax.fill_between(x, q25, q75, color=colors[group], alpha=.18)
    ax.axhline(0.0, color="black", lw=1, ls="--")
    ax.set(title=title, xlabel="Injected current (pA)", ylabel="Surrogate-centered PI (nats)")
    ax.legend(frameon=False)


def _plots(output, config, curves, auc_tests, current_tests, stability, fidelity):
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    groups = [str(x) for x in config["analysis"]["groups"]]
    primary_lag = int(config["surrogates"]["primary_code_pair_lag_ms"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _plot_band(axes[0], curves, "shuffle_lag_{}ms_centered_PI_nats".format(primary_lag),
               "Primary: non-overlap PI centered on shuffle null", groups)
    _plot_band(axes[1], curves, "iaaft_lag_{}ms_centered_PI_nats".format(primary_lag),
               "Sensitivity: non-overlap PI centered on IAAFT null", groups)
    fig.tight_layout()
    fig.savefig(figures / "primary_nonoverlap_PI_curves.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, family in zip(axes, ["shuffle", "iaaft"]):
        sub = auc_tests[auc_tests.surrogate_family == family].sort_values("lag_ms")
        ax.plot(sub.lag_ms, sub.observed_expected_direction_AUC_difference_WT_minus_SCA3,
                marker="o", color="#4c78a8")
        ax.axhline(0.0, color="black", lw=1)
        ax.set(title="{}-centered".format(family.upper()), xlabel="Non-overlap code lag (ms)")
    axes[0].set_ylabel("WT minus SCA3 cell-AUC difference (nats)")
    fig.tight_layout()
    fig.savefig(figures / "lag_sensitivity_cell_AUC.png", dpi=180)
    plt.close(fig)

    primary_current = current_tests[
        (current_tests.surrogate_family == "shuffle") & (current_tests.lag_ms == primary_lag)
    ].sort_values("current_pA")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(primary_current.current_pA, primary_current.observed_expected_direction_difference_WT_minus_SCA3,
            marker="o", color="#4c78a8")
    significant = primary_current.exact_p_expected_direction_maxT_across_lags_and_currents <= .05
    ax.scatter(primary_current.loc[significant, "current_pA"],
               primary_current.loc[significant, "observed_expected_direction_difference_WT_minus_SCA3"],
               color="#be2334", s=65, zorder=3)
    ax.axhline(0.0, color="black", lw=1)
    ax.set(xlabel="Injected current (pA)", ylabel="WT minus SCA3 difference (nats)",
           title="Primary 4-ms shuffle-centered PI by current")
    fig.tight_layout()
    fig.savefig(figures / "primary_currentwise_effect.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    sub = stability[stability.lag_ms == primary_lag]
    for kind, values in sub.groupby("surrogate_type"):
        values = values.sort_values("n_surrogates")
        ax.plot(values.n_surrogates, values.observed_expected_direction_AUC_difference_WT_minus_SCA3,
                marker="o", label=kind)
    ax.set(xlabel="Surrogates per sweep", ylabel="WT minus SCA3 cell-AUC difference (nats)",
           title="Primary-lag surrogate-count stability")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / "surrogate_count_stability.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    axes[0].hist(fidelity.spectral_amplitude_nrmse, bins=40, color="#4c78a8")
    axes[0].set(title="IAAFT spectrum fidelity", xlabel="Spectral-amplitude NRMSE", ylabel="Surrogates")
    axes[1].hist(fidelity.acf_rmse, bins=40, color="#72b7b2")
    axes[1].set(title="IAAFT autocorrelation fidelity", xlabel="ACF RMSE", ylabel="Surrogates")
    fig.tight_layout()
    fig.savefig(figures / "iaaft_fidelity.png", dpi=180)
    plt.close(fig)


def _write_summary(output, config, auc_tests, fidelity_overall, reproduction):
    primary_lag = int(config["surrogates"]["primary_code_pair_lag_ms"])
    primary = auc_tests[(auc_tests.surrogate_family == "shuffle") & (auc_tests.lag_ms == primary_lag)].iloc[0]
    iaaft = auc_tests[(auc_tests.surrogate_family == "iaaft") & (auc_tests.lag_ms == primary_lag)].iloc[0]
    text = """# NeuroThermo v0.7.2 — non-overlapping ordinal PI and validated IAAFT

## Frozen design

All 13 WT and 7 SCA3 cells are retained. The cell is the independent unit; animal-level inference is absent. Ordinal order is 4 at 1 ms delay. The primary code-pair lag is 4 ms, so paired ordinal windows share no voltage samples. Secondary lags are 8, 16 and 32 ms.

## Primary endpoint

The primary metric is 4-ms non-overlapping ordinal PI minus the median shuffled-surrogate PI. WT-minus-SCA3 cell-AUC difference: {primary_effect:.8g} nats; exact expected-direction p={primary_p:.8g}. This is the prespecified unadjusted primary test.

## IAAFT sensitivity

The 4-ms IAAFT-centered WT-minus-SCA3 cell-AUC difference is {iaaft_effect:.8g} nats; exact unadjusted p={iaaft_p:.8g}; maxT p across the four lags={iaaft_maxt:.8g}.

All IAAFT fidelity gates passed: {fidelity_pass}. The maximum adjacent-PI reconstruction error relative to v0.3.1 was {reproduction_error:.3g} nats.

## Interpretation boundary

`surrogate_centered_PI` is a signed deviation from a surrogate null, not a non-negative quantity of extra information. Survival of the 4-ms result removes the direct shared-sample artifact of adjacent ordinal windows. Survival against validated IAAFT indicates a difference not reproduced by the amplitude distribution and approximately matched linear spectrum/autocorrelation. It does not prove entropy production, thermodynamic irreversibility, disease time, causal mechanism or a phase transition.
""".format(
        primary_effect=float(primary.observed_expected_direction_AUC_difference_WT_minus_SCA3),
        primary_p=float(primary.exact_p_AUC_expected_direction_unadjusted),
        iaaft_effect=float(iaaft.observed_expected_direction_AUC_difference_WT_minus_SCA3),
        iaaft_p=float(iaaft.exact_p_AUC_expected_direction_unadjusted),
        iaaft_maxt=float(iaaft.exact_p_AUC_expected_direction_maxT_across_lags),
        fidelity_pass=bool(fidelity_overall.passed.all()),
        reproduction_error=float(reproduction.absolute_reproduction_error_nats.max()),
    )
    (output / "RUN_SUMMARY.md").write_text(text, encoding="utf-8")


def run_pipeline(config, upstream: Path, raw_root: Path, output: Path) -> None:
    checks = validate_inputs(config, upstream, raw_root)
    if not bool(checks.passed.all()):
        raise ValueError("Input validation failed:\n" + checks.to_string(index=False))
    output.mkdir(parents=True, exist_ok=True)
    checks.to_csv(output / "input_validation.csv", index=False)

    domain = _read_domain(config, upstream)
    summary, long_table, arrays, fidelity, lag_ms_values = _extract_and_compute(config, domain, raw_root)
    tolerance = float(config["surrogates"]["reproduction_tolerance_nats"])
    reproduction = summary[[
        "group", "cell_id", "current_pA", "sweep_index", "source_path_upstream",
        "resolved_source_path", "predictive_information_upstream_nats", "adjacent_PI_recomputed_nats",
        "absolute_reproduction_error_nats", "upstream_stationary_samples",
        "raw_stationary_samples", "resampled_stationary_samples", "target_dt_ms",
    ]].copy()
    reproduction["within_tolerance"] = reproduction.absolute_reproduction_error_nats <= tolerance
    reproduction.to_csv(output / "raw_trace_reproduction_audit.csv", index=False)
    if not bool(reproduction.within_tolerance.all()):
        raise ValueError("Raw adjacent-PI reproduction failed:\n" +
                         reproduction[~reproduction.within_tolerance].to_string(index=False))

    summary.to_csv(output / "nonoverlap_PI_surrogate_summary.csv", index=False)
    long_table.to_csv(output / "nonoverlap_PI_surrogates_long.csv.gz", index=False, compression="gzip")
    fidelity.to_csv(output / "iaaft_fidelity_long.csv.gz", index=False, compression="gzip")
    fidelity_by_sweep, fidelity_overall = _fidelity_summary(config, fidelity)
    fidelity_by_sweep.to_csv(output / "iaaft_fidelity_by_sweep.csv", index=False)
    fidelity_overall.to_csv(output / "iaaft_fidelity_overall.csv", index=False)
    if not bool(fidelity_overall.passed.all()):
        raise ValueError("IAAFT fidelity gate failed:\n" + fidelity_overall.to_string(index=False))

    groups = [str(x) for x in config["analysis"]["groups"]]
    currents = np.asarray(config["analysis"]["currents_pA"], dtype=int)
    cell_ids, labels = _ordered_cells(domain, groups)
    masks = label_masks(len(cell_ids), int(labels.sum()), int(config["analysis"]["exact_max_labelings"]))
    auc_tests, current_tests, cell_aucs = exact_lag_family_tests(
        summary, _families(lag_ms_values), cell_ids, labels, currents, masks, expected_direction=-1.0
    )
    auc_tests.to_csv(output / "PI_lag_cell_AUC_exact_tests.csv", index=False)
    current_tests.to_csv(output / "PI_lag_currentwise_exact_tests.csv", index=False)
    cell_aucs.to_csv(output / "PI_cell_AUC_values.csv", index=False)
    count_stability = _count_stability(
        config, summary, arrays, cell_ids, labels, currents, masks, lag_ms_values
    )
    count_stability.to_csv(output / "surrogate_count_stability.csv", index=False)
    curves = _group_curves(summary, lag_ms_values)
    curves.to_csv(output / "group_dynamic_curves.csv", index=False)

    lag_map = _lag_samples(config)
    order = int(config["surrogates"]["permutation_order"])
    delay_samples = int(round(float(config["surrogates"]["permutation_delay_ms"]) /
                              float(config["surrogates"]["resample_dt_ms"])))
    rule = {
        "pipeline_version": __version__, "analysis_unit": "cell", "animal_level_inference": False,
        "cohort": {"WT": int((~labels).sum()), "SCA3": int(labels.sum())},
        "currents_pA": [int(x) for x in currents],
        "ordinal_PI": {
            "permutation_order": order,
            "permutation_delay_ms": float(config["surrogates"]["permutation_delay_ms"]),
            "ordinal_window_width_ms": float(((order - 1) * delay_samples + 1) *
                                                float(config["surrogates"]["resample_dt_ms"])),
            "code_pair_lags_ms": [int(x) for x in lag_ms_values],
            "code_pair_lags_samples": {str(k): int(v) for k, v in lag_map.items()},
            "primary_code_pair_lag_ms": int(config["surrogates"]["primary_code_pair_lag_ms"]),
            "primary_windows_overlap": False,
            "estimator": "empirical discrete mutual information",
        },
        "primary_endpoint": {
            "metric": "shuffle-centered non-overlap PI at 4 ms",
            "summary": "normalized cell AUC over 100-600 pA",
            "test": "exact cell-label expected-direction test",
            "multiplicity": "none; single prespecified endpoint",
        },
        "secondary": {
            "null": "IAAFT", "lags_ms": [8, 16, 32],
            "multiplicity": "maxT across lags for AUC and across lag-current combinations for currentwise tests",
        },
        "surrogate_tail_diagnostics": ["lower", "upper", "two_sided_about_surrogate_median"],
        "IAAFT": {
            "n_per_sweep": int(config["surrogates"]["n_iaaft"]),
            "max_iterations": int(config["surrogates"]["iaaft_max_iterations"]),
            "improvement_tolerance": float(config["surrogates"]["iaaft_improvement_tolerance"]),
            "patience": int(config["surrogates"]["iaaft_patience"]),
            "fidelity_outputs": ["spectral_amplitude_nrmse", "log_psd_rmse", "acf_rmse", "amplitude_sorted_exact"],
        },
        "exact_inference": {"cell_labelings": int(len(masks))},
        "interpretation": "temporal phenotype validation; not disease time, entropy production, or thermodynamic phase transition",
    }
    _json_dump(output / "frozen_nonoverlap_iaaft_rule.json", rule)
    _plots(output, config, curves, auc_tests, current_tests, count_stability, fidelity)
    _write_summary(output, config, auc_tests, fidelity_overall, reproduction)

    raw_paths = sorted(set(Path(x) for x in reproduction.resolved_source_path))
    input_hashes = {"v0.3.1/" + item: _sha256(upstream / item) for item in REQUIRED_UPSTREAM}
    raw_hashes = {str(path): _sha256(path) for path in raw_paths}
    output_files = sorted(x for x in output.rglob("*") if x.is_file() and x.name != "analysis_manifest.json")
    manifest = {
        "pipeline_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_minimum": "3.9", "upstream_dir": str(upstream), "raw_root": str(raw_root),
        "output_dir": str(output), "input_sha256": input_hashes, "raw_trace_sha256": raw_hashes,
        "output_sha256": {str(path.relative_to(output)): _sha256(path) for path in output_files},
        "n_cells": int(len(cell_ids)), "n_WT_cells": int((~labels).sum()),
        "n_SCA3_cells": int(labels.sum()), "n_currents": int(len(currents)),
        "n_exact_labelings": int(len(masks)), "animal_level_inference": False,
        "primary_code_pair_lag_ms": int(config["surrogates"]["primary_code_pair_lag_ms"]),
        "IAAFT_fidelity_all_passed": bool(fidelity_overall.passed.all()),
    }
    _json_dump(output / "analysis_manifest.json", manifest)
