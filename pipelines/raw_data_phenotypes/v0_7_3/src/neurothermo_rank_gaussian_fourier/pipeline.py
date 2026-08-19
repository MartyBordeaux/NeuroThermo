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
from .config import REQUIRED_UPSTREAM, REQUIRED_V072
from .metrics import (
    fourier_fidelity,
    fourier_phase_surrogate,
    ordinal_pi_lags,
    rank_gaussianize,
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


def validate_inputs(config: Mapping[str, Any], upstream: Path, v072: Path, raw_root: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    missing_upstream = [str(upstream / x) for x in REQUIRED_UPSTREAM if not (upstream / x).is_file()]
    missing_v072 = [str(v072 / x) for x in REQUIRED_V072 if not (v072 / x).is_file()]
    rows.append({"check": "required_v031_files", "passed": not missing_upstream,
                 "detail": "none" if not missing_upstream else "; ".join(missing_upstream)})
    rows.append({"check": "required_v072_files", "passed": not missing_v072,
                 "detail": "none" if not missing_v072 else "; ".join(missing_v072)})
    if missing_upstream or missing_v072:
        return pd.DataFrame(rows)

    domain = _read_domain(config, upstream)
    v072_summary = pd.read_csv(v072 / "nonoverlap_PI_surrogate_summary.csv")
    required_domain = {
        "group", "cell_id", "source_path", "sweep_index", "current_pA",
        "stim_start_s", "stim_end_s", "stationary_samples",
    }
    required_v072_columns = {"group", "cell_id", "current_pA", "n_shuffle"}
    for lag in config["surrogates"]["code_pair_lags_ms"]:
        required_v072_columns.update({
            "observed_PI_lag_{}ms_nats".format(int(lag)),
            "shuffle_lag_{}ms_centered_PI_nats".format(int(lag)),
            "shuffle_lag_{}ms_median_nats".format(int(lag)),
        })
    absent_domain = sorted(required_domain - set(domain.columns))
    absent_v072 = sorted(required_v072_columns - set(v072_summary.columns))
    rows.append({"check": "required_v031_columns", "passed": not absent_domain,
                 "detail": "none" if not absent_domain else ", ".join(absent_domain)})
    rows.append({"check": "required_v072_columns", "passed": not absent_v072,
                 "detail": "none" if not absent_v072 else ", ".join(absent_v072)})
    if absent_domain or absent_v072:
        return pd.DataFrame(rows)

    groups = [str(x) for x in config["analysis"]["groups"]]
    currents = [int(x) for x in config["analysis"]["currents_pA"]]
    expected = config["analysis"]["expected_cells"]
    identities = domain[["group", "cell_id"]].drop_duplicates()
    counts = identities.groupby("group").size().to_dict()
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
    domain_keys = set(map(tuple, domain[["group", "cell_id", "current_pA"]].itertuples(index=False, name=None)))
    v072_keys = set(map(tuple, v072_summary[["group", "cell_id", "current_pA"]].itertuples(index=False, name=None)))
    v072_duplicates = int(v072_summary.duplicated(["cell_id", "current_pA"]).sum())
    rows.append({"check": "one_v072_row_per_cell_current", "passed": v072_duplicates == 0,
                 "detail": "duplicates={}; rows={}".format(v072_duplicates, len(v072_summary))})
    rows.append({"check": "v072_cell_current_identity", "passed": domain_keys == v072_keys,
                 "detail": "v031={}; v072={}; symmetric_difference={}".format(
                     len(domain_keys), len(v072_keys), len(domain_keys ^ v072_keys)
                 )})

    reproduction = pd.read_csv(v072 / "raw_trace_reproduction_audit.csv")
    rows.append({"check": "v072_raw_reproduction", "passed": bool(reproduction.within_tolerance.all()),
                 "detail": "passed={}/{}; max_error={:.6g}".format(
                     int(reproduction.within_tolerance.sum()), len(reproduction),
                     float(reproduction.absolute_reproduction_error_nats.max())
                 )})
    expected_shuffle = int(config["surrogates"]["n_imported_shuffle"])
    shuffle_long = pd.read_csv(
        v072 / "nonoverlap_PI_surrogates_long.csv.gz",
        usecols=["cell_id", "current_pA", "surrogate_type", "surrogate_index"],
    )
    shuffle_long = shuffle_long[shuffle_long.surrogate_type == "shuffle"]
    count_values = shuffle_long.groupby(["cell_id", "current_pA"]).size()
    duplicate_shuffle = int(shuffle_long.duplicated(["cell_id", "current_pA", "surrogate_index"]).sum())
    rows.append({"check": "complete_imported_shuffle_surrogates",
                 "passed": len(count_values) == len(domain) and bool((count_values == expected_shuffle).all()),
                 "detail": "sweeps={}; counts={}".format(len(count_values), count_values.value_counts().to_dict())})
    rows.append({"check": "unique_imported_shuffle_keys", "passed": duplicate_shuffle == 0,
                 "detail": "duplicates={}".format(duplicate_shuffle)})

    iaaft = pd.read_csv(v072 / "iaaft_fidelity_overall.csv")
    failed = iaaft.loc[~iaaft.passed.astype(bool), "check"].astype(str).tolist()
    rows.append({"check": "v072_iaaft_excluded_from_inference", "passed": True,
                 "detail": "failed_gates={}".format(json.dumps(failed))})

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
    width = (order - 1) * delay_samples + 1
    primary_lag = int(config["surrogates"]["primary_code_pair_lag_ms"])
    rows.append({"check": "primary_ordinal_windows_do_not_overlap",
                 "passed": primary_lag in lag_map and lag_map.get(primary_lag, 0) >= width,
                 "detail": "ordinal_width_samples={}; primary_code_lag_samples={}".format(width, lag_map.get(primary_lag))})
    rows.append({"check": "cell_level_only", "passed": "animal_id" not in domain.columns,
                 "detail": "all 20 cells retained; no animal-level inference"})
    return pd.DataFrame(rows)


def _seed(master: int, cell_id: str, current_pA: int) -> int:
    key = "{}|{}|{}|rank_gaussian_fourier".format(master, cell_id, current_pA).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "little")


def _fourier_worker(payload):
    cell_id, current_pA, gaussianized, order, delay, code_lags, n_fourier, master_seed, maximum_acf_lag = payload
    rng = np.random.default_rng(_seed(master_seed, cell_id, current_pA))
    pi_values = np.empty((n_fourier, len(code_lags)), dtype=float)
    diagnostics = []
    for index in range(n_fourier):
        surrogate = fourier_phase_surrogate(gaussianized, rng)
        pi_values[index] = ordinal_pi_lags(surrogate, order, delay, code_lags)
        diagnostics.append(fourier_fidelity(gaussianized, surrogate, maximum_acf_lag))
    return pi_values, diagnostics


def _surrogate_statistics(observed: float, values: np.ndarray, prefix: str) -> Dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    median = float(np.median(finite))
    mean = float(np.mean(finite))
    sd = float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan
    q = np.quantile(finite, [0.025, 0.25, 0.75, 0.975])
    distance = abs(float(observed) - median)
    return {
        prefix + "_median_nats": median,
        prefix + "_mean_nats": mean,
        prefix + "_sd_nats": sd,
        prefix + "_q025_nats": float(q[0]),
        prefix + "_q25_nats": float(q[1]),
        prefix + "_q75_nats": float(q[2]),
        prefix + "_q975_nats": float(q[3]),
        prefix + "_centered_PI_nats": float(observed - median),
        prefix + "_z": float((observed - mean) / sd) if np.isfinite(sd) and sd > 1e-15 else np.nan,
        prefix + "_p_lower": float((1 + np.sum(finite <= observed)) / (1 + len(finite))),
        prefix + "_p_upper": float((1 + np.sum(finite >= observed)) / (1 + len(finite))),
        prefix + "_p_two_sided": float((1 + np.sum(np.abs(finite - median) >= distance)) / (1 + len(finite))),
    }


def _extract_and_compute(config, domain, v072_summary, raw_root):
    dt_ms = float(config["surrogates"]["resample_dt_ms"])
    order = int(config["surrogates"]["permutation_order"])
    delay = int(round(float(config["surrogates"]["permutation_delay_ms"]) / dt_ms))
    lag_map = _lag_samples(config)
    lag_ms_values = list(lag_map)
    code_lags = [lag_map[x] for x in lag_ms_values]
    n_fourier = int(config["surrogates"]["n_fourier"])
    maximum_acf_lag = int(round(float(config["surrogates"]["linear_acf_max_lag_ms"]) / dt_ms))
    key_columns = ["group", "cell_id", "current_pA"]
    imported = v072_summary.set_index(key_columns)

    bases = []
    gaussianized_traces = []
    audit_rows = []
    for _, row in domain.iterrows():
        key = (str(row.group), str(row.cell_id), int(row.current_pA))
        old = imported.loc[key]
        values, metadata = stationary_trace(row, raw_root, config)
        gaussianized = rank_gaussianize(values)
        gaussian_pi = ordinal_pi_lags(gaussianized, order, delay, code_lags)
        base = {
            "group": key[0], "cell_id": key[1], "current_pA": key[2],
            "sweep_index": int(row.sweep_index), "source_path_upstream": str(row.source_path),
        }
        base.update(metadata)
        for lag_ms, value in zip(lag_ms_values, gaussian_pi):
            imported_observed = float(old["observed_PI_lag_{}ms_nats".format(lag_ms)])
            base["observed_PI_lag_{}ms_nats".format(lag_ms)] = imported_observed
            audit_rows.append({
                "group": key[0], "cell_id": key[1], "current_pA": key[2], "lag_ms": int(lag_ms),
                "imported_observed_PI_nats": imported_observed,
                "rank_gaussian_observed_PI_nats": float(value),
                "absolute_difference_nats": abs(float(value) - imported_observed),
            })
        for column in v072_summary.columns:
            if column.startswith("shuffle_lag_"):
                base[column] = old[column]
        base["n_imported_shuffle"] = int(old.n_shuffle)
        bases.append(base)
        gaussianized_traces.append(gaussianized)

    payloads = [
        (
            str(row.cell_id), int(row.current_pA), values, order, delay, code_lags,
            n_fourier, int(config["analysis"]["random_seed"]), maximum_acf_lag,
        )
        for (_, row), values in zip(domain.iterrows(), gaussianized_traces)
    ]
    workers = int(config.get("runtime", {}).get("workers", 2))
    if workers <= 1:
        results = [_fourier_worker(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_fourier_worker, payloads, chunksize=1))

    summary_rows = []
    long_rows = []
    fidelity_rows = []
    arrays = {}
    diagnostic_columns = [
        "spectral_amplitude_nrmse", "spectral_amplitude_max_abs_error",
        "circular_acf_rmse", "circular_acf_max_abs_error", "linear_acf_rmse",
        "surrogate_mean", "surrogate_sd", "surrogate_skewness", "surrogate_excess_kurtosis",
    ]
    for base, (pi_values, diagnostics) in zip(bases, results):
        summary = dict(base)
        key = (str(base["cell_id"]), int(base["current_pA"]))
        arrays[key] = pi_values
        for lag_index, lag_ms in enumerate(lag_ms_values):
            observed = float(base["observed_PI_lag_{}ms_nats".format(lag_ms)])
            summary.update(_surrogate_statistics(
                observed, pi_values[:, lag_index], "fourier_lag_{}ms".format(lag_ms)
            ))
        diag = pd.DataFrame(diagnostics)
        for column in diagnostic_columns:
            values = diag[column].to_numpy(float)
            summary["fourier_{}_median".format(column)] = float(np.median(values))
            summary["fourier_{}_q95".format(column)] = float(np.quantile(values, .95))
            summary["fourier_{}_max".format(column)] = float(np.max(values))
        summary["n_fourier"] = n_fourier
        summary_rows.append(summary)
        for index in range(n_fourier):
            record = {
                "group": base["group"], "cell_id": base["cell_id"],
                "current_pA": int(base["current_pA"]), "surrogate_index": int(index),
            }
            for lag_index, lag_ms in enumerate(lag_ms_values):
                record["PI_lag_{}ms_nats".format(lag_ms)] = float(pi_values[index, lag_index])
            record.update(diagnostics[index])
            long_rows.append(record)
            fidelity = dict(record)
            for lag_ms in lag_ms_values:
                fidelity.pop("PI_lag_{}ms_nats".format(lag_ms), None)
            fidelity_rows.append(fidelity)
    return (
        pd.DataFrame(summary_rows), pd.DataFrame(long_rows), arrays,
        pd.DataFrame(fidelity_rows), pd.DataFrame(audit_rows), lag_ms_values,
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


def _fidelity_summary(config, fidelity):
    by_sweep_rows = []
    columns = [
        "spectral_amplitude_nrmse", "spectral_amplitude_max_abs_error",
        "circular_acf_rmse", "circular_acf_max_abs_error", "linear_acf_rmse",
        "surrogate_mean", "surrogate_sd", "surrogate_skewness", "surrogate_excess_kurtosis",
    ]
    for keys, sub in fidelity.groupby(["group", "cell_id", "current_pA"], sort=True):
        row = {"group": keys[0], "cell_id": keys[1], "current_pA": int(keys[2]), "n_fourier": int(len(sub))}
        for column in columns:
            values = sub[column].to_numpy(float)
            row[column + "_median"] = float(np.median(values))
            row[column + "_q95"] = float(np.quantile(values, .95))
            row[column + "_max"] = float(np.max(values))
        by_sweep_rows.append(row)
    by_sweep = pd.DataFrame(by_sweep_rows)
    tolerance = float(config["surrogates"]["exact_spectrum_tolerance"])
    maximum_spectrum = float(fidelity.spectral_amplitude_nrmse.max())
    maximum_circular = float(fidelity.circular_acf_max_abs_error.max())
    finite = bool(np.isfinite(fidelity.select_dtypes(include=[np.number]).to_numpy()).all())
    overall = pd.DataFrame([
        {"check": "Fourier_spectral_amplitude_NRMSE_max", "value": maximum_spectrum,
         "threshold": tolerance, "comparison": "<=", "passed": maximum_spectrum <= tolerance},
        {"check": "Fourier_circular_ACF_max_abs_error", "value": maximum_circular,
         "threshold": tolerance, "comparison": "<=", "passed": maximum_circular <= tolerance},
        {"check": "Fourier_all_diagnostics_finite", "value": float(finite),
         "threshold": 1.0, "comparison": "==", "passed": finite},
        {"check": "Fourier_linear_ACF_RMSE_q95_informational", "value": float(fidelity.linear_acf_rmse.quantile(.95)),
         "threshold": np.nan, "comparison": "report_only", "passed": True},
    ])
    return by_sweep, overall


def _families(lags):
    return {
        "shuffle": {int(lag): "shuffle_lag_{}ms_centered_PI_nats".format(lag) for lag in lags},
        "fourier": {int(lag): "fourier_lag_{}ms_centered_PI_nats".format(lag) for lag in lags},
    }


def _leave_one_cell_out(summary, families, cell_ids, labels, currents, maximum):
    rows = []
    for family, lag_columns in families.items():
        for lag, column in sorted(lag_columns.items()):
            matrix = summary.pivot(index="cell_id", columns="current_pA", values=column).reindex(
                index=cell_ids, columns=currents
            ).to_numpy(float)
            auc_values = cell_auc(matrix, currents)
            for index, cell in enumerate(cell_ids):
                keep = np.arange(len(cell_ids)) != index
                kept_labels = labels[keep]
                masks = label_masks(int(keep.sum()), int(kept_labels.sum()), maximum)
                observed, p_expected, p_two, _ = exact_cell_difference(
                    auc_values[keep], masks, kept_labels, -1.0
                )
                rows.append({
                    "surrogate_family": family, "lag_ms": int(lag),
                    "dropped_cell_id": str(cell),
                    "dropped_group": "SCA3" if bool(labels[index]) else "WT",
                    "observed_expected_direction_AUC_difference_WT_minus_SCA3": -float(observed),
                    "exact_p_expected_direction": float(p_expected),
                    "exact_p_two_sided": float(p_two),
                    "n_exact_labelings": int(len(masks)),
                })
    return pd.DataFrame(rows)


def _count_stability(config, summary, imported_long, fourier_arrays, cell_ids, labels, currents, masks, lags):
    rows = []
    requested = sorted(set(int(x) for x in config["surrogates"]["count_sensitivity"]))
    imported = imported_long[imported_long.surrogate_type == "shuffle"].sort_values(
        ["cell_id", "current_pA", "surrogate_index"]
    )
    shuffle_arrays = {}
    for (cell, current), sub in imported.groupby(["cell_id", "current_pA"], sort=False):
        shuffle_arrays[(str(cell), int(current))] = {
            int(lag): sub["PI_lag_{}ms_nats".format(lag)].to_numpy(float) for lag in lags
        }
    for family, maximum, source in (
        ("shuffle", int(config["surrogates"]["n_imported_shuffle"]), shuffle_arrays),
        ("fourier", int(config["surrogates"]["n_fourier"]), fourier_arrays),
    ):
        for count in [x for x in requested if x <= maximum]:
            for lag_index, lag in enumerate(lags):
                observed_column = "observed_PI_lag_{}ms_nats".format(lag)
                work = summary[["cell_id", "current_pA", observed_column]].copy()
                medians = []
                for row in work.itertuples():
                    key = (str(row.cell_id), int(row.current_pA))
                    values = source[key][int(lag)][:count] if family == "shuffle" else source[key][:count, lag_index]
                    medians.append(float(np.median(values)))
                work["centered"] = work[observed_column].to_numpy(float) - np.asarray(medians)
                matrix = work.pivot(index="cell_id", columns="current_pA", values="centered").reindex(
                    index=cell_ids, columns=currents
                ).to_numpy(float)
                auc_values = cell_auc(matrix, currents)
                observed, p_expected, p_two, _ = exact_cell_difference(auc_values, masks, labels, -1.0)
                rows.append({
                    "surrogate_family": family, "lag_ms": int(lag), "n_surrogates": int(count),
                    "observed_expected_direction_AUC_difference_WT_minus_SCA3": -float(observed),
                    "exact_p_AUC_expected_direction": float(p_expected),
                    "exact_p_AUC_two_sided": float(p_two),
                })
    return pd.DataFrame(rows)


def _group_curves(summary, lags):
    rows = []
    metrics = []
    for lag in lags:
        metrics.extend([
            "observed_PI_lag_{}ms_nats".format(lag),
            "shuffle_lag_{}ms_centered_PI_nats".format(lag),
            "fourier_lag_{}ms_centered_PI_nats".format(lag),
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
               "Primary: shuffle-centered non-overlap PI", groups)
    _plot_band(axes[1], curves, "fourier_lag_{}ms_centered_PI_nats".format(primary_lag),
               "Exact-spectrum rank-Gaussian Fourier null", groups)
    fig.tight_layout(); fig.savefig(figures / "primary_and_fourier_PI_curves.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, family in zip(axes, ["shuffle", "fourier"]):
        sub = auc_tests[auc_tests.surrogate_family == family].sort_values("lag_ms")
        ax.plot(sub.lag_ms, sub.observed_expected_direction_AUC_difference_WT_minus_SCA3,
                marker="o", color="#4c78a8")
        significant = sub.exact_p_AUC_expected_direction_maxT_across_lags <= .05
        ax.scatter(sub.loc[significant, "lag_ms"],
                   sub.loc[significant, "observed_expected_direction_AUC_difference_WT_minus_SCA3"],
                   color="#be2334", s=65, zorder=3)
        ax.axhline(0.0, color="black", lw=1)
        ax.set(title="{}-centered".format(family.upper()), xlabel="Code-pair lag (ms)")
    axes[0].set_ylabel("WT minus SCA3 cell-AUC difference (nats)")
    fig.tight_layout(); fig.savefig(figures / "lag_sensitivity_cell_AUC.png", dpi=180); plt.close(fig)

    primary = current_tests[(current_tests.surrogate_family == "shuffle") &
                            (current_tests.lag_ms == primary_lag)].sort_values("current_pA")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(primary.current_pA, primary.observed_expected_direction_difference_WT_minus_SCA3,
            marker="o", color="#4c78a8")
    significant = primary.exact_p_expected_direction_maxT_across_lags_and_currents <= .05
    ax.scatter(primary.loc[significant, "current_pA"],
               primary.loc[significant, "observed_expected_direction_difference_WT_minus_SCA3"],
               color="#be2334", s=65, zorder=3)
    ax.axhline(0.0, color="black", lw=1)
    ax.set(xlabel="Injected current (pA)", ylabel="WT minus SCA3 difference (nats)",
           title="Primary 4-ms shuffle-centered PI")
    fig.tight_layout(); fig.savefig(figures / "primary_currentwise_effect.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for family, sub in stability[stability.lag_ms == primary_lag].groupby("surrogate_family"):
        sub = sub.sort_values("n_surrogates")
        ax.plot(sub.n_surrogates, sub.observed_expected_direction_AUC_difference_WT_minus_SCA3,
                marker="o", label=family)
    ax.set(xlabel="Surrogates per sweep", ylabel="WT minus SCA3 cell-AUC difference (nats)",
           title="Primary-lag surrogate-count stability")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(figures / "surrogate_count_stability.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    axes[0].hist(fidelity.spectral_amplitude_nrmse, bins=40, color="#4c78a8")
    axes[0].set(title="Exact Fourier-spectrum fidelity", xlabel="Spectral-amplitude NRMSE", ylabel="Surrogates")
    axes[1].hist(fidelity.linear_acf_rmse, bins=40, color="#72b7b2")
    axes[1].set(title="Finite-window linear ACF diagnostic", xlabel="Linear ACF RMSE", ylabel="Surrogates")
    fig.tight_layout(); fig.savefig(figures / "fourier_fidelity.png", dpi=180); plt.close(fig)


def _write_summary(output, config, auc_tests, fidelity_overall, invariance):
    primary_lag = int(config["surrogates"]["primary_code_pair_lag_ms"])
    primary = auc_tests[(auc_tests.surrogate_family == "shuffle") & (auc_tests.lag_ms == primary_lag)].iloc[0]
    fourier = auc_tests[(auc_tests.surrogate_family == "fourier") & (auc_tests.lag_ms == primary_lag)].iloc[0]
    text = """# NeuroThermo v0.7.3 — rank-Gaussian exact-spectrum Fourier validation

## Frozen design

All 13 WT and 7 SCA3 cells are retained. The cell is the independent unit; animal-level inference is absent. The primary endpoint remains the imported v0.7.2 shuffle-centered non-overlapping PI at 4 ms. Failed v0.7.2 IAAFT results are excluded from inference.

## Primary endpoint

WT-minus-SCA3 cell-AUC difference: {primary_effect:.8g} nats; exact expected-direction p={primary_p:.8g}; maxT p across lags={primary_maxt:.8g}.

## Exact-spectrum sensitivity

The voltage trace is transformed monotonically to empirical normal scores. This leaves every ordinal code unchanged. Fourier phases are randomized while the complete rFFT magnitude spectrum of the rank-Gaussian trace is preserved.

Fourier-centered 4-ms WT-minus-SCA3 cell-AUC difference: {fourier_effect:.8g} nats; exact expected-direction p={fourier_p:.8g}; maxT p across lags={fourier_maxt:.8g}.

All exact-spectrum fidelity gates passed: {fidelity_pass}. Maximum ordinal-PI difference introduced by rank Gaussianization: {invariance_error:.3g} nats.

## Interpretation boundary

Survival against this null identifies ordinal temporal structure not reproduced by a linear Gaussian-copula process with the same rank-Gaussianized spectrum. It does not prove entropy production, thermodynamic irreversibility, causal mechanism, disease time or a phase transition.
""".format(
        primary_effect=float(primary.observed_expected_direction_AUC_difference_WT_minus_SCA3),
        primary_p=float(primary.exact_p_AUC_expected_direction_unadjusted),
        primary_maxt=float(primary.exact_p_AUC_expected_direction_maxT_across_lags),
        fourier_effect=float(fourier.observed_expected_direction_AUC_difference_WT_minus_SCA3),
        fourier_p=float(fourier.exact_p_AUC_expected_direction_unadjusted),
        fourier_maxt=float(fourier.exact_p_AUC_expected_direction_maxT_across_lags),
        fidelity_pass=bool(fidelity_overall.passed.all()),
        invariance_error=float(invariance.absolute_difference_nats.max()),
    )
    (output / "RUN_SUMMARY.md").write_text(text, encoding="utf-8")


def run_pipeline(config, upstream: Path, v072: Path, raw_root: Path, output: Path) -> None:
    checks = validate_inputs(config, upstream, v072, raw_root)
    if not bool(checks.passed.all()):
        raise ValueError("Input validation failed:\n" + checks.to_string(index=False))
    output.mkdir(parents=True, exist_ok=True)
    checks.to_csv(output / "input_validation.csv", index=False)
    pd.read_csv(v072 / "iaaft_fidelity_overall.csv").to_csv(
        output / "excluded_v072_iaaft_fidelity.csv", index=False
    )

    domain = _read_domain(config, upstream)
    v072_summary = pd.read_csv(v072 / "nonoverlap_PI_surrogate_summary.csv")
    summary, fourier_long, fourier_arrays, fidelity, invariance, lags = _extract_and_compute(
        config, domain, v072_summary, raw_root
    )
    invariance_tolerance = float(config["surrogates"]["ordinal_invariance_tolerance_nats"])
    invariance["within_tolerance"] = invariance.absolute_difference_nats <= invariance_tolerance
    invariance.to_csv(output / "rank_gaussian_invariance_audit.csv", index=False)
    if not bool(invariance.within_tolerance.all()):
        raise ValueError("Rank-Gaussian ordinal invariance failed:\n" +
                         invariance[~invariance.within_tolerance].to_string(index=False))

    summary.to_csv(output / "rank_gaussian_fourier_summary.csv", index=False)
    fourier_long.to_csv(output / "fourier_surrogates_long.csv.gz", index=False, compression="gzip")
    fidelity.to_csv(output / "fourier_fidelity_long.csv.gz", index=False, compression="gzip")
    fidelity_by_sweep, fidelity_overall = _fidelity_summary(config, fidelity)
    fidelity_by_sweep.to_csv(output / "fourier_fidelity_by_sweep.csv", index=False)
    fidelity_overall.to_csv(output / "fourier_fidelity_overall.csv", index=False)
    if not bool(fidelity_overall.passed.all()):
        raise ValueError("Exact-spectrum Fourier fidelity gate failed:\n" + fidelity_overall.to_string(index=False))

    groups = [str(x) for x in config["analysis"]["groups"]]
    currents = np.asarray(config["analysis"]["currents_pA"], dtype=int)
    cell_ids, labels = _ordered_cells(domain, groups)
    maximum = int(config["analysis"]["exact_max_labelings"])
    masks = label_masks(len(cell_ids), int(labels.sum()), maximum)
    families = _families(lags)
    auc_tests, current_tests, cell_aucs = exact_lag_family_tests(
        summary, families, cell_ids, labels, currents, masks, expected_direction=-1.0
    )
    auc_tests.to_csv(output / "PI_lag_cell_AUC_exact_tests.csv", index=False)
    current_tests.to_csv(output / "PI_lag_currentwise_exact_tests.csv", index=False)
    cell_aucs.to_csv(output / "PI_cell_AUC_values.csv", index=False)
    leave_one_out = _leave_one_cell_out(
        summary, families, cell_ids, labels, currents, maximum
    )
    leave_one_out.to_csv(output / "leave_one_cell_out_robustness.csv", index=False)

    imported_long = pd.read_csv(v072 / "nonoverlap_PI_surrogates_long.csv.gz")
    stability = _count_stability(
        config, summary, imported_long, fourier_arrays, cell_ids, labels, currents, masks, lags
    )
    stability.to_csv(output / "surrogate_count_stability.csv", index=False)
    curves = _group_curves(summary, lags)
    curves.to_csv(output / "group_dynamic_curves.csv", index=False)

    rule = {
        "pipeline_version": __version__, "analysis_unit": "cell", "animal_level_inference": False,
        "cohort": {"WT": int((~labels).sum()), "SCA3": int(labels.sum())},
        "currents_pA": [int(x) for x in currents],
        "primary_endpoint": {
            "source": "frozen complete v0.7.2 shuffle surrogates",
            "metric": "shuffle-centered non-overlapping ordinal PI at 4 ms",
            "test": "exact cell-label expected-direction AUC test",
        },
        "secondary_null": {
            "name": "rank-Gaussian exact-spectrum Fourier phase randomization",
            "monotonic_normal_score_transform": True,
            "ordinal_sequence_invariant": True,
            "exact_spectrum_space": "rank-Gaussianized stationary voltage",
            "amplitude_distribution_preserved": False,
            "amplitude_distribution_role": "not constrained; ordinal metrics are invariant under monotonic transforms",
            "n_per_sweep": int(config["surrogates"]["n_fourier"]),
        },
        "lags_ms": [int(x) for x in lags],
        "multiplicity": {
            "AUC": "maxT across four lags within surrogate family",
            "currentwise": "maxT across all lag-current combinations within surrogate family",
        },
        "excluded_analysis": "v0.7.2 IAAFT; failed full-spectrum fidelity and convergence gates",
        "exact_inference": {"cell_labelings": int(len(masks))},
        "interpretation": "short-timescale ordinal temporal phenotype; not disease time, entropy production, or phase transition",
    }
    _json_dump(output / "frozen_rank_gaussian_fourier_rule.json", rule)
    _plots(output, config, curves, auc_tests, current_tests, stability, fidelity)
    _write_summary(output, config, auc_tests, fidelity_overall, invariance)

    raw_paths = sorted(set(Path(x) for x in summary.resolved_source_path))
    input_hashes = {"v0.3.1/" + item: _sha256(upstream / item) for item in REQUIRED_UPSTREAM}
    input_hashes.update({"v0.7.2/" + item: _sha256(v072 / item) for item in REQUIRED_V072})
    raw_hashes = {str(path): _sha256(path) for path in raw_paths}
    output_files = sorted(x for x in output.rglob("*") if x.is_file() and x.name != "analysis_manifest.json")
    manifest = {
        "pipeline_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_minimum": "3.9", "upstream_dir": str(upstream), "v072_results_dir": str(v072),
        "raw_root": str(raw_root), "output_dir": str(output),
        "input_sha256": input_hashes, "raw_trace_sha256": raw_hashes,
        "output_sha256": {str(path.relative_to(output)): _sha256(path) for path in output_files},
        "n_cells": int(len(cell_ids)), "n_WT_cells": int((~labels).sum()),
        "n_SCA3_cells": int(labels.sum()), "n_currents": int(len(currents)),
        "n_exact_labelings": int(len(masks)), "animal_level_inference": False,
        "primary_code_pair_lag_ms": int(config["surrogates"]["primary_code_pair_lag_ms"]),
        "fourier_fidelity_all_passed": bool(fidelity_overall.passed.all()),
        "v072_iaaft_inference_excluded": True,
    }
    _json_dump(output / "analysis_manifest.json", manifest)
