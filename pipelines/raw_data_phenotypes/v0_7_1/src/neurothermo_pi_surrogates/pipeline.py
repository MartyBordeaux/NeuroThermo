from __future__ import annotations

import copy
import hashlib
import json
import math
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from neurothermo_vulnerability import pipeline as engine

from . import __version__
from .config import REQUIRED_UPSTREAM, REQUIRED_V070
from .metrics import aaft_surrogate, ordinal_predictive_information, shuffled_surrogate
from .raw import resolve_source_path, stationary_trace
from .statistics import (
    cell_auc,
    exact_cell_difference,
    exact_fixed_metric_tests,
    label_masks,
    prepare_covariates,
    ridge_crossfit_by_cell_current,
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


def validate_inputs(
    config: Mapping[str, Any], upstream: Path, v070: Path, raw_root: Path
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    missing_upstream = [str(upstream / x) for x in REQUIRED_UPSTREAM if not (upstream / x).is_file()]
    missing_v070 = [str(v070 / x) for x in REQUIRED_V070 if not (v070 / x).is_file()]
    rows.append({"check": "required_v031_files", "passed": not missing_upstream,
                 "detail": "none" if not missing_upstream else "; ".join(missing_upstream)})
    rows.append({"check": "required_v070_files", "passed": not missing_v070,
                 "detail": "none" if not missing_v070 else "; ".join(missing_v070)})
    if missing_upstream or missing_v070:
        return pd.DataFrame(rows)

    domain = _read_domain(config, upstream)
    required = {
        "group", "cell_id", "source_path", "sweep_index", "current_pA",
        "stim_start_s", "stim_end_s", "stationary_samples",
        "predictive_information_nats", "firing_rate_hz", "mean_isi_ms",
        "baseline_noise_mV",
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
    rows.append({"check": "expected_total_rows",
                 "passed": len(domain) == sum(int(expected[x]) for x in groups) * len(currents),
                 "detail": "rows={}".format(len(domain))})

    missing_raw = []
    resolved = set()
    for _, row in domain[["group", "source_path"]].drop_duplicates().iterrows():
        try:
            resolved.add(str(resolve_source_path(str(row.source_path), str(row.group), raw_root)))
        except FileNotFoundError as exc:
            missing_raw.append(str(exc))
    rows.append({"check": "raw_trace_files", "passed": not missing_raw,
                 "detail": "resolved_unique_files={}".format(len(resolved)) if not missing_raw else " | ".join(missing_raw)})

    v070_manifest = json.loads((v070 / "analysis_manifest.json").read_text(encoding="utf-8"))
    rows.append({"check": "v070_version", "passed": str(v070_manifest.get("pipeline_version")) == "0.7.0",
                 "detail": str(v070_manifest.get("pipeline_version"))})
    v070_cells = set(pd.read_csv(v070 / "mode_cell_vulnerability_summary.csv").cell_id.astype(str))
    domain_cells = set(domain.cell_id.astype(str))
    rows.append({"check": "v070_cell_identity", "passed": v070_cells == domain_cells,
                 "detail": "v070={}; upstream={}; symmetric_difference={}".format(
                     len(v070_cells), len(domain_cells), sorted(v070_cells ^ domain_cells)
                 )})
    rows.append({"check": "cell_level_only",
                 "passed": "animal_id" not in domain.columns and not bool(v070_manifest.get("animal_level_inference", True)),
                 "detail": "animal_id_removed=True; v070_animal_level_inference={}".format(v070_manifest.get("animal_level_inference"))})
    return pd.DataFrame(rows)


def _seed(master: int, cell_id: str, current_pA: int, kind: str) -> int:
    key = "{}|{}|{}|{}".format(master, cell_id, current_pA, kind).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "little")


def _surrogate_worker(payload):
    cell_id, current_pA, values, order, delay, n_shuffle, n_aaft, master_seed = payload
    shuffle_rng = np.random.default_rng(_seed(master_seed, cell_id, current_pA, "shuffle"))
    aaft_rng = np.random.default_rng(_seed(master_seed, cell_id, current_pA, "aaft"))
    shuffle_values = np.empty(n_shuffle, dtype=float)
    aaft_values = np.empty(n_aaft, dtype=float)
    for index in range(n_shuffle):
        surrogate = shuffled_surrogate(values, shuffle_rng)
        shuffle_values[index] = ordinal_predictive_information(surrogate, order, delay)
    for index in range(n_aaft):
        surrogate = aaft_surrogate(values, aaft_rng)
        aaft_values[index] = ordinal_predictive_information(surrogate, order, delay)
    return shuffle_values, aaft_values


def _surrogate_columns(observed: float, values: np.ndarray, prefix: str) -> Dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return {prefix + suffix: np.nan for suffix in (
            "_median_nats", "_mean_nats", "_sd_nats", "_q025_nats", "_q25_nats",
            "_q75_nats", "_q975_nats", "_excess_nats", "_z", "_p_upper",
        )}
    median = float(np.median(finite))
    mean = float(np.mean(finite))
    sd = float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan
    quantiles = np.quantile(finite, [0.025, 0.25, 0.75, 0.975])
    return {
        prefix + "_median_nats": median,
        prefix + "_mean_nats": mean,
        prefix + "_sd_nats": sd,
        prefix + "_q025_nats": float(quantiles[0]),
        prefix + "_q25_nats": float(quantiles[1]),
        prefix + "_q75_nats": float(quantiles[2]),
        prefix + "_q975_nats": float(quantiles[3]),
        prefix + "_excess_nats": float(observed - median),
        prefix + "_z": float((observed - mean) / sd) if np.isfinite(sd) and sd > 1e-15 else np.nan,
        prefix + "_p_upper": float((1 + np.sum(finite >= observed)) / (1 + len(finite))),
    }


def _extract_and_compute_surrogates(config, domain: pd.DataFrame, raw_root: Path):
    trace_rows = []
    traces = []
    order = int(config["surrogates"]["permutation_order"])
    delay = max(1, int(round(
        float(config["surrogates"]["permutation_delay_ms"]) /
        float(config["surrogates"]["resample_dt_ms"])
    )))
    for _, row in domain.iterrows():
        values, metadata = stationary_trace(row, raw_root, config)
        observed = ordinal_predictive_information(values, order, delay)
        base = {
            "group": str(row.group), "cell_id": str(row.cell_id),
            "current_pA": int(row.current_pA), "sweep_index": int(row.sweep_index),
            "source_path_upstream": str(row.source_path),
            "predictive_information_upstream_nats": float(row.predictive_information_nats),
            "pi_observed_nats": observed,
            "absolute_reproduction_error_nats": abs(observed - float(row.predictive_information_nats)),
            "upstream_stationary_samples": int(row.stationary_samples),
        }
        base.update(metadata)
        trace_rows.append(base)
        traces.append(values)

    n_shuffle = int(config["surrogates"]["n_shuffle"])
    n_aaft = int(config["surrogates"]["n_aaft"])
    master_seed = int(config["analysis"]["random_seed"])
    payloads = [
        (str(row.cell_id), int(row.current_pA), values, order, delay,
         n_shuffle, n_aaft, master_seed)
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
    for base, (shuffle_values, aaft_values) in zip(trace_rows, results):
        observed = float(base["pi_observed_nats"])
        summary = dict(base)
        summary.update(_surrogate_columns(observed, shuffle_values, "pi_shuffle"))
        summary.update(_surrogate_columns(observed, aaft_values, "pi_aaft"))
        summary["n_shuffle"] = int(len(shuffle_values))
        summary["n_aaft"] = int(len(aaft_values))
        summary_rows.append(summary)
        key = (str(base["cell_id"]), int(base["current_pA"]))
        arrays[key] = {"shuffle": shuffle_values, "aaft": aaft_values}
        for kind, values in (("shuffle", shuffle_values), ("aaft", aaft_values)):
            for index, value in enumerate(values):
                long_rows.append({
                    "group": base["group"], "cell_id": base["cell_id"],
                    "current_pA": base["current_pA"], "surrogate_type": kind,
                    "surrogate_index": int(index), "predictive_information_nats": float(value),
                })
    return pd.DataFrame(summary_rows), pd.DataFrame(long_rows), arrays


def _ordered_cells(domain: pd.DataFrame, groups: Sequence[str]):
    order = domain[["group", "cell_id"]].drop_duplicates().sort_values(
        ["group", "cell_id"],
        key=lambda column: column.map({groups[0]: 0, groups[1]: 1}) if column.name == "group" else column,
    )
    cell_ids = order.cell_id.astype(str).tolist()
    by_cell = order.set_index("cell_id").group.to_dict()
    labels = np.asarray([by_cell[cell] == groups[1] for cell in cell_ids], dtype=bool)
    return cell_ids, labels


def _count_stability(config, summary, arrays, cell_ids, labels, currents, masks):
    rows = []
    requested = sorted(set(int(x) for x in config["surrogates"]["count_sensitivity"]))
    for kind, maximum in (("shuffle", int(config["surrogates"]["n_shuffle"])),
                          ("aaft", int(config["surrogates"]["n_aaft"]))):
        for count in [x for x in requested if x <= maximum]:
            work = summary[["group", "cell_id", "current_pA", "pi_observed_nats"]].copy()
            medians = []
            for _, row in work.iterrows():
                values = arrays[(str(row.cell_id), int(row.current_pA))][kind][:count]
                medians.append(float(np.median(values)))
            metric = "pi_{}_excess_nats".format(kind)
            work[metric] = work.pi_observed_nats.to_numpy(float) - np.asarray(medians)
            matrix = work.pivot(index="cell_id", columns="current_pA", values=metric).reindex(
                index=cell_ids, columns=currents
            ).to_numpy(float)
            auc = cell_auc(matrix, currents)
            observed, p_expected, p_two, _ = exact_cell_difference(auc, masks, labels, -1.0)
            rows.append({
                "surrogate_type": kind, "n_surrogates": count,
                "observed_raw_AUC_difference_SCA3_minus_WT": observed,
                "observed_expected_direction_AUC_difference_WT_minus_SCA3": -observed,
                "exact_p_AUC_expected_direction": p_expected,
                "exact_p_AUC_two_sided": p_two,
            })
    return pd.DataFrame(rows)


def _matrix_from_table(table, value_column, cell_ids, currents):
    return table.pivot(index="cell_id", columns="current_pA", values=value_column).reindex(
        index=cell_ids, columns=currents
    ).to_numpy(float)


def _primary_adjusted_analysis(config, upstream, domain, summary, cell_ids, labels, currents, masks):
    base_config = copy.deepcopy(dict(config))
    _, engine_cells, engine_labels, engine_currents, base_matrices = engine._prepare(base_config, upstream)
    if list(engine_cells) != list(cell_ids) or not np.array_equal(engine_labels, labels) or not np.array_equal(engine_currents, currents):
        raise RuntimeError("Cell/current order differs from the frozen v0.6.1 engine order.")
    required_covariates = [
        "group", "cell_id", "current_pA", "firing_rate_hz", "mean_isi_ms",
        "baseline_noise_mV", "stationary_samples",
    ]
    merged = domain[required_covariates].merge(
        summary[["group", "cell_id", "current_pA", "pi_shuffle_excess_nats"]],
        on=["group", "cell_id", "current_pA"], how="inner", validate="one_to_one",
    )
    order = pd.MultiIndex.from_product([cell_ids, currents], names=["cell_id", "current_pA"])
    merged = merged.set_index(["cell_id", "current_pA"]).reindex(order).reset_index()
    group_by_cell = domain[["cell_id", "group"]].drop_duplicates().set_index("cell_id").group
    merged["group"] = merged.cell_id.map(group_by_cell)
    merged = prepare_covariates(config, merged)
    outcome = "pi_shuffle_excess_nats"
    covariates = list(config["adjustment"]["covariates"])
    prediction, residual, coefficients = ridge_crossfit_by_cell_current(
        merged, cell_ids, currents, outcome, covariates,
        float(config["adjustment"]["ridge_lambda"]),
    )
    merged["crossfit_prediction"] = prediction
    merged["crossfit_residual"] = residual
    residual_matrix = _matrix_from_table(merged, "crossfit_residual", cell_ids, currents)
    matrices = dict(base_matrices)
    matrices["predictive_information_nats"] = residual_matrix
    cache = engine.SubsetReferenceCache(base_config, matrices, int((~labels).sum()))
    burden, excitability, predictive, strict, observed_index, observed = engine._exact_target_specific_scores(
        base_config, labels, currents, cache, masks
    )
    scores, references = engine._observed_score_tables(
        base_config, cell_ids, labels, currents, cache, observed
    )
    scores = scores.rename(columns={
        "predictive_information_nats__oriented_z": "surrogate_corrected_PI_residual__oriented_z",
        "predictive_dynamics_z": "surrogate_corrected_predictive_dynamics_z",
    })
    cell_summary = engine._cell_summary(base_config, scores.rename(columns={
        "surrogate_corrected_predictive_dynamics_z": "predictive_dynamics_z"
    }), cell_ids, labels, currents)
    current_tables = []
    curve_tables = []
    for analysis, values in (
        ("surrogate_corrected_independent_domain_burden", burden),
        ("surrogate_corrected_predictive_dynamics", predictive),
    ):
        differences = engine._group_differences(values, masks)
        current_table, curve_table = engine._exact_difference_tables(
            currents, differences, observed_index, analysis, 1.0
        )
        current_tables.append(current_table)
        curve_tables.append(pd.DataFrame([curve_table]))
    i_exit, _ = engine._exact_i_exit(base_config, currents, strict, masks, observed_index)
    i_exit["analysis"] = "secondary_surrogate_corrected_I_exit"
    return {
        "residuals": merged, "coefficients": coefficients,
        "scores": scores, "references": references, "cell_summary": cell_summary,
        "current_tests": pd.concat(current_tables, ignore_index=True),
        "curve_tests": pd.concat(curve_tables, ignore_index=True), "i_exit": i_exit,
    }


def _comparison(v070: Path, adjusted_curve_tests: pd.DataFrame, fixed_auc_tests: pd.DataFrame):
    old = pd.read_csv(v070 / "mode_group_curve_exact_tests.csv")
    old = old[old.adjustment_mode == "activity_technical"]
    rows = []
    for _, row in old.iterrows():
        rows.append({
            "version": "0.7.0", "analysis": str(row.analysis),
            "observed_expected_direction_curve_AUC_difference": float(row.observed_expected_direction_curve_auc_difference),
            "exact_p_expected_direction": float(row.exact_p_curve_auc_expected_direction),
            "primary_or_secondary": "reference",
        })
    for _, row in adjusted_curve_tests.iterrows():
        rows.append({
            "version": "0.7.1", "analysis": str(row.analysis),
            "observed_expected_direction_curve_AUC_difference": float(row.observed_expected_direction_curve_auc_difference),
            "exact_p_expected_direction": float(row.exact_p_curve_auc_expected_direction),
            "primary_or_secondary": "secondary_burden" if "burden" in str(row.analysis) else "secondary_adjusted_domain",
        })
    primary = fixed_auc_tests[fixed_auc_tests.metric == "pi_shuffle_excess_nats"].iloc[0]
    rows.append({
        "version": "0.7.1", "analysis": "primary_shuffle_corrected_PI_AUC",
        "observed_expected_direction_curve_AUC_difference": float(primary.observed_expected_direction_AUC_difference),
        "exact_p_expected_direction": float(primary.exact_p_AUC_expected_direction),
        "primary_or_secondary": "primary",
    })
    return pd.DataFrame(rows)


def _group_curves(summary, adjusted_scores):
    rows = []
    metrics = [
        "pi_observed_nats", "pi_shuffle_median_nats", "pi_shuffle_excess_nats",
        "pi_aaft_median_nats", "pi_aaft_excess_nats",
    ]
    for (group, current), sub in summary.groupby(["group", "current_pA"], sort=True):
        for metric in metrics:
            values = sub[metric].dropna()
            rows.append({
                "source": "raw_surrogate", "group": group, "current_pA": int(current),
                "metric": metric, "n": int(len(values)), "median": float(values.median()),
                "q25": float(values.quantile(.25)), "q75": float(values.quantile(.75)),
            })
    for (group, current), sub in adjusted_scores.groupby(["group", "current_pA"], sort=True):
        for metric in ["surrogate_corrected_predictive_dynamics_z", "independent_domain_burden"]:
            values = sub[metric].dropna()
            rows.append({
                "source": "target_specific_nested_LOO", "group": group,
                "current_pA": int(current), "metric": metric, "n": int(len(values)),
                "median": float(values.median()), "q25": float(values.quantile(.25)),
                "q75": float(values.quantile(.75)),
            })
    return pd.DataFrame(rows)


def _plot_band(ax, frame, metric, title, ylabel, groups):
    colors = {groups[0]: "#2878b5", groups[1]: "#be2334"}
    for group in groups:
        sub = frame[(frame.group == group) & (frame.metric == metric)].sort_values("current_pA")
        ax.plot(sub.current_pA, sub["median"], marker="o", color=colors[group], label=group)
        ax.fill_between(sub.current_pA, sub.q25, sub.q75, color=colors[group], alpha=.18)
    ax.set(title=title, xlabel="Injected current (pA)", ylabel=ylabel)
    ax.legend(frameon=False)


def _plots(output, config, curves, summary, auc_tests, current_tests, count_stability):
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    groups = [str(x) for x in config["analysis"]["groups"]]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _plot_band(axes[0], curves, "pi_observed_nats", "Observed ordinal predictive information", "PI (nats)", groups)
    _plot_band(axes[1], curves, "pi_shuffle_median_nats", "Shuffled finite-sample/overlap null", "surrogate PI (nats)", groups)
    fig.tight_layout(); fig.savefig(figures / "observed_and_shuffle_null.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _plot_band(axes[0], curves, "pi_shuffle_excess_nats", "Primary: PI excess over shuffled null", "PI excess (nats)", groups)
    _plot_band(axes[1], curves, "pi_aaft_excess_nats", "Sensitivity: PI excess over AAFT null", "PI excess (nats)", groups)
    for ax in axes:
        ax.axhline(0.0, color="black", ls="--", lw=1)
    fig.tight_layout(); fig.savefig(figures / "surrogate_corrected_PI_curves.png", dpi=180); plt.close(fig)

    plot = auc_tests.copy()
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    x = np.arange(len(plot))
    ax.bar(x, plot.observed_expected_direction_AUC_difference, color="#4c78a8")
    ax.axhline(0.0, color="black", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(plot.metric, rotation=25, ha="right")
    ax.set(ylabel="WT minus SCA3 cell-AUC difference", title="Exact cell-level surrogate sensitivity")
    for index, row in plot.iterrows():
        ax.text(index, row.observed_expected_direction_AUC_difference,
                " p={:.3g}".format(row.exact_p_AUC_expected_direction),
                rotation=90, va="bottom", ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(figures / "PI_AUC_exact_tests.png", dpi=180); plt.close(fig)

    primary_current = current_tests[current_tests.metric == "pi_shuffle_excess_nats"].sort_values("current_pA")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(primary_current.current_pA, primary_current.observed_expected_direction_difference, marker="o", color="#4c78a8")
    significant = primary_current.exact_p_expected_direction_maxT_adjusted <= .05
    ax.scatter(primary_current.loc[significant, "current_pA"], primary_current.loc[significant, "observed_expected_direction_difference"], color="#be2334", s=70, zorder=3)
    ax.axhline(0.0, color="black", lw=1)
    ax.set(xlabel="Injected current (pA)", ylabel="WT minus SCA3 difference",
           title="Primary shuffle-corrected PI by current")
    fig.tight_layout(); fig.savefig(figures / "primary_currentwise_effect.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for kind, sub in count_stability.groupby("surrogate_type"):
        sub = sub.sort_values("n_surrogates")
        ax.plot(sub.n_surrogates, sub.observed_expected_direction_AUC_difference_WT_minus_SCA3,
                marker="o", label=kind)
    ax.set(xlabel="Number of surrogates per sweep", ylabel="WT minus SCA3 AUC difference",
           title="Surrogate-count convergence")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(figures / "surrogate_count_stability.png", dpi=180); plt.close(fig)


def _write_summary(output, auc_tests, adjusted, reproduction_tolerance):
    primary = auc_tests[auc_tests.metric == "pi_shuffle_excess_nats"].iloc[0]
    aaft = auc_tests[auc_tests.metric == "pi_aaft_excess_nats"].iloc[0]
    burden = adjusted["curve_tests"][adjusted["curve_tests"].analysis == "surrogate_corrected_independent_domain_burden"].iloc[0]
    exits = adjusted["cell_summary"][~adjusted["cell_summary"].I_exit_calibrated_censored]
    exit_text = ", ".join("{}={} pA".format(row.cell_id, int(row.I_exit_calibrated_pA)) for _, row in exits.iterrows()) or "none"
    text = """# NeuroThermo v0.7.1 — surrogate validation of ordinal predictive information

## Frozen design

The analysis unit is one cell. The frozen cohort contains all 13 WT and 7 SCA3 cells. No animal-level inference is performed. The current grid is 100–600 pA. The primary null is value shuffling, which quantifies finite-sample and overlapping-window bias of ordinal predictive information. AAFT is a secondary null preserving the amplitude distribution and approximately the power spectrum.

## Primary endpoint

Primary metric: observed ordinal PI minus the median shuffled-surrogate PI. WT-minus-SCA3 cell-AUC difference: {primary_effect:.8g}; exact p={primary_p:.8g}.

## Secondary checks

AAFT-excess WT-minus-SCA3 cell-AUC difference: {aaft_effect:.8g}; exact p={aaft_p:.8g}.

Surrogate-corrected adjusted burden curve AUC difference: {burden_effect:.8g}; exact p={burden_p:.8g}.

Secondary adjusted exits: {exit_text}.

## Reproduction requirement

PI is recomputed from the raw stationary ABF interval with the frozen v0.3.1 estimator. The maximum permitted absolute discrepancy from v0.3.1 is {tolerance:.3g} nats. Failure of this check aborts the run.

## Interpretation boundary

The primary result tests whether the WT–SCA3 difference survives correction for estimator bias caused by finite record length and overlapping ordinal windows. AAFT sensitivity asks whether the difference remains beyond structure captured by the amplitude distribution and approximately by the spectrum. Neither result establishes disease time, irreversibility, causal mechanism, or a thermodynamic phase transition. I_exit remains a secondary current-stress threshold and is not a disease-onset time.
""".format(
        primary_effect=float(primary.observed_expected_direction_AUC_difference),
        primary_p=float(primary.exact_p_AUC_expected_direction),
        aaft_effect=float(aaft.observed_expected_direction_AUC_difference),
        aaft_p=float(aaft.exact_p_AUC_expected_direction),
        burden_effect=float(burden.observed_expected_direction_curve_auc_difference),
        burden_p=float(burden.exact_p_curve_auc_expected_direction),
        exit_text=exit_text, tolerance=float(reproduction_tolerance),
    )
    (output / "RUN_SUMMARY.md").write_text(text, encoding="utf-8")


def run_pipeline(config, upstream: Path, v070: Path, raw_root: Path, output: Path) -> None:
    checks = validate_inputs(config, upstream, v070, raw_root)
    if not bool(checks.passed.all()):
        raise ValueError("Input validation failed:\n" + checks.to_string(index=False))
    output.mkdir(parents=True, exist_ok=True)
    checks.to_csv(output / "input_validation.csv", index=False)

    domain = _read_domain(config, upstream)
    summary, long_table, arrays = _extract_and_compute_surrogates(config, domain, raw_root)
    tolerance = float(config["surrogates"]["reproduction_tolerance_nats"])
    reproduction = summary[[
        "group", "cell_id", "current_pA", "sweep_index", "source_path_upstream",
        "resolved_source_path", "predictive_information_upstream_nats", "pi_observed_nats",
        "absolute_reproduction_error_nats", "upstream_stationary_samples",
        "raw_stationary_samples", "resampled_stationary_samples", "target_dt_ms",
    ]].copy()
    reproduction["within_tolerance"] = reproduction.absolute_reproduction_error_nats <= tolerance
    reproduction.to_csv(output / "raw_trace_reproduction_audit.csv", index=False)
    if not bool(reproduction.within_tolerance.all()):
        failed = reproduction[~reproduction.within_tolerance]
        raise ValueError("Raw PI reproduction failed:\n" + failed.to_string(index=False))

    summary.to_csv(output / "surrogate_predictive_information_summary.csv", index=False)
    long_table.to_csv(output / "surrogate_predictive_information_long.csv.gz", index=False, compression="gzip")

    groups = [str(x) for x in config["analysis"]["groups"]]
    currents = np.asarray(config["analysis"]["currents_pA"], dtype=int)
    cell_ids, labels = _ordered_cells(domain, groups)
    masks = label_masks(len(cell_ids), int(labels.sum()), int(config["analysis"]["exact_max_labelings"]))
    metric_names = [
        "pi_observed_nats", "pi_shuffle_excess_nats", "pi_shuffle_z",
        "pi_aaft_excess_nats", "pi_aaft_z",
    ]
    auc_tests, current_tests = exact_fixed_metric_tests(
        summary, metric_names, cell_ids, labels, currents, masks, expected_direction=-1.0
    )
    auc_tests.to_csv(output / "PI_cell_AUC_exact_tests.csv", index=False)
    current_tests.to_csv(output / "PI_currentwise_exact_tests.csv", index=False)
    count_stability = _count_stability(
        config, summary, arrays, cell_ids, labels, currents, masks
    )
    count_stability.to_csv(output / "surrogate_count_stability.csv", index=False)

    adjusted = _primary_adjusted_analysis(
        config, upstream, domain, summary, cell_ids, labels, currents, masks
    )
    adjusted["residuals"].to_csv(output / "primary_crossfit_shuffle_excess_residuals.csv", index=False)
    adjusted["coefficients"].to_csv(output / "primary_crossfit_ridge_coefficients.csv", index=False)
    adjusted["scores"].to_csv(output / "primary_cell_current_scores.csv", index=False)
    adjusted["references"].to_csv(output / "primary_target_specific_references.csv", index=False)
    adjusted["cell_summary"].to_csv(output / "primary_cell_vulnerability_summary.csv", index=False)
    adjusted["current_tests"].to_csv(output / "primary_adjusted_currentwise_exact_tests.csv", index=False)
    adjusted["curve_tests"].to_csv(output / "primary_adjusted_group_curve_exact_tests.csv", index=False)
    adjusted["i_exit"].to_csv(output / "secondary_I_exit_exact_test.csv", index=False)

    curves = _group_curves(summary, adjusted["scores"])
    curves.to_csv(output / "group_dynamic_curves.csv", index=False)
    comparison = _comparison(v070, adjusted["curve_tests"], auc_tests)
    comparison.to_csv(output / "v070_v071_comparison.csv", index=False)

    rule = {
        "pipeline_version": __version__, "analysis_unit": "cell", "animal_level_inference": False,
        "cohort": {"WT": int((~labels).sum()), "SCA3": int(labels.sum())},
        "currents_pA": [int(x) for x in currents],
        "ordinal_PI": {
            "permutation_order": int(config["surrogates"]["permutation_order"]),
            "permutation_delay_ms": float(config["surrogates"]["permutation_delay_ms"]),
            "resample_dt_ms": float(config["surrogates"]["resample_dt_ms"]),
            "adjacent_code_pairs": True, "estimator": "empirical discrete mutual information",
        },
        "primary_null": {
            "type": "value_shuffle", "n_per_sweep": int(config["surrogates"]["n_shuffle"]),
            "metric": "observed_PI_minus_surrogate_median", "primary_endpoint": "cell_AUC exact label test",
        },
        "secondary_null": {
            "type": "AAFT", "n_per_sweep": int(config["surrogates"]["n_aaft"]),
            "metric": "observed_PI_minus_surrogate_median",
        },
        "adjusted_burden": {
            "outcome": "shuffle-corrected PI excess", "label_blind": True,
            "crossfit_unit": "entire cell", "fit_separately_by_current": True,
            "ridge_lambda": float(config["adjustment"]["ridge_lambda"]),
            "covariates": list(config["adjustment"]["covariates"]),
            "target_specific_nested_LOO_WT_calibration": True,
        },
        "exact_inference": {"cell_labelings": int(len(masks))},
        "I_exit_role": "secondary descriptive current-stress threshold only",
        "interpretation": "surrogate validation of ordinal PI; not disease time or thermodynamic phase transition",
    }
    _json_dump(output / "frozen_surrogate_and_inference_rule.json", rule)
    _plots(output, config, curves, summary, auc_tests, current_tests, count_stability)
    _write_summary(output, auc_tests, adjusted, tolerance)

    raw_paths = sorted(set(Path(x) for x in reproduction.resolved_source_path))
    input_hashes = {"v0.3.1/" + item: _sha256(upstream / item) for item in REQUIRED_UPSTREAM}
    input_hashes.update({"v0.7.0/" + item: _sha256(v070 / item) for item in REQUIRED_V070})
    raw_hashes = {str(path): _sha256(path) for path in raw_paths}
    output_files = sorted(x for x in output.rglob("*") if x.is_file() and x.name != "analysis_manifest.json")
    manifest = {
        "pipeline_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_minimum": "3.9", "upstream_dir": str(upstream),
        "v070_results_dir": str(v070), "raw_root": str(raw_root), "output_dir": str(output),
        "input_sha256": input_hashes, "raw_trace_sha256": raw_hashes,
        "output_sha256": {str(path.relative_to(output)): _sha256(path) for path in output_files},
        "n_cells": int(len(cell_ids)), "n_WT_cells": int((~labels).sum()),
        "n_SCA3_cells": int(labels.sum()), "n_currents": int(len(currents)),
        "n_exact_labelings": int(len(masks)), "animal_level_inference": False,
    }
    _json_dump(output / "analysis_manifest.json", manifest)
