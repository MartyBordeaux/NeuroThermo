#!/usr/bin/env python3
"""Diagnose row-level differences between recomputed and frozen stage-1 QC tables.

This script is intentionally diagnostic-only: it never changes the frozen reference
or the recomputed table and exits successfully after reporting differences. The
release preflight performs the strict assertions separately.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_KEY = ["group", "cell_id", "sweep_index", "peak_index"]
DECISION_COLUMNS = [
    "algorithm_detected",
    "fixed_qc_detected",
    "detected",
    "qc_changed",
    "qc_action",
    "qc_note",
]


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    return series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--recomputed", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--probability-threshold", type=float, default=None)
    args = parser.parse_args()

    ref = pd.read_csv(args.reference)
    new = pd.read_csv(args.recomputed)

    missing_key = [c for c in DEFAULT_KEY if c not in ref.columns or c not in new.columns]
    if missing_key:
        raise SystemExit(f"identity columns missing from comparison tables: {missing_key}")

    for label, frame in [("reference", ref), ("recomputed", new)]:
        dup = frame.duplicated(DEFAULT_KEY, keep=False)
        if dup.any():
            print(f"{label}_duplicate_identity_rows={int(dup.sum())}")
            print(frame.loc[dup, DEFAULT_KEY].sort_values(DEFAULT_KEY).to_string(index=False))
            raise SystemExit(f"duplicate candidate identities in {label}")

    print(f"reference_rows={len(ref)}")
    print(f"recomputed_rows={len(new)}")
    for column in ["algorithm_detected", "fixed_qc_detected", "detected", "qc_changed"]:
        if column in ref.columns and column in new.columns:
            print(
                f"{column}: reference={int(as_bool(ref[column]).sum())} "
                f"recomputed={int(as_bool(new[column]).sum())} "
                f"delta={int(as_bool(new[column]).sum() - as_bool(ref[column]).sum())}"
            )

    compare_columns = [c for c in DECISION_COLUMNS if c in ref.columns and c in new.columns]
    optional = [
        c
        for c in ["spike_probability", "time_ms", "current_pA", "peak_voltage_mV", "prominence_mV"]
        if c in ref.columns and c in new.columns
    ]
    merged = ref[DEFAULT_KEY + compare_columns + optional].merge(
        new[DEFAULT_KEY + compare_columns + optional],
        on=DEFAULT_KEY,
        how="outer",
        suffixes=("_ref", "_new"),
        indicator=True,
        validate="one_to_one",
    )

    missing = merged[merged["_merge"] == "left_only"]
    extra = merged[merged["_merge"] == "right_only"]
    print(f"missing_candidate_identities={len(missing)}")
    print(f"extra_candidate_identities={len(extra)}")
    if len(missing):
        print("MISSING IDENTITIES")
        print(missing[DEFAULT_KEY].to_string(index=False))
    if len(extra):
        print("EXTRA IDENTITIES")
        print(extra[DEFAULT_KEY].to_string(index=False))

    both = merged[merged["_merge"] == "both"].copy()
    mismatch = np.zeros(len(both), dtype=bool)
    decision_mismatch_columns = []
    for column in compare_columns:
        left = both[f"{column}_ref"]
        right = both[f"{column}_new"]
        if column in {"algorithm_detected", "fixed_qc_detected", "detected", "qc_changed"}:
            left = as_bool(left)
            right = as_bool(right)
        else:
            left = left.fillna("").astype(str)
            right = right.fillna("").astype(str)
        bad = left.ne(right).to_numpy()
        if bad.any():
            decision_mismatch_columns.append((column, int(bad.sum())))
            mismatch |= bad

    print("decision_mismatch_counts=" + repr(dict(decision_mismatch_columns)))
    decision_diff = both.loc[mismatch].copy()
    print(f"decision_mismatch_rows={len(decision_diff)}")

    if "spike_probability_ref" in both.columns:
        p_ref = pd.to_numeric(both["spike_probability_ref"], errors="coerce").to_numpy(float)
        p_new = pd.to_numeric(both["spike_probability_new"], errors="coerce").to_numpy(float)
        pdiff = np.abs(p_new - p_ref)
        both["probability_abs_diff"] = pdiff
        finite = np.isfinite(pdiff)
        print(f"spike_probability_nonzero_differences={int(np.count_nonzero(pdiff[finite] > 0))}")
        print(f"spike_probability_gt_1e-12={int(np.count_nonzero(pdiff[finite] > 1e-12))}")
        print(f"spike_probability_max_abs_diff={float(np.nanmax(pdiff)) if finite.any() else float('nan')}")

        if args.probability_threshold is not None:
            threshold = float(args.probability_threshold)
            crossed = finite & ((p_ref >= threshold) != (p_new >= threshold))
            print(f"probability_threshold_crossings_at_{threshold:g}={int(crossed.sum())}")
            if crossed.any():
                cols = DEFAULT_KEY + ["spike_probability_ref", "spike_probability_new", "probability_abs_diff"]
                print("THRESHOLD CROSSINGS")
                print(both.loc[crossed, cols].sort_values("probability_abs_diff", ascending=False).to_string(index=False))

        top = both.loc[finite & (pdiff > 0)].nlargest(30, "probability_abs_diff")
        if len(top):
            cols = DEFAULT_KEY + ["spike_probability_ref", "spike_probability_new", "probability_abs_diff"]
            print("TOP PROBABILITY DIFFERENCES")
            print(top[cols].to_string(index=False))

    if len(decision_diff):
        show = DEFAULT_KEY.copy()
        for column in compare_columns:
            show.extend([f"{column}_ref", f"{column}_new"])
        if "spike_probability_ref" in decision_diff.columns:
            show.extend(["spike_probability_ref", "spike_probability_new"])
        print("DECISION MISMATCH ROWS")
        print(decision_diff[show].sort_values(DEFAULT_KEY).to_string(index=False))

        per_cell = decision_diff.groupby(["group", "cell_id"], as_index=False).size()
        print("DECISION MISMATCHES BY CELL")
        print(per_cell.to_string(index=False))

    output = args.output
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        report = pd.concat([missing, extra, decision_diff], ignore_index=True, sort=False)
        report.to_csv(output, index=False)
        print(f"wrote_mismatch_report={output}")


if __name__ == "__main__":
    main()
