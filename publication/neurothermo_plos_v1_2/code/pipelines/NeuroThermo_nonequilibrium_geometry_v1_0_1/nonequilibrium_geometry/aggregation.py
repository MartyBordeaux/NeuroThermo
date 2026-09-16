from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PAIR_SEPARATOR = "__TO__"


def load_animal_mapping(cfg):
    mapping_cfg = cfg.get("animal_mapping", {})
    path_value = mapping_cfg.get("path")
    required = bool(mapping_cfg.get("required", False))
    if not path_value:
        if required:
            raise ValueError("animal_mapping.path is required")
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"Animal mapping not found: {path}")
        return None
    frame = pd.read_csv(path)
    expected = {"cell_id", "genotype", "animal_id", "provenance_status"}
    missing = sorted(expected - set(frame.columns))
    if missing:
        raise ValueError("Animal mapping is missing columns: " + ", ".join(missing))
    if frame["cell_id"].duplicated().any():
        duplicated = sorted(frame.loc[frame["cell_id"].duplicated(False), "cell_id"].unique())
        raise ValueError(f"Animal mapping has duplicated cell IDs: {duplicated}")
    return frame


def annotate_animal_pairs(frame, mapping, required=True):
    if mapping is None:
        if required:
            raise ValueError("Animal mapping is unavailable")
        return frame.copy()
    lookup = mapping.set_index("cell_id")["animal_id"].astype(str).to_dict()
    out = frame.copy()
    split = out["biological_pair_key"].astype(str).str.split(PAIR_SEPARATOR, n=1, expand=True)
    if split.shape[1] != 2:
        raise ValueError("biological_pair_key must contain exactly one __TO__ separator")
    out["wt_cell_id"], out["sca3_cell_id"] = split[0], split[1]
    out["wt_animal_id"] = out["wt_cell_id"].map(lookup)
    out["sca3_animal_id"] = out["sca3_cell_id"].map(lookup)
    missing = sorted(set(out.loc[out[["wt_animal_id", "sca3_animal_id"]].isna().any(axis=1), "biological_pair_key"]))
    if missing and required:
        raise ValueError("Animal mapping is incomplete for biological pairs: " + ", ".join(missing))
    out["animal_pair_key"] = out["wt_animal_id"].astype(str) + PAIR_SEPARATOR + out["sca3_animal_id"].astype(str)
    return out


def balanced_mean(frame, group_columns, weight_column=None, count_name=None):
    excluded = set(group_columns)
    if weight_column:
        excluded.add(weight_column)
    numeric = [column for column in frame.select_dtypes(include=[np.number]).columns if column not in excluded]
    rows = []
    for keys, group in frame.groupby(group_columns, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_columns, keys))
        if count_name:
            row[count_name] = len(group)
        if weight_column:
            weights = group[weight_column].to_numpy(float)
            weights = weights / weights.sum()
        else:
            weights = np.full(len(group), 1.0 / len(group))
        for column in numeric:
            values = group[column].to_numpy(float)
            valid = np.isfinite(values) & np.isfinite(weights)
            row[column] = float(np.average(values[valid], weights=weights[valid])) if valid.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def distribution_summary(frame, group_columns, count_columns=()):
    excluded = set(group_columns) | set(count_columns)
    numeric = [column for column in frame.select_dtypes(include=[np.number]).columns if column not in excluded]
    rows = []
    for keys, group in frame.groupby(group_columns, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_columns, keys))
        for column in count_columns:
            row[column] = int(group[column].nunique())
        for column in numeric:
            values = group[column].replace([np.inf, -np.inf], np.nan).dropna()
            row[column + "_median"] = float(values.median()) if len(values) else np.nan
            row[column + "_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[column + "_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)

