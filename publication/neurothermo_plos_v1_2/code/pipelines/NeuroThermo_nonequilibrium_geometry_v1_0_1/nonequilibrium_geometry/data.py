from __future__ import annotations

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


REQUIRED = (
    "transition_pair_scenarios.csv",
    "biological_pair_stage_summary_v1_1.csv",
    "PRIMARY_ISI_STAGING.csv",
    "staging_boundary_definitions_v1_1.csv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_frozen_dir(explicit) -> Path:
    """Resolve exactly one user-specified frozen directory; never scan parents."""
    if explicit is None or not str(explicit).strip():
        raise ValueError("An explicit --frozen-dir is required; automatic discovery is disabled.")
    root = Path(explicit).expanduser().resolve()
    missing = [name for name in REQUIRED if not (root / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing frozen input(s) in {root}: {', '.join(missing)}")
    return root


def load_frozen(root: Path):
    return tuple(pd.read_csv(root / name) for name in REQUIRED)


def select_scenarios(scenarios: pd.DataFrame, pair_stage: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    mode = str(cfg["cohort"]["mode"])
    coupled = pair_stage.loc[pair_stage["path_family"].eq("coupled")].copy()
    if mode == "core_secure_all_support":
        pairs = set(coupled.loc[coupled["both_core_secure"].astype(bool), "biological_pair_key"])
        out = scenarios.loc[scenarios["biological_pair_key"].isin(pairs)].copy()
    elif mode == "core_secure_best_only":
        pairs = set(coupled.loc[coupled["both_core_secure"].astype(bool), "biological_pair_key"])
        out = scenarios.loc[
            scenarios["biological_pair_key"].isin(pairs)
            & scenarios["wt_source"].eq("best")
            & scenarios["sca_source"].eq("best")
        ].copy()
    elif mode == "all_pairs_best_only":
        out = scenarios.loc[scenarios["wt_source"].eq("best") & scenarios["sca_source"].eq("best")].copy()
    elif mode == "all_support":
        out = scenarios.copy()
    else:
        raise ValueError(f"Unknown cohort.mode={mode}")
    requested = cfg["cohort"].get("scenario_ids")
    if requested is not None:
        requested = {int(value) for value in requested}
        out = out.loc[out["scenario_id"].astype(int).isin(requested)].copy()
        absent = sorted(requested - set(out["scenario_id"].astype(int)))
        if absent:
            raise ValueError(f"Requested scenario IDs not in selected cohort: {absent}")
    maximum = cfg["cohort"].get("max_scenarios")
    if maximum is not None:
        out = out.sort_values("scenario_id").head(int(maximum)).copy()
    if out.empty:
        raise ValueError("Scenario selection is empty.")
    support_sum = out.groupby("biological_pair_key")["within_pair_support_weight"].transform("sum")
    out["analysis_within_pair_weight"] = out["within_pair_support_weight"] / support_sum
    out["analysis_pair_weight"] = 1.0 / out["biological_pair_key"].nunique()
    out["analysis_scenario_weight"] = out["analysis_pair_weight"] * out["analysis_within_pair_weight"]
    return out.reset_index(drop=True)


def scenario_arrays(row, p_grid):
    p = np.asarray(p_grid, dtype=float)

    def linear(left, right):
        return (1.0 - p) * float(left) + p * float(right)

    def logarithmic(left, right):
        left, right = max(float(left), 1e-15), max(float(right), 1e-15)
        return np.exp((1.0 - p) * np.log(left) + p * np.log(right))

    theta = np.empty((len(p), 4), dtype=float)
    theta[:, 0] = linear(row.wt_b, row.sca_b)
    theta[:, 1] = logarithmic(row.wt_r, row.sca_r)
    theta[:, 2] = linear(row.wt_s, row.sca_s)
    theta[:, 3] = logarithmic(row.wt_kappa_I, row.sca_kappa_I)
    current = linear(row.wt_J_q75, row.sca_J_q75)
    return theta, current


def write_input_manifest(root: Path, output: Path) -> None:
    payload = {name: {"path": str(root / name), "sha256": sha256(root / name)} for name in REQUIRED}
    (output / "INPUT_MANIFEST.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
