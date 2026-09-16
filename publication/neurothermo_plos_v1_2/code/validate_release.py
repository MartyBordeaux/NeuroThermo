#!/usr/bin/env python3
"""Validate publication structure and the numerical claims used in the paper."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def require(relative: str) -> Path:
    path = ROOT / relative
    if not path.is_file():
        raise AssertionError(f"missing release file: {relative}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    provenance = json.loads(require("results/validation/PUBLICATION_INPUT_PROVENANCE.json").read_text())
    for relative, record in provenance["files"].items():
        path = require(relative)
        assert sha256(path) == record["sha256"], relative
        assert len(pd.read_csv(path)) == record["rows"], relative

    scenarios = pd.read_csv(require("data/inputs/transition_pair_scenarios.csv"))
    assert len(scenarios) == 264
    assert scenarios["biological_pair_key"].nunique() == 32

    mapping = pd.read_csv(require("data/animal_to_cell_mapping.csv"))
    sca = mapping[mapping["genotype"].eq("SCA3")]
    assert int(sca["in_raw_archive"].sum()) == 9
    assert int(sca["in_fitted_cohort"].sum()) == 7
    assert int(sca["in_primary_multisweep_cohort"].sum()) == 6
    assert set(sca.loc[sca["in_raw_archive"], "animal_id"]) == {"SCA3_DD20", "SCA3_DD24"}

    pairs = pd.read_csv(require("data/kl_convergence_v1_0_1/pair_markers_both_orders.csv"))
    primary = pairs[
        np.isclose(pairs["dt_ms"], 0.025)
        & pairs["marker_variant"].eq("seed_median_curve_isotonic")
        & pairs["view"].eq("xyz")
    ]
    expected = {"marker_first": (-0.12080968050421201, 26), "curve_first": (-0.16881793746734436, 27)}
    for order, (median_delta, n_negative) in expected.items():
        values = primary[primary["aggregation_order"].eq(order)]["kl_minus_firing_p"]
        assert len(values) == 32
        assert np.isclose(values.median(), median_delta, atol=1e-12)
        assert int((values < 0).sum()) == n_negative

    handoff = pd.read_csv(require("data/figure_source/fig3_combined_drive_handoff_summary.csv"))
    early = handoff[(handoff["path_segment"].eq("early")) & (handoff["reference_component"].eq("Applied J"))].iloc[0]
    late = handoff[(handoff["path_segment"].eq("late")) & (handoff["reference_component"].eq("kappa_I"))].iloc[0]
    assert early["pearson_curve_correlation"] > 0.99
    assert early["normalized_residual"] < 0.05
    assert late["pearson_curve_correlation"] > 0.95
    assert late["normalized_residual"] < 0.09

    verdict = json.loads(require("results/validation/nonequilibrium_geometry_v1_0_1/FORMALISM_VERDICT.json").read_text())
    assert verdict["stationary_formalism"] == "NESS"

    forbidden = tuple("/" + name + "/" for name in ("root", "workspace", "home"))
    text_suffixes = {".md", ".py", ".R", ".sh", ".yaml", ".yml", ".json", ".txt"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix in text_suffixes:
            content = path.read_text(encoding="utf-8", errors="strict")
            if any(token in content for token in forbidden):
                raise AssertionError(f"machine-specific path in {path.relative_to(ROOT)}")

    print("Publication release validation: PASS")


if __name__ == "__main__":
    main()
