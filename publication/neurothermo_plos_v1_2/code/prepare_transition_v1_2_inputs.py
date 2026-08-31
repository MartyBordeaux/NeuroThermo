#!/usr/bin/env python3
"""Assemble transition-v1.2 inputs from recomputed v1.1 plus immutable endpoint inputs.

The historical v1.2 frozen directory mixed immutable endpoint inputs with three products
of transition v1.1 and a compact boundary-definition table. This script reconstructs
that boundary between stages explicitly and verifies the recomputed v1.1 products
against the historical frozen references before v1.2 is allowed to run.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HIST = ROOT / "data" / "transition_v1_2_frozen"
V11 = ROOT / "results" / "recomputed" / "transition_v1_1"
OUT = ROOT / "results" / "recomputed" / "transition_v1_2_inputs"

STATIC = [
    "transition_ready_endpoint_support.csv",
    "endpoint_cells_full_observable.csv",
    "cell_q75_protocol_anchors.csv",
]
V11_PRODUCTS = [
    "transition_projection_reference_v1_1.csv",
    "transition_projection_transform_v1_1.csv",
    "PRIMARY_ISI_STAGING.csv",
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def compare_csv(new: Path, old: Path, atol: float = 1e-12) -> dict:
    a, b = pd.read_csv(new), pd.read_csv(old)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        raise RuntimeError(f"Schema/row mismatch: {new.name}")
    max_abs = 0.0
    for c in a.columns:
        an = pd.to_numeric(a[c], errors="coerce")
        bn = pd.to_numeric(b[c], errors="coerce")
        numeric = an.notna() | bn.notna()
        if numeric.any() and int(numeric.sum()) == len(a):
            av, bv = an.to_numpy(float), bn.to_numpy(float)
            if not np.allclose(av, bv, rtol=1e-12, atol=atol, equal_nan=True):
                raise RuntimeError(f"Numeric mismatch in {new.name}:{c}")
            max_abs = max(max_abs, float(np.nanmax(np.abs(av - bv))) if len(av) else 0.0)
        else:
            if not a[c].fillna("<NA>").astype(str).equals(b[c].fillna("<NA>").astype(str)):
                raise RuntimeError(f"Text mismatch in {new.name}:{c}")
    return {"file": new.name, "historical_sha256": sha256(old), "recomputed_sha256": sha256(new), "max_abs_numeric_diff": max_abs}


def boundary_table(ref: pd.DataFrame) -> pd.DataFrame:
    rows = []
    labels = {
        "isi_primary_v1_0_frozen": "ISI (primary; frozen v1.0 consistent reference)",
        "active_rate_experimental_v2_1": "active rate (secondary; experimental q=.75 reference)",
    }
    for key in ("isi_primary_v1_0_frozen", "active_rate_experimental_v2_1"):
        r = ref.loc[ref["projection"].eq(key)]
        if len(r) != 1:
            raise RuntimeError(f"Expected one reference row for {key}, found {len(r)}")
        r = r.iloc[0]
        label = labels[key]
        rows.extend([
            {"projection": label, "stage": "WT-exit", "A_threshold": float(r["wt_exit_A_threshold"]), "definition": "maximum A among core-secure WT endpoints"},
            {"projection": label, "stage": "balance", "A_threshold": 0.5, "definition": "midpoint between WT centroid (A=0) and SCA3 centroid (A=1)"},
            {"projection": label, "stage": "SCA3-entry", "A_threshold": float(r["sca3_entry_A_threshold"]), "definition": "minimum A among core-secure SCA3 endpoints"},
        ])
    return pd.DataFrame(rows)


def main() -> int:
    for p in [HIST, V11]:
        if not p.is_dir():
            raise FileNotFoundError(p)
    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir(parents=True)

    provenance = {"static_inputs": [], "v1_1_products": [], "generated": []}
    for name in STATIC:
        src = HIST / name
        if not src.is_file():
            raise FileNotFoundError(src)
        dst = OUT / name
        shutil.copy2(src, dst)
        provenance["static_inputs"].append({"file": name, "sha256": sha256(dst), "source": str(src.relative_to(ROOT))})

    for name in V11_PRODUCTS:
        src = V11 / name
        hist = HIST / name
        if not src.is_file() or not hist.is_file():
            raise FileNotFoundError(src if not src.is_file() else hist)
        check = compare_csv(src, hist)
        shutil.copy2(src, OUT / name)
        check["source"] = str(src.relative_to(ROOT))
        provenance["v1_1_products"].append(check)

    ref = pd.read_csv(OUT / "transition_projection_reference_v1_1.csv")
    boundaries = boundary_table(ref)
    bpath = OUT / "staging_boundary_definitions_v1_1.csv"
    boundaries.to_csv(bpath, index=False)
    check = compare_csv(bpath, HIST / bpath.name)
    check["source"] = "generated deterministically from recomputed transition_projection_reference_v1_1.csv"
    provenance["generated"].append(check)

    provenance["status"] = "PASS"
    provenance["output_dir"] = str(OUT.relative_to(ROOT))
    (OUT / "ASSEMBLY_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2))
    print("TRANSITION_V1_1_TO_V1_2_ASSEMBLY_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
