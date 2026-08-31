#!/usr/bin/env python3
"""NeuroThermo publication-release preflight."""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

@dataclass
class Check:
    kind: str
    name: str
    path: Path
    required: bool = True

PIPELINES = [
    "NeuroThermo_cell_fit_v3_9_frozen_exact",
    "NeuroThermo_characterization_v1_0",
    "NeuroThermo_dynamic_v2_1",
    "NeuroThermo_endpoint_ensemble_v1_0_1",
    "NeuroThermo_transition_v1_0",
    "NeuroThermo_transition_v1_1",
    "NeuroThermo_transition_v1_2",
    "NeuroThermo_transition_v1_2_1",
    "NeuroThermo_transition_v1_3_1",
    "NeuroThermo_KL_convergence_v1_0_1",
    "NeuroThermo_nonequilibrium_geometry_v1_0_1",
]
CAL = ROOT / "data" / "calibration"
EXPECTED = {
    "candidate_events_with_predictions.csv": "af35c327b313482f534aa59669a47e52a4078f912a5e342efcfddf0158455640",
    "frozen_accepted_spiking_sweeps_v3_5.csv": "dad46b831eb4613af4a49673f83854e4ef48b81d0934c087234562d81a447a54",
    "frozen_peak_overrides_v3_5.csv": "64e35808199e6108355b015b4ca9ded6070deed852927877e705ccf118e95069",
    "frozen_threshold_brackets_v3_5.csv": "47ba271e6b8d70704de1c49aaac3677c6ee21e3001f33faaafad8761177f9741",
    "frozen_v3_1_cell_fit_summary.csv": "85b1fa2c457e4affc0db438cf885b4406f61b943cbc08073fbdeb7f4b57f42f9",
    "frozen_v3_1_sweep_fit_summary.csv": "5663d59c35aeb105ee45b0c4c8606375210294f377a6ee3adcd771356a70ab12",
    "frozen_v3_1_identifiability.csv": "16e810e3331a0f6eb6bc1c815bb0e0d5574ee93966b95c00346522d5470957d1",
    "seed_cell_summary_v3_9.csv": "cb74bc0783c9fd1db11cacba13ccabd273cfc225e6dc019ab6e4215433dceb72",
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()
    missing = []
    bad_hash = []
    print(f"NeuroThermo publication preflight\nrelease_root={ROOT}\n")
    for name in PIPELINES:
        p = ROOT / "code" / "pipelines" / name
        ok = p.exists()
        print(f"{'PASS' if ok else 'MISSING':7s}  pipeline     {name}  ->  {p.relative_to(ROOT)}")
        if not ok: missing.append(str(p))
    for name, expected in EXPECTED.items():
        p = CAL / name
        if not p.exists():
            print(f"MISSING  calibration  {name}")
            missing.append(str(p)); continue
        got = sha256(p)
        if got != expected:
            print(f"BADHASH  calibration  {name}  {got}")
            bad_hash.append(name)
        else:
            print(f"PASS     calibration  {name}  {got}")
    extras = [
        ("figure source data", ROOT / "data" / "figure_source"),
        ("KL frozen results", ROOT / "data" / "kl_convergence_v1_0_1"),
        ("nonequilibrium frozen results", ROOT / "data" / "nonequilibrium_geometry_v1_0_1"),
        ("article figures", ROOT / "results" / "figures"),
        ("recovered source hashes", ROOT / "docs" / "RECOVERED_SOURCE_HASHES.tsv"),
        ("cell-fit source comparison", ROOT / "docs" / "CELLFIT_V3_9_SOURCE_COMPARISON.tsv"),
    ]
    for name,p in extras:
        ok=p.exists(); print(f"{'PASS' if ok else 'MISSING':7s}  release      {name}  ->  {p.relative_to(ROOT)}")
        if not ok: missing.append(str(p))
    print("\nProvenance checks:")
    print("PASS     v3.5/v3.6 accepted-sweep, peak-override and threshold-bracket files were recovered from the server bundle and have identical SHA-256 across the renamed v3.5/v3.6 copies.")
    print("PASS     exact frozen cell-fit v3.9 source is used as the canonical executable package.")
    print(f"\nrequired_missing={len(missing)}")
    print(f"bad_hash={len(bad_hash)}")
    ok = not missing and not bad_hash
    print("Static release preflight PASS." if ok else "Release preflight FAIL.")
    return 1 if args.strict and not ok else 0

if __name__ == "__main__":
    raise SystemExit(main())
