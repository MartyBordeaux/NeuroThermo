#!/usr/bin/env python3
"""NeuroThermo publication-release preflight.

This script does not run scientific analyses. It verifies that a clean clone contains
all source packages and frozen assets required to start the documented workflow.
It intentionally fails when provenance-critical inputs are missing or ambiguous.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Check:
    kind: str
    name: str
    path: Path
    required: bool = True
    note: str = ""


PIPELINES = [
    "NeuroThermo_cell_fit_v3_9",
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

# Canonical publication location. These are intentionally not server paths.
CAL = ROOT / "data" / "calibration"
CALIBRATION = [
    ("candidate_events_with_predictions.csv", "spike-event/QC table used by final cell fit"),
    ("frozen_accepted_spiking_sweeps_v3_5.csv", "publication-frozen accepted spiking sweeps"),
    ("frozen_peak_overrides_v3_5.csv", "publication-frozen manual peak overrides"),
    ("frozen_threshold_brackets_v3_5.csv", "publication-frozen rheobase brackets"),
    ("frozen_v3_1_cell_fit_summary.csv", "baseline cell-fit summary required by v3.9"),
    ("frozen_v3_1_sweep_fit_summary.csv", "baseline sweep-fit summary required by v3.9"),
    ("frozen_v3_1_identifiability.csv", "baseline identifiability table required by v3.9"),
    ("seed_cell_summary_v3_9.csv", "optimizer seed/baseline table required by v3.9"),
]

# These names appear in the historical server v3.9 config. Presence is not accepted
# silently because the publication specification freezes v3.5-named inputs.
V36_CONFLICT_NAMES = [
    "frozen_accepted_spiking_sweeps_v3_6.csv",
    "frozen_peak_overrides_v3_6.csv",
    "frozen_threshold_brackets_v3_6.csv",
]


def build_checks() -> list[Check]:
    checks: list[Check] = []
    for name in PIPELINES:
        checks.append(Check("pipeline", name, ROOT / "code" / "pipelines" / name))
    for name, note in CALIBRATION:
        checks.append(Check("calibration", name, CAL / name, True, note))
    checks.extend([
        Check("release", "figure source data", ROOT / "data" / "figure_source"),
        Check("release", "KL frozen results", ROOT / "data" / "kl_convergence_v1_0_1"),
        Check("release", "nonequilibrium frozen results", ROOT / "data" / "nonequilibrium_geometry_v1_0_1"),
        Check("release", "article figures", ROOT / "results" / "figures"),
        Check("provenance", "recovered source hashes", ROOT / "docs" / "RECOVERED_SOURCE_HASHES.tsv"),
        Check("provenance", "computational DAG audit", ROOT / "docs" / "COMPUTATIONAL_DAG_AUDIT_2026-08-31.md"),
    ])
    return checks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true", help="return non-zero if any required item is missing")
    args = ap.parse_args()

    checks = build_checks()
    missing: list[Check] = []
    print(f"NeuroThermo publication preflight\nrelease_root={ROOT}\n")
    for c in checks:
        ok = c.path.exists()
        status = "PASS" if ok else "MISSING"
        print(f"{status:7s}  {c.kind:11s}  {c.name}  ->  {c.path.relative_to(ROOT)}")
        if not ok and c.required:
            missing.append(c)

    conflicts = [CAL / name for name in V36_CONFLICT_NAMES if (CAL / name).exists()]
    print("\nProvenance checks:")
    if conflicts:
        print("BLOCKED  v3.6-named calibration files are present while publication provenance specifies v3.5 names:")
        for p in conflicts:
            print(f"         {p.relative_to(ROOT)}")
        print("         Do not substitute or rename until content hashes/provenance are reconciled.")
    else:
        print("PASS     no silent v3.5/v3.6 substitution detected")

    print("\nSummary:")
    print(f"required_missing={len(missing)}")
    print(f"provenance_conflict={int(bool(conflicts))}")
    if missing:
        print("Upstream clean-clone rerun is NOT release-ready.")
        print("Missing required items are listed above; no path editing or manual substitution should be used as a workaround.")
    elif conflicts:
        print("Upstream clean-clone rerun is BLOCKED by unresolved calibration provenance.")
    else:
        print("Static release preflight PASS. Scientific smoke/full runs must still be executed separately.")

    return 1 if args.strict and (missing or conflicts) else 0


if __name__ == "__main__":
    raise SystemExit(main())
