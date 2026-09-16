#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PACKAGE_DIR"
RESULTS_DIR="results_nonequilibrium_geometry_v1_0_1"
ARCHIVE="neurothermo_nonequilibrium_geometry_results_v1_0_1.zip"

test -s "$RESULTS_DIR/RUN_SUMMARY.json"
test -s "$RESULTS_DIR/FORMALISM_VERDICT.json"
python - <<'PY'
import json
from pathlib import Path
import pandas as pd
summary = json.loads(Path("results_nonequilibrium_geometry_v1_0_1/RUN_SUMMARY.json").read_text())
if summary["status"] != "PASS":
    raise SystemExit(f"Refusing to pack incomplete run: {summary['status']}")
root = Path("results_nonequilibrium_geometry_v1_0_1")
required = [
    "animal_balanced_geometry.csv", "animal_pair_balanced_geometry.csv",
    "animal_mapping_used.csv", "protocol_performance_summary.csv",
    "preflight_endpoint_membership.csv", "PROTOCOL_VERDICT.json",
    "NUMERICAL_VALIDATION.json",
]
missing = [name for name in required if not (root / name).is_file()]
if missing:
    raise SystemExit("Refusing to pack missing v1.0.1 outputs: " + ", ".join(missing))
fluctuation = pd.read_csv(root / "fluctuation_relations.csv")
if not (fluctuation["n_unique_path_positions"] == 15).all():
    raise SystemExit("Protocol uniqueness gate failed")
cycles = pd.read_csv(root / "markov_cycle_affinities.csv")
if not {"scenario_id", "biological_pair_key", "seed", "p"}.issubset(cycles.columns):
    raise SystemExit("Cycle-affinity provenance gate failed")
PY

python - "$RESULTS_DIR" "$ARCHIVE" <<'PY'
from pathlib import Path
import sys
import zipfile

source, archive = Path(sys.argv[1]), Path(sys.argv[2])
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
    for path in sorted(source.rglob("*")):
        if path.is_file() and "checkpoints" not in path.parts:
            output.write(path, Path(source.name) / path.relative_to(source))
print(archive)
PY
