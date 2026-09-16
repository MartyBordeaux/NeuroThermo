#!/usr/bin/env python3
"""Build the legacy XLSX compatibility view for characterization from the canonical CSV animal map."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data" / "animal_id_recovery" / "accepted_cohort.csv"
OUT = ROOT / "results" / "recomputed" / "characterization_inputs" / "NeuroThermo_animal_id_recovery.xlsx"
META = OUT.with_suffix(".provenance.json")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not CSV.is_file():
        raise FileNotFoundError(CSV)
    df = pd.read_csv(CSV, dtype={"experiment_day_code": "string"}, keep_default_na=False)
    required = ["group", "cell_id", "animal_id", "animal_id_status", "experiment_day_code"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Animal map is missing columns: {missing}")
    if len(df) != 20:
        raise RuntimeError(f"Expected 20 accepted-cohort rows, found {len(df)}")
    if df[["group", "cell_id"]].duplicated().any():
        raise RuntimeError("Duplicate group/cell_id rows in accepted cohort map")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Accepted cohort", index=False)

    meta = {
        "source": str(CSV.relative_to(ROOT)),
        "source_sha256": sha256(CSV),
        "generated": str(OUT.relative_to(ROOT)),
        "generated_sha256": sha256(OUT),
        "sheet": "Accepted cohort",
        "rows": int(len(df)),
        "columns": list(df.columns),
        "role": "legacy input-format compatibility only; accepted_cohort.csv is the canonical scientific source",
    }
    META.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
