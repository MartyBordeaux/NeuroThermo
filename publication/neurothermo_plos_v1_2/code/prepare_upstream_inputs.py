#!/usr/bin/env python3
"""Extract and verify imported v3.9 results and dynamic-v2.1 frozen inputs."""
from __future__ import annotations

import hashlib
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data" / "upstream_server_bundle" / "neurothermo_v39_dynamic_inputs_2026-08-31.tar.gz"
ARCHIVE_SHA256 = "c359ddea915621dbf2b31ce9cde2bdc1c0bf59d06f8fcfcbad8092f1f78f9128"

PREFIXES = {
    "neurothermo/v4/hr_cell_fit_v3_9/results_cellfit_v3_9/": ROOT / "data" / "v3_9_results_full",
    "neurothermo/v4/neurothermo_dynamic_v2_1/frozen/": ROOT / "data" / "dynamic_v2_1_frozen",
}

REQUIRED = {
    ROOT / "data" / "v3_9_results_full": [
        "cell_fit_summary.csv",
        "sweep_fit_summary.csv",
        "final_identifiability_alternatives.csv",
        "RUN_SUMMARY.json",
        "IDENTIFIABILITY_SUMMARY.json",
        "resolved_config.yaml",
    ],
    ROOT / "data" / "dynamic_v2_1_frozen": [
        "primary_cell_master.csv",
        "accepted_spiking_sweeps.csv",
        "selected_spike_events.csv",
        "threshold_brackets.csv",
        "final_identifiability_alternatives.csv",
        "animal_id_map.csv",
    ],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not ARCHIVE.exists():
        raise FileNotFoundError(ARCHIVE)
    got = sha256(ARCHIVE)
    if got != ARCHIVE_SHA256:
        raise RuntimeError(f"Upstream archive SHA-256 mismatch: {got} != {ARCHIVE_SHA256}")

    for dest in PREFIXES.values():
        shutil.rmtree(dest, ignore_errors=True)
        dest.mkdir(parents=True, exist_ok=True)

    extracted = {p: 0 for p in PREFIXES}
    with tarfile.open(ARCHIVE, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = member.name
            for prefix, dest in PREFIXES.items():
                if not name.startswith(prefix):
                    continue
                rel = Path(name[len(prefix):])
                if not rel.parts or ".." in rel.parts:
                    raise RuntimeError(f"Unsafe archive member: {name}")
                out = (dest / rel).resolve()
                if dest.resolve() not in out.parents:
                    raise RuntimeError(f"Unsafe extraction path: {name}")
                out.parent.mkdir(parents=True, exist_ok=True)
                src = tf.extractfile(member)
                if src is None:
                    raise RuntimeError(f"Could not read archive member: {name}")
                with out.open("wb") as f:
                    shutil.copyfileobj(src, f)
                extracted[prefix] += 1
                break

    errors = []
    for dest, names in REQUIRED.items():
        for name in names:
            if not (dest / name).is_file():
                errors.append(str((dest / name).relative_to(ROOT)))
    if errors:
        raise RuntimeError("Missing required extracted upstream files:\n" + "\n".join(errors))

    print(f"PASS upstream archive sha256={ARCHIVE_SHA256}")
    for prefix, count in extracted.items():
        print(f"PASS extracted {count} files from {prefix}")
    for dest, names in REQUIRED.items():
        for name in names:
            p = dest / name
            print(f"PASS {p.relative_to(ROOT)} {sha256(p)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
